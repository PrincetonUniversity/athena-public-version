/*
  mppnccombine - joins together netCDF data files representing a decomposed
                 domain into a unified netCDF file.  It was originally
                 designed to be used as a postprocessor for the parallel I/O
                 programming interface "mpp_io_mod"
                 (http://www.gfdl.noaa.gov/~vb/mpp_io.html) by V. Balaji.

  V2.1.7:  Added option to initialize output variables with a missing_value
           from the variables of the first input file as suggested by
           Martin Schmidt (martin.schmidt@io-warnemuende.de) and
           Franz Tauber (franz.tauber@io-warnemuende.de).
  V2.1.6:  Bug fixes for greater than 4GB record sizes.  Does not contain
           V2.1.5 modifications which were a special case.
  V2.1.5:  Supports running in an MPI environment.  Reads arguments from a
           configuration file instead of from the command line; this is needed
           to work around a bug in Cray's aprun.
  V2.1.4:  Fixed a bug with file cleanup and the debugging messages.
  V2.1.3:  Fixed a run-time segmentation fault with some compilers; changed
           ncinfile allocation in process_file function.
  V2.1.2:  Fixed a bug with input files that have decomposed dimensions
           defined after the variables that use them.
  V2.1.1:  Added option (-64) for creating output files with 64-bit offset
           format requiring NetCDF 3.6.x.
  V2.1: Added an option (-h) to pad the output file's header with empty space.
        Added an option (-e #) to specify an ending number to a range of input
        filename extensions.  It no longer aborts on missing input files, but
        gives error messages at the end of all the processing.
  V2.0: Substantial rewrite; memory buffering increases speed several times.
  V1.2: Added support for specifying the start number in filename extensions.
  V1.1.1: Added a fix for dimensions that are not also variables.
  V1.1: Changed loop order for increased I/O efficiency; records are now the
        innermost loop then the variables loop.
  V1.0: Original release.

  Written by Hans Vahlenkamp (Hans.Vahlenkamp@noaa.gov)
  Geophysical Fluid Dynamics Laboratory / NOAA
  Princeton Forrestal Campus
  Last Updated: 05/15/08
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <netcdf.h>

/* Information structure for a file */
struct fileinfo
  {
   int ncfid;  /* ID of the input netCDF file */
   int ndims;  /* Number of dimensions */
   int nvars;  /* Number of variables */
   int ngatts;  /* Number of global attributes */
   int recdim;  /* ID of the record dimensions */
   char varname[MAX_NC_VARS][MAX_NC_NAME];  /* Names of the variables */
   nc_type datatype[MAX_NC_VARS]; /* Data types of the variables */
   int varndims[MAX_NC_VARS];  /* Number of dimensions for each variable */
   int vardim[MAX_NC_VARS][MAX_NC_DIMS];  /* Dimensions for each variable */
   int natts[MAX_NC_VARS];  /* Number of attributes for each variable */
   unsigned char vardecomp[MAX_NC_VARS];  /* Is the variable decomposed */
   char dimname[MAX_NC_DIMS][MAX_NC_NAME];  /* Names of the dimensions */
   long dimsize[MAX_NC_DIMS];  /* Sizes of the dimensions (decomposed) */
   long dimfullsize[MAX_NC_DIMS];  /* Full sizes of the dimensions */
   long dimstart[MAX_NC_DIMS];  /* Start positions within full dimensions */
   long dimend[MAX_NC_DIMS];  /* End positions within full dimensions */
   unsigned char varmiss[MAX_NC_VARS];  /* Does variable have missing_value */
   unsigned char varmissval[MAX_NC_VARS][8];  /* missing_value per variable */
  };

/* Auxiliary function prototypes */
void usage();
int process_file(char *, unsigned char, struct fileinfo *, char *, int *,
                 int *, int, int, void *[], int, unsigned char,
                 unsigned char);
int process_vars(struct fileinfo *, struct fileinfo *, unsigned char, int *,
                 int, int, int, void *[], unsigned char, unsigned char);
int flush_decomp(struct fileinfo *, int, int, void *[], unsigned char);
void print_debug(struct fileinfo *, unsigned char);
char *nc_type_to_str(nc_type);


int main(int argc, char *argv[])
  {
   unsigned char verbose=0;  /* Print some progress information? */
   unsigned char appendnc=0;  /* Append to an existing netCDF file? */
   unsigned char removein=0;  /* Remove the ".####" decomposed input files? */
   int nstart=0;  /* PE number of the first input netCDF file */
   int nend=(-1);  /* PE number of the last input netCDF file */
   int headerpad=16384;  /* Additional padding at the end of the header */
   int format=NC_NOCLOBBER;  /* Format of new netCDF output file */
   unsigned char missing=0;  /* Initialize output variables with */
                             /* "missing_value" instead of 0 value? */
   int outputarg=(-1);  /* Argument # of the output netCDF file */
   int inputarg=(-1);  /* Argument # of first input netCDF file */
   struct stat statbuf;  /* Dummy structure for file-testing "stat" call */
   struct fileinfo *ncoutfile;  /* Information for the output file */
   char outfilename[2048], *strptr;  /* Name of the output netCDF file */
   int outlen;  /* Length of the output filename */
   char infilename[2048];  /* Name of an input file */
   unsigned char infileerror=0;  /* Errors reading an input file */
   unsigned char infileerrors=0;  /* Errors reading any input files */
   int nfiles=(-1);  /* Number of files in the decomposed domain */
   int a, f, r;  /* Loop variables */
   int nrecs=1;  /* Number of records in each decomposed file */
   void *varbuf[NC_MAX_VARS];  /* Buffers for decomposed variables */

   /* Check the command-line arguments */
   if (argc < 2)
     {
      usage(); return(1);
     }
   for (a=1; a < argc; a++)
     {
      if (!strcmp(argv[a],"-v")) verbose=1;
      else if (!strcmp(argv[a],"-vv")) verbose=2;  /* Hidden debug mode */
      else if (!strcmp(argv[a],"-a")) appendnc=1;
      else if (!strcmp(argv[a],"-r")) removein=1;
      else if (!strcmp(argv[a],"-n"))
        {
         a++;
         if (a < argc) nstart=atoi(argv[a]);
         else
           {
            usage(); return(1);
           }
        }
      else if (!strcmp(argv[a],"-e"))
        {
         a++;
         if (a < argc) nend=atoi(argv[a]);
         else
           {
            usage(); return(1);
           }
        }
      else if (!strcmp(argv[a],"-h"))
        {
         a++;
         if (a < argc) headerpad=atoi(argv[a]);
         else
           {
            usage(); return(1);
           }
        }
      else if (!strcmp(argv[a],"-64"))
        format=(NC_NOCLOBBER | NC_64BIT_OFFSET);
      else if (!strcmp(argv[a], "-n4"))
	format=(NC_NOCLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL);
      else if (!strcmp(argv[a],"-m")) missing=1;
      else
        {
         outputarg=a; break;
        }
     }
   if (outputarg==(-1))
     {
      usage(); return(1);
     }
   if (argc-1 > outputarg) inputarg=outputarg+1;
   sprintf(outfilename,argv[outputarg]); outlen=strlen(outfilename);
   if (outlen > 4)
     {
      strptr=outfilename+outlen-5;
      if (!strcmp(strptr,".0000")) outfilename[outlen-5]='\0';
     }

   /* Disable fatal returns from netCDF library functions */
   ncopts=0;

   /* Create a new netCDF output file */
   if ((ncoutfile=(struct fileinfo *)malloc(sizeof(struct fileinfo)))==NULL)
     {
      fprintf(stderr,"Error: cannot allocate enough memory!\n"); return(1);
     }
   if (!appendnc)
     {
      if (stat(outfilename,&statbuf)==0)
        {
         fprintf(stderr,"Error: output file seems to exist already!\n");
         free(ncoutfile); return(1);
        }
      if ((ncoutfile->ncfid=nccreate(outfilename,format))==(-1))
        {
         fprintf(stderr,"Error: cannot create the output netCDF file!\n");
         free(ncoutfile); return(1);
        }
      ncsetfill(ncoutfile->ncfid,NC_NOFILL);
     }
   /* Open an existing netCDF file for appending */
   else
     {
      if ((ncoutfile->ncfid=ncopen(outfilename,NC_WRITE))==(-1))
        {
         fprintf(stderr,"Error: cannot open the output netCDF file for appending!\n");
         free(ncoutfile); return(1);
        }
     }

   for (f=0; f < NC_MAX_VARS; f++)
     varbuf[f]=NULL;

   /* No input files are specified on the command-line */
   if (inputarg==(-1))
     {
      if (nend > -1)
        for (r=0; r < nrecs; r++)
          {
           if (verbose) printf("record = %d\n",r);
           f=0; 
           for (a=nstart; a <= nend; a++)
             {
              sprintf(infilename,"%s.%04d",outfilename,a);
              if (verbose)
                {
                 if (a==nstart && r==0) printf("  n files to go... ");
                 else printf("  %d files to go... ",nend-nstart+1-f);
                 printf("processing \"%s\"\n",infilename);
                }
              if (stat(infilename,&statbuf)!=0) continue;
              infileerror=process_file(infilename,appendnc,ncoutfile,
                                       outfilename,&nfiles,&nrecs,r,f,varbuf,
                                       headerpad,verbose,missing);
              if (infileerror) infileerrors=1;
              appendnc=1; f++;
              if (f==nfiles || a==nend)
                {
                 if (verbose > 1)
                   printf("  Write variables from previous %d files\n",f);
                 flush_decomp(ncoutfile,nfiles,r,varbuf,verbose);
                 break;
                }
             }
          }
      else
        {
         nend=nstart+1;
         for (r=0; r < nrecs; r++)
           {
            if (verbose) printf("record = %d\n",r);
            f=0;
            for (a=nstart; a < nend; a++)
              {
               sprintf(infilename,"%s.%04d",outfilename,a);
               if (verbose)
                 {
                  if (a==nstart && r==0) printf("  n files to go... ");
                  else printf("  %d files to go... ",nend-a);
                  printf("processing \"%s\"\n",infilename);
                 }
               infileerror=process_file(infilename,appendnc,ncoutfile,
                                        outfilename,&nfiles,&nrecs,r,f,varbuf,
                                        headerpad,verbose,missing);
               if (infileerror) infileerrors=1;
               if (a==nstart && nfiles > 0) nend=nstart+nfiles;
               appendnc=1; f++;
               if (f==nfiles || a==(nend-1))
                 {
                  if (verbose > 1)
                    printf("  Write variables from previous %d files\n",f);
                  flush_decomp(ncoutfile,nfiles,r,varbuf,verbose);
                  f=0; continue;
                 }
              }
           }
        }
     }
   /* Loop over all the specified input files */
   else
     for (r=0; r < nrecs; r++)
       {
        if (verbose) printf("record = %d\n",r);
        f=0;
        for (a=inputarg; a < argc; a++)
          {
           if (verbose)
             {
              if ((argc-a)==1) printf("  1 file to go... ");
              else printf("  %d files to go... ",argc-a);
              printf("processing \"%s\"\n",argv[a]);
             }
           infileerror=process_file(argv[a],appendnc,ncoutfile,
                                    outfilename,&nfiles,&nrecs,r,f,varbuf,
                                    headerpad,verbose,missing);
           if (infileerror) infileerrors=1;
           appendnc=1; f++;
           if (f==nfiles || a==(argc-1))
             {
              if (verbose > 1)
                printf("  Write variables from previous %d files\n",f);
              flush_decomp(ncoutfile,nfiles,r,varbuf,verbose);
              f=0; continue;
             }
          }
       }

   /* Clean up... return 1 on error, otherwise 0 */
   for (f=0; f < NC_MAX_VARS; f++)
     {
      if (varbuf[f]!=NULL) free(varbuf[f]);
     }
   ncclose(ncoutfile->ncfid); free(ncoutfile);
   if (!infileerrors)
     {
      if (removein)
        {
         /* No input files are specified on the command-line */
         if (inputarg==(-1))
           {
            f=0;
            for (a=nstart; a <= nend; a++)
              {
               if (++f > nfiles) break;
               sprintf(infilename,"%s.%04d",outfilename,a);
               if (verbose) printf("Removing \"%s\"\n",infilename);
               unlink(infilename);
              }
           }
         /* Loop over all the specified input files */
         else
           for (a=inputarg; a < argc; a++)
             {
              if (verbose) printf("Removing \"%s\"\n",argv[a]);
              unlink(argv[a]);
             }
        }
     }
   else
     fprintf(stderr,"Warning: output file may be incomplete!\n");
   return(infileerrors);
  }


/* Print the usage message for mppnccombine */
void usage()
  {
   printf("mppnccombine 2.1.7 - (written by Hans.Vahlenkamp@noaa.gov)\n\n");
   printf("Usage:  mppnccombine [-v] [-a] [-r] [-n #] [-e #] [-h #] [-64] [-m]\n");
   printf("                     output.nc [input ...]\n\n");
   printf("  -v    Print some progress information.\n");
   printf("  -a    Append to an existing netCDF file (not heavily tested...).\n");
   printf("  -r    Remove the \".####\" decomposed files after a successful run.\n");
   printf("  -n #  Input filename extensions start with number #### instead of 0000.\n");
   printf("  -e #  Ending number #### of a specified range of input filename extensions.\n");
   printf("        Files within the range do not have to be consecutively numbered.\n");
   printf("  -h #  Add a specified number of bytes of padding at the end of the header.\n");
   printf("  -64   Create netCDF output files with the 64-bit offset format.\n");
   printf("  -m    Initialize output variables with a \"missing_value\" from the variables\n");
   printf("        of the first input file instead of the default 0 value.\n\n");
   printf("mppnccombine joins together an arbitrary number of netCDF input files, each\n");
   printf("containing parts of a decomposed domain, into a unified netCDF output file.\n");
   printf("An output file must be specified and it is assumed to be the first filename\n");
   printf("argument.  If the output file already exists, then it will not be modified\n");
   printf("unless the option is chosen to append to it.  If no input files are specified\n");
   printf("then their names will be based on the name of the output file plus the default\n");
   printf("numeric extension \".0000\", which will increment by 1.  There is an option for\n");
   printf("starting the filename extensions with an arbitrary number instead of 0.  There\n");
   printf("is an option for specifying an end to the range of filename extension numbers;\n");
   printf("files within the range do not have to be consecutively numbered.  If input\n");
   printf("files are specified then names will be used verbatim.\n\n");
   printf("A value of 0 is returned if execution completed successfully; a value of 1\n");
   printf("otherwise.\n");
  }


/* Open an input file and get some information about it, define the   */
/* structure of the output file if necessary, prepare to copy all the */
/* variables at the current record to memory                          */
int process_file(char *ncname, unsigned char appendnc,
                 struct fileinfo *ncoutfile, char *outncname, int *nfiles,
                 int *nrecs, int r, int f, void *varbuf[], int headerpad,
                 unsigned char verbose, unsigned char missing)
  {
   struct fileinfo *ncinfile;  /* Information about an input netCDF file */
   int nfiles2;  /* Number of files in the decomposed domain */
   int d, v, n;  /* Loop variables */
   int dimid;  /* ID of a dimension */
   int decomp[4];  /* "domain_decomposition = #0, #1, #2, #3" attribute */
                   /*  #0 starting position of original dimension   */
                   /*  #1 ending position of original dimension     */
                   /*  #2 starting position of decomposed dimension */
                   /*  #3 ending position of decomposed dimension   */
   char attname[MAX_NC_NAME];  /* Name of a global or variable attribute */
   unsigned char ncinfileerror=0;  /* Were there any file errors? */

   /* Information for netCDF input file */
   if ((ncinfile=(struct fileinfo *)malloc(sizeof(struct fileinfo)))==NULL)
     {
      fprintf(stderr,"Error: cannot allocate enough memory!\n"); return(1);
     }

   /* Open an input netCDF file */
   if ((ncinfile->ncfid=ncopen(ncname,NC_NOWRITE))==(-1))
     {
      fprintf(stderr,"Error: cannot open input file \"%s\"\n",ncname);
      free(ncinfile); return(1);
     }

   /* Determine the number of files in the decomposed domain */
   if (ncattget(ncinfile->ncfid,NC_GLOBAL,"NumFilesInSet",
                (void *)&nfiles2)==(-1))
     {
      if (*nfiles==1)
        {
         fprintf(stderr,"Error: missing the \"NumFilesInSet\" global attribute!\n");
         return(1);
        }
      else if (*nfiles==(-1))
        {
         fprintf(stderr,"Warning: missing the \"NumFilesInSet\" global attribute.\n");
        }
     }
   *nfiles=nfiles2;

   /* Get some general information about the input netCDF file */
   if (ncinquire(ncinfile->ncfid,&(ncinfile->ndims),&(ncinfile->nvars),
                 &(ncinfile->ngatts),&(ncinfile->recdim))==(-1))
     {
      fprintf(stderr,"Error: cannot read the file's metadata!\n");
      ncclose(ncinfile->ncfid); free(ncinfile); return(1);
     }

   /* Get some information about the dimensions */
   for (d=0; d < ncinfile->ndims; d++)
     {
      if ((ncdiminq(ncinfile->ncfid,d,ncinfile->dimname[d],
                    &(ncinfile->dimsize[d])))==(-1))
        {
         fprintf(stderr,"Error: cannot read dimension #%d's metadata!\n",d);
         ncclose(ncinfile->ncfid); free(ncinfile); return(1);
        }
      ncinfile->dimfullsize[d]=ncinfile->dimsize[d];
      ncinfile->dimstart[d]=1; ncinfile->dimend[d]=(-1);
     }

   /* Save some information for the output file */
   if (r==0)
     {
      ncoutfile->nvars=ncinfile->nvars; ncoutfile->recdim=ncinfile->recdim;
     }

   /* Get some information about the variables */
   for (v=0; v < ncinfile->nvars; v++)
     {
      if ((ncvarinq(ncinfile->ncfid,v,ncinfile->varname[v],
                    &(ncinfile->datatype[v]),&(ncinfile->varndims[v]),
                    ncinfile->vardim[v],&(ncinfile->natts[v])))==(-1))
        {
         fprintf(stderr,"Error: cannot read variable #%d's metadata!\n",v);
         ncclose(ncinfile->ncfid); free(ncinfile); return(1);
        }

      /* If the variable is also a dimension then get decomposition info */
      if ((dimid=ncdimid(ncinfile->ncfid,ncinfile->varname[v]))!=(-1))
        {
         if (ncattget(ncinfile->ncfid,v,"domain_decomposition",
             (void *)decomp)!=(-1))
           {
            ncinfile->dimfullsize[dimid]=decomp[1]-decomp[0]+1;
            ncinfile->dimstart[dimid]=decomp[2]-(decomp[0]-1);
            ncinfile->dimend[dimid]=decomp[3]-(decomp[0]-1);
           }
         else
           {
            ncinfile->dimfullsize[dimid]=ncinfile->dimsize[dimid];
            ncinfile->dimstart[dimid]=1; ncinfile->dimend[dimid]=(-1);
           }
        }
     }

   /* Get some additional information about the variables */
   for (v=0; v < ncinfile->nvars; v++)
     {
      /* Does the variable have a decomposed dimension? */
      ncinfile->vardecomp[v]=0;
      for (d=0; d < ncinfile->varndims[v]; d++)
        {
         if (ncinfile->dimend[ncinfile->vardim[v][d]]!=(-1))
           {
            ncinfile->vardecomp[v]=1; break;
           }
        }

      /* Save some information for the output file */
      if (r==0)
        {
         ncoutfile->varndims[v]=ncinfile->varndims[v];
         for (d=0; d < ncinfile->ndims; d++)
           ncoutfile->dimfullsize[d]=ncinfile->dimfullsize[d];
         for (d=0; d < ncinfile->varndims[v]; d++)
           ncoutfile->vardim[v][d]=ncinfile->vardim[v][d];
         ncoutfile->vardecomp[v]=ncinfile->vardecomp[v];
         strcpy(ncoutfile->varname[v],ncinfile->varname[v]);
         ncoutfile->varmiss[v]=0;
        }
     }

   /* If the output netCDF file was just created then define its structure */
   if (!appendnc)
     {
      if (verbose) printf("    Creating output \"%s\"\n",outncname);

      /* Define the dimensions */
      for (d=0; d < ncinfile->ndims; d++)
        {
         if (d==ncinfile->recdim)
           ncdimdef(ncoutfile->ncfid,ncinfile->dimname[d],NC_UNLIMITED);
         else ncdimdef(ncoutfile->ncfid,ncinfile->dimname[d],
                       ncinfile->dimfullsize[d]);
        }

      /* Define the variables and copy their attributes */
      for (v=0; v < ncinfile->nvars; v++)
        {
         ncvardef(ncoutfile->ncfid,ncinfile->varname[v],ncinfile->datatype[v],
                  ncinfile->varndims[v],ncinfile->vardim[v]);
         for (n=0; n < ncinfile->natts[v]; n++)
           {
            ncattname(ncinfile->ncfid,v,n,attname);
            if (missing)
              {
               if (!strcmp(attname,"missing_value"))
                 {
                  ncoutfile->varmiss[v]=1;
                  ncattget(ncinfile->ncfid,v,"missing_value",
                           (void *)(ncoutfile->varmissval[v]));
                 }
              }
            if (!strcmp(attname,"domain_decomposition")) continue;
            else
              {
               if (ncattcopy(ncinfile->ncfid,v,attname,ncoutfile->ncfid,v)==(-1))
                 {
                  fprintf(stderr,"Error: cannot copy variable \"%s\"'s attributes!\n",
                          ncinfile->varname[v]);
                  free(ncinfile); return(1);
                 }
              }
           }
        }

      /* Copy the global attributes */
      for (n=0; n < ncinfile->ngatts; n++)
        {
         ncattname(ncinfile->ncfid,NC_GLOBAL,n,attname);
         if (!strcmp(attname,"NumFilesInSet")) continue;
         else if (!strcmp(attname,"filename"))
           ncattput(ncoutfile->ncfid,NC_GLOBAL,attname,NC_CHAR,
                    strlen(outncname),(void *)outncname);
         else
           {
            if (ncattcopy(ncinfile->ncfid,NC_GLOBAL,attname,ncoutfile->ncfid,
                          NC_GLOBAL)==(-1))
              {
               fprintf(stderr,"Error: cannot copy the file's global attributes!\n");
               return(1);
              }
           }
        }

      /* Definitions done */
      nc__enddef(ncoutfile->ncfid,headerpad,4,0,4);
     }

   /* Copy all data values of the dimensions and variables to memory */
   ncinfileerror=process_vars(ncinfile,ncoutfile,appendnc,nrecs,r,*nfiles,
                              f,varbuf,verbose,missing);

   /* Done */
   ncclose(ncinfile->ncfid); free(ncinfile); return(ncinfileerror);
  }


/* Copy all data values in an input file at the current record to memory */
int process_vars(struct fileinfo *ncinfile, struct fileinfo *ncoutfile,
                 unsigned char appendnc, int *nrecs, int r, int nfiles,
                 int f, void *varbuf[], unsigned char verbose,
                 unsigned char missing)
  {
   int v, d, i, j, k, l, b, s;  /* Loop variables */
   int dimid;  /* ID of a dimension */
   void *values;  /* Current data values */
   long instart[MAX_NC_DIMS], outstart[MAX_NC_DIMS];  /* Data array sizes */
   long count[MAX_NC_DIMS];                           /*        "         */
   long long recsize;  /* Decomposed size of one record of a variable */
   long long recfullsize;  /* Non-decomposed size of one record of a variable */
   int varrecdim;  /* Variable's record dimension */
   static unsigned char first=1;  /* First time reading variables? */
   int imax, jmax, kmax, lmax;
   int imaxfull, jmaxfull, kmaxfull, lmaxfull;
   int imaxjmaxfull, imaxjmaxkmaxfull;
   int offset, ioffset, joffset, koffset, loffset;
   long long varbufsize;

   /* Check the number of records */
   if (*nrecs==1) *nrecs=ncinfile->dimsize[ncinfile->recdim];
   else
     if (ncinfile->dimsize[ncinfile->recdim] != *nrecs)
       {
        fprintf(stderr,"Error: different number of records than the first input file!\n");
        return(1);
       }

   /* Loop over all the variables */
   for (v=0; v < ncinfile->nvars; v++)
     {
      if (verbose > 1) printf("    variable = %s\n",ncinfile->varname[v]);

      /* Get read/write dimension sizes for the variable */
      recsize=1; recfullsize=1; varrecdim=(-1);
      outstart[0]=0; outstart[1]=0; outstart[2]=0; outstart[3]=0;
      for (d=0; d < ncinfile->varndims[v]; d++)
        {
         if (ncinfile->vardim[v][d]==ncinfile->recdim)
           {
            count[d]=1; varrecdim=d;
           }
         else
           {
            count[d]=ncinfile->dimsize[ncinfile->vardim[v][d]];
            recsize*=count[d]; instart[d]=0;
            outstart[d]=ncinfile->dimstart[ncinfile->vardim[v][d]]-1;
            recfullsize*=ncinfile->dimfullsize[ncinfile->vardim[v][d]];
           }
         if (verbose > 1)
           printf("      dim %d:  instart=%ld  outstart=%ld  count=%ld\n",d,
                  instart[d],outstart[d],count[d]);
        }

      /* Prevent unnecessary reads/writes */
      if (r > 0)
        {
         /* Prevent unnecessary reads/writes of the dimensions */
         if ((dimid=ncdimid(ncinfile->ncfid,ncinfile->varname[v]))!=(-1))
           {
            if (ncinfile->recdim==dimid)
              {
               if (f!=0) continue;
              }
            else continue;
           }
         /* Prevent unnecessary reads/writes of the variables */
         else
           {
            /* Prevent unnecessary reads/writes of non-decomposed variables
            if (ncinfile->vardecomp[v]!=1 && appendnc) continue; */

            /* Non-record variables */
            if (varrecdim==(-1)) continue;

            /* Non-decomposed record variables */
            if (ncinfile->vardecomp[v]!=1 && f > 0) continue;
           }
        }
      else
        {
         if (ncinfile->vardecomp[v]!=1 && appendnc) continue;
        }

      /* Allocate a buffer for the variable's record */
      if ((values=malloc(nctypelen(ncinfile->datatype[v])*recsize))==NULL)
        {
         fprintf(stderr,"Error: cannot allocate %lld bytes for decomposed variable \"%s\"'s values!\n",
                 nctypelen(ncinfile->datatype[v])*recsize,ncinfile->varname[v]);
         return(1);
        }

      /* Read the variable */
      if (varrecdim!=(-1)) instart[varrecdim]=outstart[varrecdim]=r;
      if (ncvarget(ncinfile->ncfid,v,instart,count,values)==(-1))
        {
         fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n",
                 ncinfile->varname[v]);
         return(1);
        }

      /* Write the buffered variable immediately if it's not decomposed */
      if (ncinfile->vardecomp[v]!=1)
        {
         if (verbose > 1)
           printf("      writing %lld bytes to file\n",
                  nctypelen(ncinfile->datatype[v])*recsize);
         if (ncvarput(ncoutfile->ncfid,v,outstart,count,values)==(-1))
           {
            fprintf(stderr,"Error: cannot write variable \"%s\"'s values!\n",
                    ncinfile->varname[v]);
            return(1);
           }
        }
      /* Save the buffer */
      else
        {
         /* Allocate a buffer for the variable's non-decomposed record size */
         if (first)
           {
            varbufsize=nctypelen(ncinfile->datatype[v])*recfullsize;
            if (verbose > 1)
              printf("      allocating %lld bytes for full domain\n",
                     varbufsize);
            if ((varbuf[v]=calloc(varbufsize,1))==NULL)
              {
               fprintf(stderr,"Error: cannot allocate %lld bytes for entire variable \"%s\"'s values!\n",
                       varbufsize,ncinfile->varname[v]); return(1);
              }
            if (missing && ncoutfile->varmiss[v])
              switch (ncinfile->datatype[v])
                {
                 case NC_BYTE:
                 case NC_CHAR:
                   for (s=0; s < recfullsize; s++)
                     *((unsigned char *)(varbuf[v])+s)=
                     *((unsigned char *)(ncoutfile->varmissval[v]));
                   break;
                 case NC_SHORT:
                   for (s=0; s < recfullsize; s++)
                     *((short *)(varbuf[v])+s)=
                     *((short *)(ncoutfile->varmissval[v]));
                   break;
                 case NC_INT:
                   for (s=0; s < recfullsize; s++)
                     *((int *)(varbuf[v])+s)=
                     *((int *)(ncoutfile->varmissval[v]));
                   break;
                 case NC_FLOAT:
                   for (s=0; s < recfullsize; s++)
                     *((float *)(varbuf[v])+s)=
                     *((float *)(ncoutfile->varmissval[v]));
                   break;
                 case NC_DOUBLE:
                   for (s=0; s < recfullsize; s++)
                     *((double *)(varbuf[v])+s)=
                     *((double *)(ncoutfile->varmissval[v]));
                   break;
                }
           }
         if (varbuf[v]==NULL)
           {
            fprintf(stderr,"Internal memory usage error!\n"); return(1);
           }
         if (verbose > 1)
           printf("      writing %lld bytes to memory\n",
                   nctypelen(ncinfile->datatype[v])*recsize);

         imax=ncinfile->dimsize[ncinfile->vardim[v][ncinfile->varndims[v]-1]];
         if (ncinfile->varndims[v] > 1)
           {
            dimid=ncinfile->vardim[v][ncinfile->varndims[v]-2];
            if (dimid==ncinfile->recdim) jmax=1;
            else jmax=ncinfile->dimsize[dimid];
           }
         else jmax=1;
         if (ncinfile->varndims[v] > 2)
           {
            dimid=ncinfile->vardim[v][ncinfile->varndims[v]-3];
            if (dimid==ncinfile->recdim) kmax=1;
            else kmax=ncinfile->dimsize[dimid];
           }
         else kmax=1;
         if (ncinfile->varndims[v] > 3)
           {
            dimid=ncinfile->vardim[v][ncinfile->varndims[v]-4];
            if (dimid==ncinfile->recdim) lmax=1;
            else lmax=ncinfile->dimsize[dimid];
           }
         else lmax=1;
         if (verbose > 1)
           printf("      imax=%d  jmax=%d  kmax=%d  lmax=%d\n",imax,jmax,
                  kmax,lmax);

         imaxfull=ncinfile->dimfullsize[ncinfile->vardim[v][ncinfile->varndims[v]-1]];
         if (ncinfile->varndims[v] > 1)
           jmaxfull=ncinfile->dimfullsize[ncinfile->vardim[v][ncinfile->varndims[v]-2]];
         else jmaxfull=1;
         if (ncinfile->varndims[v] > 2)
           kmaxfull=ncinfile->dimfullsize[ncinfile->vardim[v][ncinfile->varndims[v]-3]];
         else kmaxfull=1;
         if (ncinfile->varndims[v] > 3)
           {
            if (ncinfile->vardim[v][ncinfile->varndims[v]-4]!=ncinfile->recdim)
              lmaxfull=ncinfile->dimfullsize[ncinfile->vardim[v][ncinfile->varndims[v]-4]];
            else lmaxfull=1;
           }
         else lmaxfull=1;
         if (verbose > 1)
           printf("      imaxfull=%d  jmaxfull=%d  kmaxfull=%d  lmaxfull=%d\n",
                  imaxfull,jmaxfull,kmaxfull,lmaxfull);
         imaxjmaxfull=imaxfull*jmaxfull;
         imaxjmaxkmaxfull=imaxfull*jmaxfull*kmaxfull;

         ioffset=outstart[ncinfile->varndims[v]-0-1];
         if (ncinfile->varndims[v] > 1)
           joffset=outstart[ncinfile->varndims[v]-1-1];
         else joffset=0;
         if (ncinfile->varndims[v] > 2)
           koffset=outstart[ncinfile->varndims[v]-2-1];
         else koffset=0;
         if (ncinfile->varndims[v] > 3)
           loffset=outstart[ncinfile->varndims[v]-3-1];
         else loffset=0;
         if (varrecdim!=(-1))
           {
            switch (ncinfile->varndims[v])
              {
               case 1:
                 ioffset=0;
                 break;
               case 2:
                 joffset=0;
                 break;
               case 3:
                 koffset=0;
                 break;
               case 4:
                 loffset=0;
                 break;
              }
           }
         if (verbose > 1)
           printf("      ioffset=%d  joffset=%d  koffset=%d  loffset=%d\n",
                  ioffset,joffset,koffset,loffset);
         switch (ncinfile->datatype[v])
           {
            case NC_BYTE:
            case NC_CHAR:
              if (verbose > 1) printf("      start copying byte/char\n");
              b=0;
              for (l=0; l < lmax; l++)
                for (k=0; k < kmax; k++)
                  for (j=0; j < jmax; j++)
                    for (i=0; i < imax; i++)
                      {
                       offset=(i+ioffset)+
                              (j+joffset)*imaxfull+
                              (k+koffset)*imaxjmaxfull+
                              (l+loffset)*imaxjmaxkmaxfull;
                       *((unsigned char *)(varbuf[v])+offset)=
                       *((unsigned char *)values+(b++));
                      }
              if (verbose > 1) printf("      end copying byte/char\n");
              break;
            case NC_SHORT:
              if (verbose > 1) printf("      start copying short\n");
              b=0;
              for (l=0; l < lmax; l++)
                for (k=0; k < kmax; k++)
                  for (j=0; j < jmax; j++)
                    for (i=0; i < imax; i++)
                      {
                       offset=(i+ioffset)+
                              (j+joffset)*imaxfull+
                              (k+koffset)*imaxjmaxfull+
                              (l+loffset)*imaxjmaxkmaxfull;
                       *((short *)(varbuf[v])+offset)=
                       *((short *)values+(b++));
                      }
              if (verbose > 1) printf("      end copying short\n");
              break;
            case NC_INT:
              if (verbose > 1) printf("      start copying int\n");
              b=0;
              for (l=0; l < lmax; l++)
                for (k=0; k < kmax; k++)
                  for (j=0; j < jmax; j++)
                    for (i=0; i < imax; i++)
                      {
                       offset=(i+ioffset)+
                              (j+joffset)*imaxfull+
                              (k+koffset)*imaxjmaxfull+
                              (l+loffset)*imaxjmaxkmaxfull;
                       *((int *)(varbuf[v])+offset)=
                       *((int *)values+(b++));
                      }
              if (verbose > 1) printf("      end copying int\n");
              break;
            case NC_FLOAT:
              if (verbose > 1) printf("      start copying float\n");
              b=0;
              for (l=0; l < lmax; l++)
                for (k=0; k < kmax; k++)
                  for (j=0; j < jmax; j++)
                    for (i=0; i < imax; i++)
                      {
                       offset=(i+ioffset)+
                              (j+joffset)*imaxfull+
                              (k+koffset)*imaxjmaxfull+
                              (l+loffset)*imaxjmaxkmaxfull;
                       *((float *)(varbuf[v])+offset)=
                       *((float *)values+(b++));
                      }
              if (verbose > 1) printf("      end copying float\n");
              break;
            case NC_DOUBLE:
              if (verbose > 1) printf("      start copying double\n");
              b=0;
              for (l=0; l < lmax; l++)
                for (k=0; k < kmax; k++)
                  for (j=0; j < jmax; j++)
                    for (i=0; i < imax; i++)
                      {
                       offset=(i+ioffset)+
                              (j+joffset)*imaxfull+
                              (k+koffset)*imaxjmaxfull+
                              (l+loffset)*imaxjmaxkmaxfull;
                       *((double *)(varbuf[v])+offset)=
                       *((double *)values+(b++));
                      }
              if (verbose > 1) printf("      end copying double\n");
              break;
           }
        }

      /* Deallocate the decomposed variable's buffer */
      free(values);
     }
   first=0; return(0);
  }


/* Write all the buffered decomposed variables to the output file */
int flush_decomp(struct fileinfo *ncoutfile, int nfiles, int r,
                 void *varbuf[], unsigned char verbose)
  {
   int v, d;  /* Loop variable */
   long outstart[MAX_NC_DIMS];  /* Data array sizes */
   long count[MAX_NC_DIMS];     /*        "         */
   int varrecdim;  /* Position of a variable's record dimension */

   if (verbose > 1)
     {
      printf("    nvars=%d\n",ncoutfile->nvars);
     }

   /* Write out all the decomposed variables */
   for (v=0; v < ncoutfile->nvars; v++)
     {
      if (ncoutfile->vardecomp[v]==0) continue;
      if (verbose > 1) printf("    v=%d (%s)\n",v,ncoutfile->varname[v]);
      varrecdim=(-1);
      for (d=0; d < ncoutfile->varndims[v]; d++)
        {
         outstart[d]=0;
         if (ncoutfile->vardim[v][d]==ncoutfile->recdim)
           {
            count[d]=1; varrecdim=d;
           }
         else
           {
            count[d]=ncoutfile->dimfullsize[ncoutfile->vardim[v][d]];
           }
         if (verbose > 1)
           printf("      d=%d:  outstart=%ld  count=%ld\n",d,outstart[d],
                  count[d]);
        }
      if (varrecdim!=(-1)) outstart[varrecdim]=r;
      if (varrecdim==(-1) && r > 0) continue;
      if (verbose > 1)
        printf("      writing to disk\n");
      if (ncvarput(ncoutfile->ncfid,v,outstart,count,varbuf[v])==(-1))
        {
         fprintf(stderr,"Error: cannot write variable \"%d\"'s values!\n",
                 v); return(1);
        }
     }
   return(0);
  }


/*
  U.S. Department of Commerce (DOC) Software License for "mppnccombine"
  written at NOAA's Geophysical Fluid Dynamics Laboratory, Princeton
  Forrestal Campus

  1. Scope of License

  Subject to all the terms and conditions of this license, DOC grants USER the
  royalty-free, nonexclusive, nontransferable, and worldwide rights to
  reproduce, modify, and distribute "mppnccombine", herein referred to as the
  PRODUCT.

  2. Conditions and Limitations of Use

  Warranties.  Neither the U.S. Government, nor any agency or employee
  thereof, makes any warranties, expressed or implied, with respect to the
  PRODUCT provided under this license, including but not limited to the
  implied warranties or merchantability and fitness for any particular
  purpose.

  Liability.  In no event shall the U.S. Government, nor any agency or
  employee thereof, be liable for any direct, indirect, or consequential
  damages flowing from the use of the PRODUCT provided under this license.

  Non-Assignment.  Neither this license nor any rights granted hereunder are
  transferable or assignable without the explicit prior written consent of
  DOC.

  Names and Logos.  USER shall not substitute its name or logo for the name or
  logo of DOC, or any of its agencies, in identification of the PRODUCT.

  Export of Technology.  USER shall comply with all U.S. laws and regulations
  restricting the export of the PRODUCT to other countries.

  Governing Law.  This license shall be governed by the laws of United States
  as interpreted and applied  by the Federal courts in the District of
  Columbia.

  3. Term of License

  This license shall remain in effect as long as USER uses the PRODUCT in
  accordance with Paragraphs 1 and 2.
*/
