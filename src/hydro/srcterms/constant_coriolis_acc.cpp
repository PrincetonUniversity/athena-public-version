//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief source terms due to constant coriolis acceleration

// Athena++ headers
#include "hydro_srcterms.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::ConstantCoriolisAcceleration
//  \brief Adds source terms for constant coriolis acceleration to conserved variables

void HydroSourceTerms::ConstantCoriolisAcceleration(const Real dt,const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  // acceleration in 1-direction
  if (omega1_ != 0.0 || omega2_ != 0.0 || omega3_ != 0.0) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
#pragma omp parallel for schedule(static)
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        cons(IM1,k,j,i) += 2.*dt*(omega3_*cons(IM2,k,j,i) - omega2_*cons(IM3,k,j,i));
        cons(IM2,k,j,i) += 2.*dt*(omega1_*cons(IM3,k,j,i) - omega3_*cons(IM1,k,j,i));
        cons(IM3,k,j,i) += 2.*dt*(omega2_*cons(IM1,k,j,i) - omega1_*cons(IM2,k,j,i));
      }
    }}
  }

  return;
}
