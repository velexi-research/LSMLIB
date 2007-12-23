/*
 * File:        lsm_spatial_derivatives2d_local.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.19 $
 * Modified:    $Date: 2006/10/28 01:48:58 $
 * Description: Header file for Fortran 77 2D narrow-band ENO/WENO routines
 */

#ifndef INCLUDED_LSM_SPATIAL_DERIVATIVES_2D_LOCAL_H
#define INCLUDED_LSM_SPATIAL_DERIVATIVES_2D_LOCAL_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_spatial_derivatives2d_local.h
 *
 * \brief
 * @ref lsm_spatial_derivatives2d_local.h provides support for computing spatial
 * derivatives in two space dimensions using high-order ENO and WENO 
 * discretizations for narrow banding (local) approach. 
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define LSM2D_HJ_ENO1_LOCAL              lsm2dhjeno1local_
#define LSM2D_HJ_ENO2_LOCAL              lsm2dhjeno2local_
#define LSM2D_CENTRAL_GRAD_ORDER2_LOCAL  lsm2dcentralgradorder2local_
#define LSM2D_CENTRAL_GRAD_ORDER4_LOCAL  lsm2dcentralgradorder4local_
#define LSM2D_LAPLACIAN_ORDER2_LOCAL     lsm2dlaplacianorder2local_
#define LSM2D_COMPUTE_AVE_GRAD_PHI_LOCAL lsm2dcomputeavegradphilocal_

void LSM2D_HJ_ENO1_LOCAL(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
  const int *nlo_index0,
  const int *nhi_index0,
  const int *nlo_index1,
  const int *nhi_index1,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb,
  const unsigned char *mark_D1);

void LSM2D_HJ_ENO2_LOCAL(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
  const int *nlo_index0,
  const int *nhi_index0,
  const int *nlo_index1,
  const int *nhi_index1,
  const int *nlo_index2,
  const int *nhi_index2,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb,
  const unsigned char *mark_D1,
  const unsigned char *mark_D2);


void LSM2D_CENTRAL_GRAD_ORDER2_LOCAL( 
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb);
  
void LSM2D_CENTRAL_GRAD_ORDER4_LOCAL( 
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb);
 
void LSM2D_LAPLACIAN_ORDER2_LOCAL( 
  LSMLIB_REAL *laplacian_phi,
  const int *ilo_laplacian_phi_gb,
  const int *ihi_laplacian_phi_gb,
  const int *jlo_laplacian_phi_gb,
  const int *jhi_laplacian_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb);
      
void LSM2D_COMPUTE_AVE_GRAD_PHI_LOCAL(
  LSMLIB_REAL *grad_phi_ave,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb);  
 
#ifdef __cplusplus
}
#endif
  
#endif 
