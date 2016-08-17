/*
 * File:        lsm_reinitialization3d_local.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for Fortran 77 3D narrow-band 
 *              reinitialization routines
 */

#ifndef INCLUDED_LSM_3D_REINITIALIZATION_3D_H_LOCAL
#define INCLUDED_LSM_3D_REINITIALIZATION_3D_H_LOCAL

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_reinitialization3d.h
 *
 * \brief
 * @ref lsm_reinitialization3d.h provides support for computing the right-hand
 * side of the reinitialization and orthogonalization equations in three
 * space dimensions for narrow banding approach (localization).
 *
 */

/* Link between C/C++ and Fortran function names 
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS_LOCAL                              \
                                     lsm3dcomputereinitializationeqnrhslocal_
#define LSM3D_COMPUTE_ORTHOGONALIZATION_EQN_RHS_LOCAL                           \
                                     lsm3dcomputeorthogonalizationeqnrhslocal_				     

void LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS_LOCAL(
  LSMLIB_REAL *reinit_rhs,
  const int *ilo_rhs_gb, 
  const int *ihi_rhs_gb, 
  const int *jlo_rhs_gb, 
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb, 
  const int *khi_rhs_gb,
  const LSMLIB_REAL* phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb, 
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const LSMLIB_REAL* phi0,
  const int *ilo_phi0_gb, 
  const int *ihi_phi0_gb, 
  const int *jlo_phi0_gb, 
  const int *jhi_phi0_gb,
  const int *klo_phi0_gb, 
  const int *khi_phi0_gb,
  const LSMLIB_REAL *phi_x_plus, 
  const LSMLIB_REAL *phi_y_plus,
  const LSMLIB_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,   
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,   
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,   
  const int *khi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus, 
  const LSMLIB_REAL *phi_y_minus,
  const LSMLIB_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,   
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,   
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,   
  const int *khi_grad_phi_minus_gb,
  const LSMLIB_REAL *dx, 
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const int *use_phi0_for_sgn,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index,
  const int *nhi_index,  
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const int *klo_nb_gb,
  const int *khi_nb_gb,
  const unsigned char *mark_fb);
  
void LSM3D_COMPUTE_ORTHOGONALIZATION_EQN_RHS_LOCAL(
  LSMLIB_REAL *ortho_rhs,
  const int *ilo_rhs_gb, 
  const int *ihi_rhs_gb, 
  const int *jlo_rhs_gb, 
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb, 
  const int *khi_rhs_gb,
  const LSMLIB_REAL *phi_x_plus, 
  const LSMLIB_REAL *phi_y_plus,
  const LSMLIB_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,   
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,   
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,   
  const int *khi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus, 
  const LSMLIB_REAL *phi_y_minus,
  const LSMLIB_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,   
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,   
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,   
  const int *khi_grad_phi_minus_gb,
  const LSMLIB_REAL* psi,
  const int *ilo_psi_gb, 
  const int *ihi_psi_gb, 
  const int *jlo_psi_gb, 
  const int *jhi_psi_gb,
  const int *klo_psi_gb, 
  const int *khi_psi_gb,
  const LSMLIB_REAL *psi_x_plus, 
  const LSMLIB_REAL *psi_y_plus,
  const LSMLIB_REAL *psi_z_plus,
  const int *ilo_grad_psi_plus_gb,   
  const int *ihi_grad_psi_plus_gb,
  const int *jlo_grad_psi_plus_gb,   
  const int *jhi_grad_psi_plus_gb,
  const int *klo_grad_psi_plus_gb,   
  const int *khi_grad_psi_plus_gb,
  const LSMLIB_REAL *psi_x_minus, 
  const LSMLIB_REAL *psi_y_minus,
  const LSMLIB_REAL *psi_z_minus,
  const int *ilo_grad_psi_minus_gb,   
  const int *ihi_grad_psi_minus_gb,
  const int *jlo_grad_psi_minus_gb,   
  const int *jhi_grad_psi_minus_gb,
  const int *klo_grad_psi_minus_gb,   
  const int *khi_grad_psi_minus_gb,
  const LSMLIB_REAL *dx, 
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const int *klo_nb_gb,
  const int *khi_nb_gb,
  const unsigned char *mark_fb);

  
#ifdef __cplusplus
}
#endif

#endif
