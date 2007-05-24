/*
 * File:        lsm_field_extension3d.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.7 $
 * Modified:    $Date: 2006/05/18 23:20:10 $
 * Description: Header file for 3D Fortran 77 field extension subroutines
 */

#ifndef INCLUDED_LSM_FIELD_EXTENSION_3D_H
#define INCLUDED_LSM_FIELD_EXTENSION_3D_H

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_field_extension3d.h
 *
 * \brief
 * @ref lsm_field_extension3d.h provides support for computing the right-hand
 * side of the field extension equation in three space dimensions.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in             name in
 *      C/C++ code          Fortran code
 *      ----------          ------------
 */
#define LSM3D_COMPUTE_FIELD_EXTENSION_EQN_RHS                     \
                            lsm3dcomputefieldextensioneqnrhs_


/*!
 * LSM3D_COMPUTE_FIELD_EXTENSION_EQN_RHS() computes right-hand side 
 * of the field extension equation when it is written in the form:
 *
 * \f[
 *
 *    S_t = -sgn(\phi) \vec{N} \cdot \nabla S
 *
 * \f]
 *
 * Arguments:
 *  - rhs (out):             right-hand side of field extension equation
 *  - S (in):                field to be extended off of the zero level set
 *  - phi (in):              level set function used to compute normal vector
 *  - S_*_upwind (in):       upwind spatial derivatives for \f$ \nabla S \f$
 *  - signed_normal_* (in):  signed normal
 *  - dx, dy, dz (in):       grid spacing
 *  - *_gb (in):             index range for ghostbox
 *  - *_fb (in):             index range for fillbox
 *
 * Return value:             none
 *
 * NOTES:
 * - phi requires at least one ghost cell.
 *
 */
void LSM3D_COMPUTE_FIELD_EXTENSION_EQN_RHS(
  double *rhs,
  const int *ilo_rhs_gb, 
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb, 
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb, 
  const int *khi_rhs_gb,
  const double *S,
  const int *ilo_S_gb, 
  const int *ihi_S_gb,
  const int *jlo_S_gb, 
  const int *jhi_S_gb,
  const int *klo_S_gb, 
  const int *khi_S_gb,
  const double *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const double *S_x_upwind, 
  const double *S_y_upwind, 
  const double *S_z_upwind, 
  const int *ilo_grad_S_upwind_gb, 
  const int *ihi_grad_S_upwind_gb,
  const int *jlo_grad_S_upwind_gb, 
  const int *jhi_grad_S_upwind_gb,
  const int *klo_grad_S_upwind_gb, 
  const int *khi_grad_S_upwind_gb,
  const double *signed_normal_x,
  const double *signed_normal_y,
  const double *signed_normal_z,
  const int *ilo_signed_normal_gb, 
  const int *ihi_signed_normal_gb,
  const int *jlo_signed_normal_gb, 
  const int *jhi_signed_normal_gb,
  const int *klo_signed_normal_gb, 
  const int *khi_signed_normal_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb, 
  const int *jhi_fb,
  const int *klo_fb, 
  const int *khi_fb,
  const double *dx,
  const double *dy,
  const double *dz);

#ifdef __cplusplus
}
#endif

#endif
