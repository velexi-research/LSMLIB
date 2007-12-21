/*
 * File:        lsm_field_extension2d.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.7 $
 * Modified:    $Date: 2006/05/18 23:20:10 $
 * Description: Header file for 2D Fortran 77 field extension subroutines
 */

#ifndef INCLUDED_LSM_FIELD_EXTENSION_2D_H
#define INCLUDED_LSM_FIELD_EXTENSION_2D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_field_extension2d.h
 *
 * \brief
 * @ref lsm_field_extension2d.h provides support for computing the right-hand
 * side of the field extension equation in two space dimensions.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in             name in
 *      C/C++ code          Fortran code
 *      ----------          ------------
 */
#define LSM2D_COMPUTE_FIELD_EXTENSION_EQN_RHS                     \
                            lsm2dcomputefieldextensioneqnrhs_


/*!
 * LSM2D_COMPUTE_FIELD_EXTENSION_EQN_RHS() computes right-hand side 
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
 *  - dx, dy (in):           grid spacing
 *  - *_gb (in):             index range for ghostbox
 *  - *_fb (in):             index range for fillbox
 *
 * Return value:             none
 *
 * NOTES:
 * - phi requires at least one ghost cell.
 *
 */
void LSM2D_COMPUTE_FIELD_EXTENSION_EQN_RHS(
  LSMLIB_REAL *rhs,
  const int *ilo_rhs_gb, 
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb, 
  const int *jhi_rhs_gb,
  const LSMLIB_REAL *S,
  const int *ilo_S_gb, 
  const int *ihi_S_gb,
  const int *jlo_S_gb, 
  const int *jhi_S_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const LSMLIB_REAL *S_x_upwind, 
  const LSMLIB_REAL *S_y_upwind, 
  const int *ilo_grad_S_upwind_gb, 
  const int *ihi_grad_S_upwind_gb,
  const int *jlo_grad_S_upwind_gb, 
  const int *jhi_grad_S_upwind_gb,
  const LSMLIB_REAL *signed_normal_x,
  const LSMLIB_REAL *signed_normal_y,
  const int *ilo_signed_normal_gb, 
  const int *ihi_signed_normal_gb,
  const int *jlo_signed_normal_gb, 
  const int *jhi_signed_normal_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb, 
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);

#ifdef __cplusplus
}
#endif

#endif
