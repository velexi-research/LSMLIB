/*
 * File:        lsm_field_extension1d.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.7 $
 * Modified:    $Date$
 * Description: Header file for 1D Fortran 77 field extension subroutines
 */

#ifndef INCLUDED_LSM_FIELD_EXTENSION_1D_H
#define INCLUDED_LSM_FIELD_EXTENSION_1D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_field_extension1d.h
 *
 * \brief
 * @ref lsm_field_extension1d.h provides support for computing the right-hand
 * side of the field extension equation in one space dimension.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in             name in
 *      C/C++ code          Fortran code
 *      ----------          ------------
 */
#define LSM1D_COMPUTE_FIELD_EXTENSION_EQN_RHS                     \
                            lsm1dcomputefieldextensioneqnrhs_


/*!
 * LSM1D_COMPUTE_FIELD_EXTENSION_EQN_RHS() computes right-hand side 
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
 *  - dx (in):               grid spacing
 *  - *_gb (in):             index range for ghostbox
 *  - *_fb (in):             index range for fillbox
 *
 * Return value:             none
 *
 * NOTES:
 * - phi requires at least one ghost cell. 
 *
 */
void LSM1D_COMPUTE_FIELD_EXTENSION_EQN_RHS(
  LSMLIB_REAL *rhs,
  const int *ilo_rhs_gb, 
  const int *ihi_rhs_gb,
  const LSMLIB_REAL *S,
  const int *ilo_S_gb, 
  const int *ihi_S_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const LSMLIB_REAL *S_x_upwind, 
  const int *ilo_grad_S_upwind_gb, 
  const int *ihi_grad_S_upwind_gb,
  const LSMLIB_REAL *signed_normal_x,
  const int *ilo_signed_normal_gb, 
  const int *ihi_signed_normal_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const LSMLIB_REAL *dx);

#ifdef __cplusplus
}
#endif

#endif
