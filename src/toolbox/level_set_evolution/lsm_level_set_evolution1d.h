/*
 * File:        lsm_level_set_evolution1d.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.7 $
 * Modified:    $Date: 2006/05/18 23:20:12 $
 * Description: Header file for 1D Fortran 77 level set evolution equation
 *              subroutines
 */

#ifndef INCLUDED_LSM_LEVEL_SET_EVOLUTION_1D_H
#define INCLUDED_LSM_LEVEL_SET_EVOLUTION_1D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_level_set_evolution1d.h
 *
 * \brief
 * @ref lsm_level_set_evolution1d.h provides support for contributing to the
 * right-hand side of the level set evolution equation in one space
 * dimension.
 * 
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                                name in
 *      C/C++ code                             Fortran code
 *      ----------                             ------------
 */
#define LSM1D_ZERO_OUT_LEVEL_SET_EQN_RHS       lsm1dzerooutlevelseteqnrhs_
#define LSM1D_ADD_ADVECTION_TERM_TO_LSE_RHS    lsm1daddadvectiontermtolserhs_
#define LSM1D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS   lsm1daddnormalveltermtolserhs_

/*!
 * LSM1D_ZERO_OUT_LEVEL_SET_EQN_RHS() zeros out the right-hand side of 
 * the level set equation when it is written in the form:
 *
 * \f[
 *
 *  \phi_t = ...
 *
 * \f]
 *
 * Arguments:
 *  - lse_rhs (in/out):  right-hand of level set equation
 *  - *_gb (in):         index range for ghostbox
 *
 * Return value:         none
 */
void LSM1D_ZERO_OUT_LEVEL_SET_EQN_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb);


/*!
 * LSM1D_ADD_ADVECTION_TERM_TO_LSE_RHS() adds the contribution of an 
 * advection term (external vector velocity field) to the right-hand 
 * side of the level set equation when it is written in the form:
 *
 * \f[
 *
 *    \phi_t = -\vec{V} \cdot \nabla \phi + ...
 *
 * \f]
 *
 * Arguments:
 *  - lse_rhs (in/out):  right-hand of level set equation
 *  - phi_* (in):        components of \f$ \nabla \phi \f$ at t = t_cur
 *  - vel_* (in):        components of velocity at t = t_cur
 *  - *_gb (in):         index range for ghostbox
 *  - *_fb (in):         index range for fillbox
 *
 * Return value:         none
 */
void LSM1D_ADD_ADVECTION_TERM_TO_LSE_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const LSMLIB_REAL *vel_x,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *ilo_fb, 
  const int *ihi_fb);


/*!
 * LSM1D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS() adds the contribution of a 
 * normal (scalar) velocity term to the right-hand side of the level 
 * set equation when it is written in the form:
 *   
 * \f[
 *   
 *    \phi_t = -V_n |\nabla \phi| + ... 
 *   
 * \f]
 *   
 * Arguments:
 *  - lse_rhs (in/out):  right-hand of level set equation
 *  - phi_*_plus (in):   components of forward approx to \f$ \nabla \phi \f$ 
 *                       at t = t_cur
 *  - phi_*_minus (in):  components of backward approx to \f$ \nabla \phi \f$
 *                       at t = t_cur
 *  - vel_n (in):        normal velocity at t = t_cur
 *  - *_gb (in):         index range for ghostbox
 *  - *_fb (in):         index range for fillbox
 *
 * Return value:        none
 */
void LSM1D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const LSMLIB_REAL *vel_n,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *ilo_fb, 
  const int *ihi_fb);

#ifdef __cplusplus
}
#endif

#endif
