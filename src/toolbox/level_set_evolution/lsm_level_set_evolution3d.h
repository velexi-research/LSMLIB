/*
 * File:        lsm_level_set_evolution3d.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.8 $
 * Modified:    $Date: 2006/10/28 05:09:39 $
 * Description: Header file for 3D Fortran 77 level set evolution equation
 *              subroutines
 */

#ifndef INCLUDED_LSM_LEVEL_SET_EVOLUTION_3D_H
#define INCLUDED_LSM_LEVEL_SET_EVOLUTION_3D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_level_set_evolution3d.h
 *
 * \brief
 * @ref lsm_level_set_evolution3d.h provides support for contributing to the
 * right-hand side of the level set evolution equation in three space
 * dimensions.
 * 
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                                name in
 *      C/C++ code                             Fortran code
 *      ----------                             ------------
 */
#define LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS       lsm3dzerooutlevelseteqnrhs_
#define LSM3D_ADD_ADVECTION_TERM_TO_LSE_RHS    lsm3daddadvectiontermtolserhs_
#define LSM3D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS   lsm3daddnormalveltermtolserhs_
#define LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS   \
                                          lsm3daddconstnormalveltermtolserhs_
#define LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS    lsm3daddconstcurvtermtolserhs_				  
					  
					  

/*!
 * LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS() zeros out the right-hand side of 
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
void LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *klo_lse_rhs_gb, 
  const int *khi_lse_rhs_gb);


/*!
 * LSM3D_ADD_ADVECTION_TERM_TO_LSE_RHS() adds the contribution of an 
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
void LSM3D_ADD_ADVECTION_TERM_TO_LSE_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *klo_lse_rhs_gb, 
  const int *khi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x, 
  const LSMLIB_REAL *phi_y, 
  const LSMLIB_REAL *phi_z,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb, 
  const int *khi_grad_phi_gb,
  const LSMLIB_REAL *vel_x, 
  const LSMLIB_REAL *vel_y, 
  const LSMLIB_REAL *vel_z,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *klo_vel_gb, 
  const int *khi_vel_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb, 
  const int *jhi_fb,
  const int *klo_fb, 
  const int *khi_fb);


/*!
 * LSM3D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS() adds the contribution of a 
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
 * Return value:         none
 */
void LSM3D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *klo_lse_rhs_gb, 
  const int *khi_lse_rhs_gb,
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
  const LSMLIB_REAL *vel_n, 
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *klo_vel_gb, 
  const int *khi_vel_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb, 
  const int *jhi_fb,
  const int *klo_fb, 
  const int *khi_fb);

					
/*!
 * LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS() adds the contribution 
 * of a normal (scalar) velocity term to the right-hand side of the 
 * level set equation when it is written in the form:
 *   
 * \f[
 *   
 *    \phi_t = -V_n |\nabla \phi| + ... 
 *   
 * \f]
 *  
 * where the normal velocity, \f$ V_n \f$ is a constant.
 *
 * Arguments:
 *  - lse_rhs (in/out):  right-hand of level set equation
 *  - phi_*_plus (in):   components of forward approx to \f$ \nabla \phi \f$
 *                       at t = t_cur
 *  - phi_*_minus (in):  components of backward approx to \f$ \nabla \phi \f$
 *                       at t = t_cur
 *  - vel_n (in):        constant normal velocity at t = t_cur
 *  - *_gb (in):         index range for ghostbox
 *  - *_fb (in):         index range for fillbox
 *
 * Return value:         none
 */
void LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *klo_lse_rhs_gb, 
  const int *khi_lse_rhs_gb,
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
  const LSMLIB_REAL *vel_n, 
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb, 
  const int *jhi_fb,
  const int *klo_fb, 
  const int *khi_fb);


/*!
 * LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS() adds the contribution 
 * of a normal (scalar) velocity term to the right-hand side of the 
 * level set equation when it is written in the form:
 *   
 * \f[
 *   
 *    \phi_t = -b kappa |\nabla \phi| + ... 
 *   
 * \f]
 *  
 * where the \f$ kappa \f$ is the mean curvature and \f$ b \f$ is a 
 * constant.
 *
 * Arguments:
 *  - lse_rhs (in/out):  right-hand of level set equation
 *  - phi_* (in):        first- and second-order partial derivatives 
 *                       of \f$ \phi \f$ 
 *  - b (in):            proporationality constant relating curvature
 *                       to the normal velocity
 *  - *_gb (in):         index range for ghostbox
 *  - *_fb (in):         index range for fillbox
 *
 * Return value:         none
 */
void   LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS(
  const LSMLIB_REAL  *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *klo_lse_rhs_gb, 
  const int *khi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x, 
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const  LSMLIB_REAL *phi_xx,
  const  LSMLIB_REAL *phi_xy,
  const  LSMLIB_REAL *phi_xz,
  const  LSMLIB_REAL *phi_yy,
  const  LSMLIB_REAL *phi_yz,
  const  LSMLIB_REAL *phi_zz,
  const  int *ilo_grad2_phi_gb,
  const  int *ihi_grad2_phi_gb,
  const  int *jlo_grad2_phi_gb, 
  const  int *jhi_grad2_phi_gb,
  const  int *klo_grad2_phi_gb, 
  const  int *khi_grad2_phi_gb,
  const LSMLIB_REAL *b,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb);  

#ifdef __cplusplus
}
#endif

#endif
