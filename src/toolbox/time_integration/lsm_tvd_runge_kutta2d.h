/*
 * File:        lsm_tvd_runge_kutta2d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for Fortran 77 2D TVD Runge-Kutta routines
 */

#ifndef INCLUDED_LSM_TVD_RUNGE_KUTTA_2D_H
#define INCLUDED_LSM_TVD_RUNGE_KUTTA_2D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_tvd_runge_kutta2d.h
 *
 * \brief
 * @ref lsm_tvd_runge_kutta2d.h provides support for time integration of
 * partial differential equations in two space dimensions via 
 * total-variation diminishing Runge-Kutta methods.  Support is provided 
 * for first-, second-, and third-order time integration.
 * 
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                             name in
 *      C/C++ code                          Fortran code
 *      ----------                          ------------
 */
#define LSM2D_RK1_STEP                      lsm2drk1step_
#define LSM2D_TVD_RK2_STAGE1                lsm2dtvdrk2stage1_
#define LSM2D_TVD_RK2_STAGE2                lsm2dtvdrk2stage2_
#define LSM2D_TVD_RK3_STAGE1                lsm2dtvdrk3stage1_
#define LSM2D_TVD_RK3_STAGE2                lsm2dtvdrk3stage2_
#define LSM2D_TVD_RK3_STAGE3                lsm2dtvdrk3stage3_


/*!
 * LSM2D_RK1_STEP() takes a single first-order Runge-Kutta (i.e. Forward
 * Euler) step.
 *
 * Arguments:
 *  - u_next (out):  u(t_cur+dt)
 *  - u_cur (in):    u(t_cur)
 *  - rhs (in):      right-hand side of time evolution equation
 *  - dt (in):       step size
 *  - *_gb (in):     index range for ghostbox
 *  - *_fb (in):     index range for fillbox
 *
 * Return value:     none
 */
void LSM2D_RK1_STEP(
  LSMLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dt);


/*!
 * LSM2D_TVD_RK2_STAGE1() advances the solution through first stage of the
 * second-order TVD Runge-Kutta method.
 * 
 * Arguments:
 *  - u_stage1 (out):  u_approx(t_cur+dt)
 *  - u_cur (in):      u(t_cur)
 *  - rhs (in):        right-hand side of time evolution equation
 *  - dt (in):         step size
 *  - *_gb (in):       index range for ghostbox
 *  - *_fb (in):       index range for fillbox
 *
 * Return value:       none
 *
 * NOTES:
 *  - the first stage of TVD RK2 is identical to a single RK1 step
 *
 */
void LSM2D_TVD_RK2_STAGE1(
  LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dt);


/*!
 * LSM2D_TVD_RK2_STAGE2() completes advancing the solution through a
 * single step of the second-order TVD Runge-Kutta method
 *   
 * Arguments:
 *  - u_next (out):   u(t_cur+dt)
 *  - u_stage1 (in):  u_approx(t_cur+dt)
 *  - u_cur (in):     u(t_cur)
 *  - rhs (in):       right-hand side of time evolution equation
 *  - dt (in):        step size
 *  - *_gb (in):      index range for ghostbox
 *  - *_fb (in):      index range for fillbox
 *
 * Return value:      none
 */
void LSM2D_TVD_RK2_STAGE2(
  LSMLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dt);


/*!
 * LSM2D_TVD_RK3_STAGE1() advances the solution through first stage of the
 * third-order TVD Runge-Kutta method.
 *
 * Arguments: 
 *  - u_stage1 (out):  u_approx(t_cur+dt)
 *  - u_cur (in):      u(t_cur)
 *  - rhs (in):        right-hand side of time evolution equation
 *  - dt (in):         step size
 *  - *_gb (in):       index range for ghostbox
 *  - *_fb (in):       index range for fillbox
 * 
 * Return value:       none 
 * 
 * NOTES: 
 *  - the first stage of TVD RK3 is identical to a single RK1 step
 *
 */
void LSM2D_TVD_RK3_STAGE1(
  LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dt);


/*!
 * LSM2D_TVD_RK3_STAGE2() advances the solution through second stage of
 * the third-order TVD Runge-Kutta method.
 *
 * Arguments:
 *  - u_stage2 (out):  u_approx(t_cur+dt/2)
 *  - u_stage1 (in):   u_approx(t_cur+dt)
 *  - u_cur (in):      u(t_cur)
 *  - rhs (in):        right-hand side of time evolution equation
 *  - dt (in):         step size
 *  - *_gb (in):       index range for ghostbox
 *  - *_fb (in):       index range for fillbox
 *
 * Return value:       none
 */
void LSM2D_TVD_RK3_STAGE2(
  LSMLIB_REAL *u_stage2,
  const int *ilo_u_stage2_gb,
  const int *ihi_u_stage2_gb,
  const int *jlo_u_stage2_gb,
  const int *jhi_u_stage2_gb,
  const LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dt);

 
/*!
 * LSM2D_TVD_RK3_STAGE3() completes advancing the solution through a
 * single step of the third-order TVD Runge-Kutta method.
 *   
 * Arguments:
 *  - u_next (out):   u(t_cur+dt)
 *  - u_stage2 (in):  u_approx(t_cur+dt/2)
 *  - u_cur (in):     u(t_cur)
 *  - rhs (in):       right-hand side of time evolution equation
 *  - dt (in):        step size
 *  - *_gb (in):      index range for ghostbox
 *  - *_fb (in):      index range for fillbox
 *
 * Return value:      none
 */
void LSM2D_TVD_RK3_STAGE3(
  LSMLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const LSMLIB_REAL *u_stage2,
  const int *ilo_u_stage2_gb,
  const int *ihi_u_stage2_gb,
  const int *jlo_u_stage2_gb,
  const int *jhi_u_stage2_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dt);

#ifdef __cplusplus
}
#endif

#endif
