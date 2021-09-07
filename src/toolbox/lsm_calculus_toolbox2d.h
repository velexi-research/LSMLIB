/*
 * File:        lsm_calculus_toolbox2d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file level set method calculus toolbox functions in 2d.
 */

#ifndef INCLUDED_LSM_CALCULUS_TOOLBOX2D
#define INCLUDED_LSM_CALCULUS_TOOLBOX2D

#include "lsmlib_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_calculus_toolbox2d.h
 *
 * \brief Provides higher order implementation of delta function in 2D.
 */


#define LSM2D_DELTA_FUNCTION_ORDER1 lsm2ddeltafunctionorder1_
#define LSM2D_DELTA_FUNCTION_ORDER2 lsm2ddeltafunctionorder2_

/*!
*
*  LSM2D_DELTA_FUNCTION_ORDER1() computes first order accurate delta
*  function discretization following Smereka, "The numerical approximation of
*  a delta function with application to level set methods", JCP, 2006.
*  The function is supported at a minimal set of points around the zero
*  level set.
*
*  Arguments:
*    phi(in):           level set function
*    norm_phi_* (in):   components of grad(phi)/|grad(phi)|, 
*                       obtained by 2nd order central diff
*    *_gb (in):        index range for ghostbox
*    *_fb (in):        index range for fillbox
*
*/
void LSM2D_DELTA_FUNCTION_ORDER1(
  const LSMLIB_REAL *phi,
  const LSMLIB_REAL *delta,
  const int *ilo_gb,
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb, 
  const LSMLIB_REAL *norm_phi_x,
  const LSMLIB_REAL *norm_phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL    *dx,
  const LSMLIB_REAL    *dy);
  
/*!
*
*  LSM2D_DELTA_FUNCTION_ORDER2() computes second order accurate delta
*  function discretization of a zero level set for a given ls function
*  following P. Smereka, "The numerical approximation of a delta function 
*  with application to level set methods", JCP, 2006.
*  The function is supported at a minimal set of gridpoints near the zero
*  level set.
*
*  Arguments:
*    phi(in):           level set function
*    delta(out):        discretized delta function corresp. to zero level
*    norm_phi_* (in):   components of grad(phi)/|grad(phi)|, 
*                       obtained by 2nd order central diff
*    *_gb (in):        index range for ghostbox
*    *_fb (in):        index range for fillbox
*    dx,dy(in):        grid spacing
*/
void LSM2D_DELTA_FUNCTION_ORDER2(
  const LSMLIB_REAL *phi,
  const LSMLIB_REAL *delta,
  const int *ilo_gb,
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb, 
  const LSMLIB_REAL *norm_phi_x,
  const LSMLIB_REAL *norm_phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL    *dx,
  const LSMLIB_REAL    *dy);
  
#ifdef __cplusplus
}
#endif

#endif
