/*
 * File:        lsm_calculus_toolbox3d.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and  Kevin T. Chu
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file level set method calculus toolbox functions
 */

#ifndef INCLUDED_LSM_CALCULUS_TOOLBOX3D
#define INCLUDED_LSM_CALCULUS_TOOLBOX3D

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_calculus_toolbox3d.h
 *
 * \brief Provides higher order implementation of delta function in 3D.
 */


#define LSM3D_DELTA_FUNCTION_ORDER1 lsm3ddeltafunctionorder1_

  
/*!
*
*  LSM3D_DELTA_FUNCTION_ORDER1() computes first order accurate delta
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
*  Note: still being tested.
*
*/
void LSM3D_DELTA_FUNCTION_ORDER1(
  const LSMLIB_REAL *phi,
  const LSMLIB_REAL *delta,
  const int *ilo_gb,
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb, 
  const int *klo_gb,
  const int *khi_gb, 
  const LSMLIB_REAL *norm_phi_x,
  const LSMLIB_REAL *norm_phi_y,
  const LSMLIB_REAL *norm_phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const LSMLIB_REAL    *dx,
  const LSMLIB_REAL    *dy,
  const LSMLIB_REAL    *dz);
  
#ifdef __cplusplus
}
#endif

#endif
