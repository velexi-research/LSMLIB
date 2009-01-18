/*
 * File:        lsm_curvature2d.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for Fortran 77 2D curvature routines.
 */

#ifndef INCLUDED_LSM_CURVATURE_2D_H
#define INCLUDED_LSM_CURVATURE_2D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_curvature2d.h
 *
 * \brief
 * @ref lsm_curvature2d.h provides various ways to compute curvature
 * in 2d (differences mostly in type of stencil used).
 */


/* Link between C/C++ and Fortran function names 
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */

#define LSM2D_COMPUTE_MEAN_CURVATURE_ORDER2    lsm2dcomputemeancurvatureorder2_
#define LSM2D_COMPUTE_MEAN_CURVATURE_9STENCIL  lsm2dcomputemeancurvature9stencil_
#define LSM2D_COMPUTE_MEAN_CURVATURE_SGN_DIST  lsm2dcomputemeancurvaturesgndist_	
#define LSM2D_COMPUTE_MEAN_CURVATURE_ORDER4    lsm2dcomputemeancurvatureorder4_


  
/*!
*
*  LSM2D_COMPUTE_MEAN_CURVATURE_ORDER2() computes mean curvature
*  kappa = ( phi_xx*phi_x^2 + phi_yy*phi_y^2 - 2*phi_xy*phi_x*phi_y )/
*          ( | grad phi | ^ 3 )
*  Standard centered 9 point stencil, second order differencing used.
*  First order derivatives assumed precomputed.
*
*  Arguments:
*    kappa     (in/out): curvature data array
*    phi          (in):  level set function
*    phi_*        (in):  first order derivatives of phi
*    grad_phi_mag (in):  gradient magnitude of phi
*    *_gb        (in):   index range for ghostbox
*    *_fb   (in):        index range for fillbox
*
*/
void    LSM2D_COMPUTE_MEAN_CURVATURE_ORDER2(
  LSMLIB_REAL *kappa,
  const int *ilo_kappa_gb, 
  const int *ihi_kappa_gb,
  const int *jlo_kappa_gb, 
  const int *jhi_kappa_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,  
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb,  
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);
  
/*!
*
*  LSM2D_COMPUTE_MEAN_CURVATURE_9STENCIL() computes mean
*  kappa = div ( grad*(phi) / |grad(phi)| ) 
*  This is mean curvature for any level set w.r.t. to side phi < level, 
*  or negative curvature otherwise.
*  The computation follows 9 point stencil of the nearest neighbors 
*  (formulae from Zhao/Chan/Merriman/Osher 1996.)
*
*  Arguments:
*    kappa (in/out):     curvature array
*    phi     (in):       level set function array
*    dx,dy   (in):       spacing
*    *_gb    (in):       index range for ghostbox
*    *_fb   (in):        index range for fillbox
*
*/

void    LSM2D_COMPUTE_MEAN_CURVATURE_9STENCIL(
  LSMLIB_REAL *kappa,
  const int *ilo_kappa_gb, 
  const int *ihi_kappa_gb,
  const int *jlo_kappa_gb, 
  const int *jhi_kappa_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb,  
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);
  
  
/***********************************************************************
*
*  LSM2D_COMPUTE_MEAN_CURVATURE_SGN_DIST() computes
*  kappa = (laplacian phi) / ( 1 - phi * laplacian phi)
*
*  For a signed distance function  phi, kappa equals mean curvature 
*  on phi = 0 and is constant in direction normal to the zero level set.
*  See Zhao/Chan/Merriman/Osher 1996.
*
*  Arguments:
*    kappa (in/out):    curvature array
*    phi  (in):         level set function (assumed signed distance func)
*    *_gb (in):         index range for ghostbox
*    *_fb (in):         index range for fillbox
*  
*  Notes: 
*   - Be careful in regions where denominator above is zero, computed value
*     will default to the value of Laplacian
*/  
void  LSM2D_COMPUTE_MEAN_CURVATURE_SGN_DIST(
  LSMLIB_REAL *kappa,
  const int *ilo_kappa_gb, 
  const int *ihi_kappa_gb,
  const int *jlo_kappa_gb, 
  const int *jhi_kappa_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb,  
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);  

/*!
*
*  LSM2D_COMPUTE_MEAN_CURVATURE_ORDER4() computes mean curvature
*  kappa = ( phi_xx*phi_y^2 + phi_yy*phi_x^2 - 2*phi_xy*phi_x*phi_y  )/
*          ( | grad phi | ^ 3 )
*  using 4th order central differencing for second order derivatives 
*  (phi_x,phi_y assumed precomputed with 4th order)
*
*  Arguments:
*    kappa     (in/out): curvature data array
*    phi          (in):  level set function
*    phi_*        (in):  first order derivatives of phi
*    grad_phi_mag (in):  gradient magnitude of phi
*    *_gb        (in):   index range for ghostbox
*    *_fb        (in):   index range for fillbox
*
*/
void LSM2D_COMPUTE_MEAN_CURVATURE_ORDER4(
  LSMLIB_REAL *kappa,
  const int *ilo_kappa_gb, 
  const int *ihi_kappa_gb,
  const int *jlo_kappa_gb, 
  const int *jhi_kappa_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
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
