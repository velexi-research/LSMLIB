/*
 * File:        lsm_curvature3d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for Fortran 77 3D curvature routines.
 */
 
/*! \file lsm_curvature3d.h
 *
 * \brief
 * @ref lsm_curvature3d.h provides various ways to compute curvature
 * in 2d (differences mostly in type of stencil used).
 */
  
#ifndef INCLUDED_LSM_CURVATURE_3D_H
#define INCLUDED_LSM_CURVATURE_3D_H

#include "lsmlib_config.h"

#ifdef __cplusplus
extern "C" {
#endif




/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */


#define LSM3D_COMPUTE_MEAN_CURVATURE_ORDER2   \
                                     lsm3dcomputemeancurvatureorder2_
#define LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER2  \
                                   lsm3dcomputegaussiancurvatureorder2_		     	
#define LSM3D_COMPUTE_MEAN_CURVATURE_ORDER4   \
                                     lsm3dcomputemeancurvatureorder4_
#define LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER4  \
                                   lsm3dcomputegaussiancurvatureorder4_	
				   

/*!
*  LSM3D_COMPUTE_MEAN_CURVATURE_ORDER2) computes mean curvature
*  kappa = ( phi_xx*phi_y^2 + phi_yy*phi_x^2 - 2*phi_xy*phi_x*phi_y + 
*            phi_xx*phi_z^2 + phi_zz*phi_x^2 - 2*phi_xz*phi_x*phi_z +
*            phi_yy*phi_z^2 + phi_zz*phi_y^2 - 2*phi_yz*phi_y*phi_z )
*          ( | grad phi | ^ 3 )
*  Standard centered 27 point stencil, second order differencing used.
*  First order derivatives assumed precomputed.
*
*  Arguments:
*    kappa       (out):  curvature data array
*    phi          (in):  level set function
*    phi_*        (in):  first order derivatives of phi
*    grad_phi_mag (in):  gradient magnitude of phi
*    *_gb         (in):   index range for ghostbox
*    dx, dy, dz  (in):   grid spacing
*
*  Notes:
*   - memory for 'kappa' array assumed preallocated
*/
void LSM3D_COMPUTE_MEAN_CURVATURE_ORDER2(
  LSMLIB_REAL *kappa,
  const int *ilo_kappa_gb,
  const int *ihi_kappa_gb,
  const int *jlo_kappa_gb,
  const int *jhi_kappa_gb,
  const int *klo_kappa_gb,
  const int *khi_kappa_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *phi_z,
  const LSMLIB_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const int *ilo_kappa_fb,
  const int *ihi_kappa_fb,
  const int *jlo_kappa_fb,
  const int *jhi_kappa_fb,
  const int *klo_kappa_fb,
  const int *khi_kappa_fb,  
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz); 

/*!
*
*  lsm3dComputeGaussianCurvatureOrder2() computes Gaussian curvature
*  kappa = [  phi_x^2*(phi_yy*phi_zz - phi_yz^2) +
*             phi_y^2*(phi_xx*phi_zz - phi_xz^2) + 
*             phi_z^2*(phi_xx*phi_yy - phi_xy^2) + 
*         2*( phi_x*phi_y*(phi_xz*phi_yz - phi_xy*phi_zz) +
*             phi_y*phi_z*(phi_xy*phi_xz - phi_yz*phi_xx) + 
*             phi_x*phi_z*(phi_xy*phi_yz - phi_xz*phi_yy) ) ] / 
*          ( | grad phi | ^ 4 )
*  Standard centered 27 point stencil, second order differencing used.
*  First order derivatives assumed precomputed.
*
*  Arguments:
*    kappa     (in/out): curvature data array
*    phi          (in):  level set function
*    phi_*        (in):  first order derivatives of phi
*    grad_phi_mag (in):  gradient magnitude of phi
*    *_gb        (in):   index range for ghostbox
*    *_fb        (in):   index range for fillbox
*    dx,dy,dz    (in):   grid spacing
*
*   NOTES: Data array 'kappa' assumed pre-allocated
*
*/
void LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER2(
  LSMLIB_REAL *kappa,
  const int *ilo_kappa_gb,
  const int *ihi_kappa_gb,
  const int *jlo_kappa_gb,
  const int *jhi_kappa_gb,
  const int *klo_kappa_gb,
  const int *khi_kappa_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *phi_z,
  const LSMLIB_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,  
  const int *ilo_kappa_fb,
  const int *ihi_kappa_fb,
  const int *jlo_kappa_fb,
  const int *jhi_kappa_fb,
  const int *klo_kappa_fb,
  const int *khi_kappa_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz); 

/*!
*  LSM3D_COMPUTE_MEAN_CURVATURE_ORDER4() computes mean curvature
*  kappa = ( phi_xx*phi_y^2 + phi_yy*phi_x^2 - 2*phi_xy*phi_x*phi_y + 
*            phi_xx*phi_z^2 + phi_zz*phi_x^2 - 2*phi_xz*phi_x*phi_z +
*            phi_yy*phi_z^2 + phi_zz*phi_y^2 - 2*phi_yz*phi_y*phi_z )
*          ( | grad phi | ^ 3 )
*  using 4th order central differencing.
*  First order derivatives assumed precomputed (with 4th order).
*
*  Arguments:
*    kappa       (out):  curvature data array
*    phi          (in):  level set function
*    phi_*        (in):  first order derivatives of phi
*    grad_phi_mag (in):  gradient magnitude of phi
*    *_gb         (in):   index range for ghostbox
*    dx, dy, dz  (in):   grid spacing
*
*  Notes:
*   - memory for 'kappa' array assumed preallocated
*/
void LSM3D_COMPUTE_MEAN_CURVATURE_ORDER4(
  LSMLIB_REAL *kappa,
  const int *ilo_kappa_gb,
  const int *ihi_kappa_gb,
  const int *jlo_kappa_gb,
  const int *jhi_kappa_gb,
  const int *klo_kappa_gb,
  const int *khi_kappa_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *phi_z,
  const LSMLIB_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const int *ilo_kappa_fb,
  const int *ihi_kappa_fb,
  const int *jlo_kappa_fb,
  const int *jhi_kappa_fb,
  const int *klo_kappa_fb,
  const int *khi_kappa_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz); 


/*!
*
*  lsm3dComputeGaussianCurvatureOrder4() computes Gaussian curvature
*  kappa = [  phi_x^2*(phi_yy*phi_zz - phi_yz^2) +
*             phi_y^2*(phi_xx*phi_zz - phi_xz^2) + 
*             phi_z^2*(phi_xx*phi_yy - phi_xy^2) + 
*         2*( phi_x*phi_y*(phi_xz*phi_yz - phi_xy*phi_zz) +
*             phi_y*phi_z*(phi_xy*phi_xz - phi_yz*phi_xx) + 
*             phi_x*phi_z*(phi_xy*phi_yz - phi_xz*phi_yy) ) ] / 
*          ( | grad phi | ^ 4 )
*  using 4th order central differencing.
*  First order derivatives assumed precomputed (with 4th order).
*
*  Arguments:
*    kappa     (in/out): curvature data array
*    phi          (in):  level set function
*    phi_*        (in):  first order derivatives of phi
*    grad_phi_mag (in):  gradient magnitude of phi
*    *_gb        (in):   index range for ghostbox
*    *_fb        (in):   index range for fillbox
*    dx,dy,dz    (in):   grid spacing
*
*   NOTES: Data array 'kappa' assumed pre-allocated
*
*/
void LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER4(
  LSMLIB_REAL *kappa,
  const int *ilo_kappa_gb,
  const int *ihi_kappa_gb,
  const int *jlo_kappa_gb,
  const int *jhi_kappa_gb,
  const int *klo_kappa_gb,
  const int *khi_kappa_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *phi_z,
  const LSMLIB_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const int *ilo_kappa_fb,
  const int *ihi_kappa_fb,
  const int *jlo_kappa_fb,
  const int *jhi_kappa_fb,
  const int *klo_kappa_fb,
  const int *khi_kappa_fb,  
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz); 
  
#ifdef __cplusplus
}
#endif

#endif
