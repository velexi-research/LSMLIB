/*
 * File:        lsm_curvature2d_local.h
*  Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
*  Revision:    $Revision: 1.19 $
*  Modified:    $Date$
*  Description: F77 routines for computing 2D curvature
*/

#ifndef INCLUDED_LSM_CURVATURE_2D_LOCAL_H
#define INCLUDED_LSM_CURVATURE_2D_LOCAL_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_curvature2d_local.h
 *
 * \brief
 * @ref lsm_curvature2d_local.h provides various ways to compute curvature
 * in 2d (differences mostly in type of stencil used).
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define LSM2D_COMPUTE_MEAN_CURVATURE_ORDER2_LOCAL   \
                                   lsm2dcomputemeancurvatureorder2local_
#define LSM2D_COMPUTE_MEAN_CURVATURE_9STENCIL_LOCAL \
                                   lsm2dcomputemeancurvature9stencillocal_
#define LSM2D_COMPUTE_MEAN_CURVATURE_SGN_DIST_LOCAL \
                                   lsm2dcomputemeancurvaturesgndistlocal_
#define LSM2D_COMPUTE_MEAN_CURVATURE_ORDER4_LOCAL   \
                                   lsm2dcomputemeancurvatureorder4local_	
/*!
*
*  LSM2D_COMPUTE_MEAN_CURVATURE_ORDER2_LOCAL() computes mean curvature
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
*    dx, dy      (in):   grid spacing
*    index_[xy]  (in):  [xy] coordinates of local (narrow band) points
*    n*_index    (in):  index range of points to loop over in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*
*/  
void  LSM2D_COMPUTE_MEAN_CURVATURE_ORDER2_LOCAL(
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
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb);

/*!
*
*  LSM2D_COMPUTE_MEAN_CURVATURE_9STENCIL_LOCAL() computes  mean
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
*    index_[xy](in):    [xy] coordinates of local (narrow band) points
*    n*_index(in):      index range of points to loop over in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):       upper limit narrow band value for voxels in 
*                       fillbox
*
*/
void  LSM2D_COMPUTE_MEAN_CURVATURE_9STENCIL_LOCAL(
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
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb);    

/*!
*
*  LSM2D_COMPUTE_MEAN_CURVATURE_SGN_DIST_LOCAL() computes 
*  kappa = (laplacian (phi) ) / ( 1 - phi * laplacian (phi))
*  For a signed distance function  phi, kappa equals mean curvature 
*  on phi = 0 and is constant in direction normal to the zero level set.
*  See Zhao/Chan/Merriman/Osher 1996.
*
*  Arguments:
*    kappa (in/out):     curvature array
*    phi  (in):          level set function (assumed signed distance func)
*    *_gb (in):          index range for ghostbox
*    index_[xy](in):    [xy] coordinates of local (narrow band) points
*    n*_index(in):      index range of points to loop over in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):       upper limit narrow band value for voxels in 
*                       fillbox
*  Notes: 
*   - Be careful in regions where denominator above is zero, computed value
*     will default to the value of Laplacian
*
*/
void    LSM2D_COMPUTE_MEAN_CURVATURE_SGN_DIST_LOCAL(
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
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb);
 
/*!
*
*  LSM2D_COMPUTE_MEAN_CURVATURE_ORDER4_LOCAL() computes mean curvature
*  kappa = ( phi_xx*phi_y^2 + phi_yy*phi_x^2 - 2*phi_xy*phi_x*phi_y  )/
*          ( | grad phi | ^ 3 )
*  using 4th order central differncing for second order derivatives 
*  (phi_x,phi_y should be precomputed with 4th order too)
*
*  Arguments:
*    kappa     (in/out): curvature data array
*    phi          (in):  level set function
*    phi_*        (in):  first order derivatives of phi
*    grad_phi_mag (in):  gradient magnitude of phi
*    *_gb        (in):   index range for ghostbox
*    index_[xy](in):  [xy] coordinates of local (narrow band) points
*    n*_index(in):    index range of points to loop over in index_*
*    narrow_band(in): array that marks voxels outside desired fillbox
*    mark_fb(in):     upper limit narrow band value for voxels in 
*                     fillbox
*
*/
void LSM2D_COMPUTE_MEAN_CURVATURE_ORDER4_LOCAL(
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
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb);
       
#ifdef __cplusplus
}
#endif
  
#endif 
