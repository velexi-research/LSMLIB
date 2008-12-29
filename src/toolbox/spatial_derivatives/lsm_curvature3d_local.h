/*
 * File:        lsm_curvature3d_local.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.18 $
 * Modified:    $Date: 2006/07/25 15:00:58 $
 * Description: Header file for Fortran 77 3D curvature routines.
 */
 
/*! \file lsm_curvature3d_local.h
 *
 * \brief
 * @ref lsm_curvature3d_local.h provides various ways to compute curvature
 * in 2d (differences mostly in type of stencil used).
 */

 
#ifndef INCLUDED_LSM_CURVATURE_3D_LOCAL_H
#define INCLUDED_LSM_CURVATURE_3D_LOCAL_H

#include "LSMLIB_config.h"


/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */

#define LSM3D_COMPUTE_MEAN_CURVATURE_ORDER2_LOCAL   \
                                   lsm3dcomputemeancurvatureorder2local_
#define LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER2_LOCAL   \
                                   lsm3dcomputegaussiancurvatureorder2local_
#define LSM3D_COMPUTE_MEAN_CURVATURE_ORDER4_LOCAL   \
                                   lsm3dcomputemeancurvatureorder4local_
#define LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER4_LOCAL   \
                                   lsm3dcomputegaussiancurvatureorder4local_

#ifdef __cplusplus
extern "C" {
#endif


/*!
*
*  LSM3D_COMPUTE_MEAN_CURVATURE_ORDER2_LOCAL() computes mean curvature
*  kappa = div ( grad_phi / |grad_phi|)
*  kappa = ( phi_xx*phi_y^2 + phi_yy*phi_x^2 - 2*phi_xy*phi_x*phi_y + 
*            phi_xx*phi_z^2 + phi_zz*phi_x^2 - 2*phi_xz*phi_x*phi_z +
*            phi_yy*phi_z^2 + phi_zz*phi_y^2 - 2*phi_yz*phi_y*phi_z )/
*          ( | grad phi | ^ 3 )
*  Note that this value is technically twice the mean curvature.
*  Standard centered 27 point stencil, second order differencing used.
*  First order derivatives assumed precomputed.
c
*  Arguments:
*    kappa     (in/out): curvature data array
*    phi          (in):  level set function
*    phi_*        (in):  first order derivatives of phi
*    grad_phi_mag (in):  gradient magnitude of phi
*    *_gb        (in):   index range for ghostbox
*    dx, dy      (in):   grid spacing
*    index_[xyz]  (in):  [xyz] coordinates of local (narrow band) points
*    n*_index    (in):  index range of points to loop over in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*
*/ 
void LSM3D_COMPUTE_MEAN_CURVATURE_ORDER2_LOCAL(
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
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const int *klo_nb_gb,
  const int *khi_nb_gb,
  const unsigned char *mark_fb);
 
/*!
*
*  LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER2_LOCAL() computes Gaussian curvature
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
*    dx, dy      (in):   grid spacing
*    index_[xyz]  (in):  [xyz] coordinates of local (narrow band) points
*    n*_index    (in):  index range of points to loop over in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*
*/  
void LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER2_LOCAL(
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
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const int *klo_nb_gb,
  const int *khi_nb_gb,
  const unsigned char *mark_fb);   


/*!
*
*  LSM3D_COMPUTE_MEAN_CURVATURE_ORDER4_LOCAL() computes mean curvature
*  kappa = div ( grad_phi / |grad_phi|)
*  kappa = ( phi_xx*phi_y^2 + phi_yy*phi_x^2 - 2*phi_xy*phi_x*phi_y + 
*            phi_xx*phi_z^2 + phi_zz*phi_x^2 - 2*phi_xz*phi_x*phi_z +
*            phi_yy*phi_z^2 + phi_zz*phi_y^2 - 2*phi_yz*phi_y*phi_z )/
*          ( | grad phi | ^ 3 )
*  Note that this value is technically twice the mean curvature.
*  4th order central differencing used.
*  First order derivatives assumed precomputed(with 4th order).
c
*  Arguments:
*    kappa     (in/out): curvature data array
*    phi          (in):  level set function
*    phi_*        (in):  first order derivatives of phi
*    grad_phi_mag (in):  gradient magnitude of phi
*    *_gb        (in):   index range for ghostbox
*    dx, dy      (in):   grid spacing
*    index_[xyz]  (in):  [xyz] coordinates of local (narrow band) points
*    n*_index    (in):  index range of points to loop over in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*
*/ 
void LSM3D_COMPUTE_MEAN_CURVATURE_ORDER4_LOCAL(
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
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const int *klo_nb_gb,
  const int *khi_nb_gb,
  const unsigned char *mark_fb);
 
/*!
*
*  LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER4_LOCAL() computes Gaussian curvature
*  kappa = [  phi_x^2*(phi_yy*phi_zz - phi_yz^2) +
*             phi_y^2*(phi_xx*phi_zz - phi_xz^2) + 
*             phi_z^2*(phi_xx*phi_yy - phi_xy^2) + 
*         2*( phi_x*phi_y*(phi_xz*phi_yz - phi_xy*phi_zz) +
*             phi_y*phi_z*(phi_xy*phi_xz - phi_yz*phi_xx) + 
*             phi_x*phi_z*(phi_xy*phi_yz - phi_xz*phi_yy) ) ] / 
*          ( | grad phi | ^ 4 )
*  4th order central differencing used.
*  First order derivatives assumed precomputed (with 4th order).
*
*  Arguments:
*    kappa     (in/out): curvature data array
*    phi          (in):  level set function
*    phi_*        (in):  first order derivatives of phi
*    grad_phi_mag (in):  gradient magnitude of phi
*    *_gb        (in):   index range for ghostbox
*    dx, dy      (in):   grid spacing
*    index_[xyz]  (in):  [xyz] coordinates of local (narrow band) points
*    n*_index    (in):  index range of points to loop over in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*
*/  
void LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER4_LOCAL(
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
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const int *klo_nb_gb,
  const int *khi_nb_gb,
  const unsigned char *mark_fb); 
#ifdef __cplusplus
}
#endif
  
#endif 
