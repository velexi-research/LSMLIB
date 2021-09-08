/*
 * File:        lsm_spatial_derivatives3d_local.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for Fortran 77 3D narrow-band ENO/WENO routines
 */

#ifndef INCLUDED_LSM_SPATIAL_DERIVATIVES_3D_LOCAL_H
#define INCLUDED_LSM_SPATIAL_DERIVATIVES_3D_LOCAL_H

#include "lsmlib_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_spatial_derivatives3d_local.h
 *
 * \brief
 * @ref lsm_spatial_derivatives3d_local.h provides support for computing spatial
 * derivatives in three space dimensions using high-order ENO and WENO 
 * discretizations for narrow banding (local) approach. 
 *
 */
 

/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define LSM3D_HJ_ENO1_LOCAL              lsm3dhjeno1local_
#define LSM3D_HJ_ENO2_LOCAL              lsm3dhjeno2local_

#define LSM3D_CENTRAL_GRAD_ORDER2_LOCAL  lsm3dcentralgradorder2local_
#define LSM3D_CENTRAL_GRAD_ORDER4_LOCAL  lsm3dcentralgradorder4local_

#define LSM3D_LAPLACIAN_ORDER2_LOCAL     lsm3dlaplacianorder2local_
#define LSM3D_COMPUTE_AVE_GRAD_PHI_LOCAL lsm3dcomputeavegradphilocal_
#define LSM3D_GRADIENT_MAGNITUDE_LOCAL   lsm3dgradientmagnitudelocal_


#ifdef __cplusplus
extern "C" {
#endif

/*!
*
*  LSM3D_HJ_ENO1_LOCAL() computes the forward (plus) and backward (minus)
*  first-order Hamilton-Jacobi ENO approximations to the gradient of 
*  phi.
*  The routine loops only over local (narrow band) points.
*
*  Arguments:
*    phi_*_plus (out):   components of grad(phi) in plus direction 
*    phi_*_minus (out):  components of grad(phi) in minus direction
*    phi (in):           phi
*    D1 (in):            scratch space for holding undivided first-differences
*    dx, dy, dz (in):    grid spacing
*    *_gb (in):          index range for ghostbox
*    index_*(in):        coordinates of local (narrow band) points
*    n*_index[01](in):   index range of points in index_* that are in
*                        level [01] of the narrow band
*    narrow_band(in):    array that marks voxels outside desired fillbox
*    mark_*(in):         upper limit narrow band value for voxels in 
*                        the appropriate fillbox
*
*  NOTES:
*   - it is assumed that BOTH the plus AND minus derivatives have
*     the same fillbox
*   - index_[xyz] arrays range at minimum from nlo_index0 to nhi_index1
*/
				   
void LSM3D_HJ_ENO1_LOCAL(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  LSMLIB_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,
  const int *khi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  LSMLIB_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,
  const int *khi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index0,
  const int *nhi_index0,
  const int *nlo_index1,
  const int *nhi_index1,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const int *klo_nb_gb,
  const int *khi_nb_gb,
  const unsigned char *mark_fb,
  const unsigned char *mark_D1);

/*!
*
*  LSM3D_HJ_ENO2_LOCAL() computes the forward (plus) and backward (minus)
*  second-order Hamilton-Jacobi ENO approximations to the gradient of 
*  phi.
*  The routine loops only over local (narrow band) points.
*
*  Arguments:
*    phi_*_plus (out):   components of grad(phi) in plus direction 
*    phi_*_minus (out):  components of grad(phi) in minus direction
*    phi (in):           phi
*    D1 (in):            scratch space for holding undivided first-differences
*    D2 (in):            scratch space for holding undivided second-differences
*    dx, dy, dz (in):    grid spacing
*    *_gb (in):          index range for ghostbox
*    index_*(in):        coordinates of local (narrow band) points
*    n*_index[012](in):  index range of points in index_* that are in
*                        level [012] of the narrow band
*    narrow_band(in):    array that marks voxels outside desired fillbox
*    mark_*(in):         upper limit narrow band value for voxels in 
*                        the appropriate fillbox
*
*  NOTES:
*   - it is assumed that BOTH the plus and minus derivatives have
*     the same fillbox
*   - index_* arrays range at minimum from nlo_index0 to nhi_index2
*/

void LSM3D_HJ_ENO2_LOCAL(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  LSMLIB_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,
  const int *khi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  LSMLIB_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,
  const int *khi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const int *klo_D2_gb,
  const int *khi_D2_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index0,
  const int *nhi_index0,
  const int *nlo_index1,
  const int *nhi_index1,
  const int *nlo_index2,
  const int *nhi_index2,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const int *klo_nb_gb,
  const int *khi_nb_gb,
  const unsigned char *mark_fb,
  const unsigned char *mark_D1,
  const unsigned char *mark_D2);

/*!
*
*  LSM3D_CENTRAL_GRAD_ORDER2_LOCAL() computes the second-order central 
*  approximation to the gradient of phi.
*  The routine loops only over local (narrow band) points.
*
*  Arguments:
*    phi_* (out):      components of grad(phi) 
*    phi (in):         phi
*    dx, dy, dz (in):  grid spacing
*    *_gb (in):        index range for ghostbox
*    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
*    n*_index(in):     index range of points to loop over in index_*
*    narrow_band(in):  array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*
*/
void LSM3D_CENTRAL_GRAD_ORDER2_LOCAL( 
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  LSMLIB_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
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
*  LSM3D_CENTRAL_GRAD_ORDER4_LOCAL() computes the fourth-order, central,
*  finite difference approximation to the gradient of phi.
*  The routine loops only over local (narrow band) points.
*
*  Arguments:
*    phi_* (out):      components of grad(phi) 
*    phi (in):         phi
*    dx, dy, dz (in):  grid spacing
*    *_gb (in):        index range for ghostbox
*    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
*    n*_index(in):     index range of points to loop over in index_*
*    narrow_band(in):  array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*
*/
void LSM3D_CENTRAL_GRAD_ORDER4_LOCAL( 
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  LSMLIB_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
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
*  LSM3D_LAPLACIAN_ORDER2_LOCAL() computes the second-order, central, 
*  finite difference approximation to the Laplacian of phi.
*  The routine loops only over local (narrow band) points.
*
*  Arguments:
*    laplacian_phi (out):  Laplacian of phi
*    phi (in):             phi
*    dx (in):              grid spacing
*    *_gb (in):            index range for ghostbox
*    index_[xyz](in):     [xyz] coordinates of local (narrow band) points
*    n*_index(in):        index range of points to loop over in index_*
*    narrow_band(in):     array that marks voxels outside desired fillbox
*    mark_fb(in):         upper limit narrow band value for voxels in 
*                         fillbox
*
*/
void LSM3D_LAPLACIAN_ORDER2_LOCAL( 
  LSMLIB_REAL *laplacian_phi,
  const int *ilo_laplacian_phi_gb,
  const int *ihi_laplacian_phi_gb,
  const int *jlo_laplacian_phi_gb,
  const int *jhi_laplacian_phi_gb,
  const int *klo_laplacian_phi_gb,
  const int *khi_laplacian_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
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
*  LSM3D_COMPUTE_AVE_GRAD_PHI_LOCAL() computes the average of the 
*  second-order, central, finite difference approximation to the 
*  gradient of phi within the narrow band.
*  The routine loops only over local (narrow band) points. 
*
*  Arguments:
*    phi (in):          phi
*    grad_phi_ave(out): average of the gradient 
*    dx, dy, dz (in):   grid spacing
*    index_*(in):       coordinates of local (narrow band) points
*    n*_index(in):     index range of points to loop over in index_*
*    narrow_band(in):  array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*    *_gb (in):        index range for ghostbox
*
*/     
void LSM3D_COMPUTE_AVE_GRAD_PHI_LOCAL(
  LSMLIB_REAL *grad_phi_ave,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
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
*  LSM3D_GRADIENT_MAGNITUDE_LOCAL() computes magnitude of the gradient of phi.
*
*  Arguments:
*    phi_* (in):         components of grad(phi) 
*    grad_phi_mag (out): gradient magnitude
*    *_gb (in):          index range for ghostbox
*    index_*(in):        coordinates of local (narrow band) points
*    n*_index(in):       index range of points to loop over in index_*
*    narrow_band(in):    array that marks voxels outside desired fillbox
*    mark_fb(in):        upper limit narrow band value for voxels in 
*                        fillbox
*
*/
void LSM3D_GRADIENT_MAGNITUDE_LOCAL(
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
