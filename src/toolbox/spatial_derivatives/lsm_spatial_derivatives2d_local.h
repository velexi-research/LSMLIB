/*
 * File:        lsm_spatial_derivatives2d_local.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.19 $
 * Modified:    $Date$
 * Description: Header file for Fortran 77 2D narrow-band ENO/WENO routines
 */

#ifndef INCLUDED_LSM_SPATIAL_DERIVATIVES_2D_LOCAL_H
#define INCLUDED_LSM_SPATIAL_DERIVATIVES_2D_LOCAL_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_spatial_derivatives2d_local.h
 *
 * \brief
 * @ref lsm_spatial_derivatives2d_local.h provides support for computing spatial
 * derivatives in two space dimensions using high-order ENO and WENO 
 * discretizations for narrow banding (local) approach. 
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define LSM2D_HJ_ENO1_LOCAL              lsm2dhjeno1local_
#define LSM2D_HJ_ENO2_LOCAL              lsm2dhjeno2local_
#define LSM2D_HJ_ENO3_LOCAL              lsm2dhjeno3local_
#define LSM2D_HJ_WENO5_LOCAL             lsm2dhjweno5local_

#define LSM2D_UPWIND_HJ_ENO2_LOCAL         lsm2dupwindhjeno2local_

#define LSM2D_CENTRAL_GRAD_ORDER2_LOCAL  lsm2dcentralgradorder2local_
#define LSM2D_CENTRAL_GRAD_ORDER4_LOCAL  lsm2dcentralgradorder4local_
#define LSM2D_LAPLACIAN_ORDER2_LOCAL     lsm2dlaplacianorder2local_
#define LSM2D_COMPUTE_AVE_GRAD_PHI_LOCAL lsm2dcomputeavegradphilocal_
#define LSM2D_GRADIENT_MAGNITUDE_LOCAL   lsm2dgradientmagnitudelocal_
#define LSM2D_DIVERGENCE_CENTRAL_LOCAL   lsm2ddivergencecentrallocal_

/*!
*
*  LSM2D_HJ_ENO1_LOCAL() computes the forward (plus) and backward (minus)
*  first-order Hamilton-Jacobi ENO approximations to the gradient of 
*  phi.
*  The routine loops only over local (narrow band) points.
*
*  Arguments:
*    phi_*_plus (out):   components of grad(phi) in plus direction 
*    phi_*_minus (out):  components of grad(phi) in minus direction
*    phi (in):           phi
*    D1 (in):            scratch space for holding undivided first-differences
*    dx, dy(in):         grid spacing
*    *_gb (in):          index range for ghostbox
*    index_[xy](in):     [xy] coordinates of local (narrow band) points
*    n*_index[01](in):   index range of points in index_* that are in
*                        level [01] of the narrow band
*    narrow_band(in):    array that marks voxels outside desired fillbox
*    mark_*(in):         upper limit narrow band value for voxels in 
*                        the appropriate fillbox
*
*  NOTES:
*   - it is assumed that BOTH the plus AND minus derivatives have
*     the same fillbox
*   - index_[xy] arrays range at minimum from nlo_index0 to nhi_index1
*
*/
void LSM2D_HJ_ENO1_LOCAL(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
  const int *nlo_index0,
  const int *nhi_index0,
  const int *nlo_index1,
  const int *nhi_index1,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb,
  const unsigned char *mark_D1);
  
  
/*!
*
*  LSM2D_HJ_ENO2_LOCAL() computes the forward (plus) and backward (minus)
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
*    dx, dy (in):        grid spacing
*    *_gb (in):          index range for ghostbox
*    index_[xy](in):     [xy] coordinates of local (narrow band) points
*    n*_index[012](in):  index range of points in index_* that are in
*                        level [012] of the narrow band
*    narrow_band(in):    array that marks voxels outside desired fillbox
*    mark_*(in):         upper limit narrow band value for voxels in 
*                        the appropriate fillbox
*
*  NOTES:
*   - it is assumed that BOTH the plus and minus derivatives have
*     the same fillbox
*   - index_[xy] arrays range at minimum from nlo_index0 to nhi_index2
*
*/
void LSM2D_HJ_ENO2_LOCAL(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
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
  const unsigned char *mark_fb,
  const unsigned char *mark_D1,
  const unsigned char *mark_D2);
  
  
/*!
*
*  LSM2D_HJ_ENO3_LOCAL() computes the forward (plus) and backward (minus)
*  third-order Hamilton-Jacobi ENO approximations to the gradient of 
*  phi.
*
*  Arguments:
*    phi_*_plus (out):   components of grad(phi) in plus direction
*    phi_*_minus (out):  components of grad(phi) in minus direction
*    phi (in):           phi
*    D1 (in):            scratch space for holding undivided first-differences
*    D2 (in):            scratch space for holding undivided second-differences
*    D3 (in):            scratch space for holding undivided third-differences
*    dx, dy (in):        grid spacing
*    *_gb (in):          index range for ghostbox
*    index_[xy](in):     [xy] coordinates of local (narrow band) points
*    n*_index[0123](in):  index range of points in index_* that are in
*                        level [0123] of the narrow band
*    narrow_band(in):    array that marks voxels outside desired fillbox
*    mark_*(in):         upper limit narrow band value for voxels in 
*                        the appropriate fillbox
*
*  NOTES:
*   - it is assumed that BOTH the plus AND minus derivatives have
*     the same fillbox
*
*/
void LSM2D_HJ_ENO3_LOCAL(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  LSMLIB_REAL *D3,
  const int *ilo_D3_gb,
  const int *ihi_D3_gb,
  const int *jlo_D3_gb,
  const int *jhi_D3_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
  const int *nlo_index0,
  const int *nhi_index0,
  const int *nlo_index1,
  const int *nhi_index1,
  const int *nlo_index2,
  const int *nhi_index2,
  const int *nlo_index3,
  const int *nhi_index3,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb,
  const unsigned char *mark_D1,
  const unsigned char *mark_D2,
  const unsigned char *mark_D3);
  
/*!
*
*  lsm2dHJWENO5LOCAL() computes the forward (plus) and backward (minus)
*  fifth-order Hamilton-Jacobi ENO approximations to the gradient of 
*  phi.
*
*  Arguments:
*    phi_*_plus (out):   components of grad(phi) in plus direction
*    phi_*_minus (out):  components of grad(phi) in minus direction
*    phi (in):           phi
*    D1 (in):            scratch space for holding undivided first-differences
*    dx, dy (in):        grid spacing
*    *_gb (in):          index range for ghostbox
*    *_fb (in):          index range for fillbox
*    index_[xy](in):     [xy] coordinates of local (narrow band) points
*    n*_index[0123](in):  index range of points in index_* that are in
*                        level [0123] of the narrow band
*    narrow_band(in):    array that marks voxels outside desired fillbox
*    mark_*(in):         upper limit narrow band value for voxels in 
*                        the appropriate fillbox
*
*  NOTES:
*   - it is assumed that BOTH the plus AND minus derivatives have
*     the same fillbox
*
*/
void LSM2D_HJ_WENO5_LOCAL(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
  const int *nlo_index0,
  const int *nhi_index0,
  const int *nlo_index1,
  const int *nhi_index1,
  const int *nlo_index2,
  const int *nhi_index2,
  const int *nlo_index3,
  const int *nhi_index3,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb,
  const unsigned char *mark_D1);

/*!
*
*
* LSM2D_UPWIND_HJ_ENO2_LOCAL() computes the second-order Hamilton-Jacobi ENO 
*  upwind approximation to the gradient of phi.
*
*  Arguments:
*    phi_* (out):  components of grad(phi)
*    phi (in):     phi
*    vel_* (in):   components of the velocity
*    D1 (in):      scratch space for holding undivided first-differences
*    D2 (in):      scratch space for holding undivided second-differences
*    dx, dy (in):  grid spacing
*    *_gb (in):    index range for ghostbox
*    *_fb (in):    index range for fillbox
*
*/
void LSM2D_UPWIND_HJ_ENO2_LOCAL(
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const LSMLIB_REAL *vel_x,
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  const int *jlo_vel_gb,
  const int *jhi_vel_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const int *index_x,
  const int *index_y, 
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
  const unsigned char *mark_fb,
  const unsigned char *mark_D1,
  const unsigned char *mark_D2);

/*!
*
*  LSM2D_CENTRAL_GRAD_ORDER2_LOCAL() computes the second-order central 
*  approximation to the gradient of phi.
*  The routine loops only over local (narrow band) points.
*
*  Arguments:
*    phi_* (out):      components of grad(phi) 
*    phi (in):         phi
*    dx, dy (in):      grid spacing
*    *_gb (in):        index range for ghostbox
*    index_[xy](in):   [xy] coordinates of local (narrow band) points
*    n*_index(in):     index range of points to loop over in index_*
*    narrow_band(in):  array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*
*/
void LSM2D_CENTRAL_GRAD_ORDER2_LOCAL( 
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
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
*  LSM2D_CENTRAL_GRAD_ORDER4_LOCAL() computes the second-order, central, 
*  finite difference approximation to the gradient of phi.
*  The routine loops only over local (narrow band) points.
*
*  Arguments:
*    phi_* (out):      components of grad(phi) 
*    phi (in):         phi
*    dx, dy (in):      grid spacing
*    *_gb (in):        index range for ghostbox
*    index_[xy](in):   [xy] coordinates of local (narrow band) points
*    n*_index(in):     index range of points to loop over in index_*
*    narrow_band(in):  array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*
*/  
void LSM2D_CENTRAL_GRAD_ORDER4_LOCAL( 
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
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
*  LSM2D_LAPLACIAN_ORDER2_LOCAL() computes the second-order, central, 
*  finite difference approximation to the Laplacian of phi.
*  The routine loops only over local (narrow band) points.
*
*  Arguments:
*    laplacian_phi (out):  Laplacian of phi
*    phi (in):             phi
*    dx (in):              grid spacing
*    *_gb (in):            index range for ghostbox
*    index_[xy](in):       [xy] coordinates of local (narrow band) points
*    n*_index(in):         index range of points to loop over in index_*
*    narrow_band(in):      array that marks voxels outside desired fillbox
*    mark_fb(in):          upper limit narrow band value for voxels in 
*                          fillbox
*
*/
void LSM2D_LAPLACIAN_ORDER2_LOCAL( 
  LSMLIB_REAL *laplacian_phi,
  const int *ilo_laplacian_phi_gb,
  const int *ihi_laplacian_phi_gb,
  const int *jlo_laplacian_phi_gb,
  const int *jhi_laplacian_phi_gb,
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
*  LSM2D_COMPUTE_AVE_GRAD_PHI_LOCAL() computes the average of the 
*  second-order, central, finite difference approximation to the 
*  magnitude of gradient of phi within the narrow band.
*  The routine loops only over local (narrow band) points. 
*
*  Arguments:
*    phi (in):          phi
*    grad_phi_ave(out): average of the gradient 
*    dx, dy(in):        grid spacing
*    index_[xy](in):    [xy] coordinates of local (narrow band) points
*    n*_index(in):      index range of points to loop over in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):       upper limit narrow band value for voxels in 
*                       fillbox
*    *_gb (in):         index range for ghostbox
*
*/
void LSM2D_COMPUTE_AVE_GRAD_PHI_LOCAL(
  LSMLIB_REAL *grad_phi_ave,
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
*  LSM2D_GRADIENT_MAGNITUDE_LOCAL() computes magnitude of the gradient of phi.
*
*  Arguments:
*    phi_* (in):           components of grad(phi) 
*    grad_phi_mag (out):   gradient magnitude
*    *_gb (in):           index range for ghostbox
*    index_[xy](in):     [xy] coordinates of local (narrow band) points
*    n*_index(in):       index range of points to loop over in index_*
*    narrow_band(in):    array that marks voxels outside desired fillbox
*    mark_fb(in):        upper limit narrow band value for voxels in 
*                        fillbox
*
*/
 void LSM2D_GRADIENT_MAGNITUDE_LOCAL(   
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  LSMLIB_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
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
*  LSM2D_DIVERGENCE_CENTRAL_LOCAL() computes the second-order, central, 
*  finite difference approximation to the divergence of a vector field.
*
*  Arguments:
*    divF* (out):  divergence of F
*    FX, FY(in):   x and y components of vector field F 
*    dx, dy (in):  grid spacing
*    *_gb (in):    index range for ghostbox
*    *_fb (in):    index range for fillbox
*    index_[xy](in):  [xy] coordinates of local (narrow band) points
*    n*_index(in):    index range of points to loop over in index_*
*    narrow_band(in): array that marks voxels outside desired fillbox
*    mark_fb(in):     upper limit narrow band value for voxels in 
*                     fillbox
*
*/
 void  LSM2D_DIVERGENCE_CENTRAL_LOCAL(
  LSMLIB_REAL *divF,
  const int *ilo_divf_gb, 
  const int *ihi_divf_gb,
  const int *jlo_divf_gb, 
  const int *jhi_divf_gb,
  const LSMLIB_REAL *FX,
  const LSMLIB_REAL *FY,
  const int *ilo_gb, 
  const int *ihi_gb,
  const int *jlo_gb,  
  const int *jhi_gb,
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
