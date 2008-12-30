/*
 * File:        lsm_utilities3d_local.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.23 $
 * Modified:    $Date: 2006/10/28 04:54:37 $
 * Description: Header file for 3D Fortran 77 level set method narrow-band
 *              utility subroutines
 */

#ifndef INCLUDED_LSM_UTILITIES_3D_LOCAL_H
#define INCLUDED_LSM_UTILITIES_3D_LOCAL_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_utilities3d.h
 *
 * \brief 
 * @ref lsm_utilities3d.h provides several utility functions that support
 * level set method calculations in three space dimensions.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                        name in
 *      C/C++ code                     Fortran code
 *      ----------                     ------------
 */
#define LSM3D_MAX_NORM_DIFF_LOCAL      lsm3dmaxnormdifflocal_

#define LSM3D_COMPUTE_STABLE_ADVECTION_DT_LOCAL                           \
                                       lsm3dcomputestableadvectiondtlocal_
#define LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_LOCAL                                  \
                                       lsm3dcomputestablenormalveldtlocal_
#define LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME_LOCAL            \
                          lsm3dcomputestablenormalveldtcontrolvolumelocal_
#define LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_LOCAL                     \
                                  lsm3dcomputestableconstnormalveldtlocal_

/*!
*
*  LSM3D_MAX_NORM_DIFF_LOCAL() computes the max norm of the difference 
*  between the two specified scalar fields.
*  The routine loops only over local (narrow band) points. 
*
*  Arguments:
*    max_norm_diff (out):  max norm of the difference between the fields
*    field1 (in):          scalar field 1
*    field2 (in):          scalar field 2
*    index_*](in):         coordinates of local (narrow band) points
*    n*_index(in):         index range of points in index_*
*    *_gb (in):            index range for ghostbox
*    narrow_band(in):      array that marks voxels outside desired fillbox
*    mark_fb(in):          upper limit narrow band value for voxels in 
*                          fillbox
*
*/
void LSM3D_MAX_NORM_DIFF_LOCAL(
  LSMLIB_REAL *max_norm_diff,
  const LSMLIB_REAL *field1,
  const int *ilo_field1_gb, 
  const int *ihi_field1_gb,
  const int *jlo_field1_gb, 
  const int *jhi_field1_gb,
  const int *klo_field1_gb, 
  const int *khi_field1_gb,
  const LSMLIB_REAL *field2,
  const int *ilo_field2_gb, 
  const int *ihi_field2_gb,
  const int *jlo_field2_gb, 
  const int *jhi_field2_gb,
  const int *klo_field2_gb, 
  const int *khi_field2_gb,
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
*  LSM3D_COMPUTE_STABLE_ADVECTION_DT_LOCAL() computes the stable time step size 
*  for an advection term based on a CFL criterion.
*  
*  Arguments:
*    dt (out):           step size
*    vel_* (in):         components of velocity at t = t_cur
*    *_gb (in):          index range for ghostbox
*    dx, dy (in):        grid spacing
*    index_*(in):        coordinates of local (narrow band) points
*    n*_index(in):       index range of points in index_*
*    narrow_band(in):    array that marks voxels outside desired fillbox
*    mark_fb(in):        upper limit narrow band value for voxels in 
*                        fillbox
*
*/
void LSM3D_COMPUTE_STABLE_ADVECTION_DT_LOCAL(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_x,
  const LSMLIB_REAL *vel_y,
  const LSMLIB_REAL *vel_z,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *klo_vel_gb, 
  const int *khi_vel_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const LSMLIB_REAL *cfl_number,
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
*  LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_LOCAL() computes the stable time step 
*  size for a normal velocity term based on a CFL criterion.
*  
*  Arguments:
*    dt (out):             step size
*    vel_n (in):           normal velocity at t = t_cur
*    phi_*_plus (in):      components of forward approx to grad(phi) at 
*                          t = t_cur
*    phi_*_minus (in):     components of backward approx to grad(phi) at
*                          t = t_cur
*    *_gb (in):            index range for ghostbox
*    d*(in):               grid spacing
*    index_*(in):          coordinates of local (narrow band) points
*    n*_index(in):         index range of points in index_*
*    narrow_band(in):      array that marks voxels outside desired fillbox
*    mark_fb(in):          upper limit narrow band value for voxels in 
*                          fillbox
*
*  NOTES:
*   - max(phi_*_plus , phi_*_minus) is the value of phi_* that is 
*     used in the time step size calculation.  This may be more 
*     conservative than necessary for Godunov's method, but it is 
*     cheaper to compute.
*
*/
void LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_LOCAL(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *klo_vel_gb, 
  const int *khi_vel_gb,
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
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const LSMLIB_REAL *cfl_number,
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
*  LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME_LOCAL() computes the stable time step 
*  size for a normal velocity term based on a CFL criterion.
*  
*  Arguments:
*    dt (out):             step size
*    vel_n (in):           normal velocity at t = t_cur
*    phi_*_plus (in):      components of forward approx to grad(phi) at 
*                          t = t_cur
*    phi_*_minus (in):     components of backward approx to grad(phi) at
*                          t = t_cur
*    *_gb (in):            index range for ghostbox
*    d*(in):               grid spacing
*    index_*(in):          coordinates of local (narrow band) points
*    n*_index(in):         index range of points in index_*
*    narrow_band(in):      array that marks voxels outside desired fillbox
*    mark_fb(in):          upper limit narrow band value for voxels in 
*                          fillbox
*
*  NOTES:
*   - max(phi_*_plus , phi_*_minus) is the value of phi_* that is 
*     used in the time step size calculation.  This may be more 
*     conservative than necessary for Godunov's method, but it is 
*     cheaper to compute.
*
*/
void LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME_LOCAL(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *klo_vel_gb, 
  const int *khi_vel_gb,
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
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *klo_control_vol_gb, 
  const int *khi_control_vol_gb,
  const int *control_vol_sgn,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const LSMLIB_REAL *cfl_number,
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
* LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_LOCAL() computes the stable time 
* step size for a normal velocity term based on a CFL criterion.
*  
*  Arguments:
*    dt (out):          step size
*    vel_n (in):        normal velocity at t = t_cur, constant for all pts
*    phi_*_plus (in):   components of forward approx to grad(phi) at 
*                       t = t_cur
*    phi_*_minus (in):  components of backward approx to grad(phi) at
*                       t = t_cur
*    index_*(in):       coordinates of local (narrow band) points
*    n*_index(in):      index range of points in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):       upper limit narrow band value for voxels in 
*                       fillbox
*    *_gb (in):         index range for ghostbox
*    d* (in):           grid spacing
*
*  NOTES:
*   - max(phi_*_plus , phi_*_minus) is the value of phi_* that is 
*     used in the time step size calculation.  This may be more 
*     conservative than necessary for Godunov's method, but it is 
*     cheaper to compute.
*
*/     
void LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_LOCAL(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
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
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const LSMLIB_REAL *cfl_number,  
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
