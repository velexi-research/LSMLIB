/*
 * File:        lsm_utilities2d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for 2D Fortran 77 level set method utility 
 *              subroutines
 */

#ifndef INCLUDED_LSM_UTILITIES_2D_H
#define INCLUDED_LSM_UTILITIES_2D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_utilities2d.h
 *
 * \brief 
 * @ref lsm_utilities2d.h provides several utility functions that support
 * level set method calculations in two space dimensions.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                        name in
 *      C/C++ code                     Fortran code
 *      ----------                     ------------
 */
#define LSM2D_MAX_NORM_DIFF                   lsm2dmaxnormdiff_
#define LSM2D_AVE_ABS_DIFF                    lsm2daveabsdiff_
#define LSM2D_VOXEL_COUNT_GREATER_THAN_ZERO   lsm2dvoxelcountgreaterthanzero_
#define LSM2D_VOXEL_COUNT_LESS_THAN_ZERO      lsm2dvoxelcountlessthanzero_

#define LSM2D_COMPUTE_STABLE_ADVECTION_DT                                     \
                                       lsm2dcomputestableadvectiondt_
#define LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT                                    \
                                       lsm2dcomputestablenormalveldt_
#define LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT                              \
                                       lsm2dcomputestableconstnormalveldt_
#define LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO                              \
                                       lsm2dvolumeintegralphilessthanzero_
#define LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO                           \
                                       lsm2dvolumeintegralphigreaterthanzero_
#define LSM2D_SURFACE_INTEGRAL         lsm2dsurfaceintegral_
#define LSM2D_SURFACE_INTEGRAL_DELTA                            \
                       lsm2dsurfaceintegralprecomputeddelta_

#define LSM2D_MAX_NORM_DIFF_CONTROL_VOLUME                                    \
                       lsm2dmaxnormdiffcontrolvolume_
#define LSM2D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME                      \
                       lsm2dcomputestableadvectiondtcontrolvolume_
#define LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME                    \
                       lsm2dcomputestablenormalveldtcontrolvolume_
#define LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_CONTROL_VOLUME               \
                       lsm2dcomputestableconstnormalveldtcontrolvolume_
#define LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME               \
                       lsm2dvolumeintegralphilessthanzerocontrolvolume_
#define LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME            \
                       lsm2dvolumeintegralphigreaterthanzerocontrolvolume_
#define LSM2D_SURFACE_INTEGRAL_CONTROL_VOLUME                                 \
                       lsm2dsurfaceintegralcontrolvolume_
#define LSM2D_VOXEL_COUNT_GREATER_THAN_ZERO_CONTROL_VOLUME \
                                   lsm2dvoxelcountgreaterthanzerocontrolvolume_
#define LSM2D_VOXEL_COUNT_LESS_THAN_ZERO_CONTROL_VOLUME  \
                                   lsm2dvoxelcountlessthanzerocontrolvolume_
#define LSM2D_SURFACE_INTEGRAL_DELTA_CONTROL_VOLUME              \
                       lsm2dsurfaceintegralprecomputeddeltacontrolvolume_
		       
/*!
 * LSM2D_MAX_NORM_DIFF() computes the max norm of the difference
 * between the two specified scalar fields.
 *      
 * Arguments:
 *  - max_norm_diff (out):  max norm of the difference between the fields
 *  - field1 (in):          scalar field 1
 *  - field2 (in):          scalar field 2
 *  - *_gb (in):            index range for ghostbox
 *  - *_ib (in):            index range for box to include in norm
 *                          calculation
 *
 * Return value:            none
 */
void LSM2D_MAX_NORM_DIFF(
  LSMLIB_REAL *max_norm_diff,
  const LSMLIB_REAL *field1,
  const int *ilo_field1_gb, 
  const int *ihi_field1_gb,
  const int *jlo_field1_gb, 
  const int *jhi_field1_gb,
  const LSMLIB_REAL *field2,
  const int *ilo_field2_gb, 
  const int *ihi_field2_gb,
  const int *jlo_field2_gb, 
  const int *jhi_field2_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib);

/*!
*
*  LSM2D_AVE_ABS_DIFF() computes the average pointwise abs. difference 
*  between two scalar fields. 
*
*  Arguments:
*    ave_abs_diff (out):    average of the difference between the fields
*    field1 (in):           scalar field 1
*    field2 (in):           scalar field 2
*    *_gb (in):             index range for ghostbox
*    *_ib (in):             index range for box to include in 
*                           calculation
*
*/
void LSM2D_AVE_ABS_DIFF(
  LSMLIB_REAL *ave_abs_diff,
  const LSMLIB_REAL *field1,
  const int *ilo_field1_gb, 
  const int *ihi_field1_gb,
  const int *jlo_field1_gb, 
  const int *jhi_field1_gb,
  const LSMLIB_REAL *field2,
  const int *ilo_field2_gb, 
  const int *ihi_field2_gb,
  const int *jlo_field2_gb, 
  const int *jhi_field2_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib);
  
  
/*!
*
*  LSM2D_VOXEL_COUNT_GREATER_THAN_ZERO() computes number of voxels whose phi
*  value is greater than zero.
*
*  Arguments:
*    count (out):           number of voxels where phi < 0
*    phi (in):              scalar field
*    *_gb (in):             index range for ghostbox
*    *_ib (in):             index range for box to include in
*                           calculation
*
*/  
void LSM2D_VOXEL_COUNT_GREATER_THAN_ZERO(
  int *count,
  const LSMLIB_REAL *phi,
  const int *ilo_gb,
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb);
  
  
/*!
*
*  LSM2D_VOXEL_COUNT_LESS_THAN_ZERO() computes number of voxels whose phi
*  value is less than zero.
*
*  Arguments:
*    count (out):           number of voxels where phi < 0
*    phi (in):              scalar field
*    *_gb (in):             index range for ghostbox
*    *_ib (in):             index range for box to include in
*                           calculation
*
*/  
void LSM2D_VOXEL_COUNT_LESS_THAN_ZERO(
  int *count,
  const LSMLIB_REAL *phi,
  const int *ilo_gb,
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb);  
  

/*!
 * LSM2D_COMPUTE_STABLE_ADVECTION_DT() computes the stable time step size
 * for an advection term based on a CFL criterion.
 *
 * Arguments:
 *  - dt (out):              step size
 *  - vel_* (in):            components of velocity at t = t_cur
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for box to include dt calculation
 *  - dx (in):               grid spacing
 * 
 * Return value:             none 
 */
void LSM2D_COMPUTE_STABLE_ADVECTION_DT(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_x,
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *cfl_number);


/*!
 * LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT() computes the stable time step
 * size for a normal velocity term based on a CFL criterion.
 *  
 * Arguments:
 *  - dt (out):              step size
 *  - vel_n (in):            normal velocity at t = t_cur
 *  - phi_*_plus (in):       components of forward approx to 
 *                           \f$ \nabla \phi \f$ at t = t_cur
 *  - phi_*_minus (in):      components of backward approx to 
 *                           \f$ \nabla \phi \f$ at t = t_cur
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for box to include dt calculation
 *  - dx (in):               grid spacing
 *   
 * Return value:             none   
 *   
 * NOTES:
 *  - max(phi_*_plus , phi_*_minus) is the value of phi_* that is
 *    used in the time step size calculation.  This may be more
 *    conservative than necessary for Godunov's method, but it is
 *    cheaper to compute.
 *
 */
void LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const LSMLIB_REAL *phi_x_plus,
  const LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus,
  const LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *cfl_number);


/*!
 * LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT() computes the stable time 
 * step size for a constant normal velocity term based on a CFL criterion.
 *  
 * Arguments:
 *  - dt (out):              step size
 *  - vel_n (in):            constant normal velocity at t = t_cur
 *  - phi_*_plus (in):       components of forward approx to 
 *                           \f$ \nabla \phi \f$ at t = t_cur
 *  - phi_*_minus (in):      components of backward approx to 
 *                           \f$ \nabla \phi \f$ at t = t_cur
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for box to include dt calculation
 *  - dx (in):               grid spacing
 *   
 * Return value:             none   
 *   
 * NOTES:
 *  - max(phi_*_plus , phi_*_minus) is the value of phi_* that is
 *    used in the time step size calculation.  This may be more
 *    conservative than necessary for Godunov's method, but it is
 *    cheaper to compute.
 *
 */
void LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
  const LSMLIB_REAL *phi_x_plus,
  const LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus,
  const LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *cfl_number);
  
  
/*!
 * LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO() computes the volume integral of
 *  the specified function over the region where the level set function
 *  is less than 0.  
 *    
 * Arguments:
 *  - int_F (out):           value of integral of F over the region where 
 *                           \f$ \phi < 0 \f$
 *  - F (in):                function to be integrated
 *  - phi (in):              level set function
 *                           points should be used
 *  - dx (in):               grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 */
void LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const int *jlo_F_gb, 
  const int *jhi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO() computes the volume integral
 * of the specified function over the region where the level set
 * function is greater than 0.  
 *
 * Arguments:
 *  - int_F (out):           value of integral of F over the region where 
 *                           \f$ \phi > 0 \f$
 *  - F (in):                function to be integrated
 *  - phi (in):              level set function
 *  - dx (in):               grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 */
void LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const int *jlo_F_gb, 
  const int *jhi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM2D_SURFACE_INTEGRAL() computes the surface integral of the specified
 * function over the region where the level set function equals 0.
 *     
 * Arguments:
 *  - int_F (out):           value of integral of F over the region where 
 *                           \f$ \phi = 0 \f$
 *  - F (in):                function to be integrated
 *  - phi (in):              level set function
 *  - phi_* (in):            components of \f$ \nabla \phi \f$
 *  - dx (in):               grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           delta-function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 * 
 * Return value:             none
 */
void LSM2D_SURFACE_INTEGRAL(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const int *jlo_F_gb, 
  const int *jhi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb, 
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *epsilon);

/*!
 * LSM2D_SURFACE_INTEGRAL_DELTA() computes the surface integral 
 * of the specified function over the region of the computational domain
 * where the level set function equals 0, assuming delta function has been
 * precomputed. 
 *     
 * Arguments:
 *  - int_F (out):           value of integral of F over the region where 
 *                           \f$ \phi = 0 \f$
 *  - F (in):                function to be integrated
 *  - delta_phi (in):        delta function for the zero level set
 *  - grad_phi_mag (in):     gradient magnitude
 *  - d* (in):               grid spacing
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 * 
 * Return value:             none
 */  
void LSM2D_SURFACE_INTEGRAL_DELTA(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const int *jlo_F_gb, 
  const int *jhi_F_gb,
  const LSMLIB_REAL *delta_phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb, 
  const LSMLIB_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);
  
/*!
 * LSM2D_MAX_NORM_DIFF_CONTROL_VOLUME() computes the max norm of the 
 * difference between the two specified scalar fields in the region
 * of the computational domain included by the control volume data.
 *      
 * Arguments:
 *  - max_norm_diff (out):  max norm of the difference between the fields
 *  - field1 (in):          scalar field 1
 *  - field2 (in):          scalar field 2
 *  - control_vol (in):     control volume data (used to exclude cells
 *                          from the max norm calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - *_gb (in):            index range for ghostbox
 *  - *_ib (in):            index range for box to include in norm
 *                          calculation
 *
 * Return value:            none
 */
void LSM2D_MAX_NORM_DIFF_CONTROL_VOLUME(
  LSMLIB_REAL *max_norm_diff,
  const LSMLIB_REAL *field1,
  const int *ilo_field1_gb, 
  const int *ihi_field1_gb,
  const int *jlo_field1_gb, 
  const int *jhi_field1_gb,
  const LSMLIB_REAL *field2,
  const int *ilo_field2_gb, 
  const int *ihi_field2_gb,
  const int *jlo_field2_gb, 
  const int *jhi_field2_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib);


/*!
 * LSM2D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME() computes the stable 
 * time step size for an advection term based on a CFL criterion for
 * grid cells within the computational domain included by the control 
 * volume data.
 *
 * Arguments:
 *  - dt (out):              step size
 *  - vel_* (in):            components of velocity at t = t_cur
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for box to include dt calculation
 *  - dx (in):               grid spacing
 * 
 * Return value:             none 
 */
void LSM2D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_x,
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *cfl_number);


/*!
 * LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME() computes the 
 * stable time step size for a normal velocity term based on a CFL 
 * criterion for grid cells within the computational domain included
 * by the control volume data.
 *  
 * Arguments:
 *  - dt (out):              step size
 *  - vel_n (in):            normal velocity at t = t_cur
 *  - phi_*_plus (in):       components of forward approx to 
 *                           \f$ \nabla \phi \f$ at t = t_cur
 *  - phi_*_minus (in):      components of backward approx to 
 *                           \f$ \nabla \phi \f$ at t = t_cur
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for box to include dt calculation
 *  - dx (in):               grid spacing
 *   
 * Return value:             none   
 *   
 * NOTES:
 *  - max(phi_*_plus , phi_*_minus) is the value of phi_* that is
 *    used in the time step size calculation.  This may be more
 *    conservative than necessary for Godunov's method, but it is
 *    cheaper to compute.
 *
 */
void LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const LSMLIB_REAL *phi_x_plus,
  const LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus,
  const LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *cfl_number);


/*!
 * LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_CONTROL_VOLUME() computes the 
 * stable time step size for a constant normal velocity term based on 
 * a CFL criterion for grid cells within the computational domain 
 * included by the control volume data.
 *  
 * Arguments:
 *  - dt (out):              step size
 *  - vel_n (in):            constant normal velocity at t = t_cur
 *  - phi_*_plus (in):       components of forward approx to 
 *                           \f$ \nabla \phi \f$ at t = t_cur
 *  - phi_*_minus (in):      components of backward approx to 
 *                           \f$ \nabla \phi \f$ at t = t_cur
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for box to include dt calculation
 *  - dx (in):               grid spacing
 *   
 * Return value:             none   
 *   
 * NOTES:
 *  - max(phi_*_plus , phi_*_minus) is the value of phi_* that is
 *    used in the time step size calculation.  This may be more
 *    conservative than necessary for Godunov's method, but it is
 *    cheaper to compute.
 *
 */
void LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_CONTROL_VOLUME(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
  const LSMLIB_REAL *phi_x_plus,
  const LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus,
  const LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *cfl_number);
  
  
/*!
 * LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME() computes the 
 * volume integral of the specified function over the region of the 
 * computational domain where the level set function is less than 0.  
 * The computational domain contains only those cells that are included
 * by the control volume data.
 *    
 * Arguments:
 *  - int_F (out):           value of integral of F over the region where 
 *                           \f$ \phi < 0 \f$
 *  - F (in):                function to be integrated
 *  - phi (in):              level set function
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx (in):               grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 */
void LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const int *jlo_F_gb, 
  const int *jhi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME() computes 
 * the volume integral of the specified function over the region of
 * computational domain where the level set function is greater than 0.  
 * The computational domain contains only those cells that are included
 * by the control volume data.
 *
 * Arguments:
 *  - int_F (out):           value of integral of F over the region where 
 *                           \f$ \phi > 0 \f$
 *  - F (in):                function to be integrated
 *  - phi (in):              level set function
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx (in):               grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 */
void LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const int *jlo_F_gb, 
  const int *jhi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM2D_SURFACE_INTEGRAL_CONTROL_VOLUME() computes the surface integral 
 * of the specified function over the region of the computational domain
 * where the level set function equals 0.  The computational domain 
 * contains only those cells that are included by the control volume data.
 *     
 * Arguments:
 *  - int_F (out):           value of integral of F over the region where 
 *                           \f$ \phi = 0 \f$
 *  - F (in):                function to be integrated
 *  - phi (in):              level set function
 *  - phi_* (in):            components of \f$ \nabla \phi \f$
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx (in):               grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           delta-function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 * 
 * Return value:             none
 */
void LSM2D_SURFACE_INTEGRAL_CONTROL_VOLUME(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const int *jlo_F_gb, 
  const int *jhi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb, 
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *epsilon);
  
/*!
*
*  LSM2D_VOXEL_COUNT_GREATER_THAN_ZERO_CONTROL_VOLUME() computes number of voxels whose phi
*  value is less than zero in the given control volume
*
*  Arguments:
*    count (out):           number of voxels where phi < 0
*    phi (in):              scalar field
*    *_gb (in):             index range for ghostbox
*    *_ib (in):             index range for box to include in
*                           calculation
*
*/  
void LSM2D_VOXEL_COUNT_GREATER_THAN_ZERO_CONTROL_VOLUME(
  int *count,
  const LSMLIB_REAL *phi,
  const int *ilo_gb,
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb);
  
/*!
*
*  LSM2D_VOXEL_COUNT_LESS_THAN_ZERO_CONTROL_VOLUME() computes number of voxels whose phi
*  value is less than zero in the given control volume
*
*  Arguments:
*    count (out):           number of voxels where phi < 0
*    phi (in):              scalar field
*    *_gb (in):             index range for ghostbox
*    *_ib (in):             index range for box to include in
*                           calculation
*
*/  
void LSM2D_VOXEL_COUNT_LESS_THAN_ZERO_CONTROL_VOLUME(
  int *count,
  const LSMLIB_REAL *phi,
  const int *ilo_gb,
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb);
  
/*!
 * LSM2D_SURFACE_INTEGRAL_DELTA_CONTROL_VOLUME() computes the surface integral 
 * of the specified function over the region of the computational domain
 * where the level set function equals 0, assuming delta function has been
 * precomputed.  The computational domain contains only those cells that 
 * are included by the control volume data.
 *     
 * Arguments:
 *  - int_F (out):           value of integral of F over the region where 
 *                           \f$ \phi = 0 \f$
 *  - F (in):                function to be integrated
 *  - delta_phi (in):        delta function for the zero level set
 *  - grad_phi_mag (in):     gradient magnitude
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - d* (in):               grid spacing
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 * 
 * Return value:             none
 */  
void LSM2D_SURFACE_INTEGRAL_DELTA_CONTROL_VOLUME(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const int *jlo_F_gb, 
  const int *jhi_F_gb,
  const LSMLIB_REAL *delta_phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb, 
  const LSMLIB_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);
  
#ifdef __cplusplus
}
#endif

#endif
