/*
 * File:        lsm_utilities1d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for 1D Fortran 77 level set method utility 
 *              subroutines
 */

#ifndef INCLUDED_LSM_UTILITIES_1D_H
#define INCLUDED_LSM_UTILITIES_1D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_utilities1d.h
 *
 * \brief 
 * @ref lsm_utilities1d.h provides several utility functions that support
 * level set method calculations in one space dimension.
 *
 */


/*
 * Link between C/C++ and Fortran function names
 *
 *      name in                        name in
 *      C/C++ code                     Fortran code
 *      ----------                     ------------
 */
#define LSM1D_MAX_NORM_DIFF            lsm1dmaxnormdiff_
#define LSM1D_COMPUTE_STABLE_ADVECTION_DT                                     \
                                       lsm1dcomputestableadvectiondt_
#define LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT                                    \
                                       lsm1dcomputestablenormalveldt_
#define LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO                              \
                                       lsm1dvolumeintegralphilessthanzero_
#define LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO                           \
                                       lsm1dvolumeintegralphigreaterthanzero_
#define LSM1D_SURFACE_INTEGRAL         lsm1dsurfaceintegral_
#define LSM1D_MAX_NORM_DIFF_CONTROL_VOLUME                                    \
                             lsm1dmaxnormdiffcontrolvolume_
#define LSM1D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME                      \
                             lsm1dcomputestableadvectiondtcontrolvolume_
#define LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME                     \
                             lsm1dcomputestablenormalveldtcontrolvolume_
#define LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME               \
                             lsm1dvolumeintegralphilessthanzerocontrolvolume_
#define LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME            \
                             lsm1dvolumeintegralphigreaterthanzerocontrolvolume_
#define LSM1D_SURFACE_INTEGRAL_CONTROL_VOLUME                                 \
                             lsm1dsurfaceintegralcontrolvolume_

/*!
 * LSM1D_MAX_NORM_DIFF() computes the max norm of the difference
 * between the two specified scalar fields.
 *      
 * Arguments:
 *  - max_norm_diff (out):   max norm of the difference between the fields
 *  - field1 (in):           scalar field 1
 *  - field2 (in):           scalar field 2
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for box to include in norm
 *                           calculation
 *
 * Return value:             none
 * 
 */
void LSM1D_MAX_NORM_DIFF(
  LSMLIB_REAL *max_norm_diff,
  const LSMLIB_REAL *field1,
  const int *ilo_field1_gb, 
  const int *ihi_field1_gb,
  const LSMLIB_REAL *field2,
  const int *ilo_field2_gb, 
  const int *ihi_field2_gb,
  const int *ilo_ib, 
  const int *ihi_ib);


/*!
 * LSM1D_COMPUTE_STABLE_ADVECTION_DT() computes the stable time step size
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
 * 
 */
void LSM1D_COMPUTE_STABLE_ADVECTION_DT(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_x,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *cfl_number);


/*!
 * LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT() computes the stable time step
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
void LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const LSMLIB_REAL *phi_x_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *cfl_number);


/*!
 * LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO() computes the volume integral of
 *  the specified function over the region where the level set function
 *  is less than 0.  
 *    
 * Arguments:
 *  - int_F (out):           value of integral of F over the region where 
 *                           \f$ \phi < 0 \f$
 *  - F (in):                function to be integrated
 *  - phi (in):              level set function
 *  - dx (in):               grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 * 
 */
void LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO() computes the volume integral
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
 *
 */
void LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM1D_SURFACE_INTEGRAL() computes the surface integral of the specified
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
 * 
 */
void LSM1D_SURFACE_INTEGRAL(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *epsilon);

/*!
 * LSM1D_MAX_NORM_DIFF_CONTROL_VOLUME() computes the max norm of the 
 * difference between the two specified scalar fields in the region
 * of the computational domain included by the control volume data.
 *      
 * Arguments:
 *  - max_norm_diff (out):   max norm of the difference between the fields
 *  - field1 (in):           scalar field 1
 *  - field2 (in):           scalar field 2
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the max norm calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for box to include in norm
 *                           calculation
 *
 * Return value:             none
 * 
 */
void LSM1D_MAX_NORM_DIFF_CONTROL_VOLUME(
  LSMLIB_REAL *max_norm_diff,
  const LSMLIB_REAL *field1,
  const int *ilo_field1_gb, 
  const int *ihi_field1_gb,
  const LSMLIB_REAL *field2,
  const int *ilo_field2_gb, 
  const int *ihi_field2_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib);


/*!
 * LSM1D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME() computes the stable 
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
 * 
 */
void LSM1D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_x,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *cfl_number);


/*!
 * LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME() computes the 
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
void LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const LSMLIB_REAL *phi_x_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *cfl_number);


/*!
 * LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME() computes the 
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
 * 
 */
void LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME() computes 
 * the volume integral of the specified function over the region of the
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
 *
 */
void LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM1D_SURFACE_INTEGRAL_CONTROL_VOLUME() computes the surface integral 
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
 * 
 */
void LSM1D_SURFACE_INTEGRAL_CONTROL_VOLUME(
  LSMLIB_REAL *int_F,
  const LSMLIB_REAL *F,
  const int *ilo_F_gb, 
  const int *ihi_F_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *epsilon);




#ifdef __cplusplus
}
#endif

#endif
