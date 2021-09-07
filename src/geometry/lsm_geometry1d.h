/*
 * File:        lsm_geometry1d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for 1D Fortran 77 level set method geometry
 *              subroutines
 */

#ifndef INCLUDED_LSM_GEOMETRY_1D_H
#define INCLUDED_LSM_GEOMETRY_1D_H

#include "lsmlib_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_geometry1d.h
 *
 * \brief
 * @ref lsm_geometry1d.h provides support for computing various geometric
 * quantities in one space dimension.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                        name in
 *      C/C++ code                     Fortran code
 *      ----------                     ------------
 */
#define LSM1D_COMPUTE_UNIT_NORMAL      lsm1dcomputeunitnormal_
#define LSM1D_COMPUTE_SIGNED_UNIT_NORMAL                                     \
                                       lsm1dcomputesignedunitnormal_

#define LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO                               \
                                       lsm1dlengthregionphilessthanzero_
#define LSM1D_LENGTH_REGION_PHI_GREATER_THAN_ZERO                            \
                                       lsm1dlengthregionphigreaterthanzero_
#define LSM1D_SIZE_ZERO_LEVEL_SET      lsm1dsizezerolevelset_
					    
#define LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME                \
                              lsm1dlengthregionphilessthanzerocontrolvolume_
#define LSM1D_LENGTH_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME             \
                              lsm1dlengthregionphigreaterthanzerocontrolvolume_
#define LSM1D_SIZE_ZERO_LEVEL_SET_CONTROL_VOLUME                             \
                              lsm1dsizezerolevelsetcontrolvolume_


/*!
 * LSM1D_COMPUTE_UNIT_NORMAL() computes the unit normal vector to the
 * interface from \f$ \nabla \phi \f$.
 *
 * Arguments:
 *  - normal (out):  unit normal vector
 *  - phi_* (in):    components of \f$ \nabla \phi \f$
 *  - *_gb (in):     index range for ghostbox
 *  - *_fb (in):     index range for fillbox
 *
 * Return value:     none
 *
 * NOTES:
 * - When \f$ | \nabla \phi | \f$ is close to zero, the unit normal is
 *   arbitrarily set to be 1.0.
 *
 */
void LSM1D_COMPUTE_UNIT_NORMAL(
  LSMLIB_REAL *normal,
  const int *ilo_normal_gb, 
  const int *ihi_normal_gb,
  const LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *ilo_fb, 
  const int *ihi_fb);


/*!
 * LSM1D_COMPUTE_SIGNED_UNIT_NORMAL() computes the signed unit normal
 * vector (sgn(phi)*normal) to the interface from \f$ \nabla \phi \f$ 
 * using the following smoothed sgn function 
 *   
 * \f[
 *
 *   sgn(\phi) = \phi / \sqrt{ \phi^2 + |\nabla \phi|^2 * dx^2 }
 *   
 * \f]
 *
 * Arguments:
 *  - normal_* (out):     components of unit normal vector
 *  - phi_* (in):         components of \f$ \nabla \phi \f$
 *  - phi (in):           level set function
 *  - dx (in):            grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *
 * Return value:          none
 *
 * NOTES:
 * - When \f$ | \nabla \phi | \f$ is close to zero, the unit normal is
 *   arbitrarily set to be 1.0.
 *
 */
void LSM1D_COMPUTE_SIGNED_UNIT_NORMAL(
  LSMLIB_REAL *normal,
  const int *ilo_normal_gb, 
  const int *ihi_normal_gb,
  const LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const LSMLIB_REAL *dx);

/*!
 * LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO() computes the length of the
 * region where the level set function is less than 0.
 *   
 * Arguments:
 *  - length (out):          length of the region where \f$ \phi < 0 \f$
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
void LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO(
  LSMLIB_REAL *length,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM1D_LENGTH_REGION_PHI_GREATER_THAN_ZERO() computes the length of the
 * region where the level set function is greater than 0.
 *
 * Arguments:
 *  - length (out):          length of the region where \f$ \phi > 0 \f$
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
void LSM1D_LENGTH_REGION_PHI_GREATER_THAN_ZERO(
  LSMLIB_REAL *volume,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM1D_SIZE_ZERO_LEVEL_SET() computes the size of the zero level set.
 *
 * Arguments:
 *  - size (out):            number of points contained in the zero level
 *                           set
 *  - phi (in):              level set function
 *  - phi_* (in):            components of \f$ \nabla \phi \f$
 *  - dx (in):               grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for Heaviside
 *                           function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM1D_SIZE_ZERO_LEVEL_SET(
  LSMLIB_REAL *size,
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
 * LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME() computes the 
 * length of the region of the computational domain where the level set 
 * function is less than 0.  The computational domain contains only 
 * those cells that are included by the control volume data.
 *   
 * Arguments:
 *  - length (out):          length of the region where \f$ \phi < 0 \f$
 *  - phi (in):              level set function
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral calculation)
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
void LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
  LSMLIB_REAL *length,
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
 * LSM1D_LENGTH_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME() computes the 
 * length of the region of the computational domain where the level set 
 * function is greater than 0.  The computational domain contains only 
 * those cells that are included by the control volume data.
 *
 * Arguments:
 *  - length (out):          length of the region where \f$ \phi > 0 \f$
 *  - phi (in):              level set function
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral calculation)
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
void LSM1D_LENGTH_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME(
  LSMLIB_REAL *volume,
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
 * LSM1D_SIZE_ZERO_LEVEL_SET_CONTROL_VOLUME() computes the size of the zero 
 * level set within the computational domain.  The computational domain 
 * contains only those cells that are included by the control volume data.
 *
 * Arguments:
 *  - size (out):            number of points contained in the zero level
 *                           set
 *  - phi (in):              level set function
 *  - phi_* (in):            components of \f$ \nabla \phi \f$
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx (in):               grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for Heaviside
 *                           function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM1D_SIZE_ZERO_LEVEL_SET_CONTROL_VOLUME(
  LSMLIB_REAL *size,
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
