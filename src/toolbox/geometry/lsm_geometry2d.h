/*
 * File:        lsm_geometry2d.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.13 $
 * Modified:    $Date: 2006/11/01 00:25:18 $
 * Description: Header file for 2D Fortran 77 level set method geometry
 *              subroutines
 */

#ifndef INCLUDED_LSM_GEOMETRY_2D_H
#define INCLUDED_LSM_GEOMETRY_2D_H

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_geometry2d.h
 *
 * \brief
 * @ref lsm_geometry2d.h provides support for computing various geometric
 * quantities in two space dimensions.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                        name in
 *      C/C++ code                     Fortran code
 *      ----------                     ------------
 */
#define LSM2D_COMPUTE_UNIT_NORMAL     lsm2dcomputeunitnormal_
#define LSM2D_COMPUTE_SIGNED_UNIT_NORMAL                                     \
                                       lsm2dcomputesignedunitnormal_

#define LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO                               \
                                       lsm2darearegionphilessthanzero_
#define LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO                            \
                                       lsm2darearegionphigreaterthanzero_
#define LSM2D_PERIMETER_ZERO_LEVEL_SET                                     \
                                       lsm2dperimeterzerolevelset_

#define LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME                \
                           lsm2darearegionphilessthanzerocontrolvolume_
#define LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME             \
                           lsm2darearegionphigreaterthanzerocontrolvolume_
#define LSM2D_PERIMETER_ZERO_LEVEL_SET_CONTROL_VOLUME                      \
                           lsm2dperimeterzerolevelsetcontrolvolume_

/*!
 * LSM2D_COMPUTE_UNIT_NORMAL() computes the unit normal vector to the
 * interface from \f$ \nabla \phi \f$ using the slightly modified 
 * expression for the norm:
 *
 * \f[
 *
 *   norm = \sqrt{ |\nabla \phi|^2 + dx^2 }
 *
 * \f]
 *
 * This expression avoids division by zero in computing the unit
 * normal vector.
 *
 * Arguments:
 *  - normal (out):  unit normal vector
 *  - phi_* (in):    components of \f$ \nabla \phi \f$
 *  - dx (in):       grid spacing
 *  - *_gb (in):     index range for ghostbox
 *  - *_fb (in):     index range for fillbox
 *
 * Return value:     none
 *
 */
void LSM2D_COMPUTE_UNIT_NORMAL(
  double *normal_x,
  double *normal_y,
  const int *ilo_normal_gb, 
  const int *ihi_normal_gb,
  const int *jlo_normal_gb, 
  const int *jhi_normal_gb,
  const double *phi_x,
  const double *phi_y,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb, 
  const int *jhi_fb,
  const double *dx,
  const double *dy);


/*!
 * LSM2D_COMPUTE_SIGNED_UNIT_NORMAL() computes the signed unit normal
 * vector (sgn(phi)*normal) to the interface from \f$ \nabla \phi \f$ 
 * using the following smoothed sgn function 
 *   
 * \f[
 *
 *   sgn(\phi) = \phi / \sqrt{ \phi^2 + |\nabla \phi|^2 * dx^2 }
 *   
 * \f]
 *
 * and the following slightly modified expression for the norm
 *   
 * \f[
 *
 *   norm = \sqrt{ |\nabla \phi|^2 + dx^2}
 *    
 * \f]
 *
 * This expression avoids division by zero in computing the unit
 * normal vector.
 *   
 * Arguments:
 *  - normal_* (out):     components of unit normal vector
 *  - phi_* (in):         components of \f$ \nabla \phi \f$
 *  - phi (in):           level set function
 *  - dx, dy (in):        grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *
 * Return value:          none
 *
 */
void LSM2D_COMPUTE_SIGNED_UNIT_NORMAL(
  double *normal_x,
  double *normal_y,
  const int *ilo_normal_gb, 
  const int *ihi_normal_gb,
  const int *jlo_normal_gb, 
  const int *jhi_normal_gb,
  const double *phi_x,
  const double *phi_y,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const double *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb, 
  const int *jhi_fb,
  const double *dx,
  const double *dy);

/*! 
 * LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO() computes the area of the
 * region where the level set function is less than 0.
 *
 * Arguments:
 *  - area (out):            area of the region where \f$ \phi < 0 \f$
 *  - phi (in):              level set function
 *  - dx, dy (in):           grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO(
  double *area,
  const double *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const double *dx,
  const double *dy,
  const double *epsilon);


/*!
 * LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO() computes the area of the
 * region where the level set function is greater than 0.
 *
 * Arguments:
 *  - area (out):            area of the region where \f$ \phi > 0 \f$
 *  - phi (in):              level set function
 *  - dx, dy (in):           grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO(
  double *area,
  const double *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const double *dx,
  const double *dy,
  const double *epsilon);


/*! 
 * LSM2D_PERIMETER_ZERO_LEVEL_SET() computes the perimeter of the zero
 *  level set.
 *
 * Arguments:
 *  - perimeter (out):       perimeter of curve defined by the zero level
 *                           set
 *  - phi (in):              level set function
 *  - phi_* (in):            components of \f$ \nabla \phi \f$
 *  - dx, dy (in):           grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox 
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM2D_PERIMETER_ZERO_LEVEL_SET(
  double *perimeter,
  const double *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const double *phi_x,
  const double *phi_y,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const double *dx,
  const double *dy,
  const double *epsilon);



/*! 
 * LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME() computes the 
 * area of the region of the computational domain where the level set 
 * function is less than 0.  The computational domain contains only 
 * those cells that are included by the control volume data.
 *
 * Arguments:
 *  - area (out):            area of the region where \f$ \phi < 0 \f$
 *  - phi (in):              level set function
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx, dy (in):           grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
  double *area,
  const double *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const double *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn, 
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const double *dx,
  const double *dy,
  const double *epsilon);


/*!
 * LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME() computes the 
 * area of the region of the computational domain where the level set 
 * function is greater than 0.  The computational domain contains only 
 * those cells that are included by the control volume data.
 *
 * Arguments:
 *  - area (out):            area of the region where \f$ \phi > 0 \f$
 *  - phi (in):              level set function
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx, dy (in):           grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME(
  double *area,
  const double *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const double *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb,
  const int *control_vol_sgn, 
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const double *dx,
  const double *dy,
  const double *epsilon);


/*! 
 * LSM2D_PERIMETER_ZERO_LEVEL_SET_CONTROL_VOLUME() computes the perimeter 
 * of the zero level set within the computational domain.  The computational 
 * domain contains only those cells that are included by the control volume 
 * data.
 *
 * Arguments:
 *  - perimeter (out):       perimeter of curve defined by the zero level
 *                           set
 *  - phi (in):              level set function
 *  - phi_* (in):            components of \f$ \nabla \phi \f$
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx, dy (in):           grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox 
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM2D_PERIMETER_ZERO_LEVEL_SET_CONTROL_VOLUME(
  double *perimeter,
  const double *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const double *phi_x,
  const double *phi_y,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const double *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb, 
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const double *dx,
  const double *dy,
  const double *epsilon);

#ifdef __cplusplus
}
#endif

#endif
