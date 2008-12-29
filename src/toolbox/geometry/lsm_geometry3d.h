/*
 * File:        lsm_geometry3d.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.16 $
 * Modified:    $Date: 2007/04/20 16:39:27 $
 * Description: Header file for 3D Fortran 77 level set method geometry
 *              subroutines
 */

#ifndef INCLUDED_LSM_GEOMETRY_3D_H
#define INCLUDED_LSM_GEOMETRY_3D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_geometry3d.h
 *
 * \brief
 * @ref lsm_geometry3d.h provides support for computing various geometric
 * quantities in three space dimensions.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                        name in
 *      C/C++ code                     Fortran code
 *      ----------                     ------------
 */
#define LSM3D_COMPUTE_UNIT_NORMAL      lsm3dcomputeunitnormal_
#define LSM3D_COMPUTE_SIGNED_UNIT_NORMAL                                     \
                                       lsm3dcomputesignedunitnormal_
#define LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO                               \
                                       lsm3dvolumeregionphilessthanzero_
#define LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO                            \
                                       lsm3dvolumeregionphigreaterthanzero_
#define LSM3D_SURFACE_AREA_ZERO_LEVEL_SET                                    \
                                       lsm3dsurfaceareazerolevelset_
#define LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_DELTA                                    \
                                       lsm3dsurfaceareazerolevelsetdelta_
				       
#define LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME                \
                          lsm3dvolumeregionphilessthanzerocontrolvolume_
#define LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME             \
                          lsm3dvolumeregionphigreaterthanzerocontrolvolume_
#define LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_CONTROL_VOLUME                     \
                          lsm3dsurfaceareazerolevelsetcontrolvolume_
#define LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_DELTA_CONTROL_VOLUME                     \
                          lsm3dsurfaceareazerolevelsetdeltacontrolvolume_
			  
/*!
 * LSM3D_COMPUTE_UNIT_NORMAL() computes the unit normal vector to the
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
 *   arbitrarily set to be (1.0, 0.0, 0.0).
 *
 */
void LSM3D_COMPUTE_UNIT_NORMAL(
  LSMLIB_REAL *normal_x,
  LSMLIB_REAL *normal_y,
  LSMLIB_REAL *normal_z,
  const int *ilo_normal_gb, 
  const int *ihi_normal_gb,
  const int *jlo_normal_gb, 
  const int *jhi_normal_gb,
  const int *klo_normal_gb, 
  const int *khi_normal_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *phi_z,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb, 
  const int *khi_grad_phi_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb, 
  const int *jhi_fb,
  const int *klo_fb, 
  const int *khi_fb);


/*!
 * LSM3D_COMPUTE_SIGNED_UNIT_NORMAL() computes the signed unit normal
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
 *  - dx, dy, dz (in):    grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *
 * Return value:          none
 *
 * NOTES:
 * - When \f$ | \nabla \phi | \f$ is close to zero, the unit normal is
 *   arbitrarily set to be (1.0, 0.0, 0.0).
 *
 */
void LSM3D_COMPUTE_SIGNED_UNIT_NORMAL(
  LSMLIB_REAL *normal_x,
  LSMLIB_REAL *normal_y,
  LSMLIB_REAL *normal_z,
  const int *ilo_normal_gb, 
  const int *ihi_normal_gb,
  const int *jlo_normal_gb, 
  const int *jhi_normal_gb,
  const int *klo_normal_gb, 
  const int *khi_normal_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *phi_z,
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
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb, 
  const int *jhi_fb,
  const int *klo_fb, 
  const int *khi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz);


/*!
 * LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO() computes the volume of the
 * region where the level set function is less than 0.
 *
 * Arguments:
 *  - volume (out):          volume of the region where \f$ \phi < 0 \f$
 *  - phi (in):              level set function
 *  - dx, dy, dz (in):       grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(
  LSMLIB_REAL *volume,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const LSMLIB_REAL *epsilon);


/*! 
 * LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO() computes the volume of 
 * the region where the level set function is greater than 0.
 *
 * Arguments:
 *  - volume (out):          volume of the region where \f$ \phi > 0 \f$
 *  - phi (in):              level set function
 *  - dx, dy, dz (in):       grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO(
  LSMLIB_REAL *volume,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM3D_SURFACE_AREA_ZERO_LEVEL_SET() computes the surface area of the
 * surface defined by the zero level set.
 *
 * Arguments:
 *  - area (out):            area of the surface defined by the zero level
 *                           set
 *  - phi (in):              level set function
 *  - phi_* (in):            components of \f$ \nabla \phi \f$
 *  - dx, dy, dz (in):       grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:         none
 *
 */
void LSM3D_SURFACE_AREA_ZERO_LEVEL_SET(
  LSMLIB_REAL *surface_area,
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
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb, 
  const int *khi_grad_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const LSMLIB_REAL *epsilon);

/*!
 * LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_DELTA() computes the surface area of the
 * surface defined by the zero level set.
 *
 * Arguments:
 *  - area (out):            area of the surface defined by the zero level
 *                           set
 *  - delta_phi (in):        delta function (precomputed)
 *  - phi_* (in):            components of \f$ \nabla \phi \f$
 *  - dx, dy, dz (in):       grid spacing
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:         none
 *
 */
void LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_DELTA(
  LSMLIB_REAL *surface_area,
  const LSMLIB_REAL *delta_phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *phi_z,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb, 
  const int *khi_grad_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz);

/*!
 * LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME() computes the 
 * volume of the region of the computational domain where the level set 
 * function is less than 0.  The computational domain contains only 
 * those cells that are included by the control volume data.
 *
 * Arguments:
 *  - volume (out):          volume of the region where \f$ \phi < 0 \f$
 *  - phi (in):              level set function
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx, dy, dz (in):       grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
  LSMLIB_REAL *volume,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb, 
  const int *klo_control_vol_gb, 
  const int *khi_control_vol_gb,
  const int *control_vol_sgn, 
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const LSMLIB_REAL *epsilon);


/*! 
 * LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME() computes the 
 * volume of the region of the computational domain where the level set 
 * function is greater than 0.  The computational domain contains only 
 * those cells that are included by the control volume data.
 *
 * Arguments:
 *  - volume (out):          volume of the region where \f$ \phi > 0 \f$
 *  - phi (in):              level set function
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx, dy, dz (in):       grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME(
  LSMLIB_REAL *volume,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb, 
  const int *klo_control_vol_gb, 
  const int *khi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_CONTROL_VOLUME() computes the surface 
 * area of the surface defined by the zero level set within the computational
 * domain.  The computational domain contains only those cells that are 
 * included by the control volume data.
 *
 * Arguments:
 *  - area (out):            area of the surface defined by the zero level
 *                           set
 *  - phi (in):              level set function
 *  - phi_* (in):            components of \f$ \nabla \phi \f$
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx, dy, dz (in):       grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:         none
 *
 */
void LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_CONTROL_VOLUME(
  LSMLIB_REAL *surface_area,
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
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb, 
  const int *khi_grad_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb, 
  const int *klo_control_vol_gb, 
  const int *khi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz,
  const LSMLIB_REAL *epsilon);


/*!
 * LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_DELTA_CONTROL_VOLUME() computes the surface area of the
 * surface defined by the zero level set.
 *
 * Arguments:
 *  - area (out):            area of the surface defined by the zero level
 *                           set
 *  - delta_phi (in):        externally computed delta function
 *  - phi_* (in):            components of \f$ \nabla \phi \f$
 *  - control_vol (in):      control volume data (used to exclude cells
 *                           from the integral calculation)
 *  - control_vol_sgn (in):  1 (-1) if positive (negative) control volume
 *                           points should be used
 *  - dx, dy, dz (in):       grid spacing
 *  - delta_phi (in):        externally computed delta function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:         none
 *
 */
void LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_DELTA_CONTROL_VOLUME(
  LSMLIB_REAL *surface_area,
  const LSMLIB_REAL *delta_phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *phi_z,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb, 
  const int *khi_grad_phi_gb,
  const LSMLIB_REAL *control_vol,
  const int *ilo_control_vol_gb, 
  const int *ihi_control_vol_gb,
  const int *jlo_control_vol_gb, 
  const int *jhi_control_vol_gb, 
  const int *klo_control_vol_gb, 
  const int *khi_control_vol_gb,
  const int *control_vol_sgn,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *dz);

/*!
 * LSM3D_findLineInTetrahedron() determines if there is an intersection
 * between the specified tetrahedron and the line defined by the zero-level
 * sets of \f$ \phi \f$ and \f$ \psi \f$ .  If an intersection exists, then
 * the two endpoints of the \f$ \{ \phi=0,\psi=0 \} \f$ line segment
 * contained by the tetrahedron are returned by the endpt1 and endpt2 
 * arguments.  
 *
 * Arguments:
 *  - endpt1 (out):  coordinate of one end point of \f$ \{ \phi=0,\psi=0 \} \f$
 *                   line in specified tetrahedron
 *  - endpt2 (out):  coordinate of other end point of
 *                   \f$ \{ \phi=0,\psi=0 \} \f$ line in specified tetrahedron
 *  - x1 (in):       coordinate corner #1 of tetrahedron
 *  - x2 (in):       coordinate corner #2 of tetrahedron
 *  - x3 (in):       coordinate corner #3 of tetrahedron
 *  - x4 (in):       coordinate corner #4 of tetrahedron
 *  - phi (in):      value of \f$ \phi \f$ at x1, x2, x3, and x4 (in that order)
 *  - psi (in):      value of \f$ \psi \f$ at x1, x2, x3, and x4 (in that order)
 *
 * Return value:    error code
 *  - nonnegative integer:  function exited with no errors.
 *                          value indicates number of single-point
 *                          intersections with faces of the tetrahedron
 *  - negative integer:     indicates point of failure within the function
 * 
 * NOTES:
 *  - If endpt1 and endpt2 are distinct, then (endpt2 - endpt1) points in 
 *    the direction [ grad(phi) cross grad(psi) ].
 *
 *  - It is the user's responsibility to ensure that there is sufficient
 *    memory allocated for the endpoints of the \f$ \{ \phi=0,\psi=0 \} \f$
 *    line: endpt1 and endpt2.
 *
 *  - It is the user's responsibility to ensure that x1, x2, x3, and x4
 *    are arrays of size at least 3.
 *
 *  - If the \f$ \{ \phi=0,\psi=0 \} \f$ line lies in the plane of one 
 *    of the faces of the tetrahedron, only the intersections of the 
 *    line with with the three other faces of the tetrahedron are 
 *    recorded.  
 *
 */
int LSM3D_findLineInTetrahedron(
  LSMLIB_REAL *endpt1,
  LSMLIB_REAL *endpt2,
  const LSMLIB_REAL *x1,
  const LSMLIB_REAL *x2,
  const LSMLIB_REAL *x3,
  const LSMLIB_REAL *x4,
  const LSMLIB_REAL *phi,
  const LSMLIB_REAL *psi);

#ifdef __cplusplus
}
#endif

#endif
