/*
 * File:        lsm_geometry2d_local.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.11 $
 * Modified:    $Date: 2006/09/30 18:17:22 $
 * Description: Header file for 2D Fortran 77 level set method geometry
 *              subroutines
 */

#ifndef INCLUDED_LSM_GEOMETRY_2D_LOCAL_H
#define INCLUDED_LSM_GEOMETRY_2D_LOCAL_H

#include "LSMLIB_Config.h"

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
#define LSM2D_COMPUTE_UNIT_NORMAL_LOCAL     lsm2dcomputeunitnormallocal_
#define LSM2D_COMPUTE_SIGNED_UNIT_NORMAL_LOCAL                             \
                                       lsm2dcomputesignedunitnormallocal_

/*!
 * LSM2D_COMPUTE_UNIT_NORMAL_LOCAL() computes the unit normal vector to the
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
 *  - normal (out):      unit normal vector
 *  - phi_* (in):        components of \f$ \nabla \phi \f$
 *  - dx (in):           grid spacing
 *  - *_gb (in):         index range for ghostbox
 *  - index_[xy]  (in):  [xy] coordinates of local (narrow band) points
 *  - n*_index    (in):  index range of points to loop over in index_*
 *  - narrow_band(in):   array that marks voxels outside desired fillbox
 *  - mark_fb(in):       upper limit narrow band value for voxels in 
 *                       fillbox
 *
 * Return value:         none
 *
 */
void LSM2D_COMPUTE_UNIT_NORMAL_LOCAL(
  LSMLIB_REAL *normal_x,
  LSMLIB_REAL *normal_y,
  const int *ilo_normal_gb, 
  const int *ihi_normal_gb,
  const int *jlo_normal_gb, 
  const int *jhi_normal_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
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
 * LSM2D_COMPUTE_SIGNED_UNIT_NORMAL_LOCAL() computes the signed unit normal
 * vector to the interface from \f$ \nabla \phi \f$ using the following 
 * smoothed sgn function 
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
 *  - normal_* (out):    components of unit normal vector
 *  - phi_* (in):        components of \f$ \nabla \phi \f$
 *  - phi (in):          level set function
 *  - dx, dy (in):       grid spacing
 *  - *_gb (in):         index range for ghostbox
 *  - index_[xy]  (in):  [xy] coordinates of local (narrow band) points
 *  - n*_index    (in):  index range of points to loop over in index_*
 *  - narrow_band(in):   array that marks voxels outside desired fillbox
 *  - mark_fb(in):       upper limit narrow band value for voxels in 
 *                       fillbox
 *
 * Return value:         none
 *
 */
void LSM2D_COMPUTE_SIGNED_UNIT_NORMAL_LOCAL(
  LSMLIB_REAL *normal_x,
  LSMLIB_REAL *normal_y,
  const int *ilo_normal_gb, 
  const int *ihi_normal_gb,
  const int *jlo_normal_gb, 
  const int *jhi_normal_gb,
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
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

#ifdef __cplusplus
}
#endif

#endif
