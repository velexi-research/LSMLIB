/*
 * File:        lsm_localization3d.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.27 $
 * Modified:    $Date: 2006/10/05 15:03:45 $
 * Description: Header file for 3D Fortran 77 narrow-band level set 
 *              subroutines
 */

#ifndef INCLUDED_LSM_LOCALIZATION_3D_H
#define INCLUDED_LSM_LOCALIZATION_3D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_localization3d.h
 *
 * \brief 
 * @ref lsm_localization3d.h provides several functions that support
 * localization in 3D - building the narrow band and the cut off function 
 * following Peng at al. " A PDE-based fast local level set method", 1999.
 *
 */

/* Link between C/C++ and Fortran function names
 *
 *      name in                             name in
 *      C/C++ code                          Fortran code
 *      ----------                          ------------
 */
 
 #define LSM3D_DETERMINE_NARROW_BAND           lsm3ddeterminenarrowband_
 #define LSM3D_MARK_NARROW_BAND_BOUNDARY_LAYER lsm3dmarknarrowbandboundarylayer_
 #define LSM3D_DETERMINE_NARROW_BAND_FROM_MASK lsm3ddeterminenarrowbandfrommask_
 #define LSM3D_MULTIPLY_CUT_OFF_LSE_RHS_LOCAL  lsm3dmultiplycutofflserhslocal_
 #define LSM3D_CHECK_OUTER_NARROW_BAND_LAYER   lsm3dcheckouternarrowbandlayer_ 
 
 #define LSM3D_IMPOSE_MASK_LOCAL               lsm3dimposemasklocal_
 #define LSM3D_COPY_DATA_LOCAL                 lsm3dcopydatalocal_
 
#ifdef __cplusplus
extern "C" {
#endif
 
/*!
*
*  LSM3D_DETERMINE_NARROW_BAND() finds the narrow band voxels around the zero
*  level set of the specified width. Narrow band neighbors (up to the desired
*  level) are marked as well.
*  Narrow band of level 0 - actual narrow band voxels
*  Narrow band of level L - voxels +/-L voxels in each coordinate direction
*                         (needed for correct computation of derivatives
*                         at the actual narrow band voxels) 
*
*  Arguments:
*    phi(in):          level set functions (assumed signed distance function)
*    narrow_band(out): array with values L+1 for narrow band level L voxels
*                      and 0 otherwise
*    index_*(out):     array with coordinates of narrow band voxels
*                      indices of level L narrow band stored consecutively
*    n*_index(in):     (allocated) index range of index_* arrays 
*    n_lo(out):        array, n_lo[L] is starting index of the level L narrow
*                      band voxels
*    n_hi(out):        array, n_hi[L] is ending index of the level L narrow
*                      band voxels  
*    level(in):        number of narrow band levels to mark
*    width(in):        narrow band width (distance to the zero level set)
*    width_inner(in):  inner narrow band width
*    index_outer(out): indices of the narrow band voxels such that 
*                      width_inner <= abs(phi) < width
*    n*_plus(out):     index range of 'index_outer'  elements for which
*                      phi values satisfy  0 < width_inner <= phi < width  
*    n*_minus(out):    index range of 'index_outer'  elements for which
*                      phi values satisfy  0> -width_inner >= phi > -width 
*    *_gb (in):        index range for ghostbox
*
*    Notes:
*    - narrow_band, index_*, index_outer arrays assumed allocated beforehand
*    - index_outer stores indices of points in narrow band with positive and
*     negative phi values separately in order to be able to identify change
*     of signs for phi (i.e. when zero level set crosses into the outer layer),
*     see lsm3dCheckOuterNarrowBandLayer() 
*/ 
 void LSM3D_DETERMINE_NARROW_BAND(
 const LSMLIB_REAL *phi,
 const int *ilo_gb, 
 const int *ihi_gb,
 const int *jlo_gb, 
 const int *jhi_gb,
 const int *klo_gb, 
 const int *khi_gb,
 unsigned char *narrow_band,
 const int *ilo_nb_gb, 
 const int *ihi_nb_gb,
 const int *jlo_nb_gb, 
 const int *jhi_nb_gb,
 const int *klo_nb_gb, 
 const int *khi_nb_gb,
 int *index_x,
 int *index_y, 
 int *index_z,
 const int *nlo_index, 
 const int *nhi_index,
 int *n_lo,
 int *n_hi,
 int  *index_outer,
 const int *nlo_index_outer, 
 const int *nhi_index_outer,
 int *nlo_index_outer_plus, 
 int *nhi_index_outer_plus,
 int *nlo_index_outer_minus, 
 int *nhi_index_outer_minus,
 const LSMLIB_REAL *width,
 const LSMLIB_REAL *width_inner,
 const int *level);
 
 
/*!
* LSM3D_MARK_NARROW_BAND_BOUNDARY_LAYER() marks planes x = ilo_fb, x = ihi_fb,
* y = jlo_fb, y = jhi_fb, z = klo_fb and z = khi_fb in narrow band array 
* with the specified mark. Potentially used for identification of voxel
* layers near the volume boundary.
*
*  Arguments:
*    narrow_band(out): array with values L+1 for narrow band level L voxels
*                      and 0 otherwise
*    *_gb (in):        index range for ghostbox
*    *_fb (in):        index range for fillbox
*    mark_boundary_layer  distinctive mark (>=120 or so) for boundary layer
*                      i.e. ghost box voxels that do not belong to fill box
*
*/
 void LSM3D_MARK_NARROW_BAND_BOUNDARY_LAYER(
 unsigned char *narrow_band,
 const int *ilo_gb, 
 const int *ihi_gb,
 const int *jlo_gb, 
 const int *jhi_gb,
 const int *klo_gb, 
 const int *khi_gb,
 const int *ilo_fb, 
 const int *ihi_fb,
 const int *jlo_fb, 
 const int *jhi_fb,
 const int *klo_fb, 
 const int *khi_fb,
 const unsigned char *mark_boundary_layer);

 
/*!
*
*  LSM3D_DETERMINE_NARROW_BAND_FROM_MASK() initializes the narrow band voxels 
*  as either positive of negative phase of array mask. Narrow band neighbors 
*  (up to the desired level) are marked as well.
*  Narrow band of level 0 - actual narrow band voxels
*  Narrow band of level L - voxels +/-L voxels in each coordinate direction
*                         (needed for correct computation of derivatives
*                         at the actual narrow band voxels) 
*
*  Arguments:
*    mask(in):         level set function one phase of which will determine
*                      narrow band
*    narrow_band(out): array with values L+1 for narrow band level L voxels
*                      and 0 otherwise
*    index_*(out):     array with coordinates of narrow band voxels
*                      indices of level L narrow band stored consecutively
*    n*_index:         index range of index_* arrays
*    n_lo(out):        array, n_lo[L] is starting index of the level L narrow
*                      band voxels
*    n_hi(out):        array, n_hi[L] is ending index of the level L narrow
*                      band voxels
*    level(in):        number of narrow band levels to mark
*    width(in):        narrow band width (distance to the zero level set)
*    *_gb (in):        index range for ghostbox
*    use_mask_sign     1 or -1 depending if positive or negative phase of 
*                      mask to be used to initialize the narrow band
*
*    Notes:
*    - narrow_band and index_* arrays assumed allocated beforehand
*
*/
 void LSM3D_DETERMINE_NARROW_BAND_FROM_MASK(
 const LSMLIB_REAL *mask,
 const int *ilo_mask_gb, 
 const int *ihi_mask_gb,
 const int *jlo_mask_gb, 
 const int *jhi_mask_gb,
 const int *klo_mask_gb, 
 const int *khi_mask_gb, 
 unsigned char *narrow_band,
 const int *ilo_nb_gb, 
 const int *ihi_nb_gb,
 const int *jlo_nb_gb, 
 const int *jhi_nb_gb,
 const int *klo_nb_gb, 
 const int *khi_nb_gb,
 int *index_x,
 int *index_y, 
 int *index_z,
 int *nlo_index, 
 int *nhi_index, 
 int *n_lo,
 int *n_hi,
 const int *level,
 const int *use_mask_sign);
 
/*!
*
*  LSM3D_MULTIPLY_CUT_OFF_LSE_RHS_LOCAL() multiplies the right hand side of
*  the level set method equation with the cut off function described in
*  Peng at al. '99, "A PDE-Based Fast Local Level Set Method"
*
*  Arguments:
*    phi(in):           level set method function
*    lse_rhs (in/out):  right-hand of level set equation         
*    index_*(out):      array with coordinates of narrow band voxels
*                       indices
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):       upper limit narrow band value for voxels in 
*                       fillbox
*    n*_index:          index range of index_* arrays
*    *_gb (in):         index range for ghostbox
*
*/
 void LSM3D_MULTIPLY_CUT_OFF_LSE_RHS_LOCAL(
 const LSMLIB_REAL *phi,
 LSMLIB_REAL *lse_rhs,
 const int *ilo_gb, 
 const int *ihi_gb,
 const int *jlo_gb, 
 const int *jhi_gb,
 const int *klo_gb, 
 const int *khi_gb,
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
 const unsigned char *mark_fb,
 const LSMLIB_REAL *beta,
 const LSMLIB_REAL *gamma);

/*!
*
*  LSM3D_CHECK_OUTER_NARROW_BAND_LAYER() checks outer narrow band voxels for 
*  changes in the sign of level set function phi values (indicates that the 
*  interface represented by the zero level set has crossed over into the 
*  outer layer).
*
*  Arguments:
*    phi(in):          level set functions (assumed signed distance function)
*    index_*(in):      array with coordinates of narrow band voxels
*                      indices of level L narrow band stored consecutively
*    n*_index(in):     index range of points to loop over in index_* arrays 
*    *_gb (in):        index range for ghostbox
*
*/
 void LSM3D_CHECK_OUTER_NARROW_BAND_LAYER(
 int  *change_sign,
 const LSMLIB_REAL *phi,
 const int *ilo_gb, 
 const int *ihi_gb,
 const int *jlo_gb, 
 const int *jhi_gb,
 const int *klo_gb, 
 const int *khi_gb,
 int *index_x,
 int *index_y, 
 int *index_z,
 int *nlo_index, 
 int *nhi_index,
 const int  *index_outer,
 const int *nlo_index_outer, 
 const int *nhi_index_outer,
 const int *nlo_index_outer_plus, 
 const int *nhi_index_outer_plus,
 const int *nlo_index_outer_minus, 
 const int *nhi_index_outer_minus);
 

/*!
*
*  LSM3D_IMPOSE_MASK_LOCAL() replaces phi with the maximum of
*  phi and mask only at the local (narrow band) points.
*
*  Arguments:
*    phi(in/out):     level set function
*    mask(in):        masking level set function
*    index_*(in):     array(s) with coordinates of narrow band voxels
*                     indices of level L narrow band stored consecutively
*    n*_index:        index range of index_* arrays 
*    *_gb (in):       index range for ghostbox
*
*/
 void  LSM3D_IMPOSE_MASK_LOCAL(
 LSMLIB_REAL *dest,
 const int *ilo_dest_gb, 
 const int *ihi_dest_gb,
 const int *jlo_dest_gb, 
 const int *jhi_dest_gb,
 const int *klo_dest_gb, 
 const int *khi_dest_gb,
 LSMLIB_REAL *src,
 const int *ilo_src_gb, 
 const int *ihi_src_gb,
 const int *jlo_src_gb, 
 const int *jhi_src_gb,
 const int *klo_src_gb, 
 const int *khi_src_gb,
 LSMLIB_REAL *mask,
 const int *ilo_mask_gb, 
 const int *ihi_mask_gb,
 const int *jlo_mask_gb, 
 const int *jhi_mask_gb,
 const int *klo_mask_gb, 
 const int *khi_mask_gb,
 const int *index_x,
 const int *index_y, 
 const int *index_z,
 const int *nlo_index, 
 const int *nhi_index);

/*!
*
*  LSM3D_COPY_DATA_LOCAL() copies array data from source to destination 
*  only from locations specified by narrow band.
*  Arguments:
*    dest (out):         destination scalar field
*    src (in):           source scalar field
*    *_gb (in):          index range for ghostbox
*    index_*(in):        array with coordinates of narrow band voxels
*    n*_index(in):       index range of index_* arrays
*
*/
 void  LSM3D_COPY_DATA_LOCAL(
  const LSMLIB_REAL *dest,
  const int *ilo_dest_gb, 
  const int *ihi_dest_gb,
  const int *jlo_dest_gb, 
  const int *jhi_dest_gb,
  const int *klo_dest_gb, 
  const int *khi_dest_gb,
  const LSMLIB_REAL *src,
  const int *ilo_src_gb, 
  const int *ihi_src_gb,
  const int *jlo_src_gb, 
  const int *jhi_src_gb,
  const int *klo_src_gb, 
  const int *khi_src_gb, 
  const int *index_x,
  const int *index_y,
  const int *index_z,
  const int *nlo_index,
  const int *nhi_index);  
 
#ifdef __cplusplus
}
#endif
 
 #endif
