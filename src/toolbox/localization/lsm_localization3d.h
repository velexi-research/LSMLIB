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
 * localization in 3D - building the narrow band and the cut off function following
 * Peng at al. " A PDE-based fast local level set method", 1999.
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
 #define LSM3D_IMPOSE_MAX_PHI_MASK_LOCAL       lsm3dimposemaxphimasklocal_
 #define LSM3D_CHECK_OUTER_NARROW_BAND_LAYER   lsm3dcheckouternarrowbandlayer_ 
 
#ifdef __cplusplus
extern "C" {
#endif
 
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
 
 void LSM3D_DETERMINE_NARROW_BAND_FROM_MASK(
 const LSMLIB_REAL *phi,
 const LSMLIB_REAL *mask,
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
 int *nlo_index, 
 int *nhi_index, 
 int *n_lo,
 int *n_hi,
 const int *level,
 const int *use_mask_sign);
 
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
 
 void  LSM3D_IMPOSE_MAX_PHI_MASK_LOCAL(
 LSMLIB_REAL *phi,
 const LSMLIB_REAL *mask,
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
 const int *nhi_index);
 
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
 
#ifdef __cplusplus
}
#endif
 
 #endif
