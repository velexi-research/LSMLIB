/*
 * File:        lsm_localization2d.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.27 $
 * Modified:    $Date: 2006/10/05 15:03:45 $
 * Description: Header file for 2D Fortran 77 narrow-band level set 
 *              subroutines
 */

#ifndef INCLUDED_LSM_LOCALIZATION_2D_H
#define INCLUDED_LSM_LOCALIZATION_2D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_localization2d.h
 *
 * \brief 
 * @ref lsm_localization2d.h provides several functions that support
 * localization in 2D - building the narrow band and the cut off function following
 * Peng at al. " A PDE-based fast local level set method", 1999.
 *
 */

/* Link between C/C++ and Fortran function names
 *
 *      name in                             name in
 *      C/C++ code                          Fortran code
 *      ----------                          ------------
 */
 
 #define LSM2D_DETERMINE_NARROW_BAND           lsm2ddeterminenarrowband_
 #define LSM2D_MARK_NARROW_BAND_BOUNDARY_LAYER lsm2dmarknarrowbandboundarylayer_
 #define LSM2D_DETERMINE_NARROW_BAND_FROM_MASK lsm2ddeterminenarrowbandfrommask_
 #define LSM2D_MULTIPLY_CUT_OFF_LSE_RHS_LOCAL  lsm2dmultiplycutofflserhslocal_
 #define LSM2D_IMPOSE_MAX_PHI_MASK_LOCAL       lsm2dimposemaxphimasklocal_
 #define LSM2D_CHECK_OUTER_NARROW_BAND_LAYER   lsm2dcheckouternarrowbandlayer_ 
 
#ifdef __cplusplus
extern "C" {
#endif
 
 void LSM2D_DETERMINE_NARROW_BAND(
 const LSMLIB_REAL *phi,
 const int *ilo_gb, 
 const int *ihi_gb,
 const int *jlo_gb, 
 const int *jhi_gb,
 unsigned char *narrow_band,
 const int *ilo_nb_gb, 
 const int *ihi_nb_gb,
 const int *jlo_nb_gb, 
 const int *jhi_nb_gb,
 int *index_x,
 int *index_y, 
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
 
 void LSM2D_MARK_NARROW_BAND_BOUNDARY_LAYER(
 unsigned char *narrow_band,
 const int *ilo_gb, 
 const int *ihi_gb,
 const int *jlo_gb, 
 const int *jhi_gb,
 const int *ilo_fb, 
 const int *ihi_fb,
 const int *jlo_fb, 
 const int *jhi_fb,
 const unsigned char *mark_boundary_layer);
 
 void LSM2D_DETERMINE_NARROW_BAND_FROM_MASK(
 const LSMLIB_REAL *phi,
 const LSMLIB_REAL *mask,
 const int *ilo_gb, 
 const int *ihi_gb,
 const int *jlo_gb, 
 const int *jhi_gb,
 unsigned char *narrow_band,
 const int *ilo_nb_gb, 
 const int *ihi_nb_gb,
 const int *jlo_nb_gb, 
 const int *jhi_nb_gb,
 int *index_x,
 int *index_y, 
 int *nlo_index, 
 int *nhi_index, 
 int *n_lo,
 int *n_hi,
 const int *level,
 const int *use_mask_sign);
 
 void LSM2D_MULTIPLY_CUT_OFF_LSE_RHS_LOCAL(
 const LSMLIB_REAL *phi,
 LSMLIB_REAL *lse_rhs,
 const int *ilo_gb, 
 const int *ihi_gb,
 const int *jlo_gb, 
 const int *jhi_gb,
 const int *index_x,
 const int *index_y, 
 const int *nlo_index, 
 const int *nhi_index,
 const unsigned char *narrow_band,
 const int *ilo_nb_gb, 
 const int *ihi_nb_gb,
 const int *jlo_nb_gb, 
 const int *jhi_nb_gb,
 const unsigned char *mark_fb,
 const LSMLIB_REAL *beta,
 const LSMLIB_REAL *gamma);
 
 void  LSM2D_IMPOSE_MAX_PHI_MASK_LOCAL(
 LSMLIB_REAL *phi,
 const LSMLIB_REAL *mask,
 const int *ilo_gb, 
 const int *ihi_gb,
 const int *jlo_gb, 
 const int *jhi_gb,
 const int *index_x,
 const int *index_y,
 const int *nlo_index, 
 const int *nhi_index);
 
 void LSM2D_CHECK_OUTER_NARROW_BAND_LAYER(
 int  *change_sign,
 const LSMLIB_REAL *phi,
 const int *ilo_gb, 
 const int *ihi_gb,
 const int *jlo_gb, 
 const int *jhi_gb,
 int *index_x,
 int *index_y, 
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
