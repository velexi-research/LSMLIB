/*
 * File:        lsm_samrai_f77_utilities.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for Fortran 77 3D ENO/WENO routines.
 */

#ifndef INCLUDED_LSM_SAMRAI_F77_UTILITIES
#define INCLUDED_LSM_SAMRAI_F77_UTILITIES

#include "LSMLIB_config.h"

/* Link between C/C++ and Fortran function names
 *
 *      name in                            name in
 *      C/C++ code                         Fortran code
 *      ----------                         ------------
 */
#define LSM1D_SAMRAI_UTILITIES_COPY_DATA   lsm1dsamraiutilitiescopydata_
#define LSM2D_SAMRAI_UTILITIES_COPY_DATA   lsm2dsamraiutilitiescopydata_
#define LSM3D_SAMRAI_UTILITIES_COPY_DATA   lsm3dsamraiutilitiescopydata_

void LSM1D_SAMRAI_UTILITIES_COPY_DATA(
  LSMLIB_REAL *dst_data,
  const int *ilo_dst_gb,
  const int *ihi_dst_gb,
  LSMLIB_REAL *src_data,
  const int *ilo_src_gb,
  const int *ihi_src_gb,
  const int *ilo_fb,
  const int *ihi_fb);

void LSM2D_SAMRAI_UTILITIES_COPY_DATA(
  LSMLIB_REAL *dst_data,
  const int *ilo_dst_gb,
  const int *ihi_dst_gb,
  const int *jlo_dst_gb,
  const int *jhi_dst_gb,
  LSMLIB_REAL *src_data,
  const int *ilo_src_gb,
  const int *ihi_src_gb,
  const int *jlo_src_gb,
  const int *jhi_src_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb);

void LSM3D_SAMRAI_UTILITIES_COPY_DATA(
  LSMLIB_REAL *dst_data,
  const int *ilo_dst_gb,
  const int *ihi_dst_gb,
  const int *jlo_dst_gb,
  const int *jhi_dst_gb,
  const int *klo_dst_gb,
  const int *khi_dst_gb,
  LSMLIB_REAL *src_data,
  const int *ilo_src_gb,
  const int *ihi_src_gb,
  const int *jlo_src_gb,
  const int *jhi_src_gb,
  const int *klo_src_gb,
  const int *khi_src_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb);

#endif
