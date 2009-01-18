/*
 * File:        lsm_tvd_runge_kutta2d_local.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for Fortran 77 2D narrow-band TVD Runge-Kutta 
 *              routines
 */

#ifndef INCLUDED_LSM_TVD_RUNGE_KUTTA_2D_LOCAL_H
#define INCLUDED_LSM_TVD_RUNGE_KUTTA_2D_LOCAL_H


/* Link between C/C++ and Fortran function names
 *
 *      name in                             name in
 *      C/C++ code                          Fortran code
 *      ----------                          ------------
 */
#define LSM2D_RK1_STEP_LOCAL                      lsm2drk1steplocal_
#define LSM2D_TVD_RK2_STAGE1_LOCAL                lsm2dtvdrk2stage1local_
#define LSM2D_TVD_RK2_STAGE2_LOCAL                lsm2dtvdrk2stage2local_
#define LSM2D_TVD_RK3_STAGE1_LOCAL                lsm2dtvdrk3stage1local_
#define LSM2D_TVD_RK3_STAGE2_LOCAL                lsm2dtvdrk3stage2local_
#define LSM2D_TVD_RK3_STAGE3_LOCAL                lsm2dtvdrk3stage3local_

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

void LSM2D_RK1_STEP_LOCAL(
  LSMLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const LSMLIB_REAL *dt,
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


void LSM2D_TVD_RK2_STAGE1_LOCAL(
  LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const LSMLIB_REAL *dt,
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


void LSM2D_TVD_RK2_STAGE2_LOCAL(
  LSMLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const LSMLIB_REAL *dt,
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
  
  
  void LSM2D_TVD_RK3_STAGE1_LOCAL(
  LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const LSMLIB_REAL *dt,
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
  
  
  void LSM2D_TVD_RK3_STAGE2_LOCAL(
  LSMLIB_REAL *u_stage2,
  const int *ilo_u_stage2_gb,
  const int *ihi_u_stage2_gb,
  const int *jlo_u_stage2_gb,
  const int *jhi_u_stage2_gb,
  const LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const LSMLIB_REAL *dt,
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
  
  void LSM2D_TVD_RK3_STAGE3_LOCAL(
  LSMLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const LSMLIB_REAL *u_stage2,
  const int *ilo_u_stage2_gb,
  const int *ihi_u_stage2_gb,
  const int *jlo_u_stage2_gb,
  const int *jhi_u_stage2_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const LSMLIB_REAL *dt,
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
