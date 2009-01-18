/*
 * File:        lsm_tvd_runge_kutta3d_local.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for Fortran 77 3D narrow-band TVD Runge-Kutta 
 *              routines
 */

#ifndef INCLUDED_LSM_TVD_RUNGE_KUTTA_3D_LOCAL_H
#define INCLUDED_LSM_TVD_RUNGE_KUTTA_3D_LOCAL_H


/* Link between C/C++ and Fortran function names
 *
 *      name in                             name in
 *      C/C++ code                          Fortran code
 *      ----------                          ------------
 */
#define LSM3D_RK1_STEP_LOCAL                      lsm3drk1steplocal_
#define LSM3D_TVD_RK2_STAGE1_LOCAL                lsm3dtvdrk2stage1local_
#define LSM3D_TVD_RK2_STAGE2_LOCAL                lsm3dtvdrk2stage2local_
#define LSM3D_TVD_RK3_STAGE1_LOCAL                lsm3dtvdrk3stage1local_
#define LSM3D_TVD_RK3_STAGE2_LOCAL                lsm3dtvdrk3stage2local_
#define LSM3D_TVD_RK3_STAGE3_LOCAL                lsm3dtvdrk3stage3local_

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

void LSM3D_RK1_STEP_LOCAL(
  LSMLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const int *klo_u_next_gb,
  const int *khi_u_next_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const int *klo_u_cur_gb,
  const int *khi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const LSMLIB_REAL *dt,
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
  const unsigned char *mark_fb);


void LSM3D_TVD_RK2_STAGE1_LOCAL(
  LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const int *klo_u_stage1_gb,
  const int *khi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const int *klo_u_cur_gb,
  const int *khi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const LSMLIB_REAL *dt,
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
  const unsigned char *mark_fb);


void LSM3D_TVD_RK2_STAGE2_LOCAL(
  LSMLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const int *klo_u_next_gb,
  const int *khi_u_next_gb,
  const LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const int *klo_u_stage1_gb,
  const int *khi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const int *klo_u_cur_gb,
  const int *khi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const LSMLIB_REAL *dt,
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
  const unsigned char *mark_fb);
  
void LSM3D_TVD_RK3_STAGE1_LOCAL(
  LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const int *klo_u_stage1_gb,
  const int *khi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const int *klo_u_cur_gb,
  const int *khi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const LSMLIB_REAL *dt,
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
  const unsigned char *mark_fb);  
      
 void LSM3D_TVD_RK3_STAGE2_LOCAL(
  LSMLIB_REAL *u_stage2,
  const int *ilo_u_stage2_gb,
  const int *ihi_u_stage2_gb,
  const int *jlo_u_stage2_gb,
  const int *jhi_u_stage2_gb,
  const int *klo_u_stage2_gb,
  const int *khi_u_stage2_gb,
  const LSMLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const int *klo_u_stage1_gb,
  const int *khi_u_stage1_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const int *klo_u_cur_gb,
  const int *khi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const LSMLIB_REAL *dt,
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
  const unsigned char *mark_fb);
  
 
 void LSM3D_TVD_RK3_STAGE3_LOCAL(
  LSMLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const int *klo_u_next_gb,
  const int *khi_u_next_gb,
  const LSMLIB_REAL *u_stage2,
  const int *ilo_u_stage2_gb,
  const int *ihi_u_stage2_gb,
  const int *jlo_u_stage2_gb,
  const int *jhi_u_stage2_gb,
  const int *klo_u_stage2_gb,
  const int *khi_u_stage2_gb,
  const LSMLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const int *klo_u_cur_gb,
  const int *khi_u_cur_gb,
  const LSMLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const LSMLIB_REAL *dt,
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
  const unsigned char *mark_fb);
       
#ifdef __cplusplus
}
#endif
  
#endif
