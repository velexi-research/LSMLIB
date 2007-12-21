/*
 * File:        lsm_utilities2d_local.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.23 $
 * Modified:    $Date: 2006/10/28 04:54:37 $
 * Description: Header file for 2D Fortran 77 level set method narrow-band
 *              utility subroutines
 */

#ifndef INCLUDED_LSM_UTILITIES_2D_LOCAL_H
#define INCLUDED_LSM_UTILITIES_2D_LOCAL_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_utilities2d.h
 *
 * \brief 
 * @ref lsm_utilities2d.h provides several utility functions that support
 * level set method calculations in three space dimensions.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                        name in
 *      C/C++ code                     Fortran code
 *      ----------                     ------------
 */
#define LSM2D_MAX_NORM_DIFF_LOCAL            lsm2dmaxnormdifflocal_
#define LSM2D_COMPUTE_STABLE_ADVECTION_DT_LOCAL                         \
                                       lsm2dcomputestableadvectiondtlocal_
#define LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT_LOCAL                        \
                                       lsm2dcomputestablenormalveldtlocal_
#define LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_LOCAL                  \
                                        lsm2dcomputestableconstnormalveldtlocal_

void LSM2D_MAX_NORM_DIFF_LOCAL(
  LSMLIB_REAL *max_norm_diff,
  const LSMLIB_REAL *field1,
  const int *ilo_field1_gb, 
  const int *ihi_field1_gb,
  const int *jlo_field1_gb, 
  const int *jhi_field1_gb,
  const LSMLIB_REAL *field2,
  const int *ilo_field2_gb, 
  const int *ihi_field2_gb,
  const int *jlo_field2_gb, 
  const int *jhi_field2_gb,
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

void LSM2D_COMPUTE_STABLE_ADVECTION_DT_LOCAL(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_x,
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *cfl_number,
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

void LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT_LOCAL(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const LSMLIB_REAL *phi_x_plus,
  const LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus,
  const LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *cfl_number,
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
    
void LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_LOCAL(
  LSMLIB_REAL *dt,
  const LSMLIB_REAL *vel_n,
  const LSMLIB_REAL *phi_x_plus,
  const LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus,
  const LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy,
  const LSMLIB_REAL *cfl_number,  
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
