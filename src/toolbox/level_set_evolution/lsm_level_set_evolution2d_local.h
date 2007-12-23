/*
 * File:        lsm_level_set_evolution2d_local.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.8 $
 * Modified:    $Date: 2006/10/28 05:09:39 $
 * Description: Header file for 2D Fortran 77 narrow-band, level set 
 *              evolution equation subroutines
 */

#ifndef INCLUDED_LSM_LEVEL_SET_EVOLUTION_2D_LOCAL_H
#define INCLUDED_LSM_LEVEL_SET_EVOLUTION_2D_LOCAL_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_level_set_evolution2d_local.h
 *
 * \brief
 * @ref lsm_level_set_evolution2d.h provides support for contributing to the
 * right-hand side of the level set evolution equation in three space
 * dimensions for narrow banding approach (localization).
 * 
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                                name in
 *      C/C++ code                             Fortran code
 *      ----------                             ------------
 */
#define LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS_LOCAL lsm2dzerooutlevelseteqnrhslocal_
#define LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS_LOCAL                             \
                                        lsm2daddadvectiontermtolserhslocal_
#define LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL                            \
                                          lsm2daddnormalveltermtolserhslocal_					
#define LSM2D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL                      \
                                        lsm2daddconstnormalveltermtolserhslocal_				       
#define LSM2D_ADD_CONST_CURV_TERM_TO_LSE_RHS_LOCAL                            \
                                        lsm2daddconstcurvtermtolserhslocal_				  
					

void LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS_LOCAL(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *index_x,
  const int *index_y,
  const int *nlo_index,
  const int *nhi_index);
  
void LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS_LOCAL(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x, 
  const LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *vel_x, 
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
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


void LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
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
  const LSMLIB_REAL *vel_n,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
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
  
  					
void LSM2D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
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
  const LSMLIB_REAL *vel_n, 
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
  			  

void   LSM2D_ADD_CONST_CURV_TERM_TO_LSE_RHS_LOCAL(
  LSMLIB_REAL  *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x, 
  const LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const  LSMLIB_REAL *phi_xx,
  const  LSMLIB_REAL *phi_xy,
  const  LSMLIB_REAL *phi_yy,
  const  int *ilo_grad2_phi_gb,
  const  int *ihi_grad2_phi_gb,
  const  int *jlo_grad2_phi_gb, 
  const  int *jhi_grad2_phi_gb,
  const LSMLIB_REAL *b,
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
