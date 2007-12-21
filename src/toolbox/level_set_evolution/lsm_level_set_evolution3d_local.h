#ifndef INCLUDED_LSM_LEVEL_SET_EVOLUTION_3D_LOCAL_H
#define INCLUDED_LSM_LEVEL_SET_EVOLUTION_3D_LOCAL_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_level_set_evolution3d_local.h
 *
 * \brief
 * @ref lsm_level_set_evolution3d.h provides support for contributing to the
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
#define LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS_LOCAL lsm3dzerooutlevelseteqnrhslocal_ 
#define LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL   \
                                        lsm3daddconstnormalveltermtolserhslocal_
#define LSM3D_ADD_ADVECTION_TERM_TO_LSE_RHS_LOCAL          \
                                        lsm3daddadvectiontermtolserhslocal_
#define LSM3D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL         \
                                        lsm3daddnormalveltermtolserhslocal_					
#define LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS_LOCAL         \
                                        lsm3daddconstcurvtermtolserhslocal_				  
					

void LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS_LOCAL(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *klo_lse_rhs_gb, 
  const int *khi_lse_rhs_gb,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index,
  const int *nhi_index);


void LSM3D_ADD_ADVECTION_TERM_TO_LSE_RHS_LOCAL(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *klo_lse_rhs_gb, 
  const int *khi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x, 
  const LSMLIB_REAL *phi_y, 
  const LSMLIB_REAL *phi_z,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb, 
  const int *khi_grad_phi_gb,
  const LSMLIB_REAL *vel_x, 
  const LSMLIB_REAL *vel_y, 
  const LSMLIB_REAL *vel_z,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *klo_vel_gb, 
  const int *khi_vel_gb, 
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
  

void LSM3D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *klo_lse_rhs_gb, 
  const int *khi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x_plus, 
  const LSMLIB_REAL *phi_y_plus, 
  const LSMLIB_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb, 
  const int *khi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus, 
  const LSMLIB_REAL *phi_y_minus, 
  const LSMLIB_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb, 
  const int *khi_grad_phi_minus_gb,
  const LSMLIB_REAL *vel_n, 
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *klo_vel_gb, 
  const int *khi_vel_gb,
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
			
					
void LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *klo_lse_rhs_gb, 
  const int *khi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x_plus, 
  const LSMLIB_REAL *phi_y_plus, 
  const LSMLIB_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb, 
  const int *khi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus, 
  const LSMLIB_REAL *phi_y_minus, 
  const LSMLIB_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb, 
  const int *khi_grad_phi_minus_gb,
  const LSMLIB_REAL *vel_n, 
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


void   LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS_LOCAL(
  LSMLIB_REAL  *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *klo_lse_rhs_gb, 
  const int *khi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x, 
  const LSMLIB_REAL *phi_y,
  const LSMLIB_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const  LSMLIB_REAL *phi_xx,
  const  LSMLIB_REAL *phi_xy,
  const  LSMLIB_REAL *phi_xz,
  const  LSMLIB_REAL *phi_yy,
  const  LSMLIB_REAL *phi_yz,
  const  LSMLIB_REAL *phi_zz,
  const  int *ilo_grad2_phi_gb,
  const  int *ihi_grad2_phi_gb,
  const  int *jlo_grad2_phi_gb, 
  const  int *jhi_grad2_phi_gb,
  const  int *klo_grad2_phi_gb, 
  const  int *khi_grad2_phi_gb,
  const LSMLIB_REAL *b,
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
