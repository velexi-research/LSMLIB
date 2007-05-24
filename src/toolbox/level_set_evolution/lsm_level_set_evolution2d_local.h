#ifndef INCLUDED_LSM_LEVEL_SET_EVOLUTION_2D_LOCAL_H
#define INCLUDED_LSM_LEVEL_SET_EVOLUTION_2D_LOCAL_H

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
  double *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const int *index_x,
  const int *index_y,
  const int *nlo_index,
  const int *nhi_index);
  
void LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS_LOCAL(
  double *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const double *phi_x, 
  const double *phi_y,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const double *vel_x, 
  const double *vel_y,
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
  double *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const double *phi_x_plus, 
  const double *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const double *phi_x_minus, 
  const double *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const double *vel_n,
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
  double *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const double *phi_x_plus, 
  const double *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const double *phi_x_minus, 
  const double *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,  
  const double *vel_n, 
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
  double  *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const double *phi_x, 
  const double *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const  double *phi_xx,
  const  double *phi_xy,
  const  double *phi_yy,
  const  int *ilo_grad2_phi_gb,
  const  int *ihi_grad2_phi_gb,
  const  int *jlo_grad2_phi_gb, 
  const  int *jhi_grad2_phi_gb,
  const double *b,
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
