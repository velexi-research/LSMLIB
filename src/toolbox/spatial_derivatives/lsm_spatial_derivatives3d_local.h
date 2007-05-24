#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_spatial_derivatives3d_local.h
 *
 * \brief
 * @ref lsm_spatial_derivatives3d_local.h provides support for computing spatial
 * derivatives in three space dimensions using high-order ENO and WENO 
 * discretizations for narrow banding (local) approach. 
 *
 */
 
#ifndef INCLUDED_LSM_SPATIAL_DERIVATIVES_3D_LOCAL_H
#define INCLUDED_LSM_SPATIAL_DERIVATIVES_3D_LOCAL_H



/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define LSM3D_HJ_ENO1_LOCAL              lsm3dhjeno1local_
#define LSM3D_HJ_ENO2_LOCAL              lsm3dhjeno2local_
#define LSM3D_CENTRAL_GRAD_ORDER2_LOCAL  lsm3dcentralgradorder2local_
#define LSM3D_CENTRAL_GRAD_ORDER4_LOCAL  lsm3dcentralgradorder4local_
#define LSM3D_LAPLACIAN_ORDER2_LOCAL     lsm3dlaplacianorder2local_
#define LSM3D_COMPUTE_AVE_GRAD_PHI_LOCAL lsm3dcomputeavegradphilocal_

void LSM3D_HJ_ENO1_LOCAL(
  double *phi_x_plus,
  double *phi_y_plus,
  double *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,
  const int *khi_grad_phi_plus_gb,
  double *phi_x_minus,
  double *phi_y_minus,
  double *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,
  const int *khi_grad_phi_minus_gb,
  const double *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  double *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  const double *dx,
  const double *dy,
  const double *dz,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index0,
  const int *nhi_index0,
  const int *nlo_index1,
  const int *nhi_index1,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const int *klo_nb_gb,
  const int *khi_nb_gb,
  const unsigned char *mark_fb,
  const unsigned char *mark_D1);

void LSM3D_HJ_ENO2_LOCAL(
  double *phi_x_plus,
  double *phi_y_plus,
  double *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,
  const int *khi_grad_phi_plus_gb,
  double *phi_x_minus,
  double *phi_y_minus,
  double *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,
  const int *khi_grad_phi_minus_gb,
  const double *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  double *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  double *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const int *klo_D2_gb,
  const int *khi_D2_gb,
  const double *dx,
  const double *dy,
  const double *dz,
  const int *index_x,
  const int *index_y, 
  const int *index_z,
  const int *nlo_index0,
  const int *nhi_index0,
  const int *nlo_index1,
  const int *nhi_index1,
  const int *nlo_index2,
  const int *nhi_index2,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const int *klo_nb_gb,
  const int *khi_nb_gb,
  const unsigned char *mark_fb,
  const unsigned char *mark_D1,
  const unsigned char *mark_D2);


void LSM3D_CENTRAL_GRAD_ORDER2_LOCAL( 
  double *phi_x,
  double *phi_y,
  double *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const double *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const double *dx,
  const double *dy,
  const double *dz,
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
  
void LSM3D_CENTRAL_GRAD_ORDER4_LOCAL( 
  double *phi_x,
  double *phi_y,
  double *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const double *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const double *dx,
  const double *dy,
  const double *dz,
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
 
void LSM3D_LAPLACIAN_ORDER2_LOCAL( 
  double *laplacian_phi,
  const int *ilo_laplacian_phi_gb,
  const int *ihi_laplacian_phi_gb,
  const int *jlo_laplacian_phi_gb,
  const int *jhi_laplacian_phi_gb,
  const int *klo_laplacian_phi_gb,
  const int *khi_laplacian_phi_gb,
  const double *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const double *dx,
  const double *dy,
  const double *dz,
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
      
void LSM3D_COMPUTE_AVE_GRAD_PHI_LOCAL(
  double *grad_phi_ave,
  const double *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const double *dx,
  const double *dy,
  const double *dz,
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
