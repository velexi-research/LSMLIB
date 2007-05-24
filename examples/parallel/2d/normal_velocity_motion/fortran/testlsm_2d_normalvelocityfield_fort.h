/*
 * File:        lsmtest_2d_normalvelocityfield.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Description: Header file for F77 normal velocity field routines for 2d 
 *              LSM test problem
 */


#ifndef included_lsmtest_2d_normalvelocityfield
#define included_lsmtest_2d_normalvelocityfield

/* Link between C/C++ and Fortran function names
 *
 *      name in                                  name in
 *      C/C++ code                               Fortran code
 *      ----------                               ------------
 */
#define COMPUTE_OSCILLATING_EXPANSION_VELOCITY   oscillatingexpansionvelocity_
#define COMPUTE_MEAN_CURVATURE_VELOCITY          meancurvaturevelocity_

void COMPUTE_OSCILLATING_EXPANSION_VELOCITY (
  double* normal_vel,
  const int* ilo_vel_gb, 
  const int* ihi_vel_gb, 
  const int* jlo_vel_gb, 
  const int* jhi_vel_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const double* dx,
  const double* x_lower,
  const double* speed,
  const double* omega,
  const double* time);

void COMPUTE_MEAN_CURVATURE_VELOCITY (
  double* normal_vel,
  const int* ilo_vel_gb, 
  const int* ihi_vel_gb, 
  const int* jlo_vel_gb, 
  const int* jhi_vel_gb,
  const double* normal_phi,
  const int* ilo_phi_gb, 
  const int* ihi_phi_gb, 
  const int* jlo_phi_gb, 
  const int* jhi_phi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const double* dx,
  const double* dy,
  const double* x_lower,
  const double* y_lower,
  const double* speed);

#endif
