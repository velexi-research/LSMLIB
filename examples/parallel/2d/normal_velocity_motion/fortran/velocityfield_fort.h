/*
 * File:        velocityfield_fort.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Description: Header file for F77 normal velocity field routines for 2d 
 *              LSM example problem
 */


#ifndef included_velocityfield_fort
#define included_velocityfield_fort

#include "LSMLIB_config.h"

/* Link between C/C++ and Fortran function names
 *
 *      name in                                  name in
 *      C/C++ code                               Fortran code
 *      ----------                               ------------
 */
#define COMPUTE_OSCILLATING_EXPANSION_VELOCITY   oscillatingexpansionvelocity_
#define COMPUTE_MEAN_CURVATURE_VELOCITY          meancurvaturevelocity_

void COMPUTE_OSCILLATING_EXPANSION_VELOCITY (
  LSMLIB_REAL* normal_vel,
  const int* ilo_vel_gb, 
  const int* ihi_vel_gb, 
  const int* jlo_vel_gb, 
  const int* jhi_vel_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const LSMLIB_REAL* x_lower,
  const LSMLIB_REAL* speed,
  const LSMLIB_REAL* omega,
  const LSMLIB_REAL* time);

void COMPUTE_MEAN_CURVATURE_VELOCITY (
  LSMLIB_REAL* normal_vel,
  const int* ilo_vel_gb, 
  const int* ihi_vel_gb, 
  const int* jlo_vel_gb, 
  const int* jhi_vel_gb,
  const LSMLIB_REAL* normal_phi,
  const int* ilo_phi_gb, 
  const int* ihi_phi_gb, 
  const int* jlo_phi_gb, 
  const int* jhi_phi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const LSMLIB_REAL* dx,
  const LSMLIB_REAL* dy,
  const LSMLIB_REAL* speed);

#endif
