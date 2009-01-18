/*
 * File:        velocityfield_fort.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date$
 * Description: Header file for F77 velocity field routines for 3d LSM 
 *              example problem
 */


#ifndef included_velocityfield_fort
#define included_velocityfield_fort

#include "LSMLIB_config.h"

/* Link between C/C++ and Fortran function names
 *
 *      name in               name in
 *      C/C++ code            Fortran code
 *      ----------            ------------
 */
#define UNIFORM_VELOCITY_X    uniformvelx_
#define UNIFORM_VELOCITY_Y    uniformvely_
#define UNIFORM_VELOCITY_XY   uniformvelxy_
#define ROTATING_VELOCITY     rotatingvel_
#define EXPANDING_VELOCITY    expandingvel_

void UNIFORM_VELOCITY_X(
  LSMLIB_REAL* u,
  LSMLIB_REAL* v,
  LSMLIB_REAL* w,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* klo_gb, 
  const int* khi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const int* klo_fb, 
  const int* khi_fb);

void UNIFORM_VELOCITY_Y(
  LSMLIB_REAL* u,
  LSMLIB_REAL* v,
  LSMLIB_REAL* w,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* klo_gb, 
  const int* khi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const int* klo_fb, 
  const int* khi_fb);

void UNIFORM_VELOCITY_XY(
  LSMLIB_REAL* u,
  LSMLIB_REAL* v,
  LSMLIB_REAL* w,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* klo_gb, 
  const int* khi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const int* klo_fb, 
  const int* khi_fb);

void ROTATING_VELOCITY(
  LSMLIB_REAL* u,
  LSMLIB_REAL* v,
  LSMLIB_REAL* w,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* klo_gb, 
  const int* khi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const int* klo_fb, 
  const int* khi_fb,
  const LSMLIB_REAL* dx,
  const LSMLIB_REAL* x_lower);

void EXPANDING_VELOCITY(
  LSMLIB_REAL* u,
  LSMLIB_REAL* v,
  LSMLIB_REAL* w,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* klo_gb, 
  const int* khi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const int* klo_fb, 
  const int* khi_fb,
  const LSMLIB_REAL* dx,
  const LSMLIB_REAL* x_lower,
  const LSMLIB_REAL* speed,
  const LSMLIB_REAL* omega,
  const LSMLIB_REAL* time);

#endif
