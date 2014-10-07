/*
 * File:        velocityfield_fort.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Description: Header file for F77 velocity field routines for 2d LSM 
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
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb);

void UNIFORM_VELOCITY_Y(
  LSMLIB_REAL* u,
  LSMLIB_REAL* v,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb);

void UNIFORM_VELOCITY_XY(
  LSMLIB_REAL* u,
  LSMLIB_REAL* v,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb);

void ROTATING_VELOCITY(
  LSMLIB_REAL* u,
  LSMLIB_REAL* v,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const LSMLIB_REAL* dx,
  const LSMLIB_REAL* x_lower);

void EXPANDING_VELOCITY(
  LSMLIB_REAL* u,
  LSMLIB_REAL* v,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const LSMLIB_REAL* dx,
  const LSMLIB_REAL* x_lower,
  const LSMLIB_REAL* speed,
  const LSMLIB_REAL* omega,
  const LSMLIB_REAL* time);

#endif
