/*
 * File:        advection2d_velocityfield.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Description: Header file for F77 velocity field routines for 2d LSM example
 *              problem
 */


#ifndef included_advection2d_velocityfield
#define included_advection2d_velocityfield

#include "LSMLIB_config.h"

void uniformvelx_(
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

void uniformvely_(
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

void uniformvelxy_(
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

void rotatingvel_(
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

void expandingvel_(
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
