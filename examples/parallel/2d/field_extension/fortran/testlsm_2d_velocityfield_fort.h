/*
 * File:        lsmtest_2d_velocityfield.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date: 2006/05/18 01:09:37 $
 * Description: Header file for F77 velocity field routines for 2d LSM test 
 *              problem
 */


#ifndef included_lsmtest_2d_velocityfield
#define included_lsmtest_2d_velocityfield

void uniformvelx_(
  double* u,
  double* v,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb);

void uniformvely_(
  double* u,
  double* v,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb);

void uniformvelxy_(
  double* u,
  double* v,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb);

void rotatingvel_(
  double* u,
  double* v,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const double* dx,
  const double* x_lower);

void expandingvel_(
  double* u,
  double* v,
  const int* ilo_gb, 
  const int* ihi_gb, 
  const int* jlo_gb, 
  const int* jhi_gb,
  const int* ilo_fb, 
  const int* ihi_fb, 
  const int* jlo_fb, 
  const int* jhi_fb,
  const double* dx,
  const double* x_lower,
  const double* speed,
  const double* omega,
  const double* time);

#endif
