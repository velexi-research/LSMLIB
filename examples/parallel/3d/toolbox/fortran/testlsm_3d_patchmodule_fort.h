/*
 * File:        testlsm_3d_patchmodule.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.2 $
 * Modified:    $Date: 2006/01/24 21:46:15 $
 * Description: Header file for patch module routines for 3d LSM test problem
 */

#ifndef included_lsmtest_3d_patchmodule
#define included_lsmtest_3d_patchmodule

/* Link between C/C++ and Fortran function names
 *
 *      name in               name in
 *      C/C++ code            Fortran code
 *      ----------            ------------
 */
#define INIT_SPHERE           initsphere_

void INIT_SPHERE(
  const double* level_set,
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
  const double* x_lower,
  const double* dx,
  const double* center,
  const double* radius);

#endif
