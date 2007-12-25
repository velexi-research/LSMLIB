/*
 * File:        advection2d_patchmodule.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.2 $
 * Modified:    $Date: 2006/01/24 21:45:57 $
 * Description: Header file for patch module routines for 2d LSM example 
 *              problem
 */

#ifndef included_advection2d_patchmodule
#define included_advection2d_patchmodule

#include "LSMLIB_config.h"

/* Link between C/C++ and Fortran function names
 *
 *      name in               name in
 *      C/C++ code            Fortran code
 *      ----------            ------------
 */
#define INIT_CIRCLE           initcircle_

void INIT_CIRCLE(
  const LSMLIB_REAL* level_set,
  const int* ilo_gb,
  const int* ihi_gb,
  const int* jlo_gb,
  const int* jhi_gb,
  const int* ilo_fb,
  const int* ihi_fb,
  const int* jlo_fb,
  const int* jhi_fb,
  const LSMLIB_REAL* x_lower,
  const LSMLIB_REAL* dx,
  const LSMLIB_REAL* center,
  const LSMLIB_REAL* radius);

#endif
