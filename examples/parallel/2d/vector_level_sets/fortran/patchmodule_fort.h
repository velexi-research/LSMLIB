/*
 * File:        patchmodule_fort.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for patch module routines for 2d LSM 
 *              example problem
 */

#ifndef included_patchmodule_fort
#define included_patchmodule_fort

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
