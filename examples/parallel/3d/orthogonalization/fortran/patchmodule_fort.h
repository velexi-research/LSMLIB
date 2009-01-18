/*
 * File:        patchmodule_fort.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for patch module routines for
 *              LSMLIB example program
 */

#ifndef included_patchmodule_fort
#define included_patchmodule_fort

#include "LSMLIB_config.h"

/* Link between C/C++ and Fortran function names
 *
 *      name in                              name in
 *      C/C++ code                           Fortran code
 *      ----------                           ------------
 */
#define INITIALIZE_PERIODIC_ARRAY_OF_LINES   initializeperiodicarrayoflines_

void INITIALIZE_PERIODIC_ARRAY_OF_LINES(
  const LSMLIB_REAL* phi,
  const int* ilo_phi_gb,
  const int* ihi_phi_gb,
  const int* jlo_phi_gb,
  const int* jhi_phi_gb,
  const int* klo_phi_gb,
  const int* khi_phi_gb,
  const LSMLIB_REAL* psi,
  const int* ilo_psi_gb,
  const int* ihi_psi_gb,
  const int* jlo_psi_gb,
  const int* jhi_psi_gb,
  const int* klo_psi_gb,
  const int* khi_psi_gb,
  const int* ilo_fb,
  const int* ihi_fb,
  const int* jlo_fb,
  const int* jhi_fb,
  const int* klo_fb,
  const int* khi_fb,
  const LSMLIB_REAL* x_lower,
  const LSMLIB_REAL* dx);

#endif
