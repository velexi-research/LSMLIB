/*
 * File:        lsm_reinitialization1d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for Fortran 77 1D reinitialization routines
 */

#ifndef INCLUDED_LSM_REINITIALIZATION_1D_H
#define INCLUDED_LSM_REINITIALIZATION_1D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_reinitialization1d.h
 *
 * \brief
 * @ref lsm_reinitialization1d.h provides support for computing the right-hand
 * side of the reinitialization equation in one space dimension.
 *
 */


/* Link between C/C++ and Fortran function names 
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define LSM1D_COMPUTE_REINITIALIZATION_EQN_RHS                              \
                                     lsm1dcomputereinitializationeqnrhs_


/*!
 * LSM1D_COMPUTE_REINITIALIZATION_EQN_RHS() computes the right-hand side 
 * of the reinitialization equation using a Godunov scheme to select the
 * numerical discretization of the \f$ sgn(\phi) |\nabla \phi| \f$ term.
 * Forward (plus) and backward (minus) spatial derivatives used in
 * the Godunov calculation must be supplied by the user.
 *
 * Arguments:
 *  - reinit_rhs (out):       right-hand side of reinitialization
 *                            equation
 *  - phi (in):               level set function at current iteration
 *                            of reinitialization process
 *  - phi0 (in):              level set function at initial iteration
 *                            iteration of reinitialization process
 *  - phi_*_plus (in):        forward spatial derivatives for 
 *                            \f$ \nabla \phi \f$
 *  - phi_*_minus (in):       backward spatial derivatives for 
 *                            \f$ \nabla \phi \f$
 *  - use_phi0_for_sgn (in):  flag to specify whether \f$ \phi_0 \f$ should 
 *                            be used in the computation of \f$ sgn(\phi) \f$.
 *                              0 = use \f$ \phi \f$ (do NOT use 
 *                                  \f$ \phi_0 \f$ );
 *                              1 = use \f$ \phi_0 \f$.
 *  - *_gb (in):              index range for ghostbox
 *  - *_fb (in):              index range for fillbox
 *
 * NOTES:
 * - if use_phi0_for_sgn is not equal to 0 or 1, the default
 *   behavior of LSM1D_COMPUTE_REINITIALIZATION_EQN_RHS is to use
 *   \f$ \phi \f$ (i.e. equivalent to setting use_phi0_for_sgn to 0)
 *
 */
void LSM1D_COMPUTE_REINITIALIZATION_EQN_RHS(
  LSMLIB_REAL *reinit_rhs,
  const int *ilo_rhs_gb, 
  const int *ihi_rhs_gb, 
  const LSMLIB_REAL* phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb, 
  const LSMLIB_REAL* phi0,
  const int *ilo_phi0_gb, 
  const int *ihi_phi0_gb, 
  const LSMLIB_REAL *phi_x_plus, 
  const int *ilo_grad_phi_plus_gb,   
  const int *ihi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus, 
  const int *ilo_grad_phi_minus_gb,   
  const int *ihi_grad_phi_minus_gb,
  const int *ilo_fb,   
  const int *ihi_fb,
  const LSMLIB_REAL *dx, 
  const int *use_phi0_for_sgn);

#ifdef __cplusplus
}
#endif

#endif
