/*
 * File:        lsm_reinitialization2d.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.12 $
 * Modified:    $Date: 2006/05/18 23:20:13 $
 * Description: Header file for Fortran 77 2D reinitialization routines
 */

#ifndef INCLUDED_LSM_REINITIALIZATION_2D_H
#define INCLUDED_LSM_REINITIALIZATION_2D_H

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_reinitialization2d.h
 *
 * \brief
 * @ref lsm_reinitialization2d.h provides support for computing the right-hand
 * side of the reinitialization and orthogonalization equations in two 
 * space dimensions.
 *
 */


/* Link between C/C++ and Fortran function names 
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS                               \
                                     lsm2dcomputereinitializationeqnrhs_
#define LSM2D_COMPUTE_ORTHOGONALIZATION_EQN_RHS                             \
                                     lsm2dcomputeorthogonalizationeqnrhs_


/*!
 * LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS() computes the right-hand side 
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
 *   behavior of LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS() is to use
 *   \f$ \phi \f$ (i.e. equivalent to setting use_phi0_for_sgn to 0)
 *
 */
void LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS(
  double *reinit_rhs,
  const int *ilo_rhs_gb, 
  const int *ihi_rhs_gb, 
  const int *jlo_rhs_gb, 
  const int *jhi_rhs_gb,
  const double* phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb, 
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const double* phi0,
  const int *ilo_phi0_gb, 
  const int *ihi_phi0_gb, 
  const int *jlo_phi0_gb, 
  const int *jhi_phi0_gb,
  const double *phi_x_plus, 
  const double *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,   
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,   
  const int *jhi_grad_phi_plus_gb,
  const double *phi_x_minus, 
  const double *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,   
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,   
  const int *jhi_grad_phi_minus_gb,
  const int *ilo_fb,   
  const int *ihi_fb,
  const int *jlo_fb,   
  const int *jhi_fb,
  const double *dx, 
  const double *dy,
  const int *use_phi0_for_sgn);


/*!
 * LSM2D_COMPUTE_ORTHOGONALIZATION_EQN_RHS() computes the right-hand side
 * of the orthogonalization equation:
 *
 * \f[
 *
 *   \phi_t + \nabla \phi \cdot 
 *      \left ({ \frac{sgn(\psi)}{|\nabla \psi|} \nabla \psi }\right) = 0
 *
 * \f]
 *
 * Upwinding is used to select whether the forward (plus) or backward
 * (minus) spatial derivative should be used for \f$ \nabla \phi\f$.  
 * \f$ \nabla \psi \f$ is computed by averaging the forward and backward 
 * spatial derivatives for \f$ \nabla \psi \f$.  Forward and backward 
 * spatial derivatives used in the calculation must be supplied by the user.
 *
 * Arguments:
 *  - othro_rhs (out):        right-hand side of orthogonalization
 *                            equation
 *  - psi (in):               data array for \f$ \psi \f$
 *  - phi_*_plus (in):        forward spatial derivatives for 
 *                            \f$ \nabla \phi \f$
 *  - phi_*_minus (in):       backward spatial derivatives for
 *                            \f$ \nabla \phi \f$
 *  - psi_*_plus (in):        forward spatial derivatives for 
 *                            \f$ \nabla \psi \f$
 *  - psi_*_minus (in):       backward spatial derivatives for
 *                            \f$ \nabla \psi \f$
 *  - *_gb (in):              index range for ghostbox
 *  - *_fb (in):              index range for fillbox
 *
 *  Return value:             none
 */
void LSM2D_COMPUTE_ORTHOGONALIZATION_EQN_RHS(
  double *ortho_rhs,
  const int *ilo_rhs_gb, 
  const int *ihi_rhs_gb, 
  const int *jlo_rhs_gb, 
  const int *jhi_rhs_gb,
  const double *phi_x_plus, 
  const double *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,   
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,   
  const int *jhi_grad_phi_plus_gb,
  const double *phi_x_minus, 
  const double *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,   
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,   
  const int *jhi_grad_phi_minus_gb,
  const double* psi,
  const int *ilo_psi_gb, 
  const int *ihi_psi_gb, 
  const int *jlo_psi_gb, 
  const int *jhi_psi_gb,
  const double *psi_x_plus, 
  const double *psi_y_plus,
  const int *ilo_grad_psi_plus_gb,   
  const int *ihi_grad_psi_plus_gb,
  const int *jlo_grad_psi_plus_gb,   
  const int *jhi_grad_psi_plus_gb,
  const double *psi_x_minus, 
  const double *psi_y_minus,
  const int *ilo_grad_psi_minus_gb,   
  const int *ihi_grad_psi_minus_gb,
  const int *jlo_grad_psi_minus_gb,   
  const int *jhi_grad_psi_minus_gb,
  const int *ilo_fb,   
  const int *ihi_fb,
  const int *jlo_fb,   
  const int *jhi_fb,
  const double *dx, 
  const double *dy);

#ifdef __cplusplus
}
#endif

#endif
