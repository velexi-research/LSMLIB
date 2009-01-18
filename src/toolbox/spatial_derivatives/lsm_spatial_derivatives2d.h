/*
 * File:        lsm_spatial_derivatives2d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for Fortran 77 2D ENO/WENO routines.
 */

#ifndef INCLUDED_LSM_SPATIAL_DERIVATIVES_2D_H
#define INCLUDED_LSM_SPATIAL_DERIVATIVES_2D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_spatial_derivatives2d.h
 *
 * \brief
 * @ref lsm_spatial_derivatives2d.h provides support for computing spatial
 * derivatives in two space dimensions using high-order ENO and WENO 
 * discretizations.  
 *
 */


/* Link between C/C++ and Fortran function names 
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define LSM2D_HJ_ENO1                lsm2dhjeno1_
#define LSM2D_HJ_ENO2                lsm2dhjeno2_
#define LSM2D_HJ_ENO3                lsm2dhjeno3_
#define LSM2D_HJ_WENO5               lsm2dhjweno5_
#define LSM2D_UPWIND_HJ_ENO1         lsm2dupwindhjeno1_
#define LSM2D_UPWIND_HJ_ENO2         lsm2dupwindhjeno2_
#define LSM2D_UPWIND_HJ_ENO3         lsm2dupwindhjeno3_
#define LSM2D_UPWIND_HJ_WENO5        lsm2dupwindhjweno5_
#define LSM2D_CENTRAL_GRAD_ORDER2    lsm2dcentralgradorder2_
#define LSM2D_CENTRAL_GRAD_ORDER4    lsm2dcentralgradorder4_
#define LSM2D_LAPLACIAN_ORDER2       lsm2dlaplacianorder2_
#define LSM2D_PHI_UPWIND_GRAD_F      lsm2dphiupwindgradf_
#define LSM2D_AVERAGE_GRAD_PHI       lsm2daveragegradphi_
#define LSM2D_GRADIENT_MAGNITUDE     lsm2dgradientmagnitude_
#define LSM2D_DIVERGENCE_CENTRAL     lsm2ddivergencecentral_


/*!
 * LSM2D_HJ_ENO1() computes the forward (plus) and backward (minus)
 * first-order Hamilton-Jacobi ENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - dx, dy (in):        grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *      
 * Return value:          none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 *  - phi_x_plus and phi_x_minus are face-centered data (i.e. their
 *    indices are of the form (i+1/2) and (i-1/2)).  For phi_x_plus,
 *    the data array position corresponding to the (i+1/2) is i (i.e. 
 *    the index shifted down to the nearest integer index).  For 
 *    phi_x_minus, the data array position corresponding to the (i-1/2) 
 *    is i (i.e.  the index shifted up to the nearest integer index).  
 *    Analogous conventions hold for phi_y_plus and phi_y_minus.
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void LSM2D_HJ_ENO1(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*!
 * LSM2D_HJ_ENO2() computes the forward (plus) and backward (minus)
 * second-order Hamilton-Jacobi ENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - D2 (in):            scratch space for holding undivided 
 *                        second-differences
 *  - dx, dy (in):        grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *      
 * Return value:          none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 *  - phi_x_plus and phi_x_minus are face-centered data (i.e. their
 *    indices are of the form (i+1/2) and (i-1/2)).  For phi_x_plus,
 *    the data array position corresponding to the (i+1/2) is i (i.e. 
 *    the index shifted down to the nearest integer index).  For 
 *    phi_x_minus, the data array position corresponding to the (i-1/2) 
 *    is i (i.e.  the index shifted up to the nearest integer index).  
 *    Analogous conventions hold for phi_y_plus and phi_y_minus.
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void LSM2D_HJ_ENO2(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*!
 * LSM2D_HJ_ENO3() computes the forward (plus) and backward (minus)
 * third-order Hamilton-Jacobi ENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - D2 (in):            scratch space for holding undivided 
 *                        second-differences
 *  - D3 (in):            scratch space for holding undivided third-differences
 *  - dx, dy (in):        grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *      
 * Return value:          none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 *  - phi_x_plus and phi_x_minus are face-centered data (i.e. their
 *    indices are of the form (i+1/2) and (i-1/2)).  For phi_x_plus,
 *    the data array position corresponding to the (i+1/2) is i (i.e. 
 *    the index shifted down to the nearest integer index).  For 
 *    phi_x_minus, the data array position corresponding to the (i-1/2) 
 *    is i (i.e.  the index shifted up to the nearest integer index).  
 *    Analogous conventions hold for phi_y_plus and phi_y_minus.
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void LSM2D_HJ_ENO3(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  LSMLIB_REAL *D3,
  const int *ilo_D3_gb,
  const int *ihi_D3_gb,
  const int *jlo_D3_gb,
  const int *jhi_D3_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*!
 * LSM2D_HJ_WENO5() computes the forward (plus) and backward (minus)
 * fifth-order Hamilton-Jacobi WENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - dx, dy (in):        grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *      
 * Return value:          none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 *  - phi_x_plus and phi_x_minus are face-centered data (i.e. their
 *    indices are of the form (i+1/2) and (i-1/2)).  For phi_x_plus,
 *    the data array position corresponding to the (i+1/2) is i (i.e. 
 *    the index shifted down to the nearest integer index).  For 
 *    phi_x_minus, the data array position corresponding to the (i-1/2) 
 *    is i (i.e.  the index shifted up to the nearest integer index).  
 *    Analogous conventions hold for phi_y_plus and phi_y_minus.
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void LSM2D_HJ_WENO5(
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*! 
 * LSM2D_UPWIND_HJ_ENO1() computes the first-order Hamilton-Jacobi ENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_* (out):  components of \f$ \nabla \phi \f$
 *  - phi (in):     \f$ \phi \f$
 *  - vel_* (in):   components of the velocity 
 *  - D1 (in):      scratch space for holding undivided first-differences
 *  - dx, dy (in):  grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:    none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void LSM2D_UPWIND_HJ_ENO1(
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const LSMLIB_REAL *vel_x,
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  const int *jlo_vel_gb,
  const int *jhi_vel_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*! 
 * LSM2D_UPWIND_HJ_ENO2() computes the second-order Hamilton-Jacobi ENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_* (out):  components of \f$ \nabla \phi \f$
 *  - phi (in):     \f$ \phi \f$
 *  - vel_* (in):   components of the velocity 
 *  - D1 (in):      scratch space for holding undivided first-differences
 *  - D2 (in):      scratch space for holding undivided second-differences
 *  - dx, dy (in):  grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:    none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void LSM2D_UPWIND_HJ_ENO2(
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const LSMLIB_REAL *vel_x,
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  const int *jlo_vel_gb,
  const int *jhi_vel_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*! 
 * LSM2D_UPWIND_HJ_ENO3() computes the third-order Hamilton-Jacobi ENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_* (out):  components of \f$ \nabla \phi \f$
 *  - phi (in):     \f$ \phi \f$
 *  - vel_* (in):   components of the velocity 
 *  - D1 (in):      scratch space for holding undivided first-differences
 *  - D2 (in):      scratch space for holding undivided second-differences
 *  - D3 (in):      scratch space for holding undivided third-differences
 *  - dx, dy (in):  grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:    none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void LSM2D_UPWIND_HJ_ENO3(
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const LSMLIB_REAL *vel_x,
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  const int *jlo_vel_gb,
  const int *jhi_vel_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  LSMLIB_REAL *D3,
  const int *ilo_D3_gb,
  const int *ihi_D3_gb,
  const int *jlo_D3_gb,
  const int *jhi_D3_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*! 
 * LSM2D_UPWIND_HJ_WENO5() computes the fifth-order Hamilton-Jacobi WENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_* (out):  components of \f$ \nabla \phi \f$
 *  - phi (in):     \f$ \phi \f$
 *  - vel_* (in):   components of the velocity 
 *  - D1 (in):      scratch space for holding undivided first-differences
 *  - dx, dy (in):  grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:    none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void LSM2D_UPWIND_HJ_WENO5(
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const LSMLIB_REAL *vel_x,
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  const int *jlo_vel_gb,
  const int *jhi_vel_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*! 
 * LSM2D_CENTRAL_GRAD_ORDER2() computes the second-order, central,
 * finite difference approximation to the gradient of \f$ \phi \f$ 
 * using the formula:
 *
 *    \f[
 *
 *      \left( \frac{\partial \phi}{\partial x} \right)_i \approx
 *        \frac{ \phi_{i+1} - \phi_{i-1} }{ 2 dx }
 *
 *    \f]
 * 
 * Arguments:
 *  - phi_* (out):  components of \f$ \nabla \phi \f$
 *  - phi (in):     \f$ \phi \f$
 *  - dx, dy (in):  grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:        none
 */
void LSM2D_CENTRAL_GRAD_ORDER2( 
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*! 
 * LSM2D_CENTRAL_GRAD_ORDER4() computes the fourth-order, central, 
 * finite difference approximation to the gradient of \f$ \phi \f$ 
 * using the formula:
 *
 *    \f[
 *
 *      \left( \frac{\partial \phi}{\partial x} \right)_i \approx
 *         \frac{ -\phi_{i+2} + 8 \phi_{i+1} - 8 \phi_{i-1} + \phi_{i-2} }
 *              { 12 dx }
 *
 *    \f]
 *
 * Arguments:
 *  - phi_* (out):  components of \f$ \nabla \phi \f$
 *  - phi (in):     \f$ \phi \f$
 *  - dx, dy (in):  grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:        none
 */
void LSM2D_CENTRAL_GRAD_ORDER4( 
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*! 
 * LSM2D_LAPLACIAN_ORDER2() computes the second-order, central, 
 * finite difference approximation to the Laplacian of \f$ \phi \f$ 
 * using the formula:
 *
 *    \f[
 *
 *      \nabla^2 \phi \approx
 *         \frac{ \phi_{i+1,j,k} - 2 \phi_{i,j,k} + \phi_{i-1,j,k} }
 *              { dx^2 }
 *       + \frac{ \phi_{i,j+1,k} - 2 \phi_{i,j,k} + \phi_{i,j-1,k} }
 *              { dy^2 }
 *
 *    \f]
 *
 * Arguments:
 *  - laplacian_phi (out):  Laplacian of \f$ phi \f$
 *  - phi (in):             \f$ \phi \f$
 *  - dx (in):              grid cell size
 *  - *_gb (in):            index range for ghostbox
 *  - *_fb (in):            index range for fillbox
 *
 * Return value:            none
 */
void LSM2D_LAPLACIAN_ORDER2( 
  LSMLIB_REAL *laplacian_phi,
  const int *ilo_laplacian_phi_gb,
  const int *ihi_laplacian_phi_gb,
  const int *jlo_laplacian_phi_gb,
  const int *jhi_laplacian_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);


/*!
 * LSM2D_PHI_UPWIND_GRAD_F() computes the \f$ \phi \f$-upwind gradient of a
 * function, F, using the following "upwinding" scheme to compute
 * the normal:
 *
 *   if \f$ \phi > 0 \f$:  upwind direction is direction where 
 *                         \f$ \phi \f$ is smaller
 *
 *   if \f$ \phi < 0 \f$:  upwind direction is direction where 
 *                         \f$ \phi \f$ is larger
 *
 * Arguments:
 *  - F_* (out):        components of \f$ \phi \f$-upwinded \f$ \nabla F \f$
 *  - F_*_plus (in):    components of \f$ \nabla F \f$ in plus direction
 *  - F_*_minus (in):   components of \f$ \nabla F \f$ in minus direction
 *  - phi (in):         level set function
 *  - *_gb (in):        index range for ghostbox
 *  - *_fb (in):        index range for fillbox
 *
 * Return value:        none
 *
 * NOTES:
 *  - \f$ \phi \f$ is REQUIRED to have at least one ghost cell in each
 *    coordinate direction for upwinding
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 *
 */
void LSM2D_PHI_UPWIND_GRAD_F(
  LSMLIB_REAL *F_x,
  LSMLIB_REAL *F_y,
  const int *ilo_grad_F_gb,
  const int *ihi_grad_F_gb,
  const int *jlo_grad_F_gb,
  const int *jhi_grad_F_gb,
  LSMLIB_REAL *F_x_plus,
  LSMLIB_REAL *F_y_plus,
  const int *ilo_grad_F_plus_gb,
  const int *ihi_grad_F_plus_gb,
  const int *jlo_grad_F_plus_gb,
  const int *jhi_grad_F_plus_gb,
  LSMLIB_REAL *F_x_minus,
  LSMLIB_REAL *F_y_minus,
  const int *ilo_grad_F_minus_gb,
  const int *ihi_grad_F_minus_gb,
  const int *jlo_grad_F_minus_gb,
  const int *jhi_grad_F_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb);


/*!
 * LSM2D_AVERAGE_GRAD_PHI() computes the average of the plus and minus
 * derivatives:
 *
 * \f[
 *
 *   \nabla \phi = (\nabla \phi_{plus} + \nabla \phi_{minus}) / 2
 *
 * \f]
 *
 * Arguments:
 *  - phi_* (out):       components of average \f$ \nabla \phi \f$
 *  - phi_*_plus (in):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (in):  components of \f$ \nabla \phi \f$ in minus direction
 *  - *_gb (in):         index range for ghostbox
 *  - *_fb (in):         index range for fillbox
 *
 * Return value:         none
 */
void LSM2D_AVERAGE_GRAD_PHI(
  LSMLIB_REAL *phi_x,
  LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  LSMLIB_REAL *phi_x_plus,
  LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb);

/*!
 *
 *  LSM2D_GRADIENT_MAGNITUDE() computes magnitude of the gradient of phi.
 *
 *  Arguments:
 *    phi_* (in):          components of grad(phi)
 *    grad_phi_mag (out):  gradient magnitude
 *    *_gb (in):           index range for ghostbox
 *    *_fb (in):           index range for fillbox
 * 
 */
void LSM2D_GRADIENT_MAGNITUDE(
  const LSMLIB_REAL *phi_x,
  const LSMLIB_REAL *phi_y,
  LSMLIB_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb);
 
/*! 
 *  LSM2D_DIVERGENCE_CENTRAL() computes the second-order, central,  
 *  finite difference approximation to the divergence of a vector field.
 *
 *  Arguments:
 *    divF* (out):  divergence of F
 *    FX, FY(in):   x and y components of vector field F
 *    dx, dy (in):  grid spacing
 *    *_gb (in):    index range for ghostbox
 *    *_fb (in):    index range for fillbox
 * 
 */
void  LSM2D_DIVERGENCE_CENTRAL(
  LSMLIB_REAL *divF,
  const int *ilo_divf_gb, 
  const int *ihi_divf_gb,
  const int *jlo_divf_gb,
  const int *jhi_divf_gb,
  const LSMLIB_REAL *FX,
  const LSMLIB_REAL *FY,
  const int *ilo_gb, 
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb,  
  const int *jhi_fb,
  const LSMLIB_REAL *dx,
  const LSMLIB_REAL *dy);

#ifdef __cplusplus
}
#endif

#endif
