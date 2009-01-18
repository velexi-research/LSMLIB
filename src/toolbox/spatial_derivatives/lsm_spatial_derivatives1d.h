/*
 * File:        lsm_spatial_derivatives1d.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.17 $
 * Modified:    $Date$
 * Description: Header file for Fortran 77 1D ENO/WENO routines.
 */

#ifndef INCLUDED_LSM_SPATIAL_DERIVATIVES_1D_H
#define INCLUDED_LSM_SPATIAL_DERIVATIVES_1D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_spatial_derivatives1d.h
 *
 * \brief
 * @ref lsm_spatial_derivatives1d.h provides support for computing spatial
 * derivatives in one space dimension using high-order ENO and WENO 
 * discretizations.  
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define LSM1D_HJ_ENO1                lsm1dhjeno1_
#define LSM1D_HJ_ENO2                lsm1dhjeno2_
#define LSM1D_HJ_ENO3                lsm1dhjeno3_
#define LSM1D_HJ_WENO5               lsm1dhjweno5_
#define LSM1D_UPWIND_HJ_ENO1         lsm1dupwindhjeno1_
#define LSM1D_UPWIND_HJ_ENO2         lsm1dupwindhjeno2_
#define LSM1D_UPWIND_HJ_ENO3         lsm1dupwindhjeno3_
#define LSM1D_UPWIND_HJ_WENO5        lsm1dupwindhjweno5_
#define LSM1D_CENTRAL_GRAD_ORDER2    lsm1dcentralgradorder2_
#define LSM1D_CENTRAL_GRAD_ORDER4    lsm1dcentralgradorder4_
#define LSM1D_LAPLACIAN_ORDER2       lsm1dlaplacianorder2_
#define LSM1D_PHI_UPWIND_GRAD_F      lsm1dphiupwindgradf_
#define LSM1D_AVERAGE_GRAD_PHI       lsm1daveragegradphi_


/*!
 * LSM1D_HJ_ENO1() computes the forward (plus) and backward (minus)
 * first-order Hamilton-Jacobi ENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - dx (in):            grid spacing
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
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void LSM1D_HJ_ENO1(
  LSMLIB_REAL *phi_x_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  LSMLIB_REAL *D1_x,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*!
 * LSM1D_HJ_ENO2() computes the forward (plus) and backward (minus)
 * second-order Hamilton-Jacobi ENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - D2 (in):            scratch space for holding undivided second-differences
 *  - dx (in):            grid spacing
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
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void LSM1D_HJ_ENO2(
  LSMLIB_REAL *phi_x_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*!
 * LSM1D_HJ_ENO3() computes the forward (plus) and backward (minus)
 * third-order Hamilton-Jacobi ENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - D2 (in):            scratch space for holding undivided second-differences
 *  - D3 (in):            scratch space for holding undivided third-differences
 *  - dx (in):            grid spacing
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
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void LSM1D_HJ_ENO3(
  LSMLIB_REAL *phi_x_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  LSMLIB_REAL *D3,
  const int *ilo_D3_gb,
  const int *ihi_D3_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*!
 * LSM1D_HJ_WENO5() computes the forward (plus) and backward (minus)
 * fifth-order Hamilton-Jacobi WENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - dx (in):            grid spacing
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
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void LSM1D_HJ_WENO5(
  LSMLIB_REAL *phi_x_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*! 
 * LSM1D_UPWIND_HJ_ENO1() computes the first-order Hamilton-Jacobi ENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_x (out):  derivative of \f$ \phi \f$
 *  - phi (in):     \f$ \phi \f$
 *  - vel_x (in):   velocity in the x-direction
 *  - D1 (in):      scratch space for holding undivided first-differences
 *  - dx (in):      grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:    none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void LSM1D_UPWIND_HJ_ENO1(
  LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const LSMLIB_REAL *vel_x,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*! 
 * LSM1D_UPWIND_HJ_ENO2() computes the second-order Hamilton-Jacobi ENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_x (out):  derivative of \f$ \phi \f$
 *  - phi (in):     \f$ \phi \f$
 *  - vel_x (in):   velocity in the x-direction
 *  - D1 (in):      scratch space for holding undivided first-differences
 *  - D2 (in):      scratch space for holding undivided second-differences
 *  - dx (in):      grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:    none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void LSM1D_UPWIND_HJ_ENO2(
  LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const LSMLIB_REAL *vel_x,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*! 
 * LSM1D_UPWIND_HJ_ENO3() computes the third-order Hamilton-Jacobi ENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_x (out):  derivative of \f$ \phi \f$
 *  - phi (in):     \f$ \phi \f$
 *  - vel_x (in):   velocity in the x-direction
 *  - D1 (in):      scratch space for holding undivided first-differences
 *  - D2 (in):      scratch space for holding undivided second-differences
 *  - D3 (in):      scratch space for holding undivided third-differences
 *  - dx (in):      grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:    none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void LSM1D_UPWIND_HJ_ENO3(
  LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const LSMLIB_REAL *vel_x,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  LSMLIB_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  LSMLIB_REAL *D3,
  const int *ilo_D3_gb,
  const int *ihi_D3_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*! 
 * LSM1D_UPWIND_HJ_WENO5() computes the fifth-order Hamilton-Jacobi WENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_x (out):  derivative of \f$ \phi \f$
 *  - phi (in):     \f$ \phi \f$
 *  - vel_x (in):   velocity in the x-direction
 *  - D1 (in):      scratch space for holding undivided first-differences
 *  - dx (in):      grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:    none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void LSM1D_UPWIND_HJ_WENO5(
  LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const LSMLIB_REAL *vel_x,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  LSMLIB_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*! 
 * LSM1D_CENTRAL_GRAD_ORDER2() computes the second-order, central,
 * finite difference approximation to the gradient of \f$ \phi \f$ using 
 * the formula:
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
 *  - dx (in):      grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:        none
 */
void LSM1D_CENTRAL_GRAD_ORDER2( 
  LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*! 
 * LSM1D_CENTRAL_GRAD_ORDER4() computes the fourth-order, central,
 * finite difference approximation to the gradient of \f$ \phi \f$ using 
 * the formula:
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
 *  - dx (in):      grid cell size
 *  - *_gb (in):    index range for ghostbox
 *  - *_fb (in):    index range for fillbox
 *
 * Return value:        none
 */
void LSM1D_CENTRAL_GRAD_ORDER4( 
  LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*! 
 * LSM1D_LAPLACIAN_ORDER2() computes the second-order, central, 
 * finite difference approximation to the Laplacian of \f$ \phi \f$ 
 * using the formula:
 *
 *    \f[
 *
 *      \nabla^2 \phi \approx
 *         \frac{ \phi_{i+1} - 2 \phi_{i} + \phi_{i-1} }
 *              { dx^2 }
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
void LSM1D_LAPLACIAN_ORDER2( 
  LSMLIB_REAL *laplacian_phi,
  const int *ilo_laplacian_phi_gb,
  const int *ihi_laplacian_phi_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const LSMLIB_REAL *dx);


/*!
 * LSM1D_PHI_UPWIND_GRAD_F() computes the \f$ \phi \f$-upwind gradient of a
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
void LSM1D_PHI_UPWIND_GRAD_F(
  LSMLIB_REAL *F_x,
  const int *ilo_grad_F_gb,
  const int *ihi_grad_F_gb,
  LSMLIB_REAL *F_x_plus,
  const int *ilo_grad_F_plus_gb,
  const int *ihi_grad_F_plus_gb,
  LSMLIB_REAL *F_x_minus,
  const int *ilo_grad_F_minus_gb,
  const int *ihi_grad_F_minus_gb,
  const LSMLIB_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb);


/*!
 * LSM1D_AVERAGE_GRAD_PHI() computes the average of the plus and minus
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
void LSM1D_AVERAGE_GRAD_PHI(
  LSMLIB_REAL *phi_x,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  LSMLIB_REAL *phi_x_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  LSMLIB_REAL *phi_x_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *ilo_fb,
  const int *ihi_fb);

#ifdef __cplusplus
}
#endif

#endif
