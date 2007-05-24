/*
 * File:        UPWIND_HJ_ENO1_1D.c
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.13 $
 * Modified:    $Date: 2006/04/24 11:48:36 $
 * Description: MATLAB MEX-file for 1d, first-order upwind HJ ENO 
 */

/*=======================================================================
 *
 * UPWIND_HJ_ENO1_1D() computes the first-order upwind HJ ENO
 * approximation to phi_x.
 *
 * Usage: phi_x = UPWIND_HJ_ENO1_1D(phi, vel_x, ghostcell_width, dx)
 *
 * Arguments:
 * - phi:               function for which to compute upwind 
 *                        derivative
 * - vel_x:             velocity to use in upwinding
 * - ghostcell_width:   number of ghostcells at boundary of 
 *                        computational domain
 * - dx:                grid cell size
 *
 * Return values:
 * - phi_x:             first-order, upwind HJ ENO derivative
 *
 * NOTES:
 * - phi_x has the same ghostcell width as phi.
 *
 *=======================================================================*/

#include "mex.h"
#include "matrix.h"
#include "lsm_spatial_derivatives1d.h"

/* Input Arguments */ 
#define PHI             (prhs[0])
#define VEL_X           (prhs[1])
#define GHOSTCELL_WIDTH (prhs[2])
#define DX              (prhs[3])

/* Output Arguments */ 
#define PHI_X           (plhs[0])


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  double *phi_x;
  int ilo_grad_phi_gb, ihi_grad_phi_gb;
  double *phi; 
  int ilo_phi_gb, ihi_phi_gb;
  double *vel_x; 
  int ilo_vel_gb, ihi_vel_gb;
  double *D1; 
  int ilo_D1_gb, ihi_D1_gb;
  int ilo_fb, ihi_fb;
  double dx;
  int ghostcell_width;
  int dim_M, dim_N;
  
  /* Check for proper number of arguments */
  if (nrhs != 4) { 
    mexErrMsgTxt("Four required input arguments."); 
  } else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments."); 
  } 
    
  /* Parameter Checks */
  dim_M = mxGetM(PHI);
  dim_N = mxGetN(PHI);
  if ( (dim_M > 1) && (dim_N > 1) ) {
    mexErrMsgTxt("phi should be a 1 dimensional array."); 
  }
  dim_M = mxGetM(VEL_X);
  dim_N = mxGetN(VEL_X);
  if ( (dim_M > 1) && (dim_N > 1) ) {
    mexErrMsgTxt("vel_x should be a 1 dimensional array."); 
  }

  /* Get ghostcell_width */
  ghostcell_width = mxGetPr(GHOSTCELL_WIDTH)[0];

  /* Get dx */ 
  dx = mxGetPr(DX)[0];

  /* Assign pointers for phi and vel_x */
  phi = mxGetPr(PHI);
  vel_x = mxGetPr(VEL_X);
      
  /* Get length of phi data */
  ilo_phi_gb = 1;
  ihi_phi_gb = mxGetM(PHI);
  if (1 == ihi_phi_gb) {
    /* phi is row matrix */
    ihi_phi_gb = mxGetN(PHI);
  }

  /* Get length of vel_x data */
  ilo_vel_gb = 1;
  ihi_vel_gb = mxGetM(VEL_X);
  if (1 == ihi_vel_gb) {
    /* vel_x is row matrix */
    ihi_vel_gb = mxGetN(VEL_X);
  }

  /* if necessary, shift ghostbox for velocity to be */
  /* centered with respect to the ghostbox for phi.  */
  if (ihi_vel_gb != ihi_phi_gb) {
    int shift = (ihi_phi_gb-ihi_vel_gb)/2;
    ilo_vel_gb += shift;
    ihi_vel_gb += shift;
  }

  /* Create matrices for upwind derivatives (i.e. phi_x) */
  ilo_grad_phi_gb = ilo_phi_gb;
  ihi_grad_phi_gb = ihi_phi_gb;
  PHI_X = mxCreateDoubleMatrix(ihi_grad_phi_gb-ilo_grad_phi_gb+1,1,mxREAL);
  phi_x = mxGetPr(PHI_X); 

  /* Allocate scratch memory for undivided differences */
  ilo_D1_gb = ilo_phi_gb;
  ihi_D1_gb = ihi_phi_gb;
  D1 = (double*) malloc( sizeof(double) * (ihi_D1_gb-ilo_D1_gb+1) );
  if (!D1) {
    mexErrMsgTxt("Unable to allocate memory for scratch data...aborting....");
  }

  /* Do the actual computations in a Fortran 77 subroutine */
  ilo_fb = ilo_phi_gb+ghostcell_width;
  ihi_fb = ihi_phi_gb-ghostcell_width;
  LSM1D_UPWIND_HJ_ENO1(
    phi_x, &ilo_grad_phi_gb, &ihi_grad_phi_gb,
    phi, &ilo_phi_gb, &ihi_phi_gb, 
    vel_x, &ilo_vel_gb, &ihi_vel_gb,
    D1, &ilo_D1_gb, &ihi_D1_gb, 
    &ilo_fb, &ihi_fb,
    &dx);

  /* Deallocate scratch memory for undivided differences */
  free(D1);

  return;
}


