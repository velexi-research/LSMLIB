/*
 * File:        HJ_ENO3_1D.c
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.6 $
 * Modified:    $Date: 2006/04/24 00:52:38 $
 * Description: MATLAB MEX-file for 1d, third-order plus and minus HJ ENO 
 */

/*=======================================================================
 *
 * HJ_ENO3_1D() computes the third-order plus and minus HJ ENO
 * approximation to phi_x.
 *
 * Usage: [phi_x_plus, phi_x_minus] = HJ_ENO3_1D(phi, ghostcell_width, dx)
 *
 * Arguments:
 * - phi:               function for which to compute plus and minus
 *                        spatial derivatives
 * - ghostcell_width:   number of ghostcells at boundary of 
 *                        computational domain
 * - dx:                grid cell size
 *
 * Return values:
 * - phi_x_plus:        third-order, plus HJ ENO derivative
 * - phi_x_minus:       third-order, minus HJ ENO derivative
 *
 * NOTES:
 * - phi_x_plus and phi_x_minus have the same ghostcell width as phi.
 *
 *=======================================================================*/

#include "mex.h"
#include "matrix.h"
#include "lsm_spatial_derivatives1d.h"

/* Input Arguments */ 
#define PHI             (prhs[0])
#define GHOSTCELL_WIDTH (prhs[1])
#define DX              (prhs[2])

/* Output Arguments */ 
#define PHI_X_PLUS      (plhs[0])
#define PHI_X_MINUS     (plhs[1])


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  double *phi_x_plus, *phi_x_minus;
  int ilo_grad_phi_gb, ihi_grad_phi_gb;
  double *phi; 
  int ilo_phi_gb, ihi_phi_gb;
  double *D1; 
  int ilo_D1_gb, ihi_D1_gb;
  double *D2; 
  int ilo_D2_gb, ihi_D2_gb;
  double *D3; 
  int ilo_D3_gb, ihi_D3_gb;
  int ilo_fb, ihi_fb;
  double dx;
  int ghostcell_width;
  int dim_M, dim_N;
  
  /* Check for proper number of arguments */
  if (nrhs != 3) { 
    mexErrMsgTxt("Three required input arguments."); 
  } else if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments."); 
  } 
    
  /* Parameter Checks */
  dim_M = mxGetM(PHI);
  dim_N = mxGetN(PHI);
  if ( (dim_M > 1) && (dim_N > 1) ) {
    mexErrMsgTxt("phi should be a 1 dimensional array."); 
  }

  /* Get ghostcell_width */
  ghostcell_width = mxGetPr(GHOSTCELL_WIDTH)[0];

  /* Get dx */ 
  dx = mxGetPr(DX)[0];

  /* Assign pointers for phi */
  phi = mxGetPr(PHI);
      
  /* Get length of phi data */
  ilo_phi_gb = 1;
  ihi_phi_gb = mxGetM(PHI);
  if (1 == ihi_phi_gb) {
    /* phi is row matrix */
    ihi_phi_gb = mxGetN(PHI);
  }

  /* Create matrices for plus and minus derivatives */
  ilo_grad_phi_gb = ilo_phi_gb;
  ihi_grad_phi_gb = ihi_phi_gb;
  PHI_X_PLUS  = mxCreateDoubleMatrix(
    ihi_grad_phi_gb-ilo_grad_phi_gb+1,1,mxREAL);
  PHI_X_MINUS = mxCreateDoubleMatrix(
    ihi_grad_phi_gb-ilo_grad_phi_gb+1,1,mxREAL);
  phi_x_plus  = mxGetPr(PHI_X_PLUS); 
  phi_x_minus = mxGetPr(PHI_X_MINUS); 

  /* Allocate scratch memory for undivided differences */
  ilo_D1_gb = ilo_phi_gb;
  ihi_D1_gb = ihi_phi_gb;
  D1 = (double*) malloc( sizeof(double) * (ihi_D1_gb-ilo_D1_gb+1) );

  ilo_D2_gb = ilo_phi_gb;
  ihi_D2_gb = ihi_phi_gb;
  D2 = (double*) malloc( sizeof(double) * (ihi_D2_gb-ilo_D2_gb+1) );

  ilo_D3_gb = ilo_phi_gb;
  ihi_D3_gb = ihi_phi_gb;
  D3 = (double*) malloc( sizeof(double) * (ihi_D3_gb-ilo_D3_gb+1) );

  if ( (!D1) || (!D2) || (!D3) ) {
    if (D1) free(D1);
    if (D2) free(D2);
    if (D3) free(D3);
    mexErrMsgTxt("Unable to allocate memory for scratch data...aborting....");
  }

  /* Do the actual computations in a Fortran 77 subroutine */
  ilo_fb = ilo_phi_gb+ghostcell_width;
  ihi_fb = ihi_phi_gb-ghostcell_width;
  LSM1D_HJ_ENO3(
    phi_x_plus, &ilo_grad_phi_gb, &ihi_grad_phi_gb,
    phi_x_minus, &ilo_grad_phi_gb, &ihi_grad_phi_gb,
    phi, &ilo_phi_gb, &ihi_phi_gb, 
    D1, &ilo_D1_gb, &ihi_D1_gb, 
    D2, &ilo_D2_gb, &ihi_D2_gb, 
    D3, &ilo_D3_gb, &ihi_D3_gb, 
    &ilo_fb, &ihi_fb,
    &dx);

  /* Deallocate scratch memory for undivided differences */
  free(D1);
  free(D2);
  free(D3);

  return;
}


