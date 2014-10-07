/*
 * File:        calculus_toolbox_demo.c
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Demo code for delta and heaviside functions
 */

#include <stdio.h>
#include "LSMLIB_config.h"
#include "lsm_calculus_toolbox.h"

int main(void)
{

  LSMLIB_REAL h[20], d[20];
  LSMLIB_REAL x[20];
  LSMLIB_REAL d_phi;
  LSMLIB_REAL eps;
  int i;

  /* initialize d_phi and eps (smoothing width) */
  d_phi = 0.1;
  eps = 3.0*d_phi;

  /* initialize x */
  for (i = 0; i < 20; i++) x[i] = -1.0 + i*d_phi;

  /* compute smoothed heaviside and delta functions */
  for (i = 0; i < 20; i++) {
    h[i] = LSM_HEAVISIDE(x[i],eps);
    d[i] = LSM_DELTA_FUNCTION(x[i],eps);
  }

  /* output smoothed heaviside function */
  printf("Smoothed Heaviside Function:\n");
  for (i = 0; i < 20; i++) 
    printf("%f, ",h[i]);
  printf("\n");

  /* output smoothed delta-function */
  printf("Smoothed Delta-Function:\n");
  for (i = 0; i < 20; i++) 
    printf("%f, ",d[i]);
  printf("\n");

  return 0;  
}
