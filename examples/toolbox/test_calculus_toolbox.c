/*
 * File:        test_calculus_toolbox.c
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.3 $
 * Modified:    $Date: 2006/04/26 13:45:31 $
 * Description: Test code for delta and heaviside functions
 */

#include <stdio.h>
#include "lsm_calculus_toolbox.h"

int main(void)
{

  double h[20], d[20];
  double x[20];
  double d_phi;
  double eps;
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
