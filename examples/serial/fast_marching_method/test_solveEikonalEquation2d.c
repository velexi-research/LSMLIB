/*
 * File:        test_solveEikonalEquation2d.c
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.2 $
 * Modified:    $Date: 2006/08/13 13:35:24 $
 * Description: Test program for the fast marching method functions
 */


// System headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Header for Fast Marching Method Algorithm class
#include "lsm_fast_marching_method.h"

/************************************************************************
 *
 * Test program for the Fast Marching Method classes.
 *
 ************************************************************************
 */


int main( int argc, char *argv[])
{
  /* field variables */
  double *phi;
  double *speed;
  double *mask = 0;

  /* grid parameters */
  double X_lo[2] = {-1.0,-1.0};
  double X_hi[2] = {1.0,1.0};
  double dx[2];
  int N;
  int i,j;
  int idx;
  int num_gridpts;
  int grid_dims[2];

  /* numerical parameters */
  int spatial_derivative_order = 1;

  /* auxilliary variables */
  double x,y;
  double center[2], radius;

  /* file pointer to output results */
  FILE *data_file;

  /* set up grid */
  N = 100;
  num_gridpts = 1;
  for (i = 0; i < 2; i++) {
    dx[i] = (X_hi[i]-X_lo[i])/N;
    grid_dims[i] = N+1;
    num_gridpts *= grid_dims[i];
  }

  /* allocate memory for field data */
  phi   = (double*) malloc(num_gridpts*sizeof(double));
  speed = (double*) malloc(num_gridpts*sizeof(double));

  /* initialize data */
  center[0] = 0.0; center[1] = 0.0;
  radius = 0.5;
  for (j = 0; j < grid_dims[1]; j++) {
    for (i = 0; i < grid_dims[0]; i++) {
      idx = i+j*grid_dims[0];
      x = X_lo[0]+dx[0]*i;
      y = X_lo[1]+dx[1]*j;

      /* set boundary data for phi */
      if ( ((x-center[0])*(x-center[0])+(y-center[1])*(y-center[1])) 
         - radius*radius < 0 ) {
        phi[idx] = 0.0;
      } else {;
        phi[idx] = -1.0;
      }

      /* set speed function */
      if (x<0) {
        speed[idx] = 1;
      } else {
        speed[idx] = 2;
      }

    }
  }

  /* Carry out FMM calculation */
  solveEikonalEquation2d(
    phi,
    speed,
    mask,
    spatial_derivative_order,
    grid_dims,
    dx);

  /* write results to output file */
  data_file = fopen("test_solveEikonalEquation2d.dat","w");
  for (idx = 0; idx < num_gridpts; idx++) {
    fprintf(data_file,"%f %f\n", phi[idx],speed[idx]);
  }
  fclose(data_file);

  /* clean up memory */
  free(phi);  
  free(speed);  

  return(0);
}

