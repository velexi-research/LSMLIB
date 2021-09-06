/*
 * File:        computeDistanceFunction3d.c
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Demo program for the fast marching method functions
 */


/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Header for Fast Marching Method Algorithm class */
#include "lsmlib_config.h"
#include "lsm_fast_marching_method.h"

/************************************************************************
 *
 * Demo program for the Fast Marching Method classes.
 *
 ************************************************************************
 */


int main( int argc, char *argv[])
{
  /* field variables */
  LSMLIB_REAL *phi;
  LSMLIB_REAL *distance_function;
  LSMLIB_REAL *mask = 0;

  /* grid parameters */
  LSMLIB_REAL X_lo[3] = {-1.0,-1.0,-1.0}; 
  LSMLIB_REAL X_hi[3] = {1.0,1.0,1.0}; 
  LSMLIB_REAL dx[3];
  int N;
  int i,j,k;
  int idx;
  int num_gridpts;
  int grid_dims[3];

  /* numerical parameters */
  int spatial_derivative_order = 2;

  /* auxilliary variables */
  LSMLIB_REAL x,y,z;
  LSMLIB_REAL center[3], radius;

  /* file pointer to output results */
  FILE *data_file;

  /* set up grid */
  N = 50;
  num_gridpts = 1;
  for (i = 0; i < 3; i++) {
    dx[i] = (X_hi[i]-X_lo[i])/N;
    grid_dims[i] = N+1;
    num_gridpts *= grid_dims[i];
  }

  /* allocate memory for field data */
  phi = (LSMLIB_REAL*) malloc(num_gridpts*sizeof(LSMLIB_REAL));
  distance_function = (LSMLIB_REAL*) malloc(num_gridpts*sizeof(LSMLIB_REAL));

  /* initialize data */
  center[0] = 0.0; center[1] = 0.0; center[2] = 0.0; 
  radius = 0.5;
  for (k = 0; k < grid_dims[2]; k++) {
    for (j = 0; j < grid_dims[1]; j++) {
      for (i = 0; i < grid_dims[0]; i++) {
        idx = i+j*grid_dims[0]+k*grid_dims[0]*grid_dims[1];
        x = X_lo[0]+dx[0]*i;
        y = X_lo[1]+dx[1]*j;
        z = X_lo[2]+dx[2]*k;

        phi[idx] = ( (x-center[0])*(x-center[0])
                    +(y-center[1])*(y-center[1]) 
                    +(z-center[2])*(z-center[2]) ) - radius*radius;
      }
    }
  }
  for (i = 0; i < num_gridpts; i++) {
    distance_function[i] = 0; 
  }

  /* Carry out FMM calculation */
  computeDistanceFunction3d(
    distance_function,
    phi,
    mask,
    spatial_derivative_order,
    grid_dims,
    dx);

  /* write results to output file */
  data_file = fopen("computeDistanceFunction3d.dat","w");
  for (idx = 0; idx < num_gridpts; idx++) {
    fprintf(data_file,"%f %f\n", distance_function[idx],phi[idx]);
  }
  fclose(data_file);

  /* clean up memory */
  free(phi);  
  free(distance_function);  

  return(0);
}

