/*
 * File:        neumann_bc_ENO1_1d.c
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Demo program for the 1D homogeneous Neumann BCS using
 *              ENO1 discretization
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lsmlib_config.h"
#include "lsm_boundary_conditions1d.h"
#include "lsm_spatial_derivatives1d.h"

/************************************************************************
 * Helper function declarations
 ************************************************************************/
#define DIM (1)


/************************************************************************
 *                                                                    
 * Demo program for the 1D homogeneous Neumann BCS using ENO1 
 * discretization
 *                                                                  
 ************************************************************************/

int main( int argc, char *argv[])
{
  LSMLIB_REAL *phi;
  LSMLIB_REAL *phi_x_plus;
  LSMLIB_REAL *phi_x_minus;
  LSMLIB_REAL *D1;
  int ghostcell_width = 1;
  int box_lower[DIM];
  int box_upper[DIM];
  int box_dims[DIM];
  int ghostbox_lower[DIM];
  int ghostbox_upper[DIM];
  int ghostbox_dims[DIM];
  int i;
  LSMLIB_REAL dx;
  LSMLIB_REAL err_x_lower;
  LSMLIB_REAL err_x_upper;
  int bdry_location_idx;

  /* set index space extents */
  box_dims[0] = 10;
  box_lower[0] = 0;
  box_upper[0] = box_dims[0]-1;
  ghostbox_lower[0] = box_lower[0] - ghostcell_width;
  ghostbox_upper[0] = box_upper[0] + ghostcell_width;
  ghostbox_dims[0] = ghostbox_upper[0] - ghostbox_lower[0] + 1;

  /* set grid spacing */
  dx = 1.0/box_dims[0];

  /* allocate space for data */
  phi         = malloc(sizeof(LSMLIB_REAL)*ghostbox_dims[0]);
  phi_x_plus  = malloc(sizeof(LSMLIB_REAL)*ghostbox_dims[0]);
  phi_x_minus = malloc(sizeof(LSMLIB_REAL)*ghostbox_dims[0]);
  D1          = malloc(sizeof(LSMLIB_REAL)*ghostbox_dims[0]);

  /* set phi data in interior */
  for (i = 0; i < box_dims[0]; i++) {
    int idx = i+ghostcell_width;
    LSMLIB_REAL x = (i+0.5)*dx;
    phi[idx] = (x-0.25)*(x-0.25);
  }
 
  /* impose homogeneous Neumann boundary conditions */
  bdry_location_idx = 0;
  LSM1D_HOMOGENEOUS_NEUMANN_ENO1(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &box_lower[0], 
    &box_upper[0],
    &bdry_location_idx);

  bdry_location_idx = 1;
  LSM1D_HOMOGENEOUS_NEUMANN_ENO1(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &box_lower[0], 
    &box_upper[0],
    &bdry_location_idx);


  /* compute spatial derivatives using ENO1 */
  LSM1D_HJ_ENO1(
    phi_x_plus,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    phi_x_minus,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    D1,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &box_lower[0], 
    &box_upper[0], 
    &dx);

  /* 
   * check that derivatives normal to faces of computational domain 
   * are 0.
   */

  /* x-lower */
  err_x_lower = fabs(phi_x_minus[box_lower[0]+ghostcell_width]);
  
  /* x-upper*/
  err_x_upper = fabs(phi_x_plus[box_upper[0]+ghostcell_width]);

  /* print errors */
  printf("Error in derivative x-lower:%.15f\n", err_x_lower);
  printf("Error in derivative x-upper:%.15f\n", err_x_upper);

  /* print phi values */
  printf("phi values:\n");
  for (i = 0; i < ghostbox_dims[0]-1; i++) {
    printf("%f,",phi[i]);
  }
  printf("%f\n",phi[i]);

  /* print phi_x_minus values */
  printf("phi_x_minus values:\n");
  for (i = 0; i < box_dims[0]-1; i++) {
    printf("%f,",phi_x_minus[i+ghostcell_width]);
  }
  printf("%f\n",phi_x_minus[i+ghostcell_width]);

  /* print phi_x_plus values */
  printf("phi_x_plus values:\n");
  for (i = 0; i < box_dims[0]-1; i++) {
    printf("%f,",phi_x_plus[i+ghostcell_width]);
  }
  printf("%f\n",phi_x_plus[i+ghostcell_width]);

  return 0;
}

