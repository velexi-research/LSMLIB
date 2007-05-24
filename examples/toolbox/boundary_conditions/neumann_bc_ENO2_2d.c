/*
 * File:        neumann_bc_ENO2_2d.c
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date: 2006/10/01 17:45:06 $
 * Description: test program for the 2D homogeneous Neumann BCS using
 *              ENO2 discretization
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lsm_boundary_conditions2d.h"
#include "lsm_spatial_derivatives2d.h"

/************************************************************************
 * Helper function declarations
 ************************************************************************/
#define TEST_DIM (2)


/************************************************************************
 *                                                                    
 * Test program for the 2D homogeneous Neumann BCS using ENO2 
 * discretization
 *                                                                  
 ************************************************************************/

int main( int argc, char *argv[])
{
  double *phi;
  double *phi_x_plus, *phi_y_plus;
  double *phi_x_minus, *phi_y_minus;
  double *D1, *D2;
  int ghostcell_width = 2;
  int box_lower[TEST_DIM];
  int box_upper[TEST_DIM];
  int box_dims[TEST_DIM];
  int ghostbox_lower[TEST_DIM];
  int ghostbox_upper[TEST_DIM];
  int ghostbox_dims[TEST_DIM];
  int i,j;
  double dx, dy;
  double err_x_lower;
  double err_x_upper;
  double err_y_lower;
  double err_y_upper;
  int bdry_location_idx;

  /* set index space extents */
  box_dims[0] = 10;
  box_dims[1] = 20;
  box_lower[0] = 0; box_lower[1] = 0;
  box_upper[0] = box_dims[0]-1; box_upper[1] = box_dims[1]-1;
  ghostbox_lower[0] = box_lower[0] - ghostcell_width;
  ghostbox_lower[1] = box_lower[1] - ghostcell_width;
  ghostbox_upper[0] = box_upper[0] + ghostcell_width;
  ghostbox_upper[1] = box_upper[1] + ghostcell_width;
  ghostbox_dims[0] = ghostbox_upper[0] - ghostbox_lower[0] + 1;
  ghostbox_dims[1] = ghostbox_upper[1] - ghostbox_lower[1] + 1;

  /* set grid spacing */
  dx = 1.0/box_dims[0];
  dy = 1.0/box_dims[1];

  /* allocate space for data */
  phi         = malloc(sizeof(double)*ghostbox_dims[0]*ghostbox_dims[1]);
  phi_x_plus  = malloc(sizeof(double)*ghostbox_dims[0]*ghostbox_dims[1]);
  phi_y_plus  = malloc(sizeof(double)*ghostbox_dims[0]*ghostbox_dims[1]);
  phi_x_minus = malloc(sizeof(double)*ghostbox_dims[0]*ghostbox_dims[1]);
  phi_y_minus = malloc(sizeof(double)*ghostbox_dims[0]*ghostbox_dims[1]);
  D1          = malloc(sizeof(double)*ghostbox_dims[0]*ghostbox_dims[1]);
  D2          = malloc(sizeof(double)*ghostbox_dims[0]*ghostbox_dims[1]);

  /* set phi data in interior */
  for (j = 0; j < box_dims[1]; j++) {
    for (i = 0; i < box_dims[0]; i++) {
      int idx = (i+ghostcell_width) 
              + (j+ghostcell_width)*ghostbox_dims[0];
      double x = (i+0.5)*dx;
      double y = (j+0.5)*dy;

      phi[idx] = (x-0.15)*(x-0.15) + 4*(y-0.55)*(y-0.55) - 0.25;

    }
  }
 
  /* impose homogeneous Neumann boundary conditions */
  bdry_location_idx = 0;
  LSM2D_HOMOGENEOUS_NEUMANN_ENO2(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &bdry_location_idx);

  bdry_location_idx = 1;
  LSM2D_HOMOGENEOUS_NEUMANN_ENO2(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &bdry_location_idx);

  bdry_location_idx = 2;
  LSM2D_HOMOGENEOUS_NEUMANN_ENO2(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &bdry_location_idx);

  bdry_location_idx = 3;
  LSM2D_HOMOGENEOUS_NEUMANN_ENO2(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &bdry_location_idx);


  /* compute spatial derivatives using ENO2 */
  LSM2D_HJ_ENO2(
    phi_x_plus, phi_y_plus,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1],
    phi_x_minus, phi_y_minus,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1],
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    D1,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1],
    D2,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1],
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &dx, &dy);

  /* 
   * check that derivatives normal to faces of computational domain 
   * are 0.
   */

  /* x-lower */
  i = box_lower[0];
  err_x_lower = 0.0;
  for (j = 0; j < box_dims[1]; j++) {
    int idx = (i+ghostcell_width) 
            + (j+ghostcell_width)*ghostbox_dims[0];

    if (fabs(phi_x_minus[idx]) > err_x_lower) 
      err_x_lower = fabs(phi_x_minus[idx]);
  }
  
  /* x-upper*/
  i = box_upper[0];
  err_x_upper = 0.0;
  for (j = 0; j < box_dims[1]; j++) {
    int idx = (i+ghostcell_width) 
            + (j+ghostcell_width)*ghostbox_dims[0];

    if (fabs(phi_x_plus[idx]) > err_x_upper) 
      err_x_upper = fabs(phi_x_plus[idx]);
  }
  
  /* y-lower */
  j = box_lower[1];
  err_y_lower = 0.0;
  for (i = 0; i < box_dims[0]; i++) {
    int idx = (i+ghostcell_width) 
            + (j+ghostcell_width)*ghostbox_dims[0];

    if (fabs(phi_y_minus[idx]) > err_y_lower) 
      err_y_lower = fabs(phi_y_minus[idx]);
  }
  
  /* y-upper*/
  j = box_upper[1];
  err_y_upper = 0.0;
  for (i = 0; i < box_dims[0]; i++) {
    int idx = (i+ghostcell_width) 
            + (j+ghostcell_width)*ghostbox_dims[0];

    if (fabs(phi_y_plus[idx]) > err_y_upper) 
      err_y_upper = fabs(phi_y_plus[idx]);
  }

  /* print errors */
  printf("Error in derivative x-lower:%.15f\n", err_x_lower);
  printf("Error in derivative x-upper:%.15f\n", err_x_upper);
  printf("Error in derivative y-lower:%.15f\n", err_y_lower);
  printf("Error in derivative y-upper:%.15f\n", err_y_upper);

  /* print phi values */
/*
  printf("print phi values:\n");
  for (j = 0; j < ghostbox_dims[1]; j++) {
    int idx;
    for (i = 0; i < ghostbox_dims[0]-1; i++) {
      idx = i + j*ghostbox_dims[0];
      printf("%f,",phi[idx]);
    }
    idx = i + j*ghostbox_dims[0];
    printf("%f\n",phi[idx]);
  }
*/

  return 0;
}

