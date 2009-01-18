/*
 * File:        neumann_bc_ENO1_3d.c
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Demo program for the 2D homogeneous Neumann BCS using
 *              ENO1 discretization
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "LSMLIB_config.h"
#include "lsm_boundary_conditions3d.h"
#include "lsm_spatial_derivatives3d.h"

/************************************************************************
 * Helper function declarations
 ************************************************************************/
#define DIM (3)


/************************************************************************
 *                                                                    
 * Demo program for the 3D homogeneous Neumann BCS using ENO1 
 * discretization
 *                                                                  
 ************************************************************************/

int main( int argc, char *argv[])
{
  LSMLIB_REAL *phi;
  LSMLIB_REAL *phi_x_plus, *phi_y_plus, *phi_z_plus;
  LSMLIB_REAL *phi_x_minus, *phi_y_minus, *phi_z_minus;
  LSMLIB_REAL *D1;
  int ghostcell_width = 1;
  int box_lower[DIM];
  int box_upper[DIM];
  int box_dims[DIM];
  int ghostbox_lower[DIM];
  int ghostbox_upper[DIM];
  int ghostbox_dims[DIM];
  int grid_size;
  int i,j,k;
  LSMLIB_REAL dx, dy, dz;
  LSMLIB_REAL err_x_lower;
  LSMLIB_REAL err_x_upper;
  LSMLIB_REAL err_y_lower;
  LSMLIB_REAL err_y_upper;
  LSMLIB_REAL err_z_lower;
  LSMLIB_REAL err_z_upper;
  int bdry_location_idx;

  /* set index space extents */
  box_dims[0] = 10;
  box_dims[1] = 20;
  box_dims[2] = 30;
  box_lower[0] = 0; box_lower[1] = 0; box_lower[2] = 0;
  box_upper[0] = box_dims[0]-1; 
  box_upper[1] = box_dims[1]-1;
  box_upper[2] = box_dims[2]-1;
  ghostbox_lower[0] = box_lower[0] - ghostcell_width;
  ghostbox_lower[1] = box_lower[1] - ghostcell_width;
  ghostbox_lower[2] = box_lower[2] - ghostcell_width;
  ghostbox_upper[0] = box_upper[0] + ghostcell_width;
  ghostbox_upper[1] = box_upper[1] + ghostcell_width;
  ghostbox_upper[2] = box_upper[2] + ghostcell_width;
  ghostbox_dims[0] = ghostbox_upper[0] - ghostbox_lower[0] + 1;
  ghostbox_dims[1] = ghostbox_upper[1] - ghostbox_lower[1] + 1;
  ghostbox_dims[2] = ghostbox_upper[2] - ghostbox_lower[2] + 1;

  /* set grid spacing */
  dx = 1.0/box_dims[0];
  dy = 1.0/box_dims[1];
  dz = 1.0/box_dims[2];

  /* allocate space for data */
  grid_size = ghostbox_dims[0]*ghostbox_dims[1]*ghostbox_dims[2];
  phi         = malloc(sizeof(LSMLIB_REAL)*grid_size);
  phi_x_plus  = malloc(sizeof(LSMLIB_REAL)*grid_size);
  phi_y_plus  = malloc(sizeof(LSMLIB_REAL)*grid_size);
  phi_z_plus  = malloc(sizeof(LSMLIB_REAL)*grid_size);
  phi_x_minus = malloc(sizeof(LSMLIB_REAL)*grid_size);
  phi_y_minus = malloc(sizeof(LSMLIB_REAL)*grid_size);
  phi_z_minus = malloc(sizeof(LSMLIB_REAL)*grid_size);
  D1          = malloc(sizeof(LSMLIB_REAL)*grid_size);

  /* set phi data in interior */
  for (k = 0; k < box_dims[2]; k++) {
    for (j = 0; j < box_dims[1]; j++) {
      for (i = 0; i < box_dims[0]; i++) {
        int idx = (i+ghostcell_width) 
                + (j+ghostcell_width)*ghostbox_dims[0]
                + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];
        LSMLIB_REAL x = (i+0.5)*dx;
        LSMLIB_REAL y = (j+0.5)*dy;
        LSMLIB_REAL z = (k+0.5)*dz;

        phi[idx] = (x-0.25)*(x-0.25) + 4*(y-0.35)*(y-0.35) 
                 + 0.25*(z-0.45)*(z-0.45) - 0.25;
      }
    }
  }
 
  /* impose homogeneous Neumann boundary conditions */
  bdry_location_idx = 0;
  LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &ghostbox_lower[2], 
    &ghostbox_upper[2], 
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &box_lower[2], 
    &box_upper[2],
    &bdry_location_idx);

  bdry_location_idx = 1;
  LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &ghostbox_lower[2], 
    &ghostbox_upper[2], 
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &box_lower[2], 
    &box_upper[2],
    &bdry_location_idx);

  bdry_location_idx = 2;
  LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &ghostbox_lower[2], 
    &ghostbox_upper[2], 
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &box_lower[2], 
    &box_upper[2],
    &bdry_location_idx);

  bdry_location_idx = 3;
  LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &ghostbox_lower[2], 
    &ghostbox_upper[2], 
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &box_lower[2], 
    &box_upper[2],
    &bdry_location_idx);

  bdry_location_idx = 4;
  LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &ghostbox_lower[2], 
    &ghostbox_upper[2], 
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &box_lower[2], 
    &box_upper[2],
    &bdry_location_idx);

  bdry_location_idx = 5;
  LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &ghostbox_lower[2], 
    &ghostbox_upper[2], 
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &box_lower[2], 
    &box_upper[2],
    &bdry_location_idx);


  /* compute spatial derivatives using ENO1 */
  LSM3D_HJ_ENO1(
    phi_x_plus, phi_y_plus, phi_z_plus,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1],
    &ghostbox_lower[2], 
    &ghostbox_upper[2],
    phi_x_minus, phi_y_minus, phi_z_minus,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1],
    &ghostbox_lower[2], 
    &ghostbox_upper[2],
    phi,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1], 
    &ghostbox_lower[2], 
    &ghostbox_upper[2],
    D1,
    &ghostbox_lower[0], 
    &ghostbox_upper[0], 
    &ghostbox_lower[1], 
    &ghostbox_upper[1],
    &ghostbox_lower[2], 
    &ghostbox_upper[2],
    &box_lower[0], 
    &box_upper[0], 
    &box_lower[1], 
    &box_upper[1],
    &box_lower[2], 
    &box_upper[2],
    &dx, &dy, &dz);


  /* 
   * check that derivatives normal to faces of computational domain 
   * are 0.
   */

  /* x-lower */
  i = box_lower[0];
  err_x_lower = 0.0;
  for (k = 0; k < box_dims[2]; k++) {
    for (j = 0; j < box_dims[1]; j++) {
      int idx = (i+ghostcell_width) 
              + (j+ghostcell_width)*ghostbox_dims[0]
              + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

      if (fabs(phi_x_minus[idx]) > err_x_lower) 
        err_x_lower = fabs(phi_x_minus[idx]);
   }
 }
  
  /* x-upper*/
  i = box_upper[0];
  err_x_upper = 0.0;
  for (k = 0; k < box_dims[2]; k++) {
    for (j = 0; j < box_dims[1]; j++) {
      int idx = (i+ghostcell_width) 
              + (j+ghostcell_width)*ghostbox_dims[0]
              + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

      if (fabs(phi_x_plus[idx]) > err_x_upper) 
        err_x_upper = fabs(phi_x_plus[idx]);
    }
  }
  
  /* y-lower */
  j = box_lower[1];
  err_y_lower = 0.0;
  for (k = 0; k < box_dims[2]; k++) {
    for (i = 0; i < box_dims[0]; i++) {
      int idx = (i+ghostcell_width) 
              + (j+ghostcell_width)*ghostbox_dims[0]
              + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

      if (fabs(phi_y_minus[idx]) > err_y_lower) 
        err_y_lower = fabs(phi_y_minus[idx]);
    }
  }
  
  /* y-upper*/
  j = box_upper[1];
  err_y_upper = 0.0;
  for (k = 0; k < box_dims[2]; k++) {
    for (i = 0; i < box_dims[0]; i++) {
      int idx = (i+ghostcell_width) 
              + (j+ghostcell_width)*ghostbox_dims[0]
              + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

      if (fabs(phi_y_plus[idx]) > err_y_upper) 
        err_y_upper = fabs(phi_y_plus[idx]);
    }
  }

  /* z-lower */
  k = box_lower[2];
  err_z_lower = 0.0;
  for (j = 0; j < box_dims[1]; j++) {
    for (i = 0; i < box_dims[0]; i++) {
      int idx = (i+ghostcell_width) 
              + (j+ghostcell_width)*ghostbox_dims[0]
              + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

      if (fabs(phi_z_minus[idx]) > err_z_lower) 
        err_z_lower = fabs(phi_z_minus[idx]);
   }
 }
  
  /* z-upper*/
  k = box_upper[2];
  err_z_upper = 0.0;
  for (j = 0; j < box_dims[1]; j++) {
    for (i = 0; i < box_dims[0]; i++) {
      int idx = (i+ghostcell_width) 
              + (j+ghostcell_width)*ghostbox_dims[0]
              + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

      if (fabs(phi_z_plus[idx]) > err_z_upper) 
        err_z_upper = fabs(phi_z_plus[idx]);
    }
  }
  

  /* print errors */
  printf("Error in derivative x-lower:%.15f\n", err_x_lower);
  printf("Error in derivative x-upper:%.15f\n", err_x_upper);
  printf("Error in derivative y-lower:%.15f\n", err_y_lower);
  printf("Error in derivative y-upper:%.15f\n", err_y_upper);
  printf("Error in derivative z-lower:%.15f\n", err_z_lower);
  printf("Error in derivative z-upper:%.15f\n", err_z_upper);

  /* print phi values */
/*
  printf("print phi values:\n");
  for (k = 0; k < ghostbox_dims[2]; k++) {
    for (j = 0; j < ghostbox_dims[1]; j++) {
      int idx;
      for (i = 0; i < ghostbox_dims[0]-1; i++) {
        idx = i + j*ghostbox_dims[0] + k*ghostbox_dims[0]*ghostbox_dims[1];
        printf("%f,",phi[idx]);
      }
      idx = i + j*ghostbox_dims[0] + k*ghostbox_dims[0]*ghostbox_dims[1];
      printf("%f\n",phi[idx]);
    }
    printf("\n");
  }
*/

  return 0;
}
