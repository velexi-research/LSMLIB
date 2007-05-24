/*
 * File:        lsm_FMM_eikonal_2d.c
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date: 2006/08/13 13:30:52 $
 * Description: Implementation of 2D Fast Marching Method for Eikonal equation
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "FMM_Core.h"
#include "FMM_Heap.h"
#include "lsm_fast_marching_method.h"


/*===================== lsm_FMM_2d Data Structures ====================*/
struct FMM_FieldData {
  double *phi;                 /* solution to Eikonal equation */
  double *speed;               /* speed function               */
};


/*======================= lsm_FMM_2d Constants ==========================*/
#define LSM_FMM_2D_NDIM                   (2)
#define LSM_FMM_2D_TRUE                   (1)
#define LSM_FMM_2D_FALSE                  (0)
#define LSM_FMM_2D_DEFAULT_UPDATE_VALUE   (-1)
#define LSM_FMM_2D_ZERO_TOL               (1e-8)


/*========================== Error Codes ============================*/
#define LSM_FMM_2D_ERR_SUCCESS                             (0)
#define LSM_FMM_2D_ERR_FMM_DATA_CREATION_ERROR             (1)
#define LSM_FMM_2D_ERR_INVALID_SPATIAL_DERIVATIVE_ORDER    (2)


/*========================= lsm_FMM_2d Macros ===========================*/
#define LSM_FMM_2D_IDX(i,j)     ((i) + grid_dims[0]*(j))
#define LSM_FMM_2D_ABS(x)       ((x) > 0 ? (x) : -1.0*(x))
#define LSM_FMM_2D_IDX_OUT_OF_BOUNDS(i,j,grid_dims)            \
  ( ((i)<0) || ((i)>(grid_dims)[0]-1) || ((j)<0) || ((j)>(grid_dims)[1]-1) ) 


/*================== Helper Functions Declarations ==================*/

/*
 * FMM_initializeFront_Eikonal2d() implements the callback 
 * function required by FMM_Core::FMM_initializeFront() to find and 
 * initialize the front.  
 */
void FMM_initializeFront_Eikonal2d(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int num_dims,
  int *grid_dims,
  double *dx);

/* 
 * FMM_updateGridPoint_Eikonal2d_Order1() computes and returns the 
 * updated phi value of the specified grid point using values of 
 * neighbors that have status "KNOWN". 
 */
double FMM_updateGridPoint_Eikonal2d_Order1(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int *grid_idx,
  int num_dims,
  int *grid_dims,
  double *dx);


/*==================== Function Definitions =========================*/


int solveEikonalEquation2d(
  double *phi,
  double *speed,
  double *mask,
  int spatial_derivative_order,
  int *grid_dims,
  double *dx)
{
  /* fast marching method data */
  FMM_CoreData *fmm_core_data;

  /* pointers to callback functions */
  updateGridPointFuncPtr updateGridPoint;
  initializeFrontFuncPtr initializeFront;

  /* auxiliary variables */
  int num_gridpoints;       /* number of grid points */
  int i,j;                  /* loop variables */
  double *ptr;              /* pointer to field data */


  /******************************************************
   * set up appropriate grid point update and front
   * detection/initialization functions based on the
   * specified spatial derivative order
   ******************************************************/
  initializeFront = &FMM_initializeFront_Eikonal2d;
  if (spatial_derivative_order == 1) {
    updateGridPoint = &FMM_updateGridPoint_Eikonal2d_Order1;
  } else if (spatial_derivative_order == 2) {
/*  KTC - ADD LATER
    updateGridPoint = &FMM_updateGridPoint_Eikonal2d_Order2;
*/
    fprintf(stderr,
           "ERROR: second-order spatial derivatives currently unsupported\n");
    return LSM_FMM_2D_ERR_INVALID_SPATIAL_DERIVATIVE_ORDER;
  } else {
    fprintf(stderr,
           "ERROR: Invalid spatial derivative order.  Only first-\n");
    fprintf(stderr,
           "       and second-order finite differences supported.\n");
    return LSM_FMM_2D_ERR_INVALID_SPATIAL_DERIVATIVE_ORDER;
  }

  /********************************************
   * set up FMM Field Data
   ********************************************/
  FMM_FieldData *fmm_field_data = 
    (FMM_FieldData*) malloc(sizeof(FMM_FieldData));
  if (!fmm_field_data) return LSM_FMM_2D_ERR_FMM_DATA_CREATION_ERROR;
  fmm_field_data->phi   = phi;
  fmm_field_data->speed = speed;
   
  /********************************************
   * initialize phi 
   ********************************************/
  num_gridpoints = 1;
  for (i = 0; i < LSM_FMM_2D_NDIM; i++) {
    num_gridpoints *= grid_dims[i];
  }
  for (i = 0, ptr = phi; i < num_gridpoints; i++, ptr++) {
    if ((*ptr) < 0) { /* set negative values to default value */
      *ptr = LSM_FMM_2D_DEFAULT_UPDATE_VALUE;
    }
  }

  /********************************************
   * initialize FMM Core Data
   ********************************************/
  fmm_core_data = FMM_Core_createFMM_CoreData(
    fmm_field_data,
    LSM_FMM_2D_NDIM,
    grid_dims,
    dx,
    initializeFront,
    updateGridPoint);
  if (!fmm_core_data) return LSM_FMM_2D_ERR_FMM_DATA_CREATION_ERROR;

  /* initialize grid points around the front */ 
  FMM_Core_initializeFront(fmm_core_data); 

  /* mark grid points outside of domain */
  for (j = 0; j < grid_dims[1]; j++) {
    for (i = 0; i < grid_dims[0]; i++) {
      int idx_ij = LSM_FMM_2D_IDX(i,j);

      if ((mask) && (mask[idx_ij] < 0)) {
        int grid_idx[LSM_FMM_2D_NDIM];
        grid_idx[0] = i; grid_idx[1] = j;
        FMM_Core_markPointOutsideDomain(fmm_core_data, grid_idx);

        /* set phi to DBL_MAX (i.e. infinity) */
        phi[idx_ij] = DBL_MAX;
      }

      if (speed[idx_ij] < LSM_FMM_2D_ZERO_TOL) {

        int grid_idx[LSM_FMM_2D_NDIM];
        grid_idx[0] = i; grid_idx[1] = j;
        FMM_Core_markPointOutsideDomain(fmm_core_data, grid_idx);

        /* speed is zero, so set phi to be DBL_MAX (i.e. infinity) */
        phi[idx_ij] = DBL_MAX;
      }

    }
  } 

  /* update remaining grid points */
  while (FMM_Core_moreGridPointsToUpdate(fmm_core_data)) {
    FMM_Core_advanceFront(fmm_core_data);
  }

  /* clean up memory */
  FMM_Core_destroyFMM_CoreData(fmm_core_data);
  free(fmm_field_data);

  return LSM_FMM_2D_ERR_SUCCESS;
}

void FMM_initializeFront_Eikonal2d(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int num_dims,
  int *grid_dims,
  double *dx)
{
  /* FMM Field Data variables */
  double *phi   = fmm_field_data->phi;
  double *speed = fmm_field_data->speed;
   
  /* grid variables */
  int grid_idx[LSM_FMM_2D_NDIM];

  /* auxilliary variables */
  int i,j;  /* loop variables for grid */
  int idx_ij;

  /* unused function parameters */
  (void) num_dims;


  /*
   * loop through cells in grid and initialize points on the boundary
   * for Eikonal equation.
   */
  for (j=0; j<grid_dims[1]; j++) {
    for (i=0; i<grid_dims[0]; i++) {

      /* compute index for (i,j) grid point */
      idx_ij = LSM_FMM_2D_IDX(i,j);

      /* set grid points on the initial front */
      if (phi[idx_ij] > -LSM_FMM_2D_ZERO_TOL) {

        /* the value for phi(i,j) has already been provided */
        grid_idx[0] = i; grid_idx[1] = j;
        FMM_Core_setInitialFrontPoint(fmm_core_data, grid_idx,
                                      phi[idx_ij]);
      } 

    }
  }  /* end loop over grid */

}


double FMM_updateGridPoint_Eikonal2d_Order1(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int *grid_idx,
  int num_dims,
  int *grid_dims,
  double *dx)
{
  int *gridpoint_status = FMM_Core_getGridPointStatusDataArray(fmm_core_data);

  /* FMM Field Data variables */
  double *phi   = fmm_field_data->phi; 
  double *speed = fmm_field_data->speed; 

  /* variables used in phi update */
  PointStatus  neighbor_status;
  double phi_upwind[LSM_FMM_2D_NDIM];
  double phi_plus;
  double inv_dx_sq; 
  int offset[LSM_FMM_2D_NDIM]; 
  int neighbor[LSM_FMM_2D_NDIM];

  /* coefficients of quadratic equation for phi */
  double phi_A = 0;
  double phi_B = 0;
  double phi_C = 0;
  double discriminant;
  double phi_updated;

  /* auxilliary variables */
  int dir;  /* loop variable for spatial directions */
  int l;    /* extra loop variable */ 
  int idx_cur_gridpoint, idx_neighbor;

  /* unused function parameters */
  (void) num_dims;

  /* compute index for current grid point */
  idx_cur_gridpoint = LSM_FMM_2D_IDX(grid_idx[0],grid_idx[1]);

  /* calculate update to phi */
  for (dir = 0; dir < LSM_FMM_2D_NDIM; dir++) { /* loop over coord directions */

    /* reset offset */
    for (l = 0; l < LSM_FMM_2D_NDIM; l++) { 
      offset[l] = 0; 
    }

    /* find "upwind" direction and phi value */
    phi_upwind[dir] = DBL_MAX;

    /* check minus direction */
    offset[dir] = -1;
    neighbor[0] = grid_idx[0] + offset[0];
    neighbor[1] = grid_idx[1] + offset[1];
    if (!LSM_FMM_2D_IDX_OUT_OF_BOUNDS(neighbor[0],neighbor[1],grid_dims)) {
      idx_neighbor = LSM_FMM_2D_IDX(neighbor[0],neighbor[1]);
      neighbor_status = (PointStatus) gridpoint_status[idx_neighbor];
      if (KNOWN == neighbor_status) {
        phi_upwind[dir] = phi[idx_neighbor];
      }
    }

    /* check plus direction */
    offset[dir] = 1;
    neighbor[0] = grid_idx[0] + offset[0];
    neighbor[1] = grid_idx[1] + offset[1];
    if (!LSM_FMM_2D_IDX_OUT_OF_BOUNDS(neighbor[0],neighbor[1],grid_dims)) {
      idx_neighbor = LSM_FMM_2D_IDX(neighbor[0],neighbor[1]);
      neighbor_status = (PointStatus) gridpoint_status[idx_neighbor];
      if (KNOWN == neighbor_status) {
        phi_plus = phi[idx_neighbor];

        /* 
         * choosing the upwind direction to be the direction
         * with the smaller abs(phi) value gives a consistent 
         * solution to the "upwind" Eikonal equation.
         * NOTE: to avoid having to use the C math library,
         *       we define our own ABS macro
         */
        if (LSM_FMM_2D_ABS(phi_plus) < LSM_FMM_2D_ABS(phi_upwind[dir])) {
          phi_upwind[dir] = phi_plus;
        }
      }
    }

    /*
     * accumulate coefficients for phi if either of the neighbors are "known"
     */
    if (phi_upwind[dir] < DBL_MAX) {
      /* accumulate coefs for phi */ 
      inv_dx_sq = 1/dx[dir]/dx[dir];
      phi_A += inv_dx_sq;
      phi_B += phi_upwind[dir]*inv_dx_sq;
      phi_C += phi_upwind[dir]*phi_upwind[dir]*inv_dx_sq;
    }

  } /* loop over coordinate directions */

  /* complete computation of phi_B and phi_C */
  phi_B *= -2.0;
  phi_C -= 1/speed[idx_cur_gridpoint]/speed[idx_cur_gridpoint];

  /* compute phi by solving quadratic equation */
  discriminant = phi_B*phi_B - 4*phi_A*phi_C;
  phi_updated = DBL_MAX;
  if (discriminant >= 0) {
    phi_updated = (-phi_B + sqrt(discriminant))/2/phi_A;
  } else {
    fprintf(stderr,"ERROR: phi update - discriminant negative!!!\n");
  }

  /* set phi at current grid point */
  phi[idx_cur_gridpoint] = phi_updated;

  return phi_updated;
}
