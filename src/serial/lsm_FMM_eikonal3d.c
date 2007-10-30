/*
 * File:        lsm_FMM_eikonal_3d.c
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.5 $
 * Modified:    $Date: 2007/05/24 12:01:11 $
 * Description: Implementation of 3D Fast Marching Method for Eikonal equation
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "FMM_Core.h"
#include "FMM_Heap.h"
#include "lsm_fast_marching_method.h"


/*===================== lsm_FMM_3d Data Structures ====================*/
struct FMM_FieldData {
  double *phi;                 /* solution to Eikonal equation */
  double *speed;               /* speed function               */
};


/*======================= lsm_FMM_3d Constants ==========================*/
#define LSM_FMM_3D_NDIM                   (3)
#define LSM_FMM_3D_TRUE                   (1)
#define LSM_FMM_3D_FALSE                  (0)
#define LSM_FMM_3D_DEFAULT_UPDATE_VALUE   (-1)
#define LSM_FMM_3D_ZERO_TOL               (1e-8)


/*========================== Error Codes ============================*/
#define LSM_FMM_3D_ERR_SUCCESS                             (0)
#define LSM_FMM_3D_ERR_FMM_DATA_CREATION_ERROR             (1)
#define LSM_FMM_3D_ERR_INVALID_SPATIAL_DERIVATIVE_ORDER    (2)


/*========================= lsm_FMM_3d Macros ===========================*/
#define LSM_FMM_3D_IDX(i,j,k)     ( (i) + grid_dims[0]*(j)                 \
                                  + grid_dims[0]*grid_dims[1]*(k) )
#define LSM_FMM_3D_ABS(x)         ((x) > 0 ? (x) : -1.0*(x))          
#define LSM_FMM_3D_IDX_OUT_OF_BOUNDS(i,j,k,grid_dims)                      \
  ( ((i)<0) || ((i)>(grid_dims)[0]-1) || ((j)<0) || ((j)>(grid_dims)[1]-1) \
    || ( ((grid_dims)[2]>0) && (((k)<0) || ((k)>(grid_dims)[2]-1)) ) )


/*================== Helper Functions Declarations ==================*/

/*
 * FMM_initializeFront_Eikonal3d() implements the callback 
 * function required by FMM_Core::FMM_initializeFront() to find and 
 * initialize the front.  
 */
void FMM_initializeFront_Eikonal3d(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int num_dims,
  int *grid_dims,
  double *dx);

/* 
 * FMM_updateGridPoint_Eikonal3d_Order1() computes and returns the 
 * updated phi value of the specified grid point using values of 
 * neighbors that have status "KNOWN" and a first-order accurate 
 * discretization of the gradient operator.
 */
double FMM_updateGridPoint_Eikonal3d_Order1(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int *grid_idx,
  int num_dims,
  int *grid_dims,
  double *dx);

/*
 * FMM_updateGridPoint_Eikonal3d_Order2() computes and returns the
 * updated phi value of the specified grid point using values of
 * neighbors that have status "KNOWN" and a second-order accurate
 * discretization of the gradient operator when a sufficient number
 * of "KNOWN" neighboring grid points are available.  When there are
 * an insufficient number of "KNOWN" neighbors, the discretization
 * of the gradient drops to first-order accuracy.
 */
double FMM_updateGridPoint_Eikonal3d_Order2(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int *grid_idx,
  int num_dims,
  int *grid_dims,
  double *dx);


/*==================== Function Definitions =========================*/


int solveEikonalEquation3d(
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
  int i,j,k;                /* loop variables */
  double *ptr;              /* pointer to field data */


  /******************************************************
   * set up appropriate grid point update and front
   * detection/initialization functions based on the
   * specified spatial derivative order
   ******************************************************/
  initializeFront = &FMM_initializeFront_Eikonal3d;
  if (spatial_derivative_order == 1) {
    updateGridPoint = &FMM_updateGridPoint_Eikonal3d_Order1;
  } else if (spatial_derivative_order == 2) {
    updateGridPoint = &FMM_updateGridPoint_Eikonal3d_Order2;
  } else {
    fprintf(stderr,
           "ERROR: Invalid spatial derivative order.  Only first-\n");
    fprintf(stderr,
           "       and second-order finite differences supported.\n");
    return LSM_FMM_3D_ERR_INVALID_SPATIAL_DERIVATIVE_ORDER;
  }

  /********************************************
   * set up FMM Field Data
   ********************************************/
  FMM_FieldData *fmm_field_data = 
    (FMM_FieldData*) malloc(sizeof(FMM_FieldData));
  if (!fmm_field_data) return LSM_FMM_3D_ERR_FMM_DATA_CREATION_ERROR;
  fmm_field_data->phi   = phi;
  fmm_field_data->speed = speed;
   
  /********************************************
   * initialize phi 
   ********************************************/
  num_gridpoints = 1;
  for (i = 0; i < LSM_FMM_3D_NDIM; i++) {
    num_gridpoints *= grid_dims[i];
  }
  for (i = 0, ptr = phi; i < num_gridpoints; i++, ptr++) {
    if ((*ptr) < 0) { /* set negative values to default value */
      *ptr = LSM_FMM_3D_DEFAULT_UPDATE_VALUE;
    }
  }

  /********************************************
   * initialize FMM Core Data
   ********************************************/
  fmm_core_data = FMM_Core_createFMM_CoreData(
    fmm_field_data,
    LSM_FMM_3D_NDIM,
    grid_dims,
    dx,
    initializeFront,
    updateGridPoint);
  if (!fmm_core_data) return LSM_FMM_3D_ERR_FMM_DATA_CREATION_ERROR;

  /* mark grid points outside of domain */
  for (k = 0; k < grid_dims[2]; k++) {
    for (j = 0; j < grid_dims[1]; j++) {
      for (i = 0; i < grid_dims[0]; i++) {
        int idx_ijk = LSM_FMM_3D_IDX(i,j,k);

        if ((mask) && (mask[idx_ijk] < 0)) {
          int grid_idx[LSM_FMM_3D_NDIM];
          grid_idx[0] = i; grid_idx[1] = j; grid_idx[2] = k;
          FMM_Core_markPointOutsideDomain(fmm_core_data, grid_idx);

          /* set phi to DBL_MAX (i.e. infinity) */
          phi[idx_ijk] = DBL_MAX;
        }

        if (speed[idx_ijk] < LSM_FMM_3D_ZERO_TOL) {

          int grid_idx[LSM_FMM_3D_NDIM];
          grid_idx[0] = i; grid_idx[1] = j; grid_idx[2] = k;
          FMM_Core_markPointOutsideDomain(fmm_core_data, grid_idx);

          /* speed is zero, so set phi to be DBL_MAX (i.e. infinity) */
          phi[idx_ijk] = DBL_MAX;
        } 

      }
    }
  } /* end loop over grid to mark points outside of domain */

  /* initialize grid points around the front */ 
  FMM_Core_initializeFront(fmm_core_data); 

  /* update remaining grid points */
  while (FMM_Core_moreGridPointsToUpdate(fmm_core_data)) {
    FMM_Core_advanceFront(fmm_core_data);
  }

  /* clean up memory */
  FMM_Core_destroyFMM_CoreData(fmm_core_data);
  free(fmm_field_data);

  return LSM_FMM_3D_ERR_SUCCESS;
}

void FMM_initializeFront_Eikonal3d(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int num_dims,
  int *grid_dims,
  double *dx)
{
  /* FMM Field Data variables */
  double *phi   = fmm_field_data->phi;
   
  /* grid variables */
  int grid_idx[LSM_FMM_3D_NDIM];

  /* auxilliary variables */
  int i,j,k;  /* loop variables for grid */
  int idx_ijk;

  /* unused function parameters */
  (void) num_dims;


  /*
   * loop through cells in grid and initialize points on the boundary
   * for Eikonal equation.
   */
  for (k=0; k<grid_dims[2]; k++) {
    for (j=0; j<grid_dims[1]; j++) {
      for (i=0; i<grid_dims[0]; i++) {

        /* compute index for (i,j,k) grid point */
        idx_ijk = LSM_FMM_3D_IDX(i,j,k);

        /* set grid points on the initial front */
        if (phi[idx_ijk] > -LSM_FMM_3D_ZERO_TOL) {

          /* the value for phi(i,j,k) has already been provided */
          grid_idx[0] = i; grid_idx[1] = j; grid_idx[2] = k;
          FMM_Core_setInitialFrontPoint(fmm_core_data, grid_idx,
                                        phi[idx_ijk]);
        }

      }
    }
  }  /* end loop over grid */
  
}


double FMM_updateGridPoint_Eikonal3d_Order1(
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
  PointStatus neighbor_status;
  double phi_upwind;
  double phi_plus;
  double inv_dx_sq; 
  int offset[LSM_FMM_3D_NDIM]; 
  int neighbor[LSM_FMM_3D_NDIM];

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
  idx_cur_gridpoint = LSM_FMM_3D_IDX(grid_idx[0],grid_idx[1],grid_idx[2]);

  /* calculate update to phi */
  for (dir = 0; dir < LSM_FMM_3D_NDIM; dir++) { /* loop over coord directions */

    /* reset offset */
    for (l = 0; l < LSM_FMM_3D_NDIM; l++) { 
      offset[l] = 0; 
    }

    /* find "upwind" direction and phi value */
    phi_upwind = DBL_MAX;

    /* check minus direction */
    offset[dir] = -1;
    neighbor[0] = grid_idx[0] + offset[0];
    neighbor[1] = grid_idx[1] + offset[1];
    neighbor[2] = grid_idx[2] + offset[2];
    if (!LSM_FMM_3D_IDX_OUT_OF_BOUNDS(neighbor[0],neighbor[1],neighbor[2],
                                      grid_dims)) {
      idx_neighbor = LSM_FMM_3D_IDX(neighbor[0],neighbor[1],neighbor[2]);
      neighbor_status = (PointStatus) gridpoint_status[idx_neighbor];
      if (KNOWN == neighbor_status) {
        phi_upwind = phi[idx_neighbor];
      }
    }

    /* check plus direction */
    offset[dir] = 1;
    neighbor[0] = grid_idx[0] + offset[0];
    neighbor[1] = grid_idx[1] + offset[1];
    neighbor[2] = grid_idx[2] + offset[2];
    if (!LSM_FMM_3D_IDX_OUT_OF_BOUNDS(neighbor[0],neighbor[1],neighbor[2],
                                      grid_dims)) {
      idx_neighbor = LSM_FMM_3D_IDX(neighbor[0],neighbor[1],neighbor[2]);
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
        if (LSM_FMM_3D_ABS(phi_plus) < LSM_FMM_3D_ABS(phi_upwind)) {
          phi_upwind = phi_plus;
        }
      }
    }

    /*
     * accumulate coefficients for phi if either of the neighbors are "KNOWN"
     */
    if (phi_upwind < DBL_MAX) {
      /* accumulate coefs for phi */ 
      inv_dx_sq = 1/dx[dir]/dx[dir];
      phi_A += inv_dx_sq;
      phi_B += inv_dx_sq*phi_upwind;
      phi_C += inv_dx_sq*phi_upwind*phi_upwind;
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


double FMM_updateGridPoint_Eikonal3d_Order2(
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
  PointStatus neighbor_status;
  double phi_upwind1, phi_upwind2;
  double phi_plus;
  int second_order_switch;
  double inv_dx_sq; 
  int offset[LSM_FMM_3D_NDIM]; 
  int neighbor1[LSM_FMM_3D_NDIM];
  int neighbor2[LSM_FMM_3D_NDIM];

  /* coefficients of quadratic equation for phi */
  double phi_A = 0;
  double phi_B = 0;
  double phi_C = 0;
  double discriminant;
  double phi_updated;

  /* auxilliary variables */
  int dir;  /* loop variable for spatial directions */
  int l;    /* extra loop variable */ 
  int idx_cur_gridpoint, idx_neighbor1, idx_neighbor2;

  /* unused function parameters */
  (void) num_dims;


  /* compute index for current grid point */
  idx_cur_gridpoint = LSM_FMM_3D_IDX(grid_idx[0],grid_idx[1],grid_idx[2]);

  /* calculate update to phi */
  for (dir = 0; dir < LSM_FMM_3D_NDIM; dir++) { /* loop over coord directions */

    /* reset offset */
    for (l = 0; l < LSM_FMM_3D_NDIM; l++) { 
      offset[l] = 0; 
    }

    /* reset phi_upwind1 and phi_upwind2 to DBL_MAX */
    phi_upwind1 = DBL_MAX;
    phi_upwind2 = DBL_MAX;
    
    /* reset second_order_switch to 0 (i.e. assume there are not enough */
    /* KNOWN neighbors for second-order discretization.                 */
    second_order_switch = 0;

    /* check minus direction */
    offset[dir] = -1;
    neighbor1[0] = grid_idx[0] + offset[0];
    neighbor1[1] = grid_idx[1] + offset[1];
    neighbor1[2] = grid_idx[2] + offset[2];
    neighbor2[0] = grid_idx[0] + 2*offset[0];
    neighbor2[1] = grid_idx[1] + 2*offset[1];
    neighbor2[2] = grid_idx[2] + 2*offset[2];
    if (!LSM_FMM_3D_IDX_OUT_OF_BOUNDS(neighbor1[0],neighbor1[1],neighbor1[2],
                                      grid_dims)) {
      idx_neighbor1 = LSM_FMM_3D_IDX(neighbor1[0],neighbor1[1],neighbor1[2]);
      neighbor_status = (PointStatus) gridpoint_status[idx_neighbor1];
      if (KNOWN == neighbor_status) {
        phi_upwind1 = phi[idx_neighbor1];

        /* check for neighbor required for second-order accuracy */
        if (!LSM_FMM_3D_IDX_OUT_OF_BOUNDS(neighbor2[0], neighbor2[1],
                                          neighbor2[2], grid_dims)) {
          idx_neighbor2 = LSM_FMM_3D_IDX(neighbor2[0],neighbor2[1],
                                         neighbor2[2]);
          neighbor_status = (PointStatus) gridpoint_status[idx_neighbor2];
          if ( (KNOWN == neighbor_status) &&
               (LSM_FMM_3D_ABS(phi_upwind2) < LSM_FMM_3D_ABS(phi_upwind1)) ) {
            phi_upwind2 = phi[idx_neighbor2];
            second_order_switch = 1;
          }
        }
      
      } /* end case: first-order neighbor is KNOWN */
    }

    /* check plus direction */
    offset[dir] = 1;
    neighbor1[0] = grid_idx[0] + offset[0];
    neighbor1[1] = grid_idx[1] + offset[1];
    neighbor1[2] = grid_idx[2] + offset[2];
    neighbor2[0] = grid_idx[0] + 2*offset[0];
    neighbor2[1] = grid_idx[1] + 2*offset[1];
    neighbor2[2] = grid_idx[2] + 2*offset[2];
    if (!LSM_FMM_3D_IDX_OUT_OF_BOUNDS(neighbor1[0],neighbor1[1],neighbor1[2],
                                      grid_dims)) {
      idx_neighbor1 = LSM_FMM_3D_IDX(neighbor1[0],neighbor1[1],neighbor1[2]);
      neighbor_status = (PointStatus) gridpoint_status[idx_neighbor1];
      if (KNOWN == neighbor_status) {
        phi_plus = phi[idx_neighbor1];

        /* 
         * choosing the upwind direction to be the direction
         * with the smaller abs(phi) value gives a consistent 
         * solution to the "upwind" Eikonal equation.
         * NOTE: to avoid having to use the C math library,
         *       we define our own ABS macro
         */
        if (LSM_FMM_3D_ABS(phi_plus) < LSM_FMM_3D_ABS(phi_upwind1)) {
          phi_upwind1 = phi_plus;
          phi_upwind2 = DBL_MAX;
          second_order_switch = 0;

          /* check for neighbor required for second-order accuracy */
          if (!LSM_FMM_3D_IDX_OUT_OF_BOUNDS(neighbor2[0], neighbor2[1],
                                            neighbor2[2], grid_dims)) {
            idx_neighbor2 = LSM_FMM_3D_IDX(neighbor2[0],neighbor2[1],
                                           neighbor2[2]);
            neighbor_status = (PointStatus) gridpoint_status[idx_neighbor2];
            if ( (KNOWN == neighbor_status) &&
                 (LSM_FMM_3D_ABS(phi_upwind2) < LSM_FMM_3D_ABS(phi_upwind1)) ) {
              phi_upwind2 = phi[idx_neighbor2];
              second_order_switch = 1;
            }
          }

        } /* end case: plus is upwind direction */
      
      } /* end case: first-order neighbor is KNOWN */
    }

    /*
     * accumulate coefficients for phi if either of the neighbors are "KNOWN"
     */
    if (phi_upwind1 < DBL_MAX) {
      /* temporary variables */
      double one_plus_switch_over_two = 1.0+second_order_switch/2.0;
      double phi_upwind_contrib;

      /* set phi_upwind_contrib to be first- or second-order */
      /* contribution based on value of second_order_switch  */
      if (second_order_switch == 1) {
        phi_upwind_contrib = phi_upwind1+second_order_switch*phi_upwind1
                           - second_order_switch/2.0*phi_upwind2;
      } else {
        phi_upwind_contrib = phi_upwind1;
      }

      /* accumulate coefs for phi */
      inv_dx_sq = 1/dx[dir]/dx[dir];
      phi_A += inv_dx_sq*one_plus_switch_over_two*one_plus_switch_over_two;
      phi_B += inv_dx_sq*one_plus_switch_over_two
                        *phi_upwind_contrib;
      phi_C += inv_dx_sq*phi_upwind_contrib*phi_upwind_contrib;
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
