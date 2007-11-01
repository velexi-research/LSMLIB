/*
 * File:        lsm_FMM_eikonal.c
 * Copyright:   (c) 2005-2008 Kevin T. Chu
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date: 2006/08/13 13:30:52 $
 * Description: Implementation of Fast Marching Method for Eikonal equation
 */

/*! \file lsm_FMM_eikonal.c
 *
 * \brief
 * @ref lsm_FMM_eikonal.c provides "generic" implementations of 
 *      first- and second-order accurate grid points updates within
 *      the context of the Fast Marching Method algorithm.  The 
 *      code is "templated" on the number of dimensions through
 *      the use of macro definitions that MUST be provided by
 *      the user.  The user is also responsible for providing
 *      several macros for grid index calculations that depend
 *      on the dimension of the problem.
 *
 *
 * <h3> Usage: </h3>
 *
 * -# Provide the following macros:
 *    -# FMM_NDIM:  the number of spatial dimensions.
 *    -# FMM_EIKONAL_SOLVE_EIKONAL_EQUATION:  desired name of function 
 *       that solves the Eikonal equation.
 *    -# FMM_EIKONAL_INITIALIZE_FRONT:  desired name of function that
 *       initializes the values on the front.
 *    -# FMM_EIKONAL_UPDATE_GRID_POINT_ORDER1:  desired name of function 
 *       that updates the value of the solution at grid points using
 *       a first-order accurate discretization
 *    -# FMM_EIKONAL_UPDATE_GRID_POINT_ORDER2:  desired name of function 
 *       that updates the value of the solution at grid points using
 *       a second-order accurate discretization
 * -# Include this file at the end of the implementation file
 *    for the n-dimentsional Eikonal equation solver.
 * -# Compile code.
 *
 *
 * <h3> NOTES: </h3>
 * - Because this code depends on macros, care must be taken to 
 *   ensure that macros do not conflict.
 *
 */

#ifndef included_lsm_FMM_eikonal_c
#define included_lsm_FMM_eikonal_c
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "FMM_Core.h"
#include "FMM_Heap.h"
#include "FMM_Macros.h"


/*
 * This macro protect against misuse of the code in this file.  It will
 * cause the compiler to complain.
 */
#ifndef FMM_NDIM
#error "lsm_FMM_eikonal: required macro FMM_NDIM not defined!"
#endif FMM_NDIM
#ifndef FMM_EIKONAL_SOLVE_EIKONAL_EQUATION
#error "lsm_FMM_eikonal: required macro FMM_EIKONAL_SOLVE_EIKONAL_EQUATION not defined!"
#endif FMM_EIKONAL_SOLVE_EIKONAL_EQUATION 
#ifndef FMM_EIKONAL_INITIALIZE_FRONT
#error "lsm_FMM_eikonal: required macro FMM_EIKONAL_INITIALIZE_FRONT not defined!"
#endif FMM_EIKONAL_INITIALIZE_FRONT 
#ifndef FMM_EIKONAL_UPDATE_GRID_POINT_ORDER1
#error "lsm_FMM_eikonal: required macro FMM_EIKONAL_UPDATE_GRID_POINT_ORDER1 not defined!"
#endif FMM_EIKONAL_UPDATE_GRID_POINT_ORDER1 
#ifndef FMM_EIKONAL_UPDATE_GRID_POINT_ORDER2
#error "lsm_FMM_eikonal: required macro FMM_EIKONAL_UPDATE_GRID_POINT_ORDER2 not defined!"
#endif FMM_EIKONAL_UPDATE_GRID_POINT_ORDER2 


/*===================== lsm_FMM Data Structures ====================*/
struct FMM_FieldData {
  double *phi;                 /* solution to Eikonal equation */
  double *speed;               /* speed function               */
};


/*============= FMM Eikonal Equation Solver Functions ===============*/

/*
 * FMM_EIKONAL_INITIALIZE_FRONT() implements the callback
 * function required by FMM_Core::FMM_initializeFront() to find 
 * and initialize the front.
 */
void FMM_EIKONAL_INITIALIZE_FRONT(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int num_dims,
  int *grid_dims,
  double *dx);

/* 
 * FMM_EIKONAL_UPDATE_GRID_POINT_ORDER1() implements the callback 
 * function required by FMM_Core::FMM_Core_updateNeighbors() to
 * update the soluttion at a grid point.  It computes and returns 
 * the updated phi value of the specified grid point using values of 
 * neighbors that have status "KNOWN" and a first-order accurate 
 * discretization of the gradient operator.
 */
double FMM_EIKONAL_UPDATE_GRID_POINT_ORDER1(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int *grid_idx,
  int num_dims,
  int *grid_dims,
  double *dx);

/* 
 * FMM_EIKONAL_UPDATE_GRID_POINT_ORDER2() implements the callback  
 * function required by FMM_Core::FMM_Core_updateNeighbors() to 
 * update the soluttion at a grid point.  It computes and returns 
 * the updated phi value of the specified grid point using values of 
 * neighbors that have status "KNOWN" and a second-order accurate 
 * discretization of the gradient operator when a sufficient number 
 * of "KNOWN" neighboring grid points are available.  When there are
 * an insufficient number of "KNOWN" neighbors, the discretization
 * of the gradient drops to first-order accuracy.
 */
double FMM_EIKONAL_UPDATE_GRID_POINT_ORDER2(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int *grid_idx,
  int num_dims,
  int *grid_dims,
  double *dx);


/*==================== Function Definitions =========================*/


int FMM_EIKONAL_SOLVE_EIKONAL_EQUATION(
  double *phi,
  double *speed,
  double *mask,
  int spatial_discretization_order,
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
  int i, idx;               /* loop variables */


  /******************************************************
   * set up appropriate grid point update and front
   * detection/initialization functions based on the
   * specified spatial derivative order
   ******************************************************/
  initializeFront = &FMM_EIKONAL_INITIALIZE_FRONT;
  if (spatial_discretization_order == 1) {
    updateGridPoint = &FMM_EIKONAL_UPDATE_GRID_POINT_ORDER1;
  } else if (spatial_discretization_order == 2) {
    updateGridPoint = &FMM_EIKONAL_UPDATE_GRID_POINT_ORDER2;
  } else {
    fprintf(stderr,
           "ERROR: Invalid spatial derivative order.  Only first-\n");
    fprintf(stderr,
           "       and second-order finite differences supported.\n");
    return LSM_FMM_ERR_INVALID_SPATIAL_DISCRETIZATION_ORDER;
  }

  /********************************************
   * set up FMM Field Data
   ********************************************/
  FMM_FieldData *fmm_field_data = 
    (FMM_FieldData*) malloc(sizeof(FMM_FieldData));
  if (!fmm_field_data) return LSM_FMM_ERR_FMM_DATA_CREATION_ERROR;
  fmm_field_data->phi   = phi;
  fmm_field_data->speed = speed;
   
  /********************************************
   * initialize FMM Core Data
   ********************************************/
  fmm_core_data = FMM_Core_createFMM_CoreData(
    fmm_field_data,
    FMM_NDIM,
    grid_dims,
    dx,
    initializeFront,
    updateGridPoint);
  if (!fmm_core_data) return LSM_FMM_ERR_FMM_DATA_CREATION_ERROR;

  /********************************************
   * initialize phi and mark grid points
   * outside of the mathematical/physical 
   * domain
   ********************************************/
  num_gridpoints = 1;
  for (i = 0; i < FMM_NDIM; i++) {
    num_gridpoints *= grid_dims[i];
  }

  for (idx = 0; idx < num_gridpoints; idx++) {

    /* temporary variables */
    int grid_idx[FMM_NDIM];   /* grid index */
    int idx_remainder = idx;

    /* compute grid_idx */
    for (i = 0; i < FMM_NDIM; i++) {
      grid_idx[i] = idx_remainder%grid_dims[i];
      idx_remainder -= grid_idx[i];
      idx_remainder /= grid_dims[i];
    }

    /* grid points with a negative mask value are taken to */
    /* be outside of the mathemtatical/physical domain     */
    if ((mask) && (mask[idx] < 0)) {

      FMM_Core_markPointOutsideDomain(fmm_core_data, grid_idx);

      /* set phi to DBL_MAX (i.e. infinity) */
      phi[idx] = DBL_MAX;
    }

    /* grid points with a non-positive speed are taken to */
    /* be outside of the mathemtatical/physical domain    */
    if (speed[idx] < LSM_FMM_ZERO_TOL) {

      FMM_Core_markPointOutsideDomain(fmm_core_data, grid_idx);

      /* speed is zero, so set phi to be DBL_MAX (i.e. infinity) */
      phi[idx] = DBL_MAX;
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

  return LSM_FMM_ERR_SUCCESS;
}

void FMM_EIKONAL_INITIALIZE_FRONT(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int num_dims,
  int *grid_dims,
  double *dx)
{
  /* FMM Field Data variables */
  double *phi   = fmm_field_data->phi;

  /* grid variables */
  int grid_idx[FMM_NDIM];

  /* auxilliary variables */
  int num_gridpoints;
  int i,idx;         /* loop variables */

  /* unused function parameters */
  (void) num_dims;


  /*
   * loop through cells in grid and initialize points on the boundary
   * for Eikonal equation.
   */
  num_gridpoints = 1;
  for (i = 0; i < FMM_NDIM; i++) {
    num_gridpoints *= grid_dims[i];
  }

  for (idx = 0; idx < num_gridpoints; idx++) {

    /* temporary variables */
    int idx_remainder = idx;

    /* compute grid_idx */
    for (i = 0; i < FMM_NDIM; i++) {
      grid_idx[i] = idx_remainder%grid_dims[i];
      idx_remainder -= grid_idx[i];
      idx_remainder /= grid_dims[i];
    }

    /* set grid points on the initial front */
    if (phi[idx] > -LSM_FMM_ZERO_TOL) {

      /* the value for phi(i,j) has already been provided */
      FMM_Core_setInitialFrontPoint(fmm_core_data, grid_idx,
                                    phi[idx]);
    }

  }  /* end loop over grid */

}


double FMM_EIKONAL_UPDATE_GRID_POINT_ORDER1(
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
  int offset[FMM_NDIM]; 
  int neighbor[FMM_NDIM];

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
  int grid_idx_out_of_bounds;

  /* unused function parameters */
  (void) num_dims;

  /* compute index for current grid point */
  LSM_FMM_IDX(idx_cur_gridpoint, grid_idx, grid_dims);

  /* calculate update to phi */
  for (dir = 0; dir < FMM_NDIM; dir++) { 

    /* reset offset */
    for (l = 0; l < FMM_NDIM; l++) { 
      offset[l] = 0; 
    }

    /* find "upwind" direction and phi value */
    phi_upwind = DBL_MAX;

    /* check minus direction */
    offset[dir] = -1;
    for (l = 0; l < FMM_NDIM; l++) { 
      neighbor[l] = grid_idx[l] + offset[l];
    }
    LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,neighbor,grid_dims);
    if (!grid_idx_out_of_bounds) {
      LSM_FMM_IDX(idx_neighbor, neighbor, grid_dims);
      neighbor_status = (PointStatus) gridpoint_status[idx_neighbor];
      if (KNOWN == neighbor_status) {
        phi_upwind = phi[idx_neighbor];
      }
    }

    /* check plus direction */
    offset[dir] = 1;
    for (l = 0; l < FMM_NDIM; l++) { 
      neighbor[l] = grid_idx[l] + offset[l];
    }
    LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,neighbor,grid_dims);
    if (!grid_idx_out_of_bounds) {
      LSM_FMM_IDX(idx_neighbor, neighbor, grid_dims);
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
        if (LSM_FMM_ABS(phi_plus) < LSM_FMM_ABS(phi_upwind)) {
          phi_upwind = phi_plus;
        }
      }
    }

    /*
     * accumulate coefficients for phi if any of the neighbors are "KNOWN"
     */
    if (phi_upwind < DBL_MAX) {
      /* accumulate coefs for phi */ 
      inv_dx_sq = 1/dx[dir]; inv_dx_sq *= inv_dx_sq; 
      phi_A += inv_dx_sq;
      phi_B += inv_dx_sq*phi_upwind;
      phi_C += inv_dx_sq*phi_upwind*phi_upwind;
    }

  } /* loop over coordinate directions */

  /* complete computation of phi_B and phi_C */
  phi_B *= -2.0;
  phi_C -= 1/speed[idx_cur_gridpoint]/speed[idx_cur_gridpoint];

  /* check that phi_A is nonzero */
  if (LSM_FMM_ABS(phi_A) == 0) {
    fprintf(stderr,"ERROR: phi update - no KNOWN neighbors!!!\n");
    fprintf(stderr,"       phi set to 'infinity'.\n");
    return DBL_MAX;
  }

  /* compute phi by solving quadratic equation */
  discriminant = phi_B*phi_B - 4.0*phi_A*phi_C;
  phi_updated = DBL_MAX;
  if (discriminant >= 0) {
    phi_updated = 0.5*(-phi_B + sqrt(discriminant))/phi_A;
  } else {
    fprintf(stderr,"ERROR: phi update - discriminant negative!!!\n");
    fprintf(stderr,"       phi set to 'infinity'.\n");
    fprintf(stderr,"       discriminant = %g,", discriminant);
    fprintf(stderr," grid_idx = (");
    for (l = 0; l < FMM_NDIM-1; l++) { 
      fprintf(stderr,"%d,", grid_idx[l]);
    }
    fprintf(stderr,"%d)\n",grid_idx[l]);
  }

  /* set phi at current grid point */
  phi[idx_cur_gridpoint] = phi_updated;

  return phi_updated;
}


double FMM_EIKONAL_UPDATE_GRID_POINT_ORDER2(
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
  int offset[FMM_NDIM]; 
  int neighbor1[FMM_NDIM];
  int neighbor2[FMM_NDIM];

  /* coefficients of quadratic equation for phi */
  double phi_A = 0;
  double phi_B = 0;
  double phi_C = 0;
  double discriminant;
  double phi_updated;
  double max_dx;

  /* auxilliary variables */
  int dir;  /* loop variable for spatial directions */
  int l;    /* extra loop variable */ 
  int idx_cur_gridpoint, idx_neighbor1, idx_neighbor2;
  int grid_idx_out_of_bounds;

  /* unused function parameters */
  (void) num_dims;

  /* compute index for current grid point */
  LSM_FMM_IDX(idx_cur_gridpoint, grid_idx, grid_dims);

  /* calculate update to phi */
  for (dir = 0; dir < FMM_NDIM; dir++) { 

    /* reset offset */
    for (l = 0; l < FMM_NDIM; l++) { 
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
    for (l = 0; l < FMM_NDIM; l++) { 
      neighbor1[l] = grid_idx[l] + offset[l];
      neighbor2[l] = grid_idx[l] + 2*offset[l];
    }
    LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,neighbor1,grid_dims);
    if (!grid_idx_out_of_bounds) {
      LSM_FMM_IDX(idx_neighbor1, neighbor1, grid_dims);
      neighbor_status = (PointStatus) gridpoint_status[idx_neighbor1];
      if (KNOWN == neighbor_status) {
        phi_upwind1 = phi[idx_neighbor1];

        /* check for neighbor required for second-order accuracy */
        LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,neighbor2,grid_dims);
        if (!grid_idx_out_of_bounds) {
          LSM_FMM_IDX(idx_neighbor2, neighbor2, grid_dims);
          neighbor_status = (PointStatus) gridpoint_status[idx_neighbor2];
          if ( (KNOWN == neighbor_status) &&
               (  LSM_FMM_ABS(phi[idx_neighbor2]) 
               <= LSM_FMM_ABS(phi_upwind1)) ) {
            phi_upwind2 = phi[idx_neighbor2];
            second_order_switch = 1;
          }
        } 

      } /* end case: first-order neighbor is KNOWN */
    } 

    /* check plus direction */
    offset[dir] = 1;
    for (l = 0; l < FMM_NDIM; l++) { 
      neighbor1[l] = grid_idx[l] + offset[l];
      neighbor2[l] = grid_idx[l] + 2*offset[l];
    }
    LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,neighbor1,grid_dims);
    if (!grid_idx_out_of_bounds) {
      LSM_FMM_IDX(idx_neighbor1, neighbor1, grid_dims);
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
        if (LSM_FMM_ABS(phi_plus) < LSM_FMM_ABS(phi_upwind1)) {
          phi_upwind1 = phi_plus;
          phi_upwind2 = DBL_MAX;
          second_order_switch = 0;

          /* check for neighbor required for second-order accuracy */
          LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,neighbor2,grid_dims);
          if (!grid_idx_out_of_bounds) {
            LSM_FMM_IDX(idx_neighbor2, neighbor2, grid_dims);
            neighbor_status = (PointStatus) gridpoint_status[idx_neighbor2];
            if ( (KNOWN == neighbor_status) &&
                 (  LSM_FMM_ABS(phi[idx_neighbor2]) 
                 <= LSM_FMM_ABS(phi_upwind1)) ) {
              phi_upwind2 = phi[idx_neighbor2];
              second_order_switch = 1;
            }
          } 

        }  /* end if: plus is upwind direction */
      } /* end case: first-order neighbor is KNOWN */
    }

    /*
     * accumulate coefficients for phi if any of the neighbors are "KNOWN"
     */
    if (phi_upwind1 < DBL_MAX) {
      /* temporary variables */
      double one_plus_switch_over_two = 1.0+0.5*second_order_switch;
      double phi_upwind_contrib;

      /* set phi_upwind_contrib to be first- or second-order */
      /* contribution based on value of second_order_switch  */
      if (second_order_switch == 1) {
        phi_upwind_contrib = phi_upwind1*(1+second_order_switch)
                           - 0.5*second_order_switch*phi_upwind2;
      } else {
        phi_upwind_contrib = phi_upwind1;
      }

      /* accumulate coefs for phi */ 
      inv_dx_sq = 1/dx[dir]; inv_dx_sq *= inv_dx_sq; 
      phi_A += inv_dx_sq*one_plus_switch_over_two*one_plus_switch_over_two;
      phi_B += inv_dx_sq*one_plus_switch_over_two*phi_upwind_contrib;
      phi_C += inv_dx_sq*phi_upwind_contrib*phi_upwind_contrib;
    }

  } /* loop over coordinate directions */

  /* complete computation of phi_B and phi_C */
  phi_B *= -2.0;
  phi_C -= 1/speed[idx_cur_gridpoint]/speed[idx_cur_gridpoint];

  /* check that phi_A is nonzero */
  if (LSM_FMM_ABS(phi_A) == 0) {
    fprintf(stderr,"ERROR: phi update - no KNOWN neighbors!!!\n");
    fprintf(stderr,"       phi set to 'infinity'.\n");
    return DBL_MAX;
  }

  /* compute phi by solving quadratic equation */
  discriminant = phi_B*phi_B - 4.0*phi_A*phi_C;
  phi_updated = DBL_MAX;
  max_dx = dx[0];
  for (dir = 1; dir < FMM_NDIM; dir++) {
    max_dx = (max_dx > dx[dir]) ? max_dx : dx[dir];
  }
  if (discriminant >= 0) {
    phi_updated = 0.5*(-phi_B + sqrt(discriminant))/phi_A;
  } else if (discriminant >= -4.0*max_dx*max_dx*phi_A*phi_A) {
      phi_updated = -0.5*phi_B/phi_A;
      fprintf(stderr,"WARNING: phi update - discriminant slightly negative!!!\n");
      fprintf(stderr,"         phi updated with O(dx) error.\n");
  } else {
    fprintf(stderr,"ERROR: phi update - discriminant negative!!!\n");
    fprintf(stderr,"       phi set to 'infinity'.\n");
    fprintf(stderr,"       discriminant = %g,", discriminant);
    fprintf(stderr," grid_idx = (");
    for (l = 0; l < FMM_NDIM-1; l++) { 
      fprintf(stderr,"%d,", grid_idx[l]);
    }
    fprintf(stderr,"%d)\n",grid_idx[l]);
  }

  /* set phi at current grid point */
  phi[idx_cur_gridpoint] = phi_updated;

  return phi_updated;
}

#endif
