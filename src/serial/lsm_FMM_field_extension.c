/*
 * File:        lsm_FMM_field_extension.c
 * Copyright:   (c) 2005-2008 Kevin T. Chu
 * Revision:    $Revision: 1.2 $
 * Modified:    $Date: 2006/08/13 13:30:52 $
 * Description: Implementation of Fast Marching Method for computing
 *              signed distance functions and extension fields
 */

/*! \file lsm_FMM_field_extension.c
 *
 * \brief
 * @ref lsm_FMM_field_extension.c provides "generic" implementations 
 *      first- and second-order accurate Fast Marching Method schemes 
 *      for computing signed distance functions and extension fields.  
 *      The code is "templated" on  the number of dimensions through the 
 *      use of macro definitions that MUST be provided by the user.  
 *
 *
 * <h3> Usage: </h3>
 *
 * -# Define the following macros:
 *    -# FMM_NDIM:  the number of spatial dimensions.
 *    -# FMM_COMPUTE_DISTANCE_FUNCTION:  desired name of function
 *       that computes the distance function given a level set
 *       function
 *    -# FMM_COMPUTE_EXTENSION_FIELDS:  desired name of function
 *       that computes the extensions of fields off of the zero 
 *       level set 
 *    -# FMM_INITIALIZE_FRONT_ORDER1:  desired name of function that
 *       initializes the values on the front using a first-order scheme
 *    -# FMM_INITIALIZE_FRONT_ORDER2:  desired name of function that
 *       initializes the values on the front using a second-order scheme
 *    -# FMM_UPDATE_GRID_POINT_ORDER1:  desired name of function
 *       that updates the value of the solution at grid points using
 *       a first-order accurate discretization
 *    -# FMM_UPDATE_GRID_POINT_ORDER2:  desired name of function
 *       that updates the value of the solution at grid points using
 *       a second-order accurate discretization
 * -# Include this file at the end of the implementation file
 *    for the n-dimentsional Eikonal equation solver.
 * -# Compile code.
 *
 *
 * <h3> NOTES: </h3>
 * - This implementation currently does NOT support second-order
 *   initialization grid points adjacent to the zero level set.
 *   As a consequence, the computed signed distance function and
 *   extension fields are only first-order accurate in the L-infinity
 *   norm.  All results should, however, be second-order accurate in 
 *   the L2 norm.
 * 
 * - Because this code depends on macros, care must be taken to
 *   ensure that macros do not conflict.
 *
 */

#ifndef included_lsm_FMM_field_extension_c
#define included_lsm_FMM_field_extension_c

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
#error "lsm_FMM_field_extension: required macro FMM_NDIM not defined!"
#endif FMM_NDIM
#ifndef FMM_COMPUTE_DISTANCE_FUNCTION
#error "lsm_FMM_field_extension: required macro FMM_COMPUTE_DISTANCE_FUNCTION not defined!"
#endif FMM_COMPUTE_DISTANCE_FUNCTION
#ifndef FMM_COMPUTE_EXTENSION_FIELDS
#error "lsm_FMM_field_extension: required macro FMM_COMPUTE_EXTENSION_FIELDS not defined!"
#endif FMM_COMPUTE_EXTENSION_FIELDS
#ifndef FMM_INITIALIZE_FRONT_ORDER1
#error "lsm_FMM_field_extension: required macro FMM_INITIALIZE_FRONT_ORDER1 not defined!"
#endif FMM_INITIALIZE_FRONT_ORDER1
#ifndef FMM_INITIALIZE_FRONT_ORDER2
#error "lsm_FMM_field_extension: required macro FMM_INITIALIZE_FRONT_ORDER2 not defined!"
#endif FMM_INITIALIZE_FRONT_ORDER2
#ifndef FMM_UPDATE_GRID_POINT_ORDER1
#error "lsm_FMM_field_extension: required macro FMM_UPDATE_GRID_POINT_ORDER1 not defined!"
#endif FMM_UPDATE_GRID_POINT_ORDER1
#ifndef FMM_UPDATE_GRID_POINT_ORDER2
#error "lsm_FMM_field_extension: required macro FMM_UPDATE_GRID_POINT_ORDER2 not defined!"
#endif FMM_UPDATE_GRID_POINT_ORDER2


/*=============== lsm_FMM_field_extension Data Structures =============*/
struct FMM_FieldData {
  double *phi;                 /* original level set function (input) */
  double *distance_function;   /* distance function (output)          */

  int num_extension_fields;    /* number of extension fields          */
  double **source_fields;      /* source fields to extend off of zero */
                               /* level set (input)                   */
  double **extension_fields;   /* computed extension field (output)   */
};


/*============================ FMM Functions ===========================*/

/*
 * FMM_INITIALIZE_FRONT_ORDER1() implements the callback function 
 * required by FMM_Core::FMM_initializeFront() to find and initialize 
 * the front.  
 *
 * The approximation to the distance function is computed using
 * a first-order scheme.  The values of the extension fields on 
 * the zero level set are computed using linear interpolation.  
 * The extension fields for points adjacent to the zero level set
 * are calculated using a first-order approximation to the 
 * grad(F)*grad(dist) = 0 equation.
 */
void FMM_INITIALIZE_FRONT_ORDER1(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int num_dims,
  int *grid_dims,
  double *dx);

/*
 * FMM_INITIALIZE_FRONT_ORDER2() implements the callback function 
 * required by FMM_Core::FMM_initializeFront() to find and initialize 
 * the front.  
 *
 * The approximation to the distance function is computed using
 * a first-order scheme.  The values of the extension fields on 
 * the zero level set are computed using quadratic interpolation.  
 * The extension fields for points adjacent to the zero level set
 * are calculated using a first-order approximation to the 
 * grad(F)*grad(dist) = 0 equation.
 */
void FMM_INITIALIZE_FRONT_ORDER2(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int num_dims,
  int *grid_dims,
  double *dx);

/* 
 * FMM_UPDATE_GRID_POINT_ORDER1() implements the callback function 
 * required by FMM_Core::FMM_Core_updateNeighbors() to update the
 * distance function and extension fields at a grid point.  It 
 * computes and returns the updated distance function and extension
 * field values of the specified grid point using values of neighbors
 * that have status "KNOWN".
 *
 * The approximation to the distance function is computed
 * using a first-order finite-difference scheme.  The
 * extension fields are calculated using a first-order
 * approximation to the grad(F)*grad(dist) = 0 equation.
 */
double FMM_UPDATE_GRID_POINT_ORDER1(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int *grid_idx,
  int num_dims,
  int *grid_dims,
  double *dx);

/* 
 * FMM_UPDATE_GRID_POINT_ORDER2() implements the callback function 
 * required by FMM_Core::FMM_Core_updateNeighbors() to update the
 * distance function and extension fields at a grid point.  It 
 * computes and returns the updated distance function and extension
 * field values of the specified grid point using values of neighbors
 * that have status "KNOWN".
 *
 * The approximation to the distance function is computed
 * using a second-order finite-difference scheme.  
 * The extension fields are calculated using a second-order
 * approximation to the grad(F)*grad(dist) = 0 equation.
 */
double FMM_UPDATE_GRID_POINT_ORDER2(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int *grid_idx,
  int num_dims,
  int *grid_dims,
  double *dx);


/*==================== Function Definitions =========================*/


int FMM_COMPUTE_EXTENSION_FIELDS(
  double *distance_function,
  double **extension_fields,
  double *phi,
  double *mask,
  double **source_fields,
  int num_extension_fields,
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
  int i, j, idx;            /* loop variables */
  double *ptr;              /* pointer to field data */


  /******************************************************
   * set up appropriate grid point update and front
   * detection/initialization functions based on the
   * specified spatial derivative order
   ******************************************************/
  if (spatial_discretization_order == 1) {

    initializeFront = &FMM_INITIALIZE_FRONT_ORDER1;
    updateGridPoint = &FMM_UPDATE_GRID_POINT_ORDER1;

  } else if (spatial_discretization_order == 2) {

    /* NOTE: we use first-order initialization of values      */
    /*       on the front because higher-order initialization */
    /*       is not fully implemented yet.                    */
    initializeFront = &FMM_INITIALIZE_FRONT_ORDER1; 
    updateGridPoint = &FMM_UPDATE_GRID_POINT_ORDER2;

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
  fmm_field_data->phi = phi;
  fmm_field_data->distance_function = distance_function;
  fmm_field_data->num_extension_fields = num_extension_fields;
  fmm_field_data->source_fields = source_fields;
  fmm_field_data->extension_fields = extension_fields;
   
  /********************************************
   * initialize phi and extension fields
   ********************************************/
  num_gridpoints = 1;
  for (i = 0; i < FMM_NDIM; i++) {
    num_gridpoints *= grid_dims[i];
  }
  for (i = 0, ptr = distance_function; i < num_gridpoints; i++, ptr++) {
    *ptr = LSM_FMM_DEFAULT_UPDATE_VALUE;
  }

  for (j = 0; j < num_extension_fields; j++) {
    for (i = 0, ptr = extension_fields[j]; i < num_gridpoints; i++, ptr++) {
      *ptr = LSM_FMM_DEFAULT_UPDATE_VALUE;
    }
  }

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

  /* mark grid points outside of domain */
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

    if ((mask) && (mask[idx] < 0)) {
      FMM_Core_markPointOutsideDomain(fmm_core_data, grid_idx);

      /* set distance_function and extension fields to DBL_MAX */
      distance_function[idx] = DBL_MAX;
      for (i = 0; i < num_extension_fields; i++) {
        extension_fields[i][idx] = DBL_MAX;
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

  return LSM_FMM_ERR_SUCCESS;
}

/* 
 * FMM_COMPUTE_DISTANCE_FUNCTION() just calls FMM_COMPUTE_EXTENSION_FIELDS()
 * with no source/extension fields (i.e. NULL source/extension field
 * pointers).
 */
int FMM_COMPUTE_DISTANCE_FUNCTION(
  double *distance_function,
  double *phi,
  double *mask,
  int spatial_discretization_order,
  int *grid_dims,
  double *dx)
{
  return FMM_COMPUTE_EXTENSION_FIELDS(
           distance_function,
           0, /*  NULL extension fields pointer */
           phi,
           mask,
           0, /*  NULL source fields pointer */
           0, /*  zero extension fields to compute */
           spatial_discretization_order,
           grid_dims,
           dx);
}

void FMM_INITIALIZE_FRONT_ORDER1(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int num_dims,
  int *grid_dims,
  double *dx)
{
  int *gridpoint_status = FMM_Core_getGridPointStatusDataArray(fmm_core_data);

  /* FMM Field Data variables */
  double *phi = fmm_field_data->phi;
  double *distance_function = fmm_field_data->distance_function; 
  int num_extension_fields = fmm_field_data->num_extension_fields; 
  double **source_fields = fmm_field_data->source_fields; 
  double **extension_fields = fmm_field_data->extension_fields; 
   
  /* grid variables */
  int offset[FMM_NDIM];
  int neighbor[FMM_NDIM];
  int on_interface = LSM_FMM_FALSE;
  int borders_interface = LSM_FMM_FALSE;

  /* distance function variables */
  double phi_cur;
  double phi_minus;
  double phi_plus;
  double dist_minus;
  double dist_plus;
  double dist_dir;
  int use_plus; 

  double sum_dist_inv_sq; 
  double dist_inv_sq_dir; 

  /* variables for extension field calculations */
  double *extension_fields_cur;
  double *extension_fields_sum_div_dist_sq;
  double *extension_fields_minus;
  double *extension_fields_plus;

  /* auxilliary variables */
  int num_gridpoints;       /* number of grid points */ 
  int i,idx;  /* loop variables for grid */
  int idx_neighbor;
  int m;    /* loop variable for extension fields */
  int l;    /* extra loop variable */
  int dir;  /* loop variable over spatial dimensions */
  int grid_idx_out_of_bounds;  /* flag indicating whether grid_idx */
                               /* is out of bounds                 */


  /* unused function parameters */
  (void) num_dims;


  /* allocate memory for extension field calculations */
  extension_fields_cur = (double*) malloc(num_extension_fields*sizeof(double));
  extension_fields_sum_div_dist_sq = 
    (double*) malloc(num_extension_fields*sizeof(double));
  extension_fields_minus = 
    (double*) malloc(num_extension_fields*sizeof(double));
  extension_fields_plus = 
    (double*) malloc(num_extension_fields*sizeof(double));

  /*
   * loop through cells in grid to find those that border the
   * zero level set.
   */
  num_gridpoints = 1;
  for (i = 0; i < FMM_NDIM; i++) {
    num_gridpoints *= grid_dims[i];
  }

  for (idx = 0; idx < num_gridpoints; idx++) {

    /* temporary variables */
    int grid_idx[FMM_NDIM];
    int idx_remainder; 

    /* skip point if it is out of the mathematical/physical domain */
    if (OUTSIDE_DOMAIN == gridpoint_status[idx]) {
      continue;
    }

    /* compute grid_idx */
    idx_remainder = idx;
    for (i = 0; i < FMM_NDIM; i++) {
      grid_idx[i] = idx_remainder%grid_dims[i];
      idx_remainder -= grid_idx[i];
      idx_remainder /= grid_dims[i];
    }

    /* initialize on_interface and borders_interface to FALSE */
    on_interface = LSM_FMM_FALSE;
    borders_interface = LSM_FMM_FALSE;

    /* get data values at the current grid point */
    phi_cur = phi[idx];
    for (m = 0; m < num_extension_fields; m++) {
      extension_fields_cur[m] = source_fields[m][idx];
    }

    /* zero out accumulation variables */
    sum_dist_inv_sq = 0; 
    for (m = 0; m < num_extension_fields; m++) {
      extension_fields_sum_div_dist_sq[m] = 0;
    }
  
    /* loop over neighbors */
    for (dir = 0; dir < FMM_NDIM; dir++) {

      /* check if current grid point is on the zero level set. */
      /* if so, there is no need to check other directions,    */
      /* so break out of loop.                                 */
      if (LSM_FMM_ABS(phi_cur) < LSM_FMM_ZERO_TOL) {
        on_interface = LSM_FMM_TRUE; 
        break;
      }

      for (l = 0; l < FMM_NDIM; l++) { /* reset offset */
        offset[l] = 0; 
      }

      /* reset plus and minus distances */
      dist_plus = DBL_MAX;
      dist_minus = DBL_MAX;

      /* reset plus and minus extension field values */
      for (m = 0; m < num_extension_fields; m++) {
        extension_fields_minus[m] = 0;
        extension_fields_plus[m] = 0;
      }

      /* calculate distance to interface in minus direction */
      offset[dir] = -1;
      for (l = 0; l < FMM_NDIM; l++) {
        neighbor[l] = grid_idx[l] + offset[l];
      }
      LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,neighbor,grid_dims);
      if (!grid_idx_out_of_bounds) {
        LSM_FMM_IDX(idx_neighbor, neighbor, grid_dims);
        phi_minus = phi[idx_neighbor];
        if (phi_minus*phi_cur <= 0) {

          /* locate zero level set using linear interpolant */
          dist_minus = phi_cur/(phi_cur-phi_minus);

          /* use linear interpolation for value of source field */
          /* at interface                                       */
          for (m = 0; m < num_extension_fields; m++) {
            extension_fields_minus[m] = extension_fields_cur[m] 
              + dist_minus*(source_fields[m][idx_neighbor]
                           -extension_fields_cur[m]);
          }

          /* multiply back in the units for dist_minus */
          dist_minus *= dx[dir];
        }
      }
      
      /* calculate distance to interface in plus direction */
      offset[dir] = 1;
      for (l = 0; l < FMM_NDIM; l++) {
        neighbor[l] = grid_idx[l] + offset[l];
      }
      LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,neighbor,grid_dims);
      if (!grid_idx_out_of_bounds) {
        LSM_FMM_IDX(idx_neighbor, neighbor, grid_dims);
        phi_plus = phi[idx_neighbor];
        if (phi_plus*phi_cur <= 0) {

          /* locate zero level set using linear interpolant */
          dist_plus = phi_cur/(phi_cur-phi_plus);

          /* use linear interpolation for value of source field */
          /* at interface                                       */
          for (m = 0; m < num_extension_fields; m++) {
            extension_fields_plus[m] = extension_fields_cur[m] 
              + dist_plus*(source_fields[m][idx_neighbor]
                          -extension_fields_cur[m]);
          }

          /* multiply back in the units for dist_plus */
          dist_plus *= dx[dir];
        }
      }

      /* check if current grid point borders the interface */
      if ( (dist_plus < DBL_MAX) || (dist_minus < DBL_MAX) ) {

        borders_interface = LSM_FMM_TRUE; 

        if (dist_plus < dist_minus) {
          dist_dir = dist_plus;
          use_plus = LSM_FMM_TRUE;
        } else {
          dist_dir = dist_minus;
          use_plus = LSM_FMM_FALSE;
        }

        /* update 1/dist^2 and ext_field/dist^2 values with */
        /* information from current coordinate direction    */
        dist_inv_sq_dir = 1/dist_dir/dist_dir;
        sum_dist_inv_sq += dist_inv_sq_dir; 

        if (use_plus) {
          for (m = 0; m < num_extension_fields; m++) {
            extension_fields_sum_div_dist_sq[m] += 
              extension_fields_plus[m]*dist_inv_sq_dir;
          }
        } else {
          for (m = 0; m < num_extension_fields; m++) {
            extension_fields_sum_div_dist_sq[m] += 
              extension_fields_minus[m]*dist_inv_sq_dir;
          }
        }

      } /* end cases: current grid point on or borders interface */

    } /* end loop over neighbors */

    /* set distance function and extension field of grid points */
    /* on or bordering the zero level set                       */
    if (on_interface) { 

      /* Set distance function.  We ensure the correct sign */
      /* by introducing a very small error.                 */
      if (phi_cur > 0) 
        distance_function[idx] = DBL_MIN;
      else
        distance_function[idx] = -DBL_MIN;

      /* compute extension field value */
      for (m = 0; m < num_extension_fields; m++) {
        extension_fields[m][idx] = extension_fields_cur[m];
      }

      /* set grid point as an initial front point */
      FMM_Core_setInitialFrontPoint(fmm_core_data, grid_idx,
                                    distance_function[idx]);

    } else if (borders_interface) { 

      /* compute updated value for the signed distance function */
      if (phi_cur > 0)
        distance_function[idx] = sqrt(1.0/sum_dist_inv_sq);
      else 
        distance_function[idx] = -sqrt(1.0/sum_dist_inv_sq);

      /* compute extension field value */
      for (m = 0; m < num_extension_fields; m++) {
        extension_fields[m][idx] = 
          extension_fields_sum_div_dist_sq[m]/sum_dist_inv_sq;
      }

      /* set grid point as an initial front point */
      FMM_Core_setInitialFrontPoint(fmm_core_data, grid_idx,
                                    distance_function[idx]);

    } /* end handling grid points on or near interface */

  }  /* end loop over grid */

  /* clean up memory */
  free(extension_fields_cur); 
  free(extension_fields_sum_div_dist_sq); 
  free(extension_fields_minus);
  free(extension_fields_plus);
}


void FMM_INITIALIZE_FRONT_ORDER2(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int num_dims,
  int *grid_dims,
  double *dx)
{
  int *gridpoint_status = FMM_Core_getGridPointStatusDataArray(fmm_core_data);

  /* FMM Field Data variables */
  double *phi = fmm_field_data->phi;
  double *distance_function = fmm_field_data->distance_function; 
  int num_extension_fields = fmm_field_data->num_extension_fields; 
  double **source_fields = fmm_field_data->source_fields; 
  double **extension_fields = fmm_field_data->extension_fields; 
   
  /* grid variables */
  int neighbor_plus[FMM_NDIM], neighbor_minus[FMM_NDIM];
  int on_interface = LSM_FMM_FALSE;
  int borders_interface = LSM_FMM_FALSE;

  /* distance function variables */
  double phi_cur;
  double phi_minus;
  double phi_plus;
  double dist_minus;
  double dist_plus;
  double dist_dir;
  int use_plus; 

  double sum_dist_inv_sq; 
  double dist_inv_sq_dir; 

  /* variables for extension field calculations */
  double *extension_fields_cur;
  double *extension_fields_sum_div_dist_sq;
  double *extension_fields_minus;
  double *extension_fields_plus;

  /* auxilliary variables */
  int num_gridpoints;       /* number of grid points */ 
  int i,idx;  /* loop variables for grid */
  int idx_neighbor_plus, idx_neighbor_minus;
  int m;    /* loop variable for extension fields */
  int l;    /* extra loop variable */
  int dir;  /* loop variable over spatial dimensions */
  int grid_idx_out_of_bounds;  /* flag indicating whether grid_idx */
                               /* is out of bounds                 */


  /* unused function parameters */
  (void) num_dims;


  /* allocate memory for extension field calculations */
  extension_fields_cur = (double*) malloc(num_extension_fields*sizeof(double));
  extension_fields_sum_div_dist_sq = 
    (double*) malloc(num_extension_fields*sizeof(double));
  extension_fields_minus = 
    (double*) malloc(num_extension_fields*sizeof(double));
  extension_fields_plus = 
    (double*) malloc(num_extension_fields*sizeof(double));

  /*
   * loop through cells in grid to find those that border the
   * zero level set.
   */
  num_gridpoints = 1;
  for (i = 0; i < FMM_NDIM; i++) {
    num_gridpoints *= grid_dims[i];
  }

  for (idx = 0; idx < num_gridpoints; idx++) {

    /* temporary variables */
    int grid_idx[FMM_NDIM];
    int idx_remainder; 

    /* skip point if it is out of the mathematical/physical domain */
    if (OUTSIDE_DOMAIN == gridpoint_status[idx]) {
      continue;
    }

    /* compute grid_idx */
    idx_remainder = idx;
    for (i = 0; i < FMM_NDIM; i++) {
      grid_idx[i] = idx_remainder%grid_dims[i];
      idx_remainder -= grid_idx[i];
      idx_remainder /= grid_dims[i];
    }

    /* initialize on_interface and borders_interface to FALSE */
    on_interface = LSM_FMM_FALSE;
    borders_interface = LSM_FMM_FALSE;

    /* get data values at the current grid point */
    phi_cur = phi[idx];
    for (m = 0; m < num_extension_fields; m++) {
      extension_fields_cur[m] = source_fields[m][idx];
    }

    /* zero out accumulation variables */
    sum_dist_inv_sq = 0; 
    for (m = 0; m < num_extension_fields; m++) {
      extension_fields_sum_div_dist_sq[m] = 0;
    }
  
    /* loop over neighbors */
    for (dir = 0; dir < FMM_NDIM; dir++) {

      /* check if current grid point is on the zero level set. */
      /* if so, there is no need to check other directions,    */
      /* so break out of loop.                                 */
      if (LSM_FMM_ABS(phi_cur) < LSM_FMM_ZERO_TOL) {
        on_interface = LSM_FMM_TRUE; 
        break;
      }

      /* reset plus and minus distances */
      dist_plus = DBL_MAX;
      dist_minus = DBL_MAX;

      /* reset plus and minus extension field values */
      for (m = 0; m < num_extension_fields; m++) {
        extension_fields_minus[m] = 0;
        extension_fields_plus[m] = 0;
      }

      /* calculate neighbor indices */
      for (l = 0; l < FMM_NDIM; l++) {
        neighbor_plus[l] = grid_idx[l];
        neighbor_minus[l] = grid_idx[l];
      }
      neighbor_plus[dir]++; neighbor_minus[dir]--;
      LSM_FMM_IDX(idx_neighbor_plus, neighbor_plus, grid_dims);
      LSM_FMM_IDX(idx_neighbor_minus, neighbor_minus, grid_dims);

      /* calculate distance to interface in minus direction */
      LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,
                                neighbor_minus,grid_dims);
      if (!grid_idx_out_of_bounds) {
        phi_minus = phi[idx_neighbor_minus];
        if (phi_minus*phi_cur <= 0) {

          /* locate zero level set using linear interpolant */
          dist_minus = phi_cur/(phi_cur-phi_minus);

          /* compute value of extension fields on zero level set */
          LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,
                                    neighbor_plus,grid_dims);
          if (!grid_idx_out_of_bounds) {

            /* use quadratic interpolation for value of source field */
            /* at interface                                          */
            for (m = 0; m < num_extension_fields; m++) {
              extension_fields_minus[m] = 
                  0.5*dist_minus*(dist_minus+1.0)
                                *source_fields[m][idx_neighbor_minus]
                - (dist_minus-1.0)*(dist_minus+1.0)
                                  *extension_fields_cur[m] 
                + 0.5*dist_minus*(dist_minus-1.0)
                                *source_fields[m][idx_neighbor_plus];
            }

          } else {

            /* drop to linear interpolation for value of source field */
            /* at interface since plus neighbor is out of bounds      */
            for (m = 0; m < num_extension_fields; m++) {
              extension_fields_minus[m] = extension_fields_cur[m] 
                + dist_minus*(source_fields[m][idx_neighbor_minus]
                             -extension_fields_cur[m]);
            }

          }

          /* multiply back in the units for dist_minus */
          dist_minus *= dx[dir];
        }
      }
      
      /* calculate distance to interface in plus direction */
      LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,
                                neighbor_plus,grid_dims);
      if (!grid_idx_out_of_bounds) {
        phi_plus = phi[idx_neighbor_plus];
        if (phi_plus*phi_cur <= 0) {

          /* locate zero level set using linear interpolant */
          dist_plus = phi_cur/(phi_cur-phi_plus);

          /* compute value of extension fields on zero level set */
          LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,
                                    neighbor_minus,grid_dims);
          if (!grid_idx_out_of_bounds) {

            /* use quadratic interpolation for value of source field */
            /* at interface                                          */
            for (m = 0; m < num_extension_fields; m++) {
              extension_fields_plus[m] = 
                  0.5*dist_plus*(dist_plus-1.0)
                               *source_fields[m][idx_neighbor_minus]
                - (dist_plus-1.0)*(dist_plus+1.0)
                                  *extension_fields_cur[m] 
                + 0.5*dist_plus*(dist_plus+1.0)
                               *source_fields[m][idx_neighbor_plus];
            }

          } else {

            /* drop to linear interpolation for value of source field */
            /* at interface since minus neighbor is out of bounds     */
            for (m = 0; m < num_extension_fields; m++) {
              extension_fields_plus[m] = extension_fields_cur[m] 
                + dist_plus*(source_fields[m][idx_neighbor_plus]
                            -extension_fields_cur[m]);
            }

          }

          /* multiply back in the units for dist_plus */
          dist_plus *= dx[dir];
        }
      }

      /* check if current grid point borders the interface */
      if ( (dist_plus < DBL_MAX) || (dist_minus < DBL_MAX) ) {

        borders_interface = LSM_FMM_TRUE; 

        if (dist_plus < dist_minus) {
          dist_dir = dist_plus;
          use_plus = LSM_FMM_TRUE;
        } else {
          dist_dir = dist_minus;
          use_plus = LSM_FMM_FALSE;
        }

        /* update 1/dist^2 and ext_field/dist^2 values with */
        /* information from current coordinate direction    */
        dist_inv_sq_dir = 1/dist_dir/dist_dir;
        sum_dist_inv_sq += dist_inv_sq_dir; 

        if (use_plus) {
          for (m = 0; m < num_extension_fields; m++) {
            extension_fields_sum_div_dist_sq[m] += 
              extension_fields_plus[m]*dist_inv_sq_dir;
          }
        } else {
          for (m = 0; m < num_extension_fields; m++) {
            extension_fields_sum_div_dist_sq[m] += 
              extension_fields_minus[m]*dist_inv_sq_dir;
          }
        }

      } /* end cases: current grid point on or borders interface */

    } /* end loop over neighbors */

    /* set distance function and extension field of grid points */
    /* on or bordering the zero level set                       */
    if (on_interface) { 

      /* Set distance function.  We ensure the correct sign */
      /* by introducing a very small error.                 */
      if (phi_cur > 0) 
        distance_function[idx] = DBL_MIN;
      else
        distance_function[idx] = -DBL_MIN;

      /* compute extension field value */
      for (m = 0; m < num_extension_fields; m++) {
        extension_fields[m][idx] = extension_fields_cur[m];
      }

      /* set grid point as an initial front point */
      FMM_Core_setInitialFrontPoint(fmm_core_data, grid_idx,
                                    distance_function[idx]);

    } else if (borders_interface) { 

      /* compute updated value for the signed distance function */
      if (phi_cur > 0)
        distance_function[idx] = sqrt(1.0/sum_dist_inv_sq);
      else 
        distance_function[idx] = -sqrt(1.0/sum_dist_inv_sq);

      /* compute extension field value */
      for (m = 0; m < num_extension_fields; m++) {
        extension_fields[m][idx] = 
          extension_fields_sum_div_dist_sq[m]/sum_dist_inv_sq;
      }

      /* set grid point as an initial front point */
      FMM_Core_setInitialFrontPoint(fmm_core_data, grid_idx,
                                    distance_function[idx]);

    } /* end handling grid points on or near interface */

  }  /* end loop over grid */

  /* clean up memory */
  free(extension_fields_cur); 
  free(extension_fields_sum_div_dist_sq); 
  free(extension_fields_minus);
  free(extension_fields_plus);
}


double FMM_UPDATE_GRID_POINT_ORDER1(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int *grid_idx,
  int num_dims,
  int *grid_dims,
  double *dx)
{
  int *gridpoint_status = FMM_Core_getGridPointStatusDataArray(fmm_core_data);

  /* FMM Field Data variables */
  double *distance_function = fmm_field_data->distance_function; 
  int num_extension_fields = fmm_field_data->num_extension_fields; 
  double **extension_fields = fmm_field_data->extension_fields; 

  /* variables for extension field calculations */
  double *extension_fields_numerator;
  double *extension_fields_denominator;

  /* variables used in distance function update */
  PointStatus  neighbor_status;
  int use_plus[FMM_NDIM];
  int dir_used[FMM_NDIM];
  double phi_upwind[FMM_NDIM];
  double phi_plus;
  double inv_dx_sq; 
  int offset[FMM_NDIM]; 
  int neighbor[FMM_NDIM];

  /* coefficients of quadratic equation for the updated distance function */
  double phi_A = 0;
  double phi_B = 0;
  double phi_C = 0;
  double discriminant;
  double dist_updated;

  /* auxilliary variables */
  int dir;  /* loop variable for spatial directions */
  int k;    /* loop variable for extension fields */
  int l;    /* extra loop variable */ 
  int idx_cur_gridpoint, idx_neighbor;
  int grid_idx_out_of_bounds;

  /* unused function parameters */
  (void) num_dims;

  /* allocate memory extension field calculations */
  extension_fields_numerator = 
    (double*) malloc(num_extension_fields*sizeof(double));
  extension_fields_denominator = 
    (double*) malloc(num_extension_fields*sizeof(double));;
  for (k = 0; k < num_extension_fields; k++) {
    extension_fields_numerator[k] = 0;
    extension_fields_denominator[k] = 0;
  }

  /* calculate update to distance function */
  for (dir = 0; dir < FMM_NDIM; dir++) { /* loop over coord directions */
    for (l = 0; l < FMM_NDIM; l++) { /* reset offset */
      offset[l] = 0; 
    }

    /* changed to true if has KNOWN neighbor */
    dir_used[dir] = LSM_FMM_FALSE;  

    /* find "upwind" direction and phi value */
    phi_upwind[dir] = DBL_MAX;

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
        phi_upwind[dir] = distance_function[idx_neighbor];
        use_plus[dir] = LSM_FMM_FALSE;
        dir_used[dir] = LSM_FMM_TRUE;
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
        phi_plus = distance_function[idx_neighbor];

        /* 
         * choosing the upwind direction to be the direction
         * with the smaller abs(phi) value gives a consistent 
         * solution to the "upwind" Eikonal equation.
         * NOTE: to avoid having to use the C math library,
         *       we define our own ABS macro
         */
        if (LSM_FMM_ABS(phi_plus) < LSM_FMM_ABS(phi_upwind[dir])) {
          phi_upwind[dir] = phi_plus;
          use_plus[dir] = LSM_FMM_TRUE;
          dir_used[dir] = LSM_FMM_TRUE;
        }
      }
    }

    /*
     * accumulate coefficients for updated distance function
     * if either of the neighbors are "known"
     */
    if (phi_upwind[dir] < DBL_MAX) {
      /* accumulate coefs for updated distance function */ 
      inv_dx_sq = 1/dx[dir]; inv_dx_sq *= inv_dx_sq;
      phi_A += inv_dx_sq;
      phi_B += inv_dx_sq*phi_upwind[dir];
      phi_C += inv_dx_sq*phi_upwind[dir]*phi_upwind[dir];
    }

  } /* loop over coordinate directions */

  /* check that phi_A is nonzero */
  if (LSM_FMM_ABS(phi_A) == 0) {
    fprintf(stderr,"ERROR: distance update - no KNOWN neighbors!!!\n");
    fprintf(stderr,"       distance set to 'infinity'.\n");
    return DBL_MAX;
  }

  /* complete computation of phi_B and phi_C */
  phi_B *= -2.0;
  phi_C -= 1.0;  

  /* compute updated distance function by solving quadratic equation */
  discriminant = phi_B*phi_B - 4.0*phi_A*phi_C;
  dist_updated = DBL_MAX;
  if (discriminant >= 0) {

    if (phi_B < 0) { /* distance function of neighbors is positive */
      dist_updated = 0.5*(-phi_B + sqrt(discriminant))/phi_A;
    } else if (phi_B > 0) { /* distance function of neighbors is negative */
      dist_updated = 0.5*(-phi_B - sqrt(discriminant))/phi_A;
    } else { 
      /* grid point is ON the interface, so keep it there */
      dist_updated = 0;
    } 

  } else {

    fprintf(stderr,"ERROR: distance update - discriminant negative!!!\n");

    if (phi_B < 0) { /* distance function of neighbors is positive */
      fprintf(stderr,"       distance set to 'infinity'.\n");
    } else { /* distance function of neighbors is negative */
      dist_updated *= -1;
      fprintf(stderr,"       distance set to '-infinity'.\n");
    }
    fprintf(stderr,"       discriminant = %g,", discriminant);
    fprintf(stderr," grid_idx = ("); 
    for (l = 0; l < FMM_NDIM-1; l++) { 
      fprintf(stderr,"%d,", grid_idx[l]);
    }
    fprintf(stderr,"%d)\n",grid_idx[l]);

  } /* end switch on value of discriminant */


  /* calculate extension field values */
  for (dir = 0; dir < FMM_NDIM; dir++) { /* loop over coord directions */

    /*
     * only accumulate values from the current direction if this
     * direction was used in the update of the distance function
     */
    if (dir_used[dir]) {
      for (l = 0; l < FMM_NDIM; l++) { /* reset offset */
        offset[l] = 0; 
      }
      offset[dir] = (use_plus[dir] ? 1 : -1);
      for (l = 0; l < FMM_NDIM; l++) {
        neighbor[l] = grid_idx[l] + offset[l];
      }
      LSM_FMM_IDX(idx_neighbor, neighbor, grid_dims);
  
      inv_dx_sq = 1/dx[dir]; inv_dx_sq *= inv_dx_sq;
      for (k = 0; k < num_extension_fields; k++) {
        double dist_diff = dist_updated - phi_upwind[dir];
        extension_fields_numerator[k] += 
          inv_dx_sq*dist_diff*extension_fields[k][idx_neighbor];
        extension_fields_denominator[k] += inv_dx_sq*dist_diff;
      }
    }
  } /* loop over coordinate directions */


  /* set updated quantities */
  LSM_FMM_IDX(idx_cur_gridpoint, grid_idx, grid_dims);
  distance_function[idx_cur_gridpoint] = dist_updated;
  for (k = 0; k < num_extension_fields; k++) {
    extension_fields[k][idx_cur_gridpoint] =
      extension_fields_numerator[k]/extension_fields_denominator[k];
  }

  /* free memory allocated for extension field calculations */
  free(extension_fields_numerator);
  free(extension_fields_denominator);

  return dist_updated;
}


double FMM_UPDATE_GRID_POINT_ORDER2(
  FMM_CoreData *fmm_core_data,
  FMM_FieldData *fmm_field_data,
  int *grid_idx,
  int num_dims,
  int *grid_dims,
  double *dx)
{
  int *gridpoint_status = FMM_Core_getGridPointStatusDataArray(fmm_core_data);

  /* FMM Field Data variables */
  double *distance_function = fmm_field_data->distance_function; 
  int num_extension_fields = fmm_field_data->num_extension_fields; 
  double **extension_fields = fmm_field_data->extension_fields; 

  /* variables for extension field calculations */
  double *extension_fields_numerator;
  double *extension_fields_denominator;

  /* variables used in distance function update */
  PointStatus  neighbor_status;
  int use_plus[FMM_NDIM];
  int dir_used[FMM_NDIM];
  double phi_upwind1[FMM_NDIM], phi_upwind2[FMM_NDIM];
  int second_order_switch[FMM_NDIM];
  double phi_plus;
  double inv_dx_sq; 
  int offset[FMM_NDIM]; 
  int neighbor1[FMM_NDIM];
  int neighbor2[FMM_NDIM];

  /* coefficients of quadratic equation for the updated distance function */
  double phi_A = 0;
  double phi_B = 0;
  double phi_C = 0;
  double discriminant;
  double dist_updated;
  double max_dx;

  /* auxilliary variables */
  int dir;  /* loop variable for spatial directions */
  int k;    /* loop variable for extension fields */
  int l;    /* extra loop variable */ 
  int idx_cur_gridpoint, idx_neighbor1, idx_neighbor2;
  int grid_idx_out_of_bounds;

  /* unused function parameters */
  (void) num_dims;

  /* allocate memory extension field calculations */
  extension_fields_numerator = 
    (double*) malloc(num_extension_fields*sizeof(double));
  extension_fields_denominator = 
    (double*) malloc(num_extension_fields*sizeof(double));;
  for (k = 0; k < num_extension_fields; k++) {
    extension_fields_numerator[k] = 0;
    extension_fields_denominator[k] = 0;
  }

  /* calculate update to distance function */
  for (dir = 0; dir < FMM_NDIM; dir++) { /* loop over coord directions */
    for (l = 0; l < FMM_NDIM; l++) { /* reset offset */
      offset[l] = 0; 
    }

    /* changed to true if has KNOWN neighbor */
    dir_used[dir] = LSM_FMM_FALSE;  

    /* reset phi_upwind1 and phi_upwind2 to DBL_MAX */
    phi_upwind1[dir] = DBL_MAX;
    phi_upwind2[dir] = DBL_MAX;

    /* reset second_order_switch to 0 (i.e. assume there are not enough */
    /* KNOWN neighbors for second-order discretization.                 */
    second_order_switch[dir] = 0;

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
        phi_upwind1[dir] = distance_function[idx_neighbor1];
        use_plus[dir] = LSM_FMM_FALSE;
        dir_used[dir] = LSM_FMM_TRUE;

        /* check for neighbor required for second-order accuracy */
        LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,neighbor2,grid_dims);
        if (!grid_idx_out_of_bounds) {
          LSM_FMM_IDX(idx_neighbor2, neighbor2, grid_dims);
          neighbor_status = (PointStatus) gridpoint_status[idx_neighbor2];
          if ( (KNOWN == neighbor_status) &&
               (  LSM_FMM_ABS(distance_function[idx_neighbor2])
               <= LSM_FMM_ABS(phi_upwind1[dir])) ) {
            phi_upwind2[dir] = distance_function[idx_neighbor2];
            second_order_switch[dir] = 1;
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
        phi_plus = distance_function[idx_neighbor1];

        /* 
         * choosing the upwind direction to be the direction
         * with the smaller abs(phi) value gives a consistent 
         * solution to the "upwind" Eikonal equation.
         * NOTE: to avoid having to use the C math library,
         *       we define our own ABS macro
         */
        if (LSM_FMM_ABS(phi_plus) < LSM_FMM_ABS(phi_upwind1[dir])) {
          phi_upwind1[dir] = phi_plus;
          phi_upwind2[dir] = DBL_MAX;
          second_order_switch[dir] = 0;
          use_plus[dir] = LSM_FMM_TRUE;
          dir_used[dir] = LSM_FMM_TRUE;
          
          /* check for neighbor required for second-order accuracy */
          LSM_FMM_IDX_OUT_OF_BOUNDS(grid_idx_out_of_bounds,neighbor2,grid_dims);
          if (!grid_idx_out_of_bounds) {
            LSM_FMM_IDX(idx_neighbor2, neighbor2, grid_dims);
            neighbor_status = (PointStatus) gridpoint_status[idx_neighbor2];
            if ( (KNOWN == neighbor_status) &&
                 (  LSM_FMM_ABS(distance_function[idx_neighbor2])
                 <= LSM_FMM_ABS(phi_upwind1[dir])) ) {
              phi_upwind2[dir] = distance_function[idx_neighbor2];
              second_order_switch[dir] = 1;
            }
          }
        
        }  /* end if: plus is upwind direction */
      } /* end case: first-order neighbor is KNOWN */

    }

    /* 
     * accumulate coefficients for phi if any of the neighbors are "KNOWN"
     */
    if (phi_upwind1[dir] < DBL_MAX) {
      /* temporary variables */
      double one_plus_switch_over_two = 1.0+0.5*second_order_switch[dir];
      double phi_upwind_contrib;
      
      /* set phi_upwind_contrib to be first- or second-order */
      /* contribution based on value of second_order_switch  */
      if (second_order_switch[dir] == 1) {
        phi_upwind_contrib = 2.0*phi_upwind1[dir] - 0.5*phi_upwind2[dir];
      } else {
        phi_upwind_contrib = phi_upwind1[dir];
      }
      
      /* accumulate coefs for phi */
      inv_dx_sq = 1/dx[dir]; inv_dx_sq *= inv_dx_sq;
      phi_A += inv_dx_sq*one_plus_switch_over_two*one_plus_switch_over_two;
      phi_B += inv_dx_sq*one_plus_switch_over_two*phi_upwind_contrib;
      phi_C += inv_dx_sq*phi_upwind_contrib*phi_upwind_contrib;
    }

  } /* loop over coordinate directions */

  /* check that phi_A is nonzero */
  if (LSM_FMM_ABS(phi_A) == 0) {
    fprintf(stderr,"ERROR: distance update - no KNOWN neighbors!!!\n");
    fprintf(stderr,"       distance set to 'infinity'.\n");
    return DBL_MAX;
  }

  /* complete computation of phi_B and phi_C */
  phi_B *= -2.0;
  phi_C -= 1.0;  

  /* compute updated distance function by solving quadratic equation */
  discriminant = phi_B*phi_B - 4.0*phi_A*phi_C;
  dist_updated = DBL_MAX;
  max_dx = dx[0];
  for (dir = 1; dir < FMM_NDIM; dir++) {
    max_dx = (max_dx > dx[dir]) ? max_dx : dx[dir];
  }
  if (discriminant >= 0) {

    if (phi_B < 0) { /* distance function of neighbors is positive */
      dist_updated = 0.5*(-phi_B + sqrt(discriminant))/phi_A;
    } else if (phi_B > 0) { /* distance function of neighbors is negative */
      dist_updated = 0.5*(-phi_B - sqrt(discriminant))/phi_A;
    } else { 
      /* grid point is ON the interface, so keep it there */
      dist_updated = 0;
    } 

  } else if (discriminant >= -4.0*max_dx*max_dx*phi_A*phi_A) { 

      dist_updated = -0.5*phi_B/phi_A;

      /* KTC - comment back in when error reporting is set up */
      /*
      fprintf(stderr,"WARNING: distance update - discriminant slightly negative!!!\n");    
      fprintf(stderr,"         distance updated with O(dx) error.\n");
      */

  } else {

    fprintf(stderr,"ERROR: distance update - discriminant negative!!!\n");

    if (phi_B < 0) { /* distance function of neighbors is positive */
      fprintf(stderr,"       distance set to 'infinity'.\n");
    } else { /* distance function of neighbors is negative */
      dist_updated *= -1;
      fprintf(stderr,"       distance set to '-infinity'.\n");
    }

    fprintf(stderr,"       discriminant = %g,", discriminant);
    fprintf(stderr," grid_idx = ("); 
    for (l = 0; l < FMM_NDIM-1; l++) { 
      fprintf(stderr,"%d,", grid_idx[l]);
    }
    fprintf(stderr,"%d)\n",grid_idx[l]);

  } /* end switch on value of discriminant */


  /* calculate extension field values */
  for (dir = 0; dir < FMM_NDIM; dir++) { /* loop over coord directions */

    /*
     * only accumulate values from the current direction if this
     * direction was used in the update of the distance function
     */
    if (dir_used[dir]) {
      for (l = 0; l < FMM_NDIM; l++) { /* reset offset */
        offset[l] = 0; 
      }
      offset[dir] = (use_plus[dir] == LSM_FMM_TRUE ? 1 : -1);
      for (l = 0; l < FMM_NDIM; l++) {
        neighbor1[l] = grid_idx[l] + offset[l];
        neighbor2[l] = grid_idx[l] + 2*offset[l];
      }
      LSM_FMM_IDX(idx_neighbor1, neighbor1, grid_dims);
      LSM_FMM_IDX(idx_neighbor2, neighbor2, grid_dims);
  
      inv_dx_sq = 1/dx[dir]; inv_dx_sq *= inv_dx_sq;
      if (second_order_switch[dir] == 1) {

        double grad_dist = 1.5*dist_updated - 2.0*phi_upwind1[dir] 
                         + 0.5*phi_upwind2[dir]; 

        for (k = 0; k < num_extension_fields; k++) {

        /* KTC - second-order discretization seems to lead to  */
        /*       larger errors than first-order discretization */
        /*       Currently using first-order discretization.   */
        /*       FIX ME!!!                                     */
        /*
          extension_fields_numerator[k] +=  
             inv_dx_sq*grad_dist
            *( 2.0*extension_fields[k][idx_neighbor1]
             - 0.5*extension_fields[k][idx_neighbor2] );
          extension_fields_denominator[k] += 1.5*inv_dx_sq*grad_dist;
         */

          extension_fields_numerator[k] += 
            inv_dx_sq*grad_dist*extension_fields[k][idx_neighbor1];
          extension_fields_denominator[k] += inv_dx_sq*grad_dist;
        }

      } else {

        double grad_dist = dist_updated - phi_upwind1[dir];

        for (k = 0; k < num_extension_fields; k++) {
          extension_fields_numerator[k] += 
            inv_dx_sq*grad_dist*extension_fields[k][idx_neighbor1];
          extension_fields_denominator[k] += inv_dx_sq*grad_dist;
        } 

      } /* end switch on second_order_switch[dir] */
 
    } /* end case: current direction used */
  } /* loop over coordinate directions */


  /* set updated quantities */
  LSM_FMM_IDX(idx_cur_gridpoint, grid_idx, grid_dims);
  distance_function[idx_cur_gridpoint] = dist_updated;
  for (k = 0; k < num_extension_fields; k++) {
    extension_fields[k][idx_cur_gridpoint] =
      extension_fields_numerator[k]/extension_fields_denominator[k];
  }

  /* free memory allocated for extension field calculations */
  free(extension_fields_numerator);
  free(extension_fields_denominator);

  return dist_updated;
}

#endif
