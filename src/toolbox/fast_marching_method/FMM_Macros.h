/*
 * File:        FMM_Macros.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date: 2006/08/13 13:30:52 $
 * Description: Header file that defines several common macros used 
 *              for Fast Marching Method calculations
 */

/*! \file FMM_Macros.h
 *
 * \brief
 * @ref FMM_Macros.h provides several macros used for Fast Marching
 *      Method calculations.
 *
 */

#ifndef included_FMM_Macros_h
#define included_FMM_Macros_h
 

/*============================= Constants ===========================*/
#define LSM_FMM_TRUE                   (1)
#define LSM_FMM_FALSE                  (0)
#define LSM_FMM_DEFAULT_UPDATE_VALUE   (-1)
#define LSM_FMM_ZERO_TOL               (1e-8)


/*========================== Error Codes ============================*/
#define LSM_FMM_ERR_SUCCESS                                 (0)
#define LSM_FMM_ERR_FMM_DATA_CREATION_ERROR                 (1)
#define LSM_FMM_ERR_INVALID_SPATIAL_DISCRETIZATION_ORDER    (2)


/*======================= Helper Functions ==========================*/
#define LSM_FMM_ABS(x)            ((x) > 0 ? (x) : -1.0*(x))

/*
 * LSM_FMM_IDX() computes the array index for the specified 
 * grid index and grid dimensions.
 *
 * Arguments:
 *   idx (out):       array index
 *   grid_idx (in):   grid index
 *   grid_dims (in):  grid dimensions
 * 
 * NOTES:
 *  (1) idx MUST be a valid l-value.
 *  (2) FMM_NDIM MUST be defined by user code.
 *
 */
#define LSM_FMM_IDX(idx, grid_idx, grid_dims)                            \
{                                                                        \
  int dir;                                                               \
  int grid_size_lower_dims = 1;                                          \
  idx = 0;                                                               \
  for (dir = 0; dir < FMM_NDIM; dir++) {                                 \
    idx += grid_idx[dir]*grid_size_lower_dims;                           \
    grid_size_lower_dims *= grid_dims[dir];                              \
  }                                                                      \
}

/*
 * LSM_FMM_IDX_OUT_OF_BOUNDS() determines whether the given 
 * grid index lies in the computational domain.
 *
 * Arguments:
 *   result (out):    1 if grid index is out of bounds; 0 otherwise.
 *   grid_idx (in):   grid index
 *   grid_dims (in):  grid dimensions
 * 
 * NOTES:
 *  (1) result MUST be a valid l-value.
 *  (2) FMM_NDIM MUST be defined by user code.
 *
 */
#define LSM_FMM_IDX_OUT_OF_BOUNDS(result, grid_idx, grid_dims)           \
{                                                                        \
  int dir;                                                               \
  result = 0;                                                            \
  for (dir = 0; dir < FMM_NDIM; dir++) {                                 \
    if ( (grid_idx[dir]<0) || (grid_idx[dir]>grid_dims[dir]-1) ) {       \
      result = 1;                                                        \
      break;                                                             \
    }                                                                    \
  }                                                                      \
}

#endif
