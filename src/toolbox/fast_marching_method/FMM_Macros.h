/*
 * File:        FMM_Macros.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
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


/*========================== Error Codes ============================*/
#define LSM_FMM_ERR_SUCCESS                                 (0)
#define LSM_FMM_ERR_FMM_DATA_CREATION_ERROR                 (1)
#define LSM_FMM_ERR_INVALID_SPATIAL_DISCRETIZATION_ORDER    (2)


/*======================= Helper Functions ==========================*/
/* LSM_FMM_ABS() computes the absolute value of its argument.        */
/* This macro is defined so that LSMLIB does not have to rely on the */
/* C math library.                                                   */
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
  int lsm_fmm_dir;                                                       \
  int grid_size_lower_dims = 1;                                          \
  idx = 0;                                                               \
  for (lsm_fmm_dir = 0; lsm_fmm_dir < FMM_NDIM; lsm_fmm_dir++) {         \
    idx += grid_idx[lsm_fmm_dir]*grid_size_lower_dims;                   \
    grid_size_lower_dims *= grid_dims[lsm_fmm_dir];                      \
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
  int lsm_fmm_dir;                                                       \
  result = 0;                                                            \
  for (lsm_fmm_dir = 0; lsm_fmm_dir < FMM_NDIM; lsm_fmm_dir++) {         \
    if (  (grid_idx[lsm_fmm_dir]<0)                                      \
       || (grid_idx[lsm_fmm_dir]>grid_dims[lsm_fmm_dir]-1) ) {           \
      result = 1;                                                        \
      break;                                                             \
    }                                                                    \
  }                                                                      \
}

#endif
