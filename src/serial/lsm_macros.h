/*
 * File:        lsm_macros.h
 * Copyright:   (c) 2005-2006 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date: 2006/11/03 23:22:14 $
 * Description: Header file with helpful macros for manipulating data arrays
 */

#ifndef included_lsm_macros_h
#define included_lsm_macros_h

#ifdef __cplusplus
extern "C" {
#endif


/*! \file lsm_macros.h
 *
 * \brief
 * @ref lsm_macros.h provides some helpful macros for manipulating data 
 * arrays.
 *
 */

#include <float.h>
#include "lsm_grid.h"


/*!
 * SET_DATA_TO_CONSTANT() sets entire data array to a constant value. 
 *
 * Arguments:
 *  - data (out):  data array of size (grid->num_gridpts) 
 *  - grid (in):   pointer to Grid 
 *  - value (in):  constant value to set data to
 *
 */
#define SET_DATA_TO_CONSTANT(data, grid, value)                            \
{                                                                          \
  int idx;                                                                 \
  for (idx = 0; idx < grid->num_gridpts; idx++)                            \
  {                                                                        \
    data[idx] = value;                                                     \
  }                                                                        \
}

 
/*!
 * NEGATE_DATA() sets each data value to its negative.
 *
 * Arguments:
 *  - data (in/out):  data array of size (grid->num_gridpts) 
 *  - grid (in):      pointer to Grid 
 *
 */
#define NEGATE_DATA(data, grid)                                            \
{                                                                          \
  int idx;                                                                 \
  for (idx = 0; idx < grid->num_gridpts; idx++)                            \
  {                                                                        \
     data[idx] = -data[idx];                                               \
  }                                                                        \
}   


/*!
 * IMPOSE_MASK() imposes constraint of the type 'phi <= mask'.  That is,
 * each point of array phi_masked is set to a maximum of the values
 * phi and mask. This can be used for motion in restricted (but fixed) 
 * domains (i.e. the physical domain is determined by 'mask').
 *
 * Arguments:
 *  - phi_masked (out):  array with constraint 'phi<=mask' imposed
 *  - phi (in):          array containing original values of phi
 *  - mask (in):         array containing mask values
 *  - grid (in):         pointer to Grid 
 *
 */
#define IMPOSE_MASK(phi_masked, mask, phi, grid)                           \
{                                                                          \
  int idx;                                                                 \
  for(idx = 0; idx < grid->num_gridpts; idx++)                             \
  {                                                                        \
    phi_masked[idx] = (mask[idx] > phi[idx]) ? mask[idx] : phi[idx];       \
  }                                                                        \
}

/*!
 * IMPOSE_MASK_LOCAL() imposes constraint of the type 'phi <= mask'
 * only at local points whose coordinates are stored in index_x,
 * index_y, index_z arrays.  That is,
 * each point of array phi_masked is set to a maximum of the values
 * phi and mask. This can be used for motion in restricted (but fixed) 
 * domains (i.e. the physical domain is determined by 'mask').
 *
 * Arguments:
 *  - phi_masked (out):  array with constraint 'phi<=mask' imposed
 *  - phi (in):          array containing original values of phi
 *  - mask (in):         array containing mask values
 *  - grid (in):         pointer to Grid
 *
 */
#define IMPOSE_MASK_LOCAL(phi_masked, mask, phi, grid, p)                    \
{                                                                            \
   int idx, l;                                                               \
                                                                             \
   if(grid->num_dims == 3)                                                   \
     for(l = (p->n_lo)[0]; l < (p->n_hi)[0]; l++)                            \
     {                                                                       \
	idx = (p->index_x)[l] + (p->index_y)[l]*nx + (p->index_z)[l]*nxy;    \
        phi_masked[idx] = (mask[idx] > phi[idx]) ? mask[idx] : phi[idx];     \
     }                                                                       \
   else                                                                      \
     for(l = (p->n_lo)[0]; l < (p->n_hi)[0]; l++)                            \
     {                                                                       \
	idx = (p->index_x)[l] + (p->index_y)[l]*nx;                          \
        phi_masked[idx] = (mask[idx] > phi[idx]) ? mask[idx] : phi[idx];     \
     }                                                                       \
}


/*!
 * COPY_DATA() copies values from data_src into data_dst.
 *
 * Arguments:
 *  - data_dst (out):  destination data array
 *  - data_src (in):   source data array 
 *  - grid (in):       pointer to Grid 
 *
 */
#define COPY_DATA(data_dst, data_src, grid)                                \
{                                                                          \
  int idx;                                                                 \
  for(idx = 0; idx < grid->num_gridpts; idx++)                             \
  {                                                                        \
    data_dst[idx] = data_src[idx];                                         \
  }                                                                        \
}


/*!
 * COMPUTE_MAX_ABS_ERR() computes maximal absolute value over all points
 * for (data1 - data2).
 *
 * Arguments:
 *  - max_abs_err (out):  maximal absolute value of (data1 - data2)
 *  - data1 (in):         array containing data1
 *  - data2 (in):         array containing data2
 *  - grid (in):          pointer to Grid 
 *
 */
#define COMPUTE_MAX_ABS_ERR(max_abs_err, data1, data2, grid)               \
{                                                                          \
  int idx;                                                                 \
  double min_err, max_err, err1, err2, err;                                \
  min_err = DBL_MAX; max_err = 0.0;                                        \
  for(idx = 0; idx < grid->num_gridpts; idx++)                             \
  {                                                                        \
    err = data1[idx] - data2[idx];                                         \
      if(err < min_err) min_err = err;                                     \
      if(err > max_err) max_err = err;                                     \
   }                                                                       \
   err1 = fabs(min_err); err2 = fabs(max_err);                             \
   max_abs_err = ( err1 > err2 ) ? err1 : err2;                            \
}   
 

/*!
 * COMPUTE_MAX_ABS_DATA() computes a maximal absolute value of the
 * data values in the array.
 *
 * Arguments:
 *  - max_abs (out):  maximal absolute value of data
 *  - data (in):      array containing data
 *  - grid (in):      pointer to Grid 
 *
 */
#define COMPUTE_MAX_ABS_DATA(max_abs, data, grid)                          \
{                                                                          \
  int idx;                                                                 \
  double min_err, max_err, err1, err2, err;                                \
  min_err = DBL_MAX; max_err = 0.0;                                        \
  for(idx = 0; idx < grid->num_gridpts; idx++)                             \
  {                                                                        \
    err = data[idx];                                                       \
    if(err < min_err) min_err = err;                                       \
    if(err > max_err) max_err = err;                                       \
  }                                                                        \
  err1 = fabs(min_err); err2 = fabs(max_err);                              \
  max_abs = ( err1 > err2 ) ? err1 : err2;                                 \
}


/*!
 * PRINT_ARRAY_2D() prints a 2d array as a matrix of values. Can be used 
 * for debugging purposes.
 * 
 * Arguments:
 *  - data (in):  data array
 *  - grid (in):  pointer to Grid 
 *
 */
#define PRINT_ARRAY_2D(data, grid)                                    \
{                                                                          \
  int i,j;                                                                 \
  for(j = grid->jhi_gb; j >= grid->jlo_gb; j--)                              \
  {                                                                        \
    for(i = grid->ilo_gb; i <= grid->ihi_gb; i++)                            \
    {                                                                      \
      idx = i + j*(grid->grid_dims)[0];                                    \
      printf("%g ", data[idx]);                                            \
    }                                                                      \
    printf("\n");                                                          \
  }                                                                        \
}


#ifdef __cplusplus
}
#endif

#endif
