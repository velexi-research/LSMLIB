/*
 * File:        lsm_data_arrays.h
 * Copyright:   (c) 2005-2006 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.5 $
 * Modified:    $Date: 2006/11/03 23:15:52 $
 * Description: Header file for LSM_DataArrays data structure and functions
 *              that support serial LSMLIB calculations
 */

#ifndef included_lsm_data_arrays_h
#define included_lsm_data_arrays_h

#ifdef __cplusplus
extern "C" {
#endif


/*! \file lsm_data_arrays.h
 *
 * \brief
 * @ref lsm_data_arrays.h provides support for managing the data arrays
 * required by typical level set method calculcations.
 *
 */

#include "lsm_grid.h"

/*!
 * Structure 'LSM_DataArrays' stores pointers for all arrays needed in a
 * typical LSM computation.
 */
typedef struct _LSM_DataArrays
{
  /* level set function at different time integration steps*/
  double  *phi, *phi_stage1, *phi_stage2, *phi_next;
  
  /* extra storage in reinitialization */
  double  *phi0;      

  /* extra storage for previous step functions */
  double  *phi_prev;  
  
  /* mask is the level set function that defines restricted domains */
  double  *mask;

   /* LS Equation right hand side */
  double  *lse_rhs; 
  
  /* 1st order derivatives (upwinding and central) */
  double  *phi_x_plus, *phi_x_minus, *phi_x;
  double  *phi_y_plus, *phi_y_minus, *phi_y;
  double  *phi_z_plus, *phi_z_minus, *phi_z;
  
  /* scratch space for divided differences */
  double  *D1, *D2, *D3;
  
  /* 2nd order derivatives */
  double  *phi_xx, *phi_yy, *phi_xy, *phi_zz, *phi_xz, *phi_yz;
  
  /* normal velocity */
  double *normal_velocity; 
  
   /* external velocity field */
  double *external_velocity_x;
  double *external_velocity_y;
  double *external_velocity_z;
   
  /* arrays defining narrow band position */
  unsigned char *narrow_band;
  int    num_index_pts;
  int    *index_x, *index_y, *index_z;
  int    n_lo[6], n_hi[6]; //more?
  
  /* array for outer narrow band points storage */
  int  *index_outer_pts;
  int  num_alloc_index_outer_pts;
  int  nlo_outer_plus, nhi_outer_plus;
  int  nlo_outer_minus, nhi_outer_minus;

}  LSM_DataArrays;


/*!
 * allocateLSMDataArrays() allocates a LSM_DataArrays data structure 
 * and initializes all of its data pointers to a non-NULL dummy pointer.
 *  
 * Arguments:    none
 *
 * Return value: pointer to LSM_DataArrays structure 
 *
 */
LSM_DataArrays *allocateLSMDataArrays(void);




/*!
 * allocateMemoryForLSMDataArrays() allocates memory for the data 
 * arrays contained within the LSM_DataArrays structure.
 *    
 * Arguments:
 *  - lsm_arrays(in):  pointer to LSM_DataArrays structure
 *  - grid(in):        pointer to Grid 
 *    
 * Return value:       none
 *
 * NOTES: 
 * - The memory for the LSM_DataArrays structure MUST be allocated
 *   (e.g. using the allocateLSMDataArrays() function) before 
 *   allocateMemoryForLSMDataArrays() is called.
 *
 * - Memory is allocated only for array pointers in lsm_data_arrays that 
 *   are not already associated with allocated or that are set to NULL.
 *   If memory has already been allocated for a particular data array or
 *  the data pointer is set to NULL, it will not be reallocated.
 *
 */
void allocateMemoryForLSMDataArrays(
  LSM_DataArrays *lsm_data_arrays,
  Grid *grid);


/*!
 * freeMemoryForLSMDataArrays() frees ALL memory allocated for the data 
 * arrays contained within the LSM_DataArrays structure.
 *   
 * Arguments:
 *  - lsm_data_arrays(in):  pointer to LSM_DataArrays 
 *   
 * Return value:            none
 *   
 */
void freeMemoryForLSMDataArrays(LSM_DataArrays *lsm_arrays);


/*!
 * writeDataArray() writes the specified data array out to a binary file.
 *
 * The data is output in the following order:
 * -# grid dimensions 
 * -# values of data array at all grid points.
 *
 * Arguments:
 *  - data (in):       data array to be output to file
 *  - grid (in):       pointer to Grid 
 *  - file_name (in):  name of output file 
 *   
 * Return value:       none
 *   
 * NOTES: 
 * - writeDataArray() is used for 2d and 3d data arrays.  For 2d
 *   data, the third grid dimension MUST be set to 1 (which is 
 *   the default behavior of the createGrid() function).
 *
 * - If a file with the specified file_name already exists, it is
 *   overwritten.
 *
 */   
void writeDataArray(double *data, Grid *grid, char *file_name);


/*!
 * readDataArray() loads the data from a binary file into a double
 * array and returns it to the user.  
 *   
 * Arguments:
 *  - grid_dims (out):  dimensions of grid (read from file)
 *  - file_name (in):   name of input file 
 *   
 * Return value:        pointer to data array loaded from file
 *   
 * NOTES: 
 * - readDataArray() dynamically allocates memory for the data array 
 *   that is returned.
 *
 * - The memory for grid_dims is assumed to be allocated by the user.
 *
 * - readDataArray() is used for 2d and 3d data arrays.  For 2d
 *   data, the third grid dimension is set to 1.
 *
 */   
double *readDataArray(int *grid_dims, char *file_name);


/*!
 * writeDataArray1d() writes the specified data array out to a binary file.
 *
 * The data is output in the following order:
 * -# number of grid points
 * -# values of data array at all grid points.
 *
 * Arguments:
 *  - data (in):          data array to be output to file
 *  - num_elements (in):  number of elements in the array
 *  - file_name (in):     name of output file 
 *   
 * Return value:         none
 *
 * NOTES: 
 * - If a file with the specified file_name already exists, it is
 *   overwritten.
 *
 */   
void writeDataArray1d(double *data, int num_gridpts, char *file_name);


/*!
 * readDataArray1d() loads the data from a binary file into a double
 * array and returns it to the user.  
 *   
 * Arguments:
 *  - num_elements (out):  number of elements in the array (read from file)
 *  - file_name (in):      name of input file 
 *   
 * Return value:           none
 *
 * NOTES: 
 * - readDataArray1d() dynamically allocates memory for the data array 
 *   that is returned.
 *
 * - The memory for num_elements is assumed to be allocated by the user.
 *
 *
 */   
double *readDataArray1d(int *num_elements, char *file_name);


#ifdef __cplusplus
}
#endif

#endif
