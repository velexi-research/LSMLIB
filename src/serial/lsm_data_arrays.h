/*
 * File:        lsm_data_arrays.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for LSM_DataArrays data structure and functions
 *              that support serial LSMLIB calculations
 */

#ifndef included_lsm_data_arrays_h
#define included_lsm_data_arrays_h

#include "LSMLIB_config.h"

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
#include "lsm_file.h"

/*!
 * Structure 'LSM_DataArrays' stores pointers for all arrays needed in a
 * typical LSM computation.
 */
typedef struct _LSM_DataArrays
{
  /* level set function at different time integration steps*/
  LSMLIB_REAL  *phi, *phi_stage1, *phi_stage2, *phi_next;
  
  /* extra storage in reinitialization */
  LSMLIB_REAL  *phi0;      

  /* extra storage for previous step functions */
  LSMLIB_REAL  *phi_prev;  
   
  /* extra storage for application specific purposes */
  LSMLIB_REAL  *phi_extra;  
  
  /* mask is the level set function that defines restricted domains */
  LSMLIB_REAL  *mask;

   /* LS Equation right hand side */
  LSMLIB_REAL  *lse_rhs; 
  
  /* 1st order derivatives (upwinding and central) */
  LSMLIB_REAL  *phi_x_plus, *phi_x_minus, *phi_x;
  LSMLIB_REAL  *phi_y_plus, *phi_y_minus, *phi_y;
  LSMLIB_REAL  *phi_z_plus, *phi_z_minus, *phi_z;
  
  /* scratch space for divided differences */
  LSMLIB_REAL  *D1, *D2, *D3;
  
  /* 2nd order derivatives */
  LSMLIB_REAL  *phi_xx, *phi_yy, *phi_xy, *phi_zz, *phi_xz, *phi_yz;
  
  /* normal velocity */
  LSMLIB_REAL *normal_velocity; 
  
   /* external velocity field */
  LSMLIB_REAL *external_velocity_x;
  LSMLIB_REAL *external_velocity_y;
  LSMLIB_REAL *external_velocity_z;
   
  /* arrays defining narrow band position */
  unsigned char *narrow_band;
  int    num_index_pts;
  int    *index_x, *index_y, *index_z;
  int    n_lo[10], n_hi[10]; //10 levels should be more than enough
  
  /* array for outer narrow band points storage */
  int  *index_outer_pts;
  int  num_alloc_index_outer_pts;
  int  nlo_outer_plus, nhi_outer_plus;
  int  nlo_outer_minus, nhi_outer_minus;

  /* special narrow band type of storage for solid voxels */
  unsigned char *solid_narrow_band;
  int    solid_num_index_pts;
  int    *solid_index_x, *solid_index_y, *solid_index_z;
  int    solid_n_lo[10], solid_n_hi[10];
  
  LSMLIB_REAL *solid_normal_x, *solid_normal_y, *solid_normal_z;

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
 * destroysMemoryForLSMDataArrays() frees ALL memory allocated for the data 
 * arrays contained within the LSM_DataArrays structure as well as the structure
 * itself.
 *   
 * Arguments:
 *  - lsm_data_arrays(in):  pointer to LSM_DataArrays 
 *   
 * Return value:            none
 *   
 */
void destroyLSMDataArrays(LSM_DataArrays *lsm_data_arrays);

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
 *  - zip_status(in):  integer indicating compression of the file 
 *                     (NO_ZIP,GZIP,BZIP2) 
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
void writeDataArray(LSMLIB_REAL *data, Grid *grid, char *file_name,
                    int zip_status);

/*!
 * writeDataArrayNoGrid() writes the specified data array out to a binary file.
 *
 * The data is output in the following order:
 * -# grid/data dimensions 
 * -# values of data array at all grid points.
 *
 * Arguments:
 *  - data (in):       data array to be output to file
 *  - n (in):          array of 3 integers indicating grid/array size in each
 *                     direction 
 *  - file_name (in):  name of output file
 *  - zip_status(in):  integer indicating compression of the file 
 *                     (NO_ZIP,GZIP,BZIP2) 
 *   
 * Return value:       none
 *   
 * NOTES: 
 * - writeDataArrayNoGrid() is used for 2d and 3d data arrays.  For 2d
 *   data, the third dimension MUST be set to 1 
 *
 * - If a file with the specified file_name already exists, it is
 *   overwritten.
 *
 */   
void writeDataArrayNoGrid(LSMLIB_REAL *data, int *n, char *file_name,
                    int zip_status);
		    

/*!
 * readDataArray() loads the data from a binary file into a LSMLIB_REAL
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
 * - Function recognizes if the file name contains .gz or .bz2 extention
 *   and uncompresses the file accordingly.
 */   
LSMLIB_REAL *readDataArray(int *grid_dims, char *file_name);


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
 *  - zip_status(in):     integer indicating compression of the file 
 *                        (NO_ZIP,GZIP,BZIP2)
 *   
 * Return value:         none
 *
 * NOTES: 
 * - If a file with the specified file_name already exists, it is
 *   overwritten.
 *
 */   
void writeDataArray1d(LSMLIB_REAL *data, int num_gridpts, char *file_name,
                      int  zip_status);


/*!
 * readDataArray1d() loads the data from a binary file into a LSMLIB_REAL
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
 * - Function recognizes if the file name contains .gz or .bz2 extention
 *   and uncompresses the file accordingly.
 */   
LSMLIB_REAL *readDataArray1d(int *num_elements, char *file_name);


#ifdef __cplusplus
}
#endif

#endif
