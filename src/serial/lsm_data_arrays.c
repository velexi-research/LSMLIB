/*
 * File:        lsm_data_arrays.c
 * Copyright:   (c) 2005-2006 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.3 $
 * Modified:    $Date: 2006/05/25 19:53:02 $
 * Description: Implementation file for LSM_DataArrays structure
 */

/*****************************************************************
 * TO DO:
 *   Having the following would be nice.
 *
 *   - have endian information and bite-swap id needed
 *
 *****************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "lsm_data_arrays.h"

#define DSZ  sizeof(LSMLIB_REAL)
#define ISZ  sizeof(int)
#define UCSZ sizeof(unsigned char)

#define LSMLIB_SERIAL_dummy_pointer        ((LSMLIB_REAL*)(-1))
#define LSMLIB_SERIAL_dummy_pointer_int    ((int*)(-1))
#define LSMLIB_SERIAL_dummy_pointer_uchar  ((unsigned char*)(-1))

LSM_DataArrays *allocateLSMDataArrays(void)
{
  LSM_DataArrays *lsm_data_arrays;
  int i;
  
  lsm_data_arrays  = (LSM_DataArrays *)malloc(sizeof(LSM_DataArrays));
  
  lsm_data_arrays->phi = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_stage1 = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_stage2 = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_next = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi0 = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_prev = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_extra = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->mask = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->lse_rhs = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_x_plus = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_x_minus = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_x = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_y_plus = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_y_minus = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_y  = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_z_plus = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_z_minus = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_z = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_xx = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_yy = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_zz = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_xy = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_xz = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->phi_yz = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->normal_velocity = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->external_velocity_x = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->external_velocity_y = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->external_velocity_z = LSMLIB_SERIAL_dummy_pointer;
  
  lsm_data_arrays->narrow_band = LSMLIB_SERIAL_dummy_pointer_uchar;
  lsm_data_arrays->index_x = LSMLIB_SERIAL_dummy_pointer_int;
  lsm_data_arrays->index_y = LSMLIB_SERIAL_dummy_pointer_int;
  lsm_data_arrays->index_z = LSMLIB_SERIAL_dummy_pointer_int;
  lsm_data_arrays->num_index_pts = 0;
  for(i=0; i < 10; i++)
  {
    lsm_data_arrays->n_lo[i] = lsm_data_arrays->n_hi[i] = 0;
  }
  
  lsm_data_arrays->index_outer_pts = LSMLIB_SERIAL_dummy_pointer_int;
  lsm_data_arrays->num_alloc_index_outer_pts = 0;
  lsm_data_arrays->nlo_outer_plus = lsm_data_arrays->nhi_outer_plus = 0;
  lsm_data_arrays->nlo_outer_minus = lsm_data_arrays->nhi_outer_minus = 0;
  
  lsm_data_arrays->solid_narrow_band = LSMLIB_SERIAL_dummy_pointer_uchar;
  lsm_data_arrays->solid_index_x = LSMLIB_SERIAL_dummy_pointer_int;
  lsm_data_arrays->solid_index_y = LSMLIB_SERIAL_dummy_pointer_int;
  lsm_data_arrays->solid_index_z = LSMLIB_SERIAL_dummy_pointer_int;
  lsm_data_arrays->solid_num_index_pts = 0;
  for(i=0; i < 10; i++)
  {
    lsm_data_arrays->solid_n_lo[i] = lsm_data_arrays->solid_n_hi[i] = 0;
  }
  
  lsm_data_arrays->solid_normal_x = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->solid_normal_y = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->solid_normal_z = LSMLIB_SERIAL_dummy_pointer;
  
  lsm_data_arrays->D1 = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->D2 = LSMLIB_SERIAL_dummy_pointer;
  lsm_data_arrays->D3 = LSMLIB_SERIAL_dummy_pointer;
  
  return  lsm_data_arrays;
}


void destroyLSMDataArrays(LSM_DataArrays *lsm_data_arrays)
{
  if (lsm_data_arrays) {
    freeMemoryForLSMDataArrays(lsm_data_arrays);
    free(lsm_data_arrays);
  }
}


void  allocateMemoryForLSMDataArrays(
  LSM_DataArrays *lsm_data_arrays,
  Grid *grid)
{
 /* Only arrays that are equal to LSMLIB_SERIAL_dummy_pointer will get memory allocated.
  *   If memory allocation is to be avoided, set the pointer to NULL,
  *   Non-NULL pointers different from LSMLIB_SERIAL_dummy_pointer are assumed allocated
  *   elsewhere and that will not be overridden.
  */
        
  if( lsm_data_arrays->phi == LSMLIB_SERIAL_dummy_pointer )
    lsm_data_arrays->phi = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
    
  if( lsm_data_arrays->phi_stage1  == LSMLIB_SERIAL_dummy_pointer )
    lsm_data_arrays->phi_stage1 = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
  
  if( lsm_data_arrays->phi_stage2  == LSMLIB_SERIAL_dummy_pointer )  
    lsm_data_arrays->phi_stage2 = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
  
  if( lsm_data_arrays->phi_next  == LSMLIB_SERIAL_dummy_pointer ) 
    lsm_data_arrays->phi_next = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
  
  if( lsm_data_arrays->phi0  == LSMLIB_SERIAL_dummy_pointer )
    lsm_data_arrays->phi0 = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);  
  
  if( lsm_data_arrays->phi_prev  == LSMLIB_SERIAL_dummy_pointer ) 
    lsm_data_arrays->phi_prev = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
  
  if( lsm_data_arrays->phi_extra  == LSMLIB_SERIAL_dummy_pointer ) 
    lsm_data_arrays->phi_extra = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
   
  if( lsm_data_arrays->mask  == LSMLIB_SERIAL_dummy_pointer )
    lsm_data_arrays->mask = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);  
    
  if( lsm_data_arrays->lse_rhs  == LSMLIB_SERIAL_dummy_pointer )
    lsm_data_arrays->lse_rhs = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
        
  if( lsm_data_arrays->phi_x_plus  == LSMLIB_SERIAL_dummy_pointer ) 
    lsm_data_arrays->phi_x_plus  = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
   
  if( lsm_data_arrays->phi_x_minus  == LSMLIB_SERIAL_dummy_pointer )     
    lsm_data_arrays->phi_x_minus = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
    
  if( lsm_data_arrays->phi_x  == LSMLIB_SERIAL_dummy_pointer )     
    lsm_data_arrays->phi_x = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
     
  if( lsm_data_arrays->phi_y_plus  == LSMLIB_SERIAL_dummy_pointer )      
    lsm_data_arrays->phi_y_plus  = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
  
  if( lsm_data_arrays->phi_y_minus  == LSMLIB_SERIAL_dummy_pointer )     
    lsm_data_arrays->phi_y_minus = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
  
  if( lsm_data_arrays->phi_y  == LSMLIB_SERIAL_dummy_pointer )     
    lsm_data_arrays->phi_y = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
  
  if(grid->num_dims == 3)
  {
      if( lsm_data_arrays->phi_z_plus  == LSMLIB_SERIAL_dummy_pointer )
        lsm_data_arrays->phi_z_plus  = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
      
      if( lsm_data_arrays->phi_z_minus  == LSMLIB_SERIAL_dummy_pointer )	    
        lsm_data_arrays->phi_z_minus = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
     
      if( lsm_data_arrays->phi_z  == LSMLIB_SERIAL_dummy_pointer )	    
        lsm_data_arrays->phi_z = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
  }
  else
  {
     lsm_data_arrays->phi_z_plus = (LSMLIB_REAL *)NULL;
     lsm_data_arrays->phi_z_minus = (LSMLIB_REAL *)NULL;
     lsm_data_arrays->phi_z = (LSMLIB_REAL *)NULL;
  }
  
  if( lsm_data_arrays->phi_xx  == LSMLIB_SERIAL_dummy_pointer )  
    lsm_data_arrays->phi_xx = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
    
  if( lsm_data_arrays->phi_xy  == LSMLIB_SERIAL_dummy_pointer )  
    lsm_data_arrays->phi_xy = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
    
  if( lsm_data_arrays->phi_yy  == LSMLIB_SERIAL_dummy_pointer )  
    lsm_data_arrays->phi_yy = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
  
  if(grid->num_dims == 3)
  {
     if( lsm_data_arrays->phi_zz  == LSMLIB_SERIAL_dummy_pointer )
       lsm_data_arrays->phi_zz = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
     
     if( lsm_data_arrays->phi_xz  == LSMLIB_SERIAL_dummy_pointer )
       lsm_data_arrays->phi_xz = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
     
     if( lsm_data_arrays->phi_yz  == LSMLIB_SERIAL_dummy_pointer )
       lsm_data_arrays->phi_yz = (LSMLIB_REAL*) calloc(grid->num_gridpts,DSZ);
  }
  else
  {
     lsm_data_arrays->phi_zz = (LSMLIB_REAL *)NULL;
     lsm_data_arrays->phi_xz = (LSMLIB_REAL *)NULL;
     lsm_data_arrays->phi_yz = (LSMLIB_REAL *)NULL;
  }
  
  if( lsm_data_arrays->normal_velocity  == LSMLIB_SERIAL_dummy_pointer )  
    lsm_data_arrays->normal_velocity = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
  
  if( lsm_data_arrays->external_velocity_x  == LSMLIB_SERIAL_dummy_pointer )  
    lsm_data_arrays->external_velocity_x = 
      (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
   
  if( lsm_data_arrays->external_velocity_y  == LSMLIB_SERIAL_dummy_pointer )  
    lsm_data_arrays->external_velocity_y = 
      (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
    
  if(grid->num_dims == 3)
  {
    if( lsm_data_arrays->external_velocity_z  == LSMLIB_SERIAL_dummy_pointer )  
      lsm_data_arrays->external_velocity_z = 
        (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);	
  }
  else lsm_data_arrays->external_velocity_z = (LSMLIB_REAL *)NULL;
  
  if (lsm_data_arrays->D1 == LSMLIB_SERIAL_dummy_pointer )
    lsm_data_arrays->D1 = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
  
  if (lsm_data_arrays->D2 == LSMLIB_SERIAL_dummy_pointer )  
    lsm_data_arrays->D2 = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
   
  if (lsm_data_arrays->D3 == LSMLIB_SERIAL_dummy_pointer )
    lsm_data_arrays->D3 = (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
        
  if( lsm_data_arrays->narrow_band == LSMLIB_SERIAL_dummy_pointer_uchar )  
    lsm_data_arrays->narrow_band = (unsigned char*) malloc(grid->num_gridpts*UCSZ);
    
  if( lsm_data_arrays->index_x == LSMLIB_SERIAL_dummy_pointer_int )
    lsm_data_arrays->index_x = (int*) malloc(grid->num_gridpts*ISZ); 
  
  if( lsm_data_arrays->index_y == LSMLIB_SERIAL_dummy_pointer_int )
    lsm_data_arrays->index_y = (int*) malloc(grid->num_gridpts*ISZ);
  
  if(grid->num_dims == 3)
  { 
    if( lsm_data_arrays->index_z == LSMLIB_SERIAL_dummy_pointer_int )
      lsm_data_arrays->index_z = (int*) malloc(grid->num_gridpts*ISZ);
  }
  else  lsm_data_arrays->index_z = (int*) NULL;
  
  if( lsm_data_arrays->index_outer_pts == LSMLIB_SERIAL_dummy_pointer_int )
  {
    lsm_data_arrays->index_outer_pts = (int*) malloc(grid->num_gridpts*ISZ);
    lsm_data_arrays->num_alloc_index_outer_pts = grid->num_gridpts;
  }  
    
    
  if( lsm_data_arrays->solid_narrow_band == LSMLIB_SERIAL_dummy_pointer_uchar )  
    lsm_data_arrays->solid_narrow_band = (unsigned char*) malloc(grid->num_gridpts*UCSZ);
    
  if( lsm_data_arrays->solid_index_x == LSMLIB_SERIAL_dummy_pointer_int )
    lsm_data_arrays->solid_index_x = (int*) malloc(grid->num_gridpts*ISZ); 
  
  if( lsm_data_arrays->solid_index_y == LSMLIB_SERIAL_dummy_pointer_int )
    lsm_data_arrays->solid_index_y = (int*) malloc(grid->num_gridpts*ISZ);
  
  if(grid->num_dims == 3)
  { 
    if( lsm_data_arrays->solid_index_z == LSMLIB_SERIAL_dummy_pointer_int )
      lsm_data_arrays->solid_index_z = (int*) malloc(grid->num_gridpts*ISZ);
  }
  else  lsm_data_arrays->solid_index_z = (int*) NULL;
  
  
  
  if( lsm_data_arrays->solid_normal_x  == LSMLIB_SERIAL_dummy_pointer )  
    lsm_data_arrays->solid_normal_x = 
      (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
   
  if( lsm_data_arrays->solid_normal_y  == LSMLIB_SERIAL_dummy_pointer )  
    lsm_data_arrays->solid_normal_y = 
      (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);
    
  if(grid->num_dims == 3)
  {
    if( lsm_data_arrays->solid_normal_z  == LSMLIB_SERIAL_dummy_pointer )  
      lsm_data_arrays->solid_normal_z = 
        (LSMLIB_REAL*) malloc(grid->num_gridpts*DSZ);	
  }
  else lsm_data_arrays->solid_normal_z = (LSMLIB_REAL *)NULL;  
}     



void  freeMemoryForLSMDataArrays(LSM_DataArrays *lsm_data_arrays)
{   
  free(lsm_data_arrays->phi);
  
  free(lsm_data_arrays->phi_stage1); 
  free(lsm_data_arrays->phi_stage2);
  free(lsm_data_arrays->phi_next);
  
  free(lsm_data_arrays->phi0);      
  free(lsm_data_arrays->phi_prev);
  free(lsm_data_arrays->phi_extra);
  
  free(lsm_data_arrays->mask); 
  
  free(lsm_data_arrays->lse_rhs);
   
  free(lsm_data_arrays->phi_x_plus);
  free(lsm_data_arrays->phi_x_minus);    
  free(lsm_data_arrays->phi_x);    
  free(lsm_data_arrays->phi_y_plus);
  free(lsm_data_arrays->phi_y_minus);
  free(lsm_data_arrays->phi_y);    
  free(lsm_data_arrays->phi_z_plus);
  free(lsm_data_arrays->phi_z_minus);  
  free(lsm_data_arrays->phi_z);    
  
  free(lsm_data_arrays->phi_xx);
  free(lsm_data_arrays->phi_xy);
  free(lsm_data_arrays->phi_yy);
  free(lsm_data_arrays->phi_xz);
  free(lsm_data_arrays->phi_yz);
  free(lsm_data_arrays->phi_zz);
  
  free(lsm_data_arrays->normal_velocity);
  free(lsm_data_arrays->external_velocity_x);
  free(lsm_data_arrays->external_velocity_y);
  free(lsm_data_arrays->external_velocity_z);

  free(lsm_data_arrays->narrow_band);
  free(lsm_data_arrays->index_x);
  free(lsm_data_arrays->index_y);
  free(lsm_data_arrays->index_z);

  free(lsm_data_arrays->index_outer_pts);
  
  free(lsm_data_arrays->solid_narrow_band);
  free(lsm_data_arrays->solid_index_x);
  free(lsm_data_arrays->solid_index_y);
  free(lsm_data_arrays->solid_index_z);
  
  free(lsm_data_arrays->solid_normal_x);
  free(lsm_data_arrays->solid_normal_y);
  free(lsm_data_arrays->solid_normal_z);
  
  free(lsm_data_arrays->D1);
  free(lsm_data_arrays->D2);
  free(lsm_data_arrays->D3);  
}
   

void writeDataArray(LSMLIB_REAL *data, Grid *grid, char *file_name,int zip_status)
{
   FILE *fp;
   
   fp = fopen(file_name,"w");   

   /* write grid dimensions */
   fwrite(grid->grid_dims_ghostbox, sizeof(int), 3, fp); 

   /* write data array */
   fwrite(data, DSZ, grid->num_gridpts, fp);

   fclose(fp);
   zipFile(file_name,zip_status);
}


LSMLIB_REAL *readDataArray(int *grid_dims_ghostbox,char *file_name)
{
   FILE    *fp;
   int     zip_status;
   int     num_gridpts;
   LSMLIB_REAL    *data = NULL;
   char    *file_base;
   
   checkUnzipFile(file_name,&zip_status,&file_base);
   
   fp = fopen(file_base,"r");

   if( fp != NULL)
   {
     /* read grid dimensions */
     fread(grid_dims_ghostbox, sizeof(int), 3, fp); 
  
     /* allocate memory for data array */ 
     num_gridpts = grid_dims_ghostbox[0] * grid_dims_ghostbox[1]
               * grid_dims_ghostbox[2];
     data = (LSMLIB_REAL *) malloc(num_gridpts*DSZ);

     /* read data array */ 
     fread(data, DSZ, num_gridpts, fp);

     fclose(fp);
     
     zipFile(file_base,zip_status);
   }
   else
   {
      printf("\nCould not open file %s",file_name);
   }
   free(file_base);
   return data;
}


void writeDataArray1d(LSMLIB_REAL *data, int num_elements, char *file_name,
                      int zip_status)
{
   FILE *fp;
   
   fp = fopen(file_name,"w");

   /* write number of elements */
   fwrite(&num_elements, sizeof(int), 1, fp); 

   /* write data array */
   fwrite(data, DSZ, num_elements, fp);

   fclose(fp);
   zipFile(file_name,zip_status);
}


LSMLIB_REAL *readDataArray1d(int *num_elements, char *file_name)
{
   FILE    *fp;
   LSMLIB_REAL  *data;
   char    *file_base;
   int     zip_status;
   
   checkUnzipFile(file_name,&zip_status,&file_base);
   
   fp = fopen(file_base,"r");

   if(fp)
   {
     /* read number of elements */
     fread(num_elements, sizeof(int), 1, fp); 
   
     /* allocate memory for data*/
     data = (LSMLIB_REAL *)malloc((*num_elements)*DSZ);

     /* read data array */
     fread(data ,DSZ, *num_elements, fp);

     fclose(fp);
     zipFile(file_base,zip_status);
   }
   else
   {
      printf("\nCould not open file %s",file_name);
   }
   
   free(file_base);
   return data;
}

