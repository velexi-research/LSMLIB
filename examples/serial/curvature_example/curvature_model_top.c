/*
 * File:        curvature_model_top.c
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Support file for constant curvature flow.
 */
 
/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* LSMLIB Headers */
#include "LSMLIB_config.h"
#include "lsm_initialization3d.h"

/*LSMLIB Serial headers */
#include "lsm_macros.h"
#include "lsm_grid.h"

/* Local headers */
#include "curvature_model_top.h"
#include "curvature_model3d.h"
#include "curvature_model3d_local.h"

/* 
*  Top routine for constant curvature flow: sets input options and 
*  calls appropriate main loop function.
*/

int  curvatureModelTop(
     Options     *options,
     char        *fname_data_in,
     char        *fname_grid_in,
     char        *fname_mask_in)                   
{ 
  /* structure containing all arrays */
  LSM_DataArrays *data_arrays; 
  /* grid structure */
  Grid *grid;   
   /* time parameters */
  time_t   time0, time1;
 
  int     n1[3], n2[3], i;
  char    fname[256];
  FILE    *fp_out; 
  
  LSMLIB_REAL    normalx, normaly, normalz;
  LSMLIB_REAL    pointx, pointy, pointz;
  
  time(&time0);
  
  data_arrays = allocateLSMDataArrays();
    
  if( (fname_data_in == NULL) || (fname_grid_in == NULL) || 
      (fname_mask_in == NULL) )
  {
     /* Allocate/initialize mask and Grid structure.
        Mask is assumed to be a level set function that is negative in the
	pore space and positive elsewhere.
     */    
     grid = createMaskThroatFromSpheres3d(&(data_arrays->mask),options); 
     
     /* Allocate/initialize (planar) initial interface position */
     normalx = 1.0; normaly = 0.0; normalz = 0.0;
     pointx  = (grid->x_lo[0]) + 20 * (grid->dx)[0];
     pointy = pointz = 0.0; 
     
     data_arrays->phi = (LSMLIB_REAL *)calloc(grid->num_gridpts,sizeof(LSMLIB_REAL));
     createPlane(data_arrays->phi,normalx,normaly,normalz,pointx,pointy,
                                                                  pointz,grid);    					
     if(options->do_mask)
     {
        IMPOSE_MASK(data_arrays->phi,data_arrays->mask,data_arrays->phi,grid)
     }					
  }
  else
  {
     /* Read input data arrays and grid */
     data_arrays->phi = readDataArray(n1,fname_data_in);
     data_arrays->mask = readDataArray(n2,fname_mask_in);
     grid = readGridFromBinaryFile(fname_grid_in);
     
     /* Verify that input data files describe the same geometry */
     for(i = 0; i < 3; i++)
     {
        if( (n1[i] != (grid->grid_dims_ghostbox)[i] ) || 
	    (n2[i] != (grid->grid_dims_ghostbox)[i] ))
	{
	   printf("\nInput data dimensions n1[%d]=%d, n2[%d]=%d,",
	                                          i,n1[i],i,n2[i]);
	   printf(" (grid->grid_dims_ghostbox)[%d]=%d",
	                             i,(grid->grid_dims_ghostbox)[i]);
	   printf(" don't match\n");
	   printf("\nTerminating...");
	   
	   return 1;			  
	}
     }
  }   
  
  /* Open output file */
  fp_out = fopen(options->outfile,"w");
  if( options->print_details)
  {
    fprintf(fp_out,"Time Start %s\n",ctime(&time0));
    fprintf(fp_out,"Constant Curvature Model Level Set Method simulation\n");
      
    printOptions(options,fp_out);
    printGrid(grid,fp_out);
  }
  
  fflush(fp_out);
  
  /* If desired, output initial data and grid that can be visualized and/or 
     serve as an input for a different run
  */
  if( options->save_data )
  {
    sprintf(fname,"%s/data_init",options->path);
    writeDataArray(data_arrays->phi,grid,fname,GZIP);
    sprintf(fname,"%s/grid",options->path);
    writeGridToBinaryFile(grid,fname,GZIP);
    sprintf(fname,"%s/mask",options->path);
    writeDataArray(data_arrays->mask,grid,fname,GZIP);
  }
  
  setArrayAllocationCurvatureModel(options,data_arrays);
  allocateMemoryForLSMDataArrays(data_arrays,grid); 
  
  /* Run the curvature model, only 3d supported so far */
  if( grid->num_dims == 3 )
  {
    if(options->narrow_band)
      curvatureModelMedium3dLocalMainLoop(options,data_arrays,grid,fp_out);
    else
      curvatureModelMedium3dMainLoop(options,data_arrays,grid,fp_out);								
  }
         
  /* If desired, output initial data */
  if( options->save_data )
  {
    sprintf(fname,"%s/data_final",options->path);
    writeDataArray(data_arrays->phi,grid,fname,GZIP);    
  }  
 
  /* Clean up memory */
  freeMemoryForLSMDataArrays(data_arrays);
  destroyGrid(grid);
 
  /* Print out time statistics */
  time(&time1);
  fprintf(fp_out,"\nTime End %s",ctime(&time1));
  fprintf(fp_out,"\nExecution time %g  seconds.\n", difftime(time1,time0));
  fflush(fp_out);
  
  return 0;
}


 /*  setArrayAllocationCurvatureModel()
  *  decides which arrays not to allocate by setting their pointer to a
  *  NULL value
  *
  *  This is optional and saves some memory. 
  */
void setArrayAllocationCurvatureModel(
     Options *options,
     LSM_DataArrays *data_arrays)
{  
    if(options->narrow_band == 0)
    { /* No need for narrow band (localization) storage */
       data_arrays->narrow_band = (unsigned char *)NULL;
       data_arrays->index_x = (int *)NULL;
       data_arrays->index_y = (int *)NULL;
       data_arrays->index_z = (int *)NULL;
       data_arrays->index_outer_pts = (int *)NULL;
    }
        
    if(options->b == 0)
    { /* Second order derivatives will (presumably) not be used */
       data_arrays->phi_xx = data_arrays->phi_xy = data_arrays->phi_yy = (LSMLIB_REAL *)NULL;
       data_arrays->phi_zz = data_arrays->phi_xz = data_arrays->phi_yz = (LSMLIB_REAL *)NULL;
    }
        
    if(options->a == 0)
    { /* Upwinding derivatives will (presumably) not be used */
       data_arrays->phi_x_minus  = data_arrays->phi_x_plus = (LSMLIB_REAL *)NULL;
       data_arrays->phi_y_minus  = data_arrays->phi_y_plus = (LSMLIB_REAL *)NULL;
       data_arrays->phi_z_minus  = data_arrays->phi_z_plus = (LSMLIB_REAL *)NULL;
    }

    if( options->accuracy_id != 2 )
    { /* 3rd order differences used only for HJ ENO3 */
       data_arrays->D3 = (LSMLIB_REAL *)NULL;
    }

    if( options->accuracy_id <= MEDIUM )
    {
       data_arrays->phi_stage2 = (LSMLIB_REAL *)NULL;
    } 

    if( options->do_reinit == 0)
    {
       data_arrays->phi0 = (LSMLIB_REAL *)NULL;
    }

    /* Curvature model does not assume external velocity */
    data_arrays->external_velocity_x = (LSMLIB_REAL *)NULL;
    data_arrays->external_velocity_y = (LSMLIB_REAL *)NULL;
    data_arrays->external_velocity_z = (LSMLIB_REAL *)NULL;
}

/*  createMaskThroatFromSpheres3d()
*   Creates mask for a pore space throat enclosed by three spheres of
*   radius 1.0.
*   
*   Arguments:
*   pmask(out): pointer to the array to contain masking level set function
*           mask will be negative in pore space and positive within
*	   the spheres
*   options - Options structure, using grid spacing dx and accuracy_id elements
*   
*   Returns:
*    pointer to Grid structure
*    
*   NOTES: - memory for mask level set function is allocated within  
*	   
*/
 Grid *createMaskThroatFromSpheres3d(
      LSMLIB_REAL  **pmask,
      Options  *options)
 {
      LSMLIB_REAL  x_lo[3], x_hi[3];
      int     dim;
      Grid    *grid;
      
      LSMLIB_REAL  centerx[3], centery[3], centerz[3], radius[3];
      LSMLIB_REAL  *mask, r;
      int     n, inside_flag[3];
            
      dim = 3;
      r = 1.0;
      x_lo[0] = -3*r/4;    x_hi[0] = 3*r/4;
      x_lo[1] = -3*r/4;    x_hi[1] = 3*r/4;
      x_lo[2] =  0;        x_hi[2] = r;
      grid = createGridSetDx(dim,options->dx,x_lo,x_hi,options->accuracy_id);
      	 
      mask = (LSMLIB_REAL*) malloc((grid->num_gridpts)*sizeof(LSMLIB_REAL));   
      
      n = 3; /* number of spheres */
      
      /* set sphere centers */
      centerx[0] = 0; centery[0] = -r;  centerz[0] = 0; 
      centerx[1] = 0; centery[1] =  r;  centerz[1] = 0; 
      centerx[2] = 0; centery[2] =  0;  centerz[2] = r*sqrt(3.0);
     
      /* all sphere are of the same radius r */
      radius[0] = radius[1] = radius[2] = r;
      /* inside of spheres will be set to positive values */
      inside_flag[0] = inside_flag[1] = inside_flag[2] = 1;
      
      createIntersectionOfSpheres(mask,n,centerx,centery,centerz,
                                            radius,inside_flag,grid);
      *pmask = mask;
      return grid;
 }
