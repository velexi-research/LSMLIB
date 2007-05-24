// System headers
#include <stdio.h>
#include <time.h>

#include "lsm_data_arrays.h"
#include "lsm_macros.h"
#include "lsm_options.h"
#include "lsm_initialization2d.h"

#include "lsm_reinitialization_medium2d.h"

void reinitializeTop(Options *options,int use_input_data,char *filename_data,
                                       char *filename_grid,char *filename_mask);
      
int main(int argc, char **argv) 
{
   char filename_input[256];
   char filename_data[256], filename_grid[256], filename_mask[256];
   int  use_input_data = 0;
   
   /* structure containing input options and parameters*/
   Options *options; 
   
   
   if( argc == 1 )
   {   /* input file not provided, set all options to default */
       options =  createOptionsDefault();
       options->tmax = 5*options->dx; 
   }
   else
   {   /* set options according to input file */
       sprintf(filename_input,"%s",argv[1]);
       options = createOptionsFromInputFile(filename_input);
       
       if( argc >= 4 )
       { /* initial level set function, grid are read from input files */
          sprintf(filename_data,"%s",argv[2]); /* initial data */
          sprintf(filename_grid,"%s",argv[3]); /* grid corresponding to initial 
	                                         data */
	  if( argc == 5 ) /* providing mask is optional */					 
          {
	    sprintf(filename_mask,"%s",argv[4]);
	    options->do_mask = 1;
	  }
	                                          			      
          use_input_data = 1;
       }
   }
   
   reinitializeTop(options,use_input_data,filename_data,filename_grid,
                                        filename_mask);
   
   return 0;
}



void reinitializeTop(
      Options *options,
      int      use_input_data,
      char     *filename_data,
      char     *filename_grid,
      char     *filename_mask)                   
{ 
  /* structure containing all arrays */
  LSM_DataArrays *p; 
  /* grid structure */
  Grid *g;   
   /* time parameters */
  time_t   time0, time1;
 
  double one = 1.0;
  int    idx, tmp_grid_dims1[3], tmp_grid_dims2[3], i;
  char    filename[256];
  FILE    *fp_out;                               
  
  time(&time0);
  
  /* allocate LSM_DataArrays structure */
  p = allocateLSMDataArrays();
  
  if( use_input_data == 0 )
  {
     /* set a default problem */
     int  num_gc = Accuracy_settings_menu[options->accuracy_id].num_ghostcells;
     double x_lo[2], x_hi[2];
     double centerx, centery, radius;
     int    dim, inside_flag, i;   
      
     dim = 2;
     x_lo[0] =  0;    x_hi[0] = 1.0;
     x_lo[1] =  0;    x_hi[1] = 1.0;
     
     /* allocate and set the grid structure */
     g = createGridSetDx(dim,options->dx,x_lo,x_hi,
                (LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE) options->accuracy_id);
     
     centerx = 0.5; centery = 0.5; radius = 0.3;
     /* allocate memory for the level set function */
     p->phi = (double *)malloc(g->num_gridpts*sizeof(double));
     
     /* set the level set function */
     inside_flag = -1; /* inside of the circle we'll have negative values */
     createCircle(p->phi,centerx,centery,radius,inside_flag,g);
  
     /* The level set created is a signed distance function, i.e. "smooth"
     *  We will convert it to a segmented form, so the reinitialization
     *  effects can be seen. 
     */
     
     for(i=0; i <g->num_gridpts; i++)
     {
             if( p->phi[i] < 0 ) p->phi[i] = -options->dx;
	else if( p->phi[i] > 0 ) p->phi[i] =  options->dx;
     } 
  }
  else
  {
     /* Note that mask provided this way can be non-trivial */
     p->phi = readDataArray(tmp_grid_dims1,filename_data);
     g = readGridFromBinaryFile(filename_grid);
     
     if(options->do_mask) p->mask = readDataArray(tmp_grid_dims2,filename_mask);
     
     for(i = 0; i < 3; i++)
     {
        if( (tmp_grid_dims1[i] != (g->grid_dims_ghostbox)[i] ))
	{
	   printf("\nInput data dimensions tmp_grid_dims1[%d]=%d, ",
	                                                 i,tmp_grid_dims1[i]);
	   printf("\n(g->grid_dims_ghostbox)[%d]=%d don't match\n",
                                                i,(g->grid_dims_ghostbox)[i]);				  
	}
	
	if( options->do_mask && (tmp_grid_dims2[i] != (g->grid_dims_ghostbox)[i] ))
	{
	   printf("\nInput data dimensions tmp_grid_dims2[%d]=%d, ",
	                                                 i,tmp_grid_dims2[i]);
	   printf("\n(g->grid_dims_ghostbox)[%d]=%d don't match\n",
                                                i,(g->grid_dims_ghostbox)[i]);				  
	}
     }   
  }
    
  fp_out = fopen(options->outfile,"w");
  fprintf(fp_out,"Time Start %s\n",ctime(&time0));
  fprintf(fp_out,"Level Set Method Reinitialization\n");

  printOptions(options,fp_out);
  printGrid(g,fp_out);
  
  fflush(fp_out);
  
  /* Memory allocation for arrays. 
   * Note that all of the arrays will be allocated, though some memory
   * allocations could be prevented by setting appropriate array pointers 
   * to NULL.
   */   
  allocateMemoryForLSMDataArrays(p,g); 
    
  /* save data if desired */
  if( options->save_data )
  {
    sprintf(filename,"%sdata_init",options->path);
    writeDataArray(p->phi,g,filename);
    fprintf(fp_out,"\nInitial level set function output to binary file %s",
                                                                     filename);
    sprintf(filename,"%sgrid",options->path);
    writeGridToBinaryFile(g,filename);
    fprintf(fp_out,"\nGrid structure output to binary file named %s",
                                                                     filename);
  }  
  
  /* reinitialization available only for medium accuracy at the moment */
  lsm2dReinitializationMedium(p,g,options);
  
  if( options->save_data )
  {
    sprintf(filename,"%sdata_final",options->path);
    writeDataArray(p->phi,g,filename);
    fprintf(fp_out,"\nFinal level set function output to binary file %s",
                                                                     filename);
    if(options->do_mask)
    {
      sprintf(filename,"%smask",options->path);
      writeDataArray(p->mask,g,filename);
      fprintf(fp_out,"\nMask level set function output to binary file %s",
                                                                     filename);
    }  
  }  
 
  /* clean up memory */
  freeMemoryForLSMDataArrays(p);
  destroyGrid(g);
 
  time(&time1);
  fprintf(fp_out,"\nTime End %s",ctime(&time1));
  fprintf(fp_out,"\nExecution time %g  seconds.\n", difftime(time1,time0));
  fclose(fp_out);
}
