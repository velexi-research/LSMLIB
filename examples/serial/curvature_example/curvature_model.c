/* System headers */
#include <stdio.h>

/* LSMLIB Serial package headers */
#include "lsm_data_arrays.h"
#include "lsm_macros.h"

/* Local headers */
#include "lsm_options.h"
#include "curvature_model_top.h"

/* Main driver for constant curvature level set method model */

int main(int argc, char **argv) 
{
   /* input filename storage */
   char  *in_fname, fname[256];
   char  *fname_data_in = (char *)NULL;
   char  *fname_grid_in = (char *)NULL;
   char  *fname_mask_in = (char *)NULL;
   
   int return_status;
   
   /* Structure containing input options and parameters.
      See lsm_options.h for details on the structure elements.
   */
   Options *options; 
   
   /* Initialize the problem */   
   if( argc == 1 )
   {   /* input file not provided, set all options to default */
       options =  createOptionsDefault();
   }
   else if (argc == 2)
   {   /* set options according to input file */
       in_fname = argv[1];
       sprintf(fname,"%s",argv[1]);
       options = createOptionsFromInputFile(fname);
   }
   else if( argc >= 4 )
   { /* read data from provided input files */
     fname_data_in = argv[2];
     fname_grid_in = argv[3];
     fname_mask_in = argv[4];
   }
   else
   {
     printf("\nRunning options:");
     printf("\n\t./constantCurvatureModel");
     printf("\n\t./constantCuvatureModel input_file");
     printf("\n\t./constantCurvatureModel input_file data_init grid mask");
     printf("\n"); 
   }
   
   return_status = curvatureModelTop(options,fname_data_in,fname_grid_in,
                                                               fname_mask_in);   
   /* clean up memory */
   destroyOptions(options);
   return return_status;
}
