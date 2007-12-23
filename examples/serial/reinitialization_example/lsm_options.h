/*
 * File:        lsm_options.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.2 $
 * Modified:    $Date: 2006/09/18 20:33:55 $
 * Description: Header file for Options structure implentation functions.
 */

#ifndef INCLUDED_LSM_OPTIONS
#define INCLUDED_LSM_OPTIONS

#ifdef __cplusplus
extern "C" {
#endif


/*! \file lsm_options.h
 *
 * \brief
 * LSMLIB provides support for Options structure.
 *
 */

#include "lsm_grid.h"


/*! 
 * Structure 'Options' stores various input and running options
 * that can be helpful when coding.
 * Options can be set by user via reading an input file, and each
 * option that is not mentioned in input file will be set to default.
 *   
 * Modify this structure and accompanying functions 
 * - setDefaultOptions
 * - copyOptions
 * - createOptionsFromInputFile
 * - printOptions
 * to add any additional options you need.
 *
 */
typedef struct _Options
{
      /* Main part of the structure */
   char   outfile[256]; /* output file name */
   char   path[256];    /* output file path, used internally */
     
   double dx;           /* grid spacing - assumed same in all dimensions */
   double tmax;         /* max running time */
   char   accuracy[20]; /* accuracy options: low, medium, high, very_high */
   int    accuracy_id;  /* internal accuracy identifier */
   
   /* Assuming phi_t + a |grad_phi| = b kappa |grad_phi|
      where kappa is mean curvature
   */   
   double a;   /* normal velocity */
   double b;   /* mean curvature term */
   
   int    save_data;    /* save all data (1) or not (0); data saved in
                           the same directory as the  ouput file */
   int    do_reinit;    /* reinitialize periodically (1) or not (0) */
   int    do_mask;      /* impose mask (restrict domain of movement) (1) 
                           or not (0); if 1 then set 'mask' array in
			   LSM_DataArrays structure */
   
   /* User additions */
   
   /* end User additions */

} Options;


/*!
 * Accuracy_settings structure should make it easier to switch between 
 * numerical discretization of different accuracy.
 * 'low'       <-->  HJ ENO1  in space, TVD Runge Kutta 1 in time
 * 'medium'    <-->  HJ ENO2  in space, TVD Runge Kutta 2 in time
 * 'high'      <-->  HJ ENO3  in space, TVD Runge Kutta 3 in time
 * 'very_high' <-->  HJ WENO5 in space, TVD Runge Kutta 3 in time
 */
typedef struct _Accuracy_settings_item {
   char *name;
   int  num_ghostcells;      /* number of ghostcells to be added to each
                                volume side */
} Accuracy_settings_item;


/*!
 * Accuracy_settings_menu connects the accuracy (set by user) with 
 * appropriate number of ghostcells to be added and the pointer
 * to the function that sets Grid for the appropriate numerical scheme.
 */
static Accuracy_settings_item Accuracy_settings_menu[] =
{
   {"low",       2 },
   {"medium",    3 },
   {"high",      5 },
   {"very_high", 4 }
};


/*================= Options structure manipulation ==================*/


/*!
 * createOptionsDefault() allocates a new Options structure and sets all 
 *  variables in Options structure to default values.
 *  
 * Arguments:          none
 *  
 * Return value:       pointer to an Options structure
 *
 * NOTES: Change this function if more Options structure elements are added.
 */
Options *createOptionsDefault(void);


/*!
 *  copyOptions() allocates a new Options structure and copies all 
 *  variables in Options structure from the provided Options structure.
 *  
 * Arguments:          
 *   - options_src(in):    Options structure to be copied
 *  
 * Return value:       pointer to an Options structure
 *
 * NOTES: Change this function if more Options structure elements are added.
 */
Options *copyOptions(Options *options_src);



/*! 'createOptionsFromInputFile' allocated memory for a new Options structure and sets its
 *         variables according to input file. Variables not set by input file
 *         are set to their default values.
 *
 *   Arguments:
 *   - filename (in):   input file name
 *   
 *   Return value:      pointer to a new Options structure 
 *
 *   NOTES: 
 *   - Input file is assumed to have a parameter entered in each line.
 *   - Appropriate keyword (equal to variable name from structure 'Options' )
 *     has to be at the beginnning of the line.
 *   - Change this function if more Options structure elements are added.
 */
Options *createOptionsFromInputFile(char *options);


/*!
 * 'printOptions' prints options structure variables to an output file
 *
 *   Arguments:
 *   - options(in):    Options structure to be printed
 *   - fp(in):         output file pointer (file opened for writing beforehand)
 *
 *   Return value:    none
 *
 *   NOTES: Change this function if more Options structure elements are added.
 */ 
void     printOptions(Options *options,FILE *fp);


/*!
 * destroyOptions() frees memory for an Options structure.
 *
 * Arguments:     
 *  - Options (in):   pointer to Options
 *
 * Return value:  none
 *
 */
void destroyOptions(Options *options);


/* User additions */

/* end User additions */

#ifdef __cplusplus
}
#endif

#endif
