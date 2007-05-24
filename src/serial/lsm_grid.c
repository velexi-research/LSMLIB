/*
 * File:        lsm_grid.c
 * Copyright:   (c) 2005-2006 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.5 $
 * Modified:    $Date: 2006/09/18 20:38:08 $
 * Description: Implementation file for grid data structures that support 
 *              serial LSMLIB calculations
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lsm_grid.h"

/*======== Helper Functions for Grid structure manipulation ========*/

/*
 * allocateGrid() allocates memory for a Grid structure
 *
 * Arguments:     none
 *
 * Return value:  pointer to newly allocated Grid structure
 *
 */
Grid *allocateGrid(void);
 
   
/*================= Grid structure manipulation ==================*/

Grid *allocateGrid(void)
{
  Grid *p_grid = (Grid *)calloc(1,sizeof(Grid));
  return  p_grid;
}

/* 
*  number of ghostcells for spatial discretization corresponding to 
*  the desired accuracy 
*/

static int lsmlib_num_ghostcells[] = {2,3,5,4};


/*======== Helper Functions for setting Grid index space limits  ========*/

/*
 * setIndexSpaceLimitsCentral() sets the upper and lower limits of the 
 * specified Grid for computing spatial derivatives using central 
 * differencing. 
 *
 * Arguments:
 *  - grid (in/out):  Grid data structure containing grid configuration
 *
 * Return value:      none     
 *
 * NOTES: 
 * - grid is assumed to be allocated by the user.
 *
 * - setIndexSpaceLimitsCentral() only sets the relevant data fields within 
 *   grid (i.e. it does not set ALL of the data fields of grid).
 *
 * - setIndexSpaceLimitsCentral() handles both 2D and 3D cases.  For 2D
 *   problems, the third dimension bounds might be set to 
 *   negative numbers, but these values should never be used.
 *
 */
void  setIndexSpaceLimitsCentral(Grid *grid)
{
 /* max bounding grid ('ghost grid') in each dimension */
  grid->ilo_gb = 0;   grid->ihi_gb = (grid->grid_dims_ghostbox)[0]-1;
  grid->jlo_gb = 0;   grid->jhi_gb = (grid->grid_dims_ghostbox)[1]-1;
  grid->klo_gb = 0;   grid->khi_gb = (grid->grid_dims_ghostbox)[2]-1;
  
  /* 1st order derivatives (central differencing) fill grid */
  grid->ilo_D1_fb = 1;   grid->ihi_D1_fb = (grid->grid_dims_ghostbox)[0]-2;
  grid->jlo_D1_fb = 1;   grid->jhi_D1_fb = (grid->grid_dims_ghostbox)[1]-2;
  grid->klo_D1_fb = 1;   grid->khi_D1_fb = (grid->grid_dims_ghostbox)[2]-2;
  
  /* 2nd order  derivatives (central differencing) fill grid */
  grid->ilo_D2_fb = 2;   grid->ihi_D2_fb = (grid->grid_dims_ghostbox)[0]-3;
  grid->jlo_D2_fb = 2;   grid->jhi_D2_fb = (grid->grid_dims_ghostbox)[1]-3;
  grid->klo_D2_fb = 2;   grid->khi_D2_fb = (grid->grid_dims_ghostbox)[2]-3;
  
  /* 3rd order  derivatives (central differencing) fill grid */
  grid->ilo_D3_fb = 3;   grid->ihi_D3_fb = (grid->grid_dims_ghostbox)[0]-4;
  grid->jlo_D3_fb = 3;   grid->jhi_D3_fb = (grid->grid_dims_ghostbox)[1]-4;
  grid->klo_D3_fb = 3;   grid->khi_D3_fb = (grid->grid_dims_ghostbox)[2]-4;
}  


/*
 * setIndexSpaceLimitsENO1() sets the upper and lower limits of the 
 * specified Grid for computing spatial derivatives using the HJ ENO1 
 * scheme.
 *
 * Arguments:
 *  - grid (in/out):  Grid data structure containing grid configuration
 *
 * Return value:      none     
 *
 * NOTES: 
 * - grid is assumed to be allocated by the user.
 *
 * - setIndexSpaceLimitsENO1() handles both 2D and 3D cases.  For 2D
 *   problems, the third dimension bounds might be set to negative 
 *   numbers, but these values should never be used.
 *
 */   
void setIndexSpaceLimitsENO1(Grid *grid)
{
  int ilo_plus_fb, ihi_plus_fb, ilo_minus_fb, ihi_minus_fb;
  int jlo_plus_fb, jhi_plus_fb, jlo_minus_fb, jhi_minus_fb;
  int klo_plus_fb, khi_plus_fb, klo_minus_fb, khi_minus_fb;
  
  setIndexSpaceLimitsCentral(grid);
  
  /* upwind derivatives by HJ ENO2 */
  ilo_plus_fb  = grid->ilo_D1_fb;     ihi_plus_fb  = grid->ihi_D1_fb-1;
  jlo_plus_fb  = grid->jlo_D1_fb;     jhi_plus_fb  = grid->jhi_D1_fb-1;
  klo_plus_fb  = grid->klo_D1_fb;     khi_plus_fb  = grid->khi_D1_fb-1;
  
  /* TO DO - RECHECK */
  ilo_minus_fb = grid->ilo_D1_fb;   ihi_minus_fb = grid->ihi_D1_fb;
  jlo_minus_fb = grid->jlo_D1_fb;   jhi_minus_fb = grid->jhi_D1_fb;
  klo_minus_fb = grid->klo_D1_fb;   khi_minus_fb = grid->khi_D1_fb;
  
  /* the common fill grid for phi_plus and phi_minus */  
  grid->ilo_fb = ilo_minus_fb;   grid->ihi_fb = ihi_plus_fb;
  grid->jlo_fb = jlo_minus_fb;   grid->jhi_fb = jhi_plus_fb;
  grid->klo_fb = klo_minus_fb;   grid->khi_fb = khi_plus_fb;
}


/*
 * setIndexSpaceLimitsENO2() sets the upper and lower limits of the 
 * specified Grid for computing spatial derivatives using the HJ ENO2 
 * scheme.
 *
 * Arguments:
 *  - grid (in/out):  Grid data structure containing grid configuration
 *
 * Return value:      none     
 *
 * NOTES: 
 * - grid is assumed to be allocated by the user.
 *
 * - setIndexSpaceLimitsENO2() handles both 2D and 3D cases.  For 2D
 *   problems, the third dimension bounds might be set to negative 
 *   numbers, but these values should never be used.
 *
 */   
void setIndexSpaceLimitsENO2(Grid *grid)
{
  int ilo_plus_fb, ihi_plus_fb, ilo_minus_fb, ihi_minus_fb;
  int jlo_plus_fb, jhi_plus_fb, jlo_minus_fb, jhi_minus_fb;
  int klo_plus_fb, khi_plus_fb, klo_minus_fb, khi_minus_fb;
  
  setIndexSpaceLimitsCentral(grid);
  
  /* upwind derivatives by HJ ENO 2nd order */
  ilo_plus_fb  = grid->ilo_D2_fb;     ihi_plus_fb  = grid->ihi_D2_fb-1;
  jlo_plus_fb  = grid->jlo_D2_fb;     jhi_plus_fb  = grid->jhi_D2_fb-1;
  klo_plus_fb  = grid->klo_D2_fb;     khi_plus_fb  = grid->khi_D2_fb-1;
  
  ilo_minus_fb = grid->ilo_D2_fb+1;   ihi_minus_fb = grid->ihi_D2_fb;
  jlo_minus_fb = grid->jlo_D2_fb+1;   jhi_minus_fb = grid->jhi_D2_fb;
  klo_minus_fb = grid->klo_D2_fb+1;   khi_minus_fb = grid->khi_D2_fb;
  
  /* the common fill grid for phi_plus and phi_minus */  
  grid->ilo_fb = ilo_minus_fb;   grid->ihi_fb = ihi_plus_fb;
  grid->jlo_fb = jlo_minus_fb;   grid->jhi_fb = jhi_plus_fb;
  grid->klo_fb = klo_minus_fb;   grid->khi_fb = khi_plus_fb;
}


/*  
 * setIndexSpaceLimitsENO3() sets the upper and lower limits of the 
 * specified Grid for computing spatial derivatives using the HJ ENO3 
 * scheme.
 *
 * Arguments:
 *  - grid (in/out):  Grid data structure containing grid configuration
 *
 * Return value:      none     
 *
 * NOTES: 
 * - grid is assumed to be allocated by the user.
 *
 * - setIndexSpaceLimitsENO3() only sets the relevant data fields within grid
 *   (i.e. it does not set ALL of the data fields of grid).
 *
 * - setIndexSpaceLimitsENO3() handles both 2D and 3D cases.  For 2D
 *   problems, the third dimension bounds might be set to 
 *   negative numbers, but these values should never be used.
 *
 */   
void setIndexSpaceLimitsENO3(Grid *grid)
{
  int ilo_plus_fb, ihi_plus_fb, ilo_minus_fb, ihi_minus_fb;
  int jlo_plus_fb, jhi_plus_fb, jlo_minus_fb, jhi_minus_fb;
  int klo_plus_fb, khi_plus_fb, klo_minus_fb, khi_minus_fb;
  
  setIndexSpaceLimitsCentral(grid);
  
  /* upwind derivatives by HJ ENO3 */
  ilo_plus_fb  = grid->ilo_D3_fb;     ihi_plus_fb  = grid->ihi_D3_fb-2;
  jlo_plus_fb  = grid->jlo_D3_fb;     jhi_plus_fb  = grid->jhi_D3_fb-2;
  klo_plus_fb  = grid->klo_D3_fb;     khi_plus_fb  = grid->khi_D3_fb-2;
  
  ilo_minus_fb = grid->ilo_D3_fb+1;   ihi_minus_fb = grid->ihi_D3_fb;
  jlo_minus_fb = grid->jlo_D3_fb+1;   jhi_minus_fb = grid->jhi_D3_fb;
  klo_minus_fb = grid->klo_D3_fb+1;   khi_minus_fb = grid->khi_D3_fb;
  
  /* the common fill grid for phi_plus and phi_minus */  
  grid->ilo_fb = ilo_minus_fb;   grid->ihi_fb = ihi_plus_fb;
  grid->jlo_fb = jlo_minus_fb;   grid->jhi_fb = jhi_plus_fb;
  grid->klo_fb = klo_minus_fb;   grid->khi_fb = khi_plus_fb;
}


/*
 * setIndexSpaceLimitsWENO5() sets the upper and lower limits of the 
 * specified Grid for computing spatial derivatives using the HJ WENO5 
 * scheme.
 *
 * Arguments:
 *  - grid (in/out):  Grid data structure containing grid configuration
 *
 * Return value:      none     
 *
 * NOTES: 
 * - grid is assumed to be allocated by the user.
 *
 * - setIndexSpaceLimitsWENO5() only sets the relevant data fields within grid
 *   (i.e. it does not set ALL of the data fields of grid).
 *
 * - setIndexSpaceLimitsWENO5() handles both 2D and 3D cases.  For 2D
 *   problems, the third dimension bounds might be set to 
 *   negative numbers, but these values should never be used.
 *
 */   
void setIndexSpaceLimitsWENO5(Grid *grid)
{
  int ilo_plus_fb, ihi_plus_fb, ilo_minus_fb, ihi_minus_fb;
  int jlo_plus_fb, jhi_plus_fb, jlo_minus_fb, jhi_minus_fb;
  int klo_plus_fb, khi_plus_fb, klo_minus_fb, khi_minus_fb;
  
  setIndexSpaceLimitsCentral(grid);
  
  /* upwind derivatives for HJ WENO 5th order */
  ilo_plus_fb  = grid->ilo_D1_fb;     ihi_plus_fb  = grid->ihi_D1_fb-3;
  jlo_plus_fb  = grid->jlo_D1_fb;     jhi_plus_fb  = grid->jhi_D1_fb-3;
  klo_plus_fb  = grid->klo_D1_fb;     khi_plus_fb  = grid->khi_D1_fb-3;
  
  ilo_minus_fb = grid->ilo_D1_fb+2;   ihi_minus_fb = grid->ihi_D1_fb;
  jlo_minus_fb = grid->jlo_D1_fb+2;   jhi_minus_fb = grid->jhi_D1_fb;
  klo_minus_fb = grid->klo_D1_fb+2;   khi_minus_fb = grid->khi_D1_fb;
  
  /* the common fill grid for phi_plus and phi_minus */  
  grid->ilo_fb = ilo_minus_fb;   grid->ihi_fb = ihi_plus_fb;
  grid->jlo_fb = jlo_minus_fb;   grid->jhi_fb = jhi_plus_fb;
  grid->klo_fb = klo_minus_fb;   grid->khi_fb = khi_plus_fb;
}


/*============= Function definitions for Grid management ==============*/

Grid *createGridSetDx(
      int    num_dims,
      double dx,
      double *x_lo,
      double *x_hi,
      LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE accuracy)
{
  int     i;
  double  diff;

  Grid    *g;
  
  int num_ghostcells = lsmlib_num_ghostcells[accuracy];
  
  g = allocateGrid();
  
  if( (num_dims == 2) || (num_dims == 3) )
  {
    g->num_dims = num_dims;
  }
  else
  {
    printf("\nInvalid grid dimension %d, set to default 2.\n", num_dims);
    g->num_dims = 2;
  }  
  
  /* set up grid */
  g->num_gridpts = 1;

  for(i = 0; i < g->num_dims; i++)
  {
    (g->dx)[i] = dx;

    /* compute grid dimensions */
    diff = ( x_hi[i] - x_lo[i] );
    g->grid_dims[i] = rint( diff/(dx) );

    /* set x_lo and x_hi consistent with dx and grid_dims */
    g->x_lo[i] = x_lo[i];
    g->x_hi[i] = x_lo[i] + (g->grid_dims)[i]*dx;

    /* add 'num_ghostcells' layers of voxels on each side */
    (g->grid_dims_ghostbox)[i] =  (g->grid_dims)[i] + 2*num_ghostcells;
    (g->x_lo_ghostbox)[i]= (g->x_lo)[i] - num_ghostcells*((g->dx)[i]);
    (g->x_hi_ghostbox)[i]= (g->x_hi)[i] + num_ghostcells*((g->dx)[i]);

    (g->num_gridpts) *= (g->grid_dims_ghostbox)[i];
  } 

  if( g->num_dims == 2 )
  {  /* put fake values for the third dimension */
    (g->x_lo)[2] = 0;
    (g->x_hi)[2] = 0;
    (g->x_lo_ghostbox)[2] = 0;
    (g->x_hi_ghostbox)[2] = 0;
    (g->grid_dims)[2] = 1;
    (g->grid_dims_ghostbox)[2] = 1;
    (g->dx)[2] = 0;
  }
  
  setIndexSpaceLimits(accuracy,g);
  
  return g;
}


Grid *createGridSetGridDims(
  int     num_dims,
  int    *grid_dims,
  double *x_lo,
  double *x_hi,
  LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE accuracy)
{
  Grid    *g;
  int      i;
  
  int num_ghostcells = lsmlib_num_ghostcells[accuracy];
  
  g = allocateGrid();

  /* set number of dimensions for grid */ 
  if( (num_dims == 2) || (num_dims == 3) )
  {
    g->num_dims = num_dims;
  }
  else
  {
    fprintf(stderr, "\nInvalid number of dimensions (%d).  ", num_dims);
    fprintf(stderr, "Using to default value of 2.\n");
    g->num_dims = 2;
  }  
  
  /* set grid parameters */
  g->num_gridpts = 1;
  for(i = 0; i < g->num_dims; i++)
  {
    (g->grid_dims)[i] = grid_dims[i];

    g->x_hi[i] = x_hi[i];
    g->x_lo[i] = x_lo[i];
  
    (g->dx)[i] = ( (g->x_hi)[i] - (g->x_lo)[i] )/(g->grid_dims)[i];

    /* add 'num_ghostcells' layers of voxels on each side */
    (g->grid_dims_ghostbox)[i] = 
      (g->grid_dims)[i] + 2*num_ghostcells;
    (g->x_lo_ghostbox)[i]= (g->x_lo)[i] - num_ghostcells*((g->dx)[i]);
    (g->x_hi_ghostbox)[i]= (g->x_hi)[i] + num_ghostcells*((g->dx)[i]);

    (g->num_gridpts) *= (g->grid_dims_ghostbox)[i];
  } 

  /* for 2D problems, use fake values for the third dimension */
  if( g->num_dims == 2 )
  {  
    (g->x_lo)[2] = 0;
    (g->x_hi)[2] = 0;
    (g->x_lo_ghostbox)[2] = 0;
    (g->x_hi_ghostbox)[2] = 0;
    (g->grid_dims)[2] = 1;
    (g->grid_dims_ghostbox)[2] = 1;
    (g->dx)[2] = 0;
  }
  
  setIndexSpaceLimits(accuracy,g);
  
  return g;
}


Grid *copyGrid(Grid *grid)
{
   Grid *new_grid;
   int   i;
   
   new_grid = allocateGrid();
   
   new_grid->num_dims = grid->num_dims;
   
   for(i = 0; i < 3; i++)
   {
       new_grid->x_lo[i] = grid->x_lo[i];
       new_grid->x_hi[i] = grid->x_hi[i];
       new_grid->x_lo_ghostbox[i] = grid->x_lo_ghostbox[i];
       new_grid->x_hi_ghostbox[i] = grid->x_hi_ghostbox[i];
       new_grid->grid_dims[i] = grid->grid_dims[i];
       new_grid->grid_dims_ghostbox[i] = grid->grid_dims_ghostbox[i];
       new_grid->dx[i] = grid->dx[i];      
   }
    new_grid->num_gridpts = grid->num_gridpts;
    
   new_grid->ilo_gb = grid->ilo_gb; new_grid->ihi_gb = grid->ihi_gb;
   new_grid->jlo_gb = grid->jlo_gb; new_grid->jhi_gb = grid->jhi_gb;
   new_grid->klo_gb = grid->klo_gb; new_grid->khi_gb = grid->khi_gb;
   
   new_grid->ilo_fb = grid->ilo_fb; new_grid->ihi_fb = grid->ihi_fb;
   new_grid->jlo_fb = grid->jlo_fb; new_grid->jhi_fb = grid->jhi_fb;
   new_grid->klo_fb = grid->klo_fb; new_grid->khi_fb = grid->khi_fb;
   
   new_grid->ilo_D1_fb = grid->ilo_D1_fb; 
   new_grid->ihi_D1_fb = grid->ihi_D1_fb;
   new_grid->jlo_D1_fb = grid->jlo_D1_fb; 
   new_grid->jhi_D1_fb = grid->jhi_D1_fb;
   new_grid->klo_D1_fb = grid->klo_D1_fb;
   new_grid->khi_D1_fb = grid->khi_D1_fb;
   
   new_grid->ilo_D2_fb = grid->ilo_D2_fb;
   new_grid->ihi_D2_fb = grid->ihi_D2_fb;
   new_grid->jlo_D2_fb = grid->jlo_D2_fb;
   new_grid->jhi_D2_fb = grid->jhi_D2_fb;
   new_grid->klo_D2_fb = grid->klo_D2_fb;
   new_grid->khi_D2_fb = grid->khi_D2_fb;
   
   new_grid->ilo_D3_fb = grid->ilo_D3_fb;
   new_grid->ihi_D3_fb = grid->ihi_D3_fb;
   new_grid->jlo_D3_fb = grid->jlo_D3_fb;
   new_grid->jhi_D3_fb = grid->jhi_D3_fb;
   new_grid->klo_D3_fb = grid->klo_D3_fb;
   new_grid->khi_D3_fb = grid->khi_D3_fb;
   
   return new_grid;
}


void destroyGrid(Grid *grid)
{
  if (grid) free(grid);
}
    

void printGrid(Grid *grid, FILE *fp)
{
  fprintf(fp, "Number of dimensions: %d\n", grid->num_dims);
  if (grid->num_dims == 2) 
  {
    fprintf(fp,
            "Geometric limits (without ghostcells): [%g,%g] x [%g,%g]\n",
            (grid->x_lo)[0], (grid->x_hi)[0],
            (grid->x_lo)[1], (grid->x_hi)[1]);
    fprintf(fp,
            "Geometric limits (with ghostcells): [%g,%g] x [%g,%g]\n",
            (grid->x_lo_ghostbox)[0], (grid->x_hi_ghostbox)[0],
            (grid->x_lo_ghostbox)[1], (grid->x_hi_ghostbox)[1]);
    fprintf(fp,"Grid size (without ghostcells): %d x %d\n", 
            (grid->grid_dims)[0], 
            (grid->grid_dims)[1]);
    fprintf(fp,"Grid size (with ghostcells): %d x %d\n", 
            (grid->grid_dims_ghostbox)[0],
            (grid->grid_dims_ghostbox)[1]);
    fprintf(fp,"Grid spacing: dx = %g, dy = %g \n",
            (grid->dx)[0], (grid->dx)[1]);   	    
  } 
  else 
  {
    fprintf(fp,
            "Geometric limits (without ghostcells): [%g,%g] x [%g,%g] x [%g,%g]\n",
            (grid->x_lo)[0], (grid->x_hi)[0],
            (grid->x_lo)[1], (grid->x_hi)[1],
            (grid->x_lo)[2], (grid->x_hi)[2]);
    fprintf(fp,
            "Geometric limits (with ghostcells): [%g,%g] x [%g,%g] x [%g,%g]\n",
            (grid->x_lo_ghostbox)[0], (grid->x_hi_ghostbox)[0],
            (grid->x_lo_ghostbox)[1], (grid->x_hi_ghostbox)[1],
            (grid->x_lo_ghostbox)[2], (grid->x_hi_ghostbox)[2]);
    fprintf(fp,"Grid size (without ghostcells): %d x %d x %d\n", 
            (grid->grid_dims)[0],
            (grid->grid_dims)[1],
            (grid->grid_dims)[2]);
    fprintf(fp,"Grid size (with ghostcells): %d x %d x %d\n", 
            (grid->grid_dims_ghostbox)[0],
            (grid->grid_dims_ghostbox)[1],
            (grid->grid_dims_ghostbox)[2]);
    fprintf(fp,"Grid spacing: dx = %g, dy = %g, dz = %g \n",
            (grid->dx)[0], (grid->dx)[1], (grid->dx)[2]);
  }
  fprintf(fp, "Total number of grid points: %d\n", grid->num_gridpts);
  
  /* fprintf(fp,"Index space limits\n"); */
  
  if (grid->num_dims == 2) 
  {
    fprintf(fp,
            "Ghost box index limits: [%d,%d] x [%d,%d]\n",
            grid->ilo_gb, grid->ihi_gb,
            grid->jlo_gb, grid->jhi_gb);
	    
    fprintf(fp,
            "Fill box index limits: [%d,%d] x [%d,%d]\n",
            grid->ilo_fb, grid->ihi_fb,
            grid->jlo_fb, grid->jhi_fb);
	    
    fprintf(fp,
            "1st order differences fill box: [%d,%d] x [%d,%d]\n",
            grid->ilo_D1_fb, grid->ihi_D1_fb,
            grid->jlo_D1_fb, grid->jhi_D1_fb);
    fprintf(fp,
            "2nd order differences fill box: [%d,%d] x [%d,%d]\n",
            grid->ilo_D2_fb, grid->ihi_D2_fb,
            grid->jlo_D2_fb, grid->jhi_D2_fb);
    fprintf(fp,
            "3rd order differences fill box: [%d,%d] x [%d,%d]\n",
            grid->ilo_D3_fb, grid->ihi_D3_fb,
            grid->jlo_D3_fb, grid->jhi_D3_fb);	 	    	     	    	    	    
  }
  else
  {
    fprintf(fp,
            "Ghost box index limits: [%d,%d] x [%d,%d] x [%d,%d]\n",
            grid->ilo_gb, grid->ihi_gb,
            grid->jlo_gb, grid->jhi_gb,
	    grid->klo_gb, grid->khi_gb);
    fprintf(fp,
            "Fill box index limits: [%d,%d] x [%d,%d] x [%d,%d]\n",
            grid->ilo_fb, grid->ihi_fb,
            grid->jlo_fb, grid->jhi_fb,
	    grid->klo_fb, grid->khi_fb);
	    	    
    fprintf(fp,
            "1st order differences fill box: [%d,%d] x [%d,%d] x [%d,%d]\n",
            grid->ilo_D1_fb, grid->ihi_D1_fb,
            grid->jlo_D1_fb, grid->jhi_D1_fb,
	    grid->klo_D1_fb, grid->khi_D1_fb);
	    
    fprintf(fp,
            "2nd order differences fill box: [%d,%d] x [%d,%d] x [%d,%d]\n",
            grid->ilo_D2_fb, grid->ihi_D2_fb,
            grid->jlo_D2_fb, grid->jhi_D2_fb,
	    grid->klo_D2_fb, grid->khi_D2_fb);
	    
    fprintf(fp,
            "3rd order differences fill box: [%d,%d] x [%d,%d] x [%d,%d]\n",
            grid->ilo_D3_fb, grid->ihi_D3_fb,
            grid->jlo_D3_fb, grid->jhi_D3_fb,
	    grid->klo_D3_fb, grid->khi_D3_fb);	    	       	    	    	    
  }	    	    
}


Grid *readGridFromAsciiFile(char *file_name)
{
  Grid *grid;
  FILE *fp;
  float x_lo_float[3], x_hi_float[3];
  float x_lo_ghostbox_float[3], x_hi_ghostbox_float[3];
  float dx_float[3];
  char   *line[80];
  
  /* open file and allocate Grid */
  fp = fopen(file_name,"r");
  grid = allocateGrid();

  fscanf(fp, "Number of dimensions: %d\n", &(grid->num_dims));
  if (grid->num_dims == 2) 
  {
    fscanf(fp,
           "Geometric limits (without ghostcells): [%g,%g] x [%g,%g]\n",
           &(x_lo_float[0]), &(x_hi_float[0]),
           &(x_lo_float[1]), &(x_hi_float[1]));
    grid->x_lo[0] = x_lo_float[0]; grid->x_hi[0] = x_hi_float[0];
    grid->x_lo[1] = x_lo_float[1]; grid->x_hi[1] = x_hi_float[1];
    fscanf(fp,"Geometric limits (with ghostcells): [%g,%g] x [%g,%g]\n",
           &(x_lo_ghostbox_float[0]), &(x_hi_ghostbox_float[0]),
           &(x_lo_ghostbox_float[1]), &(x_hi_ghostbox_float[1]));
    grid->x_lo_ghostbox[0] = x_lo_ghostbox_float[0]; 
    grid->x_hi_ghostbox[0] = x_hi_ghostbox_float[0];
    grid->x_lo_ghostbox[1] = x_lo_ghostbox_float[1]; 
    grid->x_hi_ghostbox[1] = x_hi_ghostbox_float[1];
    fscanf(fp,"Grid size (without ghostcells): %d x %d\n", 
           &((grid->grid_dims)[0]), 
           &((grid->grid_dims)[1]));
    fscanf(fp,"Grid size (with ghostcells): %d x %d\n", 
           &((grid->grid_dims_ghostbox)[0]),
           &((grid->grid_dims_ghostbox)[1]));
    fscanf(fp,"Grid spacing: dx = %g, dy = %g \n",
           &(dx_float[0]), &(dx_float[1]));
    grid->dx[0] = dx_float[0];
    grid->dx[1] = dx_float[1];
  } 
  else 
  {
    fscanf(fp,
           "Geometric limits (without ghostcells): [%g,%g] x [%g,%g] x [%g,%g]\n",
           &(x_lo_float[0]), &(x_hi_float[0]),
           &(x_lo_float[1]), &(x_hi_float[1]),
           &(x_lo_float[2]), &(x_hi_float[2]));
    grid->x_lo[0] = x_lo_float[0]; grid->x_hi[0] = x_hi_float[0];
    grid->x_lo[1] = x_lo_float[1]; grid->x_hi[1] = x_hi_float[1];
    grid->x_lo[2] = x_lo_float[2]; grid->x_hi[2] = x_hi_float[2];
    fscanf(fp,"Geometric limits (with ghostcells): [%g,%g] x [%g,%g] x [%g,%g]\n",
           &(x_lo_ghostbox_float[0]), &(x_hi_ghostbox_float[0]),
           &(x_lo_ghostbox_float[1]), &(x_hi_ghostbox_float[1]),
           &(x_lo_ghostbox_float[2]), &(x_hi_ghostbox_float[2]));
    grid->x_lo_ghostbox[0] = x_lo_ghostbox_float[0]; 
    grid->x_hi_ghostbox[0] = x_hi_ghostbox_float[0];
    grid->x_lo_ghostbox[1] = x_lo_ghostbox_float[1]; 
    grid->x_hi_ghostbox[1] = x_hi_ghostbox_float[1];
    grid->x_lo_ghostbox[2] = x_lo_ghostbox_float[2]; 
    grid->x_hi_ghostbox[2] = x_hi_ghostbox_float[2];
    fscanf(fp,"Grid size (without ghostcells): %d x %d x %d\n", 
           &((grid->grid_dims)[0]), 
           &((grid->grid_dims)[1]),
           &((grid->grid_dims)[2]));
    fscanf(fp,"Grid size (with ghostcells): %d x %d x %d\n", 
           &((grid->grid_dims_ghostbox)[0]),
           &((grid->grid_dims_ghostbox)[1]),
           &((grid->grid_dims_ghostbox)[2]));
    fscanf(fp,"Grid spacing: dx = %g, dy = %g, dz = %g\n",
           &(dx_float[0]), &(dx_float[1]), &(dx_float[2]));
    grid->dx[0] = dx_float[0];
    grid->dx[1] = dx_float[1];
    grid->dx[2] = dx_float[2];
  }
  fscanf(fp, "Total number of grid points: %d\n", &(grid->num_gridpts));

  /* fscanf(fp,"%s\n",line); */
  
  if (grid->num_dims == 2) 
  {
    fscanf(fp,
            "Ghost box: [%d,%d] x [%d,%d]\n",
            &(grid->ilo_gb), &(grid->ihi_gb),
            &(grid->jlo_gb), &(grid->jhi_gb));
	    
    fscanf(fp,
            "Fill box: [%d,%d] x [%d,%d]\n",
            &(grid->ilo_fb), &(grid->ihi_fb),
            &(grid->jlo_fb), &(grid->jhi_fb));
	    
    fscanf(fp,
            "1st order differences fill box: [%d,%d] x [%d,%d]\n",
            &(grid->ilo_D1_fb), &(grid->ihi_D1_fb),
            &(grid->jlo_D1_fb), &(grid->jhi_D1_fb));
    fscanf(fp,
            "2nd order differences fill box: [%d,%d] x [%d,%d]\n",
            &(grid->ilo_D2_fb), &(grid->ihi_D2_fb),
            &(grid->jlo_D2_fb), &(grid->jhi_D2_fb));
    fscanf(fp,
            "3rd order differences fill box: [%d,%d] x [%d,%d]\n",
            &(grid->ilo_D3_fb), &(grid->ihi_D3_fb),
            &(grid->jlo_D3_fb), &(grid->jhi_D3_fb));	 	    	     	    	    	    
  }
  else
  {
    fscanf(fp,
            "Ghost box: [%d,%d] x [%d,%d] x [%d,%d]\n",
            &(grid->ilo_gb), &(grid->ihi_gb),
            &(grid->jlo_gb), &(grid->jhi_gb),
	    &(grid->klo_gb), &(grid->khi_gb));
    fscanf(fp,
            "Fill box: [%d,%d] x [%d,%d] x [%d,%d]\n",
            &(grid->ilo_fb), &(grid->ihi_fb),
            &(grid->jlo_fb), &(grid->jhi_fb),
	    &(grid->klo_fb), &(grid->khi_fb));
	    	    
    fscanf(fp,
            "1st order differences fill box: [%d,%d] x [%d,%d] x [%d,%d]\n",
            &(grid->ilo_D1_fb), &(grid->ihi_D1_fb),
            &(grid->jlo_D1_fb), &(grid->jhi_D1_fb),
	    &(grid->klo_D1_fb), &(grid->khi_D1_fb));
	    
    fscanf(fp,
            "2nd order differences fill box: [%d,%d] x [%d,%d] x [%d,%d]\n",
            &(grid->ilo_D2_fb), &(grid->ihi_D2_fb),
            &(grid->jlo_D2_fb), &(grid->jhi_D2_fb),
	    &(grid->klo_D2_fb), &(grid->khi_D2_fb));
	    
    fscanf(fp,
            "3rd order differences fill box: [%d,%d] x [%d,%d] x [%d,%d]\n",
            &(grid->ilo_D3_fb), &(grid->ihi_D3_fb),
            &(grid->jlo_D3_fb), &(grid->jhi_D3_fb),
	    &(grid->klo_D3_fb), &(grid->khi_D3_fb));	    	       	    	    	    
  }	    	    

  fclose(fp);

  return grid;
}

void writeGridToAsciiFile(Grid *grid, char *file_name)
{
  FILE *fp = fopen(file_name,"w");

  printGrid(grid, fp);
 
  fclose(fp);
}


void writeGridToBinaryFile(Grid *grid, char *file_name)
{
  FILE *fp;
    
  fp = fopen(file_name,"w");

  fwrite(&(grid->num_dims), sizeof(int), 1, fp);
  fwrite(grid->x_lo, sizeof(double), 3, fp);
  fwrite(grid->x_hi, sizeof(double), 3, fp);
  fwrite(grid->x_lo_ghostbox, sizeof(double), 3, fp);
  fwrite(grid->x_hi_ghostbox, sizeof(double), 3, fp);
  fwrite(grid->grid_dims, sizeof(int), 3, fp); 
  fwrite(grid->grid_dims_ghostbox, sizeof(int), 3, fp); 
  fwrite(grid->dx, sizeof(double), 3, fp);
  fwrite(&(grid->num_gridpts), sizeof(int), 1, fp);

  fwrite(&(grid->ilo_gb), sizeof(int), 1, fp);
  fwrite(&(grid->ihi_gb), sizeof(int), 1, fp);
  fwrite(&(grid->jlo_gb), sizeof(int), 1, fp);
  fwrite(&(grid->jhi_gb), sizeof(int), 1, fp);  
  fwrite(&(grid->klo_gb), sizeof(int), 1, fp);
  fwrite(&(grid->khi_gb), sizeof(int), 1, fp);
 
  fwrite(&(grid->ilo_fb), sizeof(int), 1, fp);
  fwrite(&(grid->ihi_fb), sizeof(int), 1, fp);
  fwrite(&(grid->jlo_fb), sizeof(int), 1, fp);
  fwrite(&(grid->jhi_fb), sizeof(int), 1, fp);  
  fwrite(&(grid->klo_fb), sizeof(int), 1, fp);
  fwrite(&(grid->khi_fb), sizeof(int), 1, fp);
  
  fwrite(&(grid->ilo_D1_fb), sizeof(int), 1, fp);
  fwrite(&(grid->ihi_D1_fb), sizeof(int), 1, fp);
  fwrite(&(grid->jlo_D1_fb), sizeof(int), 1, fp);
  fwrite(&(grid->jhi_D1_fb), sizeof(int), 1, fp);  
  fwrite(&(grid->klo_D1_fb), sizeof(int), 1, fp);
  fwrite(&(grid->khi_D1_fb), sizeof(int), 1, fp);
  
  fwrite(&(grid->ilo_D2_fb), sizeof(int), 1, fp);
  fwrite(&(grid->ihi_D2_fb), sizeof(int), 1, fp);
  fwrite(&(grid->jlo_D2_fb), sizeof(int), 1, fp);
  fwrite(&(grid->jhi_D2_fb), sizeof(int), 1, fp);  
  fwrite(&(grid->klo_D2_fb), sizeof(int), 1, fp);
  fwrite(&(grid->khi_D2_fb), sizeof(int), 1, fp);
  
  fwrite(&(grid->ilo_D3_fb), sizeof(int), 1, fp);
  fwrite(&(grid->ihi_D3_fb), sizeof(int), 1, fp);
  fwrite(&(grid->jlo_D3_fb), sizeof(int), 1, fp);
  fwrite(&(grid->jhi_D3_fb), sizeof(int), 1, fp);  
  fwrite(&(grid->klo_D3_fb), sizeof(int), 1, fp);
  fwrite(&(grid->khi_D3_fb), sizeof(int), 1, fp);
  
  fclose(fp);
}
 

Grid *readGridFromBinaryFile(char *file_name)
{
  FILE *fp;
  Grid *grid;
  
 /* open file and allocate Grid */
  fp = fopen(file_name,"r");
  grid = allocateGrid();

  fread(&(grid->num_dims), sizeof(int), 1, fp);
  fread(grid->x_lo, sizeof(double), 3, fp);
  fread(grid->x_hi, sizeof(double), 3, fp);
  fread(grid->x_lo_ghostbox, sizeof(double), 3, fp);
  fread(grid->x_hi_ghostbox, sizeof(double), 3, fp);
  fread(grid->grid_dims, sizeof(int), 3, fp); 
  fread(grid->grid_dims_ghostbox, sizeof(int), 3, fp); 
  fread(grid->dx, sizeof(double), 3, fp);
  fread(&(grid->num_gridpts), sizeof(int), 1, fp);
  
  fread(&(grid->ilo_gb), sizeof(int), 1, fp);
  fread(&(grid->ihi_gb), sizeof(int), 1, fp);
  fread(&(grid->jlo_gb), sizeof(int), 1, fp);
  fread(&(grid->jhi_gb), sizeof(int), 1, fp);  
  fread(&(grid->klo_gb), sizeof(int), 1, fp);
  fread(&(grid->khi_gb), sizeof(int), 1, fp);
 
  fread(&(grid->ilo_fb), sizeof(int), 1, fp);
  fread(&(grid->ihi_fb), sizeof(int), 1, fp);
  fread(&(grid->jlo_fb), sizeof(int), 1, fp);
  fread(&(grid->jhi_fb), sizeof(int), 1, fp);  
  fread(&(grid->klo_fb), sizeof(int), 1, fp);
  fread(&(grid->khi_fb), sizeof(int), 1, fp);
  
  fread(&(grid->ilo_D1_fb), sizeof(int), 1, fp);
  fread(&(grid->ihi_D1_fb), sizeof(int), 1, fp);
  fread(&(grid->jlo_D1_fb), sizeof(int), 1, fp);
  fread(&(grid->jhi_D1_fb), sizeof(int), 1, fp);  
  fread(&(grid->klo_D1_fb), sizeof(int), 1, fp);
  fread(&(grid->khi_D1_fb), sizeof(int), 1, fp);
  
  fread(&(grid->ilo_D2_fb), sizeof(int), 1, fp);
  fread(&(grid->ihi_D2_fb), sizeof(int), 1, fp);
  fread(&(grid->jlo_D2_fb), sizeof(int), 1, fp);
  fread(&(grid->jhi_D2_fb), sizeof(int), 1, fp);  
  fread(&(grid->klo_D2_fb), sizeof(int), 1, fp);
  fread(&(grid->khi_D2_fb), sizeof(int), 1, fp);
  
  fread(&(grid->ilo_D3_fb), sizeof(int), 1, fp);
  fread(&(grid->ihi_D3_fb), sizeof(int), 1, fp);
  fread(&(grid->jlo_D3_fb), sizeof(int), 1, fp);
  fread(&(grid->jhi_D3_fb), sizeof(int), 1, fp);  
  fread(&(grid->klo_D3_fb), sizeof(int), 1, fp);
  fread(&(grid->khi_D3_fb), sizeof(int), 1, fp);

  fclose(fp);
  
  return grid;
}


void setIndexSpaceLimits(LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE accuracy, 
  Grid *grid)
{  
   switch( accuracy )
   {
     case LOW: { /* low accuracy */
       setIndexSpaceLimitsENO1(grid);
       break;
     }

     case MEDIUM: { /* medium accuracy */
       setIndexSpaceLimitsENO2(grid);
       break;
     }

     case HIGH: { /* high accuracy */
       setIndexSpaceLimitsENO3(grid);
       break;
     }

     case VERY_HIGH: { /* very high accuracy */
       setIndexSpaceLimitsWENO5(grid);
       break; 
     }

     default: { /* medium accuracy */
       setIndexSpaceLimitsENO2(grid);
     }
  }
}


