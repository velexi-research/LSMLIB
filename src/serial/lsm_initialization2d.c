/*
 * File:        lsm_initialization2d.c
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.3 $
 * Modified:    $Date: 2006/06/02 02:36:52 $
 * Description: Implementation file for 2D initialization functions
 */

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "LSMLIB_config.h"
#include "lsm_initialization2d.h"
#include "lsm_macros.h"


void createLine(
  LSMLIB_REAL *phi,
  LSMLIB_REAL normal_x, LSMLIB_REAL normal_y,
  LSMLIB_REAL point_x, LSMLIB_REAL point_y,
  Grid *grid)
{
  /* use createIntersectionOfHalfSpaces2d() */
  createIntersectionOfHalfSpaces2d(
    phi, 1, &normal_x, &normal_y, &point_x, &point_y, grid);
}


void createIntersectionOfHalfSpaces2d(
  LSMLIB_REAL* phi, int num_half_spaces,
  LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y,
  LSMLIB_REAL *point_x, LSMLIB_REAL *point_y,
  Grid *grid)
{
  int    i,j,l;
  int    idx;
  LSMLIB_REAL x,y;
  LSMLIB_REAL dot_prod,norm;
  LSMLIB_REAL max;
  LSMLIB_REAL signed_dist_to_line;
  
  for (j = 0; j < grid->grid_dims_ghostbox[1]; j++)
  {
    for (i = 0; i < grid->grid_dims_ghostbox[0]; i++) 
    {
      idx = i+j*grid->grid_dims_ghostbox[0];
      x = (grid->x_lo_ghostbox)[0] + (grid->dx)[0]*i;
      y = (grid->x_lo_ghostbox)[1] + (grid->dx)[1]*j;

      max = -LSMLIB_REAL_MAX;

      for (l = 0; l < num_half_spaces; l++)
      { 
        dot_prod = (x - point_x[l])*normal_x[l] + (y - point_y[l])*normal_y[l];
        norm = normal_x[l]*normal_x[l] + normal_y[l]*normal_y[l];
        signed_dist_to_line = dot_prod/norm;
        if (signed_dist_to_line > max) max = signed_dist_to_line;
      }

      phi[idx] = max;  

    }
  }
 
}


void createPolyhedron2d(
  LSMLIB_REAL* phi, int num_sides,
  LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y,
  LSMLIB_REAL *point_x, LSMLIB_REAL *point_y,
  Grid *grid)
{
  /* use createIntersectionOfHalfSpaces2d() */
  createIntersectionOfHalfSpaces2d(
    phi, num_sides, normal_x, normal_y, point_x, point_y, grid);
}


void createCircle(
  LSMLIB_REAL* phi,
  LSMLIB_REAL center_x, LSMLIB_REAL center_y,
  LSMLIB_REAL radius,
  int inside_flag,
  Grid *grid)
{
  /* use createIntersectionOfCircles() */
  createIntersectionOfCircles(
    phi, 1, &center_x, &center_y, &radius, &inside_flag, grid);
}


void createIntersectionOfCircles(
  LSMLIB_REAL *phi, 
  int num_circles,
  LSMLIB_REAL *center_x, LSMLIB_REAL *center_y,
  LSMLIB_REAL *radius,
  int *inside_flag,
  Grid *grid)
{
  int    i,j,l;
  int    idx;
  LSMLIB_REAL x,y;
  LSMLIB_REAL max;
  LSMLIB_REAL signed_dist_to_circle;
    
    
  for (j = 0; j < (grid->grid_dims_ghostbox)[1]; j++)
  {
    for (i = 0; i < (grid->grid_dims_ghostbox)[0]; i++) 
    {
      idx = i+j*(grid->grid_dims_ghostbox)[0];
      x = (grid->x_lo_ghostbox)[0] + (grid->dx)[0]*i;
      y = (grid->x_lo_ghostbox)[1] + (grid->dx)[1]*j;
    
      max = -LSMLIB_REAL_MAX;

      for (l = 0; l < num_circles; l++)
      {
        signed_dist_to_circle = sqrt( (x-center_x[l])*(x-center_x[l])
                                     +(y-center_y[l])*(y-center_y[l]) ) 
                                - radius[l];
        if (inside_flag[l] >= 0) {
          signed_dist_to_circle = -signed_dist_to_circle; 
        }
        if (signed_dist_to_circle > max) max = signed_dist_to_circle;
      }
      phi[idx] = max;
    }
  }
}


void createRectangle(
  LSMLIB_REAL *phi,
  LSMLIB_REAL corner_x,
  LSMLIB_REAL corner_y,
  LSMLIB_REAL side_length_x,
  LSMLIB_REAL side_length_y,
  int inside_flag,
  Grid *grid)
{
  /* use createIntersectionOfRectangles() */
  createIntersectionOfRectangles(phi, 1,
                                 &corner_x, &corner_y,
                                 &side_length_x, &side_length_y,
                                 &inside_flag, grid);
}     


void createIntersectionOfRectangles(
  LSMLIB_REAL *phi,
  int num_rectangles,   
  LSMLIB_REAL *corner_x, LSMLIB_REAL  *corner_y,
  LSMLIB_REAL *side_length_x, LSMLIB_REAL *side_length_y,
  int *inside_flag,
  Grid *grid)
{   
  LSMLIB_REAL point_x[4], point_y[4];
  LSMLIB_REAL normal_x[4], normal_y[4];
  int i, l, num_planes;
     
  LSMLIB_REAL  *phi1 = 
    (LSMLIB_REAL *)malloc(grid->num_gridpts*sizeof(LSMLIB_REAL)); 
        
  for(l = 0; l < num_rectangles; l++) {
    /* Each rectangle is the intersection of 4 half spaces */
    for(i = 0; i < 2; i++) {
      point_x[i] = corner_x[l];
      point_y[i] = corner_y[l];
    
      if (i == 0) normal_x[i] = -1;
      else        normal_x[i] = 0;
    
      if (i == 1) normal_y[i] = -1;
      else        normal_y[i] = 0;
    }
       
    for(i = 2; i < 4; i++)
    {
      point_x[i] = corner_x[l] + side_length_x[l];
      point_y[i] = corner_y[l] + side_length_y[l];
    
      if (i == 2) normal_x[i] = 1;
      else        normal_x[i] = 0;
    
      if (i == 3) normal_y[i] = 1;
      else        normal_y[i] = 0;
    }
      
    num_planes = 4;
    createIntersectionOfHalfSpaces2d(phi1, num_planes,
                                     normal_x, normal_y,
                                     point_x, point_y,
                                     grid);
    if (inside_flag[l] >= 1) 
    { 
      NEGATE_DATA(phi1,grid)
    }
      
    if (l == 0 ) COPY_DATA(phi,phi1,grid)
    else         IMPOSE_MASK(phi,phi1,phi,grid)              
  }
    
  free(phi1);
}

