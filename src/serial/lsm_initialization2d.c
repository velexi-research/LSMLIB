/*
 * File:        lsm_initialization2d.c
 * Copyright:   (c) 2005-2006 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.3 $
 * Modified:    $Date: 2006/06/02 02:36:52 $
 * Description: Implementation file for 2D initialization functions
 */

#include <math.h>
#include <stdlib.h>
#include <float.h>

#include "lsm_initialization2d.h"
#include "lsm_macros.h"


void createLine(
  double *phi,
  double normal_x, double normal_y,
  double point_x, double point_y,
  Grid *grid)
{
  /* use createIntersectionOfHalfSpaces2d() */
  createIntersectionOfHalfSpaces2d(
    phi, 1, &normal_x, &normal_y, &point_x, &point_y, grid);
}


void createIntersectionOfHalfSpaces2d(
  double* phi, int num_half_spaces,
  double *normal_x, double *normal_y,
  double *point_x, double *point_y,
  Grid *grid)
{
  int    i,j,l;
  int    idx;
  double x,y;
  double dot_prod,norm;
  double max;
  double signed_dist_to_line;
  
  for (j = 0; j < grid->grid_dims_ghostbox[1]; j++)
  {
    for (i = 0; i < grid->grid_dims_ghostbox[0]; i++) 
    {
      idx = i+j*grid->grid_dims_ghostbox[0];
      x = (grid->x_lo_ghostbox)[0] + (grid->dx)[0]*i;
      y = (grid->x_lo_ghostbox)[1] + (grid->dx)[1]*j;

      max = -DBL_MAX;

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
  double* phi, int num_sides,
  double *normal_x, double *normal_y,
  double *point_x, double *point_y,
  Grid *grid)
{
  /* use createIntersectionOfHalfSpaces2d() */
  createIntersectionOfHalfSpaces2d(
    phi, num_sides, normal_x, normal_y, point_x, point_y, grid);
}


void createCircle(
  double* phi,
  double center_x, double center_y,
  double radius,
  int inside_flag,
  Grid *grid)
{
  /* use createIntersectionOfCircles() */
  createIntersectionOfCircles(
    phi, 1, &center_x, &center_y, &radius, &inside_flag, grid);
}


void createIntersectionOfCircles(
  double *phi, 
  int num_circles,
  double *center_x, double *center_y,
  double *radius,
  int *inside_flag,
  Grid *grid)
{
  int    i,j,l;
  int    idx;
  double x,y;
  double max;
  double signed_dist_to_circle;
    
    
  for (j = 0; j < (grid->grid_dims_ghostbox)[1]; j++)
  {
    for (i = 0; i < (grid->grid_dims_ghostbox)[0]; i++) 
    {
      idx = i+j*(grid->grid_dims_ghostbox)[0];
      x = (grid->x_lo_ghostbox)[0] + (grid->dx)[0]*i;
      y = (grid->x_lo_ghostbox)[1] + (grid->dx)[1]*j;
    
      max = -DBL_MAX;

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
  double *phi,
  double corner_x,
  double corner_y,
  double side_length_x,
  double side_length_y,
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
  double *phi,
  int num_rectangles,   
  double *corner_x, double  *corner_y,
  double *side_length_x, double *side_length_y,
  int *inside_flag,
  Grid *grid)
{   
  double point_x[6], point_y[6];
  double normal_x[6], normal_y[6];
  int i, l, num_planes;
     
  double  *phi1 = (double *)malloc(grid->num_gridpts*sizeof(double)); 
        
  for(l = 0; l < num_rectangles; l++) {
    /* Each rectangle is the intersection of 4 half spaces */
    for(i = 0; i < 3; i++) {
      point_x[i] = corner_x[l];
      point_y[i] = corner_y[l];
    
      if (i == 0) normal_x[i] = -1;
      else        normal_x[i] = 0;
    
      if (i == 1) normal_y[i] = -1;
      else        normal_y[i] = 0;
    }
       
    for(i = 3; i < 6; i++)
    {
      point_x[i] = corner_x[l] + side_length_x[l];
      point_y[i] = corner_y[l] + side_length_y[l];
    
      if (i == 3) normal_x[i] = 1;
      else        normal_x[i] = 0;
    
      if (i == 4) normal_y[i] = 1;
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
