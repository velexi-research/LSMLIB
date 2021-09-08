/*
 * File:        lsm_initialization2d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for 2D initialization functions
 */

#ifndef included_lsm_initialization2d_h
#define included_lsm_initialization2d_h

#include "lsmlib_config.h"

#ifdef __cplusplus
extern "C" {
#endif


#include "lsm_grid.h"

/*! \file lsm_initialization2d.h 
 *
 * @ref lsm_initialization2d.h provides support for creating level set 
 * functions which possess zero level sets corresponding to common 
 * geometrical objects in two-dimensions.  These can be used as initial 
 * conditions for level set method calculations or masks for computations 
 * in restricted domains.  Whenever possible, the level set function is 
 * set to be an approximate signed distance function.
 *
 */   


/*!
 * createLine() sets phi to be a level set function corresponding to a line
 * through the computational domain.  The line is assumed to be represented
 * in the form:
 *
 *   normal_x * (x - point_x) + normal_y * (y - point_y) = 0.
 *
 * Arguments
 *  - phi (out):      level set function
 *  - normal_x (in):  the x-coordinate for the normal vector
 *                    to the line
 *  - normal_y (in):  the y-coordinate for the normal vector
 *                    to the line
 *  - point_x (in):   the x-coordinate of a point that lies
 *                    on the line 
 *  - point_y (in):   the y-coordinate of a point that lies
 *                    on the line 
 *  - grid (in):      pointer to Grid data structure 
 *
 * Return value:      none
 *
 * NOTES: 
 * - Is it the user's responsibility to ensure that memory for phi
 *   has been allocated.
 *
 * - The normal vector is assumed to be an ``outside'' normal; that is, the 
 *   normal vector should point from the region where phi is negative to 
 *   the region where phi is positive. 
 *
 */
void createLine(
  LSMLIB_REAL *phi, 
  LSMLIB_REAL normal_x, LSMLIB_REAL normal_y,
  LSMLIB_REAL point_x, LSMLIB_REAL point_y,
  Grid *grid);


/*!
 * createIntersectionOfHalfSpaces2d() sets phi to be a level set function 
 * corresponding to the intersection of the specified 2D half-spaces.  The 
 * i-th half-space is assumed to be represented in the form:
 *
 *   normal_x[i] * (x - point_x[i]) + normal_y[i] * (y - point_y[i]) < 0.
 *
 * Arguments
 *  - phi (out):             level set function
 *  - num_half_spaces (in):  number of half-spaces
 *  - normal_x (in):         array containing the x-coordinates for vectors 
 *                           normal to the lines defining the half-spaces
 *  - normal_y (in):         array containing the y-coordinates for vectors 
 *                           normal to the lines defining the half-spaces
 *  - point_x (in):          array containing the x-coordinates of points that 
 *                           lie on the lines defining the half-spaces
 *  - point_y (in):          array containing the y-coordinates of points that 
 *                           lie on the lines defining the half-spaces
 *  - grid (in):             pointer to Grid data structure 
 *
 * Return value:             none
 *
 * NOTES: 
 * - phi is approximately a signed distance function.  Within the region
 *   phi < 0, it is equal to the signed distance function; in the region
 *   phi > 0, it is a signed distance function everywhere except for 
 *   regions where phi > 0 simultaneously for multiple half-spaces.
 *
 * - Is it the user's responsibility to ensure that memory for phi
 *   has been allocated.
 *
 * - Normals are assumed to be ``outside'' normals; that is, the normal
 *   vector should point from the region where phi is negative to 
 *   the region where phi is positive. 
 *
 */
void createIntersectionOfHalfSpaces2d(
  LSMLIB_REAL *phi, int num_half_spaces,
  LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y,
  LSMLIB_REAL *point_x, LSMLIB_REAL *point_y,
  Grid *grid);


/*!
 * createPolyhedron2d() sets phi to be a level set function corresponding to 
 * a 2D-polyhedron (possibly unbounded) with num_sides sides given by the 
 * intersection of num_sides half-spaces.  The i-th half-spaces is assumed 
 * to be represented in the form:
 *
 *   normal_x[i] * (x - point_x[i]) + normal_y[i] * (y - point_y[i]) < 0.
 *
 * Arguments
 *  - phi (out):       level set function
 *  - num_sides (in):  number of sides of polyhedron
 *  - normal_x (in):   array containing the x-coordinates for vectors normal 
 *                     to the lines defining the sides of the polyhedron
 *  - normal_y (in):   array containing the y-coordinates for vectors normal 
 *                     to the lines defining the sides of the polyhedron
 *  - point_x (in):    array containing the x-coordinates of points that lie 
 *                     on the lines defining the sides of the polyhedron
 *  - point_y (in):    array containing the y-coordinates of points that lie 
 *                     on the lines defining the sides of the polyhedron
 *  - grid (in):       pointer to Grid data structure 
 *
 * Return value:       none
 *
 * NOTES: 
 * - phi is approximately a signed distance function.  Within the
 *   polyhedron (phi < 0), it is equal to the signed distance function.  
 *   Outside of the polyhedron (phi > 0), it is a signed distance function 
 *   except for regions where phi > 0 simultaneously for multiple half-spaces.
 *
 * - Is it the user's responsibility to ensure that memory for phi
 *   has been allocated.
 *
 * - Normals are assumed to be ``outside'' normals; that is, the normal
 *   vector should point from the region where phi is negative to 
 *   the region where phi is positive. 
 *
 */
void createPolyhedron2d(
  LSMLIB_REAL *phi, int num_sides,
  LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y,
  LSMLIB_REAL *point_x, LSMLIB_REAL *point_y,
  Grid *grid);


/*!
 * createIntersectionOfPolyhedra2d() sets phi to be a level set function 
 * corresponding to an intersection of 2D-polyhedra, each possibly unbounded,
 * and with different number of sides (half spaces). The i-th half-space in
 * each polyhedron is assumed to be represented in the form:
 *
 *   normal_x[i] * (x - point_x[i]) + normal_y[i] * (y - point_y[i]) < 0.
 *  (if inside_flag[i] is negative, reversed if positive.)
 *
 * Arguments
 *  - phi (out):       level set function
 *  - num_polyhedra (in):  number of sides of polyhedron
 *  - idx_start(in):   start indices in below arrays for each polyhedron
 *  - idx_end(in):     end indices in below arrays for each polyhedrion
 *  - normal_x (in):   array containing the x-coordinates for vectors normal 
 *                     to the lines defining the sides of the polyhedron
 *                     l-th polyhedron is represented with
 *                     normal_x[i], where i=idx_start[l],...,idx_end[l]
 *                     Same is true for all other arrays.
 *  - normal_y (in):   array containing the y-coordinates for vectors normal 
 *                     to the lines defining the sides of the polyhedron
 *  - point_x (in):    array containing the x-coordinates of points that lie 
 *                     on the lines defining the sides of the polyhedron
 *  - point_y (in):    array containing the y-coordinates of points that lie 
 *                     on the lines defining the sides of the polyhedron
 *  - inside_flag(in): array containing the num_polyhedra flags indicating 
 *                     whether the inside or outside of each polyhedron should be the
 *                     region associated with negative values of the 
 *                     level set function.  If inside_flag[l] is negative,
 *                     then phi on the inside of the l-th polyhedron is 
 *                     negative.  The reverse is true if inside_flag[l] 
 *                     is nonnegative.
 *  - grid (in):       pointer to Grid data structure 
 *
 * Return value:       none
 *
 * NOTES: 
 * - phi is approximately a signed distance function.  Within the
 *   polyhedron (phi < 0), it is equal to the signed distance function.  
 *   Outside of the polyhedron (phi > 0), it is a signed distance function 
 *   except for regions where phi > 0 simultaneously for multiple half-spaces.
 *
 * - Is it the user's responsibility to ensure that memory for phi
 *   has been allocated.
 *
 */
void createIntersectionOfPolyhedra2d(
  LSMLIB_REAL *phi,
  int num_polyhedra,
  int *idx_start,
  int *idx_end, 
  LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y,
  LSMLIB_REAL *point_x, LSMLIB_REAL *point_y,
  int *inside_flag,
  Grid *grid);
 
  
/*!
 * createCircle() sets phi to be a level set function corresponding to a 
 * circle.
 *
 * Arguments:
 *  - phi (out):         level set function array
 *  - center_x (in):     x-coordinate of the center of the circle
 *  - center_y (in):     y-coordinate of the center of the circle
 *  - radius (in):       radius of the circle
 *  - inside_flag (in):  flag indicating whether the inside or outside 
 *                       of the circle should be the region associated 
 *                       with negative values of the level set function.  
 *                       If inside_flag is negative, then phi on the inside 
 *                       of the circle is negative.  The reverse is true if 
 *                       inside_flag is nonnegative.
 *  - grid(in):          pointer to Grid structure
 *
 * Return value:         none
 *
 * NOTES: 
 * - Is it the user's responsibility to ensure that memory for phi
 *   has been allocated.
 *
 */
void createCircle(
  LSMLIB_REAL *phi, 
  LSMLIB_REAL center_x, LSMLIB_REAL center_y,
  LSMLIB_REAL radius,
  int inside_flag,
  Grid *grid);


/*!
 * createIntersectionOfCircles() sets phi to be a level set function 
 * corresponding to the intersection of num_circles circles.  
 *
 * Arguments:
 *  - phi (out):         level set function array
 *  - num_circles (in):  number of circles 
 *  - center_x (in):     array containing the x-coordinates for centers of
 *                       the circles
 *  - center_y (in):     array containing the y-coordinates for centers of
 *                       the circles
 *  - radius (in):       array containing the radii of the circles
 *  - inside_flag (in):  array containing the flags indicating whether the
 *                       inside or outside of each circle should be the
 *                       region associated with negative values of the 
 *                       level set function.  If inside_flag[i] is negative,
 *                       then phi on the inside of the i-th circle is 
 *                       negative.  The reverse is true if inside_flag[i] 
 *                       is nonnegative.
 *  - grid(in):          pointer to Grid structure
 *
 * Return value:         none
 *
 * NOTES: 
 * - phi is approximately a signed distance function.  Within the region
 *   phi < 0, it is equal to the signed distance function.  In the region
 *   phi > 0, it is a signed distance function everywhere except for 
 *   regions where phi > 0 simultaneously for the level set functions of
 *   multiple circles.
 *
 * - Is it the user's responsibility to ensure that memory for phi
 *   has been allocated.
 *
 */
void createIntersectionOfCircles(
  LSMLIB_REAL *phi, int num_circles, 
  LSMLIB_REAL *center_x, LSMLIB_REAL *center_y,
  LSMLIB_REAL *radius,
  int *inside_flag,
  Grid *grid);


/*!
 * createRectangle() sets phi to be a level set function corresponding to 
 * a rectangle with its sides parallel to coordinate axes. The rectangle is 
 * represented by the coordinates of its lower left corner and its side 
 * lengths.
 *
 * Arguments:
 *  - phi (out):           level set function
 *  - corner_x (in):       x-coordinate for lower left corner 
 *  - corner_y (in):       y-coordinate for lower left corner
 *  - side_length_x (in):  side length of rectangle in the x-coordinate 
 *                         direction
 *  - side_length_y (in):  side length of rectangle in the y-coordinate 
 *                         direction
 *  - inside_flag (in):    flag indicating whether the inside or outside of
 *                         rectangle should be the region associated with
 *                         negative values of the level set function.  If
 *                         inside_flag is negative, then phi on the inside of
 *                         the rectangle is negative.  The reverse is true if
 *                         inside_flag is nonnegative.
 *  - grid (in):           pointer to Grid data structure
 *
 * Return value:           none
 *
 * NOTES:
 * - phi is approximately a signed distance function.  Within the region
 *   phi < 0, it is equal to the signed distance function;  in the region
 *   phi > 0, it is a signed distance function everywhere except for
 *   regions where phi > 0 simultaneously for any two of the half-spaces
 *   that define the rectangle.
 *
 * - Is it the user's responsbility to ensure that memory for phi
 *   has been allocated.
 *
 */
void createRectangle(
  LSMLIB_REAL *phi,
  LSMLIB_REAL corner_x, LSMLIB_REAL corner_y,
  LSMLIB_REAL side_length_x, LSMLIB_REAL side_length_y,
  int inside_flag,
  Grid *grid);


/*! 
 * createIntersectionOfRectangles() sets phi to be a level set function 
 * corresponding to intersection of num_rectangles rectangles with sides
 * parallel to the coordinate axes.  Each rectangle is represented by the 
 * coordinates of its lower left corner and its side lengths.
 * 
 * Arguments: 
 *  - phi (out):            level set function 
 *  - num_rectangles (in):  number of rectangles 
 *  - corner_x (in):        array containing the x-coordinates for lower 
 *                          left corner of each rectangle 
 *  - corner_y (in):        array containing the y-coordinates for lower 
 *                          left corner of each rectangle 
 *  - side_length_x (in):   array containing side lengths of the rectangles
 *                          in the x-coordinate direction 
 *  - side_length_y (in):   array containing side lengths of the rectangles
 *                          in the y-coordinate direction 
 *  - inside_flag (in):     flag indicating whether the inside or outside of 
 *                          each rectangle should be the region associated 
 *                          with negative values of the level set function.  
 *                          If inside_flag[i] is negative, then phi on the 
 *                          inside of the rectangle is negative.  The reverse 
 *                          is true if inside_flag[i] is nonnegative.
 *  - grid (in):            pointer to Grid data structure
 * 
 * Return value:            none
 * 
 * NOTES: 
 * - phi is approximately a signed distance function.  Within the region
 *   phi < 0, it is equal to the signed distance function;  in the region
 *   phi > 0, it is a signed distance function everywhere except for
 *   regions where phi > 0 simultaneously for any two of the half-spaces
 *   that define the rectangles.
 * 
 * - Is it the user's responsbility to ensure that memory for phi
 *   has been allocated.
 *
 */
void createIntersectionOfRectangles(
  LSMLIB_REAL *phi,
  int num_rectangles,   
  LSMLIB_REAL *corner_x, LSMLIB_REAL *corner_y,
  LSMLIB_REAL *side_length_x, LSMLIB_REAL *side_length_y,
  int *inside_flag,
  Grid *grid);

     
#ifdef __cplusplus
}
#endif

#endif 
