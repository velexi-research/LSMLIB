/*
 * File:        lsm_initialization3d.h
 * Copyright:   (c) 2005-2006 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.2 $
 * Modified:    $Date: 2006/05/25 19:41:22 $
 * Description: Header file for 3D initialization functions
 */

#ifndef included_lsm_initialization3d_h
#define included_lsm_initialization3d_h

#ifdef __cplusplus
extern "C" {
#endif


#include "lsm_grid.h"

/*! \file lsm_initialization3d.h
 * 
 * @ref lsm_initialization3d.h provides support for creating level set 
 * functions which possess zero level sets corresponding to common 
 * geometrical objects in three -dimensions.  These can be used as initial 
 * conditions for level set method calculations or masks for computations 
 * in restricted domains.  Whenever possible, the level set function is 
 * set to be an approximate signed distance function.
 *
 */   


/*!
 * createPlane() sets phi to be a level set function corresponding to 
 * a plane cutting through the computational domain.  The plane is assumed
 * to be represented in the form:
 *
 *     normal_x[i] * (x - point_x[i]) + normal_y[i] * (y - point_y[i])
 *   + normal_z[i] * (z - point_z[i]) = 0
 * 
 * Arguments:
 *  - phi (out):      level set function 
 *  - normal_x (in):  x-coordinate of normal vector to the plane
 *  - normal_y (in):  y-coordinate of normal vector to the plane
 *  - normal_z (in):  z-coordinate of normal vector to the plane
 *  - point_x (in):   x-coordinate of point that lies on the plane
 *  - point_y (in):   y-coordinate of point that lies on the plane
 *  - point_z (in):   z-coordinate of point that lies on the plane
 *  - grid (in):      pointer to Grid data structure 
 *
 * Return value:      none
 *
 * NOTES: 
 * - Is it the user's responsbility to ensure that memory for phi
 *   has been allocated.
 *  
 * - Normals are assumed to be ``outside'' normals; that is, the normal
 *   vector should point from the region where phi is negative to 
 *   the region where phi is positive.  
 *
 */
void createPlane(
  double *phi, 
  double normal_x, double normal_y, double normal_z,
  double point_x, double point_y, double point_z,
  Grid *grid);


/*!
 * createIntersectionOfHalfSpaces3d() sets phi to be a level set function 
 * corresponding to intersection of the specified 3D half-spaces.  The
 * i-th half-space is assumed to be represented in the form:
 *
 *     normal_x[i] * (x - point_x[i]) + normal_y[i] * (y - point_y[i])
 *   + normal_z[i] * (z - point_z[i]) < 0
 * 
 * Arguments:
 *  - phi (out):             level set function 
 *  - num_half_spaces (in):  number of half-spaces
 *  - normal_x (in):         array containing the x-coordinates for vectors 
 *                           normal to the planes defining the 3D half-spaces
 *  - normal_y (in):         array containing the y-coordinates for vectors 
 *                           normal to the planes defining the 3D half-spaces
 *  - normal_z (in):         array containing the z-coordinates for vectors 
 *                           normal to the planes defining the 3D half-spaces
 *  - point_x (in):          array containing the x-coordinates of points that 
 *                           lie on the planes defining the 3D half-spaces
 *  - point_y (in):          array containing the y-coordinates of points that 
 *                           lie on the planes defining the 3D half-spaces
 *  - point_z (in):          array containing the z-coordinates of points that 
 *                           lie on the planes defining the 3D half-spaces
 *  - grid (in):             pointer to Grid data structure 
 *
 * Return value:             none
 *
 * NOTES: 
 * - phi is approximately a signed distance function.  Within the region
 *   phi < 0, it is equal to the signed distance function;  in the region 
 *   phi > 0, it is a signed distance function everywhere except for 
 *   regions where phi > 0 simultaneously for multiple half-spaces.
 *
 * - Is it the user's responsbility to ensure that memory for phi
 *   has been allocated.
 *  
 * - Normals are assumed to be ``outside'' normals; that is, the normal
 *   vector should point from the region where phi is negative to 
 *   the region where phi is positive.  
 *
 */
void createIntersectionOfHalfSpaces3d(
  double *phi, int num_half_spaces,
  double *normal_x, double *normal_y, double *normal_z,
  double *point_x, double *point_y, double *point_z,
  Grid *grid);


/*!
 * createPolyhedron3d() sets phi to be a level set function corresponding to 
 * a polyhedron (possibly unbounded) with num_faces faces given by the 
 * intersection of num_faces half-spaces.  The i-th half-space is assumed
 * to be represented in the form:
 *
 *     normal_x[i] * (x - point_x[i]) + normal_y[i] * (y - point_y[i])
 *   + normal_z[i] * (z - point_z[i]) < 0
 * 
 * Arguments:
 *  - phi (out):       level set function
 *  - num_faces (in):  number of faces of polyhedron
 *  - normal_x (in):   array containing the x-coordinates for vectors normal
 *                     to the planes defining the faces of the polyhedron
 *  - normal_y (in):   array containing the y-coordinates for vectors normal
 *                     to the planes defining the faces of the polyhedron
 *  - normal_z (in):   array containing the z-coordinates for vectors normal
 *                     to the planes defining the faces of the polyhedron
 *  - point_x (in):    array containing the x-coordinates of points that lie
 *                     on the planes defining the faces of the polyhedron
 *  - point_y (in):    array containing the y-coordinates of points that lie
 *                     on the planes defining the faces of the polyhedron
 *  - point_z (in):    array containing the z-coordinates of points that lie
 *                     on the planes defining the faces of the polyhedron
 *  - grid (in):       pointer to Grid data structure 
 *
 * Return value:       none
 *
 * NOTES: 
 * - phi is approximately a signed distance function.  Within the 
 *   polyhedron (phi < 0), it is equal to the signed distance function.
 *   Outside of the polyhedron (phi > 0), it is a signed distance function
 *   everywhere except for regions where phi > 0 simultaneously for multiple
 *   half-spaces.
 *
 * - Is it the user's responsbility to ensure that memory for phi
 *   has been allocated.
 *  
 * - Normals are assumed to be ``outside'' normals; that is, the normal
 *   vector should point from the region where phi is negative to 
 *   the region where phi is positive.  
 *
 */
void createPolyhedron3d(
  double *phi, int num_faces,
  double *normal_x, double *normal_y, double *normal_z,
  double *point_x, double *point_y, double *point_z,
  Grid *grid);


/*!
 * createSphere() sets phi to be a level set function corresponding to 
 * a sphere.
 *
 * Arguments:
 *  - phi (out):         level set function 
 *  - center_x (in):     x-coordinate of the center of the sphere 
 *  - center_y (in):     y-coordinate of the center of the sphere 
 *  - center_z (in):     z-coordinate of the center of the sphere 
 *  - radius (in):       the radius of the sphere
 *  - inside_flag (in):  flag indicating whether the inside or outside of 
 *                       each sphere should be the region associated with 
 *                       negative values of the level set function.  If 
 *                       inside_flag is negative, then phi on the inside of 
 *                       the sphere is negative.  The reverse is true if 
 *                       inside_flag is nonnegative.
 *  - grid (in):         pointer to Grid data structure
 *
 * Return value:         none
 *
 * NOTES:
 * - Is it the user's responsbility to ensure that memory for phi
 *   has been allocated.
 *
 */
void createSphere(
  double *phi, 
  double center_x, double center_y, double center_z,
  double radius, 
  int inside_flag,
  Grid *grid);


/*!
 * createIntersectionOfSpheres() sets phi to be a level set function 
 * corresponding to the intersection of num_spheres spheres.
 *
 * Arguments:
 *  - phi (out):         level set function 
 *  - num_spheres (in):  number of spheres
 *  - center_x (in):     array containing the x-coordinates of the centers
 *                       of the spheres 
 *  - center_y (in):     array containing the y-coordinates of the centers
 *                       of the spheres 
 *  - center_z (in):     array containing the z-coordinates of the centers
 *                       of the spheres 
 *  - radius (in):       array containing the radii of the spheres
 *  - inside_flag (in):  array containing the flags indicating whether the
 *                       inside or outside of each sphere should be the
 *                       region associated with negative values of the
 *                       level set function.  If inside_flag[i] is negative,
 *                       then phi on the inside of the i-th sphere is
 *                       negative.  The reverse is true if inside_flag[i]
 *                       is nonnegative.
 *  - grid (in):         pointer to Grid data structure
 *
 * Return value:         none
 *
 * NOTES:
 * - phi is approximately a signed distance function.  Within the region
 *   phi < 0, it is equal to the signed distance function; in the region
 *   phi > 0, it is a signed distance function everywhere except for 
 *   regions where phi > 0 simultaneously for multiple spheres.
 *
 * - Is it the user's responsbility to ensure that memory for phi
 *   has been allocated.
 *
 */
void createIntersectionOfSpheres(
  double *phi, int num_spheres,
  double *center_x, double *center_y, double *center_z,
  double *radius, 
  int *inside_flag,
  Grid *grid);


/*!
 * createCylinder() sets phi to be a level set function corresponding to 
 * a cylinders with arbitrary axes.
 *
 * Arguments:
 *  - phi (out):         level set function 
 *  - tangent_x(in):     x-coordinate for the vector that define the 
 *                       direction of the axes of the cylinder
 *  - tangent_y(in):     y-coordinate for the vector that define the 
 *                       direction of the axes of the cylinder
 *  - tangent_z(in):     z-coordinate for the vector that define the 
 *                       direction of the axes of the cylinder
 *  - point_x (in):      x-coordinate of a point that lies on the axes 
 *                       of the cylinder
 *  - point_y (in):      y-coordinate of a point that lies on the axes 
 *                       of the cylinder
 *  - point_z (in):      z-coordinate of a point that lies on the axes 
 *                       of the cylinder
 *  - radius (in):       radius of the cylinder
 *  - inside_flag (in):  flag indicating whether the inside or outside 
 *                       the cylinder should be the region associated 
 *                       with negative values of the level set function.  
 *                       If inside_flag is negative, then phi on the 
 *                       inside of the cylinder is negative.  The reverse 
 *                       is true if inside_flag is nonnegative.
 *  - grid (in):         pointer to Grid data structure
 *
 * Return value:         none
 *
 */
void createCylinder(
  double *phi, 
  double tangent_x, double tangent_y, double tangent_z,
  double point_x, double point_y, double point_z,
  double radius, 
  int inside_flag,
  Grid *grid);


/*!
 * createIntersectionOfCylinders() sets phi to be a level set function 
 * corresponding to the intersection of num_cylinders cylinders with
 * arbitrary axes.
 *
 * Arguments:
 *  - phi (out):           level set function 
 *  - num_cylinders (in):  number of cylinders
 *  - tangent_x(in):       array containing the x-coordinates for vectors
 *                         that define the direction of the axes of the 
 *                         cylinders
 *  - tangent_y(in):       array containing the y-coordinates for vectors
 *                         that define the direction of the axes of the 
 *                         cylinders
 *  - tangent_z(in):       array containing the z-coordinates for vectors
 *                         that define the direction of the axes of the 
 *                         cylinders
 *  - point_x (in):        array containing the x-coordinates of points that
 *                         lie on the axes of the cylinders
 *  - point_y (in):        array containing the y-coordinates of points that
 *                         lie on the axes of the cylinders
 *  - point_z (in):        array containing the z-coordinates of points that
 *                         lie on the axes of the cylinders
 *  - radius (in):         array containing the radii of the cylinders
 *  - inside_flag (in):    array containing the flags indicating whether the
 *                         inside or outside of each cylinder should be the
 *                         region associated with negative values of the
 *                         level set function.  If inside_flag[i] is negative,
 *                         then phi on the inside of the i-th cylinder is
 *                         negative.  The reverse is true if inside_flag[i]
 *                         is nonnegative.
 *  - grid (in):           pointer to Grid data structure
 *
 * Return value:           none
 *
 */
void createIntersectionOfCylinders(
  double *phi, int num_cylinders,
  double *tangent_x, double *tangent_y, double *tangent_z,
  double *point_x, double *point_y, double *point_z,
  double *radius, 
  int *inside_flag,
  Grid *grid);


/*!
 * createHyperboloid() sets phi to be a level set function corresponding 
 * to a one-sheet hyperboloid with arbitrary axes.  In a coordinate frame where s is 
 * the coordinate along the axis measured from the center of the hyperboloid 
 * and t is the coordinate perpendicular to the axis measured from the axis, 
 * the hyperboloid satisfies the equation:
 *
 *   - s^2/alpha^2 + t^2/beta^2 = 1
 *
 * Arguments:
 *  - phi (out):         level set function 
 *  - tangent_x(in):     x-coordinate for the vector that defines the direction 
 *                       of the axes of the hyperboloid
 *  - tangent_y(in):     y-coordinate for the vector that defines the direction 
 *                       of the axes of the hyperboloid
 *  - tangent_z(in):     z-coordinate for the vector that defines the direction 
 *                       of the axes of the hyperboloid
 *  - center_x (in):     x-coordinate of the centers of the hyperboloid
 *  - center_y (in):     y-coordinate of the centers of the hyperboloid
 *  - center_z (in):     z-coordinate of the centers of the hyperboloid
 *  - alpha (in):        coefficient 'alpha' of the hyperboloid
 *  - beta (in):         coefficient 'beta' of the hyperboloid
 *  - inside_flag (in):  flag indicating whether phi should be taken to be 
 *                       (s^2/alpha^2 - t^2/beta^2 - 1) or
 *                       (-s^2/alpha^2 + t^2/beta^2 + 1).  If 
 *                       inside_flag is negative, the former value is 
 *                       used; otherwise the latter value is used.
 *  - grid (in):         pointer to Grid data structure
 *
 * Return value:         none
 *
 */
void createHyperboloid(
  double *phi, 
  double tangent_x, double tangent_y, double tangent_z,
  double center_x, double center_y, double center_z,
  double alpha, double beta,
  int inside_flag,
  Grid *grid);


/*!
 * createIntersectionOfHyperboloids() sets phi to be a level set function 
 * corresponding to the intersection of num_hyperboloids one sheet hyperboloids
 * with arbitrary axes.  In a coordinate frame where s is the coordinate along
 * the axis measured from the center of the hyperboloid and t is the 
 * coordinate perpendicular to the axis measured from the axis, each 
 * hyperboloid satisfies the equation:
 *
 *  - s^2/alpha^2 + t^2/beta^2 = 1
 *
 * Arguments:
 *  - phi (out):              level set function 
 *  - num_hyperboloids (in):  number of hyperboloids
 *  - tangent_x(in):          array containing the x-coordinates for vectors
 *                            that define the direction of the axes of the 
 *                            hyperboloids
 *  - tangent_y(in):          array containing the y-coordinates for vectors
 *                            that define the direction of the axes of the 
 *                            hyperboloids
 *  - tangent_z(in):          array containing the z-coordinates for vectors
 *                            that define the direction of the axes of the 
 *                            hyperboloids
 *  - center_x (in):          array containing the x-coordinates of the 
 *                            centers of the hyperboloids
 *  - center_y (in):          array containing the y-coordinates of the 
 *                            centers of the hyperboloids
 *  - center_z (in):          array containing the z-coordinates of the 
 *                            centers of the hyperboloids
 *  - alpha (in):             array containing the coefficient 'alpha' of 
 *                            the hyperboloids
 *  - beta (in):              array containing the coefficient 'beta' of
 *                            the hyperboloids
 *  - inside_flag (in):       array containing the flags indicating whether 
 *                            phi should be taken to be 
 *                            (s^2/alpha^2 - t^2/beta^2 - 1) or
 *                            (-s^2/alpha^2 + t^2/beta^2 + 1).  If 
 *                            inside_flag[i] is negative, the former value is 
 *                            used; otherwise the latter value is used.
 *  - grid (in):              pointer to Grid data structure
 *
 * Return value:              none
 *
 */
void createIntersectionOfHyperboloids(
  double *phi, int num_hyperboloids,
  double *tangent_x, double *tangent_y, double *tangent_z,
  double *center_x, double *center_y, double *center_z,
  double *alpha, double *beta,
  int *inside_flag,
  Grid *grid);


/*!
 * createCone() sets phi to be a level set function corresponding to a
 * cone with arbitrary axes.  In a coordinate frame where s is the coordinate 
 * along the axis measured from the center of the cone and t is the coordinate 
 * perpendicular to the axis measured from the axis, the cone satisfies 
 * the equation:
 *
 *   s^2/alpha^2 - t^2/beta^2 = 0
 *
 * Arguments:
 *  - phi (out):         level set function 
 *  - tangent_x(in):     x-coordinate for the vector that defines the 
 *                       direction of the axes of the cone
 *  - tangent_y(in):     y-coordinate for the vector that defines the 
 *                       direction of the axes of the cone
 *  - tangent_z(in):     z-coordinate for the vector that defines the 
 *                       direction of the axes of the cone
 *  - center_x (in):     x-coordinate of the center of the cone
 *  - center_y (in):     y-coordinate of the center of the cone
 *  - center_z (in):     z-coordinate of the center of the cone
 *  - alpha (in):        coefficient 'alpha' of the cone
 *  - beta (in):         coefficient 'beta' of the cone
 *  - inside_flag (in):  flag indicating whether phi should be taken to be 
 *                       (s^2/alpha^2 - t^2/beta^2) or
 *                       (-s^2/alpha^2 + t^2/beta^2).  If 
 *                       inside_flag is negative, the former value is 
 *                       used; otherwise the latter value is used.
 *  - grid (in):         pointer to Grid data structure
 *
 * Return value:         none
 *
 */
void createCone(
  double *phi, 
  double tangent_x, double tangent_y, double tangent_z,
  double center_x, double center_y, double center_z,
  double alpha, double beta,
  int inside_flag,
  Grid *grid);


/*!
 * createIntersectionOfCones() sets phi to be a level set function 
 * corresponding to the intersection of num_cones cones with
 * arbitrary axes.  In a coordinate frame where s is the coordinate along
 * the axis measured from the center of the cone and t is the coordinate 
 * perpendicular to the axis measured from the axis, each cone satisfies 
 * the equation:
 *
 *   s^2/alpha^2 - t^2/beta^2 = 0
 *
 * Arguments:
 *  - phi (out):         level set function 
 *  - num_cones (in):    number of cones
 *  - tangent_x(in):     array containing the x-coordinates for vectors
 *                       that define the direction of the axes of the cones
 *  - tangent_y(in):     array containing the y-coordinates for vectors
 *                       that define the direction of the axes of the cones
 *  - tangent_z(in):     array containing the z-coordinates for vectors
 *                       that define the direction of the axes of the cones
 *  - center_x (in):     array containing the x-coordinates of the 
 *                       centers of the cones
 *  - center_y (in):     array containing the y-coordinates of the 
 *                       centers of the cones
 *  - center_z (in):     array containing the z-coordinates of the 
 *                       centers of the cones
 *  - alpha (in):        array containing the coefficient 'alpha' of 
 *                       the cones
 *  - beta (in):         array containing the coefficient 'beta' of
 *                       the cones
 *  - inside_flag (in):  array containing the flags indicating whether 
 *                       phi should be taken to be 
 *                       (s^2/alpha^2 - t^2/beta^2) or
 *                       (-s^2/alpha^2 + t^2/beta^2).  If 
 *                       inside_flag[i] is negative, the former value is 
 *                       used; otherwise the latter value is used.
 *  - grid (in):         pointer to Grid data structure
 *
 * Return value:         none
 *
 */
void createIntersectionOfCones(
  double *phi, int num_cones,
  double *tangent_x, double *tangent_y, double *tangent_z,
  double *center_x, double *center_y, double *center_z,
  double *alpha, double *beta,
  int *inside_flag,
  Grid *grid);


/*!
 * createCuboid() sets phi to be a level set function corresponding to 
 * a cuboid with its rectangular faces parallel to the coordinate axes. 
 * The cuboid is represented by coordinates of its lower corner and 
 * lengths of its sides.
 * 
 * Arguments:
 *  - phi (out):           level set function 
 *  - corner_x (in):       x-coordinate for lower corner point
 *  - corner_y (in):       y-coordinate for lower corner point
 *  - corner_z (in):       z-coordinate for lower corner point
 *  - side_length_x (in):  side length of cuboid in the x-coordinate direction 
 *  - side_length_y (in):  side length of cuboid in the y-coordinate direction
 *  - side_length_z (in):  side length of cuboid in the z-coordinate direction
 *  - inside_flag (in):    flag indicating whether the inside or outside of 
 *                         cuboid should be the region associated with 
 *                         negative values of the level set function.  If 
 *                         inside_flag is negative, then phi on the inside of 
 *                         the cuboid is negative.  The reverse is true if 
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
 *   that define the cuboid.
 *
 * - Is it the user's responsbility to ensure that memory for phi
 *   has been allocated. 
 *
 */
void createCuboid(
  double *phi,
  double corner_x, double corner_y, double corner_z,
  double side_length_x, double side_length_y, double side_length_z,
  int inside_flag,
  Grid *grid);


/*!
 * createIntersectionOfCuboids() sets phi to be a level set function 
 * corresponding to intersection of the num_cuboids cuboids with rectangular
 * faces parallel to coordinate axes.  Each cuboid is represented 
 * by coordinates of its lower corner and length of its sides.
 * 
 * Arguments:
 *  - phi (out):           level set function 
 *  - num_cuboids (in):    number of cuboids
 *  - corner_x (in):       array containing the x-coordinates for lower 
 *                         corner of each cuboid
 *  - corner_y (in):       array containing the y-coordinatesfor lower
 *                         corner of each cuboid
 *  - corner_z (in):       array containing the z-coordinatesfor lower
 *                         corner of each cuboid
 *  - side_length_x (in):  array containing side lengths of the cuboids 
 *                         in x-coordinate direction 
 *  - side_length_y (in):  array containing side lengths of the cuboids 
 *                         in y-coordinate direction 
 *  - side_length_z (in):  array containing side lengths of the cuboids 
 *                         in z-coordinate direction 
 *  - inside_flag (in):    flag indicating whether the inside or outside of 
 *                         each cuboid should be the region associated with 
 *                         negative values of the level set function.  If 
 *                         inside_flag[i] is negative, then phi on the inside 
 *                         of the cuboid is negative.  The reverse is true if 
 *                         inside_flag[i] is nonnegative.
 *  - grid (in):           pointer to Grid data structure 
 *
 * Return value:           none
 *
 * NOTES: 
 * - phi is approximately a signed distance function.  Within the region
 *   phi < 0, it is equal to the signed distance function;  in the region 
 *   phi > 0, it is a signed distance function everywhere except for 
 *   regions where phi > 0 simultaneously for any two of the half-spaces
 *   that define the cuboids.
 *
 * - Is it the user's responsbility to ensure that memory for phi
 *   has been allocated. 
 *
 */
void  createIntersectionOfCuboids(
  double *phi,
  int num_cuboids,   
  double *corner_x, double *corner_y, double *corner_z,
  double *side_length_x, double *side_length_y, double *side_length_z,
  int *inside_flag,
  Grid *grid);     


#ifdef __cplusplus
}
#endif

#endif
