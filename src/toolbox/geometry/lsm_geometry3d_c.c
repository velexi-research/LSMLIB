/*
 * File:        lsm_geometry3d_c.c
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.10 $
 * Modified:    $Date: 2007/04/20 16:39:27 $
 * Description: Implementation of 3D C geometry functions for level set method 
 */

#include <math.h>
#include <float.h>
#include "LSMLIB_config.h"
#include "lsm_geometry3d.h"

/* MACROS */
#define  LSM_GEOM_3D_ABS(x)       ( ((x) > 0) ? (x) : -(x) )
#define  LSM_GEOM_3D_ZERO_TOL     (1.0e-12)

/*
 * LSM_GEOM_3D_CROSS() computes the cross-product of two vectors.
 *
 * Arguments:
 *   cross_* (out):  components of cross-product
 *   v1_* (in):      components of vector 1
 *   v2_* (in):      components of vector 2
 *
 * NOTES:
 *  (1) cross_* MUST be valid l-values.
 *
 */
#define  LSM_GEOM_3D_CROSS( cross_x, cross_y, cross_z,                   \
                                v1_x, v1_y, v1_z,                        \
                                v2_x, v2_y, v2_z)                        \
{                                                                        \
  (cross_x) = (v1_y)*(v2_z) - (v1_z)*(v2_y);                             \
  (cross_y) = (v1_z)*(v2_x) - (v1_x)*(v2_z);                             \
  (cross_z) = (v1_x)*(v2_y) - (v1_y)*(v2_x);                             \
}

/*
 * LSM_GEOM_3D_SAME_SIDE() determines if the two points are on the same
 * side of a line.
 *
 * Arguments:
 *   result (out):  1 if two points are on the same side; 0 otherwise
 *   p1_* (in):     components of point 1
 *   p2_* (in):     components of point 2
 *   p_ref_* (in):  components of a point on line 
 *   edge_* (in):   components of direction of line
 *
 * NOTES:
 *  (1) result MUST be a valid l-value.
 *
 */
#define  LSM_GEOM_3D_SAME_SIDE( result,                                  \
                                    p1_x, p1_y, p1_z,                    \
                                    p2_x, p2_y, p2_z,                    \
                                    p_ref_x, p_ref_y, p_ref_z,           \
                                    edge_x, edge_y, edge_z )             \
{                                                                        \
  LSMLIB_REAL p1_minus_p_ref_x = (p1_x) - (p_ref_x);                          \
  LSMLIB_REAL p1_minus_p_ref_y = (p1_y) - (p_ref_y);                          \
  LSMLIB_REAL p1_minus_p_ref_z = (p1_z) - (p_ref_z);                          \
  LSMLIB_REAL p2_minus_p_ref_x = (p2_x) - (p_ref_x);                          \
  LSMLIB_REAL p2_minus_p_ref_y = (p2_y) - (p_ref_y);                          \
  LSMLIB_REAL p2_minus_p_ref_z = (p2_z) - (p_ref_z);                          \
                                                                         \
  LSMLIB_REAL cross1_x, cross1_y, cross1_z;                                   \
  LSMLIB_REAL cross2_x, cross2_y, cross2_z;                                   \
  LSMLIB_REAL norm_p1_minus_p_ref_sq;                                         \
  LSMLIB_REAL norm_p2_minus_p_ref_sq;                                         \
                                                                         \
  LSMLIB_REAL edge_len_sq;                                                    \
                                                                         \
  norm_p1_minus_p_ref_sq = p1_minus_p_ref_x*p1_minus_p_ref_x             \
                         + p1_minus_p_ref_y*p1_minus_p_ref_y             \
                         + p1_minus_p_ref_z*p1_minus_p_ref_z;            \
                                                                         \
  norm_p2_minus_p_ref_sq = p2_minus_p_ref_x*p2_minus_p_ref_x             \
                         + p2_minus_p_ref_y*p2_minus_p_ref_y             \
                         + p2_minus_p_ref_z*p2_minus_p_ref_z;            \
                                                                         \
  edge_len_sq = edge_x*edge_x + edge_y*edge_y + edge_z*edge_z;           \
                                                                         \
  if (    (norm_p1_minus_p_ref_sq > LSM_GEOM_3D_ZERO_TOL)                \
       && (norm_p2_minus_p_ref_sq > LSM_GEOM_3D_ZERO_TOL) ) {            \
                                                                         \
    LSMLIB_REAL norm_cross_1_sq, norm_cross_2_sq;                             \
                                                                         \
    LSM_GEOM_3D_CROSS(cross1_x, cross1_y, cross1_z,                      \
      p1_minus_p_ref_x, p1_minus_p_ref_y, p1_minus_p_ref_z,              \
      (edge_x), (edge_y), (edge_z) );                                    \
                                                                         \
    LSM_GEOM_3D_CROSS(cross2_x, cross2_y, cross2_z,                      \
      p2_minus_p_ref_x, p2_minus_p_ref_y, p2_minus_p_ref_z,              \
      (edge_x), (edge_y), (edge_z) );                                    \
                                                                         \
    norm_cross_1_sq = cross1_x*cross1_x + cross1_y*cross1_y              \
                    + cross1_z*cross1_z;                                 \
    norm_cross_2_sq = cross2_x*cross2_x + cross2_y*cross2_y              \
                    + cross2_z*cross2_z;                                 \
                                                                         \
    if (  (norm_cross_1_sq >                                             \
           LSM_GEOM_3D_ZERO_TOL*edge_len_sq*norm_p1_minus_p_ref_sq)      \
       && (norm_cross_2_sq >                                             \
          LSM_GEOM_3D_ZERO_TOL*edge_len_sq*norm_p2_minus_p_ref_sq) ) {   \
                                                                         \
      if ( (cross1_x*cross2_x + cross1_y*cross2_y + cross1_z*cross2_z)   \
           < 0) {                                                        \
        (result) = 0;                                                    \
      } else {                                                           \
        (result) = 1;                                                    \
      }                                                                  \
    } else {                                                             \
        (result) = 1;                                                    \
    }                                                                    \
  } else {                                                               \
      (result) = 1;                                                      \
  }                                                                      \
}


/* LSM3D_findLineInTetrahedron() */
int LSM3D_findLineInTetrahedron(
  LSMLIB_REAL *endpt1,
  LSMLIB_REAL *endpt2,
  const LSMLIB_REAL *x1,
  const LSMLIB_REAL *x2,
  const LSMLIB_REAL *x3,
  const LSMLIB_REAL *x4,
  const LSMLIB_REAL *phi,
  const LSMLIB_REAL *psi)
{
  /* number of intersections of {phi = 0,psi = 0} line with tetrahedron */
  int count = 0;

  /* coefficients for linear approximation to phi and psi */
  LSMLIB_REAL alpha_0, alpha_1, alpha_2, alpha_3;
  LSMLIB_REAL beta_0, beta_1, beta_2, beta_3;

  /* matrix for computing linear approximation to phi and psi */
  LSMLIB_REAL X11, X12, X13, X14;
  LSMLIB_REAL X21, X22, X23, X24;
  LSMLIB_REAL X31, X32, X33, X34;
  LSMLIB_REAL X41, X42, X43, X44;

  /* RHS values for linear system used to solve for alpha and beta */
  LSMLIB_REAL rhs_alpha_1, rhs_alpha_2, rhs_alpha_3, rhs_alpha_4;
  LSMLIB_REAL rhs_beta_1, rhs_beta_2, rhs_beta_3, rhs_beta_4;

  /* Householder reflection variables */
  LSMLIB_REAL norm, inv_norm_sq;
  LSMLIB_REAL inner_product_with_V;  /* V is the Householder reflection vector */

  /* {phi = 0, psi = 0} line variables */
  LSMLIB_REAL x_on_line, y_on_line, z_on_line;
  LSMLIB_REAL line_dir_x, line_dir_y, line_dir_z;

  /* variables for computing intersection of {phi=0,psi=0} line with faces */
  LSMLIB_REAL vector_in_plane_1_x, vector_in_plane_1_y, vector_in_plane_1_z;
  LSMLIB_REAL vector_in_plane_2_x, vector_in_plane_2_y, vector_in_plane_2_z;
  LSMLIB_REAL vector_in_plane_3_x, vector_in_plane_3_y, vector_in_plane_3_z;
  LSMLIB_REAL normal_x, normal_y, normal_z, plane_constant;
  LSMLIB_REAL normal_dot_line;
  LSMLIB_REAL intersect_coef_max = -LSMLIB_REAL_MAX;
  LSMLIB_REAL intersect_coef_min =  LSMLIB_REAL_MAX;


  /* check that {phi = 0,psi = 0} line intersects the tetrahedron */
  if  ( ( (phi[0]*phi[1]>0) && (phi[0]*phi[2]>0) && (phi[0]*phi[3]>0) ) ||
        ( (psi[0]*psi[1]>0) && (psi[0]*psi[2]>0) && (psi[0]*psi[3]>0) ) ) {
    return 0; 
  }

  /* initialize endpt1 and endpt2 */
  endpt1[0] = 0.0; endpt1[1] = 0.0; endpt1[2] = 0.0; 
  endpt2[0] = 0.0; endpt2[1] = 0.0; endpt2[2] = 0.0; 


  /******************************************************************* 
   *
   * Compute the linear approximations to phi and psi within the 
   * tetrahedron using the formulae:
   *
   *   phi = alpha_0 + alpha_1 * x + alpha_2 * y + alpha_3 * z
   * 
   * and 
   * 
   *   psi = beta_0 + beta_1 * x + beta_2 * y + beta_3 * z
   *
   * where alpha and beta are vectors containing the coefficients of 
   * linear approximations, the i-th row in X is the vector [1 x_i] 
   * with x_i the coordinate of the i-th corner of the tetrahedron,
   * and phi[i] and psi[i] are the value of phi and psi at the i-th 
   * corner of the tetrahedron.
   *
   *******************************************************************/

  /* 
   * form matrix of coordinates of corners of tetrahedron used to compute
   * the coefficients of linear approximations to phi and psi within
   * the tetrahedron.
   */
  X11 = 1.0; X12 = x1[0]; X13 = x1[1]; X14 = x1[2]; 
  X21 = 1.0; X22 = x2[0]; X23 = x2[1]; X24 = x2[2]; 
  X31 = 1.0; X32 = x3[0]; X33 = x3[1]; X34 = x3[2]; 
  X41 = 1.0; X42 = x4[0]; X43 = x4[1]; X44 = x4[2]; 

  /* 
   * vector of RHS values for 
   * the coefficients of linear approximations to phi and psi within
   * the tetrahedron.
   */
  rhs_alpha_1 = phi[0]; rhs_alpha_2 = phi[1]; 
  rhs_alpha_3 = phi[2]; rhs_alpha_4 = phi[3]; 
  rhs_beta_1 = psi[0]; rhs_beta_2 = psi[1]; 
  rhs_beta_3 = psi[2]; rhs_beta_4 = psi[3]; 

  /* 
   * solve for alpha and beta using QR factorization via
   * Householder reflections
   *
   * NOTE:  To conserve memory and improve computational performance,
   *        the Householder reflection vectors, V, are computed using
   *        the same variables as the matrix of tetrahedron corner 
   *        coordinates.  Also, the lower-left triangle of the X 
   *        matrix is NOT updated.
   */

  /* 
   * first-column 
   *
   * NOTE:  since we know that the first /column of the matrix X
   *        consists only of 1's, we can carry out many calculations
   *        for the first Householder reflection by hand.
   */
  norm = 2.0;
  X11 = 3.0;
  inv_norm_sq = 1.0/12.0;
  inner_product_with_V = X11*X12 + X21*X22 + X31*X32 + X41*X42;
  X12 -= 2.0*X11*inner_product_with_V*inv_norm_sq;
  X22 -= 2.0*X21*inner_product_with_V*inv_norm_sq;
  X32 -= 2.0*X31*inner_product_with_V*inv_norm_sq;
  X42 -= 2.0*X41*inner_product_with_V*inv_norm_sq;
  inner_product_with_V = X11*X13 + X21*X23 + X31*X33 + X41*X43;
  X13 -= 2.0*X11*inner_product_with_V*inv_norm_sq;
  X23 -= 2.0*X21*inner_product_with_V*inv_norm_sq;
  X33 -= 2.0*X31*inner_product_with_V*inv_norm_sq;
  X43 -= 2.0*X41*inner_product_with_V*inv_norm_sq;
  inner_product_with_V = X11*X14 + X21*X24 + X31*X34 + X41*X44;
  X14 -= 2.0*X11*inner_product_with_V*inv_norm_sq;
  X24 -= 2.0*X21*inner_product_with_V*inv_norm_sq;
  X34 -= 2.0*X31*inner_product_with_V*inv_norm_sq;
  X44 -= 2.0*X41*inner_product_with_V*inv_norm_sq;

  inner_product_with_V = ( X11*rhs_alpha_1 + X21*rhs_alpha_2 + 
                           X31*rhs_alpha_3 + X41*rhs_alpha_4 );
  rhs_alpha_1 -= 2.0*X11*inner_product_with_V*inv_norm_sq;
  rhs_alpha_2 -= 2.0*X21*inner_product_with_V*inv_norm_sq;
  rhs_alpha_3 -= 2.0*X31*inner_product_with_V*inv_norm_sq;
  rhs_alpha_4 -= 2.0*X41*inner_product_with_V*inv_norm_sq;
  inner_product_with_V = ( X11*rhs_beta_1 + X21*rhs_beta_2 + 
                           X31*rhs_beta_3 + X41*rhs_beta_4 );
  rhs_beta_1 -= 2.0*X11*inner_product_with_V*inv_norm_sq;
  rhs_beta_2 -= 2.0*X21*inner_product_with_V*inv_norm_sq;
  rhs_beta_3 -= 2.0*X31*inner_product_with_V*inv_norm_sq;
  rhs_beta_4 -= 2.0*X41*inner_product_with_V*inv_norm_sq;

  /* update X11 */
  X11 = -2.0;


  /* 
   * second-column 
   */
  norm = sqrt(X22*X22 + X32*X32 + X42*X42);
  X22 = (X22 > 0) ? (X22+norm) : (X22-norm);
  inv_norm_sq = 1.0/(X22*X22 + X32*X32 + X42*X42);
  inner_product_with_V = (X22*X23 + X32*X33 + X42*X43);
  X23 -= 2.0*X22*inner_product_with_V*inv_norm_sq;
  X33 -= 2.0*X32*inner_product_with_V*inv_norm_sq;
  X43 -= 2.0*X42*inner_product_with_V*inv_norm_sq;
  inner_product_with_V = (X22*X24 + X32*X34 + X42*X44);
  X24 -= 2.0*X22*inner_product_with_V*inv_norm_sq;
  X34 -= 2.0*X32*inner_product_with_V*inv_norm_sq;
  X44 -= 2.0*X42*inner_product_with_V*inv_norm_sq;

  inner_product_with_V = ( X22*rhs_alpha_2 + X32*rhs_alpha_3 
                         + X42*rhs_alpha_4 );
  rhs_alpha_2 -= 2.0*X22*inner_product_with_V*inv_norm_sq;
  rhs_alpha_3 -= 2.0*X32*inner_product_with_V*inv_norm_sq;
  rhs_alpha_4 -= 2.0*X42*inner_product_with_V*inv_norm_sq;
  inner_product_with_V = ( X22*rhs_beta_2 + X32*rhs_beta_3 + X42*rhs_beta_4 );
  rhs_beta_2 -= 2.0*X22*inner_product_with_V*inv_norm_sq;
  rhs_beta_3 -= 2.0*X32*inner_product_with_V*inv_norm_sq;
  rhs_beta_4 -= 2.0*X42*inner_product_with_V*inv_norm_sq;
 
  /* update X22 */
  X22 = (X22 > 0) ? -norm : norm;  


  /* 
   * third-column 
   */
  norm = sqrt(X33*X33 + X43*X43);
  X33 = (X33 > 0) ? (X33+norm) : (X33-norm);
  inv_norm_sq = 1.0/(X33*X33 + X43*X43);
  inner_product_with_V = (X33*X34 + X43*X44);
  X34 -= 2.0*X33*inner_product_with_V*inv_norm_sq;
  X44 -= 2.0*X43*inner_product_with_V*inv_norm_sq;

  inner_product_with_V = ( X33*rhs_alpha_3 + X43*rhs_alpha_4 );
  rhs_alpha_3 -= 2.0*X33*inner_product_with_V*inv_norm_sq;
  rhs_alpha_4 -= 2.0*X43*inner_product_with_V*inv_norm_sq;
  inner_product_with_V = ( X33*rhs_beta_3 + X43*rhs_beta_4 );
  rhs_beta_3 -= 2.0*X33*inner_product_with_V*inv_norm_sq;
  rhs_beta_4 -= 2.0*X43*inner_product_with_V*inv_norm_sq;

  /* update X33 */
  X33 = (X33 > 0) ? -norm : norm;  


  /* backsolve for alpha and beta */
  alpha_3 = rhs_alpha_4/X44;
  alpha_2 = (rhs_alpha_3 - X34*alpha_3)/X33;
  alpha_1 = (rhs_alpha_2 - X23*alpha_2 - X24*alpha_3)/X22;
  alpha_0 = (rhs_alpha_1 - X12*alpha_1 - X13*alpha_2 - X14*alpha_3)/X11;

  beta_3 = rhs_beta_4/X44;
  beta_2 = (rhs_beta_3 - X34*beta_3)/X33;
  beta_1 = (rhs_beta_2 - X23*beta_2 - X24*beta_3)/X22;
  beta_0 = (rhs_beta_1 - X12*beta_1 - X13*beta_2 - X14*beta_3)/X11;

/* KTC - DEBUGGING */
/*
  printf("X matrix (upper triangle part is Q in QR decomposition): \n");
  printf("%f %f %f %f\n", X11, X12, X13, X14);  
  printf("%f %f %f %f\n", X21, X22, X23, X24);  
  printf("%f %f %f %f\n", X31, X32, X33, X34);  
  printf("%f %f %f %f\n", X41, X42, X43, X44);  
  printf("alpha: = %f %f %f %f\n", alpha_0, alpha_1, alpha_2, alpha_3);
  printf("beta: = %f %f %f %f\n", beta_0, beta_1, beta_2, beta_3);
*/


  /******************************************************************* 
   *
   * Find a point and the direction of the {phi = 0,psi = 0} line (i.e
   * intersection of the planes {phi = 0} and {psi = 0}).
   *
   * To find a point on {phi = 0, psi = 0}, we find a solution for
   * the equations:
   * 
   * alpha_0 + alpha_1*x_on_line + alpha_2*y_on_line + alpha_3*z_on_line = 0
   * beta_0 + beta_1*x_on_line  + beta_2*y_on_line  + beta_3*z_on_line = 0
   *
   * The direction of the line {phi = 0, psi = 0} is computed by taking
   * the cross product of the gradients of phi and psi:
   * 
   * line_dir = grad(phi) cross grad(psi)
   *
   *******************************************************************/

  line_dir_x = alpha_2*beta_3 - alpha_3*beta_2;
  line_dir_y = alpha_3*beta_1 - alpha_1*beta_3;
  line_dir_z = alpha_1*beta_2 - alpha_2*beta_1;

/* KTC - DEBUGGING */
/*
  printf("line_dir : = %f %f %f\n", line_dir_x, line_dir_y, line_dir_z);
*/

  if (    (LSM_GEOM_3D_ABS(alpha_1) < LSM_GEOM_3D_ZERO_TOL) 
       && (LSM_GEOM_3D_ABS(beta_1) < LSM_GEOM_3D_ZERO_TOL) ) { 
    /* case: if line exists, it is parallel to x-axis     */

    if ( LSM_GEOM_3D_ABS(line_dir_x) < LSM_GEOM_3D_ZERO_TOL ) { 

      /* {phi = 0} and {psi = 0} planes do not intersect!   */
      /* this shouldn't happen!  return error!              */
      return -1;

    } else {  /* case: line exists and is parallel to x-axis */

      LSMLIB_REAL abs_alpha_2 = LSM_GEOM_3D_ABS(alpha_2);
      LSMLIB_REAL abs_beta_2 = LSM_GEOM_3D_ABS(beta_2);

      x_on_line = 0.0;  /* arbitrarily take x_on_line to be 0.0 */

      /* solve for y_on_line and z_on_line using Gaussian  */
      /* elimination with partial pivoting                 */
      if (abs_alpha_2 > abs_beta_2) {

        LSMLIB_REAL elimination_factor = beta_2/alpha_2;

        /* beta_2 = 0.0;   no need to actually set beta_2 since it  */
        /*                 is not used.                             */
        beta_3 -= elimination_factor*alpha_3;
        beta_0 -= elimination_factor*alpha_0;

        z_on_line = -beta_0/beta_3;
        y_on_line = -(alpha_0+alpha_3*z_on_line)/alpha_2;

      } else {       

        LSMLIB_REAL elimination_factor = alpha_2/beta_2;

        /* alpha_2 = 0.0;   no need to actually set alpha_2 since it  */
        /*                  is not used.                              */
        alpha_3 -= elimination_factor*beta_3;
        alpha_0 -= elimination_factor*beta_0;

        z_on_line = -alpha_0/alpha_3;
        y_on_line = -(beta_0+beta_3*z_on_line)/beta_2;
        
      }
    }

  } else {     /* case: line not parallel to x-axis */

    LSMLIB_REAL abs_alpha_1 = LSM_GEOM_3D_ABS(alpha_1);
    LSMLIB_REAL abs_beta_1 = LSM_GEOM_3D_ABS(beta_1);

    /* solve for a point on the {phi = 0, psi = 0} line and its     */
    /* direction using Gaussian elimination with partial pivoting   */
    if (abs_alpha_1 > abs_beta_1) {   /* case: abs(alpha_1) > abs(beta_1) */

      LSMLIB_REAL elimination_factor = beta_1/alpha_1;

      /* beta_1 = 0.0;   no need to actually set beta_1 since it  */
      /*                 is not used.                             */
      beta_2 -= elimination_factor*alpha_2;
      beta_3 -= elimination_factor*alpha_3;
      beta_0 -= elimination_factor*alpha_0;

      if ( LSM_GEOM_3D_ABS(beta_2) < LSM_GEOM_3D_ZERO_TOL ) { 
        /* case: line perpendicular to z-axis */

        if ( LSM_GEOM_3D_ABS(beta_3) < LSM_GEOM_3D_ZERO_TOL ) { 

          /* {phi = 0} and {psi = 0} planes do not intersect!   */
          /* this shouldn't happen!  return error!              */
          return -2;

        } else {  /* case: line lies in xy-plane with z = -beta_0/beta_3 */

          z_on_line = -beta_0/beta_3;
          y_on_line = 0.0;  /* arbitrarily take y_on_line to be 0.0 */
          x_on_line = -(alpha_0+alpha_3*z_on_line)/alpha_1;

        } 

      } else { /* case: line not parallel to x-axis or normal to z-axis */

          z_on_line = 0.0;  /* arbitrarily take z_on_line to be 0.0 */
          y_on_line = -beta_0/beta_2;
          x_on_line = -(alpha_0+alpha_2*y_on_line)/alpha_1;
      } 

    } else {   /* case: abs(alpha_1) <= abs(beta_1) */

      LSMLIB_REAL elimination_factor = alpha_1/beta_1;

      /* alpha_1 = 0.0;   no need to actually set alpha_1 since it  */
      /*                  is not used.                              */
      alpha_2 -= elimination_factor*beta_2;
      alpha_3 -= elimination_factor*beta_3;
      alpha_0 -= elimination_factor*beta_0;

      if ( LSM_GEOM_3D_ABS(alpha_2) < LSM_GEOM_3D_ZERO_TOL ) { 
        /* case: line perpendicular to z-axis */

        if ( LSM_GEOM_3D_ABS(alpha_3) < LSM_GEOM_3D_ZERO_TOL ) { 

          /* {phi = 0} and {psi = 0} planes do not intersect!   */
          /* this shouldn't happen!  return error!              */
          return -3;

        } else {  /* case: line lies in xy-plane with z = -alpha_0/alpha_3 */

          z_on_line = -alpha_0/alpha_3;
          y_on_line = 0.0;  /* arbitrarily take y_on_line to be 0.0 */
          x_on_line = -(beta_0+beta_3*z_on_line)/beta_1;

        } 

      } else { /* case: line not parallel to x-axis or normal to z-axis */

          z_on_line = 0.0;  /* arbitrarily take z_on_line to be 0.0 */
          y_on_line = -alpha_0/alpha_2;
          x_on_line = -(beta_0+beta_2*y_on_line)/beta_1;

      }
 
    }  
  }  /* end solve for {phi=0,psi=0} line */

/* KTC - DEBUGGING */
/*
  printf("Point on line: %f %f %f\n", x_on_line, y_on_line, z_on_line);
  printf("Line direction: %f %f %f\n", line_dir_x, line_dir_y, line_dir_z);
*/


  /******************************************************************* 
   *
   * Find the endpoints of the {phi = 0, psi = 0} line within the 
   * tetrahedron by checking for intersections with each of the 
   * tetrahedron's four faces.
   * 
   *******************************************************************/

  /* 
   * compute intersection with face defined by x1, x2, x3
   */

  /* compute equation of plane containing the face defined by x1, x2, x3 */
  vector_in_plane_1_x = x2[0] - x1[0];
  vector_in_plane_1_y = x2[1] - x1[1];
  vector_in_plane_1_z = x2[2] - x1[2];
  vector_in_plane_2_x = x3[0] - x1[0];
  vector_in_plane_2_y = x3[1] - x1[1];
  vector_in_plane_2_z = x3[2] - x1[2];
  LSM_GEOM_3D_CROSS(normal_x, normal_y, normal_z,
    vector_in_plane_1_x, vector_in_plane_1_y, vector_in_plane_1_z,
    vector_in_plane_2_x, vector_in_plane_2_y, vector_in_plane_2_z);

  plane_constant = -(normal_x*x1[0] + normal_y*x1[1] + normal_z*x1[2]);

  /* compute intersection of {phi = 0, psi = 0} line with current face */
  normal_dot_line = normal_x*line_dir_x + normal_y*line_dir_y 
                  + normal_z*line_dir_z;

  if ( LSM_GEOM_3D_ABS(normal_dot_line) > LSM_GEOM_3D_ZERO_TOL ) { 
    /* case: single intersection exists */
    
    LSMLIB_REAL intersect_coef = -(plane_constant + normal_x*x_on_line 
      + normal_y*y_on_line + normal_z*z_on_line)/normal_dot_line;

    LSMLIB_REAL x_intersect = x_on_line + intersect_coef*line_dir_x;
    LSMLIB_REAL y_intersect = y_on_line + intersect_coef*line_dir_y;
    LSMLIB_REAL z_intersect = z_on_line + intersect_coef*line_dir_z;


    /* check if point is in the interior of face using point in    */
    /* triangle test.                                              */

    /* check edge defined by x1 and x2 */
    int on_same_side_1 = 0;
    LSM_GEOM_3D_SAME_SIDE(on_same_side_1,
      x_intersect, y_intersect, z_intersect,
      x3[0], x3[1], x3[2],
      x1[0], x1[1], x1[2],
      vector_in_plane_1_x, vector_in_plane_1_y, vector_in_plane_1_z);

    if (on_same_side_1) {

      /* check edge defined by x1 and x3 */
      int on_same_side_2 = 0;
      LSM_GEOM_3D_SAME_SIDE(on_same_side_2,
        x_intersect, y_intersect, z_intersect,
        x2[0], x2[1], x2[2],
        x1[0], x1[1], x1[2],
        vector_in_plane_2_x, vector_in_plane_2_y, vector_in_plane_2_z);

      if (on_same_side_2) {

        /* check edge defined by x2 and x3 */
        int on_same_side_3 = 0;
        vector_in_plane_3_x = x3[0] - x2[0];
        vector_in_plane_3_y = x3[1] - x2[1];
        vector_in_plane_3_z = x3[2] - x2[2];
        LSM_GEOM_3D_SAME_SIDE(on_same_side_3,
          x_intersect, y_intersect, z_intersect,
          x1[0], x1[1], x1[2],
          x2[0], x2[1], x2[2],
          vector_in_plane_3_x, vector_in_plane_3_y, vector_in_plane_3_z);

        if (on_same_side_3) {  /* point lies within face of tetrahedron */
          
          if (intersect_coef_max < intersect_coef) {
            intersect_coef_max = intersect_coef;
            endpt2[0] = x_intersect; 
            endpt2[1] = y_intersect; 
            endpt2[2] = z_intersect; 
          }

          if (intersect_coef_min > intersect_coef) {
            intersect_coef_min = intersect_coef;
            endpt1[0] = x_intersect; 
            endpt1[1] = y_intersect; 
            endpt1[2] = z_intersect; 
          }

          count++;  /* increment intersection counter */

/* KTC - DEBUGGING */
/*
  printf("(x1,x2,x3) face\n");
  printf("intersect_coef = %f\n", intersect_coef);
  printf("endpt1 = (%f,%f,%f)\n", endpt1[0], endpt1[1], endpt1[2]);
  printf("endpt2 = (%f,%f,%f)\n", endpt2[0], endpt2[1], endpt2[2]);
*/

        }
      }
    } /* end point in triangle test */
  } /* end case: {phi=0, psi=0} line intersects face defined by x1, x2, x3 */


  /* 
   * compute intersection with face defined by x1, x2, x4
   */

  /* compute equation of plane containing the face defined by x1, x2, x4 */
  vector_in_plane_1_x = x2[0] - x1[0];
  vector_in_plane_1_y = x2[1] - x1[1];
  vector_in_plane_1_z = x2[2] - x1[2];
  vector_in_plane_2_x = x4[0] - x1[0];
  vector_in_plane_2_y = x4[1] - x1[1];
  vector_in_plane_2_z = x4[2] - x1[2];
  LSM_GEOM_3D_CROSS(normal_x, normal_y, normal_z,
    vector_in_plane_1_x, vector_in_plane_1_y, vector_in_plane_1_z,
    vector_in_plane_2_x, vector_in_plane_2_y, vector_in_plane_2_z);

  plane_constant = -(normal_x*x1[0] + normal_y*x1[1] + normal_z*x1[2]);

  /* compute intersection of {phi = 0, psi = 0} line with current face */
  normal_dot_line = normal_x*line_dir_x + normal_y*line_dir_y 
                  + normal_z*line_dir_z;

  if ( LSM_GEOM_3D_ABS(normal_dot_line) > LSM_GEOM_3D_ZERO_TOL ) { 
    /* case: single intersection exists */
    
    LSMLIB_REAL intersect_coef = -(plane_constant + normal_x*x_on_line 
      + normal_y*y_on_line + normal_z*z_on_line)/normal_dot_line;

    LSMLIB_REAL x_intersect = x_on_line + intersect_coef*line_dir_x;
    LSMLIB_REAL y_intersect = y_on_line + intersect_coef*line_dir_y;
    LSMLIB_REAL z_intersect = z_on_line + intersect_coef*line_dir_z;


    /* check if point is in the interior of face using point in    */
    /* triangle test.                                              */

    /* check edge defined by x1 and x2 */
    int on_same_side_1 = 0;
    LSM_GEOM_3D_SAME_SIDE(on_same_side_1,
      x_intersect, y_intersect, z_intersect,
      x4[0], x4[1], x4[2],
      x1[0], x1[1], x1[2],
      vector_in_plane_1_x, vector_in_plane_1_y, vector_in_plane_1_z);

    if (on_same_side_1) {

      /* check edge defined by x1 and x4 */
      int on_same_side_2 = 0;
      LSM_GEOM_3D_SAME_SIDE(on_same_side_2,
        x_intersect, y_intersect, z_intersect,
        x2[0], x2[1], x2[2],
        x1[0], x1[1], x1[2],
        vector_in_plane_2_x, vector_in_plane_2_y, vector_in_plane_2_z);

      if (on_same_side_2) {

        /* check edge defined by x2 and x4 */
        int on_same_side_3 = 0;
        vector_in_plane_3_x = x4[0] - x2[0];
        vector_in_plane_3_y = x4[1] - x2[1];
        vector_in_plane_3_z = x4[2] - x2[2];
        LSM_GEOM_3D_SAME_SIDE(on_same_side_3,
          x_intersect, y_intersect, z_intersect,
          x1[0], x1[1], x1[2],
          x2[0], x2[1], x2[2],
          vector_in_plane_3_x, vector_in_plane_3_y, vector_in_plane_3_z);

        if (on_same_side_3) {  /* point lies within face of tetrahedron */
          
          if (intersect_coef_max < intersect_coef) {
            intersect_coef_max = intersect_coef;
            endpt2[0] = x_intersect; 
            endpt2[1] = y_intersect; 
            endpt2[2] = z_intersect; 
          }

          if (intersect_coef_min > intersect_coef) {
            intersect_coef_min = intersect_coef;
            endpt1[0] = x_intersect; 
            endpt1[1] = y_intersect; 
            endpt1[2] = z_intersect; 
          }

          count++;  /* increment intersection counter */

/* KTC - DEBUGGING */
/*
  printf("(x1,x2,x4) face\n");
  printf("intersect_coef = %f\n", intersect_coef);
  printf("endpt1 = (%f,%f,%f)\n", endpt1[0], endpt1[1], endpt1[2]);
  printf("endpt2 = (%f,%f,%f)\n", endpt2[0], endpt2[1], endpt2[2]);
*/

        }
      }
    } /* end point in triangle test */
  } /* end case: {phi=0, psi=0} line intersects face defined by x1, x2, x4 */


  /* 
   * compute intersection with face defined by x1, x3, x4
   */

  /* compute equation of plane containing the face defined by x1, x3, x4 */
  vector_in_plane_1_x = x3[0] - x1[0];
  vector_in_plane_1_y = x3[1] - x1[1];
  vector_in_plane_1_z = x3[2] - x1[2];
  vector_in_plane_2_x = x4[0] - x1[0];
  vector_in_plane_2_y = x4[1] - x1[1];
  vector_in_plane_2_z = x4[2] - x1[2];
  LSM_GEOM_3D_CROSS(normal_x, normal_y, normal_z,
    vector_in_plane_1_x, vector_in_plane_1_y, vector_in_plane_1_z,
    vector_in_plane_2_x, vector_in_plane_2_y, vector_in_plane_2_z);

  plane_constant = -(normal_x*x1[0] + normal_y*x1[1] + normal_z*x1[2]);

  /* compute intersection of {phi = 0, psi = 0} line with current face */
  normal_dot_line = normal_x*line_dir_x + normal_y*line_dir_y 
                  + normal_z*line_dir_z;

  if ( LSM_GEOM_3D_ABS(normal_dot_line) > LSM_GEOM_3D_ZERO_TOL ) { 
    /* case: single intersection exists */
    
    LSMLIB_REAL intersect_coef = -(plane_constant + normal_x*x_on_line 
      + normal_y*y_on_line + normal_z*z_on_line)/normal_dot_line;

    LSMLIB_REAL x_intersect = x_on_line + intersect_coef*line_dir_x;
    LSMLIB_REAL y_intersect = y_on_line + intersect_coef*line_dir_y;
    LSMLIB_REAL z_intersect = z_on_line + intersect_coef*line_dir_z;


    /* check if point is in the interior of face using point in    */
    /* triangle test.                                              */

    /* check edge defined by x1 and x3 */
    int on_same_side_1 = 0;
    LSM_GEOM_3D_SAME_SIDE(on_same_side_1,
      x_intersect, y_intersect, z_intersect,
      x4[0], x4[1], x4[2],
      x1[0], x1[1], x1[2],
      vector_in_plane_1_x, vector_in_plane_1_y, vector_in_plane_1_z);

    if (on_same_side_1) {

      /* check edge defined by x1 and x4 */
      int on_same_side_2 = 0;
      LSM_GEOM_3D_SAME_SIDE(on_same_side_2,
        x_intersect, y_intersect, z_intersect,
        x3[0], x3[1], x3[2],
        x1[0], x1[1], x1[2],
        vector_in_plane_2_x, vector_in_plane_2_y, vector_in_plane_2_z);

      if (on_same_side_2) {

        /* check edge defined by x3 and x4 */
        int on_same_side_3 = 0;
        vector_in_plane_3_x = x4[0] - x3[0];
        vector_in_plane_3_y = x4[1] - x3[1];
        vector_in_plane_3_z = x4[2] - x3[2];
        LSM_GEOM_3D_SAME_SIDE(on_same_side_3,
          x_intersect, y_intersect, z_intersect,
          x1[0], x1[1], x1[2],
          x3[0], x3[1], x3[2],
          vector_in_plane_3_x, vector_in_plane_3_y, vector_in_plane_3_z);

        if (on_same_side_3) {  /* point lies within face of tetrahedron */
          
          if (intersect_coef_max < intersect_coef) {
            intersect_coef_max = intersect_coef;
            endpt2[0] = x_intersect; 
            endpt2[1] = y_intersect; 
            endpt2[2] = z_intersect; 
          }

          if (intersect_coef_min > intersect_coef) {
            intersect_coef_min = intersect_coef;
            endpt1[0] = x_intersect; 
            endpt1[1] = y_intersect; 
            endpt1[2] = z_intersect; 
          }

          count++;  /* increment intersection counter */

/* KTC - DEBUGGING */
/*
  printf("(x1,x3,x4) face\n");
  printf("intersect_coef = %f\n", intersect_coef);
  printf("endpt1 = (%f,%f,%f)\n", endpt1[0], endpt1[1], endpt1[2]);
  printf("endpt2 = (%f,%f,%f)\n", endpt2[0], endpt2[1], endpt2[2]);
*/

        }
      }
    } /* end point in triangle test */
  } /* end case: {phi=0, psi=0} line intersects face defined by x1, x3, x4 */


  /* 
   * compute intersection with face defined by x2, x3, x4
   */

  /* compute equation of plane containing the face defined by x2, x3, x4 */
  vector_in_plane_1_x = x3[0] - x2[0];
  vector_in_plane_1_y = x3[1] - x2[1];
  vector_in_plane_1_z = x3[2] - x2[2];
  vector_in_plane_2_x = x4[0] - x2[0];
  vector_in_plane_2_y = x4[1] - x2[1];
  vector_in_plane_2_z = x4[2] - x2[2];
  LSM_GEOM_3D_CROSS(normal_x, normal_y, normal_z,
    vector_in_plane_1_x, vector_in_plane_1_y, vector_in_plane_1_z,
    vector_in_plane_2_x, vector_in_plane_2_y, vector_in_plane_2_z);

  plane_constant = -(normal_x*x2[0] + normal_y*x2[1] + normal_z*x2[2]);

  /* compute intersection of {phi = 0, psi = 0} line with current face */
  normal_dot_line = normal_x*line_dir_x + normal_y*line_dir_y 
                  + normal_z*line_dir_z;

  if ( LSM_GEOM_3D_ABS(normal_dot_line) > LSM_GEOM_3D_ZERO_TOL ) { 
    /* case: single intersection exists */
    
    LSMLIB_REAL intersect_coef = -(plane_constant + normal_x*x_on_line 
      + normal_y*y_on_line + normal_z*z_on_line)/normal_dot_line;

    LSMLIB_REAL x_intersect = x_on_line + intersect_coef*line_dir_x;
    LSMLIB_REAL y_intersect = y_on_line + intersect_coef*line_dir_y;
    LSMLIB_REAL z_intersect = z_on_line + intersect_coef*line_dir_z;


    /* check if point is in the interior of face using point in    */
    /* triangle test.                                              */

    /* check edge defined by x2 and x3 */
    int on_same_side_1 = 0;
    LSM_GEOM_3D_SAME_SIDE(on_same_side_1,
      x_intersect, y_intersect, z_intersect,
      x4[0], x4[1], x4[2],
      x2[0], x2[1], x2[2],
      vector_in_plane_1_x, vector_in_plane_1_y, vector_in_plane_1_z);

    if (on_same_side_1) {

      /* check edge defined by x2 and x4 */
      int on_same_side_2 = 0;
      LSM_GEOM_3D_SAME_SIDE(on_same_side_2,
        x_intersect, y_intersect, z_intersect,
        x3[0], x3[1], x3[2],
        x2[0], x2[1], x2[2],
        vector_in_plane_2_x, vector_in_plane_2_y, vector_in_plane_2_z);

      if (on_same_side_2) {

        /* check edge defined by x3 and x4 */
        int on_same_side_3 = 0;
        vector_in_plane_3_x = x4[0] - x3[0];
        vector_in_plane_3_y = x4[1] - x3[1];
        vector_in_plane_3_z = x4[2] - x3[2];
        LSM_GEOM_3D_SAME_SIDE(on_same_side_3,
          x_intersect, y_intersect, z_intersect,
          x2[0], x2[1], x2[2],
          x3[0], x3[1], x3[2],
          vector_in_plane_3_x, vector_in_plane_3_y, vector_in_plane_3_z);

        if (on_same_side_3) {  /* point lies within face of tetrahedron */

          if (intersect_coef_max < intersect_coef) {
            intersect_coef_max = intersect_coef;
            endpt2[0] = x_intersect; 
            endpt2[1] = y_intersect; 
            endpt2[2] = z_intersect; 
          }

          if (intersect_coef_min > intersect_coef) {
            intersect_coef_min = intersect_coef;
            endpt1[0] = x_intersect; 
            endpt1[1] = y_intersect; 
            endpt1[2] = z_intersect; 
          }

          count++;  /* increment intersection counter */

/* KTC - DEBUGGING */
/*
  printf("(x2,x3,x4) face\n");
  printf("intersect_coef = %f\n", intersect_coef);
  printf("endpt1 = (%f,%f,%f)\n", endpt1[0], endpt1[1], endpt1[2]);
  printf("endpt2 = (%f,%f,%f)\n", endpt2[0], endpt2[1], endpt2[2]);
*/

        }
      }
    } /* end point in triangle test */
  } /* end case: {phi=0, psi=0} line intersects face defined by x2, x3, x4 */

  return count;
}

