/*
 * This program tests that the LSM3D_findLineInTetrahedron() function 
 * correctly finds the endpoints of the {phi=0,psi=0} line in the specified
 * tetrahedron.
 *
 * The linear functions for phi and psi were arbitrarily chosen by 
 * selecting random alpha and beta vectors defined by:
 *
 * phi = alpha_0 + alpha_1*x + alpha_2*y + alpha_3*z
 * psi = beta_0 + beta_1*x + beta_2*y + beta_3*z
 *
 * The corners of the test tetrahedra were chosen to test the correctness 
 * of the endpoints computed by the LSM3D_findLineInTetrahedron() function 
 * for all possible cases of intersection of a line with a tetrahedron.
 *
 * Kevin T. Chu
 * MAE, Princeton University
 *
 * =================================
 *
 * CHANGE_LOG
 * ----------
 * 2005/12: original version of code
 * 2005/05/22: error checks more automated
 *
 */

#include <stdio.h>
#include <math.h>
#include "lsm_utilities3d.h"

/* Helper Function Declarations */
void checkCase(
  double *x1, double *x2, double *x3, double *x4, 
  double *phi, double *psi, 
  double *point_on_line, double *true_line_dir);


int main(void)
{

  double alpha[4], beta[4];
  double x1[3], x2[3], x3[3], x4[3];
  double phi[4], psi[4];
  double point_on_line[3], true_line_dir[3];

  /* set alpha and beta values */
  /* NOTE: for this choice of alpha and beta, the {phi=0,psi=0} */
  /*       line is (-0.75, -1.0, 0.0) + t (-1.625, -2.5, 1.0)   */
  alpha[0] = 0.5; alpha[1] = 1.0; alpha[2] = -0.25; alpha[3] = 1.0;
  beta[0] = 0.0; beta[1] = -2.0; beta[2] = 1.5; beta[3] = 0.5;


  /* {phi=0,psi=0} line approximately goes through single corner  */
  /* of tetrahedron                                               */
  printf("\nTesting case: {phi=0,psi=0} line goes through single corner\n");
  x1[0] = -0.75;
  x1[1] = -1.0;
  x1[2] = 0.0;
  x2[0] = -0.25;
  x2[1] = -1.0;
  x2[2] = 0.0;
  x3[0] = -0.75;
  x3[1] = -1.5;
  x3[2] = 0.0;
  x4[0] = -0.75;
  x4[1] = -1.0;
  x4[2] = -0.5;

  phi[0] = 0.0;
  phi[1] = 0.5;
  phi[2] = 0.125;
  phi[3] = -0.5;

  psi[0] = 0.0;
  psi[1] = -1.0;
  psi[2] = -0.75;
  psi[3] = -0.25;

  point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
  true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line approximately goes through two corners of */
  /* tetrahedron                                                  */
  printf("\nTesting case: {phi=0,psi=0} line goes through two corners\n");
  x1[0] = -0.75;
  x1[1] = -1.0;
  x1[2] = 0.0;
  x2[0] = -0.25;
  x2[1] = -1.0;
  x2[2] = 0.0;
  x3[0] = 0.0625;
  x3[1] = 0.25;
  x3[2] = -0.5;
  x4[0] = -0.75;
  x4[1] = -1.0;
  x4[2] = -0.5;

  phi[0] = 0.0;
  phi[1] = 0.5;
  phi[2] = 0.0;
  phi[3] = -0.5;

  psi[0] = 0.0;
  psi[1] = -1.0;
  psi[2] = 0.0;
  psi[3] = -0.25;

  point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
  true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line approximately goes through one corner and */
  /* non-adjacent edge of tetrahedron                             */
  printf("\nTesting case: {phi=0,psi=0} line goes through one corner ");
  printf("and one non-adjacent edge\n");
  x1[0] = -0.75;
  x1[1] = -1.0;
  x1[2] = 0.0;
  x2[0] = 3.5;
  x2[1] = 4.0;
  x2[2] = -2.0;
  x3[0] = 1.5;
  x3[1] = 4.0;
  x3[2] = -2.0;
  x4[0] = -0.75;
  x4[1] = -1.0;
  x4[2] = 0.25;

  phi[0] = 0.0;
  phi[1] = 1.0;
  phi[2] = -1.0;
  phi[3] = 0.25;

  psi[0] = 0.0;
  psi[1] = -2.0;
  psi[2] = 2.0;
  psi[3] = 0.125;

  point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
  true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line approximately goes through one corner and */
  /* opposing face of tetrahedron                                 */
  printf("\nTesting case: {phi=0,psi=0} line goes through a corner ");
  printf("and the opposing face\n");
  x1[0] = -0.75;
  x1[1] = -1.0;
  x1[2] = 0.0;
  x2[0] = -0.25;
  x2[1] = -1.0;
  x2[2] = 0.0;
  x3[0] = -0.75;
  x3[1] = -0.5;
  x3[2] = 0.0;
  x4[0] = -0.75;
  x4[1] = -1.0;
  x4[2] = -0.5;

  phi[0] = 0.0;
  phi[1] = 0.5;
  phi[2] = -0.125;
  phi[3] = -0.5;

  psi[0] = 0.0;
  psi[1] = -1.0;
  psi[2] = 0.75;
  psi[3] = -0.25;

  point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
  true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line approximately goes through a single point on */
  /* a single edge                                                   */
  printf("\nTesting case: {phi=0,psi=0} line goes through a single point ");
  printf("on a single edge\n");
  x1[0] = -1.0;
  x1[1] = -1.0;
  x1[2] = 0.0;
  x2[0] = -0.5;
  x2[1] = -1.0;
  x2[2] = 0.0;
  x3[0] = -0.75;
  x3[1] = -0.5;
  x3[2] = 0.0;
  x4[0] = -0.75;
  x4[1] = -1.0;
  x4[2] = 1.0;

  phi[0] = -0.25;
  phi[1] = 0.25;
  phi[2] = -0.125;
  phi[3] = 1.0;

  psi[0] = 0.5;
  psi[1] = -0.5;
  psi[2] = 0.75;
  psi[3] = 0.5;

  point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
  true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line approximately goes through two adjacent edges */
  printf("\nTesting case: {phi=0,psi=0} line goes through two adjacent ");
  printf("edges\n");
  x1[0] = -3.0;
  x1[1] = -6.0;
  x1[2] = 2.0;
  x2[0] = 3.5;
  x2[1] = 4.0;
  x2[2] = -2.0;
  x3[0] = 1.5;
  x3[1] = 4.0;
  x3[2] = -2.0;
  x4[0] = -0.75;
  x4[1] = -1.0;
  x4[2] = -2.0;

  phi[0] = 1.0;
  phi[1] = 1.0;
  phi[2] = -1.0;
  phi[3] = -2.0;

  psi[0] = -2.0;
  psi[1] = -2.0;
  psi[2] = 2.0;
  psi[3] = -1.0;

  point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
  true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line approximately goes through two non-adjacent edges */
  printf("\nTesting case: {phi=0,psi=0} line goes through two non-adjacent ");
  printf("edges\n");
  x1[0] = -3.0;
  x1[1] = -6.0;
  x1[2] = 2.0;
  x2[0] = 1.5;
  x2[1] = 4.0;
  x2[2] = -2.0;
  x3[0] = 2.5;
  x3[1] = 3.0;
  x3[2] = -2.0;
  x4[0] = 2.5;
  x4[1] = 5.0;
  x4[2] = -2.0;

  phi[0] = 1.0;
  phi[1] = -1.0;
  phi[2] = 0.25;
  phi[3] = -0.25;

  psi[0] = -2.0;
  psi[1] = 2.0;
  psi[2] = -1.5;
  psi[3] = 1.5;

  point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
  true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line approximately goes through an edge and   */
  /* a non-adjacent face                                         */
  printf("\nTesting case: {phi=0,psi=0} line goes through an edge and ");
  printf("a non-adjacent face\n");
  x1[0] = -3.0;
  x1[1] = -6.0;
  x1[2] = 2.0;
  x2[0] = 1.5;
  x2[1] = 4.0;
  x2[2] = -2.0;
  x3[0] = 2.5;
  x3[1] = 2.0;
  x3[2] = -3.0;
  x4[0] = 2.5;
  x4[1] = 5.0;
  x4[2] = -2.0;

  phi[0] = 1.0;
  phi[1] = -1.0;
  phi[2] = -0.5;
  phi[3] = -0.25;

  psi[0] = -2.0;
  psi[1] = 2.0;
  psi[2] = -3.5;
  psi[3] = 1.5;

  point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
  true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line approximately goes through two faces */
  printf("\nTesting case: {phi=0,psi=0} line goes through two faces\n");
  x1[0] = -3.0;
  x1[1] = -5.0;
  x1[2] = 2.0;
  x2[0] = 1.5;
  x2[1] = 4.0;
  x2[2] = -2.0;
  x3[0] = 2.5;
  x3[1] = 2.0;
  x3[2] = -3.0;
  x4[0] = 2.5;
  x4[1] = 5.0;
  x4[2] = -2.0;

  phi[0] = 0.75;
  phi[1] = -1.0;
  phi[2] = -0.5;
  phi[3] = -0.25;

  psi[0] = -0.5;
  psi[1] = 2.0;
  psi[2] = -3.5;
  psi[3] = 1.5;

  point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
  true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line parallel to x-axis */
  printf("\nTesting case: {phi=0,psi=0} line parallel to x-axis\n");
  x1[0] = -3.0;
  x1[1] = -6.0;
  x1[2] = 2.0;
  x2[0] = 3.5;
  x2[1] = 4.0;
  x2[2] = -2.0;
  x3[0] = 1.5;
  x3[1] = 4.0;
  x3[2] = -2.0;
  x4[0] = -0.75;
  x4[1] = -1.0;
  x4[2] = -2.0;

  phi[0] = -3.25;
  phi[1] = 2.75;
  phi[2] = 2.75;
  phi[3] = -2.25;

  psi[0] = -7.75;
  psi[1] = 6.25;
  psi[2] = 6.25;
  psi[3] = 1.25;

  point_on_line[0] = 0.0; point_on_line[1] = -0.5; point_on_line[2] = -0.25;
  true_line_dir[0] = 1.0; true_line_dir[1] = 0.0; true_line_dir[2] = 0.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line parallel to y-axis */
  printf("\nTesting case: {phi=0,psi=0} line parallel to y-axis\n");
  x1[0] = -3.0;
  x1[1] = -6.0;
  x1[2] = 2.0;
  x2[0] = 3.5;
  x2[1] = 4.0;
  x2[2] = -2.0;
  x3[0] = 1.5;
  x3[1] = 4.0;
  x3[2] = -2.0;
  x4[0] = -0.75;
  x4[1] = -1.0;
  x4[2] = -2.0;

  phi[0] = -0.25;
  phi[1] = 2.25;
  phi[2] = 0.25;
  phi[3] = -2.0;

  psi[0] = -4.75;
  psi[1] = 5.75;
  psi[2] = 3.75;
  psi[3] = 1.5;

  point_on_line[0] = -0.5; point_on_line[1] = 0.0; point_on_line[2] = -0.25;
  true_line_dir[0] = 0.0; true_line_dir[1] = 1.0; true_line_dir[2] = 0.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line parallel to z-axis */
  printf("\nTesting case: {phi=0,psi=0} line parallel to z-axis\n");
  x1[0] = -3.0;
  x1[1] = -6.0;
  x1[2] = 2.0;
  x2[0] = 3.5;
  x2[1] = 4.0;
  x2[2] = -2.0;
  x3[0] = 1.5;
  x3[1] = 4.0;
  x3[2] = -2.0;
  x4[0] = -0.75;
  x4[1] = -1.0;
  x4[2] = -2.0;

  phi[0] = -9.75;
  phi[1] = 6.75;
  phi[2] = 4.75;
  phi[3] = -2.5;

  psi[0] = 2.75;
  psi[1] = -0.75;
  psi[2] = -2.75;
  psi[3] = 0.0;

  point_on_line[0] = 0.5; point_on_line[1] = 0.25; point_on_line[2] = 0.0;
  true_line_dir[0] = 0.0; true_line_dir[1] = 0.0; true_line_dir[2] = 1.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line in xy-plane with z = 0 */
  printf("\nTesting case: {phi=0,psi=0} line in xy-plane with z = 0\n");
  x1[0] = -0.25;
  x1[1] = -0.25;
  x1[2] = -0.5;
  x2[0] = 0.75;
  x2[1] = -0.25;
  x2[2] = -0.5;
  x3[0] = -0.25;
  x3[1] = 0.75;
  x3[2] = -0.5;
  x4[0] = -0.25;
  x4[1] = -0.25;
  x4[2] = 0.5;

  phi[0] = -1.5;
  phi[1] = -0.5;
  phi[2] = -0.5;
  phi[3] =  0.5;

  psi[0] = 0.0;
  psi[1] = 1.0;
  psi[2] = 1.0;
  psi[3] = -1.0;

  point_on_line[0] = 0.0; point_on_line[1] = 0.0; point_on_line[2] = 0.0;
  true_line_dir[0] = -3.0; true_line_dir[1] = 3.0; true_line_dir[2] = 0.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  /* {phi=0,psi=0} line in xy-plane with z = 0 */
  printf("\nTesting case: {phi=0,psi=0} line in xy-plane with z = 0\n");
  x1[0] = -0.25;
  x1[1] = -0.25;
  x1[2] = -0.5;
  x2[0] = 0.75;
  x2[1] = -0.25;
  x2[2] = -0.5;
  x3[0] = -0.25;
  x3[1] = 0.75;
  x3[2] = -0.5;
  x4[0] = -0.25;
  x4[1] = -0.25;
  x4[2] = 0.5;

  phi[0] = -0.5;
  phi[1] =  0.5;
  phi[2] = -1.5;
  phi[3] =  0.5;

  psi[0] =  0.5;
  psi[1] =  1.5;
  psi[2] = -0.5;
  psi[3] = -0.5;

  point_on_line[0] = 0.0; point_on_line[1] = 0.0; point_on_line[2] = 0.0;
  true_line_dir[0] = 1.0; true_line_dir[1] = 1.0; true_line_dir[2] = 0.0;
  checkCase(x1, x2, x3, x4, phi, psi, point_on_line, true_line_dir);


  return 0;
}


void checkCase(
  double *x1, double *x2, double *x3, double *x4, 
  double *phi, double *psi, 
  double *point_on_line, double *true_line_dir)
{

  double endpt1[3], endpt2[3];
  double line_dir[3];
  double check_endpt1, check_endpt2;
  double check_line_dir;
  int same_point;
  int err_code;

  /* find {phi=0,psi=0} line */
  err_code = LSM3D_findLineInTetrahedron(
    endpt1, endpt2,
    x1, x2, x3, x4,
    phi, psi);

  printf("Error Code = %d\n", err_code);
  printf("Endpoint 1 = (%f,%f,%f)\n", endpt1[0], endpt1[1], endpt1[2]);
  printf("Endpoint 2 = (%f,%f,%f)\n", endpt2[0], endpt2[1], endpt2[2]);

  /* check if two endpoints are the same */
  same_point = ( fabs(endpt1[0] - endpt2[0])
               + fabs(endpt1[1] - endpt2[1])
               + fabs(endpt1[2] - endpt2[2]) < 1e-8) ? 1 : 0;

  /* compute direction of line defined by endpt1 and endpt2 */
  if (!same_point) {
    line_dir[0] = endpt2[0] - endpt1[0];
    line_dir[1] = endpt2[1] - endpt1[1];
    line_dir[2] = endpt2[2] - endpt1[2];

    if (fabs(line_dir[2]) > 1e-8) {
      line_dir[0] /= line_dir[2];
      line_dir[1] /= line_dir[2];
      line_dir[2] /= line_dir[2];
    } else if (fabs(line_dir[1]) > 1e-8) {
      line_dir[0] /= line_dir[1];
      line_dir[1] /= line_dir[1];
      line_dir[2] /= line_dir[1];
    } else {
      line_dir[0] /= line_dir[0];
      line_dir[1] /= line_dir[0];
      line_dir[2] /= line_dir[0];
    } 

    printf("Line direction = (%f,%f,%f)\n", 
      line_dir[0], line_dir[1], line_dir[2]);

  } else {
    printf("Two points are the same...cannot compute line direction...\n");
  }

  /* check the line direction */
  if (!same_point) {
    check_line_dir = 
      fabs(line_dir[0]*true_line_dir[1]-line_dir[1]*true_line_dir[0])
    + fabs(line_dir[0]*true_line_dir[2]-line_dir[2]*true_line_dir[0])
    + fabs(line_dir[1]*true_line_dir[2]-line_dir[2]*true_line_dir[1]); 
    if (fabs(check_line_dir) < 1e-6)
      printf("OK: check line direction: %f\n", check_line_dir);
    else
      printf("ERROR: check line direction: %f\n", check_line_dir);
  }

  /* check that points lie on intersection of two planes */
  check_endpt1 = fabs( (endpt1[0] - point_on_line[0])*true_line_dir[1] 
                      -(endpt1[1] - point_on_line[1])*true_line_dir[0] )
               + fabs( (endpt1[2] - point_on_line[2])*true_line_dir[0] 
                      -(endpt1[0] - point_on_line[0])*true_line_dir[2] ) 
               + fabs( (endpt1[1] - point_on_line[1])*true_line_dir[2] 
                      -(endpt1[2] - point_on_line[2])*true_line_dir[1] );
  check_endpt2 = fabs( (endpt2[0] - point_on_line[0])*true_line_dir[1] 
                      -(endpt2[1] - point_on_line[1])*true_line_dir[0] )
               + fabs( (endpt2[2] - point_on_line[2])*true_line_dir[0] 
                      -(endpt2[0] - point_on_line[0])*true_line_dir[2] ) 
               + fabs( (endpt2[1] - point_on_line[1])*true_line_dir[2] 
                      -(endpt2[2] - point_on_line[2])*true_line_dir[1] );

  if (fabs(check_endpt1) < 1e-6)
    printf("OK: check endpt1: %f\n", check_endpt1);
  else
    printf("ERROR: check endpt1: %f\n", check_endpt1);

  if (fabs(check_endpt2) < 1e-6)
    printf("OK: check endpt2: %f\n", check_endpt2);
  else
    printf("ERROR: check endpt2: %f\n", check_endpt2);

}
