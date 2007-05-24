/*
 * This program tests that the LSM3D_findLineInTetrahedron() function 
 * correctly computes the coefficients of the linear approximation to the 
 * level set functions and the {phi=0,psi=0} line.  It does NOT test that 
 * the endpoints of the line intersecting the tetrahedron are correctly 
 * computed.
 *
 * The corners of the tetrahedron were arbitrarily chosen.  The values of 
 * phi and psi at the corners were computed using alpha and beta vectors
 * chosen to test the different cases in the LSM3D_findLineInTetrahedron()
 * function.
 *
 * The alpha and beta vectors are defined so that: 
 *
 * phi = alpha_0 + alpha_1*x + alpha_2*y + alpha_3*z
 * psi = beta_0 + beta_1*x + beta_2*y + beta_3*z
 *
 * NOTE: Use of this test program requires commenting out the check to
 *       see if there the {phi=0,psi=0} line intersects the tetrahedron 
 *       and putting printf() statements in the LSM3D_findLineInTetrahedron() 
 *       function to print out:
 *       (1)  the coefficients of the linear approximations to phi and psi
 *       (2)  the point on the {phi=0,psi=0} line
 *       (3)  the direction of the {phi=0,psi=0} line
 *
 * Kevin T. Chu
 * MAE, Princeton University
 * December 2005
 *
 */

#include <stdio.h>
#include <math.h>
#include "lsm_utilities3d.h"

int main(void)
{

  double endpt1[3], endpt2[3];
  double x1[3], x2[3], x3[3], x4[3];
  double phi[4], psi[4];
  double line_dir[3];
  int err_code;

  /* set coordinates and values at corners of tetrahedron */
  x1[0] = 0.29441081639264; 
  x1[1] = -0.69177570170229; 
  x1[2] = -1.44096443190102;
  x2[0] = -1.33618185793780; 
  x2[1] = 0.85799667282826; 
  x2[2] = 0.57114762365818;
  x3[0] = 0.71432455181895; 
  x3[1] = 1.25400142160253; 
  x3[2] = -0.39988557771536;
  x4[0] = 1.62356206444627; 
  x4[1] = -1.59372957644748; 
  x4[2] = 0.68999737546435;

/* {phi=0} and {psi=0} parallel with normal perpendicular to x-axis */
/* alpha = (1.0,0.0,1.0,1.0) */
/*
  phi[0] = -1.13274013360331;
  phi[1] = 2.42914429648644;
  phi[2] = 1.85411584388717;
  phi[3] = 0.09626779901687;
*/

/* beta = (0.5,0.0,0.5,0.5)  */
/*
  psi[0] = -0.56637006680165;
  psi[1] = 1.21457214824322;
  psi[2] = 0.92705792194358;
  psi[3] = 0.04813389950843;
*/

/* {phi=0} and {psi=0} with normal perpendicular to x-axis   */
/*   with unparallel normals and abs(alpha_2) > abs(beta_2)  */
/* alpha = (1.0,0.0,1.0,1.0) */
/*
  phi[0] = -1.13274013360331;
  phi[1] = 2.42914429648644;
  phi[2] = 1.85411584388717;
  phi[3] = 0.09626779901687;
*/

/* beta = (1.0,0.0,0.5,0.75)  */
/*
  psi[0] = -0.42661117477691;
  psi[1] = 1.85735905415776;
  psi[2] = 1.32708652751474;
  psi[3] = 0.72063324337452;
*/

/* beta = (1.0,0.0,-0.5,0.75)  */
/*
  psi[0] = 0.26516452692538;
  psi[1] = 0.99936238132950;
  psi[2] = 0.07308510591222;
  psi[3] = 2.31436281982200;
*/

/* {phi=0} and {psi=0} with normal perpendicular to x-axis   */
/*   with unparallel normals and abs(alpha_2) < abs(beta_2)  */
/* beta = (1.0,0.0,1.0,1.0) */
/*
  psi[0] = -1.13274013360331;
  psi[1] = 2.42914429648644;
  psi[2] = 1.85411584388717;
  psi[3] = 0.09626779901687;
*/

/* alpha = (1.0,0.0,0.5,0.75)  */
/*
  phi[0] = -0.42661117477691;
  phi[1] = 1.85735905415776;
  phi[2] = 1.32708652751474;
  phi[3] = 0.72063324337452;
*/

/* alpha = (1.0,0.0,-0.5,0.75)  */
/*
  phi[0] = 0.26516452692538;
  phi[1] = 0.99936238132950;
  phi[2] = 0.07308510591222;
  phi[3] = 2.31436281982200;
*/

/* {phi=0} and {psi=0} parallel with normal perpendicular to y-axis */
/* alpha = (1.0,1.0,0.0,1.0)  */
/*
  phi[0] = -0.14655361550838;
  phi[1] = 0.23496576572038;
  phi[2] = 1.31443897410359;
  phi[3] = 3.31355943991062;
*/

/* beta = (0.5,0.5,0.0,0.5)  */
/*
  psi[0] = -0.07327680775419;
  psi[1] = 0.11748288286019;
  psi[2] = 0.65721948705180;
  psi[3] = 1.65677971995531;
*/

/* {phi=0} and {psi=0} with normal perpendicular to y-axis   */
/*   with unparallel normals and abs(alpha_1) > abs(beta_1)  */
/* alpha = (1.0,1.0,0.0,1.0)  */
  phi[0] = -0.14655361550838;
  phi[1] = 0.23496576572038;
  phi[2] = 1.31443897410359;
  phi[3] = 3.31355943991062;

/* beta = (1.0,0.5,0.0,0.0)  */
/*
  psi[0] = 1.14720540819632;
  psi[1] = 0.33190907103110;
  psi[2] = 1.35716227590948;
  psi[3] = 1.81178103222314;
*/

/* beta = (1.0,-0.5,0.0,0.75)  */
/*
  psi[0] = -0.22792873212209;
  psi[1] = 2.09645164671254;
  psi[2] = 0.34292354080401;
  psi[3] = 0.70571699937513;
*/


/* {phi=0} and {psi=0} with normal perpendicular to y-axis   */
/*   with unparallel normals and abs(alpha_1) < abs(beta_1)  */
/* beta = (1.0,1.0,0.0,1.0) */
  psi[0] = -0.14655361550838;
  psi[1] = 0.23496576572038;
  psi[2] = 1.31443897410359;
  psi[3] = 3.31355943991062;

/* alpha = (1.0,0.5,0.0,0.0) */
/*
  phi[0] = 1.14720540819632;
  phi[1] = 0.33190907103110;
  phi[2] = 1.35716227590948;
  phi[3] = 1.81178103222314;
*/

/* alpha = (1.0,-0.5,0.0,0.75) */
/*
  phi[0] = -0.22792873212209;
  phi[1] = 2.09645164671254;
  phi[2] = 0.34292354080401;
  phi[3] = 0.70571699937513;
*/

/* {phi=0} and {psi=0} parallel with normal perpendicular to z-axis */
/* alpha = (1.0,1.0,1.0,0.0) */
/*
  phi[0] = 0.60263511469035;
  phi[1] = 0.52181481489046;
  phi[2] = 2.96832597342148;
  phi[3] = 1.02983248799879;
*/

/* beta = (0.5,0.5,0.5,0.0) */
/*
  psi[0] = 0.30131755734518;
  psi[1] = 0.26090740744523;
  psi[2] = 1.48416298671074;
  psi[3] = 0.51491624399940;
*/


/* {phi=0} and {psi=0} with normal perpendicular to z-axis   */
/*   with unparallel normals and abs(alpha_1) > abs(beta_1)  */
/* alpha = (1.0,1.0,1.0,0.0) */
  phi[0] = 0.60263511469035;
  phi[1] = 0.52181481489046;
  phi[2] = 2.96832597342148;
  phi[3] = 1.02983248799879;

/* beta = (1.0,0.5,2.0,0.0) */
/*
  psi[0] = -0.23634599520826;
  psi[1] = 2.04790241668762;
  psi[2] = 3.86516511911453;
  psi[3] = -1.37567812067182;
*/

/* beta = (1.0,-0.5,1.0,0.0) */
/*
  psi[0] = 0.16101889010139;
  psi[1] = 2.52608760179716;
  psi[2] = 1.89683914569305;
  psi[3] = -1.40551060867061;
*/


/* {phi=0} and {psi=0} with normal perpendicular to z-axis   */
/*   with unparallel normals and abs(alpha_1) < abs(beta_1)  */
/* beta = (1.0,1.0,1.0,0.0) */
  psi[0] = 0.60263511469035;
  psi[1] = 0.52181481489046;
  psi[2] = 2.96832597342148;
  psi[3] = 1.02983248799879;

/* alpha = (1.0,0.5,2.0,0.0) */
/*
  phi[0] = -0.23634599520826;
  phi[1] = 2.04790241668762;
  phi[2] = 3.86516511911453;
  phi[3] = -1.37567812067182;
*/

/* alpha = (1.0,-0.5,1.0,0.0) */
/*
  phi[0] = 0.16101889010139;
  phi[1] = 2.52608760179716;
  phi[2] = 1.89683914569305;
  phi[3] = -1.40551060867061;
*/

/* {phi=0} and {psi=0} parallel with normal in arbitrary direction */
/* alpha = (1.0,1.0,2.0,3.0) */
/*
  phi[0] = -4.41203388271500;
  phi[1] = 3.09325435869326;
  phi[2] = 3.02267066187793;
  phi[3] = 1.50609503794436;
*/

/* beta = (0.5,0.5,1.0,1.5) */
/*
  psi[0] = -2.20601694135750;
  psi[1] = 1.54662717934663;
  psi[2] = 1.51133533093896;
  psi[3] = 0.75304751897218;
*/

/* {phi=0} and {psi=0} with normals in unparallel, arbitrary directions */ 
/*   and abs(alpha_1) > abs(beta_1)  */
/* alpha = (1.0,1.0,2.0,3.0) */
/*
  phi[0] = -4.41203388271500;
  phi[1] = 3.09325435869326;
  phi[2] = 3.02267066187793;
  phi[3] = 1.50609503794436;
*/

/* beta = (1.0,0.5,1.0,2.0) */
/*
  psi[0] = -2.42649915730801;
  psi[1] = 2.33220099117572;
  psi[2] = 1.81139254208128;
  psi[3] = 1.59804620670436;
*/

/* beta = (1.0,0.5,2.0,2.0) */
/*
  psi[0] = -3.11827485901030;
  psi[1] = 3.19019766400398;
  psi[2] = 3.06539396368381;
  psi[3] = 0.00431663025688;
*/

/* beta = (1.0,-0.75,1.0,0.25) */
/*
  psi[0] = -0.27282492197202;
  psi[1] = 3.00291997219616;
  psi[2] = 1.61828661330948;
  psi[3] = -1.63890178091610;
*/

/* {phi=0} and {psi=0} with normals in unparallel, arbitrary directions */ 
/*   and abs(alpha_1) < abs(beta_1)  */
/* beta = (1.0,1.0,2.0,3.0) */
  psi[0] = -4.41203388271500;
  psi[1] = 3.09325435869326;
  psi[2] = 3.02267066187793;
  psi[3] = 1.50609503794436;

/* alpha = (1.0,0.5,1.0,2.0) */
/*
  phi[0] = -2.42649915730801;
  phi[1] = 2.33220099117572;
  phi[2] = 1.81139254208128;
  phi[3] = 1.59804620670436;
*/

/* alpha = (1.0,0.5,2.0,2.0) */
  phi[0] = -3.11827485901030;
  phi[1] = 3.19019766400398;
  phi[2] = 3.06539396368381;
  phi[3] = 0.00431663025688;

/* alpha = (1.0,-0.75,1.0,0.25) */
/*
  phi[0] = -0.27282492197202;
  phi[1] = 3.00291997219616;
  phi[2] = 1.61828661330948;
  phi[3] = -1.63890178091610;
*/

/* alpha = (0.0,-0.75,1.0,0.25) */
/*
  phi[0] = -1.27282492197202;
  phi[1] = 2.00291997219616;
  phi[2] = 0.61828661330948;
  phi[3] = -2.63890178091610;
*/


  /* find {phi=0,psi=0} line */
  err_code = LSM3D_findLineInTetrahedron(
    endpt1, endpt2,
    x1, x2, x3, x4,
    phi, psi);

  /* compute direction of line defined by endpt1 and endpt2 */
  line_dir[0] = endpt2[0] - endpt1[0];
  line_dir[1] = endpt2[1] - endpt1[1];
  line_dir[2] = endpt2[2] - endpt1[2];
  line_dir[0] /= line_dir[2];
  line_dir[1] /= line_dir[2];
  line_dir[2] /= line_dir[2];

  /* output results */
  printf("Error Code = %d\n", err_code);
  printf("Endpoint 1 = (%f,%f,%f)\n", endpt1[0], endpt1[1], endpt1[2]);
  printf("Endpoint 2 = (%f,%f,%f)\n", endpt2[0], endpt2[1], endpt2[2]);
  printf("Line direction = (%f,%f,%f)\n", 
    line_dir[0], line_dir[1], line_dir[2]);

  return 0;
}
