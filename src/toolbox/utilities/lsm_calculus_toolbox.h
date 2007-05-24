/*
 * File:        lsm_calculus_toolbox.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.9 $
 * Modified:    $Date: 2006/05/19 14:55:04 $
 * Description: Header file level set method calculus toolbox functions
 */

#ifndef INCLUDED_LSM_CALCULUS_TOOLBOX
#define INCLUDED_LSM_CALCULUS_TOOLBOX

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_calculus_toolbox.h
 *
 * \brief 
 * @ref lsm_calculus_toolbox.h provides macros for computing the smoothed 
 * Heaviside function and delta-function.
 *
 */


#include <math.h>

/*! \def LSM_PI
 *
 * Value of \f$ \pi \f$
 */
#define LSM_PI     (3.14159265358979323846)

/*! 
 * LSM_HEAVISIDE() computes the value of the standard smoothed Heaviside
 * function for level set method calculations
 *
 * Arguments:
 *  - x (in):        spatial position Heaviside function evaluated at
 *  - epsilon (in):  width of numerical smoothing
 *
 * Return value:     value of Heaviside function at specified position
 */
#define LSM_HEAVISIDE(x,eps)                                            \
(                                                                       \
  ((x) < -(eps)) ? 0 : ( ((x) > (eps)) ? 1 :                            \
  0.5*( 1+(x)/(eps)+1/LSM_PI*sin(LSM_PI*(x)/(eps)) ) )                  \
)

/*! 
 * LSM_DELTA_FUNCTION() computes the value of the standard smoothed 
 * delta function for level set method calculations
 *
 * Arguments:
 *  - x (in):        spatial position delta-function evaluated at
 *  - epsilon (in):  width of numerical smoothing
 *
 * Return value:     value of delta-function at specified position
 */
#define LSM_DELTA_FUNCTION(x,eps)                                       \
(                                                                       \
  ((x) < -(eps)) ? 0 : ( ((x) > (eps)) ? 0 :                            \
  0.5/(eps)*( 1+cos(LSM_PI*(x)/(eps)) ) )                               \
)

#ifdef __cplusplus
}
#endif

#endif
