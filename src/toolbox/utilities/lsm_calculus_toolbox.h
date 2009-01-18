/*
 * File:        lsm_calculus_toolbox.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file level set method calculus toolbox functions
 */

#ifndef INCLUDED_LSM_CALCULUS_TOOLBOX
#define INCLUDED_LSM_CALCULUS_TOOLBOX

#include "LSMLIB_config.h"

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

/*! 
 * LSM_HEAVISIDE() macro computes the value of the standard Heaviside
 * function for level set method calculations (smoothed using sine function).
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
  0.5*( 1+(x)/(eps)+1/M_PI*sin(M_PI*(x)/(eps)) ) )                  \
)

/*! 
 * LSM_DELTA_FUNCTION() macro computes the value of the standard
 * delta function for level set method calculations (smoothed using 
 * cosine function).
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
  0.5/(eps)*( 1+cos(M_PI*(x)/(eps)) ) )                               \
)

/*! 
 * LSM_HEAVISIDE_HAT() macro computes the value of the standard "hat"
 *  Heaviside function for level set method calculations
 *
 * Arguments:
 *  - x (in):        spatial position Heaviside function evaluated at
 *  - epsilon (in):  width of numerical smoothing
 *
 * Return value:     value of Heaviside function at specified position
 */
#define LSM_HEAVISIDE_HAT(x,eps)                                        \
(                                                                       \
  ((x) < -(eps)) ? 0 : ( ((x) > (eps)) ? 1 :                            \
  ( ((x) > (zero)) ?  0.5*((x)+eps)*((x)+eps)/(eps*eps) :               \
      1.0 - 0.5*(eps-x)*(eps-x)/(eps*eps)  ) )                          \
)

/*! 
 * LSM_DELTA_FUNCTION_HAT() macro computes the value of the standard "hat" 
 * delta function for level set method calculations
 *
 * Arguments:
 *  - x (in):        spatial position delta-function evaluated at
 *  - epsilon (in):  width of numerical smoothing
 *
 * Return value:     value of delta-function at specified position
 */
#define LSM_DELTA_FUNCTION_HAT(x,eps)                                   \
(                                                                       \
  ((x) < -(eps)) ? 0 : ( ((x) > (eps)) ? 0 :                            \
  ( eps - fabs(x) )/(eps*eps) )                                         \
)

#ifdef __cplusplus
}
#endif

#endif
