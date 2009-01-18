/*
 * File:        LSMLIB_DefaultParameters.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date$
 * Description: Header file containing default parameter values for the
 *              level set method classes
 */

#ifndef included_LSMLIB_DefaultParameters_h
#define included_LSMLIB_DefaultParameters_h

/*
 * Default values for user-defined parameters
 */
#define LSM_DIM_MAX                                      (3)
#define LSM_DEFAULT_SPATIAL_DERIVATIVE_TYPE              "WENO"
#define LSM_DEFAULT_SPATIAL_DERIVATIVE_WENO_ORDER        (5)
#define LSM_DEFAULT_SPATIAL_DERIVATIVE_ENO_ORDER         (3)
#define LSM_DEFAULT_TVD_RUNGE_KUTTA_ORDER                (3)
#define LSM_DEFAULT_CFL_NUMBER                           (0.5)
#define LSM_DEFAULT_VERBOSE_MODE                         (false)

#endif
