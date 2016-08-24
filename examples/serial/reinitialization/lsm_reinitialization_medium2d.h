/*
 * File:        lsm_reinitialization_medium2d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Support header file for reinitialization.
 */

#ifndef INCLUDED_LSM_REINITIALIZATION_MEDIUM2D
#define INCLUDED_LSM_REINITIALIZATION_MEDIUM2D

#ifdef __cplusplus
extern "C" {
#endif


/*! \file lsm_reinitialization_medium2d.h
 *
 * \brief
 * Provides reinitialization functions.
 */

#include "LSMLIB_config.h"
#include "lsm_grid.h"
#include "lsm_data_arrays.h"
#include "lsm_options.h"

/* lsm2dReinitializationMedium() reinitializes the level set 
 *   (i.e. near the zero level set it replaces the level set with a signed
 *   distance function).
 *
 *  Arguments:
 *  - lsm_arrays(in):  pointer to LSM_DataArrays structure
 *  - grid (in):       pointer to Grid
 *  - options (in):    pointer to Options structure
 *
 *  Return value:    none
 *   
 *  NOTES
 *    - Maximum running time for this function will be set to options->tmax.
 *      In reinitialization equation, the information travels at speed 1.
 *      Let dx[i] be grid spacing in direction i.
 *      tmax = M* max(dx[0],dx[1],dx[2]) means that we want to compute signed
 *      distance function in a band of M cells around zero level set. Therefore 
 *      it is suggested to set 'options->tmax' to a multiple of grid spacing.
 * 
 */
void lsm2dReinitializationMedium(LSM_DataArrays *lsm_arrays,Grid *grid,
                                                           Options *options);

#ifdef __cplusplus
}
#endif

#endif
