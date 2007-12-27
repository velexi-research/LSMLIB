/*
 * File:        lsm_reinitialization_medium2d.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision: 1.2 $
 * Modified:    $Date: 2006/09/18 20:33:55 $
 * Description: Header file for Options structure implentation functions.
 */

#ifndef INCLUDED_LSM_REINITIALIZATION_MEDIUM2D
#define INCLUDED_LSM_REINITIALIZATION_MEDIUM2D

#ifdef __cplusplus
extern "C" {
#endif


/*! \file lsm_reinitialization_medium2d.h
 *
 * \brief
 * LSMLIB provides reinitialization functions.
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
