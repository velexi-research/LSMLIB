/*
 * File:        lsm_FMM_eikonal2d.c
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Implementation of 2D Fast Marching Method for Eikonal equation
 */


/* 
 * lsm_FMM_eikonal2d.c makes use of the generic implementation of
 * the Eikonal equation solver based on the Fast Marching Method
 * provided by lsm_FMM_eikonal.c.
 */

#include "lsm_fast_marching_method.h"


/* Define required macros */
#define FMM_NDIM                               2 
#define FMM_EIKONAL_SOLVE_EIKONAL_EQUATION     solveEikonalEquation2d
#define FMM_EIKONAL_INITIALIZE_FRONT           FMM_initializeFront_Eikonal2d
#define FMM_EIKONAL_UPDATE_GRID_POINT_ORDER1                              \
        FMM_updateGridPoint_Eikonal2d_Order1
#define FMM_EIKONAL_UPDATE_GRID_POINT_ORDER2                              \
        FMM_updateGridPoint_Eikonal2d_Order2


/* Include "templated" implementation of Eikonal equation solver. */
#include "lsm_FMM_eikonal.c"
