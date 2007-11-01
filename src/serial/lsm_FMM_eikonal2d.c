/*
 * File:        lsm_FMM_eikonal2d.c
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date: 2006/08/13 13:30:52 $
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
