/*
 * File:        lsm_FMM_field_extension3d.c
 * Copyright:   (c) 2005-2008 Kevin T. Chu
 * Revision:    $Revision: 1.2 $
 * Modified:    $Date: 2006/08/13 13:30:52 $
 * Description: Implementation of 3D Fast Marching Method for computing
 *              signed distance functions and extension fields
 */
 

/*
 * lsm_FMM_field_extension3d.c makes use of the generic implementation 
 * of the Fast Marching Method algorithm for computing signed distance
 * functions and extension fields provided by lsm_FMM_field_extension.c.
 */

#include "lsm_fast_marching_method.h"


/* Define required macros */
#define FMM_NDIM                         3
#define FMM_COMPUTE_DISTANCE_FUNCTION    computeDistanceFunction3d
#define FMM_COMPUTE_EXTENSION_FIELDS     computeExtensionFields3d
#define FMM_INITIALIZE_FRONT_ORDER1                                         \
        FMM_initializeFront_FieldExtension3d_Order1
#define FMM_INITIALIZE_FRONT_ORDER2                                         \
        FMM_initializeFront_FieldExtension3d_Order2
#define FMM_UPDATE_GRID_POINT_ORDER1                                        \
        FMM_updateGridPoint_FieldExtension3d_Order1
#define FMM_UPDATE_GRID_POINT_ORDER2                                        \
        FMM_updateGridPoint_FieldExtension3d_Order2


/* Include "templated" implementation of Fast Marching Method */
/* signed distance functions and extension fields.            */
#include "lsm_FMM_field_extension.c"

