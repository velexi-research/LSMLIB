/*
 * lsm_FMM_field_extension2d.c makes use of the generic implementation 
 * of the Fast Marching Method algorithm for computing signed distance
 * functions and extension fields provided by lsm_FMM_field_extension.c.
 */

#include "lsmlib/lsm_fast_marching_method.h"


/* Define required macros */
#define FMM_NDIM                         2
#define FMM_COMPUTE_DISTANCE_FUNCTION    computeDistanceFunction2d
#define FMM_COMPUTE_EXTENSION_FIELDS     computeExtensionFields2d
#define FMM_INITIALIZE_FRONT_ORDER1                                         \
        FMM_initializeFront_FieldExtension2d_Order1
#define FMM_INITIALIZE_FRONT_ORDER2                                         \
        FMM_initializeFront_FieldExtension2d_Order2
#define FMM_UPDATE_GRID_POINT_ORDER1                                        \
        FMM_updateGridPoint_FieldExtension2d_Order1
#define FMM_UPDATE_GRID_POINT_ORDER2                                        \
        FMM_updateGridPoint_FieldExtension2d_Order2


/* Include "templated" implementation of Fast Marching Method */
/* signed distance functions and extension fields.            */
#include "lsmlib/lsm_FMM_field_extension.h"
