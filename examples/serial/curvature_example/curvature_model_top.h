#ifndef INCLUDED_CURV_MODEL_TOP_H
#define INCLUDED_CURV_MODEL_TOP_H

#define TOL       1e-16
#define EMAX_STOP 0.05
#define TPLOT     0.1

#include "lsm_options.h"
#include "lsm_data_arrays.h"

int    curvatureModelTop(Options *,char *,char *,char *);
void   setArrayAllocationCurvatureModel(Options *, LSM_DataArrays *);
Grid  *createMaskThroatFromSpheres3d(LSMLIB_REAL**,Options *);

#endif
