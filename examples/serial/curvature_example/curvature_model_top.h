/*
 * File:        curvature_model_top.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Support header file for constant curvature flow.
 */
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
