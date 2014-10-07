/*
 * File:        curvature_model3d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Support header file for 3D constant curvature flow.
 */
 
#ifndef INCLUDED_CURV_MODEL3D_H
#define INCLUDED_CURV_MODEL3D_H

void  curvatureModelMedium3dMainLoop(Options *,LSM_DataArrays  *,Grid  *,FILE *);
void  reinitializeMedium3d(LSM_DataArrays *,Grid *,Options *,LSMLIB_REAL);

#endif
