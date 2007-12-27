#ifndef INCLUDED_CURV_MODEL3D_H
#define INCLUDED_CURV_MODEL3D_H

void  curvatureModelMedium3dMainLoop(Options *,LSM_DataArrays  *,Grid  *,FILE *);
void  reinitializeMedium3d(LSM_DataArrays *,Grid *,Options *,LSMLIB_REAL);

#endif
