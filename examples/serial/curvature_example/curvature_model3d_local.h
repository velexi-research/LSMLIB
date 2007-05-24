#ifndef INCLUDED_CURV_MODEL3D_LOCAL_H
#define INCLUDED_CURV_MODEL3D_LOCAL_H

void  curvatureModelMedium3dLocalMainLoop(Options *,LSM_DataArrays  *,Grid  *,FILE *);
void  reinitializeMedium3dLocal(LSM_DataArrays *,Grid *,Options *,double);

#endif
