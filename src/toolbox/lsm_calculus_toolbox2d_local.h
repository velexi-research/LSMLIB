/*
 * File:        lsm_calculus_toolbox2d_local.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file level set method calculus toolbox functions
 */

#ifndef INCLUDED_LSM_CALCULUS_TOOLBOX2D_LOCAL
#define INCLUDED_LSM_CALCULUS_TOOLBOX2D_LOCAL

#include "lsmlib_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_calculus_toolbox2d_local.h
 *
 * \brief Provides higher order implementation of delta function in 2D.
 */


#define LSM2D_DELTA_FUNCTION_ORDER1_LOCAL lsm2ddeltafunctionorder1local_
#define LSM2D_DELTA_FUNCTION_ORDER2_LOCAL lsm2ddeltafunctionorder2local_


/*
*  LSM2D_DELTA_FUNCTION_ORDER1_LOCAL() computes first order accurate delta
*  function discretization of a zero level set for a given LS function
*  following P. Smereka, "The numerical approximation of a delta function 
*  with application to level set methods", JCP, 2006.
*  The function is supported at a minimal set of gridpoints near the zero
*  level set.
*
*  Arguments:
*    phi(in):           level set function
*    delta(out):        discretized delta function corresp. to zero level
*    norm_phi_* (in):   components of grad(phi)/|grad(phi)|, 
*                       obtained by 2nd order central diff
*    *_gb (in):        index range for ghostbox
*    dx,dy(in):        grid spacing
*    index_[xy]  (in):  [xy] coordinates of local (narrow band) points
*    n*_index    (in):  index range of points to loop over in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*/
void LSM2D_DELTA_FUNCTION_ORDER1_LOCAL(
  const LSMLIB_REAL *phi,
  const LSMLIB_REAL *delta,
  const int *ilo_gb,
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb, 
  const LSMLIB_REAL *norm_phi_x,
  const LSMLIB_REAL *norm_phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL    *dx,
  const LSMLIB_REAL    *dy,
  const int *index_x,
  const int *index_y,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb);   

  
/*
*  LSM2D_DELTA_FUNCTION_ORDER2_LOCAL() computes second order accurate delta
*  function discretization of a zero level set for a given ls function
*  following P. Smereka, "The numerical approximation of a delta function 
*  with application to level set methods", JCP, 2006.
*  The function is supported at a minimal set of gridpoints near the zero
*  level set.
*
*  Arguments:
*    phi(in):           level set function
*    delta(out):        discretized delta function corresp. to zero level
*    norm_phi_* (in):   components of grad(phi)/|grad(phi)|, 
*                       obtained by 2nd order central diff
*    *_gb (in):        index range for ghostbox
*    dx,dy(in):        grid spacing
*    index_[xy]  (in):  [xy] coordinates of local (narrow band) points
*    n*_index    (in):  index range of points to loop over in index_*
*    narrow_band(in):   array that marks voxels outside desired fillbox
*    mark_fb(in):      upper limit narrow band value for voxels in 
*                      fillbox
*/
void LSM2D_DELTA_FUNCTION_ORDER2_LOCAL(
  const LSMLIB_REAL *phi,
  const LSMLIB_REAL *delta,
  const int *ilo_gb,
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb, 
  const LSMLIB_REAL *norm_phi_x,
  const LSMLIB_REAL *norm_phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL    *dx,
  const LSMLIB_REAL    *dy,
  const int *index_x,
  const int *index_y,
  const int *nlo_index,
  const int *nhi_index,
  const unsigned char *narrow_band,
  const int *ilo_nb_gb,
  const int *ihi_nb_gb,
  const int *jlo_nb_gb,
  const int *jhi_nb_gb,
  const unsigned char *mark_fb); 
  
  
#ifdef __cplusplus
}
#endif

#endif
