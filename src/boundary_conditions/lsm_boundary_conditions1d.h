/*
 * File:        lsm_boundary_conditions1d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for 1D boundary conditions functions 
 */

#ifndef INCLUDED_LSM_BOUNDARY_CONDITIONS_1D_H
#define INCLUDED_LSM_BOUNDARY_CONDITIONS_1D_H

#include "lsmlib_config.h"

#ifdef __cplusplus
extern "C" {
#endif


/*! \file lsm_boundary_conditions1d.h
 *
 * \brief
 * @ref lsm_boundary_conditions1d.h provides support for filling 
 * ghostcells to impose common boundary conditions.  Support is 
 * provided for extrapolation and homogeneous Neumann boundary 
 * conditions.
 *
 * The boundary location index is used to identify the location of the
 * boundary relative to the computational domain.  In 1D, the boundary
 * location index conventions are:
 *
 *   x_lo: 0
 *   x_hi: 1
 *
 */


/* Link between C/C++ and Fortran function names
 * 
 *      name in                            name in 
 *      C/C++ code                         Fortran code
 *      ----------                         ------------
 */
#define LSM1D_LINEAR_EXTRAPOLATION         lsm1dlinearextrapolation_
#define LSM1D_SIGNED_LINEAR_EXTRAPOLATION  lsm1dsignedlinearextrapolation_
#define LSM1D_COPY_EXTRAPOLATION           lsm1dcopyextrapolation_
#define LSM1D_HOMOGENEOUS_NEUMANN_ENO1     lsm1dhomogeneousneumanneno1_
#define LSM1D_HOMOGENEOUS_NEUMANN_ENO2     lsm1dhomogeneousneumanneno2_
#define LSM1D_HOMOGENEOUS_NEUMANN_ENO3     lsm1dhomogeneousneumanneno3_
#define LSM1D_HOMOGENEOUS_NEUMANN_WENO5    lsm1dhomogeneousneumannweno5_


/*!
 * LSM1D_LINEAR_EXTRAPOLATION() extrapolates data from interior grid
 * cells into the ghostcells at the specified boundary location using 
 * linear extrapolation.
 *
 * Arguments:
 *  - phi (in/out):            grid function for which to set ghostcells
 *  - bdry_location_idx (in):  boundary location index
 *  - *_gb (in):               index range for ghostbox
 *  - *_fb (in):               index range for fillbox
 *      
 * Return value:               none
 */
void LSM1D_LINEAR_EXTRAPOLATION(
  LSMLIB_REAL *phi,
  const int *ilo_gb, const int *ihi_gb,
  const int *ilo_fb, const int *ihi_fb,
  const int *bdry_location_idx);


/*!
 * LSM1D_SIGNED_LINEAR_EXTRAPOLATION() extrapolates data from interior grid
 * cells into the ghostcells at the specified boundary location using 
 * "signed linear extrapolation."  This function is almost the same as 
 * linear extrapolation except that the sign of the slope used in 
 * extrapolation is chosen to ensure that the value of the level set 
 * function in the ghostcells move "away" from the zero level set.  This 
 * choice ensures that we do not accidentally introduce spurious zero level 
 * sets when we impose boundary conditions.
 * 
 * Mathematically, the extrapolation function is 
 *
 * \f[
 *
 *   \phi(s) = \phi(s_b) + sgn(\phi) |d\phi/dn (s_b)| (s - s_b)
 *
 * \f]
 *
 * where \f$ s \f$ is the coordinate normal to the boundary and \f$ s_b \f$
 * is the value of \f$ s \f$ on the boundary.
 *
 * Arguments:
 *  - phi (in/out):            grid function for which to set ghostcells
 *  - bdry_location_idx (in):  boundary location index
 *  - *_gb (in):               index range for ghostbox
 *  - *_fb (in):               index range for fillbox
 *      
 * Return value:               none
 */
void LSM1D_SIGNED_LINEAR_EXTRAPOLATION(
  LSMLIB_REAL *phi,
  const int *ilo_gb, const int *ihi_gb,
  const int *ilo_fb, const int *ihi_fb,
  const int *bdry_location_idx);
 
   
/*!
 * LSM1D_COPY_EXTRAPOLATION() trivially extrapolates data from interior 
 * grid cells into the ghostcells at the specified boundary location by 
 * copying data from the closest grid cell in the interior.
 *
 * Arguments:
 *  - phi (in/out):            grid function for which to set ghostcells
 *  - bdry_location_idx (in):  boundary location index
 *  - *_gb (in):               index range for ghostbox
 *  - *_fb (in):               index range for fillbox
 *      
 * Return value:               none
 */
void LSM1D_COPY_EXTRAPOLATION(
  LSMLIB_REAL *phi,
  const int *ilo_gb, const int *ihi_gb,
  const int *ilo_fb, const int *ihi_fb,
  const int *bdry_location_idx);


/*!
 * LSM1D_HOMOGENEOUS_NEUMANN_ENO1() sets the values of phi in the 
 * ghostcells to impose a homogeneous Neumann boundary condition at 
 * the specified boundary location for an ENO1 discretization of the 
 * derivative.  In this case, the boundary condition reduces to copy 
 * extrapolation.
 *            
 * Arguments:
 *  - phi (in/out):            grid function for which to set ghostcells
 *  - bdry_location_idx (in):  boundary location index
 *  - *_gb (in):               index range for ghostbox
 *  - *_fb (in):               index range for fillbox
 *      
 * Return value:               none
 */
void LSM1D_HOMOGENEOUS_NEUMANN_ENO1(
  LSMLIB_REAL *phi,
  const int *ilo_gb, const int *ihi_gb,
  const int *ilo_fb, const int *ihi_fb,
  const int *bdry_location_idx);


/*!
 * LSM1D_HOMOGENEOUS_NEUMANN_ENO2() sets the values of phi in the 
 * ghostcells to impose a homogeneous Neumann boundary condition at 
 * the specified boundary location for an ENO2 discretization of the 
 * derivative.  
 *            
 * Arguments:
 *  - phi (in/out):            grid function for which to set ghostcells
 *  - bdry_location_idx (in):  boundary location index
 *  - *_gb (in):               index range for ghostbox
 *  - *_fb (in):               index range for fillbox
 *      
 * Return value:               none
 */
void LSM1D_HOMOGENEOUS_NEUMANN_ENO2(
  LSMLIB_REAL *phi,
  const int *ilo_gb, const int *ihi_gb,
  const int *ilo_fb, const int *ihi_fb,
  const int *bdry_location_idx);


/*!
 * LSM1D_HOMOGENEOUS_NEUMANN_ENO3() sets the values of phi in the 
 * ghostcells to impose a homogeneous Neumann boundary condition at 
 * the specified boundary location for an ENO3 discretization of the 
 * derivative.  
 *            
 * Arguments:
 *  - phi (in/out):            grid function for which to set ghostcells
 *  - bdry_location_idx (in):  boundary location index
 *  - *_gb (in):               index range for ghostbox
 *  - *_fb (in):               index range for fillbox
 *      
 * Return value:               none
 */
void LSM1D_HOMOGENEOUS_NEUMANN_ENO3(
  LSMLIB_REAL *phi,
  const int *ilo_gb, const int *ihi_gb,
  const int *ilo_fb, const int *ihi_fb,
  const int *bdry_location_idx);


/*!
 * LSM1D_HOMOGENEOUS_NEUMANN_WENO5() sets the values of phi in the 
 * ghostcells to impose a homogeneous Neumann boundary condition at 
 * the specified boundary location for an WENO5 discretization of the 
 * derivative.  
 *            
 * Arguments:
 *  - phi (in/out):            grid function for which to set ghostcells
 *  - bdry_location_idx (in):  boundary location index
 *  - *_gb (in):               index range for ghostbox
 *  - *_fb (in):               index range for fillbox
 *      
 * Return value:               none
 */
void LSM1D_HOMOGENEOUS_NEUMANN_WENO5(
  LSMLIB_REAL *phi,
  const int *ilo_gb, const int *ihi_gb,
  const int *ilo_fb, const int *ihi_fb,
  const int *bdry_location_idx);


#ifdef __cplusplus
}
#endif

#endif
