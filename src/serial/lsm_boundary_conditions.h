/*
 * File:        lsm_boundary_conditions.h
 * Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for functions for imposing boundary conditions
 *              for serial calculations
 */

#ifndef INCLUDED_LSM_BOUNDARY_CONDITIONS_H
#define INCLUDED_LSM_BOUNDARY_CONDITIONS_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "lsm_grid.h"


/*! \file lsm_boundary_conditions.h
 *
 * \brief
 * @ref lsm_boundary_conditions.h provides support for filling 
 * ghostcells to impose several common level set method boundary 
 * conditions in two- and three-dimensions.  Support is provided 
 * for extrapolation and homogeneous Neumann boundary conditions.
 * Boundary conditions are imposed by filling ghostcells outside 
 * of the computational domain in such a way that they produce
 * the desired boundary condition.  
 *
 */

/* @enum BOUNDARY_LOCATION_IDX
 *
 * The boundary location index is used to indicate which boundaries 
 * should have the specified boundary condition imposed.  The 
 * following conventions are used for the boundary location index:
 *
 *
 * 2D
 * --
 *  - x_lo:           0
 *  - x_hi:           1
 *  - y_lo:           2
 *  - y_hi:           3
 *  - x_lo and x_hi:  6
 *  - y_lo and y_hi:  7
 *  - all boundaries: 9
 *
 *
 * 3D
 * --
 *
 *  - x_lo:           0
 *  - x_hi:           1
 *  - y_lo:           2
 *  - y_hi:           3
 *  - z_lo:           4
 *  - z_hi:           5
 *  - x_lo and x_hi:  6
 *  - y_lo and y_hi:  7
 *  - z_lo and z_hi:  8
 *  - all boundaries: 9
 * 
 */
typedef enum {
  X_LO           = 0,
  X_HI           = 1,
  Y_LO           = 2,
  Y_HI           = 3,
  Z_LO           = 4,
  Z_HI           = 5,
  X_LO_AND_X_HI  = 6,
  Y_LO_AND_Y_HI  = 7,
  Z_LO_AND_Z_HI  = 8,
  ALL_BOUNDARIES = 9} BOUNDARY_LOCATION_IDX;

/*!
 * linearExtrapolationBC() extrapolates data from interior grid
 * cells into the ghostcells at the specified boundary location(s)
 * using linear extrapolation.
 *
 * Arguments:
 *  - phi (in/out):            grid function for which to set ghostcells
 *  - grid (in):               pointer to Grid data structure
 *  - bdry_location_idx (in):  boundary location index
 *      
 * Return value:               none
 *      
 */
void linearExtrapolationBC(
  LSMLIB_REAL *phi,
  Grid *grid,
  int bdry_location_idx);


/*!
 * signedLinearExtrapolationBC() extrapolates data from interior grid
 * cells into the ghostcells at the specified boundary location(s) 
 * using "signed linear extrapolation."  This function is almost the 
 * same as linear extrapolation except that the sign of the slope 
 * used in extrapolation is chosen to ensure that the value of the 
 * level set function in the ghostcells move "away" from the zero 
 * level set.  This choice ensures that we do not accidentally introduce 
 * spurious zero level sets when we impose boundary conditions.
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
 *  - grid (in):               pointer to Grid data structure
 *  - bdry_location_idx (in):  boundary location index
 *      
 * Return value:               none
 *      
 */
void signedLinearExtrapolationBC(
  LSMLIB_REAL *phi,
  Grid *grid,
  int bdry_location_idx);
 
   
/*!
 * copyExtrapolationBC() trivially extrapolates data from interior 
 * grid cells into the ghostcells at the specified boundary 
 * location(s) by copying data from the closest grid cell in 
 * the interior.
 *
 * Arguments:
 *  - phi (in/out):            grid function for which to set ghostcells
 *  - grid (in):               pointer to Grid data structure
 *  - bdry_location_idx (in):  boundary location index
 *      
 * Return value:               none
 *      
 */
void copyExtrapolationBC(
  LSMLIB_REAL *phi,
  Grid *grid,
  int bdry_location_idx);


/*!
 * homogeneousNeumannBC() sets the values of phi in the ghostcells to 
 * impose a homogeneous Neumann boundary condition at the specified 
 * boundary location(s) for ENO1, ENO2, ENO3, and WENO5 discretizations
 * of the spatial derivative.  In all of these cases, it is sufficient
 * to fill ghostcells with the value of phi from the nearest grid
 * cell in the interior of the computational domain (i.e. copy 
 * extrapolation).
 * 
 * Arguments:
 *  - phi (in/out):            grid function for which to set ghostcells
 *  - grid (in):               pointer to Grid data structure
 *  - bdry_location_idx (in):  boundary location index
 *      
 * Return value:               none
 *      
 */
void homogeneousNeumannBC(
  LSMLIB_REAL *phi,
  Grid *grid,
  int bdry_location_idx);


#ifdef __cplusplus
}
#endif

#endif
