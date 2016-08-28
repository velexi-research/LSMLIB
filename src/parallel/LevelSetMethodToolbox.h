/*
 * File:        LevelSetMethodToolbox.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for level set method toolbox class
 */

#ifndef included_LevelSetMethodToolbox_h
#define included_LevelSetMethodToolbox_h

/*! \class LSMLIB::LevelSetMethodToolbox
 *
 * \brief
 * LevelSetMethodToolbox is a utility class that provides support
 * for several basic level set method algorithms and computations.
 *
 * In particular, it includes functionality for the following:
 *
 *  - parallel computation of high-order ENO and WENO derivatives
 *    for arbitrary functions using the SAMRAI framework to handle
 *    the data management issues;
 *
 *  - standard first-, second-, and third-order TVD Runge-Kutta
 *    integration;
 *
 *  - computation of distance functions and extension fields
 *    (e.g. velocity) using the fast-marching method algorithm;
 *
 *  - computation of the unit normal vector (both on AND off of
 *    the zero level set);
 *
 *  - computation of volume and surface integrals over regions
 *    defined by the zero level set;
 *
 *  - computation of stable time step sizes for advection and
 *    normal velocity evolution;
 *
 *  - computation of the max norm of the difference of two fields;
 *
 *  - computation of control volumes for structured adaptive meshes; and
 *
 *  - general data management/transfer procedures.
 *
 *
 * <h3> NOTES: </h3>
 *  - Since this class is a utility class, all methods are declared
 *    static (i.e. they are essentially C function calls).  There is
 *    never be a need to instantiate an object of type
 *    LevelSetMethodToolbox.  Call the functions in this
 *    class as static functions using the syntax:
 *
 *       LevelSetMethodToolbox::function(...)
 *
 */

// Standard library headers
#include <vector>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/tbox/Array.h"
#include "boost/shared_ptr.hpp"
#include "SAMRAI/tbox/Database.h"

// LSMLIB headers
#include "LSMLIB_config.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchHierarchy ; } }

// Namespaces
using namespace std;
using namespace SAMRAI;


/******************************************************************
 *
 * LevelSetMethodToolbox "Class" Definition
 *
 ******************************************************************/

namespace LSMLIB {

/*! \enum LEVEL_SET_FCN_TYPE
 *
 * Enumerated type for the two level set functions.  Used to
 * specify which level set function to use when computing
 * extension fields.
 *
 */
typedef enum { PHI = 0, PSI = 1 } LEVEL_SET_FCN_TYPE;

/*! \enum SPATIAL_DERIVATIVE_TYPE
 *
 * Enumerated type for the different methods of computing
 * spatial derivatives.
 *
 */
typedef enum { ENO = 0, WENO = 1, UNKNOWN = 2 } SPATIAL_DERIVATIVE_TYPE;

class LevelSetMethodToolbox
{

public:

  //! @{
  /*!
   ****************************************************************
   *
   * @name Type and constant definitions
   *
   ****************************************************************/

  /*! \enum UNIT_NORMAL_TYPE
   *
   * Enumerated type for the different methods of computing the
   * spatial derivatives used to compute the unit normal vector.
   *
   * PHI_UPWIND:  grad(phi) computed by selecting upwind direction
   *              based on sign of phi.
   *              If phi > 0, upwind direction is direction where
   *              phi is smaller.
   *              If phi < 0, upwind direction is direction where
   *              phi is larger.
   *
   * AVERAGE:     grad(phi) computed by averaging plus and minus
   *              spatial derivatives
   */
  typedef enum { PHI_UPWIND = 0, AVERAGE = 1 } UNIT_NORMAL_TYPE;

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for computing spatial derivatives
   *
   ****************************************************************/

  /*!
   * computeUpwindSpatialDerivatives() computes upwind spatial derivatives
   * using the specified ENO/WENO scheme specified.
   *
   * Arguments:
   *  - hierarchy (in):                 Boost pointer to PatchHierarchy
   *                                    containing data
   *  - spatial_derivative_type (in):   type of spatial derivative calculation
   *  - spatial_derivative_order (in):  order of spatial derivative
   *  - grad_phi_handle (out):          PatchData handle for grad(phi)
   *  - phi_handle (in):                PatchData handle for phi
   *  - upwind_function_handle(in):     PatchData handle for upwinding function
   *  - phi_component (in):             component of phi for which to compute
   *                                    spatial derivatives (default = 0)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - Support is only provided for the computation of the following
   *    spatial derivatives:  ENO1 (i.e. 1st-order upwinding), ENO2,
   *    ENO3, and WENO5.
   *
   *  - It is the user's responsibility to ensure that the PatchData
   *    for phi is defined with a sufficent number of ghost cells to
   *    support calculation of the selected ENO/WENO derivative in
   *    the interior of the PatchData for grad(phi).
   *
   *  - For more details about the numerical discretizations used by
   *    computeUpwindSpatialDerivatives(), see "Level Set Methods and
   *    Dynamic Implicit Surfaces" by Osher & Fedkiw.
   *
   */
  static void computeUpwindSpatialDerivatives(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
    const int spatial_derivative_order,
    const int grad_phi_handle,
    const int phi_handle,
    const int upwind_function_handle,
    const int phi_component = 0);

  /*!
   * computePlusAndMinusSpatialDerivatives() computes the forward (plus)
   * backward (minus) approximations to the spatial derivatives using
   * the specified ENO/WENO scheme.
   *
   * Arguments:
   *  - hierarchy (in):                 Boost pointer to PatchHierarchy
   *                                    containing data
   *  - spatial_derivative_type (in):   type of spatial derivative calculation
   *  - spatial_derivative_order (in):  order of spatial derivative
   *  - grad_phi_plus_handle (out):     PatchData handle for forward
   *                                    approximation to grad(phi)
   *  - grad_phi_minus_handle (out):    PatchData handle for backward
   *                                    approximation to grad(phi)
   *  - phi_handle (in):                PatchData handle for phi
   *  - phi_component (in):             component of phi for which to compute
   *                                    spatial derivatives (default = 0)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - Support is only provided for the computation of the following
   *    spatial derivatives:  ENO1 (i.e. 1st-order upwinding), ENO2,
   *    ENO3, and WENO5.
   *
   *  - It is the user's responsibility to ensure that the PatchData
   *    for phi is defined with a sufficent number of ghost cells to
   *    support calculation of the selected ENO/WENO derivative in
   *    the interior of the PatchData for grad(phi).
   *
   *  - For more details about the numerical discretizations used by
   *    computePlusAndMinusSpatialDerivatives(), see "Level Set Methods
   *    and Dynamic Implicit Surfaces" by Osher & Fedkiw.
   *
   */
  static void computePlusAndMinusSpatialDerivatives(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
    const int spatial_derivative_order,
    const int grad_phi_plus_handle,
    const int grad_phi_minus_handle,
    const int phi_handle,
    const int phi_component = 0);

  /*!
   * computeCentralSpatialDerivatives() computes central approximations
   * to the spatial derivatives with the specified order of accuracy.
   *
   * Arguments:
   *  - hierarchy (in):                 Boost pointer to PatchHierarchy
   *                                    containing data
   *  - spatial_derivative_order (in):  order of spatial derivative
   *  - grad_phi_handle (out):          PatchData handle for grad(phi)
   *  - phi_handle (in):                PatchData handle for phi
   *  - phi_component (in):             component of phi for which to compute
   *                                    spatial derivatives (default = 0)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - Support is only provided for the computation of the following
   *    order spatial derivatives:  2nd and 4th.  If the user requests
   *    a 1st order derivative, the calculation automatically defaults
   *    2nd order.  If the user requests a 3rd order derivative, the
   *    the calculation automatically defaults 4th order.
   *
   *  - It is the user's responsibility to ensure that the PatchData
   *    for phi is defined with a sufficent number of ghost cells to
   *    support calculation of the selected central derivative in
   *    the interior of the PatchData for grad(phi).
   *
   *  - The 2nd order derivative is computed using the formula
   *
   *    \f[
   *
   *      \left( \frac{\partial \phi}{\partial x} \right)_i \approx
   *         \frac{ \phi_{i+1} - \phi_{i-1} }{ 2 dx }
   *
   *    \f]
   *
   *  - The 4th order derivative is computed using the formula
   *
   *    \f[
   *
   *      \left( \frac{\partial \phi}{\partial x} \right)_i \approx
   *         \frac{ -\phi_{i+2} + 8 \phi_{i+1} - 8 \phi_{i-1} + \phi_{i+2} }
   *              { 12 dx }
   *
   *    \f]
   *
   */
  static void computeCentralSpatialDerivatives(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int spatial_derivative_order,
    const int grad_phi_handle,
    const int phi_handle,
    const int phi_component = 0);

  //! @}


  //! @{
  /*!
   ********************************************************************
   *
   * @name Methods for time integration using TVD Runge-Kutta methods
   *
   ********************************************************************/

  /*!
   * TVDRK1Step() advances the solution through a single step of
   * the first-order Runge-Kutta method (i.e. forward euler).
   *
   * Arguments:
   *  - hierarchy (in):         Boost pointer to PatchHierarchy containing
   *                            data
   *  - u_next_handle (out):    PatchData handle for u(t+dt)
   *  - u_cur_handle (in):      PatchData handle for u(t)
   *  - rhs_handle (in):        PatchData handle for rhs(t)
   *  - dt (in):                time increment to advance u
   *  - u_next_component (in):  component of u_next to use in step
   *                            (default = 0)
   *  - u_cur_component (in):   component of u_cur to use in step
   *                            (default = 0)
   *  - rhs_component (in):     component of rhs to use in step
   *                            (default = 0)
   *
   * Return value:              none
   *
   */
  static void TVDRK1Step(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int u_next_handle,
    const int u_cur_handle,
    const int rhs_handle,
    const LSMLIB_REAL dt,
    const int u_next_component = 0,
    const int u_cur_component = 0,
    const int rhs_component = 0);

  /*!
   * TVDRK2Stage1() advances the solution through the first stage
   * of the second-order TVD Runge-Kutta method
   *
   * Arguments:
   *  - hierarchy (in):           Boost pointer to PatchHierarchy containing
   *                              data
   *  - u_stage1_handle (out):    PatchData handle for u_approx(t+dt)
   *  - u_cur_handle (in):        PatchData handle for u(t)
   *  - rhs_handle (in):          PatchData handle for rhs(t)
   *  - dt (in):                  time increment to advance u
   *  - u_stage1_component (in):  component of u_stage1 to use in step
   *                              (default = 0)
   *  - u_cur_component (in):     component of u_cur to use in step
   *                              (default = 0)
   *  - rhs_component (in):       component of rhs to use in step
   *                              (default = 0)
   *
   * Return value:                none
   *
   */
  static void TVDRK2Stage1(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int u_stage1_handle,
    const int u_cur_handle,
    const int rhs_handle,
    const LSMLIB_REAL dt,
    const int u_stage1_component = 0,
    const int u_cur_component = 0,
    const int rhs_component = 0);

  /*!
   * TVDRK2Stage2() completes advancing the solution through a
   * a single step of the second-order TVD Runge-Kutta method.
   *
   * Arguments:
   *  - hierarchy (in):           Boost pointer to PatchHierarchy containing
   *                              data
   *  - u_next_handle (out):      PatchData handle for u(t+dt)
   *  - u_stage1_handle (in):     PatchData handle for u_approx(t+dt)
   *  - u_cur_handle (in):        PatchData handle for u(t)
   *  - rhs_handle (in):          PatchData handle for rhs(t)
   *  - dt (in):                  time increment to advance u
   *  - u_next_component (in):    component of u_next to use in step
   *                              (default = 0)
   *  - u_stage1_component (in):  component of u_stage1 to use in step
   *                              (default = 0)
   *  - u_cur_component (in):     component of u_cur to use in step
   *                              (default = 0)
   *  - rhs_component (in):       component of rhs to use in step
   *                              (default = 0)
   *
   * Return value:                none
   *
   */
  static void TVDRK2Stage2(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int u_next_handle,
    const int u_stage1_handle,
    const int u_cur_handle,
    const int rhs_handle,
    const LSMLIB_REAL dt,
    const int u_next_component = 0,
    const int u_stage1_component = 0,
    const int u_cur_component = 0,
    const int rhs_component = 0);

  /*!
   * TVDRK3Stage1() advances the solution through the first stage
   * of the third-order TVD Runge-Kutta method
   *
   * Arguments:
   *  - hierarchy (in):           Boost pointer to PatchHierarchy containing
   *                              data
   *  - u_stage1_handle (out):    PatchData handle for u_approx(t+dt)
   *  - u_cur_handle (in):        PatchData handle for u(t)
   *  - rhs_handle (in):          PatchData handle for rhs(t)
   *  - dt (in):                  time increment to advance u
   *  - u_stage1_component (in):  component of u_stage1 to use in step
   *                              (default = 0)
   *  - u_cur_component (in):     component of u_cur to use in step
   *                              (default = 0)
   *  - rhs_component (in):       component of rhs to use in step
   *                              (default = 0)
   *
   * Return value:                none
   *
   */
  static void TVDRK3Stage1(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int u_stage1_handle,
    const int u_cur_handle,
    const int rhs_handle,
    const LSMLIB_REAL dt,
    const int u_stage1_component = 0,
    const int u_cur_component = 0,
    const int rhs_component = 0);

  /*!
   * TVDRK3Stage2() advances the solution through the second stage
   * of the third-order TVD Runge-Kutta method
   *
   * Arguments:
   *  - hierarchy (in):           Boost pointer to PatchHierarchy containing
   *                              data
   *  - u_stage2_handle (out):    PatchData handle for u_approx(t+dt/2)
   *  - u_stage1_handle (in):     PatchData handle for u_approx(t+dt)
   *  - u_cur_handle (in):        PatchData handle for u(t)
   *  - rhs_handle (in):          PatchData handle for rhs(t)
   *  - dt (in):                  time increment to advance u
   *  - u_stage2_component (in):  component of u_stage2 to use in step
   *                              (default = 0)
   *  - u_stage1_component (in):  component of u_stage1 to use in step
   *                              (default = 0)
   *  - u_cur_component (in):     component of u_cur to use in step
   *                              (default = 0)
   *  - rhs_component (in):       component of rhs to use in step
   *                              (default = 0)
   *
   * Return value:                none
   *
   */
  static void TVDRK3Stage2(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int u_stage2_handle,
    const int u_stage1_handle,
    const int u_cur_handle,
    const int rhs_handle,
    const LSMLIB_REAL dt,
    const int u_stage2_component = 0,
    const int u_stage1_component = 0,
    const int u_cur_component = 0,
    const int rhs_component = 0);

  /*!
   * TVDRK3Stage3() completes advancing the solution through a
   * a single step of the third-order TVD Runge-Kutta method.
   *
   * Arguments:
   *  - hierarchy (in):           Boost pointer to PatchHierarchy containing
   *                              data
   *  - u_next_handle (out):      PatchData handle for u(t+dt)
   *  - u_stage2_handle (in):     PatchData handle for u_approx(t+dt/2)
   *  - u_cur_handle (in):        PatchData handle for u(t)
   *  - rhs_handle (in):          PatchData handle for rhs(t)
   *  - dt (in):                  time increment to advance u
   *  - u_next_component (in):    component of u_next to use in step
   *                              (default = 0)
   *  - u_stage2_component (in):  component of u_stage2 to use in step
   *                              (default = 0)
   *  - u_cur_component (in):     component of u_cur to use in step
   *                              (default = 0)
   *  - rhs_component (in):       component of rhs to use in step
   *                              (default = 0)
   *
   * Return value:                none
   *
   */
  static void TVDRK3Stage3(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int u_next_handle,
    const int u_stage2_handle,
    const int u_cur_handle,
    const int rhs_handle,
    const LSMLIB_REAL dt,
    const int u_next_component = 0,
    const int u_stage2_component = 0,
    const int u_cur_component = 0,
    const int rhs_component = 0);

  //! @}

  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for computing geometric quantities
   *
   ****************************************************************/

  /*!
   * computeUnitNormalVectorFromPhi() computes the unit normal vector
   * on AND off of the zero level set using the formula
   *
   *   N = grad(phi)/|grad(phi)|
   *
   * given only data for phi.
   *
   * Arguments:
   *  - patch_hierarchy (in):           PatchHierarchy on which to
   *                                    compute stable normal velocity dt
   *  - spatial_derivative_type (in):   type of spatial derivative calculation
   *  - spatial_derivative_order (in):  order of spatial derivative
   *  - normal_vector_handle (out):     PatchData handle for normal vector
   *  - phi_handle (in):                PatchData handle for phi
   *  - unit_normal_type (in):          specifies the method to use when
   *                                    computing the spatial derivatives used
   *                                    to compute the unit normal vector
   *  - phi_component (in):             component of phi to use in computing
   *                                    normal vector (default = 0)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - It is the user's responsibility to ensure that the PatchData
   *    for phi is defined with a sufficent number of ghost cells to
   *    support calculation of the selected ENO/WENO derivative in
   *    the interior of the PatchData for the unit normal vector.
   *
   */
  static void computeUnitNormalVectorFromPhi(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
    const int spatial_derivative_order,
    const int normal_vector_handle,
    const int phi_handle,
    const UNIT_NORMAL_TYPE unit_normal_type,
    const int phi_component = 0);

  /*!
   * computeSignedUnitNormalVectorFromPhi() computes the signed unit
   * normal vector on AND off of the zero level set using the formula
   *
   *   N = sgn(phi)*grad(phi)/|grad(phi)|
   *
   * given only data for phi.
   *
   * Arguments:
   *  - patch_hierarchy (in):           PatchHierarchy on which to
   *                                    compute stable normal velocity dt
   *  - spatial_derivative_type (in):   type of spatial derivative calculation
   *  - spatial_derivative_order (in):  order of spatial derivative
   *  - normal_vector_handle (out):     PatchData handle for normal vector
   *  - phi_handle (in):                PatchData handle for phi
   *  - unit_normal_type (in):          specifies the method to use when
   *                                    computing the spatial derivatives used
   *                                    to compute the unit normal vector
   *  - phi_component (in):             component of phi to use in computing
   *                                    normal vector (default = 0)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - It is the user's responsibility to ensure that the PatchData
   *    for phi is defined with a sufficent number of ghost cells to
   *    support calculation of the selected ENO/WENO derivative in
   *    the interior of the PatchData for the unit normal vector.
   *
   *  - The following formula is used for the smoothed sgn function:
   *
   *    sgn(phi) = phi/sqrt( phi^2 + |grad(phi)|^2 * dx^2 )
   *
   */
  static void computeSignedUnitNormalVectorFromPhi(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
    const int spatial_derivative_order,
    const int normal_vector_handle,
    const int phi_handle,
    const UNIT_NORMAL_TYPE unit_normal_type,
    const int phi_component = 0);

  /*!
   * computeUnitNormalVectorFromGradPhi() computes the unit normal vector
   * on AND off of the zero level set using the formula
   *
   *   N = grad(phi)/|grad(phi)|
   *
   * given data for phi AND grad(phi).  That is, the user is free to
   * choose his/her own scheme for computing grad(phi).
   *
   * Arguments:
   *  - patch_hierarchy (in):           PatchHierarchy on which to
   *                                    compute stable normal velocity dt
   *  - normal_vector_handle (out):     PatchData handle for normal vector
   *  - grad_phi_handle (in):           PatchData handle for grad(phi)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - There are no restrictions on the number of ghostcells for any
   *    of the data.
   *
   */
  static void computeUnitNormalVectorFromGradPhi(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    const int normal_vector_handle,
    const int grad_phi_handle);

  /*!
   * computeSignedUnitNormalVectorFromGradPhi() computes the signed
   * unit normal vector on AND off of the zero level set using the
   * formula
   *
   *   N = sgn(phi)*grad(phi)/|grad(phi)|
   *
   * given data for phi AND grad(phi).  That is, the user is free to
   * choose his/her own scheme for computing grad(phi).
   *
   * Arguments:
   *  - patch_hierarchy (in):           PatchHierarchy on which to
   *                                    compute stable normal velocity dt
   *  - normal_vector_handle (out):     PatchData handle for normal vector
   *  - grad_phi_handle (in):           PatchData handle for grad(phi)
   *  - phi_handle (in):                PatchData handle for phi
   *  - phi_component (in):             component of phi that was used to
   *                                    compute grad(phi) (default = 0)
   *
   * Return value:                      none
   *
   * NOTES:
   *  - There are no restrictions on the number of ghostcells for any
   *    of the data.
   *
   *  - The following formula is used for the smoothed sgn function:
   *
   *    sgn(phi) = phi/sqrt( phi^2 + |grad(phi)|^2 * dx^2 )
   *
   */
  static void computeSignedUnitNormalVectorFromGradPhi(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    const int normal_vector_handle,
    const int grad_phi_handle,
    const int phi_handle,
    const int phi_component = 0);

  /*!
   * computeVolumeOfRegionDefinedByZeroLevelSet() computes the volume of
   * one of two regions:  region with phi < 0 or region with phi > 0.
   *
   * Arguments:
   *  - patch_hierarchy (in):        PatchHierarchy on which to compute
   *                                 the volume integral
   *  - phi_handle (in):             PatchData handle for phi
   *  - control_volume_handle (in):  PatchData handle for control volume
   *  - region_indicator(in):        integer indicating which region to
   *                                 integrate over
   *                                 - (>0):  integration region =
   *                                          {x | phi(x) > 0}
   *                                 - (<=0): integration region =
   *                                          {x | phi(x) <= 0}
   *  - phi_component (in):          component of phi to use as level set
   *                                 function (default = 0)
   *  - heaviside_width (in):        width of Heaviside function as a multiple
   *                                 of the grid spacing (default = 3)
   *
   * Return value:                   volume of the specified domain that
   *                                 is enclosed by the zero level set
   *
   */
  static LSMLIB_REAL computeVolumeOfRegionDefinedByZeroLevelSet(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    const int phi_handle,
    const int control_volume_handle,
    const int region_indicator,
    const int phi_component = 0,
    const int heaviside_width = 3);

  /*!
   * computeVolumeOfZeroLevelSet() computes the volume of the codimension-one
   * region defined by the zero level set (i.e. surface area for 3D problems,
   * perimeter for 2D problems, number of points for 1D problems).
   *
   * Arguments:
   *  - patch_hierarchy (in):        PatchHierarchy on which to compute
   *                                 the volume integral
   *  - phi_handle (in):             PatchData handle for phi
   *  - grad_phi_handle (in):        PatchData handle for grad(phi)
   *  - control_volume_handle (in):  PatchData handle for control volume
   *  - phi_component (in):          component of phi to use as level set
   *                                 function (default = 0)
   *  - delta_width (in):            width of delta-function as a multiple of
   *                                 the grid spacing (default = 3)
   *
   * Return value:                   volume of zero level set
   *
   * NOTES:
   *  - The accuracy of this calculation breaks down when zero level sets
   *    are very close together because there are not enough grid cells
   *    to compute the spatial derivatives required.
   *
   */
  static LSMLIB_REAL computeVolumeOfZeroLevelSet(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    const int phi_handle,
    const int grad_phi_handle,
    const int control_volume_handle,
    const int phi_component = 0,
    const int delta_width = 3);

  /*!
   * computeVolumeIntegral() computes the volume integral of the specified
   * function over one of two regions:  region with phi < 0 or region with
   * phi > 0.
   *
   * Arguments:
   *  - patch_hierarchy (in):        PatchHierarchy on which to compute
   *                                 the volume integral
   *  - F_handle (in):               PatchData handle for function to be
   *                                 integrated
   *  - phi_handle (in):             PatchData handle for phi
   *  - control_volume_handle (in):  PatchData handle for control volume
   *  - region_indicator(in):        integer indicating which region to
   *                                 integrate over
   *                                 - (>0):  integration region =
   *                                          {x | phi(x) > 0}
   *                                 - (<=0): integration region =
   *                                          {x | phi(x) <= 0}
   *  - F_component (in):            component of F to use as integrand
   *                                 (default = 0)
   *  - phi_component (in):          component of phi to use as level set
   *                                 function (default = 0)
   *  - heaviside_width (in):        width of Heaviside function as a multiple
   *                                 of the grid spacing (default = 3)
   *
   * Return value:                   integral of F over the specified domain
   *
   */
  static LSMLIB_REAL computeVolumeIntegral(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    const int F_handle,
    const int phi_handle,
    const int control_volume_handle,
    const int region_integrator,
    const int F_component = 0,
    const int phi_component = 0,
    const int heaviside_width = 3);

  /*!
   * computeSurfaceIntegral() computes the surface integral over the
   * specified function over the zero level set.
   *
   * Arguments:
   *  - patch_hierarchy (in):        PatchHierarchy on which to compute
   *                                 the volume integral
   *  - F_handle (in):               PatchData handle for function to be
   *                                 integrated
   *  - phi_handle (in):             PatchData handle for phi
   *  - grad_phi_handle (in):        PatchData handle for grad(phi)
   *  - control_volume_handle (in):  PatchData handle for control volume
   *  - F_component (in):            component of F to use as integrand
   *                                 (default = 0)
   *  - phi_component (in):          component of phi to use as level set
   *                                 function (default = 0)
   *  - delta_width (in):            width of delta-function as a multiple of
   *                                 the grid spacing (default = 3)
   *
   * Return value:                   integral of F over the zero level set
   *
   * NOTES:
   *  - The accuracy of this calculation breaks down when zero level sets
   *    are very close together because there are not enough grid cells
   *    to compute the spatial derivatives required.
   *
   */
  static LSMLIB_REAL computeSurfaceIntegral(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    const int F_handle,
    const int phi_handle,
    const int grad_phi_handle,
    const int control_volume_handle,
    const int F_component = 0,
    const int phi_component = 0,
    const int delta_width = 3);

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Utility methods
   *
   ****************************************************************/

  /*!
   * computeStableAdvectionDt() computes the maximum stable
   * time step size for the advective term in a time evolution
   * equation.
   *
   * Arguments:
   *  - patch_hierarchy (in):        PatchHierarchy on which to compute
   *                                 stable advection dt
   *  - velocity_handle (in):        PatchData handle for velocity field
   *                                 data
   *  - control_volume_handle (in):  PatchData handle for control volume
   *  - cfl_number (in):             CFL number
   *
   * Return value:                   none
   *
   */
  static LSMLIB_REAL computeStableAdvectionDt(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    const int velocity_handle,
    const int control_volume_handle,
    const LSMLIB_REAL cfl_number);

  /*!
   * computeStableNormalVelocityDt() computes the maximum stable
   * time step size for motion under a velocity in the normal
   * direction.
   *
   * Arguments:
   *  - patch_hierarchy (in):         PatchHierarchy on which to
   *                                  compute stable normal velocity dt
   *  - normal_velocity_handle (in):  PatchData handle for normal
   *                                  velocity field data
   *  - grad_phi_plus_handle (in):    PatchData handle for grad(phi)
   *                                  computed using forward differencing
   *  - grad_phi_minus_handle (in):   PatchData handle for grad(phi)
   *                                  computed using forward differencing
   *  - control_volume_handle (in):   PatchData handle for control volume
   *  - cfl_number (in):              CFL number
   *
   * Return value:                    none
   *
   * NOTES:
   *  - KTC - fill me in
   *
   */
  static LSMLIB_REAL computeStableNormalVelocityDt(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    const int normal_velocity_handle,
    const int grad_phi_plus_handle,
    const int grad_phi_minus_handle,
    const int control_volume_handle,
    const LSMLIB_REAL cfl_number);

  /*!
   * maxNormOfDifference() computes the max norm of the difference
   * two scalar fields.
   *
   * Arguments:
   *  - hierarchy (in):              Boost pointer to PatchHierarchy containing
   *                                 data
   *  - field1_handle (in):          PatchData handle for field1
   *  - field2_handle (in):          PatchData handle for field2
   *  - control_volume_handle (in):  PatchData handle for control volume
   *  - field1_component (in):       component of field1 to use
   *                                 (default = 0)
   *  - field2_component (in):       component of field2 to use
   *                                 (default = 0)
   *
   * Return value:                   max norm of (field1 - field2)
   *
   */
  static LSMLIB_REAL maxNormOfDifference(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int field1_handle,
    const int field2_handle,
    const int control_volume_handle,
    const int field1_component = 0,
    const int field2_component = 0);

  /*!
   * computeControlVolumes() computes the control volumes for the
   * cells in the specified PatchHierarchy.
   *
   * Arguments:
   *  - control_volume_handle (out):  PatchData handle where control volumes
   *                                  are to be stored
   *  - patch_hierarchy (in):         PatchHierarchy on which to compute
   *                                  control volumes
   *
   * Return value:                    none
   *
   */
  static void computeControlVolumes(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    int control_volume_handle);

  /*!
   * copySAMRAIData() copies data from the PatchData associated
   * with the source handle to the PatchData associated with the
   * destination handle for the specified components.  It only copies
   * data in the interior of the Patch.
   *
   * Arguments:
   *  - patch_hierarchy (in):  PatchHierarchy on which to copy data
   *  - dst_handle (out):      PatchData handle for scratch space
   *  - src_handle (in):       PatchData handle for data to be copied to
   *                           scratch space
   *  - dst_component (in):    component of scratch space data to copy to
   *                           (default = 0)
   *  - src_component (in):    component of data to be copied to scratch
   *                           space
   *                           (default = 0)
   *
   * Return value:             none
   *
   */
  static void copySAMRAIData(
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    const int dst_handle,
    const int src_handle,
    const int dst_component = 0,
    const int src_component = 0);

  //! @}
protected:


  //! @{
  /*!
   ******************************************************************
   *
   * @name Helper methods
   *
   ******************************************************************/

  /*!
   * initializeComputeSpatialDerivativesParameters() sets up the
   * parameters required for computing spatial derivatives.
   *
   * Arguments:     none
   *
   * Return value:  none
   *
   */
  static void initializeComputeSpatialDerivativesParameters(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy);

  /*!
   * initializeComputeUnitNormalParameters() sets up the parameters
   * required for computing the unit normal vector.
   *
   * Arguments:     none
   *
   * Return value:  none
   *
   */

  static void initializeComputeUnitNormalParameters(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy);

  //! @}

  /******************************************************************
   *
   * Static Data Members
   *
   ******************************************************************/

  // Parameters for computing spatial derivatives
  // NOTE:  these are set up as needed.
  static int s_D1_one_ghostcell_handle;
  static int s_D1_two_ghostcells_handle;
  static int s_D2_two_ghostcells_handle;
  static int s_D1_three_ghostcells_handle;
  static int s_D2_three_ghostcells_handle;
  static int s_D3_three_ghostcells_handle;

  // Parameters for computing unit normal vector
  // NOTE:  these are set up as needed.
  static int s_compute_normal_grad_phi_handle;
  static int s_compute_normal_grad_phi_plus_handle;
  static int s_compute_normal_grad_phi_minus_handle;


private:

  /*
   * Private default constructor to prevent use.
   *
   * Arguments: none
   *
   */
  LevelSetMethodToolbox(){}

  /*
   * Private copy constructor to prevent use.
   *
   * Arguments:
   *  - rhs (in):  LevelSetMethodToolbox object to copy
   *
   */
  LevelSetMethodToolbox(const LevelSetMethodToolbox& rhs){}

  /*
   * Private assignment operator to prevent use.
   *
   * Arguments:
   *  - rhs (in):    LevelSetMethodToolbox to copy
   *
   * Return value:   *this
   *
   */
  const LevelSetMethodToolbox& operator=(
    const LevelSetMethodToolbox& rhs){
      return *this;
  }

};

} // end LSMLIB namespace

#endif
