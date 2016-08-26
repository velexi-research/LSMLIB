/*
 * File:        LevelSetFunctionIntegrator.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:     $09/19/2014$ jrigelo- pointers replaced by boost pointer: boost::shared_ptr
 * Description: Header file for level set method integrator class
 */

#ifndef included_LevelSetFunctionIntegrator_h
#define included_LevelSetFunctionIntegrator_h

/*! \class LSMLIB::LevelSetFunctionIntegrator
 *
 * \brief
 * The LevelSetFunctionIntegrator class provides the core functionality
 * required to integrate level set equation for evolving the level set
 * function that capture the dynamics of implicit surfaces and curves.
 * Support is provided for level set method calculations in one-, two-
 * and three-dimensions.
 *
 *
 * <h3> Numerical Discretization </h3>
 *
 * The LevelSetFunctionIntegrator supports the most of the standard
 * numerical schemes used by the level set method community.
 * To compute spatial derivatives, it provides the following
 * high-order finite-difference schemes: ENO1, ENO2, ENO3, and
 * WENO5.  For time-integration, total-variation-diminishing
 * Runge-Kutta schemes of order 1, 2, and 3 are available.
 * All of these options are easily configured via the input file
 * (see "User-specified parameters" section below).
 *
 * Flexibility is provided in determining the time step size used
 * to advance the level set functions.  The value of dt used in a
 * given TVD Runge-Kutta time step is one of the following time step
 * candidates:
 *
 *  - user_specified_dt
 *  - min{physics_dt, advection_dt, normal_vel_dt}
 *
 * A detailed description of each of these time step sizes is given
 * in the comments for the computeStableDt() method.
 *
 *
 * <h3> User-specified parameters (input database fields) </h3>
 *
 * <h4> Level Set Method Parameters: </h4>
 *
 * - start_time                  = start time for calculation (default = 0.0)
 * - end_time                    = end time for calculation (default = 10.0)
 * - cfl_number                  = CFL number (default = 0.5)
 * - spatial_derivative_type     = type of spatial derivative
 *                                 (default = "WENO")
 * - spatial_derivative_order    = order of spatial derivative (default = 5)
 * - tvd_runge_kutta_order       = order of Runge-Kutta time integration
 *                                 (default = 3)
 * - reinitialization_interval   = interval between reinitialization
 *                                 (default = 10)
 *                                 (reinitialization disabled if <= 0)
 * - reinitialization_stop_tol   = stopping criterion for termination of
 *                                 evolution of reinitialization equation.
 *                                 Reinitialization stops when the max norm
 *                                 of the change in the level set function
 *                                 drops below the specified tolerance.
 * - reinitialization_stop_dist  = approximate stopping criterion for
 *                                 reinitialization of level set functions.
 *                                 Reinitialization terminates after the
 *                                 information from the zero level set has
 *                                 propogated by approximately the specified
 *                                 distance.
 * - reinitialization_max_iters  = maximum number of time steps to take
 *                                 during the reinitialization process
 *                                 (default = 20)
 * - orthogonalization_interval  = interval between orthogonalizing phi
 *                                 and psi for codimension-two problems
 *                                 (default = 10)
 *                                 (orthogonalization disabled if <= 0)
 * - orthogonalization_stop_tol  = stopping criterion for termination of
 *                                 evolution of orthogonalization equation.
 *                                 Orthogonalization stops when the max norm
 *                                 of the change in the level set function
 *                                 drops below the specified tolerance.
 * - orthogonalization_stop_dist = approximate stopping criterion for
 *                                 orthogonalization of level set functions.
 *                                 Orthogonalization terminates after the
 *                                 information from the zero level set has
 *                                 propogated by approximately the specified
 *                                 distance.
 * - orthogonalization_max_iters = maximum number of time steps to take
 *                                 during the orthogonalization process
 *                                 (default = 20)
 *
 * <h4> Boundary Condition Parameters: </h4>
 *
 * - lower_bc_phi_[i]            = boundary conditions for the lower
 *                                 face of the computational domain in
 *                                 each coordinate direction for the i-th
 *                                 component of the PHI level set function.
 *                                 lower_bc_phi_[i] is a vector of length
 *                                 DIM.  The j-th entry should contain the
 *                                 type of boundary condition to impose at
 *                                 the lower boundary in the j-th coordinate
 *                                 direction.  See NOTES section for
 *                                 boundary condition types.  For
 *                                 information about the boundary
 *                                 condition types, see documentation
 *                                 of the BoundaryConditionModule class.
 *                                 (default = vector of zeros)
 * - upper_bc_phi_[i]            = boundary conditions for the upper
 *                                 face of the computational domain in
 *                                 each coordinate direction for the i-th
 *                                 component of the PHI level set function.
 *                                 upper_bc_phi_[i] is a vector of length
 *                                 DIM.  The j-th entry should contain the
 *                                 type of boundary condition to impose at
 *                                 the upper boundary in the j-th coordinate
 *                                 direction.  See NOTES section for
 *                                 boundary condition types.  For
 *                                 information about the boundary
 *                                 condition types, see documentation
 *                                 of the BoundaryConditionModule class.
 *                                 (default = vector of zeros)
 * - lower_bc_psi_[i]            = boundary conditions for the lower
 *                                 face of the computational domain in
 *                                 each coordinate direction for the i-th
 *                                 component of the PSI level set function.
 *                                 lower_bc_psi_[i] is a vector of length
 *                                 DIM.  The j-th entry should contain the
 *                                 type of boundary condition to impose at
 *                                 the lower boundary in the j-th coordinate
 *                                 direction.  See NOTES section for
 *                                 boundary condition types.  For
 *                                 information about the boundary
 *                                 condition types, see documentation
 *                                 of the BoundaryConditionModule class.
 *                                 (default = vector of zeros)
 * - upper_bc_psi_[i]            = boundary conditions for the upper
 *                                 face of the computational domain in
 *                                 each coordinate direction for the i-th
 *                                 component of the PSI level set function.
 *                                 upper_bc_psi_[i] is a vector of length
 *                                 DIM.  The j-th entry should contain the
 *                                 type of boundary condition to impose at
 *                                 the upper boundary in the j-th coordinate
 *                                 direction.  See NOTES section for
 *                                 boundary condition types.  For
 *                                 information about the boundary
 *                                 condition types, see documentation
 *                                 of the BoundaryConditionModule class.
 *                                 (default = vector of zeros)
 *
 * <h4> AMR Parameters: </h4>
 *
 * - use_AMR                     = TRUE if AMR should be used (default = FALSE)
 * - regrid_interval             = regridding interval (default = 5)
 * - tag_buffer_width            = number of buffer cells to use around
 *                                 cells tagged for refinement
 *                                 (default = 2)
 * - refinement_cutoff_value     = cutoff value for distance function
 *                                 (default = 1.0)
 *
 * <h4> Miscellaneous Parameters: </h4>
 *
 * - verbose_mode                = TRUE if status should be output during
 *                                 integration (default = FALSE)
 *
 * When restarting a computation, the following input parameters override
 * the values from the restart file:
 *    end_time,
 *    reinitialization_interval,
 *    reinitialization_stop_tol,
 *    reinitialization_stop_dist,
 *    reinitialization_max_iters,
 *    orthogonalization_interval,
 *    orthogonalization_stop_tol,
 *    orthogonalization_stop_dist,
 *    orthogonalization_max_iters,
 *    verbose_mode
 *
 *
 * <h3> NOTES: </h3>
 *
 *  - The precedence of the three input parameters for reinitialization
 *    and orthogonalization is as follows:
 *    -#  if provided, the stop tolerance is used with a maximum
 *        number of iterations determined from the calculation in (b)
 *        or a default value of 1000 iterations (to avoid failure to
 *        terminate)
 *    -#  when both the stop distance and maximum number of iterations
 *        are specified, the one that results in the smaller number of
 *        iterations is used.  when only one is specified, it is the
 *        sole parameter used to compute the maximum number of time steps
 *        to take when advancing the reinitialization/orthogonalization
 *        equation
 *    -#  when no stopping criteria are supplied, the maximum number of
 *        time steps taken defaults to 25.
 *
 *  - The reinitialization (and orthogonalization for codimension-two
 *    problems) process uses the same order TVD Runge-Kutta as specified
 *    for the time evolution of the level set equation(s).
 *
 *  - This class takes care of making sure that the scratch spaces
 *    for the level set functions have sufficient ghost cells to
 *    carry out the spatial derivative calculations.
 *
 *  - The spatial_derivative_type field in the input file is expected
 *    to be a string value, so LSMLIB_REAL quotes (") must be used; otherwise
 *    the InputManager will throw an error while parsing the input file.
 *
 *    - The numbering for the boundary conditions begins at 0.
 *
 *    - The boundary condition types are as follows:
 *      - NONE:                         0
 *      - HOMOEGENEOUS_NEUMANN:         1
 *      - LINEAR_EXTRAPOLATION:         2
 *      - SIGNED_LINEAR_EXTRAPOLATION:  3
 *      - ANTI_PERIODIC:                4
 *
 *      By default, no boundary conditions are imposed at any
 *      non-periodic boundary of the computational domain.
 *
 *    - Non-periodic boundary conditions MAY be imposed at
 *      boundaries that are specified to be periodic directions
 *      for the GridGeometry object associated with the
 *      PatchHierarchy set in the constructor.  The non-periodic
 *      boundary conditions will override the periodic boundary
 *      conditions.
 *
 *    - In order to use custom boundary conditions at a particular
 *      boundary location, the boundary condition type for that
 *      boundary location MUST be set to NONE.  Otherwise, the
 *      boundary condition will be overwritten by specified boundary
 *      condition type.
 *
 *  - AMR is currently UNAVAILABLE.  It is still in the development
 *    stages.
 *
 */

// Standard library headers
#include <ostream>
#include <string>

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Serializable.h"

// LSMLIB headers
#include "LSMLIB_config.h"
#include "LevelSetFunctionIntegratorStrategy.h"
#include "LevelSetMethodToolbox.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations
namespace LSMLIB { class BoundaryConditionModule; }
namespace LSMLIB { class OrthogonalizationAlgorithm; }
namespace LSMLIB { class ReinitializationAlgorithm; }
namespace SAMRAI { namespace geom { class CartesianGridGeometry; } }
namespace SAMRAI { namespace hier { class Box; } }
namespace SAMRAI { namespace hier { class Patch; } }
namespace SAMRAI { namespace hier { class PatchDataRestartManager; } }
namespace SAMRAI { namespace hier { class PatchHierarchy; } }
namespace SAMRAI { namespace hier { class PatchLevel; } }
namespace SAMRAI { namespace tbox { class Database; } }
namespace SAMRAI { namespace xfer { class RefineAlgorithm; } }
namespace SAMRAI { namespace xfer { class RefineSchedule; } }


namespace LSMLIB {

// explicit declaration of LevelSetMethodPatchStrategy and
// LevelSetMethodVelocityFieldStrategy to
// avoid circular dependencies in the header files.
class LevelSetMethodPatchStrategy;
class LevelSetMethodVelocityFieldStrategy;


/******************************************************************
 *
 * LevelSetFunctionIntegrator Class Definition
 *
 ******************************************************************/

class LevelSetFunctionIntegrator:
  public LevelSetFunctionIntegratorStrategy,
  public RefinePatchStrategy,
  public CoarsenPatchStrategy,
  public Serializable
{
public:

  //! @{
  /*!
   ****************************************************************
   *
   * @name Constructor and destructor
   *
   ****************************************************************/

  /*!
   * The constructor for LevelSetFunctionIntegrator sets up the
   * level set method and AMR parameters from the input and
   * restart databases, registers the level set variables with
   * the VariableDatabase, and sets up the objects involved in
   * the communication of data between patches in the PatchHierarchy.
   *
   * Arguments:
   *  - input_db (in):                      input database containing
   *                                        user-defined parameters for
   *                                        integration of level set equation
   *  - patch_hierarchy (in):               PatchHierarchy for computation
   *  - lsm_patch_strategy (in):            LevelSetMethodPatchStrategy
   *  - lsm_velocity_field_strategy (in):   LevelSetMethodVelocityFieldStrategy
   *  - num_level_set_fcn_components (in):  number of components of level set
   *                                        functions (default = 1)
   *  - codimension (in):                   codimension of problem
   *                                        (default = 1)
   *  - object_name (in):                   string name for object (default =
   *                                        "LevelSetFunctionIntegrator")
   *
   */
  LevelSetFunctionIntegrator(
    const tbox::Dimension& dim,
    boost::shared_ptr<Database> input_db,
    boost::shared_ptr< PatchHierarchy > patch_hierarchy,
    LevelSetMethodPatchStrategy* lsm_patch_strategy,
    LevelSetMethodVelocityFieldStrategy* lsm_velocity_field_strategy,
    const int num_level_set_fcn_components = 1,
    const int codimension = 1,
    const string& object_name = "LevelSetFunctionIntegrator");

  /*!
   * The destructor for LevelSetFunctionIntegrator unregisters
   * the LevelSetFunctionIntegrator object a restart item.
   */
  virtual ~LevelSetFunctionIntegrator();

  //! @}


  //! @{

  /****************************************************************
   *
   * Overridden Methods from the LevelSetFunctionIntegratorStrategy
   * class
   *
   ****************************************************************/

  /*!
   ********************************************************************
   *
   * @name Accessor methods for level set functions and related data
   *
   ********************************************************************/

  /*!
   * getPhiPatchDataHandle() returns the patch data handle for phi.
   *
   * Arguments:     none
   *
   * Return value:  PatchData handle for phi
   *
   * NOTES:
   *  - The PatchData for phi associated with the returned PatchData
   *    handle has the number of ghostcells required for the spatial
   *    derivative type and order specified when the
   *    LevelSetFunctionIntegrator object was constructed.
   *
   */
  virtual int getPhiPatchDataHandle() const;

  /*!
   * getPsiPatchDataHandle() returns the patch data handle for phi.
   *
   * Arguments:     none
   *
   * Return value:  PatchData handle for phi
   *
   * NOTES:
   *  - The PatchData for psi associated with the returned PatchData
   *    handle has the number of ghostcells required for the spatial
   *    derivative type and order specified when the
   *    LevelSetFunctionIntegrator object was constructed.
   *
   */
  virtual int getPsiPatchDataHandle() const;

  /*!
   * getControlVolumePatchDataHandle() returns the patch data handle for
   * the control volume.
   *
   * Arguments:     none
   *
   * Return value:  PatchData handle for control volume
   *
   */
  virtual int getControlVolumePatchDataHandle() const;

  //! @}


  //! @{
  /*!
   *******************************************************************
   *
   *  @name Accessor methods for simulation state
   *
   *******************************************************************/

  /*!
   * getStartTime() returns the simulation start time.
   *
   * Arguments:      none
   *
   * Return value :  start time
   *
   */
  virtual LSMLIB_REAL getStartTime() const;

  /*!
   * getEndTime() returns the simulation end time.
   *
   * Arguments:      none
   *
   * Return value :  end time
   *
   */
  virtual LSMLIB_REAL getEndTime() const;

  /*!
   * getCurrentTime() returns the current simulation time.
   *
   * Arguments:      none
   *
   * Return value :  current time
   *
   */
  virtual LSMLIB_REAL getCurrentTime() const;

  /*!
   * endTimeReached() returns true if the end time for the
   * integration has been reached; otherwise, it returns false.
   *
   * Arguments:     none
   *
   * Return value:  true if the current time is equal to or after
   *            the end time for the integration; false otherwise
   *
   */
  virtual bool endTimeReached() const;

  /*!
   * numIntegrationStepsTaken() returns the number of integration steps
   * that have been taken during the level set method calculation.
   *
   * Arguments:      none
   *
   * Return value :  current integration step count
   *
   */
  virtual int numIntegrationStepsTaken() const;

  /*!
   * printClassData() prints the values of the data members for
   * an instance of the LevelSetFunctionIntegrator class.
   *
   * Arguments:
   *  - os (in):     output stream to write object information
   *
   * Return value:   none
   *
   */
  virtual void printClassData(ostream& os) const;

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @ Accessor methods for numerical parameters
   *
   ****************************************************************/

  /*!
   * getSpatialDerivativeType() returns type of spatial derivative
   * discretization used in the level set method calculation.
   *
   * Arguments:      none
   *
   * Return value :  spatial derivative type
   *
   * NOTES:
   *  - Meaning of return values: 0 = ENO, 1 = WENO.
   */
  virtual int getSpatialDerivativeType() const;

  /*!
   * getSpatialDerivativeOrder() returns order of the spatial derivative
   * discretization used in the level set method calculation.
   *
   * Arguments:      none
   *
   * Return value :  spatial derivative order
   *
   */
  virtual int getSpatialDerivativeOrder() const;

  /*!
   * getTVDRungeKuttaOrder() returns order of the TVD Runge-Kutta
   * integration scheme used in the level set method calculation.
   *
   * Arguments:      none
   *
   * Return value :  TVD Runge-Kutta order
   *
   */
  virtual int getTVDRungeKuttaOrder() const;

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for time advancing level set functions
   *
   ****************************************************************/

  /*!
   * computeStableDt() computes the maximum allowable dt for the
   * next time step of the level set functions.  The maximum stable
   * dt is computed using the algorithm:
   *
   * <pre>
   *   if (user_specified_dt < LSMLIB_REAL_MAX)
   *
   *     dt = user_specified_dt
   *
   *   else
   *
   *     dt = min{physics_dt, advection_dt, normal_vel_dt}
   * </pre>
   *
   *
   * In this algorithm, the various dt values are as follows:
   *
   *  - user_specified_dt: arbitrarily specified by the user and computed
   *                       on a patch-by-patch basis by the user-defined
   *                       concrete subclass of LevelSetMethodPatchStrategy
   *
   *  - physics_dt:        based on the physics that determines the velocity
   *                       field and computed by the user-defined concrete
   *                       subclass of LevelSetMethodVelocityFieldStrategy
   *
   *  - advection_dt:      computed (in this method if required) using
   *                       the current external (vector) velocity field
   *                       and the "CFL time step" equation:
   *
   *                dt * max{|u|/dx + |v|/dy + |w|/dz} = cfl_number
   *
   *                       See Osher & Fedkiw, "Level Set Methods and Dynamic
   *                       Implicit Surfaces"
   *
   *  - normal_vel_dt:     computed (in this method if required) using
   *                       current normal (scalar) velocity field and
   *                       the "CFL time step" equation:
   *
   *     dt * max{|vel_n|*(|phi_x|/dx + |phi_y|/dy + |phi_z|/dz)} = cfl_number
   *
   *                       See Osher & Fedkiw, "Level Set Methods and Dynamic
   *                       Implicit Surfaces"
   *
   * Arguments:     none
   *
   * Return value:  maximum stable dt for next time step
   *
   * NOTES:
   *  - advection_dt is only used if the LevelSetMethodVelocityFieldStrategy
   *    provides an external (vector) velocity field.
   *  - normal_vel_dt is only used if the LevelSetMethodVelocityFieldStrategy
   *    provides a normal (scalar) velocity field.
   *  - if the value of physics_dt or user_specified_dt is negative, then
   *    is it ignored.
   *
   */
  virtual LSMLIB_REAL computeStableDt();

  /*!
   * advanceLevelSetFunctions() advances the level set function
   * phi (and psi for codimension-two problems) by the specified
   * time increment, dt.
   *
   * Arguments:
   *  - dt (in):     time increment to advance the level set functions
   *
   * Return value:   true if patch hierarchy needs to be regridded after
   *                 this time step; false otherwise.
   *
   * NOTES:
   *  - AMR is NOT yet implemented, so advanceLevelSetFunctions()
   *    currently always returns false.
   *
   */
  virtual bool advanceLevelSetFunctions(const LSMLIB_REAL dt);

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Reinitialization and orthogonalization methods
   *
   ****************************************************************/

  /*!
   * getReinitializationInterval() should return the number of time
   * steps between reinitializations of the level set functions.
   *
   * Arguments:     none
   *
   * Return value:  reinitialization interval
   *
   */
  virtual int getReinitializationInterval() const;

  /*!
   * setReinitializationInterval() should set the number of time
   * steps between reinitializations of the level set functions.
   *
   * Arguments:
   *  - interval (in):  number of time steps to take between
   *                    reinitializations of the level set functions
   *
   * Return value:      none
   *
   * NOTES:
   *  - To disable reinitialization, set interval to zero.
   *
   */
  virtual void setReinitializationInterval(const int interval);

  /*!
   * reinitializeLevelSetFuntions() reinitializes the level set
   * functions to be distance functions using the reinitialization
   * equation:
   *
   *   phi_t + sgn(phi) * ( |grad(phi)| - 1 ) = 0
   *
   * This Hamilton-Jacobi equation is advanced in time towards steady-state
   * using the same TVD Runge-Kutta method selected to advance the level
   * set equation.  The number of steps taken is a function of the
   * user-specified input parameters: reinitialization_stop_dist
   * and reinitialization_max_iters.  grad(phi) is computed using the
   * same numerical discretization used to calculate the spatial
   * derivatives when advancing the level set equation.  A Godunov scheme
   * is used to select whether the appropriate spatial derivative
   * approximation for each component of grad(phi).
   *
   * Arguments:
   *  - level_set_fcn (in):   level set function (i.e. PHI or PSI) to
   *                          reinitialize (default = PHI)
   *  - max_iterations (in):  maximum number of iterations to use
   *                          for reinitialization.  Set max_iterations
   *                          to -1 to use the value specified in the
   *                          input file.
   *                          (default = -1)
   *
   * Return value:            none
   *
   * NOTES:
   *  - If max_iterations is set to a non-negative value, it overrides
   *    ALL of the stopping criteria specified in the input file.
   *
   */
  virtual void reinitializeLevelSetFunctions(
    const LEVEL_SET_FCN_TYPE level_set_fcn = LSMLIB::PHI,
    const int max_iterations = -1);

  /*!
   * getOrthogonalizationInterval() should return the number of time
   * steps between orthogonalizations of the level set functions.
   *
   * Arguments:     none
   *
   * Return value:  orthogonalization interval
   *
   */
  virtual int getOrthogonalizationInterval() const;

  /*!
   * setOrthogonalizationInterval() should set the number of time
   * steps between orthogonalizations of the level set functions.
   *
   * Arguments:
   *  - interval (in):  number of time steps to take between
   *                    orthogonalizations of the level set functions
   *
   * Return value:      none
   *
   * NOTES:
   *  - To disable orthogonalization, set interval to zero.
   *
   */
  virtual void setOrthogonalizationInterval(const int interval);

  /*!
   * orthogonalizeLevelSetFunctions() reinitializes the level set functions
   * phi and psi and evolves them so that they have orthogonal gradients.
   * This goal is achieved by solving the orthogonalization equation:
   *
   *   phi_t + sgn(psi) * ( grad(psi)/|grad(psi)| ) dot grad(phi) = 0
   *
   * This Hamilton-Jacobi equation is advanced in time towards steady-state
   * using the same TVD Runge-Kutta method selected to advance the level
   * set equation.  The number of steps taken is a function of the
   * user-specified input parameters: orthogonalization_stop_dist
   * and orthogonalization_max_iters.  grad(phi) is computed using the
   * same numerical discretization used to calculate the spatial
   * derivatives when advancing the level set equation.  grad(psi) is
   * computed using by taking the average of the forward and backward
   * spatial derivatives.  grad(phi) is computed via a simple upwinding
   * scheme that treats grad(psi) as the velocity.
   *
   * Arguments:
   *  - level_set_fcn (in):          level set function (i.e. PHI or PSI) to
   *                                 evolve to orthogonalize grad(phi) and
   *                                 grad(psi)
   *  - max_reinit_iterations (in):  maximum number of iterations to use
   *                                 for reinitialization.  Set max_iterations
   *                                 to -1 to use the value specified in the
   *                                 input file.
   *                                 (default = -1)
   *  - max_ortho_iterations (in):   maximum number of iterations to use
   *                                 for orthogonalization.  Set max_iterations
   *                                 to -1 to use the value specified in the
   *                                 input file.
   *                                 (default = -1)
   *
   * Return value:                   none
   *
   * NOTES:
   *  - This method is ONLY used for codimension-two problems.  For
   *    codimension-one problems, this method is never invoked.
   *
   *  - If max_reinit_iterations is set to a non-negative value, it overrides
   *    ALL of the reinitialization stopping criteria specified in the input
   *    file.
   *
   *  - If max_ortho_iterations is set to a non-negative value, it overrides
   *    ALL of the orthogonalization stopping criteria specified in the input
   *    file.
   *
   */
  virtual void orthogonalizeLevelSetFunctions(
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int max_reinit_iterations = -1,
    const int max_ortho_iterations = -1);

  //! @}


  //! @{
  /*!
   ********************************************************************
   *
   * @name Methods for pre-/post-processing level set function data
   *
   * These methods are used to prepare level set function data for
   * computation of the velocity field.
   *
   ********************************************************************/

  /*!
   * preprocessInitializeVelocityField() pre-processes level set
   * function data for use in initializing velocity field data.
   *
   * Arguments:
   *  - phi_handle (out):   PatchData handle for phi
   *  - psi_handle (out):   PatchData handle for psi
   *  - hierarchy (in):     PatchHierarchy containing level to
   *                        pre-process
   *  - level_number (in):  level number for level to be
   *                        pre-processed
   *
   * Return value:          none
   *
   */
  virtual void preprocessInitializeVelocityField(
    int& phi_handle,
    int& psi_handle,
    const boost::shared_ptr< PatchHierarchy > hierarchy,
    const int level_number);

  /*!
   * postprocessInitializeVelocityField() post-processes level set
   * function data after initializing velocity field data.
   *
   * Arguments:
   *  - hierarchy (in):     PatchHierarchy containing level to
   *                        post-process
   *  - level_number (in):  level number for level to be
   *                        post-processed
   *
   * Return value:          none
   *
   */
  virtual void postprocessInitializeVelocityField(
    const boost::shared_ptr< PatchHierarchy > hierarchy,
    const int level_number);

  //! @}


  //! @{
  /*!
   **************************************************************************
   *
   * @name Utility methods for managing state of LevelSetFunctionIntegrator
   *
   **************************************************************************/

  /*!
   * putToRestart() writes the state of the LevelSetFunctionIntegrator
   * object to the given database (typically for restart purposes).
   *
   * Arguments:
   *  - db (in):     database in which to write current state of
   *                 LevelSetFunctionIntegrator object
   *
   * Return value:   none
   *
   */
  virtual void putToRestart(const boost::shared_ptr<Database>& restart_db) const;

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for setting initial and boundary conditions
   *
   * NOTES:
   *  - initializeLevelData() overrides virtual methods
   *    declared in the StandardTagAndInitStrategy base class.
   *  - setPhysicalBoundaryConditions() overrides the pure virtual
   *    method from the RefinePatchStrategy base class.
   *
   ****************************************************************/

  /*!
   * initializeLevelData() allocates and initializes the level
   * set function(s) for a new patch level in the patch hierarchy.
   *
   * Arguments:
   *  - hierarchy (in):       BasePatchHierarchy on which to tag cells for
   *                          refinement
   *  - level_number (in):    BasePatchLevel number on which to tag cells for
   *                          refinement
   *  - can_be_refined (in):  true if this is NOT the finest level in the
   *                          BasePatchHierarchy; false it if is
   *  - init_data_time (in):  true if the BasePatchLevel is being introduced
   *                          for the first time; false otherwise
   *  - old_level (in):       old BasePatchLevel from which data for new
   *                          BasePatchLevel should be taken
   *  - allocate_data (in):   true if PatchData needs to be allocated before
   *                          it is initialized; false otherwise
   *
   * Return value:            none
   *
   */
  virtual void initializeLevelData (
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const boost::shared_ptr<hier::PatchLevel>& old_level =
         boost::shared_ptr<hier::PatchLevel>(),
      const bool allocate_data = true );

  /*!
   * setPhysicalBoundaryConditions() sets the data in ghost
   * cells corresponding to custom physical boundary conditions
   * implemented by the user in the method --
   * LevelSetMethodPatchStrategy::setLevelSetFunctionBoundaryConditions()
   *
   * Arguments:
   *  - patch (in):                Patch on which to set boundary conditions
   *  - fill_time (in):            time at which boundary conditions are to
   *                               be set
   *  - ghost_width_to_fill (in):  width of the ghostcells to fill
   *
   * Return value:                 none
   *
   */
  virtual void setPhysicalBoundaryConditions(
    Patch& patch,
    const double fill_time,
    const IntVector & ghost_width_to_fill);

  /*!
   * setBoundaryConditions() sets the boundary conditions to impose
   * at the outer faces of the computational domain for the
   * specified components of the vector level set function.  If no
   * component is specified, then ALL components of a vector level
   * set function will be set to have the same boundary conditions
   * imposed.
   *
   * Arguments:
   *  - lower_bc (in):       vector of integers specifying the
   *                         type of boundary conditions to impose
   *                         on the lower face of the computational
   *                         domain in each coordinate direction.
   *                         The i-th entry should contain the type
   *                         of boundary condition to impose at the lower
   *                         boundary in the i-th coordinate direction.
   *                         See NOTES for boundary condition
   *                         types.  For information about the
   *                         boundary condition types, see documentation
   *                         of the BoundaryConditionModule class.
   *  - upper_bc (in):       vector of integers specifying the
   *                         type of boundary conditions to impose
   *                         on the upper face of the computational
   *                         domain in each coordinate direction.
   *                         The i-th entry should contain the type
   *                         of boundary condition to impose at the upper
   *                         boundary in the i-th coordinate direction.
   *                         See NOTES for boundary condition
   *                         types.  For information about the
   *                         boundary condition types, see documentation
   *                         of the BoundaryConditionModule class.
   *  - level_set_fcn (in):  level set function (i.e. PHI or PSI) for
   *                         which to set boundary conditions
   *  - component (in):      component of level set function for
   *                         which to set boundary conditions
   *                         (default = -1)
   *
   * Return value:           none
   *
   * NOTES:
   *  - The boundary condition types are as follows:
   *    - NONE:                         0
   *    - HOMOEGENEOUS_NEUMANN:         1
   *    - LINEAR_EXTRAPOLATION:         2
   *    - SIGNED_LINEAR_EXTRAPOLATION:  3
   *    - ANTI_PERIODIC:                4
   *
   *    If anti-periodic boundary conditions are imposed in the i-th
   *    coordinate direction, then the i-th entry of lower_bc
   *    and upper_bc MUST both be set to ANTI_PERIODIC.
   *
   *    For more information about the various boundary conditions,
   *    see the documentation for the BoundaryConditionModule class.
   *
   *  - Anti-periodic boundary conditions are only imposed for those
   *    directions that are specified by bc AND that are periodic
   *    directions for the GridGeometry object associated with the
   *    PatchHierarchy set in the constructor.  If a direction is
   *    specified to be anti-periodic by the bc variable but
   *    is not a periodic direction for the GridGeometry object, then
   *    that direction is NOT treated as an anti-periodic direction.
   *
   *  - If component is set to a negative number, than ALL components of
   *    the level set function will be set to have the specified
   *    homogeneous Neumann boundaries.
   *
   */
  virtual void setBoundaryConditions(
    const IntVector& lower_bc,
    const IntVector& upper_bc,
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int component = -1);

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for managing grid configuration
   *
   * NOTES:
   *  - applyGradientDetector() and resetHierarchyConfiguration()
   *    override virtual methods declared in the
   *    StandardTagAndInitStrategy base class.
   *
   ****************************************************************/

  /*!
   * Set integer tags to "1" in cells where refinement of the given
   * level should occur according to the criteria that the absolute
   * value of the level set functions, phi and psi, is less than some
   * user-supplied threshold.
   *
   * Arguments:
   *  - hierarchy (in):       BasePatchHierarchy on which to tag cells for
   *                          refinement
   *  - level_number (in):    BasePatchLevel number on which to tag cells for
   *                          refinement
   *  - error_data_time (in): ignored by LevelSetFunctionIntegrator class
   *  - tag_index (in):       PatchData index of the cell-centered integer
   *                          tag data
   *  - initial_time (in):    ignored by LevelSetFunctionIntegrator class
   *  - uses_richardson_extrapolation_too (in):
   *                          ignored by LevelSetFunctionIntegrator class
   *
   * Return value:            none
   *
   */
  virtual void applyGradientDetector(
      const boost::shared_ptr< PatchHierarchy > hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

  /*!
   * resetHierarchyConfiguration() resets any internal
   * hierarchy-dependent information.
   *
   * Arguments:
   *  - new_hierarchy (in):   BasePatchHierarchy to configure
   *  - coarsest_level (in):  level number of coarsest PatchLevel to reset
   *  - finest_level (in):    level number of finest BasePatchLevel to reset
   *
   * Return value:            none
   *
   */
  virtual void resetHierarchyConfiguration(
     const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
     const int coarsest_level,
     const int finest_level );

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @ name Inherited Methods
   *
   * Additional Declarations/Implementations of Pure Virtual
   * Methods from the RefinePatchStrategy, CoarsenPatchStrategy,
   * and Base Classes.
   *
   * NOTE: understanding these methods is NOT critical for the
   *       typical use of the LevelSetFunctionIntegrator class.
   *
   ****************************************************************/

  /*!
   * getRefineOpStencilWidth() returns the maximum stencil width needed
   * for user-defined data interpolation operations.  Default is to
   * return zero, assuming no user-defined operations provided.
   *
   * Arguments:     none
   *
   * Return value:  (0,0,0)
   *
   */
  virtual IntVector getRefineOpStencilWidth(const tbox::Dimension& dim) const;

  /*!
   * preprocessRefine() does no special user-defined spatial
   * interpolation.
   *
   * Arguments:
   *  - fine (in):      ignored by LevelSetFunctionIntegrator class
   *  - coarse (in):    ignored by LevelSetFunctionIntegrator class
   *  - fine_box (in):  ignored by LevelSetFunctionIntegrator class
   *  - ratio (in):     ignored by LevelSetFunctionIntegrator class
   *
   * Return value:      none
   *
   */
  virtual void preprocessRefine(
    Patch& fine,
    const Patch& coarse,
    const Box& fine_box,
    const IntVector& ratio);

  /*!
   * postprocessRefine() does no special user-defined spatial
   * interpolation.
   *
   * Arguments:
   *  - fine (in):      ignored by LevelSetFunctionIntegrator class
   *  - coarse (in):    ignored by LevelSetFunctionIntegrator class
   *  - fine_box (in):  ignored by LevelSetFunctionIntegrator class
   *  - ratio (in):     ignored by LevelSetFunctionIntegrator class
   *
   * Return value:      none
   *
   */
  virtual void postprocessRefine(
    Patch& fine,
    const Patch& coarse,
    const Box& fine_box,
    const IntVector& ratio);

  /*!
   * getCoarsenOpStencilWidth() returns the maximum stencil width
   * needed for user-defined data coarsen operations.  Default is
   * to return zero, assuming no user-defined operations provided.
   *
   * Arguments:     none
   *
   * Return value:  (0,0,0)
   *
   */
  virtual IntVector getCoarsenOpStencilWidth(const tbox::Dimension& dim) const;

  /*!
   * preprocessCoarsen() does no special user-defined spatial
   * coarsening.
   *
   * Arguments:
   *  - coarse (in):      ignored by LevelSetFunctionIntegrator class
   *  - fine (in):        ignored by LevelSetFunctionIntegrator class
   *  - coarse_box (in):  ignored by LevelSetFunctionIntegrator class
   *  - ratio (in):       ignored by LevelSetFunctionIntegrator class
   *
   * Return value:        none
   *
   */
  virtual void preprocessCoarsen(
    Patch& coarse,
    const Patch& fine,
    const Box& coarse_box,
    const IntVector& ratio);

  /*!
   * preprocessCoarsen() does no special user-defined spatial
   * coarsening.
   *
   * Arguments:
   *  - coarse (in):      ignored by LevelSetFunctionIntegrator class
   *  - fine (in):        ignored by LevelSetFunctionIntegrator class
   *  - coarse_box (in):  ignored by LevelSetFunctionIntegrator class
   *  - ratio (in):       ignored by LevelSetFunctionIntegrator class
   *
   * Return value:        none
   *
   */
  virtual void postprocessCoarsen(
    Patch& coarse,
    const Patch& fine,
    const Box& coarse_box,
    const IntVector& ratio);

  //! @}


protected:

  //! @{
  /*!
   ****************************************************************
   *
   * @name Main stages of time-advance
   *
   ****************************************************************/

  /*!
   * advanceLevelSetEqnUsingTVDRK*() advances the level set functions
   * using the level set equation using a first-, second-, or third-order
   * TVD Runge-Kutta step.
   *
   * Arguments:
   *  - dt (in):         time increment to advance the level set functions
   *
   * Return value:       none
   *
   */
  virtual void advanceLevelSetEqnUsingTVDRK1(
    const LSMLIB_REAL dt);
  virtual void advanceLevelSetEqnUsingTVDRK2(
    const LSMLIB_REAL dt);
  virtual void advanceLevelSetEqnUsingTVDRK3(
    const LSMLIB_REAL dt);

  /*!
   * computeLevelSetEquationRHS() computes the right-hand side of
   * the level set equation when it is written in the form:
   *
   *   phi_t = ...
   *
   * Arguments:
   *  - level_set_fcn (in):  level set function to compute RHS
   *                         of evolution equation (i.e. PHI or PSI)
   *  - phi_handle (in):     PatchData handle for phi that should
   *                         be used to compute spatial derivatives
   *  - component (in):      component of level set function that for
   *                         which the RHS is being computed
   *                         (default = 0)
   *
   *
   * Return value:           none
   *
   */
  virtual void computeLevelSetEquationRHS(
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int phi_handle,
    const int component = 0);

  /*!
   * addAdvectionTermToLevelSetEquationRHS() adds the contribution of
   * an advection term to right-hand side of the level set equation when
   * it is written in the form:
   *
   *   phi_t = - velocity dot grad(phi) + ...
   *
   * Arguments:
   *  - level_set_fcn (in):  level set function to compute RHS
   *                         of evolution equation (i.e. PHI or PSI)
   *  - phi_handle (in):     PatchData handle for phi that should
   *                         be used to compute spatial derivatives
   *  - component (in):      component of level set function that for
   *                         which the RHS is being computed
   *                         (default = 0)
   *
   * Return value:           none
   *
   * NOTES:
   *  - The phi_handle that is passed into this method is for scratch
   *    data which only has one component (even for vector level set
   *    calculations).  The component argument is only used to determine
   *    which velocity field to use.
   *
   */
  virtual void addAdvectionTermToLevelSetEquationRHS(
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int phi_handle,
    const int component = 0);

  /*!
   * addNormalVelocityTermToLevelSetEquationRHS() adds the contribution
   * of a normal velocity term to right-hand side of the level set
   * equation when it is written in the form:
   *
   *   phi_t = - vel_n |grad(phi)| + ...
   *
   * Arguments:
   *  - level_set_fcn (in):  level set function to compute RHS
   *                         of evolution equation (i.e. PHI or PSI)
   *  - phi_handle (in):     PatchData handle for phi that should
   *                         be used to compute spatial derivatives
   *  - component (in):      component of level set function that for
   *                         which the RHS is being computed
   *                         (default = 0)
   *
   * Return value:           none
   *
   * NOTES:
   *  - The phi_handle that is passed into this method is for scratch
   *    data which only has one component (even for vector level set
   *    calculations).  The component argument is only used to determine
   *    which velocity field to use.
   *
   */
  virtual void addNormalVelocityTermToLevelSetEquationRHS(
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int phi_handle,
    const int component = 0);

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Utility methods
   *
   ****************************************************************/

  /*!
   * initializeVariables() sets up the variables and PatchData handles
   * for the level set method.
   *
   * Arguments:     none
   *
   * Return value:  none
   *
   */
  virtual void initializeVariables();

  /*!
   * initializeCommunicationObjects() initializes the objects
   * involved in communication between patches (including
   * inter-node data transfers).
   *
   * Arguments:     none
   *
   * Return value:  none
   *
   */
  virtual void initializeCommunicationObjects();

  /*!
   * getFromInput() configures the LevelSetFunctionIntegrator object
   * from the values in the specified input database.
   *
   * Arguments:
   *  - db (in):               input database from which to take parameter
   *                           values
   *  - is_from_restart (in):  true if object is to be configured using
   *                           parameter values in a restart file;
   *                           false otherwise
   *
   * Return value:             none
   *
   * NOTES:
   *  - An assertion results if the database boost pointer is null.
   *
   */
  virtual void getFromInput(boost::shared_ptr<Database> db,
                            bool is_from_restart);

  /*!
   * getFromRestart() configures the LevelSetFunctionIntegrator object
   * based on parameter values from the restart file.
   *
   * Arguments:      none
   *
   * Return values:  none
   *
   * NOTES:
   *  - An assertion results if the database boost pointer is null.
   *
   */
  virtual void getFromRestart();

  //! @}

  /****************************************************************
   *
   * Data Members
   *
   ****************************************************************/

  /*
   * The object name is used for error/warning reporting and also as a
   * string label for restart database entries.
   */
  string d_object_name;


  /*
   * User-defined parameters
   */

  // general level set method parameters
  int d_num_level_set_fcn_components;   // number of components for level
                                        //   set functions
  int d_codimension;                    // codimension of problem
  LSMLIB_REAL d_start_time;                  // start time for calculation
  LSMLIB_REAL d_end_time;                    // end time for calculation
  LSMLIB_REAL d_cfl_number;                  // CFL number
  SPATIAL_DERIVATIVE_TYPE               // type of spatial derivative:
       d_spatial_derivative_type;       //   default = WENO
  int d_spatial_derivative_order;       // order of spatial derivative
  int d_tvd_runge_kutta_order;          // order of TVD Runge-Kutta time
                                        //   integration
  int d_reinitialization_interval;      // interval between reinitialization
  LSMLIB_REAL d_reinitialization_stop_tol;   // stopping criterion for termination
                                        //   of evolution of reinitialization
                                        //   equation.  Reinitialization
                                        //   stops when the max norm of the
                                        //   change in the level set function
                                        //   drops below the specified
                                        //   tolerance.
  LSMLIB_REAL d_reinitialization_stop_dist;  // stopping criterion for termination
                                        //   of evolution of reinitialization
                                        //   equation.  Reinitialization
                                        //   approximately propogates
                                        //   information outward from the
                                        //   zero level set by the specified
                                        //   distance
  int d_reinitialization_max_iters;     // maximum number of time steps
                                        //   for reinitialization iteration
  int d_orthogonalization_interval;     // interval between orthogonalizing
                                        //   phi and psi for codim-two problems
  LSMLIB_REAL d_orthogonalization_stop_tol;  // stopping criterion for termination
                                        //   of evolution of orthogonalization
                                        //   equation.  Orthogonalization
                                        //   stops when the max norm of the
                                        //   change in the level set function
                                        //   drops below the specified
                                        //   tolerance.
  LSMLIB_REAL d_orthogonalization_stop_dist; // stopping criterion for termination
                                        //   of evolution of orthogonalization
                                        //   equation.  Orthogonalization
                                        //   approximately propogates
                                        //   information outward from the
                                        //   zero level set by the specified
                                        //   distance
  int d_orthogonalization_max_iters;    // maximum number of time steps
                                        //   for orthogonalization iteration

  // AMR parameters
  bool d_use_AMR;                       // true if AMR should be used
  int d_regrid_interval;                // regridding interval
  int d_tag_buffer_width;               // number of buffer cells to use around
                                        //   cells tagged for refinement
  LSMLIB_REAL d_refinement_cutoff_value;     // cutoff value for distance function

  // Miscellaneous parameters
  bool d_verbose_mode;                  // true if status information should
                                        //   be output

  /*
   * User-defined level set method strategy objects
   */

  // Boost pointer to the LevelSetMethodPatchStrategy.  This object
  // is used to initialize and set boundary conditions for the
  // level set functions.
  LevelSetMethodPatchStrategy*  d_lsm_patch_strategy;

  // Boost pointer to LevelSetMethodVelocityFieldStrategy.  This object
  // is used to set the velocity field for the time advance of the
  // level set functions.  It is also used to provide some physics-based
  // restrictions on the maximum allowable dt to use for advancing
  // the level set functions in time.
  LevelSetMethodVelocityFieldStrategy*  d_lsm_velocity_field_strategy;


  /*
   * Grid management objects
   */

  // Boost pointer to the patch hierarchy object
  boost::shared_ptr< PatchHierarchy > d_patch_hierarchy;

  // Boost pointer to the grid geometry.  The CartesianGridGeometry object
  // is used to set up initial data, set physical boundary conditions,
  // and register plot variables.
  boost::shared_ptr< CartesianGridGeometry > d_grid_geometry;

  // Boost pointer to reinitialization algorithms which manage the
  // reinitialization of level set functions to distance functions
  boost::shared_ptr< ReinitializationAlgorithm > d_phi_reinitialization_alg;
  boost::shared_ptr< ReinitializationAlgorithm > d_psi_reinitialization_alg;

  // Boost pointer to orthogonalization algorithm which manages the
  // orthogonalization of gradients of phi and psi for
  // codimension-two problems
  boost::shared_ptr< OrthogonalizationAlgorithm > d_orthogonalization_alg;

  /*
   * PatchData handles for data required by level set method
   */

  // PatchData handles for level set functions
  // NOTE:  the 0-th component of both of these vectors
  //        contains the PatchData handle for the current
  //        solution of the level set function.  The
  //        remainder are PatchData handles for scratch
  //        space associated with the TVD Runge-Kutta
  //        time integration.
  vector<int> d_phi_handles;
  vector<int> d_psi_handles;

  // forward and backward spatial derivatives
  int d_grad_phi_plus_handle;
  int d_grad_psi_plus_handle;
  int d_grad_phi_minus_handle;
  int d_grad_psi_minus_handle;

  // upwind spatial derivatives
  int d_grad_phi_upwind_handle;
  int d_grad_psi_upwind_handle;

  // right-hand side of evolution equations
  // (i.e. level set equation, reinitialization
  //  equation, and orthogonalization equation)
  int d_rhs_phi_handle;
  int d_rhs_psi_handle;

  // auxilliary variables
  int d_control_volume_handle;

  // level set ghostcell width
  IntVector d_level_set_ghostcell_width;

  /*
   * Problem dimension.
  */
  const tbox::Dimension d_dim;

  /*
   * Component selectors to organize variables into logical groups
   */
  ComponentSelector d_solution_variables;
  ComponentSelector d_time_advance_scratch_variables;
  ComponentSelector d_compute_stable_dt_scratch_variables;
  ComponentSelector d_reinitialization_scratch_variables;
  ComponentSelector d_orthogonalization_scratch_variables;
  ComponentSelector d_persistent_variables;


  /* internal state  variables */

  // flags
  bool d_use_reinitialization;
  bool d_use_reinitialization_stop_tol;
  bool d_use_reinitialization_stop_dist;
  bool d_use_reinitialization_max_iters;
  bool d_use_orthogonalization;
  bool d_use_orthogonalization_stop_tol;
  bool d_use_orthogonalization_stop_dist;
  bool d_use_orthogonalization_max_iters;

  // counter variables
  LSMLIB_REAL d_current_time;
  int d_num_integration_steps_taken;
  int d_reinitialization_count;
  int d_orthogonalization_count;
  LEVEL_SET_FCN_TYPE d_orthogonalization_evolved_field;
  int d_regrid_count;

  /*
   * Boundary condition objects
   */
  boost::shared_ptr< BoundaryConditionModule > d_bc_module;
  Array< IntVector > d_lower_bc_phi;
  Array< IntVector > d_upper_bc_phi;
  Array< IntVector > d_lower_bc_psi;
  Array< IntVector > d_upper_bc_psi;

  /*
   * Communication objects.
   */

  // for filling a new level
  boost::shared_ptr< RefineAlgorithm > d_fill_new_level;

  // for filling boundary data used during the calculation of
  // a stable dt for motion under normal velocity
  boost::shared_ptr< RefineAlgorithm > d_fill_bdry_compute_stable_dt;
  Array< boost::shared_ptr< RefineSchedule > >
    d_fill_bdry_sched_compute_stable_dt;

  // for filling bdry data before doing time advance
  Array< boost::shared_ptr< RefineAlgorithm > > d_fill_bdry_time_advance;
  Array< Array< boost::shared_ptr< RefineSchedule > > >
    d_fill_bdry_sched_time_advance;

//Variable Database
 // hier::VariableDatabase* var_db;
  hier::PatchDataRestartManager* pdrm;
private:


  /*
   * Private assignment operator to prevent use.
   *
   * Arguments:
   *  - rhs (in):    LevelSetFunctionIntegrator to copy
   *
   * Return value:   *this
   *
   */
  const LevelSetFunctionIntegrator& operator=(
    const LevelSetFunctionIntegrator& rhs){
      return *this;
  }

};

} // end LSMLIB namespace


#endif
