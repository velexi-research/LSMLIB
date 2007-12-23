/*
 * File:        ReinitializationAlgorithm.h
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.17 $
 * Modified:    $Date: 2006/10/05 15:03:46 $
 * Description: Header file for level set method reinitialization algorithm 
 */
 
#ifndef included_ReinitializationAlgorithm_h
#define included_ReinitializationAlgorithm_h

/*! \class LSMLIB::ReinitializationAlgorithm
 *
 * \brief 
 * ReinitializationAlgorithm is a utility class that provides support for 
 * reinitializing the a level set function to be a distance function.  
 *
 * It reinitializes the level set functions to be distance functions using 
 * the reinitialization equation:
 *
 * \f[
 *
 *   \phi_t + sgn(\phi) ( |\nabla \phi| - 1 ) = 0
 *
 * \f] 
 *
 * This Hamilton-Jacobi equation is advanced in time towards steady-state
 * using a user-specified TVD Runge-Kutta method.  The number of steps 
 * taken is a function of the user-specified input parameters: 
 * stop_distance, max_iterations, and iteration_stop_tolerance.
 * The numerical discretization used to compute grad(phi) is also
 * user-specified.  A Godunov scheme is used to select the appropriate 
 * spatial derivative approximation for each component of grad(phi).
 *
 *
 * <h3> USAGE: </h3>
 *
 *  -# Set up level set functions and initialize the level set calculation
 *     (e.g. using the LevelSetMethodAlgorithm class).
 *  -# Create a ReinitializationAlgorithm object.
 *  -# Invoke @ref resetHierarchyConfiguration() if PatchHierarchy 
 *     configuration of has been changed since last reinitialization 
 *     procedure. 
 *  -# Reinitialize level set functions using the 
 *     @ref reinitializeLevelSetFunctions() or the
 *     @ref reinitializeLevelSetFunctionForSingleComponent() methods.
 *
 * 
 * <h3> User-specified parameters </h3>
 * 
 * - spatial_derivative_type    = type of spatial derivative calculation
 * - spatial_derivative_order   = order of spatial derivative
 * - tvd_runge_kutta_order      = order of TVD Runge-Kutta integration
 * - cfl_number                 = CFL number (default = 0.5)
 * - stop_distance              = approximate stopping criterion for  
 *                                evolution of reinitialization equation.
 *                                Reinitialization terminates after the
 *                                information from the zero level set has
 *                                propagated by approximately the specified
 *                                distance.  (default = 0.0)
 * - max_iterations             = maximum number of time steps to take during
 *                                the reinitialization process (default = 0)
 * - iteration_stop_tolerance   = stopping criterion for termination of 
 *                                evolution of reinitialization equation.
 *                                Reinitialization stops when the max norm
 *                                of the change in the level set function
 *                                drops below the specified tolerance.
 *                                (default = 0.0)
 * - verbose_mode               = flag to activate/deactivate verbose-mode 
 *                                (default = false)
 *
 *
 * <h3> NOTES: </h3>
 *
 * - If no stopping criteria are specified, the reinitialization
 *   calculation is terminated using the stop_distance criterion
 *   with the stop_distance set to the length of the largest 
 *   dimension of the computational domain.  This choice ensures
 *   that the level set function is reinitialized throughout the 
 *   entire computational domain.
 * 
 * - Unless it is important that the reinitialized level set function
 *   is a distance function (i.e. |grad(phi)|) to very high accuracy, 
 *   it is recommended that reinitialzation is done using a low-order
 *   spatial discretization and time integration scheme (e.g. ENO1 and 
 *   TVD Runge-Kutta 1).  This will significantly improve preformance.
 *
 * - A reinitialization calculation that extends throughout the
 *   entire domain can take a rather large amount of computational
 *   time.  The calculation may be sped up by specifying a smaller
 *   stop_distance (which effectively limits the distance off of the 
 *   zero level set where the level set function is reinitialized to 
 *   be a distance function) or using a lower order spatial-and 
 *   time-discretization.
 * 
 */


#include <vector>
#include "SAMRAI_config.h"
#include "CartesianGridGeometry.h"
#include "ComponentSelector.h"
#include "PatchHierarchy.h"
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "LSMLIB_config.h"
#include "BoundaryConditionModule.h"
#include "LevelSetMethodToolbox.h"

// SAMRAI namespaces 
using namespace SAMRAI;
using namespace hier;
using namespace tbox;
using namespace xfer;


/******************************************************************
 *
 * ReinitializationAlgorithm Class Definition
 *
 ******************************************************************/

namespace LSMLIB {

template<int DIM> class ReinitializationAlgorithm
{

public:

  //! @{ 
  /*!
   ****************************************************************
   *
   * @name Constructors and destructor
   *
   ****************************************************************/

  /*!
   * This constructor sets up the required variables, PatchData, etc. for 
   * reinitializing the level set function to be a distance function via
   * the reinitialization equation using parameters provided through 
   * the specified input database.
   *
   * Arguments:      
   *  - input_db (in):               input database containing user-defined
   *                                 parameters
   *  - hierarchy (in):              Pointer to PatchHierarchy 
   *                                 containing data
   *  - phi_handle (in):             PatchData handle for the level set
   *                                 function (phi)
   *  - control_volume_handle (in):  PatchData handle for control volume
   *                                 data
   *  - object_name (in):            string name for object (default = 
   *                                 "ReinitializationAlgorithm")
   *
   * Return value:                   none
   *
   * NOTES:
   *  - Only one ReinitializationAlgorithm object may be 
   *    created for each level set function.
   *  - Linear interpolation is used to refine data from coarser to 
   *    finer levels.
   *
   */
  ReinitializationAlgorithm(
    Pointer<Database> input_db,
    Pointer< PatchHierarchy<DIM> > hierarchy,
    const int phi_handle,
    const int control_volume_handle,
    const string& object_name = "ReinitializationAlgorithm");

  /*!
   * This constructor sets up the required variables, PatchData, etc. for 
   * reinitializing the level set function to be a distance function via
   * the reinitialization equation using parameters provided as arguments 
   * to the constructor.
   *
   * Arguments:      
   *  - hierarchy (in):                 Pointer to PatchHierarchy 
   *                                    containing data
   *  - phi_handle (in):                PatchData handle for the level set
   *                                    function (phi)
   *  - control_volume_handle (in):     PatchData handle for control volume
   *                                    data
   *  - spatial_derivative_type (in):   type of spatial derivative calculation
   *  - spatial_derivative_order (in):  order of spatial derivative
   *  - tvd_runge_kutta_order (in):     order of TVD Runge-Kutta integration
   *  - cfl_number (in):                CFL number used to compute dt
   *  - stop_distance (in):             approximate stopping criterion for 
   *                                    evolution of reinitialization equation.
   *                                    Reinitialization stops after the 
   *                                    information from the zero level set 
   *                                    has propagated by approximately the 
   *                                    specified distance. (default = 0.0) 
   *  - max_iterations (in):            maximum number of time steps to take
   *                                    during the reinitialization process
   *                                    (default = 0)
   *  - iteration_stop_tolerance (in):  stopping criterion for termination of 
   *                                    evolution of reinitialization equation.
   *                                    Reinitialization stops when the max 
   *                                    norm of the change in level set function
   *                                    drops below the specified tolerance.
   *                                    (default = 0.0) 
   *  - verbose_mode (in):              flag to activate/deactivate 
   *                                    verbose-mode (default = false)
   *  - object_name (in):               string name for object (default = 
   *                                    "ReinitializationAlgorithm")
   *
   * Return value:                      none
   *
   * NOTES:
   *  - Only one ReinitializationAlgorithm object may be 
   *    created for each level set function.
   *  - Linear interpolation is used to refine data from coarser to 
   *    finer levels.
   *
   */
  ReinitializationAlgorithm(
    Pointer< PatchHierarchy<DIM> > hierarchy,
    const int phi_handle,
    const int control_volume_handle,
    const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
    const int spatial_derivative_order,
    const int tvd_runge_kutta_order,
    const LSMLIB_REAL cfl_number,
    const LSMLIB_REAL stop_distance = 0.0,
    const int max_iterations = 0,
    const LSMLIB_REAL iteration_stop_tolerance = 0.0,
    const bool verbose_mode = false,
    const string& object_name = "ReinitializationAlgorithm");

  /*!
   * The destructor does nothing.
   *
   * Arguments:  none
   *
   */
  inline virtual ~ReinitializationAlgorithm(){}

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for managing reinitialization calculation
   *
   ****************************************************************/

  /*!
   * reinitializeLevelSetFunctions() reinitializes all of the components 
   * of the level set functions corresponding to the PatchData handle
   * passed into the constructor (i.e. phi_handle).
   *
   * Arguments:     
   *  - max_iterations (in):  maximum number of iterations to use
   *                          for reinitialization.  Set max_iterations
   *                          to -1 to use the value specified in the 
   *                          input file.
   *                          (default = -1)
   *  - lower_bc (in):        vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the lower face of the computational
   *                          domain in .  The i-th entry should
   *                          contain the type of boundary condition
   *                          impose at the lower boundary in the
   *                          i-th coordinate direction.
   *                          For information about the boundary 
   *                          condition types, see documentation of 
   *                          BoundaryConditionModule class.
   *                          (default = vector of -1's)
   *  - upper_bc (in):        vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the lower face of the computational
   *                          domain in .  The i-th entry should
   *                          contain the type of boundary condition
   *                          impose at the upper boundary in the
   *                          i-th coordinate direction.
   *                          For information about the boundary 
   *                          condition types, see documentation of 
   *                          BoundaryConditionModule class.
   *                          (default = vector of -1's)
   *
   * Return value:            none
   *
   * NOTES:
   *  - If max_iterations is set to a non-negative value, it overrides
   *    ALL of the stopping criteria specified in the input file.
   * 
   *  - The default behavior when no boundary conditions are supplied
   *    (i.e. either lower_bc or upper_bc is a vector of -1's) is that
   *    homogeneous Neumann boundary conditions will be imposed on
   *    non-periodic boundaries.
   *
   *  - The boundary conditions specified by lower_bc and upper_bc
   *    are used for ALL of the components of the level set function
   *    that are being reinitialized.
   *
   */
  virtual void reinitializeLevelSetFunctions(
    const int max_iterations = -1,
    const IntVector<DIM>& lower_bc = IntVector<DIM>(-1),
    const IntVector<DIM>& upper_bc = IntVector<DIM>(-1));

  /*!
   * reinitializeLevelSetFunctionForSingleComponent() reinitializes the
   * specified component of the level set function corresponding to the 
   * PatchData handle, phi_handle, passed into the constructor.
   *
   * Arguments:      
   *  - component (in):       component of level set function to 
   *                          reinitialize
   *                          (default = 0)
   *  - max_iterations (in):  maximum number of iterations to use
   *                          for reinitialization.  Set max_iterations
   *                          to -1 to use the value specified in the 
   *                          input file.
   *                          (default = -1)
   *  - lower_bc (in):        vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the lower face of the computational
   *                          domain in each coordinate direction.  
   *                          The i-th entry should contain the type 
   *                          of boundary condition to impose at the lower 
   *                          boundary in the i-th coordinate direction.
   *                          For information about the boundary 
   *                          condition types, see documentation of 
   *                          BoundaryConditionModule class.
   *                          (default = vector of -1's)
   *  - upper_bc (in):        vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the upper face of the computational
   *                          domain in each coordinate direction.  
   *                          The i-th entry should contain the type 
   *                          of boundary condition to impose at the upper
   *                          boundary in the i-th coordinate direction.
   *                          For information about the boundary 
   *                          condition types, see documentation of 
   *                          BoundaryConditionModule class.
   *                          (default = vector of -1's)
   *
   * Return value:            none
   *
   * NOTES:
   *  - If max_iterations is set to a non-negative value, it overrides
   *    ALL of the stopping criteria specified in the input file.
   *
   *  - The default behavior when no boundary conditions are supplied
   *    (i.e. either lower_bc or upper_bc is a vector of -1's) is that
   *    homogeneous Neumann boundary conditions will be imposed on
   *    non-periodic boundaries.
   *
   */
  virtual void reinitializeLevelSetFunctionForSingleComponent(
    const int component = 0,
    const int max_iterations = -1,
    const IntVector<DIM>& lower_bc = IntVector<DIM>(-1),
    const IntVector<DIM>& upper_bc = IntVector<DIM>(-1));

  //! @}


  //! @{
  /*!
   ****************************************************************
   * 
   * @name Methods for managing grid configuration
   *
   ****************************************************************/

  /*!
   * resetHierarchyConfiguration() updates resets the communication 
   * schedules to be consistent with the new configuration of the 
   * PatchHierarchy.  This method should be invoked any time the 
   * PatchHierarchy has been changed.
   *
   * Arguments:      
   *  - hierarchy (in):       Pointer to new PatchHierarchy 
   *  - coarsest_level (in):  coarsest level in the hierarchy to be updated
   *  - finest_level (in):    finest level in the hierarchy to be updated
   *
   * Return value:            none
   *
   */
  virtual void resetHierarchyConfiguration(
    Pointer< PatchHierarchy<DIM> > hierarchy,
    const int coarsest_level,
    const int finest_level);

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
   * advanceReinitializationEqnUsingTVDRK*() advances the reinitialization
   * equation a first-, second-, or third-order TVD Runge-Kutta step.
   *
   * Arguments:
   *  - dt (in):                  time increment to advance the level set 
   *                              functions
   *  - component (in):           component of level set function to advance 
   *                              in time
   *  - lower_bc (in):            vector of integers specifying the
   *                              type of boundary conditions to impose
   *                              on the lower face of the computational
   *                              domain in each coordinate direction for
   *                              the level set function.  The i-th entry 
   *                              should contain the type of boundary 
   *                              condition to impose at the lower boundary 
   *                              in the i-th coordinate direction.
   *                              For information about the boundary
   *                              condition types, see documentation of
   *                              BoundaryConditionModule class.
   *  - upper_bc (in):            vector of integers specifying the
   *                              type of boundary conditions to impose
   *                              on the upper face of the computational
   *                              domain in each coordinate direction for
   *                              the level set function.  The i-th entry 
   *                              should contain the type of boundary 
   *                              condition to impose at the upper boundary 
   *                              in the i-th coordinate direction.
   *                              For information about the boundary
   *                              condition types, see documentation of
   *                              BoundaryConditionModule class.
   *
   * Return value:                none
   *
   * NOTES:
   *  - These methods are NOT intended to be used directly by the 
   *    user.  They are only helper methods for the 
   *    reinitializeLevelSetFunctionForSingleComponent() and
   *    reinitializeLevelSetFunctions() methods
   *
   */
  virtual void advanceReinitializationEqnUsingTVDRK1(
    const LSMLIB_REAL dt,
    const int component,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc);
  virtual void advanceReinitializationEqnUsingTVDRK2(
    const LSMLIB_REAL dt,
    const int component,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc);
  virtual void advanceReinitializationEqnUsingTVDRK3(
    const LSMLIB_REAL dt,
    const int component,
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc);

  /*!
   * computeReinitializationEqnRHS() computes the right-hand side of
   * the reinitialization equation when it is written in the form:
   *
   *   phi_t = ...
   *
   * Arguments:
   *  - phi_handle (in):  PatchData handle to use in computing RHS of 
   *                      reinitialization equation
   *
   * Return value:        none
   *
   */
  virtual void computeReinitializationEqnRHS(
    const int phi_handle);

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
   * for the reinitialization calculation.
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
   * getFromInput() configures the ReinitializationAlgorithm 
   * object from the values in the specified input database.
   *
   * Arguments:
   *  - db (in):               input database from which to take parameter
   *                           values
   *
   * Return value:             none
   *
   * NOTES:
   *  - An assertion results if the database pointer is null.
   *
   */
  virtual void getFromInput(Pointer<Database> db);

  /*!
   * checkParameters() checks to make sure that user-specified parameters
   * are acceptable.  This method throws unrecoverable errors if any of the
   * parameters are unacceptable.
   *
   * Arguments:     none
   *
   * Return value:  none
   *
   */
  virtual void checkParameters();

  //! @}


  /****************************************************************
   *
   * Data Members
   *
   ****************************************************************/

  /*
   * The object name is used for error/warning reporting. 
   */
  string d_object_name;

  /*
   * User-defined parameters
   */

  // numerical parameters
  SPATIAL_DERIVATIVE_TYPE d_spatial_derivative_type;
  int d_spatial_derivative_order;
  int d_tvd_runge_kutta_order;
  LSMLIB_REAL d_cfl_number;
  LSMLIB_REAL d_stop_distance;
  int d_max_iterations;
  LSMLIB_REAL d_iteration_stop_tol;

  // verbose mode
  bool d_verbose_mode;

  /*
   * Grid management objects 
   */

  // Pointer to PatchHierarchy object
  Pointer< PatchHierarchy<DIM> > d_patch_hierarchy;

  // Pointer to GridGeometry
  Pointer< CartesianGridGeometry<DIM> > d_grid_geometry;

  /*
   * PatchData handles for data required to solve reinitialization equation
   */

  // data field handles
  int d_phi_handle;
  int d_control_volume_handle;

  // scratch data 
  vector<int> d_phi_scr_handles;
  IntVector<DIM> d_phi_scratch_ghostcell_width;
  int d_rhs_handle;
  int d_grad_phi_plus_handle;
  int d_grad_phi_minus_handle;

  /* 
   * internal state  variables 
   */

  // iteration termination flags
  bool d_use_stop_distance;
  bool d_use_max_iterations;
  bool d_use_iteration_stop_tol;

  // flag indicating that communication schedules need to be recomputed
  bool d_hierarchy_configuration_needs_reset; 

  // level set data parameters
  int d_num_phi_components;

  // ComponentSelector to organize variables
  ComponentSelector d_scratch_data;

  /*
   * Boundary condition objects
   */
  Pointer< BoundaryConditionModule<DIM> > d_bc_module;

  /*
   * Communication objects.
   */

  // data communication parameters
  Array< Pointer< RefineAlgorithm<DIM> > > d_phi_fill_bdry_alg;
  Array< Array< Pointer< RefineSchedule<DIM> > > > d_phi_fill_bdry_sched;

private: 

  /*
   * Private copy constructor to prevent use.
   * 
   * Arguments:
   *  - rhs (in):  object to copy
   *
   */
  ReinitializationAlgorithm(
    const ReinitializationAlgorithm& rhs){}
   
  /*
   * Private assignment operator to prevent use.
   *
   * Arguments:
   *  - rhs (in):    object to copy
   *
   * Return value:   *this
   *
   */
  const ReinitializationAlgorithm& operator=(
    const ReinitializationAlgorithm& rhs){
      return *this;
  }

};

} // end LSMLIB namespace 

#endif
