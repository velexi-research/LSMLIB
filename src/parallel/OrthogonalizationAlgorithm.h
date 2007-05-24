/*
 * File:        OrthogonalizationAlgorithm.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.21 $
 * Modified:    $Date: 2006/10/05 15:03:45 $
 * Description: Header file for level set method orthogonalization algorithm 
 */
 
#ifndef included_OrthogonalizationAlgorithm_h
#define included_OrthogonalizationAlgorithm_h

/*! \class LSMLIB::OrthogonalizationAlgorithm
 *
 * \brief
 * OrthogonalizationAlgorithm is a utility class that provides support for 
 * orthogonalizing two level set functions, phi and psi. 
 *
 * This goal is achieved by solving either of the two following 
 * orthogonalization equations:
 *
 * \f[
 *
 *   \phi_t + sgn(\psi) * ( \nabla \psi / |\nabla \psi| ) \cdot \nabla \phi = 0
 *
 * \f]
 * \f[
 *
 *   \psi_t + sgn(\phi) * ( \nabla \phi / |\nabla \phi| ) \cdot \nabla \psi = 0
 *
 * \f]
 *
 * which evolve either \f$ \phi \f$ or \f$ \psi \f$ until \f$ \nabla \phi \f$ 
 * and \f$ \nabla \psi \f$ are.  
 * orthogonal.  These Hamilton-Jacobi equations are advanced in time towards 
 * steady-state using a user-specified TVD Runge-Kutta method.  The number 
 * of steps taken is a function of the user-specified input parameters: 
 * stop_distance, max_iterations, and iteration_stop_tolerance.  The 
 * numerical discretization used to compute \f$ \nabla \phi \f$ and 
 * \f$ \nabla \psi \f$ are also user-specified.  In the first of the above 
 * equation, \f$ \nabla \psi \f$ is computed using by taking the average 
 * of the forward and backward spatial derivatives, and \f$ \nabla \phi \f$ 
 * is computed via a simple upwinding scheme that treats \f$ \nabla \psi \f$ 
 * as the velocity.  In the second of the above equations, the reverse 
 * procedure is used to compute \f$ \nabla \phi \f$ and \f$ \nabla \psi \f$.
 *
 *
 * <h3> USAGE: </h3>
 *
 *  -# Set up level set functions and initialize the level set calculation
 *     (e.g. using the LevelSetMethodAlgorithm class).
 *  -# Create a OrthogonalizationAlgorithm object.
 *  -# Invoke @ref resetHierarchyConfiguration() if PatchHierarchy 
 *     configuration of has been changed since last orthogonalization 
 *     procedure. 
 *  -# Orthogonalize level set functions using the 
 *     @ref orthogonalizeLevelSetFunctions() or the
 *     @ref orthogonalizeLevelSetFunctionForSingleComponent() methods.
 *
 * 
 * <h3> User-specified parameters </h3>
 * 
 * - spatial_derivative_type    = type of spatial derivative calculation
 * - spatial_derivative_order   = order of spatial derivative
 * - tvd_runge_kutta_order      = order of TVD Runge-Kutta integration
 * - cfl_number                 = CFL number (default = 0.5)
 * - stop_distance              = approximate stopping criterion for  
 *                                evolution of orthogonalization equation.
 *                                Orthogonalization terminates after the
 *                                information from the zero level set has
 *                                propagated by approximately the specified
 *                                distance.  (default = 0.0)
 * - max_iterations             = maximum number of time steps to take during
 *                                the orthogonalization process (default = 0)
 * - iteration_stop_tolerance   = stopping criterion for termination of 
 *                                evolution of orthogonalization equation.
 *                                Orthogonalization stops when the max norm
 *                                of the change in the level set function
 *                                drops below the specified tolerance.
 *                                (default = 0.0)
 * - verbose_mode               = flag to activate/deactivate verbose-mode 
 *                                (default = false)
 *
 *
 * <h3> NOTES: </h3>
 *
 * - If no stopping criteria are specified, the orthogonalization
 *   calculation is terminated using the stop_distance criterion
 *   with the stop_distance set to the length of the largest 
 *   dimension of the computational domain.  This choice ensures
 *   that the level set functions are orthogonalized throughout the 
 *   entire computational domain.
 * 
 * - Unless it is important that the gradients of the level set 
 *   functions are orthogonal to very high accuracy, it is 
 *   recommended that orthogonalization is done using a low-order
 *   spatial discretization and time integration scheme (e.g. ENO1 and 
 *   TVD Runge-Kutta 1).  This will significantly improve preformance.
 *
 * - An orthogonalization calculation that extends throughout the
 *   entire domain can take a rather large amount of computational
 *   time.  The calculation may be sped up by specifying a smaller
 *   stop_distance (which effectively limits the distance off of
 *   the zero level sets that the level set functions are 
 *   orthogonalized) or using a lower order spatial-and time-discretization.
 * 
 */


#include "SAMRAI_config.h"
#include "CartesianGridGeometry.h"
#include "PatchHierarchy.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "LSMLIB_config.h"
#include "LevelSetMethodToolbox.h"
#include "FieldExtensionAlgorithm.h"
#include "BoundaryConditionModule.h"

// SAMRAI namespaces 
using namespace SAMRAI;
using namespace hier;
using namespace tbox;
using namespace xfer;


/******************************************************************
 *
 * OrthogonalizationAlgorithm Class Definition
 *
 ******************************************************************/

namespace LSMLIB {

template<int DIM> class OrthogonalizationAlgorithm
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
   * This constructor sets up the required variables, PatchData, etc. 
   * for orthogonalizing two level set functions when using the level 
   * set method for codimension-two problems.  Parameters for the 
   * orthogonalization are provided through the specified input database.
   *
   * Arguments:      
   *  - input_db (in):               input database containing user-defined
   *                                 parameters
   *  - hierarchy (in):              Pointer to PatchHierarchy 
   *                                 containing data
   *  - phi_handle (in):             PatchData handle for phi
   *  - psi_handle (in):             PatchData handle for psi
   *  - control_volume_handle (in):  PatchData handle for control volume
   *                                 data
   *  - object_name (in):            string name for object (default = 
   *                                 "OrthogonalizationAlgorithm")
   *  - phi_ghostcell_width (in):    ghostcell width for the data associated
   *                                 with phi_handle.  It is used to 
   *                                 determine if it is necessary to
   *                                 allocate scratch space when computing
   *                                 grad(phi) during orthogonalization
   *                                 procedures where phi is held fixed.  
   *                                 If a phi_handle that possesses a 
   *                                 sufficient number of ghostcells for 
   *                                 the grad(phi) calculation AND the 
   *                                 phi_ghostcell_width argument is used, 
   *                                 the memory required to orthogonalize
   *                                 the level set functions is reduced.
   *                                 (default = [0,0,0])
   *  - psi_ghostcell_width (in):    ghostcell width for the data associated
   *                                 with psi_handle.  It is used to 
   *                                 determine if it is necessary to
   *                                 allocate scratch space when computing
   *                                 grad(psi) during orthogonalization
   *                                 procedures where psi is held fixed.  
   *                                 If a psi_handle that possesses a 
   *                                 sufficient number of ghostcells for 
   *                                 the grad(psi) calculation AND the 
   *                                 psi_ghostcell_width argument is used, 
   *                                 the memory required to orthogonalize
   *                                 the level set functions is reduced.
   *                                 (default = [0,0,0])
   *
   * Return value:                   none
   *
   * NOTES:
   *  - Only one OrthogonalizationAlgorithm object may be created 
   *    for each level set function.
   *  - phi and psi MUST have the same number of components.
   *  - Linear interpolation is used to refine data from coarser to 
   *    finer levels.
   *
   */
  OrthogonalizationAlgorithm(
    Pointer<Database> input_db,
    Pointer< PatchHierarchy<DIM> > hierarchy,
    const int phi_handle,
    const int psi_handle,
    const int control_volume_handle,
    const string& object_name = "OrthogonalizationAlgorithm",
    const IntVector<DIM>& phi_ghostcell_width = 0,
    const IntVector<DIM>& psi_ghostcell_width = 0);

  /*!
   * This constructor sets up the required variables, PatchData, etc. 
   * for orthogonalizing two level set functions when using the level
   * set method for codimension-two problems.  Parameters for the 
   * orthogonalization procedure are provided through the specified 
   * input database.
   *
   * Arguments:      
   *  - hierarchy (in):                 Pointer to PatchHierarchy 
   *                                    containing data
   *  - phi_handle (in):                PatchData handle for phi
   *  - psi_handle (in):                PatchData handle for psi
   *  - control_volume_handle (in):     PatchData handle for control volume
   *                                    data
   *  - spatial_derivative_type (in):   type of spatial derivative calculation
   *  - spatial_derivative_order (in):  order of spatial derivative
   *  - tvd_runge_kutta_order (in):     order of TVD Runge-Kutta integration
   *  - cfl_number (in):                CFL number used to compute dt
   *  - stop_distance (in):             approximate stopping criterion for 
   *                                    evolution of orthogonalization equation.
   *                                    Orthogonalization stops after the 
   *                                    information from the zero level set 
   *                                    has propagated by approximately the 
   *                                    specified distance. (default = 0.0) 
   *  - max_iterations (in):            maximum number of time steps to take
   *                                    during the orthogonalization process
   *                                    (default = 0)
   *  - iteration_stop_tolerance (in):  stopping criterion for termination of 
   *                                    evolution of orthogonalization equation.
   *                                    Orthogonalization stops when the max 
   *                                    norm of the change in level set function
   *                                    drops below the specified tolerance.
   *                                    (default = 0.0) 
   *  - verbose_mode (in):              flag to activate/deactivate 
   *                                    verbose-mode (default = false)
   *  - object_name (in):               string name for object (default = 
   *                                    "OrthogonalizationAlgorithm")
   *  - phi_ghostcell_width (in):       ghostcell width for the data associated
   *                                    with phi_handle.  It is used to 
   *                                    determine if it is necessary to
   *                                    allocate scratch space when computing
   *                                    grad(phi) during orthogonalization
   *                                    procedures where phi is held fixed.  
   *                                    If a phi_handle that possesses a 
   *                                    sufficient number of ghostcells for 
   *                                    the grad(phi) calculation AND the 
   *                                    phi_ghostcell_width argument is used, 
   *                                    the memory required to orthogonalize
   *                                    the level set functions is reduced.
   *                                    (default = [0,0,0])
   *  - psi_ghostcell_width (in):       ghostcell width for the data associated
   *                                    with psi_handle.  It is used to 
   *                                    determine if it is necessary to
   *                                    allocate scratch space when computing
   *                                    grad(psi) during orthogonalization
   *                                    procedures where psi is held fixed.  
   *                                    If a psi_handle that possesses a 
   *                                    sufficient number of ghostcells for 
   *                                    the grad(psi) calculation AND the 
   *                                    psi_ghostcell_width argument is used, 
   *                                    the memory required to orthogonalize
   *                                    the level set functions is reduced.
   *                                    (default = [0,0,0])
   *   
   * Return value:                      none
   *
   * NOTES:
   *  - Only one OrthogonalizationAlgorithm object may be 
   *    created for each level set function.
   *  - phi and psi MUST have the same number of components.
   *  - Linear interpolation is used to refine data from coarser to 
   *    finer levels.
   *
   */
  OrthogonalizationAlgorithm(
    Pointer< PatchHierarchy<DIM> > hierarchy,
    const int phi_handle,
    const int psi_handle,
    const int control_volume_handle,
    const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
    const int spatial_derivative_order,
    const int tvd_runge_kutta_order,
    const double cfl_number,
    const double stop_distance = 0.0,
    const int max_iterations = 0,
    const double iteration_stop_tolerance = 0.0,
    const bool verbose_mode = false,
    const string& object_name = "OrthogonalizationAlgorithm",
    const IntVector<DIM>& phi_ghostcell_width = 0,
    const IntVector<DIM>& psi_ghostcell_width = 0);

  /*!
   * The destructor does nothing.
   *
   * Arguments:  none
   *
   */
  inline virtual ~OrthogonalizationAlgorithm(){}

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for managing orthogonalization calculation
   *
   ****************************************************************/

  /*!
   * orthogonalizeLevelSetFunctions() orthogonalizes grad(phi) and 
   * grad(psi) for all components of the level set functions corresponding 
   * to the PatchData handles passed into the constructor (i.e. phi_handle 
   * and psi_handle).  The level set function specified by the argument, 
   * level_set_fcn, is the one that is evolved to orthogonalize grad(phi) 
   * and grad(psi).
   * 
   * Arguments:     
   *  - level_set_fcn (in):     level set function which is evolved 
   *                            to produce orthogonal gradients
   *  - max_iterations (in):    maximum number of iterations to use
   *                            for orthogonalization.  Set 
   *                            max_iterations to -1 to use the 
   *                            value specified in the input file.
   *                            (default = -1)
   *  - lower_bc_fixed (in):    vector of integers specifying the
   *                            type of boundary conditions to impose
   *                            on the lower face of the computational
   *                            domain in each coordinate direction for 
   *                            the fixed level set function.
   *                            The i-th entry should contain the type 
   *                            of boundary condition to impose at the 
   *                            lower boundary in the i-th coordinate 
   *                            direction.
   *                            For information about the boundary
   *                            condition types, see documentation of
   *                            BoundaryConditionModule class.
   *                            (default = vector of -1's)
   *  - upper_bc_fixed (in):    vector of integers specifying the
   *                            type of boundary conditions to impose
   *                            on the upper face of the computational
   *                            domain in each coordinate direction for 
   *                            the fixed level set function.
   *                            The i-th entry should contain the type 
   *                            of boundary condition to impose at the 
   *                            upper boundary in the i-th coordinate 
   *                            direction.
   *                            For information about the boundary
   *                            condition types, see documentation of
   *                            BoundaryConditionModule class.
   *                            (default = vector of -1's)
   *  - lower_bc_evolved (in):  vector of integers specifying the
   *                            type of boundary conditions to impose
   *                            on the lower face of the computational
   *                            domain in each coordinate direction for 
   *                            the evolved level set function.
   *                            The i-th entry should contain the type 
   *                            of boundary condition to impose at the 
   *                            lower boundary in the i-th coordinate 
   *                            direction.
   *                            For information about the boundary
   *                            condition types, see documentation of
   *                            BoundaryConditionModule class.
   *                            (default = vector of -1's)
   *  - upper_bc_evolved (in):  vector of integers specifying the
   *                            type of boundary conditions to impose
   *                            on the upper face of the computational
   *                            domain in each coordinate direction for 
   *                            the evolved level set function.
   *                            The i-th entry should contain the type 
   *                            of boundary condition to impose at the 
   *                            upper boundary in the i-th coordinate 
   *                            direction.
   *                            For information about the boundary
   *                            condition types, see documentation of
   *                            BoundaryConditionModule class.
   *                            (default = vector of -1's)
   *
   * Return value:              none
   *
   * NOTES:
   *  - If max_iterations is set to a non-negative value, it overrides
   *    ALL of the stopping criteria specified in the input file.
   *
   *  - It is recommended that the level set function which is held
   *    fixed during the orthogonalization process is reinitialized 
   *    before invoking orthogonalizeLevelSetFunctions() and that
   *    the level set function which is evolved during the orthogonalization 
   *    process is reinitialized after invoking 
   *    orthogonalizeLevelSetFunctions(). 
   *
   *  - The default behavior when no boundary conditions are supplied
   *    (i.e. either lower_bc or upper_bc is a vector of -1's) is that
   *    homogeneous Neumann boundary conditions will be imposed on
   *    non-periodic boundaries.
   *
   */
  virtual void orthogonalizeLevelSetFunctions(
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int max_iterations = -1,
    const IntVector<DIM>& lower_bc_fixed = IntVector<DIM>(-1),
    const IntVector<DIM>& upper_bc_fixed = IntVector<DIM>(-1),
    const IntVector<DIM>& lower_bc_evolved = IntVector<DIM>(-1),
    const IntVector<DIM>& upper_bc_evolved = IntVector<DIM>(-1));

  /*!
   * orthogonalizeLevelSetFunctionForSingleComponent() orthogonalizes 
   * grad(phi) and grad(psi) for the specified component of the level 
   * set function corresponding to the the level set functions 
   * corresponding to the PatchData handles passed into the constructor 
   * (i.e. phi_handle and psi_handle).  The level set function specified 
   * by the argument, level_set_fcn, is the one that is evolved to 
   * orthogonalize grad(phi) and grad(psi).
   *
   * Arguments:      
   *  - level_set_fcn (in):     level set function which is evolved 
   *                            to produce orthogonal gradients
   *  - component (in):         component to level set functions to 
   *                            orthogonalize (default = 0)
   *  - max_iterations (in):    maximum number of iterations to use
   *                            for orthogonalization.  Set 
   *                            max_iterations to -1 to use the 
   *                            value specified in the input file.
   *                            (default = -1)
   *  - lower_bc_fixed (in):    vector of integers specifying the
   *                            type of boundary conditions to impose
   *                            on the lower face of the computational
   *                            domain for the fixed level set function.
   *                            The i-th entry should contain the type 
   *                            of boundary condition to impose at the 
   *                            lower boundary in the i-th coordinate 
   *                            direction.
   *                            For information about the boundary
   *                            condition types, see documentation of
   *                            BoundaryConditionModule class.
   *                            (default = vector of -1's)
   *  - upper_bc_fixed (in):    vector of integers specifying the
   *                            type of boundary conditions to impose
   *                            on the upper face of the computational
   *                            domain for the fixed level set function.
   *                            The i-th entry should contain the type 
   *                            of boundary condition to impose at the 
   *                            upper boundary in the i-th coordinate 
   *                            direction.
   *                            For information about the boundary
   *                            condition types, see documentation of
   *                            BoundaryConditionModule class.
   *                            (default = vector of -1's)
   *  - lower_bc_evolved (in):  vector of integers specifying the
   *                            type of boundary conditions to impose
   *                            on the lower face of the computational
   *                            domain for the evolved level set function.
   *                            The i-th entry should contain the type 
   *                            of boundary condition to impose at the 
   *                            lower boundary in the i-th coordinate 
   *                            direction.
   *                            For information about the boundary
   *                            condition types, see documentation of
   *                            BoundaryConditionModule class.
   *                            (default = vector of -1's)
   *  - upper_bc_evolved (in):  vector of integers specifying the
   *                            type of boundary conditions to impose
   *                            on the upper face of the computational
   *                            domain for the evolved level set function.
   *                            The i-th entry should contain the type 
   *                            of boundary condition to impose at the 
   *                            upper boundary in the i-th coordinate 
   *                            direction.
   *                            For information about the boundary
   *                            condition types, see documentation of
   *                            BoundaryConditionModule class.
   *                            (default = vector of -1's)
   *
   * Return value:              none
   *
   * NOTES:
   *  - If max_iterations is set to a non-negative value, it overrides
   *    ALL of the stopping criteria specified in the input file.
   *
   *  - It is recommended that the level set function which is held
   *    fixed during the orthogonalization process is reinitialized 
   *    before invoking orthogonalizeLevelSetFunctions() and that
   *    the level set function which is evolved during the orthogonalization 
   *    process is reinitialized after invoking 
   *    orthogonalizeLevelSetFunctions(). 
   *
   *  - The default behavior when no boundary conditions are supplied
   *    (i.e. either lower_bc or upper_bc is a vector of -1's) is that
   *    homogeneous Neumann boundary conditions will be imposed on
   *    non-periodic boundaries.
   *  
   */
  virtual void orthogonalizeLevelSetFunctionForSingleComponent(
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int component = 0,
    const int max_iterations = -1,
    const IntVector<DIM>& lower_bc_fixed = IntVector<DIM>(-1),
    const IntVector<DIM>& upper_bc_fixed = IntVector<DIM>(-1),
    const IntVector<DIM>& lower_bc_evolved = IntVector<DIM>(-1),
    const IntVector<DIM>& upper_bc_evolved = IntVector<DIM>(-1));

  // @}


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
   * @name Utility methods
   *
   ****************************************************************/

  /*!
   * getFromInput() configures the OrthogonalizationAlgorithm 
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
  double d_cfl_number;
  double d_stop_distance;
  int d_max_iterations;
  double d_iteration_stop_tol;

  // verbose mode
  bool d_verbose_mode;

  /*
   * Grid management objects 
   */

  // Pointer to PatchHierarchy and GridGeometry object
  Pointer< PatchHierarchy<DIM> > d_patch_hierarchy;
  Pointer< CartesianGridGeometry<DIM> > d_grid_geometry;

  /*
   * Pointers to FieldExtensionAlgorithm objects 
   */
  Pointer< FieldExtensionAlgorithm<DIM> > d_fixed_phi_field_ext_alg;
  Pointer< FieldExtensionAlgorithm<DIM> > d_fixed_psi_field_ext_alg;

  /*
   * PatchData handles for data required to solve orthogonalization equation
   */

  // data field handles
  int d_phi_handle;
  int d_psi_handle;
  int d_control_volume_handle;

  /* 
   * internal state  variables 
   */

  // iteration termination flags
  bool d_use_stop_distance;
  bool d_use_max_iterations;
  bool d_use_iteration_stop_tol;

  // level set data parameters
  int d_num_field_components;

private: 

  /*
   * Private copy constructor to prevent use.
   * 
   * Arguments:
   *  - rhs (in):  object to copy
   *
   */
  OrthogonalizationAlgorithm(
    const OrthogonalizationAlgorithm& rhs){}
   
  /*
   * Private assignment operator to prevent use.
   *
   * Arguments:
   *  - rhs (in):    object to copy
   *
   * Return value:   *this
   *
   */
  const OrthogonalizationAlgorithm& operator=(
    const OrthogonalizationAlgorithm& rhs){
      return *this;
  }

};

} // end LSMLIB namespace

#endif
