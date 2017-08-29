/*
 * File:        FieldExtensionAlgorithm.h
 * Description: Header file for level set method field extension algorithm
 */

#ifndef included_FieldExtensionAlgorithm_h
#define included_FieldExtensionAlgorithm_h

/*! \class LSMLIB::FieldExtensionAlgorithm
 *
 * \brief
 * The FieldExtensionAlgorithm class is a utility class that provides
 * support for extending field values off of the interface defined by the
 * zero level set.
 *
 *
 * <h3> USAGE: </h3>
 *
 *  -# Set up level set functions and initialize the level set calculation
 *     (e.g. using the LevelSetMethodAlgorithm class).
 *  -# Set up and initialize the extension field variables.
 *  -# Create a FieldExtensionAlgorithm object or request one by using
 *     LevelSetMethodAlgorithm::getFieldExtensionAlgorithm().
 *  -# Invoke @ref resetHierarchyConfiguration() if PatchHierarchy
 *     configuration of has been changed since last field extension
 *     computation.
 *  -# Set values of the field variable in grid cells adjacent to
 *     the zero level set.
 *  -# Extend fields off of the zero level set using the
 *     @ref computeExtensionField() or
 *     @ref computeExtensionFieldForSingleComponent() methods.
 *
 * NOTE: steps (1) and (2) may be reversed.
 *
 *
 * <h3> User-specified parameters </h3>
 *
 * - spatial_derivative_type    = type of spatial derivative calculation
 * - spatial_derivative_order   = order of spatial derivative
 * - tvd_runge_kutta_order      = order of TVD Runge-Kutta integration
 * - cfl_number                 = CFL number (default = 0.5)
 * - stop_distance              = approximate stopping criterion for
 *                                evolution of field extension equation.
 *                                Field extension terminates after the
 *                                information from the zero level set has
 *                                propagated by approximately the specified
 *                                distance.  (default = 0.0)
 * - max_iterations             = maximum number of time steps to take during
 *                                the field extension process (default = 0)
 * - iteration_stop_tolerance   = stopping criterion for termination of
 *                                evolution of field extension equation.
 *                                Field extension stops when the max norm
 *                                of the change in the extension field
 *                                drops below the specified tolerance.
 *                                (default = 0.0)
 * - verbose_mode               = flag to activate/deactivate verbose-mode
 *                                (default = false)
 *
 * <h3> NOTES: </h3>
 *
 * - The values of the extension field in grid cells adjacent to
 *   the zero level set in the directions of the coordinate axes
 *   are guaranteed to remain unchanged by the field extension
 *   procedure.
 *
 * - If no stopping criteria are specified, the field extension
 *   calculation is terminated using the stop_distance criterion
 *   with the stop_distance set to the length of the largest
 *   dimension of the computational domain.  This choice ensures
 *   that the extension field is extended throughout the entire
 *   computational domain.
 *
 * - Unless it is important that the extension fields be computed
 *   to very high accuracy, it is recommended that field extension
 *   is done using a low-order spatial discretization and time
 *   integration scheme (e.g. ENO1 and TVD Runge-Kutta 1).  This will
 *   significantly improve preformance.
 *
 * - A field extension calculation that extends throughout the
 *   entire domain can take a rather large amount of computational
 *   time.  The calculation may be sped up by specifying a smaller
 *   stop_distance (which effectively limits the distance that
 *   the fields are extended off of the zero level set) or using
 *   a lower order spatial-and time-discretization.
 *
 */

// Standard library headers
#include <iosfwd>
#include <vector>

// SAMRAI headers
#include "boost/smart_ptr/shared_ptr.hpp"

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Array.h"

// LSMLIB headers
#include "LSMLIB_config.h"
#include "LevelSetMethodToolbox.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations
namespace LSMLIB { class BoundaryConditionModule; }
namespace SAMRAI { namespace geom { class CartesianGridGeometry; } }
namespace SAMRAI { namespace hier { class PatchHierarchy; } }
namespace SAMRAI { namespace tbox { class Database; } }
namespace SAMRAI { namespace xfer { class RefineAlgorithm; } }
namespace SAMRAI { namespace xfer { class RefineSchedule; } }

/******************************************************************
 *
 * FieldExtensionAlgorithm Class Definition
 *
 ******************************************************************/

namespace LSMLIB {

class FieldExtensionAlgorithm
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
   * computing extension fields via the advection equation using parameters
   * provided through the specified input database.
   *
   * Arguments:
   *  - input_db (in):                input database containing user-defined
   *                                  parameters
   *  - hierarchy (in):               Boost pointer to PatchHierarchy
   *                                  containing data
   *  - field_handle (in):            PatchData handle for field to
   *                                  extend off of the interface
   *  - phi_handle (in):              PatchData handle for the level set
   *                                  function (phi)
   *  - control_volume_handle (in):   PatchData handle for control volume
   *                                  data
   *  - object_name (in):             string name for object (default =
   *                                  "FieldExtensionAlgorithm")
   *  - phi_ghostcell_width (in):     ghostcell width for the data associated
   *                                  with phi_handle.  It is used to
   *                                  determine if it is necessary to
   *                                  allocate scratch space when computing
   *                                  grad(phi).  If a phi_handle that
   *                                  possesses a sufficient number of
   *                                  ghostcells for the grad(phi) calculation
   *                                  AND the phi_ghostcell_width argument
   *                                  is used, the memory required to compute
   *                                  extension fields is reduced.
   *                                  (default = [0,0,0])
   *
   * NOTES:
   *  - Only one FieldExtensionAlgorithm object may be
   *    created for each field handle.
   *  - Linear interpolation is used to refine data from coarser to
   *    finer levels.
   *
   */

 FieldExtensionAlgorithm(
    boost::shared_ptr<tbox::Database> input_db,
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int field_handle,
    const int phi_handle,
    const int control_volume_handle,
    const hier::IntVector& phi_ghostcell_width,
    const string& object_name = "FieldExtensionAlgorithm");

   /*!
   * This constructor sets up the required variables, PatchData, etc. for
   * computing extension fields via the advection equation using parameters
   * provided as arguments to the constructor.
   *
   * Arguments:
   *  - hierarchy (in):                 Boost pointer to PatchHierarchy
   *                                    containing data
   *  - field_handle (in):              PatchData handle for field to
   *                                    extend off of the interface
   *  - phi_handle (in):                PatchData handle for the level set
   *                                    function (phi)
   *  - control_volume_handle (in):     PatchData handle for control volume
   *                                    data
   *  - spatial_derivative_type (in):   type of spatial derivative calculation
   *                                    (default = ENO)
   *  - spatial_derivative_order (in):  order of spatial derivative
   *                                    (default = 1)
   *  - tvd_runge_kutta_order (in):     order of TVD Runge-Kutta integration
   *                                    (default = 1)
   *  - cfl_number (in):                CFL number used to compute dt
   *                                    (default = 0.5)
   *  - stop_distance (in):             approximate stopping criterion for
   *                                    evolution of field extension equation.
   *                                    Field extension stops after the
   *                                    information from the zero level set
   *                                    has propagated by approximately the
   *                                    specified distance. (default = 0.0)
   *  - max_iterations (in):            maximum number of time steps to take
   *                                    during the field extension process
   *                                    (default = 0)
   *  - iteration_stop_tolerance (in):  stopping criterion for termination of
   *                                    evolution of field extension equation.
   *                                    Field extension stops when the max norm
   *                                    of the change in the extension field
   *                                    drops below the specified tolerance.
   *                                    (default = 0.0)
   *  - verbose_mode (in):              flag to activate/deactivate
   *                                    verbose-mode (default = false)
   *  - object_name (in):               string name for object (default =
   *                                    "FieldExtensionAlgorithm")
   *  - phi_ghostcell_width (in):       ghostcell width for the data associated
   *                                    with phi_handle.  It is used to
   *                                    determine if it is necessary to
   *                                    allocate scratch space when computing
   *                                    grad(phi).  If a phi_handle that
   *                                    possesses a sufficient number of
   *                                    ghostcells for the grad(phi) calculation
   *                                    AND the phi_ghostcell_width argument
   *                                    is used, the memory required to compute
   *                                    extension fields is reduced.
   *                                    (default = [0,0,0])
   *
   * NOTES:
   *  - Only one FieldExtensionAlgorithm object may be
   *    created for each field handle.
   *  - Linear interpolation is used to refine data from coarser to
   *    finer levels.
   *
   */

  FieldExtensionAlgorithm(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int field_handle,
    const int phi_handle,
    const int control_volume_handle,
    const hier::IntVector& phi_ghostcell_width,
    const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type = ENO,
    const int spatial_derivative_order = 1,
    const int tvd_runge_kutta_order = 1,
    const LSMLIB_REAL cfl_number = 0.5,
    const LSMLIB_REAL stop_distance = 0.0,
    const int max_iterations = 0,
    const LSMLIB_REAL iteration_stop_tolerance = 0.0,
    const bool verbose_mode = false,
    const string& object_name = "FieldExtensionAlgorithm");

  /*!
   * The destructor does nothing.
   *
   * Arguments:  none
   *
   */
  inline virtual ~FieldExtensionAlgorithm(){}

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for managing extension field calculation
   *
   ****************************************************************/

  /*!
   * computeExtensionField() extends all of the components of the specified
   * field off of the interface defined by the zero level set of the
   * function phi by advecting it in the direction normal to the interface.
   *
   * Arguments:
   *  - phi_component (in):   component of level set function to
   *                          use in field extension calculation
   *  - lower_bc_phi (in):    vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the lower face of the computational
   *                          domain in each coordinate direction for
   *                          the level set function.
   *                          The i-th entry should contain the type of
   *                          boundary condition to impose at the lower
   *                          boundary in the i-th coordinate direction.
   *                          For information about the boundary
   *                          condition types, see documentation of
   *                          BoundaryConditionModule class.
   *  - upper_bc_phi (in):    vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the upper face of the computational
   *                          domain in each coordinate direction for
   *                          the level set function.
   *                          The i-th entry should contain the type of
   *                          boundary condition to impose at the upper
   *                          boundary in the i-th coordinate direction.
   *                          For information about the boundary
   *                          condition types, see documentation of
   *                          BoundaryConditionModule class.
   *  - lower_bc_ext (in):    vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the lower face of the computational
   *                          domain in each coordinate direction for
   *                          the extension field.
   *                          The i-th entry should contain the type of
   *                          boundary condition to impose at the lower
   *                          boundary in the i-th coordinate direction.
   *                          For information about the boundary
   *                          condition types, see documentation of
   *                          BoundaryConditionModule class.
   *  - upper_bc_ext (in):    vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the upper face of the computational
   *                          domain in each coordinate direction for
   *                          the extension field.
   *                          The i-th entry should contain the type of
   *                          boundary condition to impose at the upper
   *                          boundary in the i-th coordinate direction.
   *                          For information about the boundary
   *                          condition types, see documentation of
   *                          BoundaryConditionModule class.
   *  - max_iterations (in):  maximum number of iterations to use
   *                          for field extension.  Set max_iterations
   *                          to -1 to use the value specified in the
   *                          input file. (default value: -1)
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
virtual void computeExtensionField(
    const int phi_component,
    const hier::IntVector& lower_bc_phi,
    const hier::IntVector& upper_bc_phi,
    const hier::IntVector& lower_bc_ext,
    const hier::IntVector& upper_bc_ext,
    const int max_iterations=-1);

/*!
   * computeExtensionFieldForSingleComponent() extends the specified
   * component of the field off of the interface defined by the zero
   * level set of the function phi by advecting it in the direction normal
   * to the interface.
   *
   * Arguments:
   *  - component (in):       component to extend off of the zero
   *                          level set
   *  - phi_component (in):   component of level set function to
   *                          use in field extension calculation
   *  - lower_bc_phi (in):    vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the lower face of the computational
   *                          domain in each coordinate direction for
   *                          the level set function.
   *                          The i-th entry should contain the type of
   *                          boundary condition to impose at the lower
   *                          boundary in the i-th coordinate direction.
   *                          For information about the boundary
   *                          condition types, see documentation of
   *                          BoundaryConditionModule class.
   *  - upper_bc_phi (in):    vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the upper face of the computational
   *                          domain in each coordinate direction for
   *                          the level set function.
   *                          The i-th entry should contain the type of
   *                          boundary condition to impose at the upper
   *                          boundary in the i-th coordinate direction.
   *                          For information about the boundary
   *                          condition types, see documentation of
   *                          BoundaryConditionModule class.
   *  - lower_bc_ext (in):    vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the lower face of the computational
   *                          domain in each coordinate direction for
   *                          the extension field.
   *                          The i-th entry should contain the type of
   *                          boundary condition to impose at the lower
   *                          boundary in the i-th coordinate direction.
   *                          For information about the boundary
   *                          condition types, see documentation of
   *                          BoundaryConditionModule class.
   *  - upper_bc_ext (in):    vector of integers specifying the
   *                          type of boundary conditions to impose
   *                          on the upper face of the computational
   *                          domain in each coordinate direction for
   *                          the extension field.
   *                          The i-th entry should contain the type of
   *                          boundary condition to impose at the upper
   *                          boundary in the i-th coordinate direction.
   *                          For information about the boundary
   *                          condition types, see documentation of
   *                          BoundaryConditionModule class.
   *  - max_iterations (in):  maximum number of iterations to use
   *                          for field extension.  Set max_iterations
   *                          to -1 to use the value specified in the
   *                          input file. (default value: -1)
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
  virtual void computeExtensionFieldForSingleComponent(
    const int component,
    const int phi_component,
    const hier::IntVector& lower_bc_phi,
    const hier::IntVector& upper_bc_phi,
    const hier::IntVector& lower_bc_ext,
    const hier::IntVector& upper_bc_ext,
    const int max_iterations=-1);

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
   *  - hierarchy (in):       Boost pointer to new PatchHierarchy
   *  - coarsest_level (in):  coarsest level in the hierarchy to be updated
   *  - finest_level (in):    finest level in the hierarchy to be updated
   *
   * Return value:            none
   *
   */
  virtual void resetHierarchyConfiguration(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
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
   * advanceFieldExtensionEqnUsingTVDRK*() advances the advection
   * equation for extending fields off of the zero level set using
   * a first-, second-, or third-order TVD Runge-Kutta step.
   *
   * Arguments:
   *  - dt (in):                      time increment to advance the level set
   *                                  functions
   *  - field_component (in):         component of field to advance in time
   *  - phi_component (in):           phi component to use time advance
   *  - lower_bc_ext (in):            vector of integers specifying the
   *                                  type of boundary conditions to impose
   *                                  on the lower face of the computational
   *                                  domain in each coordinate direction for
   *                                  the extension field.
   *                                  The i-th entry should contain the
   *                                  type of boundary condition to
   *                                  impose at the lower boundary in
   *                                  the i-th coordinate direction.
   *                                  For information about the boundary
   *                                  condition types, see documentation of
   *                                  BoundaryConditionModule class.
   *  - upper_bc_ext (in):            vector of integers specifying the
   *                                  type of boundary conditions to impose
   *                                  on the upper face of the computational
   *                                  domain in each coordinate direction for
   *                                  the extension field.
   *                                  The i-th entry should contain the
   *                                  type of boundary condition to
   *                                  impose at the upper boundary in
   *                                  the i-th coordinate direction.
   *                                  For information about the boundary
   *                                  condition types, see documentation of
   *                                  BoundaryConditionModule class.
   *
   * Return value:                    none
   *
   * NOTES:
   *  - These methods are NOT intended to be used directly by the
   *    user.  They are only helper methods for the
   *    computeExtensionFieldForSingleComponent() and
   *    computeExtensionField() methods.
   *  - phi_component is provided in case the user provided phi data
   *    is used directly for computation (rather than using a scratch
   *    copy).
   *
   */
  virtual void advanceFieldExtensionEqnUsingTVDRK1(
    const LSMLIB_REAL dt,
    const int field_component,
    const int phi_component,
    const hier::IntVector& lower_bc_ext,
    const hier::IntVector& upper_bc_ext);
  virtual void advanceFieldExtensionEqnUsingTVDRK2(
    const LSMLIB_REAL dt,
    const int field_component,
    const int phi_component,
    const hier::IntVector& lower_bc_ext,
    const hier::IntVector& upper_bc_ext);
  virtual void advanceFieldExtensionEqnUsingTVDRK3(
    const LSMLIB_REAL dt,
    const int field_component,
    const int phi_component,
    const hier::IntVector& lower_bc_ext,
    const hier::IntVector& upper_bc_ext);

  /*!
   * computeFieldExtensionEqnRHS() computes the right-hand side of
   * the field extension equation when it is written in the form:
   *
   *   S_t = ...
   *
   * Arguments:
   *  - extension_field_handle (in):  PatchData handle for field data that
   *                                  should be used to compute spatial
   *                                  derivatives
   *  - phi_component (in):           phi component to use to compute RHS
   *
   * Return value:                    none
   *
   * NOTES:
   *  - phi_component is provided in case the user provided phi data
   *    is used directly for computation (rather than using a scratch
   *    copy).
   *
   */
  virtual void computeFieldExtensionEqnRHS(
    const int extension_field_handle,
    const int phi_component);

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
   * for the field extension calculation.
   *
   * Arguments:
   *  - phi_ghostcell_width (in):  ghostcell width for PatchData associated
   *                               with phi_handle
   *
   * Return value:                 none
   *
   */
  virtual void initializeVariables(const hier::IntVector& phi_ghostcell_width);

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
   * getFromInput() configures the FieldExtensionAlgorithm
   * object from the values in the specified input database.
   *
   * Arguments:
   *  - db (in):               input database from which to take parameter
   *                           values
   *
   * Return value:             none
   *
   * NOTES:
   *  - An assertion results if the database boost pointer is null.
   *
   */
  virtual void getFromInput(boost::shared_ptr<tbox::Database> db);

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

  // Boost pointer to PatchHierarchy object
  boost::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;

  // Boost pointer to GridGeometry
  boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

  /*
   * PatchData handles for data required to solve field extension equation
   */

  // data field handles
  int d_extension_field_handle;
  int d_phi_handle;
  int d_control_volume_handle;

  // scratch data
  vector<int> d_extension_field_scr_handles;
  hier::IntVector d_ext_field_scratch_ghostcell_width;
  int d_phi_scr_handle;
  hier::IntVector d_phi_scratch_ghostcell_width;
  int d_rhs_handle;
  int d_normal_vector_handle;
  int d_grad_field_handle;
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

  // field data parameters
  int d_num_field_components;

  // ComponentSelector to organize variables
  hier::ComponentSelector d_scratch_data;

  /*
   * Boundary condition objects
   */
  boost::shared_ptr<BoundaryConditionModule> d_phi_bc_module;
  boost::shared_ptr<BoundaryConditionModule> d_ext_field_bc_module;

  /*
   * Communication objects.
   */
  tbox::Array< boost::shared_ptr<xfer::RefineAlgorithm> >
      d_extension_field_fill_bdry_alg;
  tbox::Array< tbox::Array< boost::shared_ptr<xfer::RefineSchedule> > >
      d_extension_field_fill_bdry_sched;
  boost::shared_ptr<xfer::RefineAlgorithm> d_phi_fill_bdry_alg;
  tbox::Array< boost::shared_ptr<xfer::RefineSchedule> > d_phi_fill_bdry_sched;


private:

  /*
   * Private assignment operator to prevent use.
   *
   * Arguments:
   *  - rhs (in):    object to copy
   *
   * Return value:   *this
   *
   */
  const FieldExtensionAlgorithm& operator=(
    const FieldExtensionAlgorithm& rhs){
      return *this;
  }

};

} // end LSMLIB namespace

#endif
