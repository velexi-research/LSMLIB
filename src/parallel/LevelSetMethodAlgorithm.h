/*
 * File:        LevelSetMethodAlgorithm.h
 * Description: Header file for the main level set method algorithm class
 */

#ifndef included_LevelSetMethodAlgorithm_h
#define included_LevelSetMethodAlgorithm_h

/*! \class LSMLIB::LevelSetMethodAlgorithm
 *
 * \brief
 * The LevelSetMethodAlgorithm class provides a simple interface for using
 * the LevelSetMethod classes to capture the dynamics of implicit surfaces
 * and curves in one-, two- and three-dimensions.  It also provides support
 * for vector level set method calculations.
 *
 * It provides the flexibility of (a) managing the level set method
 * calculation with the standard LevelSetFunctionIntegrator and
 * LevelSetMethodGriddingAlgorithm classes or (b) using custom integrator
 * and gridding algorithm classes (by subclassing the
 * LevelSetFunctionIntegratorStrategy and
 * LevelSetMethodGriddingStrategy classes).
 *
 * For more details about the numerical algorithms used by the standard
 * LevelSetFunctionIntegrator class, see the documentation at the beginning
 * of the header file for that class.
 *
 *
 * <h3> NOTES </h3>
 *
 *  - AMR is NOT yet supported.  It is still currently under development.
 *
 *
 * <h3> USAGE </h3>
 * There are two ways to use the LevelSetMethodAlgorithm:
 *
 * (A) use the standard LevelSetFunctionIntegrator and
 *     LevelSetMethodGriddingAlgorithm classes or
 *
 * (B) use custom level set function integrator and level set method
 *     gridding algorithm classes that conform to the interfaces defined
 *     by the LevelSetFunctionIntegratorStrategy and
 *     LevelSetMethodGriddingStrategy classes.
 *
 *
 * <h4> Method (A) </h4>
 *
 * Method (A) should be sufficient for most applications.  To use (A),
 * an application developer should use the following steps:
 *
 * -# Provide a concrete implementation of the LevelSetMethodPatchStrategy
 *    class.  This class must provide numerical routines for computations
 *    on a single patch (i.e., logically rectangular, uniform grid).
 *    In particular, methods must be provided for initialization of
 *    data, setting the physical boundary conditions for the level set
 *    functions, and computing a stable time step size for the next
 *    TVD Runge-Kutta time step.  For more details, see the header file
 *    for the LevelSetMethodPatchStrategy class.
 * -# Provide a concrete implementation of the
 *    LevelSetMethodVelocityFieldStrategy class.  This class must
 *    provide routines associated with computing the velocity field
 *    that is used to advance the level set functions.  In particular,
 *    methods must be provided for computing the velocity field used
 *    in a level set method calculation, computing a stable time step
 *    size based on the velocity field computation, initializing data
 *    required for the velocity computation, and returning the PatchData
 *    handle for the velocity field data.  For more details, see the
 *    header file for the LevelSetMethodVelocityFieldStrategy class.
 * -# Create a LevelSetMethodAlgorithm object using the constructor
 *    that takes boost pointers to LevelSetMethodVelocityFieldStrategy and
 *    LevelSetMethodPatchStrategy objects as arguments.
 * -# Begin the level set method calculation by invoking the
 *    @ref initializeLevelSetMethodCalculation() method.  This will
 *    initialize the data on the PatchHierarchy.
 * -# Integrate the level set functions in time by invoking the
 *    @ref advanceLevelSetFunctions() method for each time step.
 *    Stable time step sizes can be computed for each time step using
 *    the @ref computeStableDt() method.
 * -# Access to level set function data and other integration data
 *    is available through several accessor method defined below.
 *
 * Users interested in the details of the level set method calculation
 * used by the LevelSetFunctionIntegrator should consult the class
 * definitions in the header files for the LevelSetFunctionIntegrator and
 * LevelSetMethodToolbox classes.
 *
 *
 * <h4> Method (B) </h4>
 *
 * To use (B), an application developer should use the following steps:
 *
 * -# Implement a level set method integrator class that conforms to
 *    the interface defined by the LevelSetFunctionIntegratorStrategy
 *    class.
 * -# Implement a level set gridding algorithm class that conforms to
 *    the interface defined by the LevelSetMethodGriddingStrategy
 *    class.
 * -# Create a LevelSetMethodAlgorithm object using the constructor
 *    that takes a boost pointer to a LevelSetFunctionIntegratorStrategy
 *    and a LevelSetMethodGriddingStrategy as arguments.
 * -# Follow steps 4 through 6 from Method (A).
 *
 *
 * <h3> User-specified parameters (input database field) </h3>
 *
 * When using the standard LevelSetFunctionIntegrator, it is possible
 * to modify the behavior of the integrator through an input file.
 * The input data parameters available for/required by the user are
 * described below. In the input file, the LevelSetFunctionIntegrator
 * database and LevelSetMethodGriddingAlgorithm database must be organized
 * as separate sub-databases of a single, high-level LevelSetMethodAlgorithm
 * database.  That is, the input file should have the following form:
 *
 * <pre>
 * LevelSetMethodAlgorithm {
 *
 *   LevelSetFunctionIntegrator {
 *     ... input parameters for LevelSetFunctionIntegrator ...
 *   }
 *
 *   LevelSetMethodGriddingAlgorithm {
 *     ... input parameters for LevelSetMethodGriddingAlgorithm ...
 *   }
 *
 * } // end input database for LevelSetMethodAlgorithm
 * </pre>
 *
 *
 * <h4> LevelSetFunctionIntegrator Input Database Parameters </h4>
 *
 *   <h5> Level Set Method Parameters: </h5>
 *
 *   - start_time                  = start time for calculation
 *                                   (default = 0.0)
 *   - end_time                    = end time for calculation
 *                                   (default = 10.0)
 *   - cfl_number                  = CFL number (default = 0.5)
 *   - spatial_derivative_type     = type of spatial derivative
 *                                   (default = "WENO")
 *   - spatial_derivative_order    = order of spatial derivative
 *                                   (default = 5)
 *   - tvd_runge_kutta_order       = order of Runge-Kutta time integration
 *                                   (default = 3)
 *   - reinitialization_interval   = interval between reinitialization
 *                                   (default = 10)
 *                                   (reinitialization disabled if <= 0)
 *   - reinitialization_stop_tol   = stopping criterion for termination of
 *                                   evolution of reinitialization equation.
 *                                   Reinitialization stops when the max norm
 *                                   of the change in the level set function
 *                                   drops below the specified tolerance.
 *   - reinitialization_stop_dist  = approximate stopping criterion for
 *                                   reinitialization of level set functions.
 *                                   Reinitialization terminates after the
 *                                   information from the zero level set has
 *                                   propogated by approximately the specified
 *                                   distance.
 *   - reinitialization_max_iters  = maximum number of time steps to take
 *                                   during the reinitialization process
 *                                   (default = 25)
 *   - orthogonalization_interval  = interval between orthogonalizing phi
 *                                   and psi for codimension-two problems
 *                                   (default = 10)
 *                                   (orthogonalization disabled if <= 0)
 *   - orthogonalization_stop_tol  = stopping criterion for termination of
 *                                   evolution of orthogonalization equation.
 *                                   Orthogonalization stops when the max norm
 *                                   of the change in the level set function
 *                                   drops below the specified tolerance.
 *   - orthogonalization_stop_dist = approximate stopping criterion for
 *                                   orthogonalization of level set functions.
 *                                   Orthogonalization terminates after the
 *                                   information from the zero level set has
 *                                   propogated by approximately the specified
 *                                   distance.
 *   - orthogonalization_max_iters = maximum number of time steps to take
 *                                   during the orthogonalization process
 *                                   (default = 25)
 *
 *   <h5> Boundary Condition Parameters: </h5>
 *
 *   - lower_bc_phi_[i]            = boundary conditions for the lower
 *                                   face of the computational domain in
 *                                   each coordinate direction for the i-th
 *                                   component of the PHI level set function.
 *                                   lower_bc_phi_[i] is a vector of length
 *                                   DIM.  The j-th entry should contain the
 *                                   type of boundary condition to impose at
 *                                   the lower boundary in the j-th coordinate
 *                                   direction.  See NOTES ON INPUT PARAMETERS
 *                                   section for boundary condition types.
 *                                   For information about the boundary
 *                                   condition types, see documentation
 *                                   of the BoundaryConditionModule class.
 *                                   (default = vector of zeros)
 *   - upper_bc_phi_[i]            = boundary conditions for the upper
 *                                   face of the computational domain in
 *                                   each coordinate direction for the i-th
 *                                   component of the PHI level set function.
 *                                   upper_bc_phi_[i] is a vector of length
 *                                   DIM.  The j-th entry should contain the
 *                                   type of boundary condition to impose at
 *                                   the upper boundary in the j-th coordinate
 *                                   direction.  See NOTES ON INPUT PARAMETERS
 *                                   section for boundary condition types.
 *                                   For information about the boundary
 *                                   condition types, see documentation
 *                                   of the BoundaryConditionModule class.
 *                                   (default = vector of zeros)
 *   - lower_bc_psi_[i]            = boundary conditions for the lower
 *                                   face of the computational domain in
 *                                   each coordinate direction for the i-th
 *                                   component of the PSI level set function.
 *                                   lower_bc_psi_[i] is a vector of length
 *                                   DIM.  The j-th entry should contain the
 *                                   type of boundary condition to impose at
 *                                   the lower boundary in the j-th coordinate
 *                                   direction.  See NOTES ON INPUT PARAMETERS
 *                                   section for boundary condition types.
 *                                   For information about the boundary
 *                                   condition types, see documentation
 *                                   of the BoundaryConditionModule class.
 *                                   (default = vector of zeros)
 *   - upper_bc_psi_[i]            = boundary conditions for the upper
 *                                   face of the computational domain in
 *                                   each coordinate direction for the i-th
 *                                   component of the PSI level set function.
 *                                   upper_bc_psi_[i] is a vector of length
 *                                   DIM.  The j-th entry should contain the
 *                                   type of boundary condition to impose at
 *                                   the upper boundary in the j-th coordinate
 *                                   direction.  See NOTES ON INPUT PARAMETERS
 *                                   section for boundary condition types.
 *                                   For information about the boundary
 *                                   condition types, see documentation
 *                                   of the BoundaryConditionModule class.
 *                                   (default = vector of zeros)
 *
 *   <h5> Miscellaneous Parameters: </h5>
 *
 *   - verbose_mode                = TRUE if status should be output during
 *                                   integration (default = TRUE)
 *
 *   <h5> AMR Parameters (CURRENTLY UNUSED): </h5>
 *
 *   - use_AMR                     = TRUE if AMR should be used
 *                                   (default = FALSE)
 *   - regrid_interval             = regridding interval (default = 5)
 *   - tag_buffer_width            = number of buffer cells to use around
 *                                   cells tagged for refinement
 *                                   (default = 2)
 *   - refinement_cutoff_value     = cutoff value for distance function
 *                                   (default = 1.0)
 *
 *
 * <h4> LevelSetMethodGriddingAlgorithm Input Database Parameters </h4>
 *
 *   <h5> Hierarchy Structure Input: </h5>
 *
 *   - max_levels (REQUIRED)          =  integer value specifying maximum
 *                                       number of levels allowed in the AMR
 *                                       PatchHierarchy.
 *   - largest_patch_size  (REQUIRED) =  an array of integer vectors (each has
 *                                       length = DIM) that specify the
 *                                       dimensions of largest patch allowed
 *                                       on each level in the hierarchy.  If
 *                                       more than max_levels entries are
 *                                       given, extra entries will be ignored.
 *                                       If fewer than max_levels entries are
 *                                       given, then the last element in the
 *                                       array will be used on each level
 *                                       without a specified input value.
 *   - ratio_to_coarser (REQUIRED)    =  set of (max_levels - 1) integer
 *                                       vectors, each of which indicates the
 *                                       ratio of the index space of a
 *                                       PatchLevel to that of the next
 *                                       coarser level.  The input for each
 *                                       level must correspond to the format
 *                                       ``level_n = vector'', where n is the
 *                                       level number and each vector must
 *                                       have length DIM.
 *
 *   <h5> Adaptive Refinement Input (CURRENTLY UNUSED): </h5>
 *
 *   - tagging_method (OPTIONAL)      =  string array specification of the type
 *                                       of cell-tagging used. Valid choices
 *                                       include: ``GRADIENT_DETECTOR'' and
 *                                       ``REFINE_BOXES''.  A combination of
 *                                       any or all of the above may be placed
 *                                       in any order. If no input is given, no
 *                                       tagging will be performed.
 *   - RefineBoxes (OPTIONAL)         =  input section describing the refine
 *                                       boxes for each level.
 *     - Level<ln>                    =  input section provides the sequence of
 *                                       Box arrays describing where
 *                                       user-specified refinement is to occur
 *                                       on Level ln.
 *       - times (OPTIONAL)           =  LSMLIB_REAL array specifying times at which
 *                                       a particular box sequence is to be
 *                                       used.
 *       - cycles (OPTIONAL)          =  integer array specifying regrid cycles
 *                                       at which a particular box seqence is
 *                                       to be used.
 *       - boxes_0                    =  box array specifying refine boxes for
 *                                       sequence 0.
 *       - boxes_1                    =  box array specifying refine boxes for
 *                                       sequence 1.
 *       - boxes_n                    =  box array specifying refine boxes for
 *                                       sequence n.
 *
 *   <h5> Load Balancer Input: </h5>
 *
 *   - NO REQUIRED INPUT PARAMETERS (several OPTIONAL input parameters)
 *
 *
 *  <h4> NOTES ON INPUT PARAMETERS </h4>
 *
 *    - When restarting a computation, the following input parameters
 *      override the values from the restart file:
 *        end_time,
 *        reinitialization_interval,
 *        reinitialization_stop_tol,
 *        reinitialization_stop_dist,
 *        reinitialization_max_iters,
 *        orthogonalization_interval,
 *        orthogonalization_stop_tol,
 *        orthogonalization_stop_dist,
 *        orthogonalization_max_iters,
 *        verbose
 *
 *    - The precedence of the three input parameters for reinitialization
 *      and orthogonalization is as follows:
 *      -# if provided, the stop tolerance is used with a maximum
 *         number of iterations determined from the calculation in (b)
 *         or a default value of 1000 iterations (to avoid failure to
 *         terminate)
 *      -# when both the stop distance and maximum number of iterations
 *         are specified, the one that results in the smaller number of
 *         iterations is used.  when only one is specified, it is the
 *         sole parameter used to compute the maximum number of time steps
 *         to take when advancing the reinitialization/orthogonalization
 *         equation
 *      -# when no stopping criteria are supplied, the maximum number of
 *         time steps taken defaults to 25.
 *
 *    - The names of the sub-databases, LevelSetFunctionIntegrator and
 *      LevelSetMethodGriddingAlgorithm, MUST be named exactly as above.
 *      The name of the top-level database, LevelSetMethodAlgorithm, may
 *      be arbitrary (it is up to the user to make sure that the correct
 *      top-level input database is passed to the constructor).
 *
 *    - The numbering for the boundary conditions begins at 1.
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
 *    - When using a custom level set function integrator, all input
 *      parameters are ignored.
 *
 *    - The descriptions of the input parameters for the
 *      LevelSetGriddingAlgorithm were taken almost verbatim from the class
 *      descriptions of the SAMRAI::mesh::GriddingAlgorithm,
 *      SAMRAI::mesh::StandardTagAndInitialize, and
 *      SAMRAI::mesh::TagAndInitializeStrategy classes in the files
 *      GriddingAlgorithm.h, StandardTagAndInitialize.h, and
 *      TagAndInitializeStrategy.h.   For more details about the input
 *      parameters, please consult the SAMRAI documentation for these
 *      classes.
 *
 *    - For a list and description of optional gridding algorithm input
 *      fields, see the documentation for the SAMRAI::mesh::GriddingAlgorithm
 *      class.
 *
 *    - For a list and description of optional load balancer input
 *      fields, see the documentation for the SAMRAI::mesh::LoadBalancer
 *      class.
 *
 * <h4> Sample Input File </h4>
 *
 *  <pre>
 *  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *  LevelSetMethodAlgorithm{
 *
 *    LevelSetFunctionIntegrator {
 *      start_time  = 0.0
 *      end_time    = 0.5
 *
 *      cfl_number               = 0.5
 *      spatial_derivative_type  = "WENO"
 *      spatial_derivative_order = 5
 *      tvd_runge_kutta_order    = 3
 *
 *      reinitialization_interval = 0
 *      reinitialization_max_iters = 20
 *      reinitialization_stop_dist = 0.2
 *      orthogonalization_interval = 0
 *      orthogonalization_max_iters = 20
 *      orthogonalization_stop_dist = 0.2
 *
 *      lower_bc_phi_0 = 1, 1, 1
 *      upper_bc_phi_0 = 1, 1, 1
 *
 *      use_AMR = FALSE
 *      refinement_cutoff_value = 0.25
 *      tag_buffer = 2,2,2,2,2,2
 *
 *      verbose = false
 *
 *    } // end of LevelSetFunctionIntegrator database
 *
 *
 *    LevelSetMethodGriddingAlgorithm {
 *      max_levels = 4
 *
 *      ratio_to_coarser {
 *         level_1            = 2, 2
 *         level_2            = 2, 2
 *         level_3            = 2, 2
 *      }
 *
 *      largest_patch_size {
 *        level_0 = 50,50
 *        level_1 = 100,100
 *        // all finer levels will use same values as level_1...
 *      }
 *
 *      tagging_method = "GRADIENT_DETECTOR","REFINE_BOXES"
 *
 *      RefineBoxes {
 *        level_0 = [(15,0),(29,14)]
 *        level_1 = [(65,10),(114,40)]
 *      }
 *
 *      LoadBalancer {
 *        // load balancer input parameters
 *      }
 *
 *    } // end of LevelSetMethodGriddingAlgorithm database
 *
 *  } // end of LevelSetMethodAlgorithm database
 *
 *  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *  </pre>
 *
 */

// Standard library headers
#include <ostream>

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/tbox/Array.h"

// LSMLIB header
#include "LSMLIB_config.h"
#include "LevelSetMethodToolbox.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations
namespace LSMLIB { class FieldExtensionAlgorithm; }
namespace LSMLIB { class LevelSetFunctionIntegratorStrategy; }
namespace LSMLIB { class LevelSetMethodGriddingStrategy; }
namespace LSMLIB { class LevelSetMethodPatchStrategy; }
namespace LSMLIB { class LevelSetMethodVelocityFieldStrategy; }
namespace SAMRAI { namespace hier { class IntVector; } }
namespace SAMRAI { namespace hier { class PatchHierarchy; } }
namespace SAMRAI { namespace tbox { class Database; } }


/******************************************************************
 *
 * LevelSetMethodAlgorithm Class Definition
 *
 ******************************************************************/

namespace LSMLIB {

class LevelSetMethodAlgorithm
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
   * This constructor for LevelSetMethodAlgorithm creates a
   * LevelSetFunctionIntegrator object from the user-specified
   * LevelSetMethodPatchStrategy and LevelSetMethodVelocityFieldStrategy
   * objects and the parameters in the input and restart databases.
   * This LevelSetFunctionIntegrator object manages the time integration
   * of the level set functions.
   *
   * Arguments:
   *  - lsm_algorithm_input_db (in):        input database containing
   *                                        user-defined parameters for the
   *                                        LevelSetMethodAlgorithm
   *  - patch_hierarchy (in):               PatchHierarchy for computation
   *  - patch_strategy (in):                LevelSetMethodPatchStrategy
   *  - velocity_field_strategy (in):       LevelSetMethodVelocityFieldStrategy
   *  - num_level_set_fcn_components (in):  number of components of level set
   *                                        functions (default = 1)
   *  - codimension (in):                   codimension of problem
   *                                        (default = 1)
   *  - object_name (in):                   name for object (default =
   *                                        "LevelSetMethodAlgorithm")
   *
   */
  LevelSetMethodAlgorithm(
    boost::shared_ptr<tbox::Database> lsm_algorithm_input_db,
    boost::shared_ptr<hier::PatchHierarchy > patch_hierarchy,
    LevelSetMethodPatchStrategy* patch_strategy,
    LevelSetMethodVelocityFieldStrategy* velocity_field_strategy,
    const int num_level_set_fcn_components = 1,
    const int codimension = 1,
    const std::string& object_name = "LevelSetMethodAlgorithm");

  /*!
   * This constructor for LevelSetMethodAlgorithm uses a user-specified
   * concrete subclass of the LevelSetFunctionIntegratorStrategy class to
   * manage the time integration of the level set functions.
   *
   * Arguments:
   *  - lsm_integrator_strategy (in): boost pointer to concrete subclass of the
   *                                   LevelSetFunctionIntegratorStrategy
   *  - lsm_gridding_strategy (in):   boost pointer to concrete subclass of the
   *                                   LevelSetMethodGriddingStrategy
   *  - object_name (in):              name for object
   *
   */
  LevelSetMethodAlgorithm(
    boost::shared_ptr<LevelSetFunctionIntegratorStrategy> lsm_integrator_strategy,
    boost::shared_ptr<LevelSetMethodGriddingStrategy> lsm_gridding_strategy,
    const std::string& object_name = "LevelSetMethodAlgorithm");

  /*!
   * The destructor for LevelSetMethodAlgorithm does nothing.
   */
  virtual ~LevelSetMethodAlgorithm();

  //! @}


  //! @{
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
   *  - When the standard constructor for the LevelSetMethodAlgorithm
   *    is used, the PatchData for phi associated with the returned
   *    PatchData handle has the number of ghostcells required for the
   *    spatial derivative type and order specified when the
   *    LevelSetMethodAlgorithm object was constructed.
   *
   */
  virtual int getPhiPatchDataHandle() const;

  /*!
   * getPsiPatchDataHandle() returns the patch data handle for psi.
   *
   * Arguments:     none
   *
   * Return value:  PatchData handle for psi
   *
   * NOTES:
   *  - When the standard constructor for the LevelSetMethodAlgorithm
   *    is used, the PatchData for psi associated with the returned
   *    PatchData handle has the number of ghostcells required for the
   *    spatial derivative type and order specified when the
   *    LevelSetMethodAlgorithm object was constructed.
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
   *                the end time for the integration; false otherwise
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
   * an instance of the LevelSetMethodAlgorithm class.
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
   * @ Accessor methods for simulation parameters
   *
   ****************************************************************/

  /*!
   * getSpatialDerivativeType() returns type of spatial derivative
   * discretization used by the standard LevelSetFunctionIntegrator object.
   *
   * Arguments:      none
   *
   * Return value :  spatial derivative type
   *
   * NOTES:
   *  - Meaning of return values: 0 = ENO, 1 = WENO.
   *
   */
  virtual int getSpatialDerivativeType() const;

  /*!
   * getSpatialDerivativeOrder() returns order of the spatial derivative
   * discretization used by the standard LevelSetFunctionIntegrator object.
   *
   * Arguments:      none
   *
   * Return value :  spatial derivative order
   *
   */
  virtual int getSpatialDerivativeOrder() const;

  /*!
   * getTVDRungeKuttaOrder() returns order of the TVD Runge-Kutta
   * integration scheme used by the standard LevelSetFunctionIntegrator
   * object.
   *
   * Arguments:      none
   *
   * Return value :  TVD Runge-Kutta order
   *
   */
  virtual int getTVDRungeKuttaOrder() const;

  /*!
   * getCFLNumber() returns CFL number used by the standard
   * LevelSetFunctionIntegrator object to compute a stable time step.
   *
   * Arguments:      none
   *
   * Return value :  CFL number
   *
   */
  virtual LSMLIB_REAL getCFLNumber() const;

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @ Methods for setting boundary conditions
   *
   ****************************************************************/

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
   *    NONE, HOMOEGENEOUS_NEUMANN, LINEAR_EXTRAPOLATION,
   *    SIGNED_LINEAR_EXTRAPOLATION, ANTI_PERIODIC
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
    const hier::IntVector& lower_bc,
    const hier::IntVector& upper_bc,
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int component = -1);

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for time advancing level set functions
   *
   ****************************************************************/

  /*!
   * initializeLevelSetMethodCalculation() initializes the
   * data on the PatchHierarchy in preparation for a level
   * set method calculation.
   *
   * Arguments:     none
   *
   * Return value:  none
   *
   */
  virtual void initializeLevelSetMethodCalculation();

  /*!
   * computeStableDt() computes the maximum allowable dt for the
   * next time step of the level set functions.
   *
   * If the standard LevelSetFunctionIntegrator is used, the
   * maximum stable dt is computed using the algorithm:
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
   ****************************************************************
   *
   * @name Methods for managing the grid configuration
   *
   ****************************************************************/

  /*!
   * resetHierarchyConfiguration() resets the configuration of the
   * calculation and communication objects to be consistent with
   * the specified PatchHierarchy.
   *
   * Arguments:
   *  - hierarchy (in):       PatchHierarchy to reconfigure
   *  - coarsest_level (in):  coarsest level in hierarchy to reconfigure
   *  - finest_level (in):    finest level in hierarchy to reconfigure
   *
   * Return value:            none
   *
   */
  void resetHierarchyConfiguration(
    const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int coarsest_level,
    const int finest_level);

  /*!
   * regridPatchHierarchy() regrids the entire PatchHierarchy
   * and reinitializes the data on the PatchHierarchy using
   * interpolation and averaging, as necessary.
   *
   * Arguments:      none
   *
   * Return value:   none
   *
   */
  virtual void regridPatchHierarchy();

  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for creating FieldExtensionAlgorithm objects
   *
   ****************************************************************/

  /*!
   * getFieldExtensionAlgorithm() creates a FieldExtensionAlgorithm
   * object that can be used to extend the specified field off of the
   * zero level set of the specified level set function (PHI or PSI).
   * The PatchHierarchy used for by the FieldExtensionAlgorithm is the
   * same as the one that is used by the LevelSetMethodAlgorithm.
   * This version of getFieldExtensionAlgorithm() read the parameters
   * for the field extension calculation from specified input database.
   *
   * Arguments:
   *  - input_db (in):                  input database containing user-defined
   *                                    parameters
   *  - field_handle (in):              PatchData handle for field to
   *                                    extend off of the interface
   *  - level_set_fcn (in):             the level set function to use
   *                                    in field extension calculation
   *  - verbose_mode (in):              flag to activate/deactivate
   *                                    verbose-mode (default = false)
   *  - object_name (in):               string name for object (default =
   *                                    "FieldExtensionAlgorithm")
   *
   * Return value:                      none
   *
   * NOTES:
   *  - getFieldExtensionAlgorithm() may only be used when using the
   *    standard LevelSetFunctionIntegrator.
   *
   *  - FieldExtensionAlgorithm objects creted using this method do
   *    NOT need to have the
   *    FieldExtensionAlgorithm::resetHierarchyConfiguration() method
   *    invoked when the hierarchy configuration changes.  The
   *    LevelSetMethodAlgorithm automatically handles the reset when
   *    LevelSetMethodAlgorithm::resetHierarchyConfiguration() is invoked.
   *
   */
  virtual boost::shared_ptr<FieldExtensionAlgorithm>
    getFieldExtensionAlgorithm(
      boost::shared_ptr<tbox::Database> input_db,
      const int field_handle,
      const LEVEL_SET_FCN_TYPE level_set_fcn,
      const std::string& object_name = "FieldExtensionAlgorithm");

  /*!
   * getFieldExtensionAlgorithm() creates a FieldExtensionAlgorithm
   * object that can be used to extend the specified field off of the
   * zero level set of the specified level set function (PHI or PSI).
   * The PatchHierarchy used for by the FieldExtensionAlgorithm is the
   * same as the one that is used by the LevelSetMethodAlgorithm.
   * This version of getFieldExtensionAlgorithm() takes the parameters
   * for the field extension calculation as explicit arguments.
   *
   * Arguments:
   *  - field_handle (in):              PatchData handle for field to
   *                                    extend off of the interface
   *  - level_set_fcn (in):             the level set function to use
   *                                    in field extension calculation
   *  - spatial_derivative_type (in):   type of spatial derivative calculation
   *                                    (default = spatial derivative type
   *                                     used for evolving level set function)
   *  - spatial_derivative_order (in):  order of spatial derivative
   *                                    (default = spatial derivative order
   *                                     used for evolving level set function)
   *  - tvd_runge_kutta_order (in):     order of TVD Runge-Kutta integration
   *                                    (default = TVD Runge-Kutta order
   *                                     used for evolving level set function)
   *  - cfl_number (in):                CFL number used to compute dt
   *                                    (default = CFL number used for
   *                                     evolving level set function)
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
   *
   * Return value:                      none
   *
   * NOTES:
   *  - getFieldExtensionAlgorithm() may only be used when using the
   *    standard LevelSetFunctionIntegrator.
   *
   *  - FieldExtensionAlgorithm objects created using this method do
   *    NOT need to have the
   *    FieldExtensionAlgorithm::resetHierarchyConfiguration() method
   *    invoked when the hierarchy configuration changes.  The
   *    LevelSetMethodAlgorithm automatically handles the reset when
   *    LevelSetMethodAlgorithm::resetHierarchyConfiguration() is invoked.
   *
   */
  virtual boost::shared_ptr<FieldExtensionAlgorithm> getFieldExtensionAlgorithm(
    const int field_handle,
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    SPATIAL_DERIVATIVE_TYPE spatial_derivative_type = UNKNOWN,
    int spatial_derivative_order = 0,
    int tvd_runge_kutta_order = 0,
    LSMLIB_REAL cfl_number = 0,
    const LSMLIB_REAL stop_distance = 0.0,
    const int max_iterations = 0,
    const LSMLIB_REAL iteration_stop_tolerance = 0.0,
    const bool verbose_mode = false,
    const std::string& object_name = "FieldExtensionAlgorithm");

  //! @}

protected:

  /****************************************************************
   *
   * Data Members
   *
   ****************************************************************/

  // object name
  std::string d_object_name;

  // boost pointers to LevelSetFunctionIntegratorStrategy and
  // LevelSetMethodGriddingStrategy objects
  boost::shared_ptr<LevelSetFunctionIntegratorStrategy> d_lsm_integrator_strategy;
  boost::shared_ptr<LevelSetMethodGriddingStrategy> d_lsm_gridding_strategy;


  // lists of FieldExtensionAlgorithm objects
  // ( used in resetHierarchyConfiguration() )
  tbox::Array< boost::shared_ptr<FieldExtensionAlgorithm> >
      d_field_extension_alg_list;


  // simulation parameters when using standard LevelSetFunctionIntegrator
  // class
  boost::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
  bool d_using_standard_level_set_fcn_integrator;
  SPATIAL_DERIVATIVE_TYPE d_spatial_derivative_type;
  int d_spatial_derivative_order;
  int d_tvd_runge_kutta_order;
  LSMLIB_REAL d_cfl_number;

private:

  /*
   * Private copy constructor to prevent use.
   *
   * Arguments:
   *  - rhs (in):  LevelSetMethodAlgorithm object to copy
   *
   */
  LevelSetMethodAlgorithm(const LevelSetMethodAlgorithm& rhs){}

  /*
   * Private assignment operator to prevent use.
   *
   * Arguments:
   *  - rhs (in):    LevelSetMethodAlgorithm to copy
   *
   * Return value:   *this
   *
   */
  const LevelSetMethodAlgorithm& operator=(
    const LevelSetMethodAlgorithm& rhs){
      return *this;
  }

};

} // end LSMLIB namespace


#endif
