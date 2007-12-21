/*
 * File:        LevelSetFunctionIntegratorStrategy.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.15 $
 * Modified:    $Date: 2006/10/05 15:03:44 $
 * Description: Header file for level set method integrator strategy class
 */
 
#ifndef included_LevelSetFunctionIntegratorStrategy_h
#define included_LevelSetFunctionIntegratorStrategy_h

/*! \class LSMLIB::LevelSetFunctionIntegratorStrategy
 *
 * \brief 
 * The LevelSetFunctionIntegratorStrategy class defines the interface
 * for the time integration of the level set functions in the level set
 * method.  
 * 
 * In particular, it defines interfaces for the time integration 
 * of level set functions, initialization of data for the level set 
 * method calculation, choice of numerical schemes for computing spatial 
 * derivatives, and access to level set function data and status 
 * information (e.g. current time).
 * 
 */


#include "SAMRAI_config.h"     
#include "BasePatchHierarchy.h"
#include "PatchHierarchy.h"
#include "StandardTagAndInitStrategy.h"
#include "tbox/Pointer.h"

#include "LSMLIB_config.h"     
#include "LevelSetMethodToolbox.h"
#include "BoundaryConditionModule.h"

// namespaces
using namespace std;
using namespace SAMRAI;
using namespace hier;
using namespace mesh;
using namespace tbox;


/******************************************************************
 *
 * LevelSetFunctionIntegratorStrategy Class Definition
 *
 ******************************************************************/

namespace LSMLIB {

template<int DIM> class LevelSetFunctionIntegratorStrategy:
  public StandardTagAndInitStrategy<DIM>
{
public:

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
   */
  virtual int getPhiPatchDataHandle() const = 0;

  /*!
   * getPsiPatchDataHandle() returns the patch data handle for psi.
   * 
   * Arguments:     none
   * 
   * Return value:  PatchData handle for psi
   *
   */
  virtual int getPsiPatchDataHandle() const = 0;

  /*!
   * getControlVolumePatchDataHandle() returns the patch data handle for 
   * the control volume.
   * 
   * Arguments:     none
   * 
   * Return value:  PatchData handle for control volume
   *
   */
  virtual int getControlVolumePatchDataHandle() const = 0;

  //! @}


  //! @{
  /*!
   *******************************************************************
   *
   *  @name Accessor methods for simulation state
   *
   *******************************************************************/

  /*!
   * getStartTime() returns the start time value for the time integration.
   * 
   * Arguments:      none
   * 
   * Return value :  start time
   *
   */
  virtual LSMLIB_REAL getStartTime() const = 0;

  /*!
   * getEndTime() returns the end time value for the time integration.
   * 
   * Arguments:      none
   * 
   * Return value :  end time
   *
   */
  virtual LSMLIB_REAL getEndTime() const = 0;

  /*!
   * getCurrentTime() returns the current time value in the time integration.
   * 
   * Arguments:      none
   * 
   * Return value :  current time
   *
   */
  virtual LSMLIB_REAL getCurrentTime() const = 0;

  /*!
   * endTimeReached() returns true if the end time for the integration has 
   * been reached; otherwise, it returns false.
   * 
   * Arguments:     none
   * 
   * Return value:  true if the current time is equal to or after 
   *                the end time for the integration; false otherwise
   *
   */
  virtual bool endTimeReached() const = 0;

  /*!
   * numIntegrationStepsTaken() returns the number of integration steps 
   * that have been taken during the level set method calculation.
   *
   * Arguments:      none
   *
   * Return value :  current integration step count
   *
   */
  virtual int numIntegrationStepsTaken() const = 0;

  //! @}


  //! @{
  /*!
   ********************************************************************
   *
   * @name Methods for managing boundary conditions
   *
   ********************************************************************/

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
    const IntVector<DIM>& lower_bc,
    const IntVector<DIM>& upper_bc,
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int component = -1) = 0;

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
   * next time step of the level set functions.  
   *
   * Arguments:     none
   *
   * Return value:  maximum stable dt for next time step
   *
   */
  virtual LSMLIB_REAL computeStableDt() = 0;

  /*!
   * advanceLevelSetFunctions() advances the level set function
   * phi (and psi for codimension-two problems) by the specified
   * time increment, dt.
   * 
   * Arguments: 
   *  - dt (in):    time increment to advance the level set functions
   *
   * Return value:  true if patch hierarchy needs to be regridded after 
   *                this time step; false otherwise.
   *
   */
  virtual bool advanceLevelSetFunctions(const LSMLIB_REAL dt) = 0;

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
  virtual int getReinitializationInterval() const = 0;

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
  virtual void setReinitializationInterval(const int interval) = 0;

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
   *  - If max_iterations is set to a non-negative value, it should override
   *    ALL of the stopping criteria specified in the input file.
   *
   */
  virtual void reinitializeLevelSetFunctions(
    const LEVEL_SET_FCN_TYPE level_set_fcn = LSMLIB::PHI,
    const int max_iterations = -1) = 0;

  /*!
   * getOrthogonalizationInterval() should return the number of time 
   * steps between orthogonalizations of the level set functions.
   *
   * Arguments:     none
   *
   * Return value:  orthogonalization interval
   *
   */
  virtual int getOrthogonalizationInterval() const = 0;

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
  virtual void setOrthogonalizationInterval(const int interval) = 0;

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
    const int max_ortho_iterations = -1) = 0;
  
  //! @}


  //! @{
  /*!
   ****************************************************************
   *
   * @name Methods for managing grid configuration 
   *
   * NOTES:
   *  - initializeLevelData(), applyGradientDetector(), and
   *    resetHierarchyConfiguration() override virtual methods
   *    declared in the StandardTagAndInitStrategy base class.
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
    const Pointer< BasePatchHierarchy<DIM> > hierarchy ,
    const int level_number ,
    const double init_data_time ,
    const bool can_be_refined ,
    const bool initial_time ,
    const Pointer< BasePatchLevel<DIM> > old_level
      = Pointer< BasePatchLevel<DIM> >((0)) ,
    const bool allocate_data = true ) = 0;

  /*!
   * applyGradientDetector() sets integer tags to "1" in cells where 
   * refinement of the given level should occur according to the criteria 
   * that the absolute value of the level set functions, phi and psi, is 
   * less than some user-supplied threshold.
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
      const Pointer< BasePatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too) = 0;

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
  virtual void resetHierarchyConfiguration (
    Pointer< BasePatchHierarchy<DIM> > hierarchy ,
    int coarsest_level ,
    int finest_level ) = 0;

  //! @}


  //! @{
  /*!
   *******************************************************************
   *
   * @name Methods for pre-/post-processing level set function data 
   *
   *******************************************************************/

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
   * NOTES:
   *  - phi_handle and psi_handle are intended to be set to be handles
   *    to scratch data that possess a sufficient number of ghostcells 
   *    to compute the initial velocity field 
   *  - for codimension-one problems, psi_handle may be set arbitrarily 
   *    because psi has no meaning for codimension-one problems.
   *
   */
  virtual void preprocessInitializeVelocityField(
    int& phi_handle,
    int& psi_handle,
    const Pointer< PatchHierarchy<DIM> > hierarchy,
    const int level_number) = 0;

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
    const Pointer< PatchHierarchy<DIM> > hierarchy,
    const int level_number) = 0;

  //! @}


  //! @{
  /*!
   ****************************************************************
   * 
   * @name Destructors
   *
   ****************************************************************/

  /*!
   * Destructor does nothing but must be declared virtual so that 
   * concrete subclasses are properly destroyed.
   */
  virtual ~LevelSetFunctionIntegratorStrategy(){}

  //! @}

};

} // end LSMLIB namespace

#endif
