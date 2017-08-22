/*
 * File:        LevelSetMethodPatchStrategy.h
 * Description: Interface for single patch numerical routines for the
 *              level set algorithm
 */
 
#ifndef included_LevelSetMethodPatchStrategy_h
#define included_LevelSetMethodPatchStrategy_h

/*! \class LSMLIB::LevelSetMethodPatchStrategy
 *
 * \brief 
 * The LevelSetMethodPatchStrategy class defines the interface
 * for single patch, problem specific numerical routines. 
 *
 * In particular, it defines interfaces for initialization of data
 * on a patch, setting the boundary conditions for the level set
 * functions, and computing a stable time step size for the next
 * TVD Runge-Kutta time step.
 *
 * It is worth mentioning that the dt returned by the 
 * computeStableDtOnPatch() method is intended to give the user
 * the flexibility to specify a time step size in a completely 
 * arbitrary manner that is appropriate for his/her specific 
 * problem.
 *
 *
 * <h3> NOTES: </h3>
 *  - In order to use the LevelSetFunctionIntegrator class, a user MUST
 *    implement a concrete subclass of this ``strategy'' class.
 *
 */

// Standard library headers
#include <cfloat>  // IWYU pragma: keep

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep

// LSMLIB headers
#include "LSMLIB_config.h"

// Class/type declarations
namespace LSMLIB { class LevelSetFunctionIntegrator; }
namespace LSMLIB { class LevelSetMethodVelocityFieldStrategy; }
namespace SAMRAI { namespace hier { class IntVector; } }
namespace SAMRAI { namespace hier { class Patch; } }

// Namespaces
using namespace SAMRAI;

/******************************************************************
 *
 * LevelSetMethodPatchStrategy Class Definition
 *
 ******************************************************************/

namespace LSMLIB {

class LevelSetMethodPatchStrategy
{
public:

  //! @{
  /*!
 / ****************************************************************
   * 
   * @name Methods for setting initial and boundary conditions
   *
   ****************************************************************/


  
/* initializeLevelSetFunctionsOnPatch() initializes the level set 
   * functions on a single patch.  
   * 
   * Arguments:
   *  - patch (in):         Patch on which to initialize level set functions
   *  - current_time (in):  current time
   *  - phi_handle (in):    PatchData handle for phi
   *  - psi_handle (in):    PatchData handle for psi
   * 
   * Return value:          none
   * 
   * NOTES: 
   *  - For codimension-one problems, psi_handle is NOT guaranteed
   *    to be set to a valid PatchData handle, so it should be 
   *    ignored in the user-defined version of this method.
   *
   *  - This is a pure abstract method that the user MUST override in 
   *    order to use the LevelSetFunctionIntegrator class.
   * 
   */
  virtual void initializeLevelSetFunctionsOnPatch(hier::Patch& patch,
                                                  const LSMLIB_REAL time,
                                                  const int phi_handle,
                                                  const int psi_handle) = 0;

  /*!
   * setLevelSetFunctionBoundaryConditions() can be used to set the 
   * data in ghost cells in user-defined ways in order to impose 
   * boundary conditions for the level set functions.  The standard 
   * boundary conditions provided by BoundaryConditionModule class 
   * (which can be managed by the LevelSetMethodAlgorithm input 
   * file) are sufficient for most level set method calculations.
   * 
   * Arguments:
   *  - patch (in):               Patch on which to set boundary conditions
   *                              for the level set functions
   *  - fill_time (in):           time at which boundary data are being filled
   *  - phi_handle (in):          PatchData handle for phi
   *  - psi_handle (in):          PatchData handle for psi
   *  - ghost_width_to_fill (in): width of ghost cells to fill
   * 
   * Return value:                none
   * 
   * NOTES: 
   *  - For codimension-one problems, psi_handle is NOT guaranteed
   *    to be set to a valid PatchData handle, so it should be 
   *    ignored in the user-defined version of this method.
   *
   *  - Because this method is called BEFORE the standard boundary 
   *    conditions are imposed, the user MUST set the boundary condtion 
   *    type to NONE for boundaries where custom boundary conditions are 
   *    to be imposed.  Otherwise, any ghost cells set by this method
   *    will be overwritten. 
   *
   *  - It is virtual with an empty implementation here (rather than
   *    pure virtual) so that users are not required to provide an
   *    implementation.      
   * 
   */
  virtual void setLevelSetFunctionBoundaryConditions(
    hier::Patch& patch,
    const LSMLIB_REAL fill_time,
    const int phi_handle,
    const int psi_handle,
    const hier::IntVector& ghost_width_to_fill){}

  //! @}


  //! @{
  /*!
   ****************************************************************
   * 
   * @name Methods related to time integration 
   *
   ****************************************************************/

  /*!
   * computeStableDtOnPatch() returns a user-specified stable dt on a patch.  
   *
   * Arguments:
   *  - patch (in):                    Patch on which to compute a stable
   *                                   time step size
   *  - lsm_integrator (in):           boost pointer to LevelSetFunctionIntegrator
   *  - velocity_field_strategy (in):  bosst pointer to 
   *                                   LevelSetMethodVelocityFieldStrategy
   * 
   * Return value:                     user-specified time step dt
   * 
   * NOTES: 
   *  - The dt returned by this method need not be related to the
   *    current values of the level set functions or velocity field.
   *    This method is intended to allow the user to specify a completely 
   *    arbitrary time step size.
   *
   *  - If this method returns any value smaller than LSMLIB_REAL_MAX, 
   *    the LevelSetFunctionIntegrator will ignore the internally computed
   *    maximum stable dt and the stable dt returned by the user 
   *    implemented subclass of LevelSetMethodVelocityFieldStrategy.
   *
   *  - The boost pointer to the LevelSetMethodVelocityFieldStrategy object
   *    is provided in case the stable dt calculation requires 
   *    information about the velocity field.
   *
   *  - It is virtual with an default implementation here (rather than
   *    pure virtual) so that users are not required to provide an
   *    implementation. 
   * 
   */
  inline virtual LSMLIB_REAL computeStableDtOnPatch(
    hier::Patch& patch,
    LevelSetFunctionIntegrator* lsm_integrator,
    LevelSetMethodVelocityFieldStrategy* velocity_field_strategy) {
      return LSMLIB_REAL_MAX;
  } 

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
  virtual ~LevelSetMethodPatchStrategy(){}

  //! @}

};

} // end LSMLIB namespace 

#endif
