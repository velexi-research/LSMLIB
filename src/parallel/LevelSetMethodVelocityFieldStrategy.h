/*
 * File:        LevelSetMethodVelocityFieldStrategy.h
 * Description: Header for strategy for the velocity field for the level
 *              set method
 */

#ifndef included_LevelSetMethodVelocityFieldStrategy_h
#define included_LevelSetMethodVelocityFieldStrategy_h

/*! \class LSMLIB::LevelSetMethodVelocityFieldStrategy
 *
 * \brief
 * The LevelSetMethodVelocityFieldStrategy class defines the interface for
 * supplying the external velocity field and/or normal velocity field for
 * a level set method calculation.
 *
 * In particular, it defines interfaces for computing the velocity field
 * used in a level set method calculation, computing a stable time step
 * size based on the velocity field computation, and accessing the PatchData
 * handles for the velocity field data.  The interface also defines methods
 * for initializing the data required for the velocity field calculation,
 * resetting the PatchHierarchy configuation (after a regridding operation),
 * and tagging cells for refinement.  These methods emulate similarly named
 * methods in the SAMRAI::mesh::StandardTagAndInitStrategy class (in
 * the SAMRAI library).
 *
 *
 * <h3> NOTES: </h3>
 *  - In order to use the LevelSetFunctionIntegrator class, a user MUST
 *    implement a concrete subclass of this ``strategy'' class.
 *
 *  - It is NOT required that the user-defined subclass provides
 *    support for both external (vector) and normal (scalar) velocity
 *    fields.  However, at least one of these velocity fields MUST be
 *    provided by the user-defined subclass.  Depending on the
 *    application either one (or both) velocity fields may be required
 *    in the user-defined subclass.
 *
 *  - If restart capabilities are desired for the application, the
 *    concrete subclass of LevelSetMethodVelocityFieldStrategy
 *    MUST register all PatchData for required restart using the
 *    SAMRAI::VariableDatabase::registerPatchDataForRestart() method.
 *
 *  - Default (empty) implementations of resetHierarchyConfiguration()
 *    and tagCellsForRefinement() are provided for convenience if
      AMR is not required by the user's application.
 *
 */

// Boost headers
#include "boost/smart_ptr/shared_ptr.hpp"

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep

// LSMLIB headers
#include "LSMLIB_config.h"
#include "LevelSetMethodToolbox.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchHierarchy; } }
namespace SAMRAI { namespace hier { class PatchLevel; } }

// Namespaces
using namespace SAMRAI;


/******************************************************************
 *
 * LevelSetMethodVelocityFieldStrategy Class Definition
 *
 ******************************************************************/

namespace LSMLIB {

class LevelSetMethodVelocityFieldStrategy
{
public:

  //! @{
  /*!
   *******************************************************************
   *
   * @name Methods related to calculation of velocity field
   *
   *******************************************************************/

  /*!
   * providesExternalVelocityField() indicates whether the concrete
   * subclass of LevelSetMethodVelocityFieldStrategy provides an
   * external (vector) velocity field for the level set method
   * calculation.
   *
   * Arguments:     none
   *
   * Return value:  true if the object provides calculation of an
   *                external (vector) velocity field for a level set
   *                method calculation; false otherwise
   *
   * NOTES:
   *  - This is a pure abstract method that the user MUST override in
   *    order to use the LevelSetFunctionIntegrator class.
   *
   */
  virtual bool providesExternalVelocityField() const = 0;

  /*!
   * providesNormalVelocityField() indicates whether the concrete
   * subclass of LevelSetMethodVelocityFieldStrategy provides an
   * normal (scalar) velocity field for the level set method
   * calculation.
   *
   * Arguments:     none
   *
   * Return value:  true if the object provides calculation of an normal
   *                (scalar) velocity field for a level set method
   *                calculation; false otherwise
   *
   * NOTES:
   *  - This is a pure abstract method that the user MUST override in
   *    order to use the LevelSetFunctionIntegrator class.
   *
   */
  virtual bool providesNormalVelocityField() const = 0;

  /*!
   * Accessor method for the external (vector) velocity field PatchData
   * handle.
   *
   * Arguments:
   *  - component (in):  component of vector level set function that the
   *                     velocity field handle is being requested for
   *
   * Return value:       PatchData handle for the external velocity field
   *                     data
   *
   * NOTES:
   *  - This is a pure abstract method that the user MUST override in
   *    order to use the LevelSetFunctionIntegrator class.
   *
   */
  virtual int getExternalVelocityFieldPatchDataHandle(
    const int component) const = 0;

  /*!
   * Accessor method for the normal (scalar) velocity field PatchData
   * handle.
   *
   * Arguments:
   *  - level_set_fcn (in):  level set function for which to get
   *                         normal velocity field PatchData handle
   *  - component (in):      component of vector level set function that the
   *                         normal velocity field handle is being requested
   *                         for
   *
   * Return value:           PatchData handle for the normal velocity
   *                         field data
   *
   * NOTES:
   *  - The level set function is required as an argument because the
   *    the normal velocity depends on the orientation of the zero level
   *    set (which is a function of the level set function).
   *
   *  - This is a pure abstract method that the user MUST override in
   *    order to use the LevelSetFunctionIntegrator class.
   *
   */
  virtual int getNormalVelocityFieldPatchDataHandle(
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int component) const = 0;

  /*!
   * computeVelocityField() computes all necessary level set method
   * velocity fields on the entire hierarchy.
   *
   * Arguments:
   *  - time (in):           time that velocity field is to be computed
   *  - phi_handle (in):     PatchData handle for phi
   *  - psi_handle (in):     PatchData handle for psi
   *  - component (in):      component of level set functions for which to
   *                         compute velocity field
   *
   * Return value:           none
   *
   * NOTES:
   *  - Whether the external (vector) velocity field or just a normal
   *    (scalar) velocity is computed depends on the specific application.
   *
   *  - phi_handle, psi_handle, phi_component, and psi_component are provided
   *    for applications that use information about the level set functions
   *    to compute the velocity field.  They are NOT guaranteed to be the
   *    same between calls to computeVelocityField().
   *
   *  - The number of ghostcells for the PatchData associated with
   *    phi_handle and psi_handle is equal to the number required to
   *    compute spatial derivatives using the type and order specified
   *    in the input database for the LevelSetFunctionIntegrator.
   *    Specifically, the number of ghostcells are one, two, and three
   *    for first-, second-, and third-order ENO spatial derivatives,
   *    respectively.  The number of ghostcells is equal to three for
   *    fifth-order WENO spatial derivatives.
   *
   *  - For codimension-one problems, psi_handle is NOT guaranteed
   *    to be set to a valid PatchData handle, so it should be
   *    ignored in the user-defined version of this method.
   *
   *  - This is a pure abstract method that the user MUST override in
   *    order to use the LevelSetFunctionIntegrator class.
   *
   */
  virtual void computeVelocityField(
    const LSMLIB_REAL time,
    const int phi_handle,
    const int psi_handle,
    const int component) = 0;

  //! @}


  //! @{
  /*!
   *******************************************************************
   *
   * @name Methods related to time integration
   *
   *******************************************************************/

  /*!
   * setCurrentTime() sets the current time for the
   * VelocityFieldStrategy class so that the simulation
   * time for the velocity field calculation can be synchronized with
   * the simulation time for the level set method calculation.
   *
   * Arguments:
   *  - time (in):   new current time
   *
   * Return value:   none
   *
   * NOTES:
   *  - This is a pure abstract method that the user MUST override in
   *    order to use the LevelSetFunctionIntegrator class.
   *
   */
  virtual void setCurrentTime(const LSMLIB_REAL time) = 0;

  /*!
   * computeStableDt() returns the maximum acceptable (stable) time
   * step size allowed in order to maintain numerical stability of
   * the calculation used to compute the velocity field.  This often
   * comes from the physics of the problem (e.g. inherent time-variation
   * in an external velocity field).
   *
   * Arguments:     none
   *
   * Return value:  maximum acceptable time step
   *
   * NOTES:
   *  - This is a pure abstract method that the user MUST override in
   *    order to use the LevelSetFunctionIntegrator class.
   *
   */
  virtual LSMLIB_REAL computeStableDt() = 0;

  //! @}


  //! @{
  /*!
   ************************************************************************
   *
   * @name Methods for initializing data
   *
   * NOTES:
   *  - initializeLevelData() essentially emulate the method of the same
   *    name declared in the SAMRAI::mesh::StandardTagAndInitStrategy class
   *    but with access to level set function data
   *
   ************************************************************************/

  /*!
   * initializeLevelData() allocates and initializes any PatchData
   * related to the computation of the velocity field for a PatchLevel
   * that has been newly added to the PatchHierarchy.
   *
   * Arguments:
   *  - hierarchy (in):       PatchHierarchy on which to tag cells for
   *                          refinement
   *  - level_number (in):    PatchLevel number on which to tag cells for
   *                          refinement
   *  - phi_handle (in):      PatchData handle for phi
   *  - psi_handle (in):      PatchData handle for psi
   *  - can_be_refined (in):  true if this is NOT the finest level in the
   *                          PatchHierarchy; false it if is
   *  - init_data_time (in):  true if the PatchLevel is being introduced for
   *                          the first time; false otherwise
   *  - old_level (in):       old PatchLevel from which data for new
   *                          PatchLevel should be taken
   *  - allocate_data (in):   true if PatchData needs to be allocated before
   *                          it is initialized; false otherwise
   *
   * Return value:            none
   *
   * NOTES:
   *  - The number of ghostcells for phi and psi data depends on the
   *    type and order of spatial derivative selected for the
   *    LevelSetFunctionIntegrator:
   *
   *    - ENO-1:  1
   *    - ENO-2:  2
   *    - ENO-3:  3
   *    - WENO-5: 3
   *
   *  - For codimension-one problems, psi_handle is NOT guaranteed
   *    to be set to a valid PatchData handle, so it should be
   *    ignored in the user-defined version of this method.
   *
   *  - This is a pure abstract method that the user MUST override in
   *    order to use the LevelSetFunctionIntegrator class.
   *
   */
  virtual void initializeLevelData(
    const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int level_number,
    const LSMLIB_REAL init_data_time,
    const int phi_handle,
    const int psi_handle,
    const bool can_be_refined,
    const bool initial_time,
    const boost::shared_ptr<hier::PatchLevel> old_level
    = boost::shared_ptr<hier::PatchLevel>(),
    const bool allocate_data = true) = 0;


  //! @}


  //! @{
  /*!
   ************************************************************************
   *
   * @name Methods for managing grid configuration
   *
   * NOTES:
   *  - resetHierarchyConfiguration() exactly mimics the method of the
   *    same name declared in the SAMRAI::mesh::StandardTagAndInitStrategy
   *    class
   *
   *  - tagCellsForRefinement() essentially emulates applyGradientDetector()
   *    in the the SAMRAI::mesh::StandardTagAndInitStrategy class
   ************************************************************************/

  /*!
   * resetHierarchyConfiguration() resets any internal
   * hierarchy-dependent information.
   *
   * Arguments:
   *  - new_hierarchy (in):   PatchHierarchy to configure
   *  - coarsest_level (in):  level number of coarsest PatchLevel to reset
   *  - finest_level (in):    level number of finest Patchlevel to reset
   *
   * Return value:            none
   *
   * NOTES:
   *   - This method is virtual with an empty implementation here
   *     (rather than pure virtual) so that users do not need to
   *     provide an implementation when the method is not needed
   *     (e.g. when parallel computation and AMR is not required).
   *
   */
  virtual void resetHierarchyConfiguration(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    int coarsest_level,
    int finest_level){}

  /*!
   * tagCellsForRefinement() sets integer tags to "1" in cells where
   * refinement of the given level should occur.
   *
   * Arguments:
   *  - hierarchy (in):       PatchHierarchy on which to tag cells for
   *                          refinement
   *  - level_number (in):    PatchLevel number on which to tag cells for
   *                          refinement
   *  - tag_handle (in):      PatchData handle for the cell-centered integer
   *                          tag data
   *
   * Return value:            none
   *
   * NOTES:
   *   - This method is virtual with an empty implementation here
   *     (rather than pure virtual) so that users do not need to
   *     provide an implementation when the method is not needed
   *     (e.g. when AMR is not required).
   *
   */
  virtual void tagCellsForRefinement(
      const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
      const int level_number,
      const int tag_handle){}

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
  virtual ~LevelSetMethodVelocityFieldStrategy(){}

  //! @}

};

} // end LSMLIB namespace

#endif
