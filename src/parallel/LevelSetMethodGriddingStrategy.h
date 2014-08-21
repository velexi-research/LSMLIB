/*
 * File:        LevelSetMethodGriddingStrategy.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Header file for the level set method gridding strategy class
 */
 
#ifndef included_LevelSetMethodGriddingStrategy_h
#define included_LevelSetMethodGriddingStrategy_h

/*! \class LSMLIB::LevelSetMethodGriddingStrategy
 *
 * \brief
 * The LevelSetMethodGriddingStrategy class defines the interface for
 * data initialization and management of grid refinement for the level
 * set method computation.
 *
 */

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/StandardTagAndInitStrategy.h"

#include "LSMLIB_config.h"
#include "LevelSetMethodVelocityFieldStrategy.h"

// namespaces
using namespace SAMRAI;
using namespace mesh;


/******************************************************************
 *
 * LevelSetMethodGriddingStrategy Class Definition
 *
 ******************************************************************/

namespace LSMLIB {

template<int DIM> class LevelSetMethodGriddingStrategy
{
public:

  //! @{
  /*!
   ****************************************************************
   * 
   * @name Method for registering velocity field 
   *
   ****************************************************************/

  /*!
   * registerVelocityFieldStrategy() registers the specified instance of 
   * a concrete subclass of the LevelSetMethodVelocityFieldStrategy class
   * with LevelSetMethodGriddingStrategy object.  
   *
   * Arguments:     
   *  - velocity_field_strategy (in):  pointer to instance of subclass of 
   *                                   LevelSetMethodVelocityFieldStrategy
   *                                   used to manage the variables involved 
   *                                   in the calculation of the velocity field
   *
   * Return value:                     none
   *
   */
  virtual void registerVelocityFieldStrategy(
     LevelSetMethodVelocityFieldStrategy<DIM>* velocity_field_strategy) = 0;

  //! @{
  /*!
   ****************************************************************
   * 
   * @name Methods for managing grid configuration
   *
   ****************************************************************/

  /*!
   * initializePatchHierarchy() constructs the PatchHierarchy and 
   * initializes the level set functions and variables involved in 
   * the computation of the velocity field.
   *
   * Arguments:     
   *  - time (in):   simulation time that PatchHierarchy is being 
   *                 initialized
   *
   * Return value:   none
   *
   * NOTES:
   *  - all variables required to calculate the velocity field MUST
   *    be registered using registerTagAndInitStrategy() before the 
   *    PatchHierarchy is constructed and initialized.
   */
  virtual void initializePatchHierarchy(const LSMLIB_REAL time) = 0;

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
   * NOTES:
   *  - Overrides StandardTagAndInitialize::resetHierarchyConfiguration()
   *
   */
  virtual void resetHierarchyConfiguration(
    const Pointer< BasePatchHierarchy<DIM> > hierarchy,
    const int coarsest_level,
    const int finest_level) = 0;

  /*!
   * regridPatchHierarchy() regrids the entire PatchHierarchy and 
   * reinitializes the data on the PatchHierarchy using interpolation
   * and averaging, as necessary.
   *
   * Arguments:     
   *  - time (in):   simulation time when PatchHierarchy is regrid
   *
   * Return value:   none
   *
   */
  virtual void regridPatchHierarchy(LSMLIB_REAL time) = 0;

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
  virtual ~LevelSetMethodGriddingStrategy(){}

  //! @}

};

} // end LSMLIB namespace

#endif
