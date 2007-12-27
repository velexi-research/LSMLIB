/*
 * File:        PatchModule.h
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date: 2006/04/10 15:25:50 $
 * Description: Header for concrete subclass of LevelSetMethodPatchStrategy 
 *              that computes the single patch numerical routines for the 
 *              3d level set method example problem
 */

#ifndef included_PatchModule
#define included_PatchModule

/*************************************************************************
 *
 * This PatchModule class provides routines for initializing
 * the level set function.  It does not implement any boundary conditions
 * because it is assumed that periodic boundary conditions are used for
 * the level set function.
 *
 *************************************************************************/

   
// SAMRAI configuration header must be included
// before any other SAMRAI header files
#include "SAMRAI_config.h"

#include <string>
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Pointer.h"

// LevelSetMethod configuration header must be included
// before any other LevelSetMethod header files
#include "LSMLIB_config.h"
#include "LevelSetMethodPatchStrategy.h"
#include "LevelSetMethodVelocityFieldStrategy.h"

// SAMRAI namespaces
using namespace SAMRAI;
using namespace hier;
using namespace tbox;
using namespace LSMLIB;

class PatchModule:
  public LevelSetMethodPatchStrategy<3>
{
public:

  /*!
   * This constructor sets the object name and reads in the 
   * user input from the specified input database.
   *
   * Arguments:
   *  - object_name (in):      string name for object
   *
   * Return value:             none
   *
   */
  PatchModule(
    const string& object_name = "PatchModule");

  /*!
   * Empty destructor.
   */
  virtual ~PatchModule() {};


  /****************************************************************
   *
   * Methods Inherited from LevelSetMethodPatchStrategy
   *
   ****************************************************************/

  /*!
   * initializeLevelSetFunctionsOnPatch() initializes the level set
   * function on the patch based on the value of 
   * "initial_level_set" in the input database.
   *
   * Arguments:
   *  - 
   *
   * Return value:             none
   *
   */
  virtual void initializeLevelSetFunctionsOnPatch(Patch<3>& patch,
                                                  const LSMLIB_REAL data_time,
                                                  const int phi_handle,
                                                  const int psi_handle);

  /*!
   * setLevelSetFunctionBoundaryConditions() sets the data in ghost cells 
   * corresponding to boundary conditions.  
   */
  virtual void setLevelSetFunctionBoundaryConditions(
    Patch<3>& patch,
    const LSMLIB_REAL fill_time,
    const int phi_handle,
    const int psi_handle,
    const IntVector<3>& ghost_width_to_fill);

  /*!
   * Print all data members for FluidSolver class.
   */
  void printClassData(ostream& os) const;

protected:

  /*
   * The object name is used for error/warning reporting and also as a
   * string label for restart database entries.
   */
  string d_object_name;

};

#endif
