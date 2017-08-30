/*
 * File:        PatchModule.h
 * Description: Header for concrete subclass of LevelSetMethodPatchStrategy 
 *              that computes the single patch numerical routines for the 
 *              level set method example problem
 */

#ifndef included_PatchModule
#define included_PatchModule

/*************************************************************************
 *
 * The PatchModule class provides routines for initializing the level 
 * set function.  It does not implement any boundary conditions because 
 * it is assumed that periodic boundary conditions are used for the level 
 * set function.
 *
 *************************************************************************/

// Standard headers
#include <iosfwd>

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/tbox/Dimension.h"

// LevelSetMethod configuration header must be included
// before any other LevelSetMethod header files
#include "LSMLIB_config.h"
#include "LevelSetMethodPatchStrategy.h"

// Class/type declarations
namespace SAMRAI { namespace tbox { class Database; } }
namespace SAMRAI { namespace hier { class IntVector; } }
namespace SAMRAI { namespace hier { class Patch; } }

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace LSMLIB;

class PatchModule: public LevelSetMethodPatchStrategy
{
public:

  /*!
   * Enumeration of Initial Level Set Configurations
   */
  typedef enum {
      // 2D
      CIRCLE = 0,
      STAR = 1,

      // 3D
      SPHERE = 0,
      BUMPY_SPHERE = 1
  } INITIAL_SET_TYPE;

  /*!
   * This constructor sets the object name and reads in the 
   * user input from the specified input database.
   *
   * Arguments:
   *  - input_db (in):         pointer to database containing user input
   *  - dim (in):              problem dimension
   *  - object_name (in):      string name for object
   *
   * Return value:             none
   *
   */
  PatchModule(
    boost::shared_ptr<tbox::Database> input_db,
    const tbox::Dimension& dim,
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
   *  - input_db (in):         pointer to database containing user input
   *  - patch_hierarchy (in):  PatchHierarchy on which to compute velocity field
   *  - grid_geometry (in):    geometry of the computational grid
   *  - object_name (in):      string name for object
   *
   * Return value:             none
   *
   */
  virtual void initializeLevelSetFunctionsOnPatch(hier::Patch& patch,
                                                  const LSMLIB_REAL data_time,
                                                  const int phi_handle,
                                                  const int psi_handle);

  /*!
   * setLevelSetFunctionBoundaryConditions() sets the data in ghost cells 
   * corresponding to physical boundary conditions.  
   */
  virtual void setLevelSetFunctionBoundaryConditions(
    hier::Patch& patch,
    const LSMLIB_REAL fill_time,
    const int phi_handle,
    const int psi_handle,
    const hier::IntVector& ghost_width_to_fill);

  /*!
   * Print all data members for FluidSolver class.
   */
  void printClassData(ostream& os) const;

protected:

  /****************************************************************
   *
   * Utility Methods
   *
   ****************************************************************/

  void getFromInput(boost::shared_ptr<tbox::Database> db);

  /*
   * The object name is used for error/warning reporting and also as a
   * string label for restart database entries.
   */
  string d_object_name;

  /*
   * problem dimension
   */
  tbox::Dimension d_dim;

  /*
   * d_initial_level_set is set by the initial_level_set field in the
   * input database.
   */
  int d_initial_level_set;

  /*
   * other initial level set parameters
   */
  LSMLIB_REAL d_radius;
  LSMLIB_REAL *d_center;

  /*
   * class constants
   */
  static const LSMLIB_REAL s_default_radius;
};

#endif
