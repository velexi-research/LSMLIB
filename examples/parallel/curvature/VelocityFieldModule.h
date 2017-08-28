/*
 * File:        VelocityFieldModule.h
 * Description: Header for class that computes the velocity field for
 *              the level set method
 */

#ifndef included_VelocityFieldModule
#define included_VelocityFieldModule

/*************************************************************************
 *
 * The VelocityFieldModule class provides several simple
 * 2D external velocity fields to be used by the LevelSetMethodAlgorithm.
 *
 *************************************************************************/
   
// Standard headers
#include <iosfwd>
   
// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/math/PatchCellDataBasicOps.h"

// Level set method velocity field interface definition 
// LevelSetMethod configuration header must be included
// before any other LevelSetMethod header files
#include "LSMLIB_config.h"
#include "LevelSetMethodToolbox.h"
#include "LevelSetMethodVelocityFieldStrategy.h"

// Class/type declarations
namespace SAMRAI { namespace geom { class CartesianGridGeometry; } }
namespace SAMRAI { namespace hier { class PatchHierarchy; } }
namespace SAMRAI { namespace hier { class PatchLevel; } }
namespace SAMRAI { namespace tbox { class Database; } }

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace LSMLIB;

class VelocityFieldModule:
    public LevelSetMethodVelocityFieldStrategy
{
public:

  /*!
   * The constructor for VelocityFieldModule caches pointers
   * to the hierarchy and geometry objects that are to be used in the 
   * level set method computation and sets up the velocity field 
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
  VelocityFieldModule(
    boost::shared_ptr<tbox::Database> input_db,
    boost::shared_ptr< hier::PatchHierarchy > patch_hierarchy,
    boost::shared_ptr< geom::CartesianGridGeometry > grid_geometry,
    const string& object_name = "VelocityFieldModule");

  /*!
   * The destructor for VelocityFieldModule does nothing.
   */
  virtual ~VelocityFieldModule(){}


  /****************************************************************
   *
   * Methods Inherited from LevelSetMethodVelocityFieldStrategy
   *
   ****************************************************************/

  /*!
   * providesExternalVelocityField() always returns false because
   * curvature-driven motion provides a velocity in the direction
   * normal to the surface.
   *
   * Arguments:     none
   *
   * Return value:  returns false
   *
   */
  virtual inline bool providesExternalVelocityField() const {
    return false;
  }

  /*!
   * providesNormalVelocityField() always returns true because
   * curvature-driven motion provides a velocity in the direction
   * normal to the surface.
   *
   * Arguments:     none
   *
   * Return value:  returns true
   *
   */
  virtual inline bool providesNormalVelocityField() const {
    return true;
  }

  /*!
   * getExternalVelocityFieldPatchDataHandle() returns -1 (a bogus
   * PatchData handle value) because this example module does not 
   * provide an external velocity field.
   * 
   * Arguments:     
   *  - component (in):  component of vector level set function that the
   *                     velocity field handle is being requested for
   *
   * Return value:       PatchData handle for the velocity field data
   * 
   */
  virtual inline int getExternalVelocityFieldPatchDataHandle(
    const int component) const
  {
    (void) component;

    return -1;
  }

  /*!
   * getNormalVelocityFieldPatchDataHandle() returns the PatchData handle for
   * the normal velocity field.
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
   */
  virtual inline int getNormalVelocityFieldPatchDataHandle(
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int component) const
  {
    (void) level_set_fcn;
    (void) component;

    return d_velocity_handle;
  }

  /*!
   *
   * setCurrentTime() sets the current time so that the simulation
   * time for the velocity field calculation can be synchronized with
   * the simulation time for the level set method calculation.
   *   
   * Arguments:     
   *  - time (in):   new current time
   *
   * Return value:  none
   *
   */
  virtual inline void setCurrentTime(const LSMLIB_REAL time)
  {
    d_current_time = time;
  }

  /*!
   * computeStableDt() returns the the time step specified by the 
   * user in the input file or the default value of 1.0.
   *
   * Arguments:     none
   *
   * Return value:  maximum acceptable (stable) time step
   *
   */
  virtual inline LSMLIB_REAL computeStableDt()
  {
    return d_min_dt;
  }

  /*!
   * computeVelocityField() sets the velocity field on the entire 
   * hierarchy based on the time and the velocity_field set in the 
   * input database.
   *
   * Arguments:
   *  - time (in):        time that velocity field is to be computed
   *  - phi_handle (in):  PatchData handle for phi
   *  - psi_handle (in):  PatchData handle for psi
   *  - component (in):   component of level set functions for which to
   *                      compute velocity field
   *
   * Return value:        none
   *
   */
  virtual void computeVelocityField(
    const LSMLIB_REAL time,
    const int phi_handle,
    const int psi_handle,
    const int component);


  /*!
   * Allocate and initialize data for a new level in the patch hierarchy.
   */
  virtual void initializeLevelData (
    const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
    const int level_number,
    const LSMLIB_REAL init_data_time,
    const int phi_handle,
    const int psi_handle,
    const bool can_be_refined,
    const bool initial_time,
    const boost::shared_ptr< hier::PatchLevel > old_level
      = boost::shared_ptr< hier::PatchLevel >(),
    const bool allocate_data = true);

  /*!
   * Print all data members for VelocityFieldModule class.
   */
  void printClassData(ostream& os) const;

protected:

  /*
   * computeVelocityFieldOnLevel() computes the velocity field on an
   * entire PatchLevel based on the time and the velocity_field set 
   * in the input database.
   */
  void computeVelocityFieldOnLevel(
    const boost::shared_ptr< hier::PatchLevel > level, 
    const LSMLIB_REAL time,
    const int phi_handle);

  /*
   * These private member functions read data from input. 
   *
   * An assertion results if the database pointer is null.
   */
  void getFromInput(boost::shared_ptr<tbox::Database> db);

  /*
   * The object name is used for error/warning reporting and also as a 
   * string label for restart database entries. 
   */
  string d_object_name;

  /*
   * Pointer to the patch hierarchy.
   */
  boost::shared_ptr< hier::PatchHierarchy > d_patch_hierarchy;

  /*
   * Cache pointer to the grid geometry object to set up initial data, 
   * set physical boundary conditions, and register plot variables.
   */
  boost::shared_ptr< geom::CartesianGridGeometry > d_grid_geometry;

  /*
   * Pointer to the math operations object.
   */
  math::PatchCellDataBasicOps<LSMLIB_REAL> d_math_ops;

  /*
   * current time
   */
  LSMLIB_REAL d_current_time;

  /*
   * flag indicating if velocity has ever been computed
   */
  bool d_velocity_never_computed;

  /*
   * Physical parameters
   */

  // proportionality constant between mean curvature and normal velocity
  double d_b;

  /*
   * minimum time step size (read in from input file or set to default value)
   */
  LSMLIB_REAL d_min_dt;

  /*
   * PatchData handles.
   */
  int d_velocity_handle; // velocity
  int d_grad_phi_handle; // grad(phi)
  int d_hessian_phi_handle; // Hessian of phi
  int d_kappa_handle; // mean curvature
};

#endif
