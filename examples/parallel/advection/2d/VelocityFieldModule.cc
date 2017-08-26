/*
 * File:        VelocityFieldModule.cc
 * Description: Implementation of class that computes the velocity field
 *              for the level set method
 */

// Class header
#include "VelocityFieldModule.h"

// Standard headers
#include <assert.h>
#include <cfloat>
#include <sstream>

// Boost headers
// IWYU pragma: no_include <boost/smart_ptr/detail/operator_bool.hpp>
#include <boost/smart_ptr/make_shared_object.hpp>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Utilities.h"

extern "C" {
  #include "velocityfield_fort.h"
}

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchData; } }
namespace SAMRAI { namespace hier { class PatchGeometry; } }
namespace SAMRAI { namespace hier { class VariableContext; } }

// SAMRAI namespaces
using namespace pdat;


/* Constructor */
VelocityFieldModule::VelocityFieldModule(
  boost::shared_ptr<tbox::Database> input_db,
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geom,
  const string& object_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(input_db);
  assert(patch_hierarchy);
  assert(grid_geom);
  assert(!object_name.empty());
#endif

  // set object name, patch hierarchy, and  grid geometry
  d_object_name = object_name;
  d_patch_hierarchy = patch_hierarchy;
  d_grid_geometry = grid_geom;

  // read in input data
  getFromInput(input_db);

  // Allocate velocity variable
  const tbox::Dimension dim = patch_hierarchy->getDim();
  const int depth = dim.getValue();
  boost::shared_ptr< pdat::CellVariable<LSMLIB_REAL> > velocity =
      boost::shared_ptr< pdat::CellVariable<LSMLIB_REAL> >(
          new pdat::CellVariable<LSMLIB_REAL>(dim, "velocity field", depth));

  // Register velocity variable with VariableDatabase.
  hier::VariableDatabase *vdb = hier::VariableDatabase::getDatabase();
  boost::shared_ptr<hier::VariableContext> cur_ctxt = vdb->getContext("CURRENT");
  d_velocity_handle = vdb->registerVariableAndContext(
    velocity, cur_ctxt, hier::IntVector(d_patch_hierarchy->getDim(), 0));

  hier::PatchDataRestartManager::getManager()->
    registerPatchDataForRestart(d_velocity_handle);

  // set d_velocity_never_computed to true to ensure that velocity is
  // computed on first call to computeVelocityField()
  d_velocity_never_computed = true;
}


/* computeVelocityField() */
void VelocityFieldModule::computeVelocityField(
  const LSMLIB_REAL time,
  const int phi_handle,
  const int psi_handle,
  const int component)
{
  (void) psi_handle; // psi is meaningless for co-dimension one problems
  (void) component;  // component is not used because this example problem
                     // only has one component for level set function

  // only carry out computation if the time has changed
  if (!d_velocity_never_computed && (d_current_time == time)) return;

  // set d_velocity_never_computed to false
  d_velocity_never_computed = false;

  // update the current time
  d_current_time = time;

  // set velocity on all levels of hierarchy
  const int finest_level = d_patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    boost::shared_ptr< hier::PatchLevel > level =
        d_patch_hierarchy->getPatchLevel(ln);
    computeVelocityFieldOnLevel(level, time, phi_handle);

  } // end loop over hierarchy
}


/* initializeLevelData() */
void VelocityFieldModule::initializeLevelData (
  const boost::shared_ptr<hier::PatchHierarchy> hierarchy ,
  const int level_number ,
  const LSMLIB_REAL init_data_time ,
  const int phi_handle,
  const int psi_handle,
  const bool can_be_refined ,
  const bool initial_time ,
  const boost::shared_ptr< hier::PatchLevel > old_level,
  const bool allocate_data)
{
  (void) psi_handle;  // psi is meaningless for co-dimension one problems

  boost::shared_ptr< hier::PatchLevel > level =
      hierarchy->getPatchLevel(level_number);
  if (allocate_data) {
    level->allocatePatchData(d_velocity_handle);
  }

  /*
   * Initialize data on all patches in the level.
   */
  computeVelocityFieldOnLevel(level,init_data_time,phi_handle);

}

/* computeVelocityFieldOnLevel() */
void VelocityFieldModule::computeVelocityFieldOnLevel(
  const boost::shared_ptr< hier::PatchLevel > level,
  const LSMLIB_REAL time,
  const int phi_handle)
{
  for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
    // loop over patches
    boost::shared_ptr<hier::Patch> patch = *pi;
    if (!patch) {
      TBOX_ERROR(d_object_name << ": Cannot find patch. Null patch pointer.");
    }

    boost::shared_ptr< CellData<LSMLIB_REAL> > velocity_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
            patch->getPatchData(d_velocity_handle));

    boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch->getPatchGeometry());

#ifdef LSMLIB_DOUBLE_PRECISION
  const double* dx = patch_geom->getDx();
  const double* x_lower = patch_geom->getXLower();
#else
  const double* dx_double = patch_geom->getDx();
  const double* x_lower_double = patch_geom->getXLower();
  float dx[2], x_lower[2];
  dx[0] = dx_double[0]; dx[1] = dx_double[1];
  x_lower[0] = x_lower_double[0]; x_lower[1] = x_lower_double[1];
#endif

    hier::Box vel_ghostbox = velocity_data->getGhostBox();
    const hier::IntVector vel_ghostbox_lower = vel_ghostbox.lower();
    const hier::IntVector vel_ghostbox_upper = vel_ghostbox.upper();

    hier::Box vel_box = velocity_data->getBox();
    const hier::IntVector vel_lower = vel_box.lower();
    const hier::IntVector vel_upper = vel_box.upper();

    // get velocity data pointers
    LSMLIB_REAL* vel_x_data_ptr = velocity_data->getPointer(0);
    LSMLIB_REAL* vel_y_data_ptr = velocity_data->getPointer(1);

    switch (d_velocity_field_selector) {
      case 0: { // uniform velocity field (1,0)
        UNIFORM_VELOCITY_X(
          vel_x_data_ptr,
          vel_y_data_ptr,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_lower[0],
          &vel_upper[0],
          &vel_lower[1],
          &vel_upper[1]);
        break;
      }
      case 1: { // uniform velocity field (0,1)
        UNIFORM_VELOCITY_Y(
          vel_x_data_ptr,
          vel_y_data_ptr,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_lower[0],
          &vel_upper[0],
          &vel_lower[1],
          &vel_upper[1]);
        break;
      }
      case 2: { // uniform velocity field (1,1)
        UNIFORM_VELOCITY_XY(
          vel_x_data_ptr,
          vel_y_data_ptr,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_lower[0],
          &vel_upper[0],
          &vel_lower[1],
          &vel_upper[1]);
        break;
      }
      case 3: { // rotating velocity field
        ROTATING_VELOCITY(
          vel_x_data_ptr,
          vel_y_data_ptr,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_lower[0],
          &vel_upper[0],
          &vel_lower[1],
          &vel_upper[1],
          dx,
          x_lower);
        break;
      }
      case 4: { // oscillating expanding/contracting velocity field
        // (u,v) = speed*cos(omega*time) * (x/r,y/r)
        LSMLIB_REAL speed = 0.1;
        LSMLIB_REAL omega = 1.0;
        EXPANDING_VELOCITY(
          vel_x_data_ptr,
          vel_y_data_ptr,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_lower[0],
          &vel_upper[0],
          &vel_lower[1],
          &vel_upper[1],
          dx,
          x_lower,
          &speed,
          &omega,
          &time);
        break;
      }
      default: {}
    }; // switch statement

  }  // loop over patches

}

void VelocityFieldModule::printClassData(ostream& os) const
{
  os << "\nVelocityFieldModule::printClassData..." << endl;
  os << "VelocityFieldModule: this = " <<
     (VelocityFieldModule*)this << endl;
  os << "d_object_name = " << d_object_name << endl;
  os << "d_velocity_field = " << d_velocity_field_selector << endl;

  // KTC - put more here...
  os << endl;
}

void VelocityFieldModule::getFromInput(
  boost::shared_ptr<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(db);
#endif

  // set d_min_dt (declared in parent class)
#ifdef LSMLIB_DOUBLE_PRECISION
  d_min_dt = db->getDoubleWithDefault("min_dt", DBL_MAX);
#else
  d_min_dt = db->getFloatWithDefault("min_dt", FLT_MAX);
#endif

  // set velocity field type
  d_velocity_field_selector =
    db->getIntegerWithDefault("velocity_field", 0);

}
