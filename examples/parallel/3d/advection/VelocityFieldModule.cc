/*
 * File:        VelocityFieldModule.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Implementation of class that computes the velocity field
 *              for the 3d level set method example program
 */

#include "VelocityFieldModule.h" 

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "IntVector.h"
#include "Patch.h"
#include "VariableContext.h"
#include "VariableDatabase.h"

#include <float.h>

extern "C" {
  #include "velocityfield_fort.h"
}


// SAMRAI namespaces
using namespace pdat;
using namespace LSMLIB;


/* Constructor */
VelocityFieldModule::VelocityFieldModule(
  Pointer<Database> input_db,
  Pointer< PatchHierarchy<3> > patch_hierarchy,
  Pointer< CartesianGridGeometry<3> > grid_geom,
  const string& object_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!input_db.isNull());
  assert(!patch_hierarchy.isNull());
  assert(!grid_geom.isNull());
  assert(!object_name.empty());
#endif

  // set object name, patch hierarchy, and  grid geometry 
  d_object_name = object_name;
  d_patch_hierarchy = patch_hierarchy;
  d_grid_geometry = grid_geom;

  // read in input data
  getFromInput(input_db);

  // Allocate velocity variable
  Pointer< CellVariable<3,LSMLIB_REAL> > velocity = 
    new CellVariable<3,LSMLIB_REAL>("velocity field",3); 
 
  // Register velocity variable with VariableDatabase.
  VariableDatabase<3> *vdb = VariableDatabase<3>::getDatabase();
  Pointer<VariableContext> cur_ctxt = vdb->getContext("CURRENT");
  d_velocity_handle = vdb->registerVariableAndContext(
    velocity, cur_ctxt, IntVector<3>(0));
  vdb->registerPatchDataForRestart(d_velocity_handle);

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
  (void) psi_handle; // psi is meaningless for codimension-one problems
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

    Pointer< PatchLevel<3> > level = d_patch_hierarchy->getPatchLevel(ln);
    computeVelocityFieldOnLevel(level,time);

  } // end loop over hierarchy
}


/* initializeLevelData() */
void VelocityFieldModule::initializeLevelData (
  const Pointer< PatchHierarchy<3> > hierarchy ,
  const int level_number ,
  const LSMLIB_REAL init_data_time ,
  const int phi_handle,
  const int psi_handle,
  const bool can_be_refined ,
  const bool initial_time ,
  const Pointer< PatchLevel<3> > old_level,
  const bool allocate_data)
{

  Pointer< PatchLevel<3> > level = hierarchy->getPatchLevel(level_number);
  if (allocate_data) {
    level->allocatePatchData(d_velocity_handle);
  }

  /*
   * Initialize data on all patches in the level.
   */ 
  computeVelocityFieldOnLevel(level,init_data_time);

}

/* computeVelocityFieldOnLevel() */
void VelocityFieldModule::computeVelocityFieldOnLevel(
  const Pointer< PatchLevel<3> > level,
  const LSMLIB_REAL time) 
{
  for (PatchLevelIterator<3> pi(level); pi; pi++) { // loop over patches
    const int pn = *pi;
    Pointer< Patch<3> > patch = level->getPatch(pn);
    if ( patch.isNull() ) {
      TBOX_ERROR(d_object_name << ": Cannot find patch. Null patch pointer.");
    }

    Pointer< CellData<3,LSMLIB_REAL> > velocity_data = 
      patch->getPatchData( d_velocity_handle );

    Pointer< CartesianPatchGeometry<3> > patch_geom 
      = patch->getPatchGeometry();
#ifdef LSMLIB_DOUBLE_PRECISION
  const double* dx = patch_geom->getDx();
  const double* x_lower = patch_geom->getXLower();
#else
  const double* dx_double = patch_geom->getDx();
  const double* x_lower_double = patch_geom->getXLower();
  float dx[3], x_lower[3];
  dx[0] = dx_double[0]; 
  dx[1] = dx_double[1];
  dx[2] = dx_double[2];
  x_lower[0] = x_lower_double[0]; 
  x_lower[1] = x_lower_double[1];
  x_lower[2] = x_lower_double[2];
#endif

    Box<3> vel_ghostbox = velocity_data->getGhostBox();
    const IntVector<3> vel_ghostbox_lower = vel_ghostbox.lower();
    const IntVector<3> vel_ghostbox_upper = vel_ghostbox.upper();

    Box<3> vel_box = velocity_data->getBox();
    const IntVector<3> vel_lower = vel_box.lower();
    const IntVector<3> vel_upper = vel_box.upper();

    // get velocity data pointers
    LSMLIB_REAL* vel_x_data_ptr = velocity_data->getPointer(0);
    LSMLIB_REAL* vel_y_data_ptr = velocity_data->getPointer(1);
    LSMLIB_REAL* vel_z_data_ptr = velocity_data->getPointer(2);

    switch (d_velocity_field_selector) {
      case 0: { // uniform velocity field (1,0)
        UNIFORM_VELOCITY_X(
          vel_x_data_ptr,
          vel_y_data_ptr,
          vel_z_data_ptr,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_ghostbox_lower[2],
          &vel_ghostbox_upper[2],
          &vel_lower[0],
          &vel_upper[0],
          &vel_lower[1],
          &vel_upper[1],
          &vel_lower[2],
          &vel_upper[2]);
        break;
      }
      case 1: { // uniform velocity field (0,1)
        UNIFORM_VELOCITY_Y(
          vel_x_data_ptr,
          vel_y_data_ptr,
          vel_z_data_ptr,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_ghostbox_lower[2],
          &vel_ghostbox_upper[2],
          &vel_lower[0],
          &vel_upper[0],
          &vel_lower[1],
          &vel_upper[1],
          &vel_lower[2],
          &vel_upper[2]);
        break;
      }
      case 2: { // uniform velocity field (1,1)
        UNIFORM_VELOCITY_XY(
          vel_x_data_ptr,
          vel_y_data_ptr,
          vel_z_data_ptr,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_ghostbox_lower[2],
          &vel_ghostbox_upper[2],
          &vel_lower[0],
          &vel_upper[0],
          &vel_lower[1],
          &vel_upper[1],
          &vel_lower[2],
          &vel_upper[2]);
        break;
      }
      case 3: { // rotating velocity field
        ROTATING_VELOCITY(
          vel_x_data_ptr,
          vel_y_data_ptr,
          vel_z_data_ptr,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_ghostbox_lower[2],
          &vel_ghostbox_upper[2],
          &vel_lower[0],
          &vel_upper[0],
          &vel_lower[1],
          &vel_upper[1],
          &vel_lower[2],
          &vel_upper[2],
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
          vel_z_data_ptr,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_ghostbox_lower[2],
          &vel_ghostbox_upper[2],
          &vel_lower[0],
          &vel_upper[0],
          &vel_lower[1],
          &vel_upper[1],
          &vel_lower[2],
          &vel_upper[2],
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
  Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!db.isNull());
#endif

  // set d_min_dt (declared in parent class) 
#ifdef LSMLIB_DOUBLE_PRECISION
  d_min_dt = db->getDoubleWithDefault("min_dt", DBL_MAX);
#else
  d_min_dt = db->getFloatWithDefault("min_dt", FLT_MAX);
#endif

  // set velocity field selector
  d_velocity_field_selector = 
    db->getIntegerWithDefault("velocity_field", 0);

}

