/*
 * File:        TestLSM_2d_VelocityFieldModule.cc
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date: 2006/05/18 01:09:35 $
 * Description: Implementation of class that computes the velocity field
 *              for the 2d level set method test program
 */

#include "TestLSM_2d_VelocityFieldModule.h" 

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
  #include "testlsm_2d_velocityfield_fort.h"
}

// SAMRAI namespaces
using namespace pdat;
using namespace LSMLIB;


/* Constructor */
TestLSM_2d_VelocityFieldModule::TestLSM_2d_VelocityFieldModule(
  Pointer<Database> input_db,
  Pointer< PatchHierarchy<2> > patch_hierarchy,
  Pointer< CartesianGridGeometry<2> > grid_geom,
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
  Pointer< CellVariable<2,double> > velocity = 
    new CellVariable<2,double>("velocity field",2); 
 
  // Register velocity variable with VariableDatabase.
  VariableDatabase<2> *vdb = VariableDatabase<2>::getDatabase();
  Pointer<VariableContext> cur_ctxt = vdb->getContext("CURRENT");
  d_velocity_handle = vdb->registerVariableAndContext(
    velocity, cur_ctxt, IntVector<2>(0));
  vdb->registerPatchDataForRestart(d_velocity_handle);

  // set d_velocity_never_computed to true to ensure that velocity is
  // computed on first call to computeVelocityField()
  d_velocity_never_computed = true;
}


/* computeVelocityField() */
void TestLSM_2d_VelocityFieldModule::computeVelocityField(
  const double time,
  const int phi_handle,
  const int psi_handle,
  const int component)
{
  (void) psi_handle; // psi is meaningless for codimension-one problems
  (void) component;  // component is not used because this test problem
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

    Pointer< PatchLevel<2> > level = d_patch_hierarchy->getPatchLevel(ln);
    computeVelocityFieldOnLevel(level,time);

  } // end loop over hierarchy
}


/* initializeLevelData() */
void TestLSM_2d_VelocityFieldModule::initializeLevelData (
  const Pointer< PatchHierarchy<2> > hierarchy ,
  const int level_number ,
  const double init_data_time ,
  const int phi_handle,
  const int psi_handle,
  const bool can_be_refined ,
  const bool initial_time ,
  const Pointer< PatchLevel<2> > old_level,
  const bool allocate_data)
{

  Pointer< PatchLevel<2> > level = hierarchy->getPatchLevel(level_number);
  if (allocate_data) {
    level->allocatePatchData(d_velocity_handle);
  }

  /*
   * Initialize data on all patches in the level.
   */ 
  computeVelocityFieldOnLevel(level,init_data_time);

}

/* computeVelocityFieldOnLevel() */
void TestLSM_2d_VelocityFieldModule::computeVelocityFieldOnLevel(
  const Pointer< PatchLevel<2> > level,
  const double time) 
{
  for (PatchLevelIterator<2> pi(level); pi; pi++) { // loop over patches
    const int pn = *pi;
    Pointer< Patch<2> > patch = level->getPatch(pn);
    if ( patch.isNull() ) {
      TBOX_ERROR(d_object_name << ": Cannot find patch. Null patch pointer.");
    }

    Pointer< CellData<2,double> > velocity_data = 
      patch->getPatchData( d_velocity_handle );

    Pointer< CartesianPatchGeometry<2> > patch_geom 
      = patch->getPatchGeometry();
    const double* dx = patch_geom->getDx();
    const double* x_lower = patch_geom->getXLower();

    Box<2> vel_ghostbox = velocity_data->getGhostBox();
    const IntVector<2> vel_ghostbox_lower = vel_ghostbox.lower();
    const IntVector<2> vel_ghostbox_upper = vel_ghostbox.upper();

    Box<2> vel_box = velocity_data->getBox();
    const IntVector<2> vel_lower = vel_box.lower();
    const IntVector<2> vel_upper = vel_box.upper();

    // get velocity data pointers
    double* vel_x_data_ptr = velocity_data->getPointer(0);
    double* vel_y_data_ptr = velocity_data->getPointer(1);

    switch (d_velocity_field_selector) {
      case 0: { // uniform velocity field (1,0)
        uniformvelx_(
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
        uniformvely_(
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
        uniformvelxy_(
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
        rotatingvel_(
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
        double speed = 0.1;
        double omega = 1.0;
        expandingvel_(
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

void TestLSM_2d_VelocityFieldModule::printClassData(ostream& os) const
{
  os << "\nTestLSM_2d_VelocityFieldModule::printClassData..." << endl;
  os << "TestLSM_2d_VelocityFieldModule: this = " << 
     (TestLSM_2d_VelocityFieldModule*)this << endl;
  os << "d_object_name = " << d_object_name << endl;
  os << "d_velocity_field = " << d_velocity_field_selector << endl;

  // KTC - put more here...
  os << endl;
}

void TestLSM_2d_VelocityFieldModule::getFromInput(
  Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!db.isNull());
#endif

  // set d_min_dt (declared in parent class) 
  d_min_dt = db->getDoubleWithDefault("min_dt", DBL_MAX);

  // set test_problem_id
  d_velocity_field_selector = 
    db->getIntegerWithDefault("velocity_field", 0);

}

