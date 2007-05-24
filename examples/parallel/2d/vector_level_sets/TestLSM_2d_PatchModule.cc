/*
 * File:        TestLSM_2d_PatchModule.cc
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.3 $
 * Modified:    $Date: 2006/03/23 14:00:04 $
 * Description: Implementation for concrete subclass of 
 *              LevelSetMethodPatchStrategy that computes the single patch 
 *              numerical routines for the level set method test problem
 */


#include "TestLSM_2d_PatchModule.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"

// headers for level set method numerical kernels
extern "C" {
  #include "testlsm_2d_patchmodule_fort.h"
}

// SAMRAI namespaces
using namespace geom; 
using namespace pdat; 

// CONSTANTS
const double TestLSM_2d_PatchModule::s_default_radius = 0.25;

TestLSM_2d_PatchModule::TestLSM_2d_PatchModule(
  Pointer<Database> input_db,
  const string& object_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!input_db.isNull());
  assert(!object_name.empty());
#endif

  // set object name and grid geometry
  d_object_name = object_name;

  // read in input data
  getFromInput(input_db);

}

void TestLSM_2d_PatchModule::initializeLevelSetFunctionsOnPatch(
  Patch<2>& patch,
  const double data_time,
  const int phi_handle,
  const int psi_handle)
{
  Pointer< CellData<2,double> > level_set_data =
    patch.getPatchData( phi_handle );

  Pointer< CartesianPatchGeometry<2> > patch_geom 
    = patch.getPatchGeometry();
  const double* dx = patch_geom->getDx();
  const double* x_lower = patch_geom->getXLower();

  Box<2> box = level_set_data->getBox();
  Box<2> ghostbox = level_set_data->getGhostBox();
  const IntVector<2> box_lower = box.lower();
  const IntVector<2> box_upper = box.upper();
  const IntVector<2> ghostbox_lower = ghostbox.lower();
  const IntVector<2> ghostbox_upper = ghostbox.upper();

  // initialize component 0 of phi
  double* level_set_data_ptr = level_set_data->getPointer(0);
  switch (d_initial_level_set_0) 
  {
    case CIRCLE: // circle
    {
      INIT_CIRCLE(
        level_set_data_ptr,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        x_lower,
        dx,
        d_center_0,
        &d_radius_0);
      break;
    }
    default: {}
  }; 

  // initialize component 1 of phi
  level_set_data_ptr = level_set_data->getPointer(1);
  switch (d_initial_level_set_1) 
  {
    case CIRCLE: // circle
    {
      INIT_CIRCLE(
        level_set_data_ptr,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        x_lower,
        dx,
        d_center_1,
        &d_radius_1);
      break;
    }
    default: {}
  }; 

}

void TestLSM_2d_PatchModule::setLevelSetFunctionBoundaryConditions(
    Patch<2>& patch,
    const double fill_time,
    const int phi_handle,
    const int psi_handle,
    const IntVector<2>& ghost_width_to_fill)
{
}

void TestLSM_2d_PatchModule::printClassData(ostream &os) const
{
  os << "\nTestLSM_2d_PatchModule::printClassData..." << endl;
  os << "TestLSM_2d_PatchModule: this = " << (TestLSM_2d_PatchModule*)this 
     << endl;
  os << "d_object_name = " << d_object_name << endl;
  os << "d_initial_level_set_0 = " << d_initial_level_set_0 << endl;
  os << "d_initial_level_set_1 = " << d_initial_level_set_1 << endl;

 // KTC - PUT MORE HERE
  os << endl;
}


void TestLSM_2d_PatchModule::getFromInput(
  Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!db.isNull());
#endif

  // set initial level_set_selector
  d_initial_level_set_0 = db->getIntegerWithDefault("initial_level_set_0", 0);
  d_initial_level_set_1 = db->getIntegerWithDefault("initial_level_set_1", 0);

  // get auxilliary parameters for initial level set 0
  switch (d_initial_level_set_0) {
    case CIRCLE: {
      d_radius_0 = db->getDoubleWithDefault("radius_0", s_default_radius);
      if (db->keyExists("center_0")) {
        db->getDoubleArray("center_0", d_center_0, 2);
      } else {
        d_center_0[0] = 0.0;
        d_center_0[1] = 0.0;
      }
      break;
    }

    default: { 
      TBOX_ERROR(  "TestLSM_2d_PatchModule"
                << "::getFromInput()"
                << ":Invalid type of initial level set for component 0"
                << endl );
    } 
  };

  // get auxilliary parameters for initial level set 1
  switch (d_initial_level_set_1) {
    case CIRCLE: {
      d_radius_1 = db->getDoubleWithDefault("radius_1", s_default_radius);
      if (db->keyExists("center_1")) {
        db->getDoubleArray("center_1", d_center_1, 2);
      } else {
        d_center_1[0] = 0.0;
        d_center_1[1] = 0.0;
      }
      break;
    }

    default: { 
      TBOX_ERROR(  "TestLSM_2d_PatchModule"
                << "::getFromInput()"
                << ":Invalid type of initial level set for component 0"
                << endl );
    } 
  };

}

