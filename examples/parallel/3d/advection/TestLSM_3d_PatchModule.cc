/*
 * File:        TestLSM_3d_PatchModule.cc
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date: 2006/03/23 14:00:05 $
 * Description: Implementation for concrete subclass of 
 *              LevelSetMethodPatchStrategy that computes the single patch 
 *              numerical routines for the 3d level set method test problem
 */


#include "TestLSM_3d_PatchModule.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"

// headers for level set method numerical kernels
extern "C" {
  #include "testlsm_3d_patchmodule_fort.h"
}

// SAMRAI namespaces
using namespace geom; 
using namespace pdat; 
using namespace LSMLIB;

// CONSTANTS
const double TestLSM_3d_PatchModule::s_default_radius = 0.25;

TestLSM_3d_PatchModule::TestLSM_3d_PatchModule(
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

void TestLSM_3d_PatchModule::initializeLevelSetFunctionsOnPatch(
  Patch<3>& patch,
  const double data_time,
  const int phi_handle,
  const int psi_handle)
{
  Pointer< CellData<3,double> > level_set_data =
    patch.getPatchData( phi_handle );

  double* level_set_data_ptr = level_set_data->getPointer();

  Pointer< CartesianPatchGeometry<3> > patch_geom 
    = patch.getPatchGeometry();
  const double* dx = patch_geom->getDx();
  const double* x_lower = patch_geom->getXLower();

  Box<3> box = level_set_data->getBox();
  Box<3> ghostbox = level_set_data->getGhostBox();
  const IntVector<3> box_lower = box.lower();
  const IntVector<3> box_upper = box.upper();
  const IntVector<3> ghostbox_lower = ghostbox.lower();
  const IntVector<3> ghostbox_upper = ghostbox.upper();

  switch (d_initial_level_set) 
  {
    case SPHERE: // circle
    {
      INIT_SPHERE(
        level_set_data_ptr,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        x_lower,
        dx,
        d_center,
        &d_radius);
      break;
    }
    default: {}
  }; 

}

void TestLSM_3d_PatchModule::setLevelSetFunctionBoundaryConditions(
    Patch<3>& patch,
    const double fill_time,
    const int phi_handle,
    const int psi_handle,
    const IntVector<3>& ghost_width_to_fill)
{
}

void TestLSM_3d_PatchModule::printClassData(ostream &os) const
{
  os << "\nTestLSM_3d_PatchModule::printClassData..." << endl;
  os << "TestLSM_3d_PatchModule: this = " << (TestLSM_3d_PatchModule*)this 
     << endl;
  os << "d_object_name = " << d_object_name << endl;
  os << "d_initial_level_set = " << d_initial_level_set << endl;

   // KTC - PUT MORE HERE
  os << endl;
}


void TestLSM_3d_PatchModule::getFromInput(
  Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!db.isNull());
#endif

  // set initial level_set_selector
  d_initial_level_set = db->getIntegerWithDefault("initial_level_set", 0);

  // get auxilliary parameters for initial level set
  switch (d_initial_level_set) {
    case SPHERE: {
      d_radius = db->getDoubleWithDefault("radius", s_default_radius);
      if (db->keyExists("center")) {
        db->getDoubleArray("center", d_center, 3);
      } else {
        d_center[0] = 0.0;
        d_center[1] = 0.0;
        d_center[2] = 0.0;
      }
      break;
    }

    default: {      
      TBOX_ERROR(  "TestLSM_3d_PatchModule"
                << "::getFromInput()"
                << ":Invalid type of initial level set"
                << endl );
    }  

  };

}

