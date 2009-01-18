/*
 * File:        PatchModule.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.3 $
 * Modified:    $Date$
 * Description: Implementation for concrete subclass of 
 *              LevelSetMethodPatchStrategy that computes the single patch 
 *              numerical routines for the level set method example problem
 */


#include "PatchModule.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"

// headers for level set method numerical kernels
extern "C" {
  #include "patchmodule_fort.h"
}

// SAMRAI namespaces
using namespace geom; 
using namespace pdat; 

// CONSTANTS
const LSMLIB_REAL PatchModule::s_default_radius = 0.25;

PatchModule::PatchModule(
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

void PatchModule::initializeLevelSetFunctionsOnPatch(
  Patch<2>& patch,
  const LSMLIB_REAL data_time,
  const int phi_handle,
  const int psi_handle)
{
  Pointer< CellData<2,LSMLIB_REAL> > level_set_data =
    patch.getPatchData( phi_handle );

  Pointer< CartesianPatchGeometry<2> > patch_geom 
    = patch.getPatchGeometry();

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

  Box<2> box = level_set_data->getBox();
  Box<2> ghostbox = level_set_data->getGhostBox();
  const IntVector<2> box_lower = box.lower();
  const IntVector<2> box_upper = box.upper();
  const IntVector<2> ghostbox_lower = ghostbox.lower();
  const IntVector<2> ghostbox_upper = ghostbox.upper();

  // initialize component 0 of phi
  LSMLIB_REAL* level_set_data_ptr = level_set_data->getPointer(0);
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

void PatchModule::setLevelSetFunctionBoundaryConditions(
    Patch<2>& patch,
    const LSMLIB_REAL fill_time,
    const int phi_handle,
    const int psi_handle,
    const IntVector<2>& ghost_width_to_fill)
{
}

void PatchModule::printClassData(ostream &os) const
{
  os << "\nPatchModule::printClassData..." << endl;
  os << "PatchModule: this = " << (PatchModule*)this 
     << endl;
  os << "d_object_name = " << d_object_name << endl;
  os << "d_initial_level_set_0 = " << d_initial_level_set_0 << endl;
  os << "d_initial_level_set_1 = " << d_initial_level_set_1 << endl;

 // KTC - PUT MORE HERE
  os << endl;
}


void PatchModule::getFromInput(
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

#ifdef LSMLIB_DOUBLE_PRECISION
      d_radius_0 = db->getDoubleWithDefault("radius_0", s_default_radius);
#else
      d_radius_0 = db->getFloatWithDefault("radius_0", s_default_radius);
#endif

      if (db->keyExists("center_0")) {
#ifdef LSMLIB_DOUBLE_PRECISION
        db->getDoubleArray("center_0", d_center_0, 2);
#else
        db->getFloatArray("center_0", d_center_0, 2);
#endif

      } else {
        d_center_0[0] = 0.0;
        d_center_0[1] = 0.0;
      }
      break;
    }

    default: { 
      TBOX_ERROR(  "PatchModule"
                << "::getFromInput()"
                << ":Invalid type of initial level set for component 0"
                << endl );
    } 
  };

  // get auxilliary parameters for initial level set 1
  switch (d_initial_level_set_1) {
    case CIRCLE: {

#ifdef LSMLIB_DOUBLE_PRECISION
      d_radius_1 = db->getDoubleWithDefault("radius_1", s_default_radius);
#else
      d_radius_1 = db->getFloatWithDefault("radius_1", s_default_radius);
#endif

      if (db->keyExists("center_1")) {
#ifdef LSMLIB_DOUBLE_PRECISION
        db->getDoubleArray("center_1", d_center_1, 2);
#else
        db->getFloatArray("center_1", d_center_1, 2);
#endif

      } else {
        d_center_1[0] = 0.0;
        d_center_1[1] = 0.0;
      }
      break;
    }

    default: { 
      TBOX_ERROR(  "PatchModule"
                << "::getFromInput()"
                << ":Invalid type of initial level set for component 0"
                << endl );
    } 
  };

}

