/*
 * File:        VelocityFieldModule.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.4 $
 * Modified:    $Date$
 * Description: Implementation of class that computes the normal velocity 
 *              field for the level set method
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

/*
 * class constants
 */
const int VelocityFieldModule::s_max_scratch_ghostcell_width = 3;


/* Constructor */
VelocityFieldModule::VelocityFieldModule(
  Pointer<Database> input_db,
  Pointer< PatchHierarchy<2> > patch_hierarchy,
  Pointer< CartesianGridGeometry<2> > grid_geom,
  const string& object_name) :
d_min_dt_computed(false)
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

  // Allocate normal velocity variable
  Pointer< CellVariable<2,LSMLIB_REAL> > normal_velocity = 
    new CellVariable<2,LSMLIB_REAL>("normal velocity field",1); 
 
  // Register normal velocity variable with VariableDatabase.
  VariableDatabase<2> *vdb = VariableDatabase<2>::getDatabase();
  Pointer<VariableContext> cur_ctxt = vdb->getContext("CURRENT");
  d_normal_velocity_handle = vdb->registerVariableAndContext(
    normal_velocity, cur_ctxt, IntVector<2>(0));
  vdb->registerPatchDataForRestart(d_normal_velocity_handle);

  // initialize d_velocity_never_computed to true to ensure that velocity
  // is computed on first call to computeVelocityField()
  d_velocity_never_computed = true;
}


/* computeStableDt() */
LSMLIB_REAL VelocityFieldModule::computeStableDt()
{
  if (!d_min_dt_computed) {
    switch (d_normal_velocity_field_selector) {
      case EXPANSION_VELOCITY: { 
        // make sure there are approximately 10 time steps per period
        if (2.0/d_omega < d_min_dt) d_min_dt = 2.0/d_omega;
        break;
      }
      case MEAN_CURVATURE_VELOCITY: {
        // use approximate stable time step for diffusion
        const int finest_level_num = d_patch_hierarchy->getFinestLevelNumber();
        Pointer< PatchLevel<2> > level = 
          d_patch_hierarchy->getPatchLevel(finest_level_num);
        IntVector<2> ratio_to_coarsest = level->getRatio();
        const double* dx_coarsest = d_grid_geometry->getDx();
        LSMLIB_REAL dx[2];
        dx[0] = dx_coarsest[0]/ratio_to_coarsest[0];
        dx[1] = dx_coarsest[1]/ratio_to_coarsest[1];
        LSMLIB_REAL one_over_dx_sq = 1.0/dx[0]/dx[0] + 1.0/dx[1]/dx[1];
        if (0.5/d_speed/one_over_dx_sq < d_min_dt) 
          d_min_dt = 0.5/d_speed/one_over_dx_sq; 
      }
  
      default: {}
    }; // switch statement

    // set d_min_dt_computed to true
    d_min_dt_computed = true;
  }
  return d_min_dt;
}


/* computeVelocityField() */
void VelocityFieldModule::computeVelocityField(
  const LSMLIB_REAL time,
  const int phi_handle,
  const int psi_handle,
  const int component)
{
  (void) psi_handle;  // psi not relevant in 2D
  (void) component;  // component is not used because this example problem
                     // only has one component for level set function

  // only carry out computation if the time has changed
  if (!d_velocity_never_computed && (d_current_time == time)) return;
  
  // set d_velocity_never_computed to false
  d_velocity_never_computed = false;

  // update the current time
  d_current_time = time;

  // set normal velocity on all levels of hierarchy
  const int finest_level = d_patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    Pointer< PatchLevel<2> > level = d_patch_hierarchy->getPatchLevel(ln);
    computeVelocityFieldOnLevel(level,time,phi_handle);

  } // end loop over hierarchy
}


/* initializeLevelData() */
void VelocityFieldModule::initializeLevelData (
  const Pointer< PatchHierarchy<2> > hierarchy ,
  const int level_number ,
  const LSMLIB_REAL init_data_time ,
  const int phi_handle,
  const int psi_handle,
  const bool can_be_refined ,
  const bool initial_time ,
  const Pointer< PatchLevel<2> > old_level,
  const bool allocate_data)
{
//  (void) psi_handle;  // psi is meaningless for 2D problems

  Pointer< PatchLevel<2> > level = hierarchy->getPatchLevel(level_number);
  if (allocate_data) {
    level->allocatePatchData(d_normal_velocity_handle);
  }

  /*
   * Initialize data on all patches in the level.
   */ 
  computeVelocityFieldOnLevel(level,init_data_time,phi_handle);

}

/* computeVelocityFieldOnLevel() */
void VelocityFieldModule::computeVelocityFieldOnLevel(
  const Pointer< PatchLevel<2> > level,
  const LSMLIB_REAL time,
  const int phi_handle)
{
  for (PatchLevelIterator<2> pi(level); pi; pi++) { // loop over patches
    const int pn = *pi;
    Pointer< Patch<2> > patch = level->getPatch(pn);
    if ( patch.isNull() ) {
      TBOX_ERROR(d_object_name << ": Cannot find patch. Null patch pointer.");
    }

    // get geometry 
    Pointer< CartesianPatchGeometry<2> > patch_geom 
      = patch->getPatchGeometry();
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

    // get field data
    Pointer< CellData<2,LSMLIB_REAL> > normal_velocity_data = 
      patch->getPatchData( d_normal_velocity_handle );
    Pointer< CellData<2,LSMLIB_REAL> > phi_data = 
      patch->getPatchData( phi_handle );

    Box<2> normal_vel_ghostbox = normal_velocity_data->getGhostBox();
    const IntVector<2> normal_vel_ghostbox_lower = normal_vel_ghostbox.lower();
    const IntVector<2> normal_vel_ghostbox_upper = normal_vel_ghostbox.upper();

    Box<2> phi_ghostbox = phi_data->getGhostBox();
    const IntVector<2> phi_ghostbox_lower = phi_ghostbox.lower();
    const IntVector<2> phi_ghostbox_upper = phi_ghostbox.upper();

    Box<2> normal_vel_box = normal_velocity_data->getBox();
    const IntVector<2> normal_vel_lower = normal_vel_box.lower();
    const IntVector<2> normal_vel_upper = normal_vel_box.upper();

    // get normal velocity data pointer
    LSMLIB_REAL* normal_vel_data_ptr = normal_velocity_data->getPointer();
    LSMLIB_REAL* phi_data_ptr = phi_data->getPointer();

    switch (d_normal_velocity_field_selector) {
      case EXPANSION_VELOCITY: { 
        // oscillating expanding/contracting normal velocity field 
        // V_n = speed*cos(omega*time) 
        COMPUTE_OSCILLATING_EXPANSION_VELOCITY(
          normal_vel_data_ptr,
          &normal_vel_ghostbox_lower[0],
          &normal_vel_ghostbox_upper[0],
          &normal_vel_ghostbox_lower[1],
          &normal_vel_ghostbox_upper[1],
          &normal_vel_lower[0],
          &normal_vel_upper[0],
          &normal_vel_lower[1],
          &normal_vel_upper[1],
          x_lower,
          &d_speed,
          &d_omega,
          &time);
        break;
      }

      case MEAN_CURVATURE_VELOCITY: {
        // mean curvature velocity
        // V_n = -speed*kappa
        LSMLIB_REAL speed = 0.1;
        COMPUTE_MEAN_CURVATURE_VELOCITY(
          normal_vel_data_ptr,
          &normal_vel_ghostbox_lower[0],
          &normal_vel_ghostbox_upper[0],
          &normal_vel_ghostbox_lower[1],
          &normal_vel_ghostbox_upper[1],
          phi_data_ptr,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          &normal_vel_lower[0],
          &normal_vel_upper[0],
          &normal_vel_lower[1],
          &normal_vel_upper[1],
          &dx[0],
          &dx[1],
          &speed);
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
  os << "d_normal_velocity_field = " 
     << d_normal_velocity_field_selector << endl;

  // KTC - put more here...
  os << endl;
}

void VelocityFieldModule::getFromInput(
  Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!db.isNull());
#endif

  // set velocity type
  d_normal_velocity_field_selector = 
    db->getIntegerWithDefault("normal_velocity_field", 0);

  // set d_speed 
#ifdef LSMLIB_DOUBLE_PRECISION
  d_speed = db->getDoubleWithDefault("speed", 1.0);
#else
  d_speed = db->getFloatWithDefault("speed", 1.0);
#endif

  // set d_omega
#ifdef LSMLIB_DOUBLE_PRECISION
  d_omega = db->getDoubleWithDefault("omega", 1.0);
#else
  d_omega = db->getFloatWithDefault("omega", 1.0);
#endif


  // get min_dt from input
#ifdef LSMLIB_DOUBLE_PRECISION
  d_min_dt = db->getDoubleWithDefault("min_dt", DBL_MAX);
#else
  d_min_dt = db->getFloatWithDefault("min_dt", FLT_MAX);
#endif

}

