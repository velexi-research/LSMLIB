/*
 * File:        FieldExtensionAlgorithm.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Implementation file for level set method field extension class
 */

#ifndef included_FieldExtensionAlgorithm_cc
#define included_FieldExtensionAlgorithm_cc

// System Headers
#include <sstream>
extern "C" {
  #include <limits.h>
  #include <math.h>
}

#include "FieldExtensionAlgorithm.h" 
#include "LSMLIB_DefaultParameters.h"

// SAMRAI Headers
#include "SAMRAI/geom/CartesianPatchGeometry.h" 
#include "SAMRAI/pdat/CellData.h" 
#include "SAMRAI/pdat/CellVariable.h" 
#include "SAMRAI/hier/IntVector.h" 
#include "SAMRAI/hier/Patch.h" 
#include "SAMRAI/hier/PatchLevel.h" 
#include "SAMRAI/hier/VariableContext.h" 
#include "SAMRAI/hier/VariableDatabase.h" 

// Headers for Fortran kernels
extern "C" {
  #include "lsm_field_extension1d.h"
  #include "lsm_field_extension2d.h"
  #include "lsm_field_extension3d.h"
  #include "lsm_samrai_f77_utilities.h"
}

// SAMRAI namespaces
using namespace std;
using namespace pdat;
using namespace geom;
using namespace hier;
// Constant
#define LSM_FEA_STOP_TOLERANCE_MAX_ITERATIONS                   (1000)

namespace LSMLIB {

/* Constructor - parameters from input database */
 
FieldExtensionAlgorithm::FieldExtensionAlgorithm(
  boost::shared_ptr<Database> input_db,
  boost::shared_ptr< PatchHierarchy > hierarchy,
  const int field_handle,
  const int phi_handle,
  const int control_volume_handle,
  const IntVector& phi_ghostcell_width,
  const string& object_name)
{
  // set object_name
  d_object_name = object_name;

  // set d_patch_hierarchy and d_grid_geometry
  d_patch_hierarchy = hierarchy;
  d_grid_geometry = BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(d_patch_hierarchy->getGridGeometry());

  // set data field handles
  d_extension_field_handle = field_handle;
  d_phi_handle = phi_handle;
  d_control_volume_handle = control_volume_handle;

  // get input parameters
  getFromInput(input_db);

  // check input parameters
  checkParameters();

  // create empty BoundaryConditionModule
  boost::shared_ptr< BoundaryConditionModule >d_phi_bc_module = 
  boost::shared_ptr< BoundaryConditionModule > (new BoundaryConditionModule);
  boost::shared_ptr< BoundaryConditionModule >d_ext_field_bc_module = 
  boost::shared_ptr< BoundaryConditionModule > (new BoundaryConditionModule);

  // initialize variables and communication objects
  initializeVariables(phi_ghostcell_width);
  initializeCommunicationObjects();

}


/* Constructor - parameters from arguments */
FieldExtensionAlgorithm::FieldExtensionAlgorithm(
  boost::shared_ptr< PatchHierarchy > hierarchy,
  const tbox::Dimension& dim,
  const int field_handle,
  const int phi_handle,
  const int control_volume_handle,
  const IntVector& phi_ghostcell_width,
  const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  const int spatial_derivative_order,
  const int tvd_runge_kutta_order,
  const LSMLIB_REAL cfl_number,
  const LSMLIB_REAL stop_distance, 
  const int max_iterations,
  const LSMLIB_REAL iteration_stop_tolerance,
  const bool verbose_mode,
  const string& object_name)
{
  // set object_name
  d_object_name = object_name;

  // set d_patch_hierarchy and d_grid_geometry
  d_patch_hierarchy = hierarchy;
  d_grid_geometry = hierarchy->getGridGeometry();

  // set data field handles
  d_extension_field_handle = field_handle;
  d_phi_handle = phi_handle;
  d_control_volume_handle = control_volume_handle;

  // set integration parameters
  d_spatial_derivative_type = spatial_derivative_type;
  d_spatial_derivative_order = spatial_derivative_order;
  d_tvd_runge_kutta_order = tvd_runge_kutta_order;
  d_cfl_number = cfl_number;

  // set iteration termination parameters
  d_stop_distance = stop_distance;
  d_max_iterations = max_iterations;
  d_iteration_stop_tol = iteration_stop_tolerance;

  // set termination criteria flags 
  d_use_stop_distance = (d_stop_distance > 0.0);
  d_use_max_iterations = (d_max_iterations > 0);
  d_use_iteration_stop_tol = (d_iteration_stop_tol > 0.0);

  // if no stopping criteria are specified, use the length of the 
  // largest dimension of the computational domain as stop distance
  if ( !( d_use_iteration_stop_tol || d_use_stop_distance ||
          d_use_max_iterations) ) {
    d_use_stop_distance = true;
    
#ifdef LSMLIB_DOUBLE_PRECISION
    const double *X_lower = d_grid_geometry->getXLower();
    const double *X_upper = d_grid_geometry->getXUpper();
#else 
    const double *X_lower_double = d_grid_geometry->getXLower();
    const double *X_upper_double = d_grid_geometry->getXUpper();
    int DIM = hierarchy->getDim().getValue();
    float X_lower[DIM];
    float X_upper[DIM];
    for (int i = 0; i < DIM; i++) {
      X_lower[i] = (float) X_lower_double[i];
      X_upper[i] = (float) X_upper_double[i];
    }
#endif
    d_stop_distance = X_upper[0]-X_lower[0];
    for (int dim = 1; dim < DIM; dim++) {
      if ( d_stop_distance < X_upper[dim]-X_lower[dim] ) {
        d_stop_distance = X_upper[dim]-X_lower[dim]; 
      }
    }
  }

  // set verbose-mode
  d_verbose_mode = verbose_mode;

  // check that the user-specifeid parameters are acceptable
  checkParameters();

  // create empty BoundaryConditionModule
  boost::shared_ptr< BoundaryConditionModule >d_phi_bc_module = 
  boost::shared_ptr< BoundaryConditionModule >( new BoundaryConditionModule);
  boost::shared_ptr< BoundaryConditionModule >d_ext_field_bc_module = 
  boost::shared_ptr< BoundaryConditionModule >( new BoundaryConditionModule);

  // initialize variables and communication objects
  initializeVariables(phi_ghostcell_width);
  initializeCommunicationObjects();
}


/* computeExtensionField() */
void FieldExtensionAlgorithm::computeExtensionField(
  const int phi_component,
  const int max_iterations,
  const IntVector& lower_bc_phi,
  const IntVector& upper_bc_phi,
  const IntVector& lower_bc_ext,
  const IntVector& upper_bc_ext)
{

  // reset hierarchy configuration if necessary
  if (d_hierarchy_configuration_needs_reset) {
    resetHierarchyConfiguration(d_patch_hierarchy, 
      0, d_patch_hierarchy->getFinestLevelNumber());
  }

  // allocate patch data for requird to compute extension field
  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    boost::shared_ptr< PatchLevel > level = 
      d_patch_hierarchy->getPatchLevel(ln);
    level->allocatePatchData(d_scratch_data);
  }

  /*
   * compute dt for extension field calculation
   */

  // get dx
  int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();
  boost::shared_ptr< PatchLevel > patch_level = 
    d_patch_hierarchy->getPatchLevel(finest_level_number);
  IntVector ratio_to_coarsest = patch_level->getRatioToCoarserLevel();
  const double* coarsest_dx = d_grid_geometry->getDx();
  int DIM = d_patch_hierarchy->getDim().getValue();
  LSMLIB_REAL finest_dx[DIM];
  for (int i = 0; i < DIM; i++) {
    finest_dx[i] = coarsest_dx[i]/ratio_to_coarsest[i];
  }
  LSMLIB_REAL min_dx = finest_dx[0];
  for (int i = 1; i < DIM; i++) {
    if (finest_dx[i] < min_dx) min_dx = finest_dx[i];
  }

  // compute dt
  const LSMLIB_REAL dt = d_cfl_number*min_dx;  

  // compute the maximum number of iterations
  int num_steps = LSM_FEA_STOP_TOLERANCE_MAX_ITERATIONS;
  if (max_iterations >= 0) {
    num_steps = max_iterations;
  } else {
    if (d_use_max_iterations) num_steps = d_max_iterations; 
    if (d_use_stop_distance) {
      int stop_dist_num_steps = (int) (d_stop_distance/dt);
      if (d_use_max_iterations) {
        if (stop_dist_num_steps < num_steps) num_steps = stop_dist_num_steps;
      } else {
        num_steps = stop_dist_num_steps;
      }
    } 
  } 


  /*
   *  compute signed normal vector
   */
  // fill phi scratch space for computing signed normal vector
  if (d_phi_scr_handle != d_phi_handle) {
    LevelSetMethodToolbox::copySAMRAIData(
      d_patch_hierarchy,
      d_phi_scr_handle, d_phi_handle, 
      0, phi_component);
  }
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical 
    //       boundary conditions should be set.
    d_phi_fill_bdry_sched[ln]->fillData(0.0,true);
  }
  if (d_phi_scr_handle != d_phi_handle) {
    d_phi_bc_module->imposeBoundaryConditions(
      d_phi_scr_handle, 
      lower_bc_phi,
      upper_bc_phi,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      0);
  } else {
    d_phi_bc_module->imposeBoundaryConditions(
      d_phi_scr_handle, 
      lower_bc_phi,
      upper_bc_phi,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      phi_component);
  }

  if (d_phi_scr_handle != d_phi_handle) {  

    // case:  use FieldExtensionAlgorithm's scratch space for phi
    LevelSetMethodToolbox::computeSignedUnitNormalVectorFromPhi(
      d_patch_hierarchy,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      d_normal_vector_handle,
      d_phi_scr_handle,
      LevelSetMethodToolbox::PHI_UPWIND,
      0);  // phi scratch data is stored in component 0 

  } else {

    // case:  use user allocated data for phi
    LevelSetMethodToolbox::computeSignedUnitNormalVectorFromPhi(
      d_patch_hierarchy,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      d_normal_vector_handle,
      d_phi_handle,
      LevelSetMethodToolbox::PHI_UPWIND,
      phi_component);  // phi data to use for field extension calculation
                       // stored in phi_component 

  }
  
  /*
   *  compute the number of components in the extension field data
   *  (if it has not already been computed) 
   */
  if (d_num_field_components == 0) {
    boost::shared_ptr< PatchLevel> level = d_patch_hierarchy->getPatchLevel(0);
    for (PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) { // loop over patches
      boost::shared_ptr< Patch > patch = *pi;//returns second patch in line.
      if ( patch==NULL ) {
        TBOX_ERROR(  d_object_name
                  << "::computeExtensionField(): " 
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }
  
      boost::shared_ptr< CellData<LSMLIB_REAL> > field_data= 
        BOOST_CAST<CellData<LSMLIB_REAL>, PatchData>(
        patch->getPatchData( d_extension_field_handle ));
      d_num_field_components = field_data->getDepth();
  
      break;  // only need PatchData from one patch for computation
    }
  }

  /*
   *  main field extension loop
   */
  int count = 0;
  LSMLIB_REAL delta = 1.0;
  const int field_handle_after_step = d_extension_field_handle;
  const int field_handle_before_step = d_extension_field_scr_handles[0];
  while ( (count < num_steps) &&
          (!d_use_iteration_stop_tol || (delta > d_iteration_stop_tol)) ) {

    // reset delta to zero
    delta = 0.0;

    // loop over components in extension field
    for (int component = 0; component < d_num_field_components; component++) {

      // advance extension field equation using TVD Runge-Kutta 
      switch(d_tvd_runge_kutta_order) {
        case 1: { // first-order TVD RK (e.g. Forward Euler)
          if (d_phi_scr_handle != d_phi_handle) {  
            advanceFieldExtensionEqnUsingTVDRK1(
              dt, component, 0, lower_bc_ext, upper_bc_ext);
          } else {
            advanceFieldExtensionEqnUsingTVDRK1(
              dt, component, phi_component, lower_bc_ext, upper_bc_ext);
          }
          break;
        }
        case 2: { // second-order TVD RK 
          if (d_phi_scr_handle != d_phi_handle) {  
            advanceFieldExtensionEqnUsingTVDRK2(
              dt, component, 0, lower_bc_ext, upper_bc_ext);
          } else {
            advanceFieldExtensionEqnUsingTVDRK2(
              dt, component, phi_component, lower_bc_ext, upper_bc_ext);
          }
          break;
        }
        case 3: { // third-order TVD RK 
          if (d_phi_scr_handle != d_phi_handle) {  
            advanceFieldExtensionEqnUsingTVDRK3(
              dt, component, 0, lower_bc_ext, upper_bc_ext);
          } else {
            advanceFieldExtensionEqnUsingTVDRK3(
              dt, component, phi_component, lower_bc_ext, upper_bc_ext);
          } 
          break;
        }
        default: { // UNSUPPORTED ORDER
          TBOX_ERROR(  d_object_name
                    << "::computeExtensionField(): " 
                    << "Unsupported TVD Runge-Kutta order.  "
                    << "Only TVD-RK1, TVD-RK2, and TVD-RK3 supported."
                    << endl);
        }
      } // end switch on TVD Runge-Kutta order

      // update count and delta
      if (d_use_iteration_stop_tol) {
        delta += LevelSetMethodToolbox::maxNormOfDifference(
          d_patch_hierarchy, field_handle_after_step, field_handle_before_step, 
          d_control_volume_handle, component, 0);  // 0 is component of field
                                                   // before the time step
                                                   // which is just a single
                                                   // component scratch space
      }
    } // end loop over components of extension field

    // VERBOSE MODE
    if (d_verbose_mode) {
      pout << endl;
      pout << d_object_name << " iteration count: " << count << endl;
      if (d_use_stop_distance) { 
        pout << "  Fields on zero level set extended to a distance " 
             << "of approximately " << dt*count << endl;
      } 
      if (d_use_iteration_stop_tol) {
        pout << "  Max norm of change in extension fields: " 
             << delta << endl;
      }
    }
  
    count++;

  } // end loop over evolution of extension field equation

  // warn if iteration terminated before stop_tol reached
  if ( d_use_iteration_stop_tol && (delta > d_iteration_stop_tol) ) {
    TBOX_WARNING(  d_object_name
                << "::computeExtensionField(): "
                << "target stop tolerance (" 
                << d_iteration_stop_tol << ") NOT reached after "
                << count << " time steps. "
                << "delta = " << delta
                << endl );
  }

  // VERBOSE MODE
  if (d_verbose_mode) {
    pout << endl;
    pout << "Total number of iterations: " << count << endl;
    if (d_use_stop_distance) { 
      pout << "  Fields on zero level set extended to a distance " 
           << "of approximately " << dt*count << endl;
    } 
    if (d_use_iteration_stop_tol)
      pout << "  Last change in max norm of extension field: " 
           << delta << endl;
  }

  // deallocate patch data that was allocated to compute extension fields
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    boost::shared_ptr< PatchLevel > level 
      = d_patch_hierarchy->getPatchLevel(ln);
    level->deallocatePatchData(d_scratch_data);
  }
}


/* computeExtensionFieldForSingleComponent() */ 
void FieldExtensionAlgorithm::computeExtensionFieldForSingleComponent(
  const int component,
  const int phi_component,
  const int max_iterations,
  const IntVector& lower_bc_phi,
  const IntVector& upper_bc_phi,
  const IntVector& lower_bc_ext,
  const IntVector& upper_bc_ext)
{

  // reset hierarchy configuration if necessary
  if (d_hierarchy_configuration_needs_reset) {
    resetHierarchyConfiguration(d_patch_hierarchy, 
      0, d_patch_hierarchy->getFinestLevelNumber());
  }

  // allocate patch data for required to compute extension field
  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    boost::shared_ptr< PatchLevel > level = 
      d_patch_hierarchy->getPatchLevel(ln);
    level->allocatePatchData(d_scratch_data);
  }

  /*
   * compute dt for extension field calculation
   */

  // get dx
  int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();
  boost::shared_ptr< PatchLevel > patch_level = 
    d_patch_hierarchy->getPatchLevel(finest_level_number);
  IntVector ratio_to_coarsest = patch_level->getRatioToCoarserLevel();
  const double* coarsest_dx = d_grid_geometry->getDx();
  int DIM = d_patch_hierarchy->getDim().getValue(); 
  LSMLIB_REAL finest_dx[DIM];
  for (int i = 0; i < DIM; i++) {
    finest_dx[i] = coarsest_dx[i]/ratio_to_coarsest[i];
  }
  LSMLIB_REAL min_dx = finest_dx[0];
  for (int i = 1; i < DIM; i++) {
    if (finest_dx[i] < min_dx) min_dx = finest_dx[i];
  }

  // compute dt 
  const LSMLIB_REAL dt = d_cfl_number*min_dx;  

  // compute the maximum number of iterations
  int num_steps = LSM_FEA_STOP_TOLERANCE_MAX_ITERATIONS;
  if (max_iterations >= 0) {
    num_steps = max_iterations;
  } else {
    if (d_use_max_iterations) num_steps = d_max_iterations; 
    if (d_use_stop_distance) {
      int stop_dist_num_steps = (int) (d_stop_distance/dt);
      if (d_use_max_iterations) {
        if (stop_dist_num_steps < num_steps) num_steps = stop_dist_num_steps;
      } else {
        num_steps = stop_dist_num_steps;
      }
    } 
  } 


  /*
   *  compute signed normal vector
   */
  // fill phi scratch space for computing signed normal vector
  if (d_phi_scr_handle != d_phi_handle) {
    LevelSetMethodToolbox::copySAMRAIData(
      d_patch_hierarchy,
      d_phi_scr_handle, d_phi_handle, 
      0, phi_component);
  }
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical 
    //       boundary conditions should be set.
    d_phi_fill_bdry_sched[ln]->fillData(0.0,true);
  }
  if (d_phi_scr_handle != d_phi_handle) {
    d_phi_bc_module->imposeBoundaryConditions(
      d_phi_scr_handle, 
      lower_bc_phi,
      upper_bc_phi,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      0);
  } else {
    d_phi_bc_module->imposeBoundaryConditions(
      d_phi_scr_handle, 
      lower_bc_phi,
      upper_bc_phi,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      phi_component);
  }

  if (d_phi_scr_handle != d_phi_handle) {  

    // case:  use FieldExtensionAlgorithm's scratch space for phi
    LevelSetMethodToolbox::computeSignedUnitNormalVectorFromPhi(
      d_patch_hierarchy,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      d_normal_vector_handle,
      d_phi_scr_handle,
      LevelSetMethodToolbox::PHI_UPWIND,
      0);  // phi scratch data is stored in component 0 

  } else {

    // case:  use user allocated data for phi
    LevelSetMethodToolbox::computeSignedUnitNormalVectorFromPhi(
      d_patch_hierarchy,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      d_normal_vector_handle,
      d_phi_handle,
      LevelSetMethodToolbox::PHI_UPWIND,
      phi_component);  // phi data to use for field extension calculation
                       // stored in phi_component 

  }
  
  /*
   *  main field extension loop
   */
  int count = 0;
  LSMLIB_REAL delta = 1.0;
  const int field_handle_after_step = d_extension_field_handle;
  const int field_handle_before_step = d_extension_field_scr_handles[0];
  while ( (count < num_steps) &&
          (!d_use_iteration_stop_tol || (delta > d_iteration_stop_tol)) ) {

    // advance extension field equation using TVD Runge-Kutta 
    switch(d_tvd_runge_kutta_order) {
      case 1: { // first-order TVD RK (e.g. Forward Euler)
        if (d_phi_scr_handle != d_phi_handle) {  
          advanceFieldExtensionEqnUsingTVDRK1(
              dt, component, 0, lower_bc_ext, upper_bc_ext);
        } else {
          advanceFieldExtensionEqnUsingTVDRK1(
              dt, component, phi_component, lower_bc_ext, upper_bc_ext);
        }
        break;
      }
      case 2: { // second-order TVD RK 
        if (d_phi_scr_handle != d_phi_handle) {  
          advanceFieldExtensionEqnUsingTVDRK2(
              dt, component, 0, lower_bc_ext, upper_bc_ext);
        } else {
          advanceFieldExtensionEqnUsingTVDRK2(
              dt, component, phi_component, lower_bc_ext, upper_bc_ext);
        }
        break;
      }
      case 3: { // third-order TVD RK 
        if (d_phi_scr_handle != d_phi_handle) {  
          advanceFieldExtensionEqnUsingTVDRK3(
              dt, component, 0, lower_bc_ext, upper_bc_ext);
        } else {
          advanceFieldExtensionEqnUsingTVDRK3(
              dt, component, phi_component, lower_bc_ext, upper_bc_ext);
        } 
        break;
      }
      default: { // UNSUPPORTED ORDER
        TBOX_ERROR(  d_object_name
                  << "::computeExtensionFieldForSingleComponent(): " 
                  << "Unsupported TVD Runge-Kutta order.  "
                  << "Only TVD-RK1, TVD-RK2, and TVD-RK3 supported."
                  << endl);
      }
    } // end switch on TVD Runge-Kutta order

    // update count and delta
    if (d_use_iteration_stop_tol) {
      delta = LevelSetMethodToolbox::maxNormOfDifference(
        d_patch_hierarchy, field_handle_after_step, field_handle_before_step, 
        d_control_volume_handle, component, 0);  // 0 is component of field
                                                 // before the time step
                                                 // which is just a single
                                                 // component scratch space
    }

    // VERBOSE MODE
    if (d_verbose_mode) {
      pout << endl;
      pout << d_object_name << " iteration count: " << count << endl;
      if (d_use_stop_distance) { 
        pout << "  Fields on zero level set extended to a distance " 
             << "of approximately " << dt*count << endl;
      } 
      if (d_use_iteration_stop_tol) {
        pout << "  Max norm of change in extension field: " 
             << delta << endl;
      }
    }
  
    count++;

  } // end loop over evolution of extension field equation

  // warn if iteration terminated before stop_tol reached
  if ( d_use_iteration_stop_tol && (delta > d_iteration_stop_tol) ) {
    TBOX_WARNING(  d_object_name
                << "::computeExtensionField(): "
                << "target stop tolerance (" 
                << d_iteration_stop_tol << ") NOT reached after "
                << count << " time steps. "
                << "delta = " << delta
                << endl );
  }

  // VERBOSE MODE
  if (d_verbose_mode) {
    pout << endl;
    pout << "Total number of iterations: " << count << endl;
    if (d_use_stop_distance) { 
      pout << "  Fields on zero level set extended to a distance " 
           << "of approximately " << dt*count << endl;
    } 
    if (d_use_iteration_stop_tol)
      pout << "  Last max norm of change in extension field: " 
           << delta << endl;
  }

  // deallocate patch data that was allocated to compute extension fields
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    boost::shared_ptr< PatchLevel > level 
      = d_patch_hierarchy->getPatchLevel(ln);
    level->deallocatePatchData(d_scratch_data);
  }
}


/* resetHierarchyConfiguration() */
void FieldExtensionAlgorithm::resetHierarchyConfiguration(
  boost::shared_ptr< PatchHierarchy > hierarchy,
  const int coarsest_level,
  const int finest_level)
{
  // reset d_patch_hierarchy
  d_patch_hierarchy = hierarchy;

  // reset d_grid_geometry
  d_grid_geometry = d_patch_hierarchy->getGridGeometry();

  // compute RefineSchedules for filling extension field boundary data
  int num_levels = hierarchy->getNumberOfLevels();
  for (int k = 0; k < d_tvd_runge_kutta_order; k++) {
    d_extension_field_fill_bdry_sched[k].resizeArray(num_levels);

    for (int ln = coarsest_level; ln <= finest_level; ln++) {
      boost::shared_ptr< PatchLevel > level = hierarchy->getPatchLevel(ln);
 
      // reset data transfer configuration for boundary filling
      // before time advance
      d_extension_field_fill_bdry_sched[k][ln] =
        d_extension_field_fill_bdry_alg[k]->createSchedule(
          level, ln-1, hierarchy, 0);  // NULL RefinePatchStrategy
 
    } // end loop over levels
  } // end loop over TVD Runge-Kutta stages

  // compute RefineSchedules for filling phi boundary data 
  // (required for calculating the signed normal vector)
  d_phi_fill_bdry_sched.resizeArray(num_levels);

  for (int ln = coarsest_level; ln <= finest_level; ln++) {
    boost::shared_ptr< PatchLevel > level = hierarchy->getPatchLevel(ln);

    // reset data transfer configuration for filling phi boundary data
    // before computing signed normal vector
    d_phi_fill_bdry_sched[ln] =
      d_phi_fill_bdry_alg->createSchedule(
        level, ln-1, hierarchy, 0);  // NULL RefinePatchStrategy

  } // end loop over levels

  // create anti-periodic boundary condition module
  d_phi_bc_module->resetHierarchyConfiguration(
    d_patch_hierarchy,
    coarsest_level,
    finest_level,
    d_phi_scratch_ghostcell_width);
  d_ext_field_bc_module->resetHierarchyConfiguration(
    d_patch_hierarchy,
    coarsest_level,
    finest_level,
    d_ext_field_scratch_ghostcell_width);

  // set d_hierarchy_configuration_needs_reset to (finest_level < 0)
  d_hierarchy_configuration_needs_reset = (finest_level < 0);
}


/* advanceFieldExtensionEqnUsingTVDRK1() */
void FieldExtensionAlgorithm::advanceFieldExtensionEqnUsingTVDRK1(
  const LSMLIB_REAL dt,
  const int field_component,
  const int phi_component,
  const IntVector& lower_bc_ext,
  const IntVector& upper_bc_ext)
{
  // initialize counter for current stage of TVD RK step
  // NOTE: the rk_stage begins at 0 for convenience
  int rk_stage = 0;

  /*
   * fill scratch space for time advance
   */

  // copy component of field data to scratch space
  LevelSetMethodToolbox::copySAMRAIData(
    d_patch_hierarchy,
    d_extension_field_scr_handles[0], d_extension_field_handle, 
    0, field_component);

  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical 
    //       boundary conditions should be set.
    d_extension_field_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_ext_field_bc_module->imposeBoundaryConditions(
    d_extension_field_scr_handles[rk_stage],
    lower_bc_ext,
    upper_bc_ext,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance extension field through TVD-RK1 step
  computeFieldExtensionEqnRHS(d_extension_field_scr_handles[rk_stage],
                              phi_component);
  LevelSetMethodToolbox::TVDRK1Step(
    d_patch_hierarchy,
    d_extension_field_handle,
    d_extension_field_scr_handles[rk_stage], 
    d_rhs_handle, dt,
    field_component, 0, 0); // components of PatchData to use in TVD-RK1 step
}


/* advanceFieldExtensionEqnUsingTVDRK2() */
void FieldExtensionAlgorithm::advanceFieldExtensionEqnUsingTVDRK2(
  const LSMLIB_REAL dt,
  const int field_component,
  const int phi_component,
  const IntVector& lower_bc_ext,
  const IntVector& upper_bc_ext)
{
  // { begin Stage 1

  // initialize counter for current stage of TVD RK step
  // NOTE: the rk_stage begins at 0 for convenience
  int rk_stage = 0;

  /*
   * fill scratch space for first stage of time advance
   */

  // copy component of field data to scratch space
  LevelSetMethodToolbox::copySAMRAIData(
    d_patch_hierarchy,
    d_extension_field_scr_handles[0], d_extension_field_handle, 
    0, field_component);

  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical 
    //       boundary conditions should be set.
    d_extension_field_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_ext_field_bc_module->imposeBoundaryConditions(
    d_extension_field_scr_handles[rk_stage],
    lower_bc_ext,
    upper_bc_ext,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance extension field through the first stage of TVD-RK2
  computeFieldExtensionEqnRHS(d_extension_field_scr_handles[rk_stage],
                              phi_component);
  LevelSetMethodToolbox::TVDRK2Stage1(
    d_patch_hierarchy,
    d_extension_field_scr_handles[rk_stage+1],
    d_extension_field_scr_handles[rk_stage],
    d_rhs_handle, dt,
    0, 0, 0);  // components of PatchData to use in TVD-RK2 step

  // } end Stage 1


  // { begin Stage 2

  // advance TVD RK2 stage counter
  rk_stage = 1;

  // fill scratch space for second stage of time advance
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical 
    //       boundary conditions should be set.
    d_extension_field_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_ext_field_bc_module->imposeBoundaryConditions(
    d_extension_field_scr_handles[rk_stage],
    lower_bc_ext,
    upper_bc_ext,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance extension field through the second stage of TVD-RK2
  computeFieldExtensionEqnRHS(d_extension_field_scr_handles[rk_stage],
                              phi_component);
  LevelSetMethodToolbox::TVDRK2Stage2(
    d_patch_hierarchy,
    d_extension_field_handle,
    d_extension_field_scr_handles[rk_stage],
    d_extension_field_scr_handles[0],
    d_rhs_handle, dt,
    field_component, 0, 0, 0);  // components of PatchData to use in 
                                // TVD-RK2 step

  // } end Stage 2
}


/* advanceFieldExtensionEqnUsingTVDRK3() */
void FieldExtensionAlgorithm::advanceFieldExtensionEqnUsingTVDRK3(
  const LSMLIB_REAL dt,
  const int field_component,
  const int phi_component,
  const IntVector& lower_bc_ext,
  const IntVector& upper_bc_ext)
{
  // { begin Stage 1

  // initialize counter for current stage of TVD RK step
  // NOTE: the rk_stage begins at 0 for convenience
  int rk_stage = 0;

  /*
   * fill scratch space for first stage of time advance
   */

  // copy component of field data to scratch space
  LevelSetMethodToolbox::copySAMRAIData(
    d_patch_hierarchy,
    d_extension_field_scr_handles[0], d_extension_field_handle, 
    0, field_component);

  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical 
    //       boundary conditions should be set.
    d_extension_field_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_ext_field_bc_module->imposeBoundaryConditions(
    d_extension_field_scr_handles[rk_stage],
    lower_bc_ext,
    upper_bc_ext,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance extension field through the first stage of TVD-RK3
  computeFieldExtensionEqnRHS(d_extension_field_scr_handles[rk_stage],
                              phi_component);
  LevelSetMethodToolbox::TVDRK3Stage1(
    d_patch_hierarchy,
    d_extension_field_scr_handles[rk_stage+1],
    d_extension_field_scr_handles[rk_stage],
    d_rhs_handle, dt,
    0, 0, 0);  // components of PatchData to use in TVD-RK3 step

  // } end Stage 1


  // { begin Stage 2

  // advance TVD RK3 stage counter
  rk_stage = 1;

  // fill scratch space for second stage of time advance
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical 
    //       boundary conditions should be set.
    d_extension_field_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_ext_field_bc_module->imposeBoundaryConditions(
    d_extension_field_scr_handles[rk_stage],
    lower_bc_ext,
    upper_bc_ext,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance extension field through the second stage of TVD-RK3
  computeFieldExtensionEqnRHS(d_extension_field_scr_handles[rk_stage],
                              phi_component);
  LevelSetMethodToolbox::TVDRK3Stage2(
    d_patch_hierarchy,
    d_extension_field_scr_handles[rk_stage+1],
    d_extension_field_scr_handles[rk_stage],
    d_extension_field_scr_handles[rk_stage-1],
    d_rhs_handle, dt,
    0, 0, 0, 0);  // components of PatchData to use in TVD-RK3 step

  // } end Stage 2


  // { begin Stage 3

  // advance TVD RK3 stage counter
  rk_stage = 2;

  // fill scratch space for second stage of time advance
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical 
    //       boundary conditions should be set.
    d_extension_field_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_ext_field_bc_module->imposeBoundaryConditions(
    d_extension_field_scr_handles[rk_stage],
    lower_bc_ext,
    upper_bc_ext,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance extension field through the third stage of TVD-RK3
  computeFieldExtensionEqnRHS(d_extension_field_scr_handles[rk_stage],
                              phi_component);
  LevelSetMethodToolbox::TVDRK3Stage3(
    d_patch_hierarchy,
    d_extension_field_handle,
    d_extension_field_scr_handles[rk_stage],
    d_extension_field_scr_handles[0],
    d_rhs_handle, dt,
    field_component, 0, 0, 0);  // components of PatchData to use in 
                                // TVD-RK3 step

  // } end Stage 3
}


/* computeFieldExtensionEqnRHS() */
template <int DIM> 
void FieldExtensionAlgorithm<DIM>::computeFieldExtensionEqnRHS(
  const int extension_field_handle,
  const int phi_component)
{
  // compute spatial derivatives of the extension field for 
  // the current stage
  LevelSetMethodToolbox<DIM>::computeUpwindSpatialDerivatives(
    d_patch_hierarchy,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    d_grad_field_handle,
    extension_field_handle,
    d_normal_vector_handle);

  // loop over PatchHierarchy and compute RHS for level set equation
  // by calling Fortran routines
  const int num_levels = d_patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = d_patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  d_object_name
                  << "::computeFieldExtensionRHS(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get grid spacing
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geom->getDx();
#else 
      const double* dx_double = patch_geom->getDx();
      float dx[DIM]; 
      for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,LSMLIB_REAL> > rhs_data =
        patch->getPatchData( d_rhs_handle );
      Pointer< CellData<DIM,LSMLIB_REAL> > field_data =
        patch->getPatchData( extension_field_handle );
      Pointer< CellData<DIM,LSMLIB_REAL> > phi_data =
        patch->getPatchData( d_phi_scr_handle );
      Pointer< CellData<DIM,LSMLIB_REAL> > normal_vector_data =
        patch->getPatchData( d_normal_vector_handle );
      Pointer< CellData<DIM,LSMLIB_REAL> > grad_field_data =
        patch->getPatchData( d_grad_field_handle );

      Box<DIM> rhs_ghostbox = rhs_data->getGhostBox();
      const IntVector<DIM> rhs_ghostbox_lower = rhs_ghostbox.lower();
      const IntVector<DIM> rhs_ghostbox_upper = rhs_ghostbox.upper();

      Box<DIM> field_ghostbox = field_data->getGhostBox();
      const IntVector<DIM> field_ghostbox_lower = field_ghostbox.lower();
      const IntVector<DIM> field_ghostbox_upper = field_ghostbox.upper();

      Box<DIM> phi_ghostbox = phi_data->getGhostBox();
      const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
      const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();

      Box<DIM> normal_vector_ghostbox = normal_vector_data->getGhostBox();
      const IntVector<DIM> normal_vector_ghostbox_lower = 
        normal_vector_ghostbox.lower();
      const IntVector<DIM> normal_vector_ghostbox_upper = 
        normal_vector_ghostbox.upper();

      Box<DIM> grad_field_ghostbox = grad_field_data->getGhostBox();
      const IntVector<DIM> grad_field_ghostbox_lower = 
        grad_field_ghostbox.lower();
      const IntVector<DIM> grad_field_ghostbox_upper = 
        grad_field_ghostbox.upper();

      // fill box
      Box<DIM> fillbox = rhs_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      LSMLIB_REAL* rhs = rhs_data->getPointer();
      LSMLIB_REAL* field = field_data->getPointer();
      LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
      LSMLIB_REAL* normal_vector[LSM_DIM_MAX];
      LSMLIB_REAL* upwind_grad_field[LSM_DIM_MAX];
      for (int dim = 0; dim < DIM; dim++) {
        upwind_grad_field[dim] = grad_field_data->getPointer(dim);
        normal_vector[dim] = normal_vector_data->getPointer(dim);
      }

      if (DIM == 3) {

        LSM3D_COMPUTE_FIELD_EXTENSION_EQN_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2],
          field,
          &field_ghostbox_lower[0],
          &field_ghostbox_upper[0],
          &field_ghostbox_lower[1],
          &field_ghostbox_upper[1],
          &field_ghostbox_lower[2],
          &field_ghostbox_upper[2],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          &phi_ghostbox_lower[2],
          &phi_ghostbox_upper[2],
          upwind_grad_field[0], upwind_grad_field[1], upwind_grad_field[2],
          &grad_field_ghostbox_lower[0],
          &grad_field_ghostbox_upper[0],
          &grad_field_ghostbox_lower[1],
          &grad_field_ghostbox_upper[1],
          &grad_field_ghostbox_lower[2],
          &grad_field_ghostbox_upper[2],
          normal_vector[0], normal_vector[1], normal_vector[2],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &normal_vector_ghostbox_lower[1],
          &normal_vector_ghostbox_upper[1],
          &normal_vector_ghostbox_lower[2],
          &normal_vector_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2],
          &dx[0], &dx[1], &dx[2]);

      } else if (DIM == 2) {

        LSM2D_COMPUTE_FIELD_EXTENSION_EQN_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          field,
          &field_ghostbox_lower[0],
          &field_ghostbox_upper[0],
          &field_ghostbox_lower[1],
          &field_ghostbox_upper[1],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          upwind_grad_field[0], upwind_grad_field[1], 
          &grad_field_ghostbox_lower[0],
          &grad_field_ghostbox_upper[0],
          &grad_field_ghostbox_lower[1],
          &grad_field_ghostbox_upper[1],
          normal_vector[0], normal_vector[1], 
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &normal_vector_ghostbox_lower[1],
          &normal_vector_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &dx[0], &dx[1]);

      } else if (DIM == 1) {

        LSM1D_COMPUTE_FIELD_EXTENSION_EQN_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          field,
          &field_ghostbox_lower[0],
          &field_ghostbox_upper[0],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          upwind_grad_field[0], 
          &grad_field_ghostbox_lower[0],
          &grad_field_ghostbox_upper[0],
          normal_vector[0],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &dx[0]);

      } else {  // Unsupported dimension
        TBOX_ERROR(  d_object_name
                  << "::computeFieldExtensionRHS(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

    } // end loop over patches in level
  } // end loop over levels in hierarchy

}


/* initializeVariables() */
template <int DIM>
void FieldExtensionAlgorithm<DIM>::initializeVariables(
  const IntVector<DIM>& phi_ghostcell_width)
{
  // initialize d_num_field_components to zero
  d_num_field_components = 0;

  // compute ghost cell widths
  int scratch_ghostcell_width_for_grad = -1; // bogus value which is reset in 
                                             // switch statement
  switch (d_spatial_derivative_type) { 
    case ENO: {
      scratch_ghostcell_width_for_grad = d_spatial_derivative_order;
      break;
    } 
    case WENO: {
      scratch_ghostcell_width_for_grad = d_spatial_derivative_order/2 + 1;
      break;
    }
    default:
      TBOX_ERROR(  d_object_name
                << "::initializeVariables(): "
                << "Unsupported spatial derivative type.  "
                << "Only ENO and WENO derivatives are supported."
                << endl );
  }

  IntVector<DIM> scratch_ghostcell_width(scratch_ghostcell_width_for_grad);
  IntVector<DIM> zero_ghostcell_width(0);

  /*
   * create variables and PatchData for scratch data
   */
  VariableDatabase<DIM> *var_db = VariableDatabase<DIM>::getDatabase();
  
  // compute minimum ghostcell width for phi data 
  int min_ghostcell_width = phi_ghostcell_width(0);
  for (int k = 1; k < DIM; k++ ) {
    if (min_ghostcell_width > phi_ghostcell_width(k)) 
      min_ghostcell_width = phi_ghostcell_width(k);
  }
  
  // get variable associated with phi_handle
  Pointer< Variable<DIM> > tmp_variable;
  Pointer<VariableContext> tmp_context;
  Pointer< CellVariable<DIM,LSMLIB_REAL> > phi_variable;
  if (var_db->mapIndexToVariableAndContext(d_phi_handle, 
                                           tmp_variable, tmp_context)) {
    phi_variable = tmp_variable;
  } else {
    TBOX_ERROR(  d_object_name
              << "::initializeVariables(): "
              << "Specified phi handle does not exist or does "
              << "not correspond to cell-centered data"
              << endl);
  }

  // create scratch context for extension field calculation
  Pointer<VariableContext> scratch_context = 
    var_db->getContext("EXTENSION_FIELD_SCRATCH");

  // reserve space for extension field scratch PatchData handles
  d_extension_field_scr_handles.reserve(d_tvd_runge_kutta_order);

  // clear out ComponentSelector for scratch data
  d_scratch_data.clrAllFlags();

  // get CellVariable associated with d_extension_field_handle
  Pointer< CellVariable<DIM,LSMLIB_REAL> > field_variable;
  if (var_db->mapIndexToVariableAndContext(d_extension_field_handle, 
                                           tmp_variable, tmp_context)) {
    field_variable = tmp_variable;
  } else {
    TBOX_ERROR(  d_object_name
              << "::initializeVariables(): "
              << "Specified field handle does not exist or does "
              << "not correspond to cell-centered data"
              << endl);
  }

  // create "SCRATCH" context for extension field 
  d_ext_field_scratch_ghostcell_width = scratch_ghostcell_width;
  stringstream ext_field_scratch_name("");
  ext_field_scratch_name << d_object_name << "::" << field_variable->getName() 
                   << "::EXTENSION_FIELD_SCRATCH"; 
  Pointer< CellVariable<DIM,LSMLIB_REAL> > ext_field_scratch_variable;
  if (var_db->checkVariableExists(ext_field_scratch_name.str())) {
    ext_field_scratch_variable = 
      var_db->getVariable(ext_field_scratch_name.str());
  } else {
    ext_field_scratch_variable = 
      new CellVariable<DIM,LSMLIB_REAL>(ext_field_scratch_name.str(), 1);
  }
  for (int k=0; k < d_tvd_runge_kutta_order; k++) {
    stringstream context_name("");
    context_name << "EXTENSION_FIELD_TVDRK::" 
                 << d_extension_field_handle
                 << "::" << k;
    d_extension_field_scr_handles[k] = 
      var_db->registerVariableAndContext(
        ext_field_scratch_variable,
        var_db->getContext(context_name.str()),
        d_ext_field_scratch_ghostcell_width);
    d_scratch_data.setFlag(d_extension_field_scr_handles[k]);
  }

  // create "SCRATCH" context for phi if there are insufficient ghostcells 
  // for grad(phi) calculation
  if ( (min_ghostcell_width < scratch_ghostcell_width_for_grad) ||
       (min_ghostcell_width < 1) ) {  // at least 1 ghostcell is required
                                      // to compute upwind normal vector

    stringstream phi_scratch_name("");
    phi_scratch_name << d_object_name << "::" << phi_variable->getName() 
                     << "::EXTENSION_FIELD_PHI_SCRATCH"; 
    Pointer< CellVariable<DIM,LSMLIB_REAL> > phi_scratch_variable;
    if (var_db->checkVariableExists(phi_scratch_name.str())) {
      phi_scratch_variable = var_db->getVariable(phi_scratch_name.str());
    } else {
      phi_scratch_variable = 
        new CellVariable<DIM,LSMLIB_REAL>(phi_scratch_name.str(), 1);
    }
    d_phi_scr_handle = var_db->registerVariableAndContext(
      phi_scratch_variable, scratch_context, scratch_ghostcell_width);
    d_scratch_data.setFlag(d_phi_scr_handle);
    d_phi_scratch_ghostcell_width = scratch_ghostcell_width;
  } else { // case: phi_handle has enough ghostcells for grad(phi) calculation
    d_phi_scr_handle = d_phi_handle;
    d_phi_scratch_ghostcell_width = phi_ghostcell_width;
  }

  // create RHS variables
  stringstream rhs_name("");
  rhs_name << field_variable->getName() << "::EXTENSION_FIELD_RHS";
  Pointer< CellVariable<DIM,LSMLIB_REAL> > rhs_variable;
  if (var_db->checkVariableExists(rhs_name.str())) {
   rhs_variable = var_db->getVariable(rhs_name.str());
  } else {
   rhs_variable = new CellVariable<DIM,LSMLIB_REAL>(rhs_name.str(), 1);
  }
  d_rhs_handle = var_db->registerVariableAndContext(
    rhs_variable, scratch_context, zero_ghostcell_width);
  d_scratch_data.setFlag(d_rhs_handle);

  // create variables for unit normal vector ( grad(phi)/|grad(phi)| )
  stringstream normal_vector_name("");
  normal_vector_name << phi_variable->getName() 
                     << "::EXTENSION_FIELD_NORMAL_VECTOR";
  Pointer< CellVariable<DIM,LSMLIB_REAL> > normal_vector_variable;
  if (var_db->checkVariableExists(normal_vector_name.str())) {
   normal_vector_variable = var_db->getVariable(normal_vector_name.str());
  } else {
   normal_vector_variable = 
     new CellVariable<DIM,LSMLIB_REAL>(normal_vector_name.str(), DIM);
  }
  d_normal_vector_handle = var_db->registerVariableAndContext(
    normal_vector_variable, scratch_context, zero_ghostcell_width);
  d_scratch_data.setFlag(d_normal_vector_handle);

  // create variables for grad(phi)
  stringstream grad_phi_name("");
  grad_phi_name << phi_variable->getName() 
                << "::EXTENSION_FIELD_GRAD_PHI";
  Pointer< CellVariable<DIM,LSMLIB_REAL> > grad_phi_variable;
  if (var_db->checkVariableExists(grad_phi_name.str())) {
   grad_phi_variable = var_db->getVariable(grad_phi_name.str());
  } else {
   grad_phi_variable = new CellVariable<DIM,LSMLIB_REAL>(grad_phi_name.str(), DIM);
  }
  Pointer<VariableContext> grad_phi_plus_context =  
    var_db->getContext("EXTENSION_FIELD_GRAD_PHI_PLUS");
  Pointer<VariableContext> grad_phi_minus_context =  
    var_db->getContext("EXTENSION_FIELD_GRAD_PHI_MINUS");
  d_grad_phi_plus_handle = var_db->registerVariableAndContext(
    grad_phi_variable, grad_phi_plus_context, zero_ghostcell_width);
  d_grad_phi_minus_handle = var_db->registerVariableAndContext(
    grad_phi_variable, grad_phi_minus_context, zero_ghostcell_width);
  d_scratch_data.setFlag(d_grad_phi_plus_handle);
  d_scratch_data.setFlag(d_grad_phi_minus_handle);

  // set up grad(field) handle
  stringstream grad_field_name("");
  grad_field_name << field_variable->getName() << "::" 
                  << phi_variable->getName()
                  << "::EXTENSION_FIELD_GRAD_FIELD";
  Pointer< CellVariable<DIM,LSMLIB_REAL> > grad_field_variable;
  if (var_db->checkVariableExists(grad_field_name.str())) {
   grad_field_variable = var_db->getVariable(grad_field_name.str());
  } else {
   grad_field_variable = 
     new CellVariable<DIM,LSMLIB_REAL>(grad_field_name.str(), DIM);
  }
  d_grad_field_handle = var_db->registerVariableAndContext(
    grad_field_variable, 
    var_db->getContext("EXTENSION_FIELD_UPWIND_GRAD_FIELD"), 
    zero_ghostcell_width);
  d_scratch_data.setFlag(d_grad_field_handle);

}


/* initializeCommunicationObjects() */
template <int DIM>
void FieldExtensionAlgorithm<DIM>::initializeCommunicationObjects()
{

  // initialize d_hierarchy_configuration_needs_reset to true
  d_hierarchy_configuration_needs_reset = true;

  // get pointer to VariableDatabase
  VariableDatabase<DIM> *var_db = VariableDatabase<DIM>::getDatabase();

  /*
   * Look up refine operation
   */
  // get CellVariable associated with d_extension_field_handle
  Pointer< CellVariable<DIM,LSMLIB_REAL> > field_variable;
  Pointer< Variable<DIM> > tmp_variable;
  Pointer<VariableContext> tmp_context;
  if (var_db->mapIndexToVariableAndContext(d_extension_field_handle,
                                           tmp_variable, tmp_context)) {
    field_variable = tmp_variable;
  } else {
    TBOX_ERROR(  d_object_name
              << "::initializeCommunicationObjects(): "
              << "Specified field handle does not exist or does "
              << "not correspond to cell-centered data"
              << endl);
  }

  // lookup refine operations
  Pointer< RefineOperator<DIM> > refine_op =
    d_grid_geometry->lookupRefineOperator(field_variable, "LINEAR_REFINE");


  /*
   * create RefineAlgorithms for filling boundary data for extension fields
   */
  d_extension_field_fill_bdry_alg.resizeArray(d_tvd_runge_kutta_order);
  d_extension_field_fill_bdry_sched.resizeArray(d_tvd_runge_kutta_order);
  for (int k = 0; k < d_tvd_runge_kutta_order; k++) {
    d_extension_field_fill_bdry_alg[k] = new RefineAlgorithm<DIM>;

    // empty out the boundary bdry fill schedules
    d_extension_field_fill_bdry_sched[k].setNull();

    // register data transfer to fill boundary data before time advance
    d_extension_field_fill_bdry_alg[k]->registerRefine(
      d_extension_field_scr_handles[k],
      d_extension_field_scr_handles[k],
      d_extension_field_scr_handles[k],
      refine_op);

  } // end loop over TVD-Runge-Kutta stages

  // create RefineAlgorithms for filling boundary data for phi 
  // (required to calculate the signed normal vector)
  d_phi_fill_bdry_alg = new RefineAlgorithm<DIM>;

  // empty out the boundary bdry fill schedules
  d_phi_fill_bdry_sched.setNull();

  // fill algorithm before computing extension fields
  d_phi_fill_bdry_alg->registerRefine(
    d_phi_scr_handle,
    d_phi_scr_handle,
    d_phi_scr_handle,
    refine_op);

  // configure communications schedules for ALL levels using the 
  // specified PatchHierarchy
  resetHierarchyConfiguration(d_patch_hierarchy, 
    0, d_patch_hierarchy->getFinestLevelNumber());

}


/* getFromInput() */
template <int DIM> 
void FieldExtensionAlgorithm<DIM>::getFromInput(
  Pointer<Database> db)
{

  // get numerical parameters
  string spatial_derivative_type = db->getStringWithDefault(
    "spatial_derivative_type", LSM_DEFAULT_SPATIAL_DERIVATIVE_TYPE);
  if (spatial_derivative_type == "WENO") {
    d_spatial_derivative_type = WENO;
    d_spatial_derivative_order = db->getIntegerWithDefault(
      "spatial_derivative_order", LSM_DEFAULT_SPATIAL_DERIVATIVE_WENO_ORDER);
  } else if (spatial_derivative_type == "ENO") {
    d_spatial_derivative_type = ENO;
    d_spatial_derivative_order = db->getIntegerWithDefault(
      "spatial_derivative_order", LSM_DEFAULT_SPATIAL_DERIVATIVE_ENO_ORDER);
  } else {
    d_spatial_derivative_type = UNKNOWN;
  }
  d_tvd_runge_kutta_order = db->getIntegerWithDefault(
    "tvd_runge_kutta_order", LSM_DEFAULT_TVD_RUNGE_KUTTA_ORDER);
  d_cfl_number = db->getDoubleWithDefault("cfl_number",
                                          LSM_DEFAULT_CFL_NUMBER);

  // set iteration termination parameters
  d_stop_distance = db->getDoubleWithDefault("stop_distance", 0.0);
  d_max_iterations = db->getIntegerWithDefault("max_iterations", 0);
  d_iteration_stop_tol = db->getDoubleWithDefault(
      "iteration_stop_tolerance", 0.0);

  // set termination criteria flags 
  d_use_stop_distance = (d_stop_distance > 0.0);
  d_use_max_iterations = (d_max_iterations > 0);
  d_use_iteration_stop_tol = (d_iteration_stop_tol > 0.0);

  // if no stopping criteria are specified, use the length of the 
  // largest dimension of the computational domain as stop distance
  if ( !( d_use_iteration_stop_tol || d_use_stop_distance ||
          d_use_max_iterations) ) {
    d_use_stop_distance = true;
    
#ifdef LSMLIB_DOUBLE_PRECISION
    const double *X_lower = d_grid_geometry->getXLower();
    const double *X_upper = d_grid_geometry->getXUpper();
#else
    const double *X_lower_double = d_grid_geometry->getXLower();
    const double *X_upper_double = d_grid_geometry->getXUpper();
    float X_lower[DIM];
    float X_upper[DIM];
    for (int i = 0; i < DIM; i++) {
      X_lower[i] = (float) X_lower_double[i];
      X_upper[i] = (float) X_upper_double[i];
    }
#endif

    d_stop_distance = X_upper[0]-X_lower[0];
    for (int dim = 1; dim < DIM; dim++) {
      if ( d_stop_distance < X_upper[dim]-X_lower[dim] ) {
        d_stop_distance = X_upper[dim]-X_lower[dim]; 
      }
    }
  }

  // get verbose mode
  d_verbose_mode = db->getBoolWithDefault(
    "verbose_mode", LSM_DEFAULT_VERBOSE_MODE);

}


/* checkParameters() */
template <int DIM> 
void FieldExtensionAlgorithm<DIM>::checkParameters()
{
  // check that spatial derivative type, spatial derivative order,
  // and TVD Runge-Kutta order are valid.
  if (d_spatial_derivative_type == ENO) {
    if ( (d_spatial_derivative_order < 1) ||
         (d_spatial_derivative_order > 3) ) {
      TBOX_ERROR(  d_object_name
                << "::checkParameters(): "
                << "Unsupported order for ENO derivative.  "
                << "Only ENO1, ENO2, and ENO3 supported."
                << endl );
    }
  } else if (d_spatial_derivative_type == WENO) {
    if ( (d_spatial_derivative_order != 5) )
      TBOX_ERROR(  d_object_name
                << "::checkParameters(): "
                << "Unsupported order for WENO derivative.  "
                << "Only WENO5 supported."
                << endl );
  }
  if ( (d_tvd_runge_kutta_order < 1) ||
       (d_tvd_runge_kutta_order > 3) ) {
    TBOX_ERROR(  d_object_name
              << "::checkParameters(): "
              << "Unsupported TVD Runge-Kutta order.  "
              << "Only TVD-RK1, TVD-RK2, and TVD-RK3 supported."
              << endl );
  }
}

} // end LSMLIB namespace

#endif
