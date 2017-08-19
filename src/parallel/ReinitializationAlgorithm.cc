/*
 * File:        ReinitializationAlgorithm.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Implementation file for level set method reinitialization class
 */

#ifndef included_ReinitializationAlgorithm_cc
#define included_ReinitializationAlgorithm_cc

#include "ReinitializationAlgorithm.h"

// Standard library headers
#include <cstddef>
#include <sstream>

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

// LSMLIB headers
#include "LSMLIB_config.h"
#include "LSMLIB_DefaultParameters.h"
#include "BoundaryConditionModule.h"

extern "C" {
  #include "lsm_reinitialization1d.h"
  #include "lsm_reinitialization2d.h"
  #include "lsm_reinitialization3d.h"
}

// Namespaces
using namespace SAMRAI;


// Constant
#define LSM_REINIT_ALG_STOP_TOLERANCE_MAX_ITERATIONS                   (1000)

namespace LSMLIB {

/* Constructor - parameters from input database */
ReinitializationAlgorithm::ReinitializationAlgorithm(
  boost::shared_ptr<tbox::Database> input_db,
  boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const int phi_handle,
  const int control_volume_handle,
  const string& object_name)
:
d_phi_scratch_ghostcell_width(hierarchy->getDim(),3)
{
  // set object_name
  d_object_name = object_name;

  // set d_patch_hierarchy and d_grid_geometry
  d_patch_hierarchy = hierarchy;
  d_grid_geometry =
      BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
          d_patch_hierarchy->getGridGeometry());

  // set data field handles
  d_phi_handle = phi_handle;
  d_control_volume_handle = control_volume_handle;

  // get input parameters
  getFromInput(input_db);

  // check input parameters
  checkParameters();

  // create empty BoundaryConditionModule
  d_bc_module = boost::shared_ptr<BoundaryConditionModule> (
      new BoundaryConditionModule(d_patch_hierarchy,
                                  d_phi_scratch_ghostcell_width));

  // initialize variables and communication objects
  initializeVariables();
  initializeCommunicationObjects();

}


/* Constructor - parameters from arguments */
ReinitializationAlgorithm::ReinitializationAlgorithm(
  boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const int phi_handle,
  const int control_volume_handle,
  const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  const int spatial_derivative_order,
  const int tvd_runge_kutta_order,
  const LSMLIB_REAL cfl_number,
  const LSMLIB_REAL stop_distance,
  const int max_iterations,
  const LSMLIB_REAL iteration_stop_tolerance,
  const bool verbose_mode,
  const string& object_name)
:
d_phi_scratch_ghostcell_width(hierarchy->getDim(),0)
{
  // set object_name
  d_object_name = object_name;

  // set d_patch_hierarchy and d_grid_geometry
  d_patch_hierarchy = hierarchy;
  d_grid_geometry = BOOST_CAST<geom::CartesianGridGeometry,
                               hier::BaseGridGeometry>(
                                 d_patch_hierarchy->getGridGeometry());
  // set data field handles
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

    int DIM = d_patch_hierarchy->getDim().getValue();
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


  // set verbose-mode
  d_verbose_mode = verbose_mode;

  // check that the user-specifeid parameters are acceptable
  checkParameters();

  // create empty BoundaryConditionModule
  d_bc_module = boost::shared_ptr<BoundaryConditionModule> (
    new BoundaryConditionModule(d_patch_hierarchy,
                                d_phi_scratch_ghostcell_width));

  // initialize variables and communication objects
  initializeVariables();
  initializeCommunicationObjects();
}


/* reinitializeLevelSetFunctions() */
void ReinitializationAlgorithm::reinitializeLevelSetFunctions(
  const hier::IntVector& lower_bc,
  const hier::IntVector& upper_bc,
  const int max_iterations)
{

  // reset hierarchy configuration if necessary
  if (d_hierarchy_configuration_needs_reset) {
    resetHierarchyConfiguration(d_patch_hierarchy,
      0, d_patch_hierarchy->getFinestLevelNumber());
  }

  // allocate patch data for required to reinitialize level set function
  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    boost::shared_ptr<hier::PatchLevel> level =
        d_patch_hierarchy->getPatchLevel(ln);
    level->allocatePatchData(d_scratch_data);
  }

  /*
   * compute dt for reinitialization calculation
   */

  // get dx
  int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();
  boost::shared_ptr<hier::PatchLevel> patch_level =
    d_patch_hierarchy->getPatchLevel(finest_level_number);
  hier::IntVector ratio_to_coarsest = patch_level->getRatioToCoarserLevel();
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
  int num_steps = LSM_REINIT_ALG_STOP_TOLERANCE_MAX_ITERATIONS;
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
   *  compute the number of components in the level set data
   *  (if it has not already been computed)
   */
  if (d_num_phi_components == 0) {
    boost::shared_ptr<hier::PatchLevel> level =
        d_patch_hierarchy->getPatchLevel(0);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  d_object_name
                  << "::reinitializeLevelSetFunctions(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
            patch->getPatchData( d_phi_handle ));
      d_num_phi_components = phi_data->getDepth();

      break;  // only need PatchData from one patch for computation
    }
  }


  /*
   *  main reinitialization loop
   */
  int count = 0;
  LSMLIB_REAL delta = 1.0;
  const int phi_handle_after_step = d_phi_handle;
  const int phi_handle_before_step = d_phi_scr_handles[0];
  while ( (count < num_steps) &&
          (!d_use_iteration_stop_tol || (delta > d_iteration_stop_tol)) ) {

    // reset delta to zero
    delta = 0.0;

    // loop over components in level set function
    for (int component = 0; component < d_num_phi_components; component++) {

      // advance reinitialization equation using TVD Runge-Kutta
      switch(d_tvd_runge_kutta_order) {
        case 1: { // first-order TVD RK (e.g. Forward Euler)
          advanceReinitializationEqnUsingTVDRK1(
            dt, component, lower_bc, upper_bc);
          break;
        }
        case 2: { // second-order TVD RK
          advanceReinitializationEqnUsingTVDRK2(
            dt, component, lower_bc, upper_bc);
          break;
        }
        case 3: { // third-order TVD RK
          advanceReinitializationEqnUsingTVDRK3(
            dt, component, lower_bc, upper_bc);
          break;
        }
        default: { // UNSUPPORTED ORDER
          TBOX_ERROR(  d_object_name
                    << "::reinitializeLevelSetFunctions(): "
                    << "Unsupported TVD Runge-Kutta order.  "
                    << "Only TVD-RK1, TVD-RK2, and TVD-RK3 supported."
                    << endl);
        }
      } // end switch on TVD Runge-Kutta order

      // update count and delta
      if (d_use_iteration_stop_tol) {
        delta += LevelSetMethodToolbox::maxNormOfDifference(
          d_patch_hierarchy, phi_handle_after_step, phi_handle_before_step,
          d_control_volume_handle, component, 0);  // 0 is component of field
                                                   // before the time step
                                                   // which is just a single
                                                   // component scratch space
      }
    } // end loop over components of level set function

    // VERBOSE MODE
    if (d_verbose_mode) {
        tbox::pout << endl;
        tbox::pout << d_object_name << " iteration count: " << count
                   << endl;
      if (d_use_stop_distance) {
          tbox::pout << "  Level set functions reinitialized to a distance "
                     << "of approximately " << dt*count << endl;
      }
      if (d_use_iteration_stop_tol) {
          tbox::pout << "  Max norm of change in level set functions: "
                     << delta << endl;
      }
    }

    count++;

  } // end loop over evolution of reinitialization equation

  // warn if iteration terminated before stop_tol reached
  if ( d_use_iteration_stop_tol && (delta > d_iteration_stop_tol) ) {
    TBOX_WARNING(  d_object_name
                << "::reinitializeLevelSetFunctions(): "
                << "target stop tolerance ("
                << d_iteration_stop_tol << ") NOT reached after "
                << count << " time steps. "
                << "delta = " << delta
                << endl );
  }

  // VERBOSE MODE
  if (d_verbose_mode) {
      tbox::pout << endl;
      tbox::pout << "Total number of iterations: " << count << endl;
    if (d_use_stop_distance) {
        tbox::pout << "  Level set functions reinitialized to a distance "
                   << "of approximately " << dt*count << endl;
    }
    if (d_use_iteration_stop_tol)
      tbox::pout << "  Last max norm of change in level set function: "
                   << delta << endl;
  }

  // deallocate patch data that was allocated for reinitialization
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    boost::shared_ptr<hier::PatchLevel> level =
        d_patch_hierarchy->getPatchLevel(ln);
    level->deallocatePatchData(d_scratch_data);
  }
}


/* reinitializeLevelSetFunctionForSingleComponent() */
void ReinitializationAlgorithm::
  reinitializeLevelSetFunctionForSingleComponent(
    const hier::IntVector& lower_bc,
    const hier::IntVector& upper_bc,
    const int component,
    const int max_iterations)
{

  // reset hierarchy configuration if necessary
  if (d_hierarchy_configuration_needs_reset) {
    resetHierarchyConfiguration(d_patch_hierarchy,
      0, d_patch_hierarchy->getFinestLevelNumber());
  }

  // allocate patch data for required for reinitialization calculation
  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    boost::shared_ptr<hier::PatchLevel> level =
        d_patch_hierarchy->getPatchLevel(ln);
    level->allocatePatchData(d_scratch_data);
  }

  /*
   * compute dt for reinitialization calculation
   */

  // get dx
  int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();
  boost::shared_ptr<hier::PatchLevel> patch_level =
    d_patch_hierarchy->getPatchLevel(finest_level_number);
  hier::IntVector ratio_to_coarsest = patch_level->getRatioToCoarserLevel();
  int DIM = d_patch_hierarchy->getDim().getValue();
  const double* coarsest_dx = d_grid_geometry->getDx();
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
  int num_steps = LSM_REINIT_ALG_STOP_TOLERANCE_MAX_ITERATIONS;
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
   *  main reinitialization loop
   */
  int count = 0;
  LSMLIB_REAL delta = 1.0;
  const int phi_handle_after_step = d_phi_handle;
  const int phi_handle_before_step = d_phi_scr_handles[0];
  while ( (count < num_steps) &&
          (!d_use_iteration_stop_tol || (delta > d_iteration_stop_tol)) ) {

    // advance reinitialization equation using TVD Runge-Kutta
    switch(d_tvd_runge_kutta_order) {
      case 1: { // first-order TVD RK (e.g. Forward Euler)
        advanceReinitializationEqnUsingTVDRK1(
          dt, component, lower_bc, upper_bc);
        break;
      }
      case 2: { // second-order TVD RK
        advanceReinitializationEqnUsingTVDRK2(
          dt, component, lower_bc, upper_bc);
        break;
      }
      case 3: { // third-order TVD RK
        advanceReinitializationEqnUsingTVDRK3(
          dt, component, lower_bc, upper_bc);
        break;
      }
      default: { // UNSUPPORTED ORDER
        TBOX_ERROR(  d_object_name
                  << "::reinitializeLevelSetFunctionForSingleComponent(): "
                  << "Unsupported TVD Runge-Kutta order.  "
                  << "Only TVD-RK1, TVD-RK2, and TVD-RK3 supported."
                  << endl);
      }
    } // end switch on TVD Runge-Kutta order

    // update count and delta
    if (d_use_iteration_stop_tol) {
      delta = LevelSetMethodToolbox::maxNormOfDifference(
        d_patch_hierarchy, phi_handle_after_step, phi_handle_before_step,
        d_control_volume_handle, component, 0);  // 0 is component of field
                                                 // before the time step
                                                 // which is just a single
                                                 // component scratch space
    }

    // VERBOSE MODE
    if (d_verbose_mode) {
        tbox::pout << endl;
        tbox::pout << d_object_name << " iteration count: " << count
                   << endl;
      if (d_use_stop_distance) {
          tbox::pout << "  Level set functions reinitialized to a distance "
                     << "of approximately " << dt*count << endl;
      }
      if (d_use_iteration_stop_tol) {
          tbox::pout << "  Max norm of change in level set function: "
                     << delta << endl;
      }
    }

    count++;

  } // end loop over evolution of reinitialization equation

  // warn if iteration terminated before stop_tol reached
  if ( d_use_iteration_stop_tol && (delta > d_iteration_stop_tol) ) {
    TBOX_WARNING(  d_object_name
                << "::reinitializeLevelSetFunctionForSingleComponent(): "
                << "target stop tolerance ("
                << d_iteration_stop_tol << ") NOT reached after "
                << count << " time steps. "
                << "delta = " << delta
                << endl );
  }

  // VERBOSE MODE
  if (d_verbose_mode) {
      tbox::pout << endl;
      tbox::pout << "Total number of iterations: " << count << endl;
    if (d_use_stop_distance) {
        tbox::pout << "  Level set functions reinitialized to a distance "
                   << "of approximately " << dt*count << endl;
    }
    if (d_use_iteration_stop_tol)
      tbox::pout << "  Last change in max norm of level set function: "
                   << delta << endl;
  }

  // deallocate patch data that was allocated for reinitialization calculation
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    boost::shared_ptr<hier::PatchLevel> level =
      d_patch_hierarchy->getPatchLevel(ln);
    level->deallocatePatchData(d_scratch_data);
  }
}


/* resetHierarchyConfiguration() */
void ReinitializationAlgorithm::resetHierarchyConfiguration(
  boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const int coarsest_level,
  const int finest_level)
{
  // reset d_patch_hierarchy
  d_patch_hierarchy = hierarchy;

  // reset d_grid_geometry
  d_grid_geometry = BOOST_CAST<geom::CartesianGridGeometry,
                               hier::BaseGridGeometry>(
                                    d_patch_hierarchy->getGridGeometry());
  // compute RefineSchedules for filling phi boundary data
  int num_levels = hierarchy->getNumberOfLevels();
  for (int k = 0; k < d_tvd_runge_kutta_order; k++) {
    d_phi_fill_bdry_sched[k].resizeArray(num_levels);

    for (int ln = coarsest_level; ln <= finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      // reset data transfer configuration for boundary filling
      // before time advance
      d_phi_fill_bdry_sched[k][ln] =
        d_phi_fill_bdry_alg[k]->createSchedule(
          level, ln-1, hierarchy, 0);  // NULL RefinePatchStrategy

    } // end loop over levels
  } // end loop over TVD Runge-Kutta stages


  // create anti-periodic boundary condition module
  d_bc_module->resetHierarchyConfiguration(
    d_patch_hierarchy,
    coarsest_level,
    finest_level,
    d_phi_scratch_ghostcell_width);

  // set d_hierarchy_configuration_needs_reset to (finest_level < 0)
  d_hierarchy_configuration_needs_reset = (finest_level < 0);
}


/* advanceReinitializationEqnUsingTVDRK1() */
void ReinitializationAlgorithm::advanceReinitializationEqnUsingTVDRK1(
  const LSMLIB_REAL dt,
  const int phi_component,
  const hier::IntVector& lower_bc,
  const hier::IntVector& upper_bc)
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
    d_phi_scr_handles[0], d_phi_handle,
    0, phi_component);

  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical
    //       boundary conditions should be set.
    d_phi_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_bc_module->imposeBoundaryConditions(
    d_phi_scr_handles[rk_stage],
    lower_bc,
    upper_bc,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance reinitialization equation through TVD-RK1 step
  computeReinitializationEqnRHS(d_phi_scr_handles[rk_stage]);
  LevelSetMethodToolbox::TVDRK1Step(
    d_patch_hierarchy,
    d_phi_handle,
    d_phi_scr_handles[rk_stage],
    d_rhs_handle, dt,
    phi_component, 0, 0); // components of PatchData to use in TVD-RK1 step
}


/* advanceReinitializationEqnUsingTVDRK2() */
void ReinitializationAlgorithm::advanceReinitializationEqnUsingTVDRK2(
  const LSMLIB_REAL dt,
  const int phi_component,
  const hier::IntVector& lower_bc,
  const hier::IntVector& upper_bc)
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
    d_phi_scr_handles[0], d_phi_handle,
    0, phi_component);

  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical
    //       boundary conditions should be set.
    d_phi_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_bc_module->imposeBoundaryConditions(
    d_phi_scr_handles[rk_stage],
    lower_bc,
    upper_bc,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance reinitialization equation through the first stage of TVD-RK2
  computeReinitializationEqnRHS(d_phi_scr_handles[rk_stage]);
  LevelSetMethodToolbox::TVDRK2Stage1(
    d_patch_hierarchy,
    d_phi_scr_handles[rk_stage+1],
    d_phi_scr_handles[rk_stage],
    d_rhs_handle, dt,
    0, 0, 0);  // components of PatchData to use in TVD-RK2 step

  // } end Stage 1


  // { begin Stage 2

  // advance TVD RK2 stage counter
  rk_stage = 1;

  // fill scratch space for secont stage of time advance
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical
    //       boundary conditions should be set.
    d_phi_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_bc_module->imposeBoundaryConditions(
    d_phi_scr_handles[rk_stage],
    lower_bc,
    upper_bc,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance reinitialization equation through the second stage of TVD-RK2
  computeReinitializationEqnRHS(d_phi_scr_handles[rk_stage]);
  LevelSetMethodToolbox::TVDRK2Stage2(
    d_patch_hierarchy,
    d_phi_handle,
    d_phi_scr_handles[rk_stage],
    d_phi_scr_handles[0],
    d_rhs_handle, dt,
    phi_component, 0, 0, 0);  // components of PatchData to use in TVD-RK2 step

  // } end Stage 2
}


/* advanceReinitializationEqnUsingTVDRK3() */
void ReinitializationAlgorithm::advanceReinitializationEqnUsingTVDRK3(
  const LSMLIB_REAL dt,
  const int phi_component,
  const hier::IntVector& lower_bc,
  const hier::IntVector& upper_bc)
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
    d_phi_scr_handles[0], d_phi_handle,
    0, phi_component);

  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical
    //       boundary conditions should be set.
    d_phi_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_bc_module->imposeBoundaryConditions(
    d_phi_scr_handles[rk_stage],
    lower_bc,
    upper_bc,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance reinitialization equation through the first stage of TVD-RK3
  computeReinitializationEqnRHS(d_phi_scr_handles[rk_stage]);
  LevelSetMethodToolbox::TVDRK3Stage1(
    d_patch_hierarchy,
    d_phi_scr_handles[rk_stage+1],
    d_phi_scr_handles[rk_stage],
    d_rhs_handle, dt,
    0, 0, 0);  // components of PatchData to use in TVD-RK3 step

  // } end Stage 1


  // { begin Stage 2

  // advance TVD RK3 stage counter
  rk_stage = 1;

  // fill scratch space for secont stage of time advance
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical
    //       boundary conditions should be set.
    d_phi_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_bc_module->imposeBoundaryConditions(
    d_phi_scr_handles[rk_stage],
    lower_bc,
    upper_bc,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance reinitialization equation through the second stage of TVD-RK3
  computeReinitializationEqnRHS(d_phi_scr_handles[rk_stage]);
  LevelSetMethodToolbox::TVDRK3Stage2(
    d_patch_hierarchy,
    d_phi_scr_handles[rk_stage+1],
    d_phi_scr_handles[rk_stage],
    d_phi_scr_handles[rk_stage-1],
    d_rhs_handle, dt,
    0, 0, 0, 0);  // components of PatchData to use in TVD-RK3 step

  // } end Stage 2


  // { begin Stage 3

  // advance TVD RK3 stage counter
  rk_stage = 2;

  // fill scratch space for secont stage of time advance
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: 0.0 is "current time" and true indicates that physical
    //       boundary conditions should be set.
    d_phi_fill_bdry_sched[rk_stage][ln]->fillData(0.0,true);
  }
  d_bc_module->imposeBoundaryConditions(
    d_phi_scr_handles[rk_stage],
    lower_bc,
    upper_bc,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    0);

  // advance reinitialization equation through the third stage of TVD-RK3
  computeReinitializationEqnRHS(d_phi_scr_handles[rk_stage]);
  LevelSetMethodToolbox::TVDRK3Stage3(
    d_patch_hierarchy,
    d_phi_handle,
    d_phi_scr_handles[rk_stage],
    d_phi_scr_handles[0],
    d_rhs_handle, dt,
    phi_component, 0, 0, 0);  // components of PatchData to use in TVD-RK3 step

  // } end Stage 3
}


/* computeReinitializationEqnRHS() */
void ReinitializationAlgorithm::computeReinitializationEqnRHS(
  const int phi_handle)
{

  // compute spatial derivatives for the current stage
  LevelSetMethodToolbox::computePlusAndMinusSpatialDerivatives(
    d_patch_hierarchy,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    d_grad_phi_plus_handle,
    d_grad_phi_minus_handle,
    phi_handle);

  // loop over PatchHierarchy and compute RHS for level set equation
  // by calling Fortran routines
  const int num_levels = d_patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

      boost::shared_ptr<hier::PatchLevel> level =
          d_patch_hierarchy->getPatchLevel(0);
    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  d_object_name
                  << "::computeReinitializationEqnRHS(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> rhs_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( d_rhs_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_plus_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( d_grad_phi_plus_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_minus_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( d_grad_phi_minus_handle ));

      hier::Box rhs_ghostbox = rhs_data->getGhostBox();
      const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
      const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

      hier::Box phi_ghostbox = phi_data->getGhostBox();
      const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
      const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

      hier::Box grad_phi_plus_ghostbox = grad_phi_plus_data->getGhostBox();
      const hier::IntVector grad_phi_plus_ghostbox_lower =
        grad_phi_plus_ghostbox.lower();
      const hier::IntVector grad_phi_plus_ghostbox_upper =
        grad_phi_plus_ghostbox.upper();

      hier::Box grad_phi_minus_ghostbox = grad_phi_minus_data->getGhostBox();
      const hier::IntVector grad_phi_minus_ghostbox_lower =
        grad_phi_minus_ghostbox.lower();
      const hier::IntVector grad_phi_minus_ghostbox_upper =
        grad_phi_minus_ghostbox.upper();

      // fill box
      hier::Box fillbox = rhs_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      LSMLIB_REAL* rhs = rhs_data->getPointer();
      LSMLIB_REAL* phi = phi_data->getPointer();
      LSMLIB_REAL* grad_phi_plus[LSM_DIM_MAX];
      LSMLIB_REAL* grad_phi_minus[LSM_DIM_MAX];

     int DIM = d_patch_hierarchy->getDim().getValue();
       for (int dim = 0; dim < DIM; dim++) {
        grad_phi_plus[dim] = grad_phi_plus_data->getPointer(dim);
        grad_phi_minus[dim] = grad_phi_minus_data->getPointer(dim);
      }

      // get dx
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST <geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry());
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geom->getDx();
#else
      const double* dx_double = patch_geom->getDx();
      float dx[DIM];
      for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif

      // flag for whether or not to use phi0 in computing sgn(phi)
      int use_phi0 = 0; // KTC do NOT use phi0 for sgn(phi) calculation

      if (DIM == 3) {

        LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          &phi_ghostbox_lower[2],
          &phi_ghostbox_upper[2],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          &phi_ghostbox_lower[2],
          &phi_ghostbox_upper[2],
          grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
          &grad_phi_plus_ghostbox_lower[0],
          &grad_phi_plus_ghostbox_upper[0],
          &grad_phi_plus_ghostbox_lower[1],
          &grad_phi_plus_ghostbox_upper[1],
          &grad_phi_plus_ghostbox_lower[2],
          &grad_phi_plus_ghostbox_upper[2],
          grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
          &grad_phi_minus_ghostbox_lower[0],
          &grad_phi_minus_ghostbox_upper[0],
          &grad_phi_minus_ghostbox_lower[1],
          &grad_phi_minus_ghostbox_upper[1],
          &grad_phi_minus_ghostbox_lower[2],
          &grad_phi_minus_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2],
          &dx[0], &dx[1], &dx[2],
          &use_phi0);

      } else if (DIM == 2) {

        LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          grad_phi_plus[0], grad_phi_plus[1],
          &grad_phi_plus_ghostbox_lower[0],
          &grad_phi_plus_ghostbox_upper[0],
          &grad_phi_plus_ghostbox_lower[1],
          &grad_phi_plus_ghostbox_upper[1],
          grad_phi_minus[0], grad_phi_minus[1],
          &grad_phi_minus_ghostbox_lower[0],
          &grad_phi_minus_ghostbox_upper[0],
          &grad_phi_minus_ghostbox_lower[1],
          &grad_phi_minus_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &dx[0], &dx[1],
          &use_phi0);

      } else if (DIM == 1) {

        LSM1D_COMPUTE_REINITIALIZATION_EQN_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          grad_phi_plus[0],
          &grad_phi_plus_ghostbox_lower[0],
          &grad_phi_plus_ghostbox_upper[0],
          grad_phi_minus[0],
          &grad_phi_minus_ghostbox_lower[0],
          &grad_phi_minus_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &dx[0],
          &use_phi0);

      } else {  // Unsupported dimension
        TBOX_ERROR(  d_object_name
                  << "::computeReinitializationEqnRHS(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

    } // end loop over patches in level
  } // end loop over levels in hierarchy

}


/* initializeVariables() */
void ReinitializationAlgorithm::initializeVariables()
{
  // initialize d_num_phi_components to zero
  d_num_phi_components = 0;

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

  d_phi_scratch_ghostcell_width =
    hier::IntVector(d_patch_hierarchy->getDim(),
                    scratch_ghostcell_width_for_grad);
  hier::IntVector zero_ghostcell_width(d_patch_hierarchy->getDim(),0);

  /*
   * create variables and PatchData for scratch data
   */
  hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

  // get variable associated with phi_handle
  boost::shared_ptr<hier::Variable> tmp_variable;
  boost::shared_ptr<hier::VariableContext> tmp_context;
  boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>> phi_variable;
  if (var_db->mapIndexToVariableAndContext(d_phi_handle,
                                           tmp_variable, tmp_context)) {
    phi_variable = BOOST_CAST<pdat::CellVariable<LSMLIB_REAL>, hier::Variable>
     (tmp_variable);
  } else {
    TBOX_ERROR(  d_object_name
              << "::initializeVariables(): "
              << "Specified phi handle does not exist or does "
              << "not correspond to cell-centered data"
              << endl);
  }

  // create scratch context for reinitialization
  boost::shared_ptr<hier::VariableContext> scratch_context =
    var_db->getContext("REINITIALIZATION_SCRATCH");

  // reserve space for scratch PatchData handles
  d_phi_scr_handles.reserve(d_tvd_runge_kutta_order);

  // clear out ComponentSelector for scratch data
  d_scratch_data.clrAllFlags();

  // create "SCRATCH" context for phi
  stringstream phi_scratch_variable_name("");
  phi_scratch_variable_name << phi_variable->getName()
                            << "::REINITIALIZATION_PHI_SCRATCH";
  boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>> phi_scratch_variable;
  if (var_db->checkVariableExists(phi_scratch_variable_name.str())) {
   phi_scratch_variable = BOOST_CAST<pdat::CellVariable<LSMLIB_REAL>,
                                     hier::Variable >
     (var_db->getVariable(phi_scratch_variable_name.str()));
  } else {
   phi_scratch_variable = boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>>(
     new pdat::CellVariable<LSMLIB_REAL>(d_patch_hierarchy->getDim(),
                                         phi_scratch_variable_name.str(), 1));
  }
  for (int k=0; k < d_tvd_runge_kutta_order; k++) {
    stringstream context_name("");
    context_name << "REINITIALIZATION_SCRATCH_"
                 << d_phi_handle
                 << "::" << k;
    d_phi_scr_handles[k] =
      var_db->registerVariableAndContext(
        phi_scratch_variable,
        var_db->getContext(context_name.str()),
        d_phi_scratch_ghostcell_width);
    d_scratch_data.setFlag(d_phi_scr_handles[k]);
  }

  // create RHS variables
  stringstream rhs_name("");
  rhs_name << phi_variable->getName() << "::REINITIALIZATION_RHS";
  boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>> rhs_variable;
  if (var_db->checkVariableExists(rhs_name.str())) {
   rhs_variable = BOOST_CAST<pdat::CellVariable<LSMLIB_REAL>, hier::Variable>
     (var_db->getVariable(rhs_name.str()));
  } else {
   rhs_variable = boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>>
    (new pdat::CellVariable<LSMLIB_REAL>(d_patch_hierarchy->getDim(),
                                         rhs_name.str(), 1));
  }
  d_rhs_handle = var_db->registerVariableAndContext(
    rhs_variable, scratch_context, zero_ghostcell_width);
  d_scratch_data.setFlag(d_rhs_handle);

  // create variables for grad(phi)
  stringstream grad_phi_name("");
  grad_phi_name << phi_variable->getName()
                << "::REINITIALIZATION_GRAD_PHI";
  boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>> grad_phi_variable;
  if (var_db->checkVariableExists(grad_phi_name.str())) {
   grad_phi_variable = BOOST_CAST<pdat::CellVariable<LSMLIB_REAL>,
                                  hier::Variable>
     (var_db->getVariable(grad_phi_name.str()));
  } else {
   grad_phi_variable = boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>>(
    new pdat::CellVariable<LSMLIB_REAL>(d_patch_hierarchy->getDim(),
                                        grad_phi_name.str()));
  }
  boost::shared_ptr<hier::VariableContext> grad_phi_plus_context =
    var_db->getContext("REINITIALIZATION_GRAD_PHI_PLUS");
  boost::shared_ptr<hier::VariableContext> grad_phi_minus_context =
    var_db->getContext("REINITIALIZATION_GRAD_PHI_MINUS");
  d_grad_phi_plus_handle = var_db->registerVariableAndContext(
    grad_phi_variable, grad_phi_plus_context, zero_ghostcell_width);
  d_grad_phi_minus_handle = var_db->registerVariableAndContext(
    grad_phi_variable, grad_phi_minus_context, zero_ghostcell_width);
  d_scratch_data.setFlag(d_grad_phi_plus_handle);
  d_scratch_data.setFlag(d_grad_phi_minus_handle);

}


/* initializeCommunicationObjects() */
void ReinitializationAlgorithm::initializeCommunicationObjects()
{

  // initialize d_hierarchy_configuration_needs_reset to true
  d_hierarchy_configuration_needs_reset = true;

  // get pointer to VariableDatabase
  hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

  /*
   * Lookup refine operations
   */
  // get CellVariable associated with d_phi_handle
  boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>> phi_variable;
  boost::shared_ptr<hier::Variable> tmp_variable;
  boost::shared_ptr<hier::VariableContext> tmp_context;
  if (var_db->mapIndexToVariableAndContext(d_phi_handle,
                                           tmp_variable, tmp_context)) {
    phi_variable = BOOST_CAST<pdat::CellVariable<LSMLIB_REAL>, hier::Variable>
     (tmp_variable);
  } else {
    TBOX_ERROR(  d_object_name
              << "::initializeCommunicationObjects(): "
              << "Specified phi handle does not exist or does "
              << "not correspond to cell-centered data"
              << endl);
  }
  // lookup refine operations
  boost::shared_ptr<hier::RefineOperator> refine_op =
    d_grid_geometry->lookupRefineOperator(phi_variable, "LINEAR_REFINE");


  /*
   * create RefineAlgorithms for filling boundary data for phi
   */
  d_phi_fill_bdry_alg.resizeArray(d_tvd_runge_kutta_order);
  d_phi_fill_bdry_sched.resizeArray(d_tvd_runge_kutta_order);
  for (int k = 0; k < d_tvd_runge_kutta_order; k++) {
    d_phi_fill_bdry_alg[k] = boost::shared_ptr<xfer::RefineAlgorithm>
     (new xfer::RefineAlgorithm);

    // empty out the boundary bdry fill schedules
    d_phi_fill_bdry_sched[k].setNull();

      // register data transfer to fill boundary data before time advance
      d_phi_fill_bdry_alg[k]->registerRefine(
        d_phi_scr_handles[k],
        d_phi_scr_handles[k],
        d_phi_scr_handles[k],
        refine_op);
  } // end loop over TVD-Runge-Kutta stages
// configure communications schedules for ALL levels using the
  // specified PatchHierarchy
    resetHierarchyConfiguration(d_patch_hierarchy, 0, d_patch_hierarchy->getFinestLevelNumber());
}


/* getFromInput() */
void ReinitializationAlgorithm::getFromInput(
  boost::shared_ptr<tbox::Database> db)
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

    int DIM = d_patch_hierarchy->getDim().getValue();
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
void ReinitializationAlgorithm::checkParameters()
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
