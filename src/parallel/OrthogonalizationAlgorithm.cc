/*
 * File:        OrthogonalizationAlgorithm.cc
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.27 $
 * Modified:    $Date: 2006/10/05 15:03:45 $
 * Description: Implementation file for level set method orthogonalization 
 *              class
 */

#ifndef included_OrthogonalizationAlgorithm_cc
#define included_OrthogonalizationAlgorithm_cc

// System Headers
#include <sstream>
extern "C" {
  #include <limits.h>
}

#include "LSMLIB_DefaultParameters.h"
#include "OrthogonalizationAlgorithm.h" 

// SAMRAI Headers
#include "CartesianPatchGeometry.h" 
#include "CellData.h" 
#include "IntVector.h" 
#include "Patch.h" 
#include "PatchLevel.h" 

// Headers for Fortran kernels
extern "C" {
  #include "lsm_reinitialization2d.h"
  #include "lsm_reinitialization3d.h"
  #include "lsm_samrai_f77_utilities.h"
}

// SAMRAI namespaces
using namespace std;
using namespace pdat;


namespace LSMLIB {

/* Constructor - parameters from input database */
template <int DIM>
OrthogonalizationAlgorithm<DIM>::OrthogonalizationAlgorithm(
  Pointer<Database> input_db,
  Pointer< PatchHierarchy<DIM> > hierarchy,
  const int phi_handle,
  const int psi_handle,
  const int control_volume_handle,
  const string& object_name,
  const IntVector<DIM>& phi_ghostcell_width,
  const IntVector<DIM>& psi_ghostcell_width)
{
  // set object_name
  d_object_name = object_name;

  // set d_patch_hierarchy and d_grid_geometry
  d_patch_hierarchy = hierarchy;
  d_grid_geometry = hierarchy->getGridGeometry();

  // set data field handles
  d_phi_handle = phi_handle;
  d_psi_handle = psi_handle;
  d_control_volume_handle = control_volume_handle;

  // get input parameters
  getFromInput(input_db);

  // check input parameters
  checkParameters();

  // initialize d_num_field_components to zero
  d_num_field_components = 0;

  // create FieldExtensionAlgorithm objects
  d_fixed_phi_field_ext_alg = new FieldExtensionAlgorithm<DIM>(
    d_patch_hierarchy,
    psi_handle,
    phi_handle,
    control_volume_handle,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    d_tvd_runge_kutta_order,
    d_cfl_number,
    d_stop_distance,
    d_max_iterations,
    d_iteration_stop_tol,
    d_verbose_mode,
    object_name + "::FIXED_PHI",
    phi_ghostcell_width);

  d_fixed_psi_field_ext_alg = new FieldExtensionAlgorithm<DIM>(
    d_patch_hierarchy,
    phi_handle,
    psi_handle,
    control_volume_handle,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    d_tvd_runge_kutta_order,
    d_cfl_number,
    d_stop_distance,
    d_max_iterations,
    d_iteration_stop_tol,
    d_verbose_mode,
    object_name + "::FIXED_PSI",
    psi_ghostcell_width);

}


/* Constructor - parameters from arguments */
template <int DIM>
OrthogonalizationAlgorithm<DIM>::OrthogonalizationAlgorithm(
  Pointer< PatchHierarchy<DIM> > hierarchy,
  const int phi_handle,
  const int psi_handle,
  const int control_volume_handle,
  const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  const int spatial_derivative_order,
  const int tvd_runge_kutta_order,
  const double cfl_number,
  const double stop_distance, 
  const int max_iterations,
  const double iteration_stop_tolerance,
  const bool verbose_mode,
  const string& object_name,
  const IntVector<DIM>& phi_ghostcell_width,
  const IntVector<DIM>& psi_ghostcell_width)
{
  // set object_name
  d_object_name = object_name;

  // set d_patch_hierarchy and d_grid_geometry
  d_patch_hierarchy = hierarchy;
  d_grid_geometry = hierarchy->getGridGeometry();

  // set data field handles
  d_phi_handle = phi_handle;
  d_psi_handle = psi_handle;
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
    
    const double *X_lower = d_grid_geometry->getXLower();
    const double *X_upper = d_grid_geometry->getXUpper();
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

  // initialize d_num_field_components to zero
  d_num_field_components = 0;

  // create FieldExtensionAlgorithm objects
  d_fixed_phi_field_ext_alg = new FieldExtensionAlgorithm<DIM>(
    d_patch_hierarchy,
    psi_handle,
    phi_handle,
    control_volume_handle,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    d_tvd_runge_kutta_order,
    d_cfl_number,
    d_stop_distance,
    d_max_iterations,
    d_iteration_stop_tol,
    d_verbose_mode,
    object_name + "::FIXED_PHI",
    phi_ghostcell_width);

  d_fixed_psi_field_ext_alg = new FieldExtensionAlgorithm<DIM>(
    d_patch_hierarchy,
    phi_handle,
    psi_handle,
    control_volume_handle,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    d_tvd_runge_kutta_order,
    d_cfl_number,
    d_stop_distance,
    d_max_iterations,
    d_iteration_stop_tol,
    d_verbose_mode,
    object_name + "::FIXED_PSI",
    psi_ghostcell_width);

}


/* orthogonalizeLevelSetFunctions() */
template <int DIM> 
void OrthogonalizationAlgorithm<DIM>::orthogonalizeLevelSetFunctions(
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const int max_iterations,
  const IntVector<DIM>& lower_bc_fixed,
  const IntVector<DIM>& upper_bc_fixed,
  const IntVector<DIM>& lower_bc_evolved,
  const IntVector<DIM>& upper_bc_evolved)
{
  // compute the number of components of level set function
  if (d_num_field_components == 0) {
    Pointer< PatchLevel<DIM> > level = d_patch_hierarchy->getPatchLevel(0);
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  d_object_name
                  << "::orthogonalizeLevelSetFunctions(): " 
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }
  
      Pointer< CellData<DIM,double> > phi_data = 
        patch->getPatchData( d_phi_handle );
      d_num_field_components = phi_data->getDepth();
  
      break;  // only need PatchData from one patch for computation
    }
  }

  for (int component = 0; component < d_num_field_components; component++) {
    if (level_set_fcn == PHI) {
      d_fixed_psi_field_ext_alg->computeExtensionFieldForSingleComponent(
        component, component, 
        max_iterations,
        lower_bc_fixed, upper_bc_fixed,
        lower_bc_evolved, upper_bc_evolved);
    } else {
      d_fixed_phi_field_ext_alg->computeExtensionFieldForSingleComponent(
        component, component,
        max_iterations,
        lower_bc_fixed, upper_bc_fixed,
        lower_bc_evolved, upper_bc_evolved);
    }
  }

}


/* orthogonalizeLevelSetFunctionForSingleComponent() */
template <int DIM> 
void OrthogonalizationAlgorithm<DIM>::
  orthogonalizeLevelSetFunctionForSingleComponent(
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const int component,
    const int max_iterations,
    const IntVector<DIM>& lower_bc_fixed,
    const IntVector<DIM>& upper_bc_fixed,
    const IntVector<DIM>& lower_bc_evolved,
    const IntVector<DIM>& upper_bc_evolved)
{

  if (level_set_fcn == PHI) {
    d_fixed_psi_field_ext_alg->computeExtensionFieldForSingleComponent(
      component, component,
      max_iterations,
      lower_bc_fixed, upper_bc_fixed,
      lower_bc_evolved, upper_bc_evolved);
  } else {
    d_fixed_phi_field_ext_alg->computeExtensionFieldForSingleComponent(
      component, component,
      max_iterations,
      lower_bc_fixed, upper_bc_fixed,
      lower_bc_evolved, upper_bc_evolved);
  }

}


/* resetHierarchyConfiguration() */
template <int DIM>
void OrthogonalizationAlgorithm<DIM>::resetHierarchyConfiguration(
  Pointer< PatchHierarchy<DIM> > hierarchy,
  const int coarsest_level,
  const int finest_level)
{
    d_fixed_psi_field_ext_alg->resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level);

    d_fixed_phi_field_ext_alg->resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level);

}


/* getFromInput() */
template <int DIM> 
void OrthogonalizationAlgorithm<DIM>::getFromInput(
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
    
    const double *X_lower = d_grid_geometry->getXLower();
    const double *X_upper = d_grid_geometry->getXUpper();
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
void OrthogonalizationAlgorithm<DIM>::checkParameters()
{
  // check that (DIM >= 2)
  if (DIM < 2) {
    TBOX_ERROR(  d_object_name
              << "::checkParameters(): "
              << "Invalid value of DIM.  "
              << "Orthogonalization procedure only applies for DIM > 1."
              << endl );
  }

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
