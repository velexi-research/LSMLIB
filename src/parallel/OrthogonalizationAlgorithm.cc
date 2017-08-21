/*
 * File:        OrthogonalizationAlgorithm.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Implementation file for level set method orthogonalization
 *              class
 */

#ifndef included_OrthogonalizationAlgorithm_cc
#define included_OrthogonalizationAlgorithm_cc

// Class header file
#include "OrthogonalizationAlgorithm.h"

// Standard library headers
#include <cstddef>
#include <ostream>
#include <string>

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI Headers
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

// LSMLIB headers
#include "LSMLIB_config.h"
#include "LSMLIB_DefaultParameters.h"
#include "FieldExtensionAlgorithm.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class BaseGridGeometry; } }
namespace SAMRAI { namespace hier { class PatchData; } }

// Namespaces
using namespace SAMRAI;

namespace LSMLIB {

/* Constructor - parameters from input database */
OrthogonalizationAlgorithm::OrthogonalizationAlgorithm(
  boost::shared_ptr<tbox::Database> input_db,
  boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const int phi_handle,
  const int psi_handle,
  const int control_volume_handle,
  const hier::IntVector& phi_ghostcell_width,
  const hier::IntVector& psi_ghostcell_width,
  const string& object_name)
{
  // set object_name
  d_object_name = object_name;

  // set d_patch_hierarchy and d_grid_geometry
  d_patch_hierarchy = hierarchy;
  d_grid_geometry =
      BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
        hierarchy->getGridGeometry());

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
  boost::shared_ptr<FieldExtensionAlgorithm> d_fixed_phi_field_ext_alg =
    boost::shared_ptr<FieldExtensionAlgorithm>(new FieldExtensionAlgorithm(
      d_patch_hierarchy,
      psi_handle,
      phi_handle,
      control_volume_handle,
      phi_ghostcell_width,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      d_tvd_runge_kutta_order,
      d_cfl_number,
      d_stop_distance,
      d_max_iterations,
      d_iteration_stop_tol,
      d_verbose_mode,
      object_name + "::FIXED_PHI"));

  boost::shared_ptr<FieldExtensionAlgorithm> d_fixed_psi_field_ext_alg =
    boost::shared_ptr<FieldExtensionAlgorithm>(new FieldExtensionAlgorithm(
      d_patch_hierarchy,
      phi_handle,
      psi_handle,
      control_volume_handle,
      phi_ghostcell_width,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      d_tvd_runge_kutta_order,
      d_cfl_number,
      d_stop_distance,
      d_max_iterations,
      d_iteration_stop_tol,
      d_verbose_mode,
      object_name + "::FIXED_PSI"));

}


/* Constructor - parameters from arguments */
OrthogonalizationAlgorithm::OrthogonalizationAlgorithm(
  boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const int phi_handle,
  const int psi_handle,
  const int control_volume_handle,
  const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  const int spatial_derivative_order,
  const int tvd_runge_kutta_order,
  const LSMLIB_REAL cfl_number,
  const hier::IntVector& phi_ghostcell_width,
  const hier::IntVector& psi_ghostcell_width,
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
  d_grid_geometry =
      BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
        hierarchy->getGridGeometry());

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
      int DIM = d_patch_hierarchy->getDim().getValue();
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
  boost::shared_ptr<FieldExtensionAlgorithm> d_fixed_phi_field_ext_alg =
    boost::shared_ptr<FieldExtensionAlgorithm>(new FieldExtensionAlgorithm(
      d_patch_hierarchy,
      psi_handle,
      phi_handle,
      control_volume_handle,
      phi_ghostcell_width,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      d_tvd_runge_kutta_order,
      d_cfl_number,
      d_stop_distance,
      d_max_iterations,
      d_iteration_stop_tol,
      d_verbose_mode,
      object_name + "::FIXED_PHI"));

  boost::shared_ptr<FieldExtensionAlgorithm> d_fixed_psi_field_ext_alg =
    boost::shared_ptr<FieldExtensionAlgorithm>(new FieldExtensionAlgorithm(
      d_patch_hierarchy,
      phi_handle,
      psi_handle,
      control_volume_handle,
      phi_ghostcell_width,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      d_tvd_runge_kutta_order,
      d_cfl_number,
      d_stop_distance,
      d_max_iterations,
      d_iteration_stop_tol,
      d_verbose_mode,
      object_name + "::FIXED_PSI"));

}


/* orthogonalizeLevelSetFunctions() */
void OrthogonalizationAlgorithm::orthogonalizeLevelSetFunctions(
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const hier::IntVector& lower_bc_fixed,
  const hier::IntVector& upper_bc_fixed,
  const hier::IntVector& lower_bc_evolved,
  const hier::IntVector& upper_bc_evolved,
  const int max_iterations)
{
  // compute the number of components of level set function
  if (d_num_field_components == 0) {
    boost::shared_ptr<hier::PatchLevel> level =
        d_patch_hierarchy->getPatchLevel(0);
    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  d_object_name
                  << "::orthogonalizeLevelSetFunctions(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
            patch->getPatchData(d_phi_handle));
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
void OrthogonalizationAlgorithm::
  orthogonalizeLevelSetFunctionForSingleComponent(
    const LEVEL_SET_FCN_TYPE level_set_fcn,
    const hier::IntVector& lower_bc_fixed,
    const hier::IntVector& upper_bc_fixed,
    const hier::IntVector& lower_bc_evolved,
    const hier::IntVector& upper_bc_evolved,
    const int component,
    const int max_iterations)
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
void OrthogonalizationAlgorithm::resetHierarchyConfiguration(
  boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const int coarsest_level,
  const int finest_level)
{
    d_fixed_psi_field_ext_alg->resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level);

    d_fixed_phi_field_ext_alg->resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level);

}


/* getFromInput() */
void OrthogonalizationAlgorithm::getFromInput(
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
void OrthogonalizationAlgorithm::checkParameters()
{
   int DIM = d_patch_hierarchy->getDim().getValue();
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
