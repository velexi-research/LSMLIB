/*
 * File:        LevelSetMethodAlgorithm.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Implementation file for level set method integrator class
 */

#ifndef included_LevelSetMethodAlgorithm_cc
#define included_LevelSetMethodAlgorithm_cc

#include "LevelSetMethodAlgorithm.h"

// Standard library headers
#include <cstddef>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Utilities.h"

// LSMLIB headers
#include "LSMLIB_config.h"
#include "LSMLIB_DefaultParameters.h"
#include "FieldExtensionAlgorithm.h"
#include "LevelSetMethodGriddingAlgorithm.h"
#include "LevelSetFunctionIntegrator.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <cassert>
#endif
#endif

/****************************************************************
 *
 * Non-inline methods for LevelSetMethodAlgorithm class
 *
 ****************************************************************/

namespace LSMLIB {

/* Constructor for standard integrator and gridding strategy */
LevelSetMethodAlgorithm::LevelSetMethodAlgorithm(
  const tbox::Dimension& dim,
  boost::shared_ptr<Database> lsm_algorithm_input_db,
  boost::shared_ptr< PatchHierarchy > patch_hierarchy,
  LevelSetMethodPatchStrategy* patch_strategy,
  LevelSetMethodVelocityFieldStrategy* velocity_field_strategy,
  const int num_level_set_fcn_components,
  const int codimension,
  const string& object_name)
{

#ifdef DEBUG_CHECK_ASSERTIONS
  assert(lsm_algorithm_input_db != NULL);
  assert(patch_hierarchy != NULL);
  assert(patch_strategy != NULL);
  assert(velocity_field_strategy != NULL);
  assert(!object_name.empty());
#endif
  d_patch_hierarchy = patch_hierarchy;
  d_object_name = object_name;

  // create new LevelSetFunctionIntegrator object
  boost::shared_ptr<Database> level_set_fcn_integrator_db =
    lsm_algorithm_input_db->getDatabase("LevelSetFunctionIntegrator");
  d_lsm_integrator_strategy =
    boost::shared_ptr<LevelSetFunctionIntegrator> (
      new LevelSetFunctionIntegrator(
        dim,
        level_set_fcn_integrator_db,
        patch_hierarchy,
        patch_strategy,
        velocity_field_strategy,
        num_level_set_fcn_components,
        codimension,
        object_name + "::LevelSetFunctionIntegrator"));

  // cache simulation parameters for standard LevelSetFunctionIntegrator
  d_using_standard_level_set_fcn_integrator = true;
  string spatial_derivative_type =
    level_set_fcn_integrator_db->getStringWithDefault(
      "spatial_derivative_type", LSM_DEFAULT_SPATIAL_DERIVATIVE_TYPE);
  if (spatial_derivative_type == "WENO") {
    d_spatial_derivative_type = WENO;
    d_spatial_derivative_order =
      level_set_fcn_integrator_db->getIntegerWithDefault(
        "spatial_derivative_order", LSM_DEFAULT_SPATIAL_DERIVATIVE_WENO_ORDER);
  } else if (spatial_derivative_type == "ENO") {
    d_spatial_derivative_type = ENO;
    d_spatial_derivative_order =
      level_set_fcn_integrator_db->getIntegerWithDefault(
        "spatial_derivative_order", LSM_DEFAULT_SPATIAL_DERIVATIVE_ENO_ORDER);
  } else {
    d_spatial_derivative_type = UNKNOWN;
  }
  d_tvd_runge_kutta_order =
    level_set_fcn_integrator_db->getIntegerWithDefault(
      "tvd_runge_kutta_order", LSM_DEFAULT_TVD_RUNGE_KUTTA_ORDER);
  d_cfl_number =
    level_set_fcn_integrator_db->getDoubleWithDefault("cfl_number",
                                                       LSM_DEFAULT_CFL_NUMBER);

  // create new LevelSetMethodGriddingAlgorithm object
  d_lsm_gridding_strategy = boost::shared_ptr<LevelSetMethodGriddingAlgorithm> (new LevelSetMethodGriddingAlgorithm(
    lsm_algorithm_input_db->getDatabase("LevelSetMethodGriddingAlgorithm"),
    patch_hierarchy,
    d_lsm_integrator_strategy,
    object_name + "::LevelSetMethodGriddingAlgorithm"));

  // register velocity_field_strategy with d_lsm_gridding_strategy
  d_lsm_gridding_strategy->registerVelocityFieldStrategy(
    velocity_field_strategy);

}
/* Constructor for custom integrator and gridding strategy */
LevelSetMethodAlgorithm::LevelSetMethodAlgorithm(
  boost::shared_ptr< LevelSetFunctionIntegratorStrategy > lsm_integrator_strategy,
  boost::shared_ptr< LevelSetMethodGriddingStrategy > lsm_gridding_strategy,
  const string& object_name)
{

#ifdef DEBUG_CHECK_ASSERTIONS
  assert(lsm_integrator_strategy);
  assert(lsm_gridding_strategy);
  assert(!object_name.empty());
#endif

  d_lsm_integrator_strategy = lsm_integrator_strategy;
  d_lsm_gridding_strategy = lsm_gridding_strategy;
  d_object_name = object_name;

  // set simulation parameters used by standard LevelSetFunctionIntegrator
  // class to invalid values
  d_using_standard_level_set_fcn_integrator = false;
  d_patch_hierarchy = boost::shared_ptr<PatchHierarchy>();
  d_spatial_derivative_type = UNKNOWN;
  d_spatial_derivative_order = 0;
  d_tvd_runge_kutta_order = 0;
  d_cfl_number = 0.0;

}


/* printClassData() */
void LevelSetMethodAlgorithm::printClassData(ostream& os) const
{
  os << "\n===================================" << endl;
  os << "LevelSetMethodAlgorithm" << endl;

  os << "Object Pointers" << endl;
  os << "---------------" << endl;
  os << "(LevelSetMethodAlgorithm*) this = "
     << (LevelSetMethodAlgorithm*) this << endl;
  os << "d_lsm_integrator_strategy = "
     << d_lsm_integrator_strategy.get() << endl;
  os << "d_lsm_gridding_strategy = "
     << d_lsm_gridding_strategy.get() << endl;
  os << "===================================" << endl << endl;
}


/* resetHierarchyConfiguration() */
void LevelSetMethodAlgorithm::resetHierarchyConfiguration(
  const boost::shared_ptr< PatchHierarchy > hierarchy,
  const int coarsest_level,
  const int finest_level)
{
  d_lsm_gridding_strategy->resetHierarchyConfiguration(
    hierarchy, coarsest_level, finest_level);

  // invoke resetHierarchyConfiguration() for FieldExtensionAlgorithm
  // objects created by this LevelSetMethodAlgorithm object
  int num_field_extension_algs = d_field_extension_alg_list.size();
  for (int k = 0; k < num_field_extension_algs; k++) {
    d_field_extension_alg_list[k]->resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level);
  }

}


/* getFieldExtensionAlgorithm() */
boost::shared_ptr< FieldExtensionAlgorithm >
LevelSetMethodAlgorithm::getFieldExtensionAlgorithm(
  boost::shared_ptr<Database> input_db,
  const int field_handle,
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const string& object_name)
{
  if (!d_using_standard_level_set_fcn_integrator) {
    TBOX_ERROR( d_object_name
              << "::getFieldExtensionAlgorithm(): "
              << "this method may only be used when using the standard "
              << "LevelSetFunctionIntegrator."
              << endl );

    boost::shared_ptr< FieldExtensionAlgorithm > field_extension_alg=boost::shared_ptr< FieldExtensionAlgorithm > ();
    return field_extension_alg;
  }

  int level_set_fcn_handle;
  if (level_set_fcn == PHI) {
    level_set_fcn_handle = d_lsm_integrator_strategy->getPhiPatchDataHandle();
  } else {
    level_set_fcn_handle = d_lsm_integrator_strategy->getPsiPatchDataHandle();
  }

  boost::shared_ptr< FieldExtensionAlgorithm > field_extension_alg =
   boost::shared_ptr< FieldExtensionAlgorithm > ( new FieldExtensionAlgorithm(
      input_db,
      d_patch_hierarchy,
      field_handle,
      level_set_fcn_handle,
      getControlVolumePatchDataHandle(),
      IntVector(d_patch_hierarchy->getDim(), 0),
      object_name));

  // add new object to list of FieldExtensionAlgorithms to be
  // reset by resetHierarchyConfiguration
  int length = d_field_extension_alg_list.size();
  d_field_extension_alg_list.resizeArray(length+1);
  d_field_extension_alg_list[length] = field_extension_alg;

  return field_extension_alg;
}


/* getFieldExtensionAlgorithm() */
boost::shared_ptr< FieldExtensionAlgorithm >
LevelSetMethodAlgorithm::getFieldExtensionAlgorithm(
  const tbox::Dimension& dim,
  const int field_handle,
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  int spatial_derivative_order,
  int tvd_runge_kutta_order,
  LSMLIB_REAL cfl_number,
  const LSMLIB_REAL stop_distance,
  const int max_iterations,
  const LSMLIB_REAL iteration_stop_tolerance,
  const bool verbose_mode,
  const string& object_name)
{
  if (!d_using_standard_level_set_fcn_integrator) {
    TBOX_ERROR( d_object_name
              << "::getFieldExtensionAlgorithm(): "
              << "this method may only be used when using the standard "
              << "LevelSetFunctionIntegrator."
              << endl );

    boost::shared_ptr< FieldExtensionAlgorithm > field_extension_alg = boost::shared_ptr< FieldExtensionAlgorithm> ();
    return field_extension_alg;
  }

  int level_set_fcn_handle;
  if (level_set_fcn == PHI) {
    level_set_fcn_handle = d_lsm_integrator_strategy->getPhiPatchDataHandle();
  } else {
    level_set_fcn_handle = d_lsm_integrator_strategy->getPsiPatchDataHandle();
  }

  // get default parameters if necesary
  if (spatial_derivative_type == UNKNOWN) {
    spatial_derivative_type = d_spatial_derivative_type;
  }
  if (spatial_derivative_order == 0) {
    spatial_derivative_order = d_spatial_derivative_order;
  }
  if (tvd_runge_kutta_order == 0) {
    tvd_runge_kutta_order = d_tvd_runge_kutta_order;
  }
  if (cfl_number == 0) {
    cfl_number = d_cfl_number;
  }

  boost::shared_ptr< FieldExtensionAlgorithm > field_extension_alg =
    boost::shared_ptr< FieldExtensionAlgorithm > (new FieldExtensionAlgorithm(
      d_patch_hierarchy,
      field_handle,
      level_set_fcn_handle,
      getControlVolumePatchDataHandle(),
      IntVector(d_patch_hierarchy->getDim(), 0),
      (SPATIAL_DERIVATIVE_TYPE) spatial_derivative_type,
      spatial_derivative_order,
      tvd_runge_kutta_order,
      cfl_number,
      stop_distance,
      max_iterations,
      iteration_stop_tolerance,
      verbose_mode,
      object_name));

  // add new object to list of FieldExtensionAlgorithms to be
  // reset by resetHierarchyConfiguration
  int length = d_field_extension_alg_list.size();
  d_field_extension_alg_list.resizeArray(length+1);
  d_field_extension_alg_list[length] = field_extension_alg;

  return field_extension_alg;
}





/* Destructor */
LevelSetMethodAlgorithm::~LevelSetMethodAlgorithm(){}


/* getPhiPatchDataHandle() */
int LevelSetMethodAlgorithm::getPhiPatchDataHandle() const
{
  return d_lsm_integrator_strategy->getPhiPatchDataHandle();
}


/* getPsiPatchDataHandle() */
int LevelSetMethodAlgorithm::getPsiPatchDataHandle() const
{
  return d_lsm_integrator_strategy->getPsiPatchDataHandle();
}


/* getControlVolumePatchDataHandle() */
int LevelSetMethodAlgorithm::getControlVolumePatchDataHandle() const
{
  return d_lsm_integrator_strategy->getControlVolumePatchDataHandle();
}


/* getStartTime() */
LSMLIB_REAL LevelSetMethodAlgorithm::getStartTime() const
{
  return d_lsm_integrator_strategy->getStartTime();
}


/* getEndTime() */
LSMLIB_REAL LevelSetMethodAlgorithm::getEndTime() const
{
  return d_lsm_integrator_strategy->getEndTime();
}


/* getCurrentTime() */
LSMLIB_REAL LevelSetMethodAlgorithm::getCurrentTime() const
{
  return d_lsm_integrator_strategy->getCurrentTime();
}


/* endTimeReached() */
bool LevelSetMethodAlgorithm::endTimeReached() const
{
  return d_lsm_integrator_strategy->endTimeReached();
}


/* numIntegrationStepsTaken() */
int LevelSetMethodAlgorithm::numIntegrationStepsTaken() const
{
  return d_lsm_integrator_strategy->numIntegrationStepsTaken();
}


/* getSpatialDerivativeType() */
int LevelSetMethodAlgorithm::getSpatialDerivativeType() const
{
  if (!d_using_standard_level_set_fcn_integrator) {
    TBOX_WARNING( d_object_name
              << "::getSpatialDerivativeType(): "
              << "not using standard LevelSetFunctionIntegrator..."
              << "return value meaningless."
              << endl );
  }

  return d_spatial_derivative_type;
}


/* getSpatialDerivativeOrder() */
int LevelSetMethodAlgorithm::getSpatialDerivativeOrder() const
{
  if (!d_using_standard_level_set_fcn_integrator) {
    TBOX_WARNING( d_object_name
              << "::getSpatialDerivativeOrder(): "
              << "not using standard LevelSetFunctionIntegrator..."
              << "return value meaningless."
              << endl );
  }

  return d_spatial_derivative_order;
}


/* getTVDRungeKuttaOrder() */
int LevelSetMethodAlgorithm::getTVDRungeKuttaOrder() const
{
  if (!d_using_standard_level_set_fcn_integrator) {
    TBOX_WARNING( d_object_name
              << "::getTVDRungeKuttaOrder(): "
              << "not using standard LevelSetFunctionIntegrator..."
              << "return value meaningless."
              << endl );
  }

  return d_tvd_runge_kutta_order;
}


/* getCFLNumber() */
LSMLIB_REAL LevelSetMethodAlgorithm::getCFLNumber() const
{
  if (!d_using_standard_level_set_fcn_integrator) {
    TBOX_WARNING( d_object_name
              << "::getCFLNumber(): "
              << "not using standard LevelSetFunctionIntegrator..."
              << "return value meaningless."
              << endl );
  }

  return d_cfl_number;
}


/* setBoundaryConditions() */
void LevelSetMethodAlgorithm::setBoundaryConditions(
  const IntVector& lower_bc,
  const IntVector& upper_bc,
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const int component)
{
  d_lsm_integrator_strategy->setBoundaryConditions(
    lower_bc, upper_bc, level_set_fcn, component);
}


/* initializeLevelSetMethodCalculation() */
void LevelSetMethodAlgorithm::initializeLevelSetMethodCalculation()
{
  d_lsm_gridding_strategy->initializePatchHierarchy(
    d_lsm_integrator_strategy->getCurrentTime());
}


/* computeStableDt() */
LSMLIB_REAL LevelSetMethodAlgorithm::computeStableDt()
{
  return d_lsm_integrator_strategy->computeStableDt();
}


/* advanceLevelSetFunctions() */
bool LevelSetMethodAlgorithm::advanceLevelSetFunctions(const LSMLIB_REAL dt)
{
  return d_lsm_integrator_strategy->advanceLevelSetFunctions(dt);
}


/* getReinitializationInterval() */
int LevelSetMethodAlgorithm::getReinitializationInterval() const
{
  return d_lsm_integrator_strategy->getReinitializationInterval();
}


/* setReinitializationInterval() */
void LevelSetMethodAlgorithm::setReinitializationInterval(
  const int interval)
{
  d_lsm_integrator_strategy->setReinitializationInterval(interval);
}


/* reinitializeLevelSetFunctions() */
void LevelSetMethodAlgorithm::reinitializeLevelSetFunctions(
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const int max_iterations)
{
  return d_lsm_integrator_strategy->reinitializeLevelSetFunctions(
           level_set_fcn, max_iterations);
}


/* getOrthogonalizationInterval() */
int LevelSetMethodAlgorithm::getOrthogonalizationInterval() const
{
  return d_lsm_integrator_strategy->getOrthogonalizationInterval();
}


/* setOrthogonalizationInterval() */
void LevelSetMethodAlgorithm::setOrthogonalizationInterval(
  const int interval)
{
  d_lsm_integrator_strategy->setOrthogonalizationInterval(interval);
}


/* orthogonalizeLevelSetFunctions() */
void LevelSetMethodAlgorithm::orthogonalizeLevelSetFunctions(
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const int max_reinit_iterations,
  const int max_ortho_iterations)
{
  return d_lsm_integrator_strategy->orthogonalizeLevelSetFunctions(
           level_set_fcn, max_reinit_iterations, max_ortho_iterations);
}


/* regridPatchHierarchy() */
void LevelSetMethodAlgorithm::regridPatchHierarchy()
{
  d_lsm_gridding_strategy->regridPatchHierarchy(
    d_lsm_integrator_strategy->getCurrentTime());
}

} // end LSMLIB namespace

#endif
