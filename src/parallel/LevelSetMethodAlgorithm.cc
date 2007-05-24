/*
 * File:        LevelSetMethodAlgorithm.cc
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.9 $
 * Modified:    $Date: 2006/08/02 19:05:13 $
 * Description: Implementation file for level set method integrator class
 */

#ifndef included_LevelSetMethodAlgorithm_cc
#define included_LevelSetMethodAlgorithm_cc

#include "LevelSetMethodAlgorithm.h" 

#ifdef LSMLIB_DEBUG_NO_INLINE
#include "LevelSetMethodAlgorithm.inline"
#endif

/****************************************************************
 *
 * Non-inline methods for LevelSetMethodAlgorithm class
 *
 ****************************************************************/

namespace LSMLIB {

/* Constructor for standard integrator and gridding strategy */
template<int DIM> 
LevelSetMethodAlgorithm<DIM>::LevelSetMethodAlgorithm(
  Pointer<Database> lsm_algorithm_input_db,
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  LevelSetMethodPatchStrategy<DIM>* patch_strategy,
  LevelSetMethodVelocityFieldStrategy<DIM>* velocity_field_strategy,
  const int num_level_set_fcn_components,
  const int codimension,
  const string& object_name)
{

#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!lsm_algorithm_input_db.isNull());
  assert(!patch_hierarchy.isNull());
  assert(patch_strategy);
  assert(velocity_field_strategy);
  assert(!object_name.empty());
#endif

  d_patch_hierarchy = patch_hierarchy;
  d_object_name = object_name;

  // create new LevelSetFunctionIntegrator object
  Pointer<Database> level_set_fcn_integrator_db =
    lsm_algorithm_input_db->getDatabase("LevelSetFunctionIntegrator");
  d_lsm_integrator_strategy = new LevelSetFunctionIntegrator<DIM>(
    level_set_fcn_integrator_db,
    patch_hierarchy,
    patch_strategy,
    velocity_field_strategy,
    num_level_set_fcn_components,
    codimension,
    object_name + "::LevelSetFunctionIntegrator");

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
  d_lsm_gridding_strategy = new LevelSetMethodGriddingAlgorithm<DIM>(
    lsm_algorithm_input_db->getDatabase("LevelSetMethodGriddingAlgorithm"),
    patch_hierarchy,
    d_lsm_integrator_strategy,
    object_name + "::LevelSetMethodGriddingAlgorithm");

  // register velocity_field_strategy with d_lsm_gridding_strategy
  d_lsm_gridding_strategy->registerVelocityFieldStrategy(
    velocity_field_strategy);

}


/* Constructor for custom integrator and gridding strategy */
template<int DIM> 
LevelSetMethodAlgorithm<DIM>::LevelSetMethodAlgorithm(
  Pointer< LevelSetFunctionIntegratorStrategy<DIM> > lsm_integrator_strategy,
  Pointer< LevelSetMethodGriddingStrategy<DIM> > lsm_gridding_strategy,
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
  d_patch_hierarchy.setNull(); 
  d_spatial_derivative_type = UNKNOWN;
  d_spatial_derivative_order = 0;
  d_tvd_runge_kutta_order = 0;
  d_cfl_number = 0.0;

}


/* printClassData() */
template<int DIM> 
void LevelSetMethodAlgorithm<DIM>::printClassData(ostream& os) const
{
  os << "\n===================================" << endl;
  os << "LevelSetMethodAlgorithm<DIM>" << endl;

  os << "Object Pointers" << endl;
  os << "---------------" << endl;
  os << "(LevelSetMethodAlgorithm*) this = "
     << (LevelSetMethodAlgorithm*) this << endl;
  os << "d_lsm_integrator_strategy = "
     << d_lsm_integrator_strategy.getPointer() << endl;
  os << "d_lsm_gridding_strategy = " 
     << d_lsm_gridding_strategy.getPointer() << endl;
  os << "===================================" << endl << endl;
}


/* resetHierarchyConfiguration() */
template<int DIM> 
void LevelSetMethodAlgorithm<DIM>::resetHierarchyConfiguration(
  const Pointer< PatchHierarchy<DIM> > hierarchy,
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
template<int DIM> 
Pointer< FieldExtensionAlgorithm<DIM> > 
LevelSetMethodAlgorithm<DIM>::getFieldExtensionAlgorithm(
  Pointer<Database> input_db,
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

    Pointer< FieldExtensionAlgorithm<DIM> > field_extension_alg;
    field_extension_alg.setNull();
    return field_extension_alg;
  }

  int level_set_fcn_handle;
  if (level_set_fcn == PHI) {
    level_set_fcn_handle = d_lsm_integrator_strategy->getPhiPatchDataHandle();
  } else {
    level_set_fcn_handle = d_lsm_integrator_strategy->getPsiPatchDataHandle();
  }

  Pointer< FieldExtensionAlgorithm<DIM> > field_extension_alg =
    new FieldExtensionAlgorithm<DIM>(
      input_db,
      d_patch_hierarchy,
      field_handle,
      level_set_fcn_handle,  
      getControlVolumePatchDataHandle(),
      object_name);

  // add new object to list of FieldExtensionAlgorithms to be 
  // reset by resetHierarchyConfiguration
  int length = d_field_extension_alg_list.size();
  d_field_extension_alg_list.resizeArray(length+1);
  d_field_extension_alg_list[length] = field_extension_alg;

  return field_extension_alg;
}


/* getFieldExtensionAlgorithm() */
template<int DIM> 
Pointer< FieldExtensionAlgorithm<DIM> > 
LevelSetMethodAlgorithm<DIM>::getFieldExtensionAlgorithm(
  const int field_handle,
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  int spatial_derivative_order,
  int tvd_runge_kutta_order,
  double cfl_number,
  const double stop_distance,
  const int max_iterations,
  const double iteration_stop_tolerance,
  const bool verbose_mode,
  const string& object_name)
{
  if (!d_using_standard_level_set_fcn_integrator) {
    TBOX_ERROR( d_object_name
              << "::getFieldExtensionAlgorithm(): "
              << "this method may only be used when using the standard "
              << "LevelSetFunctionIntegrator."
              << endl );

    Pointer< FieldExtensionAlgorithm<DIM> > field_extension_alg;
    field_extension_alg.setNull();
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

  Pointer< FieldExtensionAlgorithm<DIM> > field_extension_alg =
    new FieldExtensionAlgorithm<DIM>(
      d_patch_hierarchy,
      field_handle,
      level_set_fcn_handle,  
      getControlVolumePatchDataHandle(),
      (SPATIAL_DERIVATIVE_TYPE) spatial_derivative_type,
      spatial_derivative_order,
      tvd_runge_kutta_order,
      cfl_number,
      stop_distance,
      max_iterations,
      iteration_stop_tolerance,
      verbose_mode,
      object_name);

  // add new object to list of FieldExtensionAlgorithms to be 
  // reset by resetHierarchyConfiguration
  int length = d_field_extension_alg_list.size();
  d_field_extension_alg_list.resizeArray(length+1);
  d_field_extension_alg_list[length] = field_extension_alg;

  return field_extension_alg;
}


} // end LSMLIB namespace

#endif
