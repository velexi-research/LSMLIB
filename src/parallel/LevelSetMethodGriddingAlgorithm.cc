/*
 * File:        LevelSetMethodGriddingAlgorithm.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.6 $
 * Modified:    $Date: 2007/03/24 01:28:10 $
 * Description: Implementation of the level set method grid management class
 */

#ifndef included_LevelSetMethodGriddingAlgorithm_cc
#define included_LevelSetMethodGriddingAlgorithm_cc

#include "LevelSetMethodGriddingAlgorithm.h" 
#include "BergerRigoutsos.h" 
#include "CellData.h" 
#include "LoadBalancer.h" 
#include "tbox/RestartManager.h" 

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif


/****************************************************************
 *
 * Implementation for LevelSetMethodGriddingAlgorithm methods.
 *
 ****************************************************************/

namespace LSMLIB {

/* Constructor */
template <int DIM>
LevelSetMethodGriddingAlgorithm<DIM>::LevelSetMethodGriddingAlgorithm(
  Pointer<Database> input_db,
  Pointer< BasePatchHierarchy<DIM> > patch_hierarchy,
  Pointer< LevelSetFunctionIntegratorStrategy<DIM> > lsm_integrator_strategy,
  const string& object_name) :
StandardTagAndInitialize<DIM>(object_name, 
  lsm_integrator_strategy, input_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(input_db);
  assert(patch_hierarchy);
  assert(lsm_integrator_strategy);
  assert(!object_name.empty());
#endif

  // set d_object_name
  d_object_name = object_name;

  // set d_patch_hierarchy
  d_patch_hierarchy = patch_hierarchy;

  // set level set method integrator and initialize 
  d_lsm_integrator_strategy = lsm_integrator_strategy;

  // make sure that the d_velocity_field_strategies
  // is empty.
  d_velocity_field_strategies.setNull();

  // read input parameters
  getFromInput(input_db);

  /*
   * Set up LevelSetMethodGriddingAlgorithm
   */ 

  // construct box generator and load balancer objects
  Pointer< BergerRigoutsos<DIM> > box_generator = new BergerRigoutsos<DIM>();
  Pointer<Database> load_balancer_input_db;
  Pointer< LoadBalancer<DIM> > load_balancer;
  if (input_db->isDatabase("LoadBalancer")) {
    load_balancer_input_db = input_db->getDatabase("LoadBalancer");
    load_balancer = new 
      LoadBalancer<DIM> ("load balancer", load_balancer_input_db);
  } else {
    load_balancer = new LoadBalancer<DIM> ("load balancer");
  }

  // construct gridding algorithm using "this" as the 
  // TagAndInitializeStrategy.  
  // NOTE: "this" is passed to the SAMRAI::mesh::GriddingAlgorithm as an 
  // unmanaged smart-pointer to avoid having two distinct smart-pointers 
  // reference the same SAMRAI::mesh::TagAndInitializeStrategy<DIM>*.
  d_gridding_alg = new SAMRAI::mesh::GriddingAlgorithm<DIM>(
    "LevelSetMethodGriddingAlgorithm", 
    input_db, 
    Pointer< LevelSetMethodGriddingAlgorithm<DIM> >(this, false),
    box_generator, 
    load_balancer);
}


/* registerVelocityFieldStrategy() */
template <int DIM>
void LevelSetMethodGriddingAlgorithm<DIM>::registerVelocityFieldStrategy(
  LevelSetMethodVelocityFieldStrategy<DIM>* velocity_field_strategy)
{
  int size = d_velocity_field_strategies.size();
  d_velocity_field_strategies.resizeArray(size+1);

  // create an "unmanaged" smart-pointer to wrap the pointer to the
  // velocity field strategy
  d_velocity_field_strategies[size] = 
    Pointer< LevelSetMethodVelocityFieldStrategy<DIM> >(
      velocity_field_strategy,false);
}


/* initializePatchHierarchy() */
template <int DIM>
void LevelSetMethodGriddingAlgorithm<DIM>::initializePatchHierarchy(
  const LSMLIB_REAL time)
{
  /*
   * Construct and initialize the levels of hierarchy.
   */
  if (RestartManager::getManager()->isFromRestart()) { 
    // from restart

    d_patch_hierarchy->getFromRestart(d_gridding_alg->getMaxLevels());

    resetHierarchyConfiguration(d_patch_hierarchy,
                                0,
                                d_patch_hierarchy->getFinestLevelNumber());

  } else {  
    // not from restart

    d_gridding_alg->makeCoarsestLevel(d_patch_hierarchy,time);
    bool done = false;
    for (int level_num = 1; 
         d_gridding_alg->levelCanBeRefined(level_num) && !done;
         level_num++) {

      plog << "Adding finer level with level_num = " << level_num << endl;
      d_gridding_alg->makeFinerLevel(d_patch_hierarchy, time, true, 0);

      plog << "Just added finer level with level_num = " << level_num << endl;
      done = !(d_patch_hierarchy->finerLevelExists(level_num));
    }


    //KTC - temporarily left out - need to fix for AMR
    /*
     * After data on each level is initialized at simulation start time,
     * coarser levels are synchronized with finer levels that didn't exist
     * when the coarser level initial data was set.  This synchronization
     * process is defined by the integration algorithm.
     */
    if (d_patch_hierarchy->getFinestLevelNumber() > 0) {
      // "true" argument: const bool initial_time = true;
      /*
          d_level_set_alg->synchronizeNewLevels(
            0,
            d_patch_hierarchy->
            getFinestLevelNumber(),
            true);
          d_velocity_field_alg->synchronizeNewLevels(
            0,
            d_patch_hierarchy->
            getFinestLevelNumber(),
            true);
      */
    }

  }
}


/* 
 * regridPatchHierarchy() invokes the method
 * SAMRAI::LevelSetMethodGriddingAlgorithm::regridAllFinerLevels()
 * to carry out the regridding.   
 */
template<int DIM> 
void LevelSetMethodGriddingAlgorithm<DIM>::regridPatchHierarchy(
  LSMLIB_REAL time)
{
  int num_levels = d_patch_hierarchy->getNumberLevels();
  Array<int> tag_buffer(num_levels, true);
  for (int ln=0; ln < num_levels ; ln++) 
    tag_buffer[ln] = d_gridding_alg->getProperNestingBuffer(ln);

  d_gridding_alg->regridAllFinerLevels(
    d_patch_hierarchy, 
    0,    // regrid all levels finer than the coarsest level
    time,
    tag_buffer);
}


/* initializeLevelData() */
template<int DIM> 
void LevelSetMethodGriddingAlgorithm<DIM>::initializeLevelData(
  const Pointer< BasePatchHierarchy<DIM> > hierarchy,
  const int level_number,
  const LSMLIB_REAL init_data_time,
  const bool can_be_refined,
  const bool initial_time,
  const Pointer< BasePatchLevel<DIM> > old_level,
  const bool allocate_data) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!hierarchy.isNull());
  assert( (level_number >= 0)
          && (level_number <= hierarchy->getFinestLevelNumber()) );
  if ( !(old_level.isNull()) ) {
    assert( level_number == old_level->getLevelNumber() );            
  }
  assert(!(hierarchy->getPatchLevel(level_number)).isNull()); 
#endif

  if (!d_lsm_integrator_strategy.isNull()) { 
    d_lsm_integrator_strategy->initializeLevelData(hierarchy, 
      level_number, init_data_time, can_be_refined, initial_time, old_level, 
      allocate_data);
  }

  int phi_handle;
  int psi_handle;
  d_lsm_integrator_strategy->preprocessInitializeVelocityField(
    phi_handle,psi_handle, hierarchy,level_number);
  int num_velocity_field_strategies = 
    d_velocity_field_strategies.size();
  for (int k = 0; k < num_velocity_field_strategies; k++) {
    if (!d_velocity_field_strategies[k].isNull()) {
      d_velocity_field_strategies[k]
        ->initializeLevelData(hierarchy,
                              level_number,
                              init_data_time,
                              phi_handle,
                              psi_handle,
                              can_be_refined,
                              initial_time,
                              old_level,
                              allocate_data);
    }
  } // end loop over the elements of d_velocity_field_strategies
  d_lsm_integrator_strategy->postprocessInitializeVelocityField(
    hierarchy, level_number);

}


/* resetHierarchyConfiguration() */
template<int DIM> 
void LevelSetMethodGriddingAlgorithm<DIM>::resetHierarchyConfiguration(
  const Pointer< BasePatchHierarchy<DIM> > hierarchy,
  const int coarsest_level,
  const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!hierarchy.isNull());
  assert( (coarsest_level >= 0)
          && (coarsest_level <= finest_level)
          && (finest_level <= hierarchy->getFinestLevelNumber()) );
  for (int ln0 = 0; ln0 <= finest_level; ln0++) {
    assert(!(hierarchy->getPatchLevel(ln0)).isNull());
  }
#endif

  if (!d_lsm_integrator_strategy.isNull()) { 
    d_lsm_integrator_strategy->resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level);
  }

  int num_velocity_field_strategies = 
    d_velocity_field_strategies.size();
  for (int k = 0; k < num_velocity_field_strategies; k++) {
    if (!d_velocity_field_strategies[k].isNull()) {
      d_velocity_field_strategies[k]
        ->resetHierarchyConfiguration(hierarchy,
                                      coarsest_level,
                                      finest_level);
    }
  } // end loop over the elements of d_velocity_field_strategies
}


/* tagCellsForRefinements() 
 *
 * The implementation of this method is essentially copied from 
 * SAMRAI::StandardTagAndInitialize::tagCellsForRefinement().  The
 * main difference is that Richardson extrapolation is not a valid
 * way to set the refinement cells.
 */
template<int DIM> 
void LevelSetMethodGriddingAlgorithm<DIM>::tagCellsForRefinement(
  const Pointer< BasePatchHierarchy<DIM> > hierarchy,
  const int level_number,
  const LSMLIB_REAL regrid_time,
  const int tag_index, 
  const bool initial_time,
  const bool coarsest_sync_level,
  const bool can_be_refined,
  const LSMLIB_REAL regrid_start_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!(hierarchy.isNull()));
  assert( (level_number>=0) 
          && (level_number <= hierarchy->getFinestLevelNumber()) );
  assert(!(hierarchy->getPatchLevel(level_number).isNull()));
  assert(tag_index>=0);
#endif

  if (d_use_richardson_extrapolation)  {
    TBOX_WARNING(   d_object_name 
                 << "::tagCellsForRefinement(): "
                 << "Richardson extrapolation is NOT implemented in "
                 << "the LevelSetMethodGriddingAlgorithm.");
  }

  if (d_use_gradient_detector) {

    NULL_USE(regrid_start_time);
    NULL_USE(can_be_refined);
    NULL_USE(coarsest_sync_level);

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!d_lsm_integrator_strategy.isNull());
#endif
    d_lsm_integrator_strategy->applyGradientDetector(hierarchy, 
      level_number, regrid_time, tag_index, initial_time, 
      d_use_richardson_extrapolation);

    int num_velocity_field_strategies = 
      d_velocity_field_strategies.size();
    for (int k = 0; k < num_velocity_field_strategies; k++) {
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!d_velocity_field_strategies[k].isNull());
#endif

      d_velocity_field_strategies[k]
        ->tagCellsForRefinement(hierarchy,
                                level_number,
                                tag_index);
    }
  }


  /*
   * If user-supplied refine boxes are to be used, get refine box information
   * from the TagAndInitializeStrategy class, from which this class is
   * derived.
   */
  if (d_use_refine_boxes) {

    hier::BoxArray<DIM> refine_boxes;
    getUserSuppliedRefineBoxes(refine_boxes, level_number, regrid_time);

    tbox::Pointer< hier::PatchLevel<DIM> > level =
       hierarchy->getPatchLevel(level_number);

    for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip());

      tbox::Pointer< pdat::CellData<DIM,int> > tag_data =
        patch->getPatchData(tag_index);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert( !(tag_data.isNull()) );
#endif

      for (int ib = 0; ib < refine_boxes.getNumberOfBoxes(); ib++) {
        hier::Box<DIM> intersection =
          refine_boxes.getBox(ib) * tag_data->getBox();
        if ( !(intersection.empty()) ) {
            tag_data->fill(1, intersection);
        }
      }
    }
  }

}


/* getFromInput() 
 *
 * This method is essentially copied from 
 * StandardTagAndInitialize::getFromInput(). 
 * It just does not read in the input related to 
 * the explicitly-specified refinement boxes.
 */
template<int DIM> 
void LevelSetMethodGriddingAlgorithm<DIM>::getFromInput(
  Pointer<Database> input_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!input_db.isNull());
#endif

  tbox::Array<string> tagging_method;
  if (input_db->keyExists("tagging_method")) {
    tagging_method = input_db->getStringArray("tagging_method");
  }

  if (tagging_method.getSize() > 3) {
    TBOX_ERROR(  d_object_name 
              << "::getFromInput(): "
              << tagging_method.getSize() 
              << " entries specified in `tagging_method' input.  "
              << "Maximum allowable is 3."
              << endl);
  }

  d_use_gradient_detector = false;
  d_use_richardson_extrapolation = false;
  d_use_refine_boxes = false;

  /*
   * Check tagging method input.
   */

  bool found_method = false;
  for (int i = 0; i < tagging_method.getSize(); i++) {

    if ( tagging_method[i] == "GRADIENT_DETECTOR" ) {

      d_use_gradient_detector = true;
      found_method = true;

    }

    if (tagging_method[i] == "RICHARDSON_EXTRAPOLATION") {

      d_use_richardson_extrapolation = true;
      found_method = true;

    }

    if (tagging_method[i] == "REFINE_BOXES") {

      d_use_refine_boxes = true;
      found_method = true;

    }
  }

  /*
   * Check for valid entries
   */
  if (!found_method) {
    TBOX_WARNING(d_object_name << "::getFromInput(): "
      << "No `tagging_method' entry specified, so cell tagging \n"
      << "will NOT be performed.  If you wish to invoke cell \n"
      << "tagging, you must enter one or more valid tagging \n"
      << "methods, of type GRADIENT_DETECTOR, "
      << "RICHARDSON_EXTRAPOLATION, or REFINE_BOXES\n"
      << "See class header for details.\n");
  }
}

/* Copy Constructor */
template <int DIM>
LevelSetMethodGriddingAlgorithm<DIM>::LevelSetMethodGriddingAlgorithm(
  const LevelSetMethodGriddingAlgorithm& lsm_gridding_algorithm) :
StandardTagAndInitialize<DIM>(
  string("copy constructor"),
  Pointer< LevelSetFunctionIntegratorStrategy<DIM> >(0),
  Pointer<Database>(0) )
{}

} // end LSMLIB namespace

#endif
