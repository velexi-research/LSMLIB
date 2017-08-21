/*
 * File:        LevelSetMethodGriddingAlgorithm.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Implementation of the level set method grid management class
 */

#ifndef included_LevelSetMethodGriddingAlgorithm_cc
#define included_LevelSetMethodGriddingAlgorithm_cc

// Class header
#include "LevelSetMethodGriddingAlgorithm.h"

// Standard library headers
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// SAMRAI headers
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"

// LSMLIB headers
#include "LevelSetFunctionIntegratorStrategy.h"
#include "LevelSetMethodVelocityFieldStrategy.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <cassert>
#endif
#endif

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchData; } }

/****************************************************************
 *
 * Implementation for LevelSetMethodGriddingAlgorithm methods.
 *
 ****************************************************************/

namespace LSMLIB {

/* Constructor */
LevelSetMethodGriddingAlgorithm::LevelSetMethodGriddingAlgorithm(
  boost::shared_ptr<tbox::Database> input_db,
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  boost::shared_ptr<LevelSetFunctionIntegratorStrategy> lsm_integrator_strategy,
  const string& object_name) :
  StandardTagAndInitialize(object_name,
                           lsm_integrator_strategy.get(),
                           input_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(input_db);
  assert(patch_hierarchy);
  assert(lsm_integrator_strategy);
  assert(!object_name.empty());
#endif
  cout << "LevelSetMethodGriddingAlgorithm::constructor "
       << lsm_integrator_strategy.get() << endl;
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
  boost::shared_ptr<mesh::BergerRigoutsos> box_generator =
    boost::shared_ptr<mesh::BergerRigoutsos>(
      new mesh::BergerRigoutsos(d_patch_hierarchy->getDim(),input_db));

  boost::shared_ptr<tbox::Database> load_balancer_input_db;

  boost::shared_ptr<mesh::ChopAndPackLoadBalancer> load_balancer;
  if (input_db->isDatabase("LoadBalancer")) {
    load_balancer_input_db = input_db->getDatabase("LoadBalancer");
    load_balancer = boost::shared_ptr<mesh::ChopAndPackLoadBalancer> (
       new mesh::ChopAndPackLoadBalancer(d_patch_hierarchy->getDim(),
                                         "LoadBalancer",
                                         load_balancer_input_db));

  } else {
    throw runtime_error(
      "LoadBalancer' section not found in input database");
  }

  // construct gridding algorithm using "this" as the TagAndInitializeStrategy
  cout << this << endl;
  d_gridding_alg = boost::shared_ptr<mesh::GriddingAlgorithm> (
    new mesh::GriddingAlgorithm(
      patch_hierarchy,
      "LevelSetMethodGriddingAlgorithm",
      input_db,
//      boost::shared_ptr<TagAndInitializeStrategy>(this),
      boost::shared_ptr<LevelSetMethodGriddingAlgorithm>(this),
      box_generator,
      load_balancer));
}


/* registerVelocityFieldStrategy() */
void LevelSetMethodGriddingAlgorithm::registerVelocityFieldStrategy(
  LevelSetMethodVelocityFieldStrategy* velocity_field_strategy)
{
  int size = d_velocity_field_strategies.size();
  d_velocity_field_strategies.resizeArray(size+1);

  // velocity field strategy
 d_velocity_field_strategies[size] =
   boost::shared_ptr<LevelSetMethodVelocityFieldStrategy>(
      velocity_field_strategy);
}


/* initializePatchHierarchy() */
void LevelSetMethodGriddingAlgorithm::initializePatchHierarchy(
  const LSMLIB_REAL time)
{
  /*
   * Construct and initialize the levels of hierarchy.
   */
  if (tbox::RestartManager::getManager()->isFromRestart()) {
    cout << "initializePatchHierarchy from restart" << endl;

    // from restart
    d_patch_hierarchy->getMaxNumberOfLevels();

    resetHierarchyConfiguration(d_patch_hierarchy,
                                0,
                                d_patch_hierarchy->getFinestLevelNumber());

  } else {
    cout << "initializePatchHierarchy not from restart 1" << endl;
    // not from restart
    d_gridding_alg->makeCoarsestLevel(time);
    cout << "initializePatchHierarchy not from restart 2" << endl;

    bool done = false;
    for (int level_num = 1;
         d_patch_hierarchy->levelCanBeRefined(level_num) && !done;
         level_num++) {
      tbox::plog << "Adding finer level with level_num = "
                 << level_num << endl;
      d_gridding_alg->makeFinerLevel(time, true, 0,
                                     d_patch_hierarchy->getFinestLevelNumber());

      tbox::plog << "Just added finer level with level_num = "
                 << level_num << endl;
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
void LevelSetMethodGriddingAlgorithm::regridPatchHierarchy(
  const double time)
{
  int num_levels = d_patch_hierarchy->getNumberOfLevels();
  vector<int> tag_buffer(num_levels, true);
  for (int ln=0; ln < num_levels ; ln++)
    tag_buffer[ln] = d_patch_hierarchy->getProperNestingBuffer(ln);

  d_gridding_alg->regridAllFinerLevels(
    0, // regrid all levels finer than the coarsest level
    tag_buffer,
    0, //unused, look for iteration_number if iteration is needed
    time);
}


/* initializeLevelData() */
void LevelSetMethodGriddingAlgorithm::initializeLevelData(
  const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const int level_number,
  const double init_data_time,
  const bool can_be_refined,
  const bool initial_time,
  const boost::shared_ptr<hier::PatchLevel> old_level,
  const bool allocate_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(hierarchy!=NULL);
  assert( (level_number >= 0)
          && (level_number <= hierarchy->getFinestLevelNumber()) );
  if ( old_level!=NULL ) {
    assert( level_number == old_level->getLevelNumber() );
  }
  assert((hierarchy->getPatchLevel(level_number))!=NULL);
#endif
  cout << "LevelSetMethodGriddingAlgorithm::initializeLevelData"
       << d_lsm_integrator_strategy.get()
       << endl;
  if (d_lsm_integrator_strategy!=NULL) {
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
  cout << "LevelSetMethodGriddingAlgorithm::initializeLevelData "
       << num_velocity_field_strategies << endl;
  for (int k = 0; k < num_velocity_field_strategies; k++) {
    cout << "LevelSetMethodGriddingAlgorithm::initializeLevelData "
         << d_velocity_field_strategies[k] << endl;
    if (d_velocity_field_strategies[k]!=NULL) {
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
void LevelSetMethodGriddingAlgorithm::resetHierarchyConfiguration(
  const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const int coarsest_level,
  const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(hierarchy!=NULL);
  assert( (coarsest_level >= 0)
          && (coarsest_level <= finest_level)
          && (finest_level <= hierarchy->getFinestLevelNumber()) );
  for (int ln0 = 0; ln0 <= finest_level; ln0++) {
    assert((hierarchy->getPatchLevel(ln0))!=NULL);
  }
#endif

  if (d_lsm_integrator_strategy!=NULL) {
    d_lsm_integrator_strategy->resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level);
  }

  int num_velocity_field_strategies =
    d_velocity_field_strategies.size();
  for (int k = 0; k < num_velocity_field_strategies; k++) {
    if (d_velocity_field_strategies[k]!=NULL) {
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
void LevelSetMethodGriddingAlgorithm::tagCellsForRefinement(
  const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const int level_number,
  const double regrid_time,
  const int tag_index,
  const bool initial_time,
  const bool coarsest_sync_level,
  const bool can_be_refined,
  const double regrid_start_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(hierarchy!=NULL);
  assert( (level_number>=0)
          && (level_number <= hierarchy->getFinestLevelNumber()) );
  assert((hierarchy->getPatchLevel(level_number)!=NULL));
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
    assert(d_lsm_integrator_strategy!=NULL);
#endif
    d_lsm_integrator_strategy->applyGradientDetector(hierarchy,
      level_number, regrid_time, tag_index, initial_time,
      d_use_richardson_extrapolation);

    int num_velocity_field_strategies =
      d_velocity_field_strategies.size();
    for (int k = 0; k < num_velocity_field_strategies; k++) {
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(d_velocity_field_strategies[k]!=NULL);
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
    hier::BoxContainer  refine_boxes;
    getUserSuppliedRefineBoxes(refine_boxes, level_number,0, regrid_time);
// 0 is an unused regrid_cycle
    boost::shared_ptr<hier::PatchLevel> level =
       hierarchy->getPatchLevel(level_number);

    for (hier::PatchLevel::Iterator ip(level->begin()); ip!=level->end(); ip++) {
      boost::shared_ptr<hier::Patch> patch = *ip;

   boost::shared_ptr<pdat::CellData<int>> tag_data(
      BOOST_CAST<pdat::CellData<int>, hier::PatchData>(
         patch->getPatchData(tag_index)));
#ifdef DEBUG_CHECK_ASSERTIONS
      assert( tag_data!=NULL );
#endif
hier::Box pbox = patch->getBox();
      for (hier::BoxContainer::iterator ib = refine_boxes.begin();
              ib != refine_boxes.end(); ++ib) {
            hier::Box  intersection = pbox* (*ib);

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

const int nsize = 3; // size of arrays
void LevelSetMethodGriddingAlgorithm::getFromInput(
  boost::shared_ptr<tbox::Database> input_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(input_db!=NULL);
#endif
  vector<string> tagging_method;
  tagging_method.resize(nsize);
  if (input_db->keyExists("tagging_method")) {
    tagging_method = input_db->getStringVector("tagging_method");
  }

 /* if (tagging_method.getSize() > 3) {
    TBOX_ERROR(  d_object_name
              << "::getFromInput(): "
              << tagging_method.getSize()
              << " entries specified in `tagging_method' input.  "
              << "Maximum allowable is 3."
              << endl);
  }*/
  d_use_gradient_detector = false;
  d_use_richardson_extrapolation = false;
  d_use_refine_boxes = false;

  /*
   * Check tagging method input.
   */

  bool found_method = false;
  for (int i = 0; i < nsize; i++) {

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


} // end LSMLIB namespace

#endif
