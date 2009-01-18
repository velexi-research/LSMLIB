/*
 * File:        main_integrals.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.5 $
 * Modified:    $Date$
 * Description: 2d example program for Level Set Method computing integrals
 */

/************************************************************************
 *
 * This example program demonstrates how to use the LSMLIB C++ classes
 * to for a 2D problem to compute area and line integrals over regions 
 * defined by the zero level set.
 *
 ************************************************************************/

// SAMRAI Configuration 
#include "SAMRAI_config.h"

/* 
 * Headers for basic SAMRAI objects
 */

// variables and variable management
#include "CellVariable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"

// geometry and patch hierarchy
#include "CartesianGridGeometry.h" 
#include "PatchHierarchy.h"

// basic SAMRAI classes
#include "IntVector.h" 
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "tbox/Database.h" 
#include "tbox/InputDatabase.h" 
#include "tbox/InputManager.h" 
#include "tbox/MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAIManager.h"
#include "tbox/Utilities.h"


// Headers for level set method
// LevelSetMethod configuration header must be included
// before any other LevelSetMethod header files
#include "LSMLIB_config.h"
#include "LevelSetMethodAlgorithm.h"
#include "LevelSetMethodToolbox.h"
#include "VelocityFieldModule.h"
#include "PatchModule.h"


// namespaces
using namespace std;
using namespace SAMRAI;
using namespace geom;
using namespace hier;
using namespace pdat;
using namespace tbox;
using namespace LSMLIB;


int main(int argc, char *argv[])
{

  /*
   * Initialize MPI and SAMRAI, enable logging, and process command line.
   */
  tbox::MPI::init(&argc, &argv);
  tbox::MPI::initialize();
  SAMRAIManager::startup();

  string input_filename;
  string restart_read_dirname;
  int restore_num = 0;

  bool is_from_restart = false;

  if ( (argc != 2) && (argc != 4) ) {
    pout << "USAGE:  " << argv[0] << " <input filename> "
         << "\n"
         << "<restart dir> <restore number> [options]\n"
         << "  options:\n"
         << "  none at this time"
         << endl;
    tbox::MPI::abort();
    return (-1);
  } else {
    input_filename = argv[1];
    if (argc == 4) {
      restart_read_dirname = argv[2];
      restore_num = atoi(argv[3]);
      is_from_restart = true;
    }
  }

  /*
   * Create input database and parse all data in input file.  
   */
  Pointer<Database> input_db = new InputDatabase("input_db");
  InputManager::getManager()->parseInputFile(input_filename, input_db);

  /*
   * Read in the input from the "Main" section of the input database.  
   */
  Pointer<Database> main_db = input_db->getDatabase("Main");

  /* 
   * The base_name variable is a base name for all name strings in 
   * this program.
   */
   string base_name = "unnamed";
   base_name = main_db->getStringWithDefault("base_name", base_name);

  /*
   * Start logging.
   */
   const string log_file_name = base_name + ".log";
   bool log_all_nodes = false;
   log_all_nodes = main_db->getBoolWithDefault("log_all_nodes", log_all_nodes);
   if (log_all_nodes) {
      PIO::logAllNodes(log_file_name);
   } else {
      PIO::logOnlyNodeZero(log_file_name);
   }

  /*
   * Get the restart manager and root restart database.  If run is from 
   * restart, open the restart file.
   */

  RestartManager* restart_manager = RestartManager::getManager();
  if (is_from_restart) {
    restart_manager->
       openRestartFile(restart_read_dirname, restore_num, 
                       tbox::MPI::getNodes() );
  }

  // log the command-line args
  plog << "input_filename = " << input_filename << endl;
  plog << "restart_read_dirname = " << restart_read_dirname << endl;
  plog << "restore_num = " << restore_num << endl;

  /*
   *  Create major algorithm and data objects. 
   */

  Pointer< CartesianGridGeometry<2> > grid_geometry =
    new CartesianGridGeometry<2>(
      base_name+"::CartesianGeometry",
      input_db->getDatabase("CartesianGeometry"));
  plog << "CartesianGridGeometry:" << endl;
  grid_geometry->printClassData(plog);

  Pointer< PatchHierarchy<2> > patch_hierarchy =
    new PatchHierarchy<2>(base_name+"::PatchHierarchy",
                             grid_geometry);

  VelocityFieldModule* velocity_field_module = 
    new VelocityFieldModule( 
      input_db->getDatabase("VelocityFieldModule"),
      patch_hierarchy,
      grid_geometry,
      base_name+"::VelocityFieldModule");
  plog << "VelocityFieldModule:" << endl;
  velocity_field_module->printClassData(plog);

  PatchModule* patch_module = 
    new PatchModule(
      input_db->getDatabase("PatchModule"),
      base_name+"::PatchModule");
  plog << "PatchModule:" << endl;
  patch_module->printClassData(plog);

  int num_level_set_fcn_components = 1;
  int codimension = 1;
  Pointer< LevelSetMethodAlgorithm<2> > lsm_algorithm = 
    new LevelSetMethodAlgorithm<2>( 
      input_db->getDatabase("LevelSetMethodAlgorithm"),
      patch_hierarchy,
      patch_module,
      velocity_field_module,
      num_level_set_fcn_components,
      codimension,
      base_name+"::LevelSetMethodAlgorithm");
  plog << "LevelSetMethodAlgorithm:" << endl;
  lsm_algorithm->printClassData(plog);


  /*
   * After creating all objects and initializing their state, 
   * print the input database and variable database contents to the 
   * log file.
   */
  plog << "\nCheck input data and variables before simulation:" << endl;
  plog << "Input database..." << endl;
  input_db->printClassData(plog);
  plog << "\nVariable database..." << endl;
  VariableDatabase<2>::getDatabase()->printClassData(plog);


  // get PatchData handles
  int phi_handle = lsm_algorithm->getPhiPatchDataHandle();
  int velocity_handle = velocity_field_module
      ->getExternalVelocityFieldPatchDataHandle(0);
  int control_volume_handle = 
    lsm_algorithm->getControlVolumePatchDataHandle();

  /*
   * Initialize level set method calculation
   */ 
  lsm_algorithm->initializeLevelSetMethodCalculation();

  /*
   * Set up SAMRAI variables for grad(phi)
   */ 
  VariableDatabase<2> *var_db = VariableDatabase<2>::getDatabase();
  Pointer<VariableContext> plus_context = var_db->getContext("PLUS");
  Pointer<VariableContext> minus_context = var_db->getContext("MINUS");
  Pointer<VariableContext> scratch_context = 
    var_db->getContext("SCRATCH");
  Pointer< CellVariable<2,LSMLIB_REAL> > phi_variable =
    new CellVariable<2,LSMLIB_REAL>("phi",1);
  Pointer< CellVariable<2,LSMLIB_REAL> > grad_phi_variable =
    new CellVariable<2,LSMLIB_REAL>("grad(phi)",2);
 
  int phi_scratch_handle = var_db->registerVariableAndContext(
    phi_variable, scratch_context, IntVector<2>(3));
 int grad_phi_plus_handle = var_db->registerVariableAndContext(
    grad_phi_variable, plus_context, IntVector<2>(0));
  int grad_phi_minus_handle = var_db->registerVariableAndContext(
    grad_phi_variable, minus_context, IntVector<2>(0));
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    Pointer< PatchLevel<2> > level =
      patch_hierarchy->getPatchLevel(ln);
    level->allocatePatchData( phi_scratch_handle );
    level->allocatePatchData( grad_phi_plus_handle );
    level->allocatePatchData( grad_phi_minus_handle );
  }

  /*
   * Compute grad(phi) 
   */ 

  // create and fill scratch space 
  RefineAlgorithm<2> refine_alg;
  Pointer< RefineOperator<2> > refine_op =  
    grid_geometry->lookupRefineOperator(phi_variable, "LINEAR_REFINE");
  refine_alg.registerRefine(
    phi_scratch_handle, phi_handle, phi_scratch_handle, refine_op);
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    Pointer< PatchLevel<2> > level =
      patch_hierarchy->getPatchLevel(ln);
    Pointer< RefineSchedule<2> > sched = refine_alg.createSchedule(
      level, ln-1, patch_hierarchy, 0);
    sched->fillData(0.0,true); // physical boundary conditions set
  }

  // compute plus and minus spatial derivatives
  LevelSetMethodToolbox<2>::computePlusAndMinusSpatialDerivatives(
    patch_hierarchy,
    LSMLIB::WENO,
    5, // fifth-order
    grad_phi_plus_handle,
    grad_phi_minus_handle,
    phi_scratch_handle);
/*
  LevelSetMethodToolbox<2>::computePlusAndMinusSpatialDerivatives(
    patch_hierarchy,
    LSMLIB::ENO,
    1, // first-order
    grad_phi_plus_handle,
    grad_phi_minus_handle,
    phi_scratch_handle);
*/

  /*
   * Compute and output volume integral of the x-component
   * of the velocity.  For the  velocity field 0, this integral
   * is the volume of the sphere bounded by the zero level set.
   */
  LSMLIB_REAL volume_integral = LevelSetMethodToolbox<2>::computeVolumeIntegral(
    patch_hierarchy,
    velocity_handle, 
    phi_handle,
    control_volume_handle,
    -1); // -1 indicates that integral should be over region with phi <= 0
  pout << "Volume integral of x-component of velocity = " 
       << volume_integral << endl;
 
  /*
   * Compute and output surface integral of the x-component
   * of the velocity.  For the velocity field 0, this integral
   * is the surface area of the sphere defined by the zero level set.
   */
  LSMLIB_REAL surface_integral = LevelSetMethodToolbox<2>::computeSurfaceIntegral(
    patch_hierarchy,
    velocity_handle, 
    phi_handle,
    grad_phi_plus_handle,
    control_volume_handle);
  pout << "Surface integral of x-component of velocity = " 
       << surface_integral << endl;
 
  /*
   * At conclusion of simulation, deallocate objects and free memory.
   */
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    Pointer< PatchLevel<2> > level =
      patch_hierarchy->getPatchLevel(ln);
    level->deallocatePatchData( phi_scratch_handle );
    level->deallocatePatchData( grad_phi_plus_handle );
    level->deallocatePatchData( grad_phi_minus_handle );
  }
  delete patch_module;
  delete velocity_field_module;

  SAMRAIManager::shutdown();
  tbox::MPI::finalize();

  return(0);
}
