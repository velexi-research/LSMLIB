/*
 * File:        main.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.7 $
 * Modified:    $Date: 2006/08/01 16:13:43 $
 * Description: 3D example program for Level Set Method Classes
 */

/************************************************************************
 *
 * This example program demonstrates how to use the LSMLIB C++ classes
 * to for a 3D problem where the motion of the level sets is determined
 * by an external velocity field (defined throughout the entire
 * computational domain) is specified.
 *
 ************************************************************************/

// SAMRAI Configuration 
#include "SAMRAI_config.h"

/* 
 * Headers for basic SAMRAI objects
 */

// variables and variable management
#include "CellVariable.h"
#include "FaceVariable.h"
#include "VariableDatabase.h"

// geometry and patch hierarchy
#include "CartesianGridGeometry.h" 
#include "PatchHierarchy.h"

// basic SAMRAI classes
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
#include "VelocityFieldModule.h"
#include "PatchModule.h"

// Header for data writers
#include "VisItDataWriter.h"


// namespaces
using namespace std;
using namespace SAMRAI;
using namespace appu;
using namespace geom;
using namespace hier;
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

  int restart_interval = 0;
  if (main_db->keyExists("restart_interval")) {
    restart_interval = main_db->getInteger("restart_interval");
  } 
  string restart_write_dirname = base_name + ".restart";
  const bool write_restart = (restart_interval > 0); 

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

  Pointer< CartesianGridGeometry<3> > grid_geometry =
    new CartesianGridGeometry<3>(
      base_name+"::CartesianGeometry",
      input_db->getDatabase("CartesianGeometry"));
  plog << "CartesianGridGeometry:" << endl;
  grid_geometry->printClassData(plog);

  Pointer< PatchHierarchy<3> > patch_hierarchy =
    new PatchHierarchy<3>(base_name+"::PatchHierarchy",
                             grid_geometry);

  VelocityFieldModule* velocity_field_module = new VelocityFieldModule( 
      input_db->getDatabase("VelocityFieldModule"),
      patch_hierarchy,
      grid_geometry,
      base_name+"::VelocityFieldModule");
  plog << "VelocityFieldModule:" << endl;
  velocity_field_module->printClassData(plog);

  PatchModule* patch_module = new PatchModule(
      input_db->getDatabase("PatchModule"),
      base_name+"::PatchModule");
  plog << "PatchModule:" << endl;
  patch_module->printClassData(plog);

  int num_level_set_fcn_components = 1;
  int codimension = 1;
  Pointer< LevelSetMethodAlgorithm<3> > lsm_algorithm = 
    new LevelSetMethodAlgorithm<3>( 
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
  VariableDatabase<3>::getDatabase()->printClassData(plog);


  /*
   * Set up visualization data writers
   */

  bool use_visit = false;
  if (main_db->keyExists("use_visit")) {
    use_visit = main_db->getBool("use_visit");
  }

  // set up viz write interval  
  int viz_write_interval = -1;
  if (use_visit) {
    if (main_db->keyExists("viz_write_interval")) {
      viz_write_interval =
        main_db->getInteger("viz_write_interval");
    }
  }

  // set up extra VisIt parameters
  int visit_number_procs_per_file = 1;
  if (use_visit) {
    if (main_db->keyExists("visit_number_procs_per_file")) {
      visit_number_procs_per_file =
        main_db->getInteger("visit_number_procs_per_file");
    } 
  }

  // get PatchData handles
  int phi_patch_data_handle = lsm_algorithm->getPhiPatchDataHandle();
  int psi_patch_data_handle = lsm_algorithm->getPsiPatchDataHandle();
  int velocity_patch_data_handle = 
    velocity_field_module
      ->getExternalVelocityFieldPatchDataHandle(0);

  Pointer<VisItDataWriter<3> > visit_data_writer = 0;
  if (use_visit) {
    string visit_data_dirname = base_name + ".visit";
    visit_data_writer = new VisItDataWriter<3>("VisIt Writer",
                                               visit_data_dirname,
                                               visit_number_procs_per_file);

    // register level set functions and velocity fields for plotting
    visit_data_writer->registerPlotQuantity(
      "phi", "SCALAR", 
      phi_patch_data_handle, 0, 1.0, "CELL");
  
    if (psi_patch_data_handle >= 0) {
      visit_data_writer->registerPlotQuantity(
        "psi", "SCALAR", 
        psi_patch_data_handle, 0, 1.0, "CELL");
    }
  
    visit_data_writer->registerPlotQuantity(
      "velocity", "VECTOR", 
      velocity_patch_data_handle, 0, 1.0, "CELL");
  }  

  /*
   * Initialize level set method calculation
   */ 
  lsm_algorithm->initializeLevelSetMethodCalculation();


  /*
   * Close restart file before starting main time-stepping loop.
   */ 
  restart_manager->closeRestartFile();


  /* 
   * Set up loop variables
   */
  int count = 0;
  int max_num_time_steps = main_db->getInteger("max_num_time_steps");
  LSMLIB_REAL dt = 0;
  LSMLIB_REAL current_time = lsm_algorithm->getCurrentTime();
  int cur_integrator_step = lsm_algorithm->numIntegrationStepsTaken();

  /* 
   * Output initial conditions (if this run is not from restart).
   */
  if ( write_restart && (!is_from_restart) ) {
    restart_manager->writeRestartFile(restart_write_dirname,
                                      cur_integrator_step);
  }

  // write VisIt data for first time step
  if ( use_visit && (!is_from_restart) ) {
    visit_data_writer->writePlotData(patch_hierarchy, cur_integrator_step,
                                     current_time);
  }


  /*
   * Main time loop
   */ 
  while ( !lsm_algorithm->endTimeReached() && 
          ((max_num_time_steps <= 0) || (count < max_num_time_steps)) ) {

    pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;
    pout << "  Time step (in current run): " << count << endl;
    pout << "  Integrator time step: " << cur_integrator_step << endl;
    pout << "  Current time:  " << current_time << endl;

    // compute next time step
    dt = lsm_algorithm->computeStableDt();
    LSMLIB_REAL end_time = lsm_algorithm->getEndTime(); 
    if (end_time - current_time < dt) dt = end_time - current_time;
    pout << "  dt:  " << dt << endl;

    // advance level set functions
    lsm_algorithm->advanceLevelSetFunctions(dt);
 
    // add an extra line to output for aesthetic reasons
    pout << endl;

    /* 
     * output data for current time step if this is the
     * initial time step or if the next write interval has
     * been reached
     */
    cur_integrator_step = lsm_algorithm->numIntegrationStepsTaken();

    // write restart file
    if ( write_restart && (0==cur_integrator_step%restart_interval) ) {
      restart_manager->writeRestartFile(restart_write_dirname,
                                        cur_integrator_step);
    }

    // write VisIt data
    if ( use_visit && (0==cur_integrator_step%viz_write_interval) ) {
      visit_data_writer->writePlotData(patch_hierarchy, cur_integrator_step, 
                                       lsm_algorithm->getCurrentTime());
    }

    // update counter and current time
    count++; 
    current_time = lsm_algorithm->getCurrentTime();

  }

  // output information for final time step 
  // (if it hasn't already been output)
  current_time = lsm_algorithm->getCurrentTime();
  pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;
  pout << "  Final time step (in current run): " << count << endl;
  pout << "  Final integrator time step: " << cur_integrator_step << endl;
  pout << "  Current time:  " << current_time << endl;
  pout << endl;
  pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;

  // write restart file for final time step
  if ( write_restart && (0!=cur_integrator_step%restart_interval) ) {
    restart_manager->writeRestartFile(restart_write_dirname,
                                      cur_integrator_step);
  }

  // write VisIt data for final time step
  if ( use_visit && (0!=cur_integrator_step%viz_write_interval) ) {
    visit_data_writer->writePlotData(patch_hierarchy, cur_integrator_step,
                                     lsm_algorithm->getCurrentTime());
                                       
  }


  /*
   * At conclusion of simulation, deallocate objects.
   */
  delete patch_module;
  delete velocity_field_module;

  SAMRAIManager::shutdown();
  tbox::MPI::finalize();

  return(0);
}
