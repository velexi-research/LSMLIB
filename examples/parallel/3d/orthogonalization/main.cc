/*
 * File:        main.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: example program for the OrthogonalizationAlgorithm
 */

/************************************************************************
 *
 * This program example the OrthogonalizationAlgorithm.
 *
 ************************************************************************/

// SAMRAI Headers
#include "SAMRAI_config.h"

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
#include "tbox/SAMRAIManager.h"
#include "tbox/Utilities.h"
#include "VisItDataWriter.h"

// Headers for level set method
// LevelSetMethod configuration header must be included
// before any other LevelSetMethod header files
#include "LSMLIB_config.h"
#include "LevelSetMethodAlgorithm.h"
#include "VelocityFieldModule.h"
#include "PatchModule.h"


// namespaces
using namespace std;
using namespace SAMRAI;
using namespace appu;
using namespace geom;
using namespace hier;
using namespace tbox;


// CONSTANTS
#define DIM   (3)

int main(int argc, char *argv[])
{

  /*
   * Initialize MPI and SAMRAI, enable logging, and process command line.
   */
  tbox::MPI::init(&argc, &argv);
  tbox::MPI::initialize();
  SAMRAIManager::startup();

 
  string input_filename;

  // process command-line arguments
  if (argc != 2) {
    pout << "USAGE:  " << argv[0] << " <input filename> "
         << "\n"
         << "  options:\n"
         << "  none at this time"
         << endl;
    tbox::MPI::abort();
    return (-1);
  } else {
    input_filename = argv[1];
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

  // log the command-line args
  plog << "input_filename = " << input_filename << endl;

  /*
   * Create PatchHierarchy and CartesianGridGeometry 
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
      base_name+"::PatchModule");
  plog << "PatchModule:" << endl;
  patch_module->printClassData(plog);

  int num_level_set_fcn_components = 1;
  int codimension = 2;
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
   * Set up VisIt data writer
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
  int phi_handle = lsm_algorithm->getPhiPatchDataHandle();
  int psi_handle = lsm_algorithm->getPsiPatchDataHandle();

  // set up VisIt data writer
  Pointer<VisItDataWriter<DIM> > visit_data_writer = 0;
  if ( use_visit ) {
    string visit_data_dirname = base_name + ".visit";
    visit_data_writer = new VisItDataWriter<DIM>("VisIt Writer",
                                                  visit_data_dirname,
                                                  visit_number_procs_per_file);

    // register phi and psi for plotting
    visit_data_writer->registerPlotQuantity(
      "phi", "SCALAR", phi_handle, 0, 1.0, "CELL");
    visit_data_writer->registerPlotQuantity(
      "psi", "SCALAR", psi_handle, 0, 1.0, "CELL");

  }  

  /*
   * Initialize level set method calculation
   */
  lsm_algorithm->initializeLevelSetMethodCalculation();

  /*
   * output phi and psi before orthogonalization
   */
  pout << "=============================" << endl;
  pout << "Writing initial level set functions to file...";

  // write out stress field to VisIt data format
  if ( use_visit ) {
    visit_data_writer->writePlotData(patch_hierarchy, 0, 0);
  }

  pout << "done" << endl;
  pout << "=============================" << endl;


  /*
   * Orthogonalize level set functions
   */
  pout << "=============================" << endl;
  pout << "Orthogonalizing level set functions..." << endl << endl;

  lsm_algorithm->orthogonalizeLevelSetFunctions(LSMLIB::PHI);

  pout << endl << "=============================" << endl;

  /*
   * output results
   */
  pout << "=============================" << endl;
  pout << "Writing results to file...";

  // write out stress field to VisIt data format
  if ( use_visit ) {
    visit_data_writer->writePlotData(patch_hierarchy, 1, 1);
  }

  pout << "done" << endl;
  pout << "=============================" << endl;

  /*
   * At conclusion of simulation, deallocate objects.
   */
  delete patch_module;
  delete velocity_field_module;

  SAMRAIManager::shutdown();
  tbox::MPI::finalize();

  return(0);
}
