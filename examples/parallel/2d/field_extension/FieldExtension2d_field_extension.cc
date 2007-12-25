/*
 * File:        FieldExtension2d_field_extension.cc
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date: 2006/05/18 01:09:35 $
 * Description: 2D example program for Level Set Method Classes
 */

/*************************************************************************
 *
 * This example program demonstrates how to use the LSMLIB C++ classes
 * to for a 2D problem to extend a field variable off of the zero
 * level set.
 *
 **************************************************************************/


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
#include "VisItDataWriter.h"

// Headers for level set method
// LevelSetMethod configuration header must be included
// before any other LevelSetMethod header files
#include "LSMLIB_config.h"
#include "LevelSetMethodAlgorithm.h"
#include "FieldExtensionAlgorithm.h"
#include "FieldExtension2d_VelocityFieldModule.h"
#include "FieldExtension2d_PatchModule.h"


// namespaces
using namespace std;
using namespace SAMRAI;
using namespace appu;
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

  if ( argc != 2 ) {
    pout << "USAGE:  " << argv[0] << " <input filename> " << endl;
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

  FieldExtension2d_VelocityFieldModule* fieldextension2d_velocity_field_module = 
    new FieldExtension2d_VelocityFieldModule( 
      input_db->getDatabase("FieldExtension2d_VelocityFieldModule"),
      patch_hierarchy,
      grid_geometry,
      base_name+"::FieldExtension2d_VelocityFieldModule");
  plog << "FieldExtension2d_VelocityFieldModule:" << endl;
  fieldextension2d_velocity_field_module->printClassData(plog);

  FieldExtension2d_PatchModule* fieldextension2d_patch_module = 
    new FieldExtension2d_PatchModule(
      input_db->getDatabase("FieldExtension2d_PatchModule"),
      base_name+"::FieldExtension2d_PatchModule");
  plog << "FieldExtension2d_PatchModule:" << endl;
  fieldextension2d_patch_module->printClassData(plog);

  int num_level_set_fcn_components = 1;
  int codimension = 1;
  Pointer< LevelSetMethodAlgorithm<2> > lsm_algorithm = 
    new LevelSetMethodAlgorithm<2>( 
      input_db->getDatabase("LevelSetMethodAlgorithm"),
      patch_hierarchy,
      fieldextension2d_patch_module,
      fieldextension2d_velocity_field_module,
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


  /*
   * get PatchData handles
   */
  int phi_handle = lsm_algorithm->getPhiPatchDataHandle();
  int velocity_handle = fieldextension2d_velocity_field_module
      ->getExternalVelocityFieldPatchDataHandle(0);


  /*
   * Set up VisIt data writer
   */
  Pointer<VisItDataWriter<2> > visit_data_writer = 0;
  bool use_visit = false;
  int visit_number_procs_per_file = 1;
  if (main_db->keyExists("use_visit")) {
    use_visit = main_db->getBool("use_visit");
    if (main_db->keyExists("visit_number_procs_per_file")) {
      visit_number_procs_per_file =
        main_db->getInteger("visit_number_procs_per_file");
    }

    string visit_data_dirname = base_name + ".visit";
    visit_data_writer = new VisItDataWriter<2>("VisIt Writer",
                                               visit_data_dirname,
                                               visit_number_procs_per_file);

    // register level set functions and velocity fields for plotting
    visit_data_writer->registerPlotQuantity(
      "phi", "SCALAR",
      phi_handle, 0, 1.0, "CELL");

    visit_data_writer->registerPlotQuantity(
      "velocity-x", "SCALAR",
      velocity_handle, 0, 1.0, "CELL");
    visit_data_writer->registerPlotQuantity(
      "velocity-y", "SCALAR",
      velocity_handle, 1, 1.0, "CELL");
  }


  /*
   * Create FieldExtensionAlgorithm object
   */
  Pointer< FieldExtensionAlgorithm<2> > field_ext_alg = 
    lsm_algorithm->getFieldExtensionAlgorithm(
      input_db->getDatabase("FieldExtensionAlgorithm"),
      velocity_handle, 
      LSMLIB::PHI);
 
  /*
   * Initialize level set method calculation
   */ 
  lsm_algorithm->initializeLevelSetMethodCalculation();

  /*
   * Extend the x and y components of the velocity off of the 
   * zero level set.
   */
  field_ext_alg->computeExtensionField();

  /* 
   * write results to VisIt data 
   */
  if ( use_visit ) {
    visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
  }

  /*
   * At conclusion of simulation, deallocate objects and free memory.
   */
  delete fieldextension2d_patch_module;
  delete fieldextension2d_velocity_field_module;

  SAMRAIManager::shutdown();
  tbox::MPI::finalize();

  return(0);
}
