/*
 * File:        main.cc
 * Description: Example program for LSMLIB Parallel Package
 */

// Standard headers
#include <iosfwd>
#include <iostream>
#include <stdlib.h>
#include <string>

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

// LSMLIB headers
#include "LSMLIB_config.h" // LSMLIB configuration header must be included
                           // before any other LSMLIB header files
#include "LevelSetMethodAlgorithm.h"
#include "LevelSetMethodToolbox.h"

// Application headers
#include "VelocityFieldModule.h"
#include "PatchModule.h"

// Class/type declarations

// namespaces
using namespace std;
using namespace SAMRAI;
using namespace LSMLIB;

int main(int argc, char *argv[])
{

    // --- Initialize MPI and SAMRAI, enable logging, and process command line

    tbox::SAMRAI_MPI::init(&argc, &argv);
    tbox::SAMRAIManager::initialize();
    tbox::SAMRAIManager::startup();
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

    string input_filename;
    string restart_read_dirname;
    int restore_num = 0;

    bool is_from_restart = false;

    if ( (argc != 2) && (argc != 4) ) {
        tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
             << "\n"
             << "<restart dir> <restore number> [options]\n"
             << "  options:\n"
             << "  none at this time"
             << endl;
        tbox::SAMRAI_MPI::abort();
        return (-1);
    } else {
        input_filename = argv[1];
        if (argc == 4) {
            restart_read_dirname = argv[2];
            restore_num = atoi(argv[3]);
            is_from_restart = true;
        }
    }

    // --- Load parameters from input file

    boost::shared_ptr<tbox::InputDatabase> input_db =
        boost::shared_ptr<tbox::InputDatabase>(
            new tbox::InputDatabase("input_db"));
    tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

    // Process "Main" section of the input database.
    boost::shared_ptr<tbox::Database> main_db = input_db->getDatabase("Main");

    // Get dimension
    const tbox::Dimension dim(
        static_cast<unsigned short>(main_db->getInteger("dim")));

    // Get the base name for all name strings in this program
    string base_name = "unnamed";
    base_name = main_db->getStringWithDefault("base_name", base_name);

    // --- Set up restart manager

    int restart_interval = 0;
    if (main_db->keyExists("restart_interval")) {
        restart_interval = main_db->getInteger("restart_interval");
    }
    string restart_write_dirname = base_name + ".restart";
    const bool write_restart = (restart_interval > 0);

    tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

    // If run is from restart, open the restart file
    if (is_from_restart) {
        restart_manager->openRestartFile(restart_read_dirname, restore_num,
                                         mpi.getSize());
    }

    // --- Start logging

    const string log_file_name = base_name + ".log";
    bool log_all_nodes = false;
    log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
                                                log_all_nodes);
    if (log_all_nodes) {
        tbox::PIO::logAllNodes(log_file_name);
    } else {
        tbox::PIO::logOnlyNodeZero(log_file_name);
    }

    // Log the command-line args
    tbox::plog << "input_filename = " << input_filename << endl;
    tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
    tbox::plog << "restore_num = " << restore_num << endl;

    // --- Create major algorithm and data objects

    boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry =
        boost::shared_ptr<geom::CartesianGridGeometry>(
            new geom::CartesianGridGeometry(
                dim,
                base_name+"::CartesianGeometry",
                input_db->getDatabase("CartesianGeometry")));
    tbox::plog << "CartesianGridGeometry:" << endl;
    grid_geometry->printClassData(tbox::plog);

    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy =
        boost::shared_ptr<hier::PatchHierarchy>(
            new hier::PatchHierarchy(base_name+"::PatchHierarchy",
                                     grid_geometry));

    VelocityFieldModule* velocity_field_module = new VelocityFieldModule(
        input_db->getDatabase("VelocityFieldModule"),
        patch_hierarchy,
        grid_geometry,
        base_name+"::VelocityFieldModule");
    tbox::plog << "VelocityFieldModule:" << endl;
    velocity_field_module->printClassData(tbox::plog);

    PatchModule* patch_module = new PatchModule(
        input_db->getDatabase("PatchModule"),
        patch_hierarchy->getDim(),
        base_name+"::PatchModule");
    tbox::plog << "PatchModule:" << endl;
    patch_module->printClassData(tbox::plog);

    int num_level_set_fcn_components = 1;
    int codimension = 1;
    boost::shared_ptr<LevelSetMethodAlgorithm> lsm_algorithm =
        boost::shared_ptr<LevelSetMethodAlgorithm>(
            new LevelSetMethodAlgorithm(
                input_db->getDatabase("LevelSetMethodAlgorithm"),
                patch_hierarchy,
                patch_module,
                velocity_field_module,
                num_level_set_fcn_components,
                codimension,
                base_name+"::LevelSetMethodAlgorithm"));
    tbox::plog << "LevelSetMethodAlgorithm:" << endl;
    lsm_algorithm->printClassData(tbox::plog);

    // Initialize field extension algorithm
    int phi_handle = lsm_algorithm->getPhiPatchDataHandle();
    int control_volume_handle =
        lsm_algorithm->getControlVolumePatchDataHandle();
    velocity_field_module->initializeFieldExtensionAlgorithm(
        input_db->getDatabase("VelocityFieldModule")
                ->getDatabase("FieldExtensionAlgorithm"),
        phi_handle, control_volume_handle);

    // Emit contents of input database and variable database to log file.
    tbox::plog << "\nCheck input data and variables before simulation:"
        << endl;
    tbox::plog << "Input database..." << endl;
    input_db->printClassData(tbox::plog);
    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

    // --- Set up visualization data writers

    bool use_visit = false;
    if (main_db->keyExists("use_visit")) {
        use_visit = main_db->getBool("use_visit");
    }

    // Set up viz write interval
    int viz_write_interval = -1;
    if (use_visit) {
        if (main_db->keyExists("viz_write_interval")) {
            viz_write_interval =
                main_db->getInteger("viz_write_interval");
        }
    }

    // Set up extra VisIt parameters
    int visit_number_procs_per_file = 1;
    if (use_visit) {
        if (main_db->keyExists("visit_number_procs_per_file")) {
            visit_number_procs_per_file =
                main_db->getInteger("visit_number_procs_per_file");
        }
    }

    // Set up Visit
    boost::shared_ptr<appu::VisItDataWriter> visit_data_writer = 0;
    if (use_visit) {
        string visit_data_dirname = base_name + ".visit";
        visit_data_writer = boost::shared_ptr<appu::VisItDataWriter>
           (new appu::VisItDataWriter(
                dim, "VisIt Writer", visit_data_dirname,
                visit_number_procs_per_file));

      // Register level set functions and velocity fields for plotting
      visit_data_writer->registerPlotQuantity(
          "phi", "SCALAR", phi_handle, 0, 1.0, "CELL");

      visit_data_writer->registerPlotQuantity(
          "velocity", "SCALAR",
          velocity_field_module->getNormalVelocityFieldPatchDataHandle(
              LEVEL_SET_FCN_TYPE::PHI, 0),
          0, 1.0, "CELL");
    }

    // --- Initialize level set method calculation

    lsm_algorithm->initializeLevelSetMethodCalculation();

    // Close restart file before starting main time-stepping loop
    restart_manager->closeRestartFile();

    // Set up loop variables
    int count = 0;
    int max_num_time_steps = main_db->getInteger("max_num_time_steps");
    LSMLIB_REAL dt = 0;
    LSMLIB_REAL current_time = lsm_algorithm->getCurrentTime();
    int cur_integrator_step = lsm_algorithm->numIntegrationStepsTaken();

    // Reinitialize level set functions
    lsm_algorithm->reinitializeLevelSetFunctions();

    // Get grad(phi) PatchData handle
    int grad_phi_handle = velocity_field_module->getGradPhiHandle();

    // Output initial conditions (if this run is not from restart)
    if ( write_restart && (!is_from_restart) ) {
        restart_manager->writeRestartFile(restart_write_dirname,
                                          cur_integrator_step);
    }

    // Write VisIt data for initial time step
    if ( use_visit && (!is_from_restart) ) {
        visit_data_writer->writePlotData(patch_hierarchy, cur_integrator_step,
                                         current_time);
    }

    // --- Main time loop

    while ( !lsm_algorithm->endTimeReached() &&
            ((max_num_time_steps <= 0) || (count < max_num_time_steps)) ) {

        // Get current time
        current_time = lsm_algorithm->getCurrentTime();

        tbox::pout << "++++++++++++++++++++++++++++++++++++++++++"
            << endl;
        tbox::pout << "  Time step (in current run): " << count << endl;
        tbox::pout << "  Integrator time step: " << cur_integrator_step
            << endl;
        tbox::pout << "  Current time:  " << current_time << endl;

        // Compute volume and surface area
        double volume =
            LevelSetMethodToolbox::computeVolumeOfRegionDefinedByZeroLevelSet(
                patch_hierarchy,
                phi_handle,
                control_volume_handle,
                -1); // compute volume of region with phi < 0
        double surface_area =
            LevelSetMethodToolbox::computeVolumeOfZeroLevelSet(
                patch_hierarchy,
                phi_handle,
                grad_phi_handle,
                control_volume_handle);

        if (dim.getValue() == 3) {
            tbox::pout << "  Current volume:  " << volume << endl;
        } else if (dim.getValue() == 2) {
            tbox::pout << "  Current area:  " << volume << endl;
        }
        if (count == 0) {
            tbox::pout << "  Current surface area: N/A" << endl;
        } else {
            if (dim.getValue() == 3) {
                tbox::pout << "  Current surface area:  "
                           << surface_area << endl;
            } else if (dim.getValue() == 2) {
                tbox::pout << "  Current perimeter:  "
                           << surface_area << endl;
            }
        }

        // Compute next time step
        dt = lsm_algorithm->computeStableDt();
        LSMLIB_REAL end_time = lsm_algorithm->getEndTime();
        if (end_time - current_time < dt) dt = end_time - current_time;
        tbox::pout << "  dt:  " << dt << endl;

        // Advance level set functions
        lsm_algorithm->advanceLevelSetFunctions(dt);

        // Add an extra line to output for aesthetic reasons
        tbox::pout << endl;

        // Output data for current time step if this is the
        // initial time step or if the next write interval has
        // been reached
        cur_integrator_step = lsm_algorithm->numIntegrationStepsTaken();

        // Write restart file
        if ( write_restart && (0==cur_integrator_step%restart_interval) ) {
            restart_manager->writeRestartFile(restart_write_dirname,
                                          cur_integrator_step);
        }

        // Write VisIt data
        if ( use_visit && (0==cur_integrator_step%viz_write_interval) ) {
            visit_data_writer->writePlotData(patch_hierarchy,
                                             cur_integrator_step,
                                             lsm_algorithm->getCurrentTime());
      }

      // Update count
      count++;
    }

    // Output information for final time step
    current_time = lsm_algorithm->getCurrentTime();
    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;
    tbox::pout << "  Final time step (in current run): " << count << endl;
    tbox::pout << "  Final integrator time step: " << cur_integrator_step
        << endl;
    tbox::pout << "  Current time:  " << current_time << endl;

    // Compute volume and surface area
    double volume =
        LevelSetMethodToolbox::computeVolumeOfRegionDefinedByZeroLevelSet(
            patch_hierarchy,
            phi_handle,
            control_volume_handle,
            -1); // compute volume of region with phi < 0
    double surface_area =
        LevelSetMethodToolbox::computeVolumeOfZeroLevelSet(
            patch_hierarchy,
            phi_handle,
            grad_phi_handle,
            control_volume_handle);

    if (dim.getValue() == 3) {
        tbox::pout << "  Current volume:  " << volume << endl;
        tbox::pout << "  Current surface area:  " << surface_area << endl;
    } else if (dim.getValue() == 2) {
        tbox::pout << "  Current area:  " << volume << endl;
        tbox::pout << "  Current perimeter:  " << surface_area << endl;
    }

    tbox::pout << endl;
    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;

    // Write restart file for final time step
    if ( write_restart && (0!=cur_integrator_step%restart_interval) ) {
        restart_manager->writeRestartFile(restart_write_dirname,
                                          cur_integrator_step);
    }

    // Write VisIt data for final time step
    if ( use_visit && (0!=cur_integrator_step%viz_write_interval) ) {
        visit_data_writer->writePlotData(patch_hierarchy, cur_integrator_step,
                                         lsm_algorithm->getCurrentTime());

    }

    // --- At conclusion of simulation, clean up

    delete patch_module;
    delete velocity_field_module;

    tbox::SAMRAIManager::shutdown();
    tbox::SAMRAIManager::finalize();
    tbox::SAMRAI_MPI::finalize();

    return(0);
}
