/*
 * File:        VelocityFieldModule.cc
 * Description: Implementation of class that computes the velocity field
 *              for the level set method
 */

// Class header
#include "VelocityFieldModule.h"

// Standard headers
#include <assert.h>
#include <cmath>
#include <sstream>

// Boost headers
// IWYU pragma: no_include <boost/smart_ptr/detail/operator_bool.hpp>
#include <boost/smart_ptr/make_shared_object.hpp>

// SAMRAI headers
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Utilities.h"

// LSMLIB headers
#include "FieldExtensionAlgorithm.h"

extern "C" {
  #include "lsm_geometry2d.h"
  #include "lsm_geometry3d.h"
  #include "lsm_spatial_derivatives2d.h"
  #include "lsm_spatial_derivatives3d.h"
}

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchData; } }
namespace SAMRAI { namespace hier { class PatchGeometry; } }
namespace SAMRAI { namespace hier { class VariableContext; } }

// SAMRAI namespaces
using namespace pdat;


/* Constructor */
VelocityFieldModule::VelocityFieldModule(
  boost::shared_ptr<tbox::Database> input_db,
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geom,
  const string& object_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(input_db);
  assert(patch_hierarchy);
  assert(grid_geom);
  assert(!object_name.empty());
#endif

  // set object name, patch hierarchy, and  grid geometry
  d_object_name = object_name;
  d_patch_hierarchy = patch_hierarchy;
  d_grid_geometry = grid_geom;

  // read in input data
  getFromInput(input_db);

  // Allocate velocity variable
  const tbox::Dimension dim = patch_hierarchy->getDim();
  const int num_dims = dim.getValue();
  boost::shared_ptr< pdat::CellVariable<LSMLIB_REAL> > velocity =
      boost::shared_ptr< pdat::CellVariable<LSMLIB_REAL> >(
          new pdat::CellVariable<LSMLIB_REAL>(dim, "velocity field", 1));

  boost::shared_ptr< pdat::CellVariable<LSMLIB_REAL> > grad_phi =
      boost::shared_ptr< pdat::CellVariable<LSMLIB_REAL> >(
          new pdat::CellVariable<LSMLIB_REAL>(dim, "grad phi", num_dims));

  boost::shared_ptr< pdat::CellVariable<LSMLIB_REAL> > hessian_phi =
      boost::shared_ptr< pdat::CellVariable<LSMLIB_REAL> >(
          new pdat::CellVariable<LSMLIB_REAL>(dim, "hessian phi",
                                              num_dims * (num_dims + 1) / 2));

  boost::shared_ptr< pdat::CellVariable<LSMLIB_REAL> > kappa =
      boost::shared_ptr< pdat::CellVariable<LSMLIB_REAL> >(
          new pdat::CellVariable<LSMLIB_REAL>(dim, "mean curvature", 1));

  // Register velocity variable with VariableDatabase.
  hier::VariableDatabase *vdb = hier::VariableDatabase::getDatabase();
  boost::shared_ptr<hier::VariableContext> cur_ctxt = vdb->getContext("CURRENT");
  d_grad_phi_handle = vdb->registerVariableAndContext(
    grad_phi, cur_ctxt, hier::IntVector(dim, 0));

  d_hessian_phi_handle = vdb->registerVariableAndContext(
    hessian_phi, cur_ctxt, hier::IntVector(dim, 0));

  d_kappa_handle = vdb->registerVariableAndContext(
    kappa, cur_ctxt, hier::IntVector(dim, 0));

  d_velocity_handle = vdb->registerVariableAndContext(
    velocity, cur_ctxt, hier::IntVector(dim, 0));
}

/* initializeFieldExtensionAlgorithm() */
void VelocityFieldModule::initializeFieldExtensionAlgorithm(
    boost::shared_ptr<tbox::Database> input_db,
    int phi_handle, int control_volume_handle)
{
    if (d_use_field_extension) {
        tbox::Dimension dim = d_patch_hierarchy->getDim();
        d_field_extension_alg = boost::shared_ptr<FieldExtensionAlgorithm>(
            new FieldExtensionAlgorithm(
                input_db,
                d_patch_hierarchy,
                d_velocity_handle,
                phi_handle,
                control_volume_handle,
                hier::IntVector(dim, 0)));
    }
}

/* computeStableDt() */
LSMLIB_REAL VelocityFieldModule::computeStableDt()
{
    // get problem dimension and grid spacing
    const tbox::Dimension dim = d_patch_hierarchy->getDim();
    const int num_dims = dim.getValue();
    const double* dx = d_grid_geometry->getDx();

    // compute contribution to stability factor from constant normal velocity
    // Note: the calculation used overestimates the contributions in 2D and 3D
    //       by bounding H_i by 1.
    //
    //     - 2D: a * max{|H_1|/dx + |H_2|/dy}
    //
    //     - 3D: a * max{|H_1|/dx + |H_2|/dy + |H_3|/dz}
    double const_normal_velocity_term = 0.0;
    if (d_a != 0) {
      for (int i = 0; i < num_dims; i++) {
          const_normal_velocity_term += 1.0 / dx[i];
      }
    }
    const_normal_velocity_term *= abs(d_a);

    // compute contribution to stability factor from curvature term
    double curvature_term = 0.0;
    for (int i = 0; i < num_dims; i++) {
        curvature_term += 1.0 / (dx[i] * dx[i]);
    }
    curvature_term *= 2 * d_b;

    double dt = 1.0 / (const_normal_velocity_term + curvature_term);

    return dt;
}

/* computeVelocityField() */
void VelocityFieldModule::computeVelocityField(
  const LSMLIB_REAL time,
  const int phi_handle,
  const int psi_handle,
  const int component)
{
  (void) psi_handle; // psi is meaningless for co-dimension one problems
  (void) component;  // component is not used because this example problem
                     // only has one component for level set function

  // update the current time
  d_current_time = time;

  // set velocity on all levels of hierarchy
  const int finest_level = d_patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    boost::shared_ptr< hier::PatchLevel > level =
        d_patch_hierarchy->getPatchLevel(ln);
    computeVelocityFieldOnLevel(level, time, phi_handle);

  } // end loop over hierarchy

  // Extend velocity field off of interface
  if (d_use_field_extension) {
      const tbox::Dimension dim = d_patch_hierarchy->getDim();
      const hier::IntVector bc_type(dim, -1);
      const int phi_component = 0;
      d_field_extension_alg->computeExtensionField(phi_component,
        bc_type, bc_type, bc_type, bc_type);
  }
}


/* initializeLevelData() */
void VelocityFieldModule::initializeLevelData (
  const boost::shared_ptr<hier::PatchHierarchy> hierarchy ,
  const int level_number ,
  const LSMLIB_REAL init_data_time ,
  const int phi_handle,
  const int psi_handle,
  const bool can_be_refined ,
  const bool initial_time ,
  const boost::shared_ptr< hier::PatchLevel > old_level,
  const bool allocate_data)
{
  (void) psi_handle;  // psi is meaningless for co-dimension one problems

  boost::shared_ptr< hier::PatchLevel > level =
      hierarchy->getPatchLevel(level_number);
  if (allocate_data) {
    level->allocatePatchData(d_velocity_handle);
    level->allocatePatchData(d_grad_phi_handle);
    level->allocatePatchData(d_hessian_phi_handle);
    level->allocatePatchData(d_kappa_handle);
  }

  /*
   * Initialize data on all patches in the level.
   */
  computeVelocityFieldOnLevel(level,init_data_time,phi_handle);

}

/* computeVelocityFieldOnLevel() */
void VelocityFieldModule::computeVelocityFieldOnLevel(
    const boost::shared_ptr< hier::PatchLevel > level,
    const LSMLIB_REAL time,
    const int phi_handle)
{
    // get problem dimension
    const tbox::Dimension dim = d_patch_hierarchy->getDim();
    const int num_dims = dim.getValue();

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end();
         pi++) {

        // loop over patches
        boost::shared_ptr<hier::Patch> patch = *pi;
        if (!patch) {
            TBOX_ERROR(d_object_name
                       << ": Cannot find patch. Null patch pointer.");
        }

        boost::shared_ptr< CellData<LSMLIB_REAL> > phi_data =
            BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData(phi_handle));

        boost::shared_ptr< CellData<LSMLIB_REAL> > grad_phi_data =
            BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData(d_grad_phi_handle));

        boost::shared_ptr< CellData<LSMLIB_REAL> > hessian_phi_data =
            BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData(d_hessian_phi_handle));

        boost::shared_ptr< CellData<LSMLIB_REAL> > kappa_data =
            BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData(d_kappa_handle));

        boost::shared_ptr< CellData<LSMLIB_REAL> > velocity_data =
            BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData(d_velocity_handle));

        boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch->getPatchGeometry());

#ifdef LSMLIB_DOUBLE_PRECISION
        const double* dx = patch_geom->getDx();
        const double* x_lower = patch_geom->getXLower();
#else
        const double* dx_double = patch_geom->getDx();
        const double* x_lower_double = patch_geom->getXLower();
        float *dx = new float[num_dims];
        float *x_lower = new float[num_dims];
        for (int i = 0; i < num_dims; i++) {
            dx[i] = dx_double[i];
            x_lower[i] = x_lower_double[i];
        }
#endif

        hier::Box phi_ghostbox = phi_data->getGhostBox();
        const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
        const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

        hier::Box grad_phi_ghostbox = grad_phi_data->getGhostBox();
        const hier::IntVector grad_phi_ghostbox_lower =
            grad_phi_ghostbox.lower();
        const hier::IntVector grad_phi_ghostbox_upper =
            grad_phi_ghostbox.upper();

        hier::Box hessian_phi_ghostbox = hessian_phi_data->getGhostBox();
        const hier::IntVector hessian_phi_ghostbox_lower =
            hessian_phi_ghostbox.lower();
        const hier::IntVector hessian_phi_ghostbox_upper =
            hessian_phi_ghostbox.upper();

        hier::Box kappa_ghostbox = kappa_data->getGhostBox();
        const hier::IntVector kappa_ghostbox_lower = kappa_ghostbox.lower();
        const hier::IntVector kappa_ghostbox_upper = kappa_ghostbox.upper();

        hier::Box vel_ghostbox = velocity_data->getGhostBox();
        const hier::IntVector vel_ghostbox_lower = vel_ghostbox.lower();
        const hier::IntVector vel_ghostbox_upper = vel_ghostbox.upper();

        hier::Box fillbox = velocity_data->getBox();
        const hier::IntVector fillbox_lower = fillbox.lower();
        const hier::IntVector fillbox_upper = fillbox.upper();

        if (num_dims == 3) {

            // get data pointers
            LSMLIB_REAL* phi = phi_data->getPointer(0);

            LSMLIB_REAL* phi_x = grad_phi_data->getPointer(0);
            LSMLIB_REAL* phi_y = grad_phi_data->getPointer(1);
            LSMLIB_REAL* phi_z = grad_phi_data->getPointer(2);

            LSMLIB_REAL* phi_xx = hessian_phi_data->getPointer(0);
            LSMLIB_REAL* phi_xy = hessian_phi_data->getPointer(1);
            LSMLIB_REAL* phi_xz = hessian_phi_data->getPointer(2);
            LSMLIB_REAL* phi_yy = hessian_phi_data->getPointer(3);
            LSMLIB_REAL* phi_yz = hessian_phi_data->getPointer(4);
            LSMLIB_REAL* phi_zz = hessian_phi_data->getPointer(5);

            LSMLIB_REAL* kappa = kappa_data->getPointer(0);

            // compute gradient of phi
            LSM3D_CENTRAL_GRAD_ORDER2(
                phi_x,
                phi_y,
                phi_z,
                &grad_phi_ghostbox_lower[0],
                &grad_phi_ghostbox_upper[0],
                &grad_phi_ghostbox_lower[1],
                &grad_phi_ghostbox_upper[1],
                &grad_phi_ghostbox_lower[2],
                &grad_phi_ghostbox_upper[2],
                phi,
                &phi_ghostbox_lower[0],
                &phi_ghostbox_upper[0],
                &phi_ghostbox_lower[1],
                &phi_ghostbox_upper[1],
                &phi_ghostbox_lower[2],
                &phi_ghostbox_upper[2],
                &fillbox_lower[0],
                &fillbox_upper[0],
                &fillbox_lower[1],
                &fillbox_upper[1],
                &fillbox_lower[2],
                &fillbox_upper[2],
                &dx[0],
                &dx[1],
                &dx[2]);

            // compute Hessian of phi
            LSM3D_CENTRAL_HESSIAN(
                phi_xx,
                phi_xy,
                phi_xz,
                phi_yy,
                phi_yz,
                phi_zz,
                &hessian_phi_ghostbox_lower[0],
                &hessian_phi_ghostbox_upper[0],
                &hessian_phi_ghostbox_lower[1],
                &hessian_phi_ghostbox_upper[1],
                &hessian_phi_ghostbox_lower[2],
                &hessian_phi_ghostbox_upper[2],
                phi,
                &phi_ghostbox_lower[0],
                &phi_ghostbox_upper[0],
                &phi_ghostbox_lower[1],
                &phi_ghostbox_upper[1],
                &phi_ghostbox_lower[2],
                &phi_ghostbox_upper[2],
                &fillbox_lower[0],
                &fillbox_upper[0],
                &fillbox_lower[1],
                &fillbox_upper[1],
                &fillbox_lower[2],
                &fillbox_upper[2],
                &dx[0],
                &dx[1],
                &dx[2]);

            // compute mean curvature
            LSM3D_COMPUTE_MEAN_CURVATURE(
                kappa,
                &kappa_ghostbox_lower[0],
                &kappa_ghostbox_upper[0],
                &kappa_ghostbox_lower[1],
                &kappa_ghostbox_upper[1],
                &kappa_ghostbox_lower[2],
                &kappa_ghostbox_upper[2],
                phi_x,
                phi_y,
                phi_z,
                &grad_phi_ghostbox_lower[0],
                &grad_phi_ghostbox_upper[0],
                &grad_phi_ghostbox_lower[1],
                &grad_phi_ghostbox_upper[1],
                &grad_phi_ghostbox_lower[2],
                &grad_phi_ghostbox_upper[2],
                phi_xx,
                phi_xy,
                phi_xz,
                phi_yy,
                phi_yz,
                phi_zz,
                &hessian_phi_ghostbox_lower[0],
                &hessian_phi_ghostbox_upper[0],
                &hessian_phi_ghostbox_lower[1],
                &hessian_phi_ghostbox_upper[1],
                &hessian_phi_ghostbox_lower[2],
                &hessian_phi_ghostbox_upper[2],
                &fillbox_lower[0],
                &fillbox_upper[0],
                &fillbox_lower[1],
                &fillbox_upper[1],
                &fillbox_lower[2],
                &fillbox_upper[2]);

            // compute normal velocity = a - b * kappa
            d_math_ops.scale(velocity_data, -d_b, kappa_data, fillbox);
            if (d_a != 0) {
                d_math_ops.addScalar(velocity_data, velocity_data,
                                     d_a, fillbox);
            }

        } else if (num_dims == 2) {
            // get data pointers
            LSMLIB_REAL* phi = phi_data->getPointer(0);

            LSMLIB_REAL* phi_x = grad_phi_data->getPointer(0);
            LSMLIB_REAL* phi_y = grad_phi_data->getPointer(1);

            LSMLIB_REAL* phi_xx = hessian_phi_data->getPointer(0);
            LSMLIB_REAL* phi_xy = hessian_phi_data->getPointer(1);
            LSMLIB_REAL* phi_yy = hessian_phi_data->getPointer(2);

            LSMLIB_REAL* kappa = kappa_data->getPointer(0);

            // compute gradient of phi
            LSM2D_CENTRAL_GRAD_ORDER2(
                phi_x,
                phi_y,
                &grad_phi_ghostbox_lower[0],
                &grad_phi_ghostbox_upper[0],
                &grad_phi_ghostbox_lower[1],
                &grad_phi_ghostbox_upper[1],
                phi,
                &phi_ghostbox_lower[0],
                &phi_ghostbox_upper[0],
                &phi_ghostbox_lower[1],
                &phi_ghostbox_upper[1],
                &fillbox_lower[0],
                &fillbox_upper[0],
                &fillbox_lower[1],
                &fillbox_upper[1],
                &dx[0],
                &dx[1]);

            // compute Hessian of phi
            LSM2D_CENTRAL_HESSIAN(
                phi_xx,
                phi_xy,
                phi_yy,
                &hessian_phi_ghostbox_lower[0],
                &hessian_phi_ghostbox_upper[0],
                &hessian_phi_ghostbox_lower[1],
                &hessian_phi_ghostbox_upper[1],
                phi,
                &phi_ghostbox_lower[0],
                &phi_ghostbox_upper[0],
                &phi_ghostbox_lower[1],
                &phi_ghostbox_upper[1],
                &fillbox_lower[0],
                &fillbox_upper[0],
                &fillbox_lower[1],
                &fillbox_upper[1],
                &dx[0],
                &dx[1]);

            // compute mean curvature
            LSM2D_COMPUTE_MEAN_CURVATURE(
                kappa,
                &kappa_ghostbox_lower[0],
                &kappa_ghostbox_upper[0],
                &kappa_ghostbox_lower[1],
                &kappa_ghostbox_upper[1],
                phi_x,
                phi_y,
                &grad_phi_ghostbox_lower[0],
                &grad_phi_ghostbox_upper[0],
                &grad_phi_ghostbox_lower[1],
                &grad_phi_ghostbox_upper[1],
                phi_xx,
                phi_xy,
                phi_yy,
                &hessian_phi_ghostbox_lower[0],
                &hessian_phi_ghostbox_upper[0],
                &hessian_phi_ghostbox_lower[1],
                &hessian_phi_ghostbox_upper[1],
                &fillbox_lower[0],
                &fillbox_upper[0],
                &fillbox_lower[1],
                &fillbox_upper[1]);

            // compute normal velocity = a - b * kappa
            d_math_ops.scale(velocity_data, -d_b, kappa_data, fillbox);
            if (d_a != 0) {
                d_math_ops.addScalar(velocity_data, velocity_data,
                                     d_a, fillbox);
            }

        } else {
            TBOX_ERROR("VelocityFieldModule::computeVelocityFieldOnLevel():"
                       << "Invalid number of dimensions. Valid values: 2 and 3."
                       << endl);
        } // end switch over dimensions
    } // loop over patches

#ifndef LSMLIB_DOUBLE_PRECISION
    // Clean up memory
    delete [] dx;
    delete [] x_lower;
#endif
}

void VelocityFieldModule::printClassData(ostream& os) const
{
    os << "\nVelocityFieldModule::printClassData..." << endl;
    os << "VelocityFieldModule: this = "
       << (VelocityFieldModule*)this << endl;
    os << "d_object_name = " << d_object_name << endl;

    // KTC - put more here...
    os << endl;
}

void VelocityFieldModule::getFromInput(
    boost::shared_ptr<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(db);
#endif

    d_a = db->getDoubleWithDefault("a", 0);
    d_b = db->getDouble("b");
    d_use_field_extension = db->getBoolWithDefault("use_field_extension",
                                                   true);
}
