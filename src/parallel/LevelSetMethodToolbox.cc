/*
 * File:        LevelSetMethodToolbox.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Implementation file for level set method toolbox class
 */

#ifndef included_LevelSetMethodToolbox_cc
#define included_LevelSetMethodToolbox_cc

#include "LevelSetMethodToolbox.h"

// Standard library headers
#include <cstddef>
#include <cfloat>
//#include <vector>

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// MPI headers
#include "mpi.h"

// SAMRAI Headers
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

// LSMLIB headers
#include "LSMLIB_config.h"
#include "LSMLIB_DefaultParameters.h"

extern "C" {
  #include "lsm_geometry1d.h"
  #include "lsm_geometry2d.h"
  #include "lsm_geometry3d.h"
  #include "lsm_spatial_derivatives1d.h"
  #include "lsm_spatial_derivatives2d.h"
  #include "lsm_spatial_derivatives3d.h"
  #include "lsm_samrai_f77_utilities.h"
  #include "lsm_tvd_runge_kutta1d.h"
  #include "lsm_tvd_runge_kutta2d.h"
  #include "lsm_tvd_runge_kutta3d.h"
  #include "lsm_utilities1d.h"
  #include "lsm_utilities2d.h"
  #include "lsm_utilities3d.h"
}

// Class/type declarations
namespace SAMRAI { namespace hier { class Patch; } }
namespace SAMRAI { namespace hier { class VariableContext; } }

// Namespaces
using namespace std;
using namespace SAMRAI;

/****************************************************************
 *
 * Initialize static data members
 *
 * NOTE:  all of the initialization values correspond to
 *        invalid states.
 *
 ****************************************************************/

namespace LSMLIB {

// parameters for computing spatial derivatives
int
LevelSetMethodToolbox::s_D1_one_ghostcell_handle = -1;

int
LevelSetMethodToolbox::s_D1_two_ghostcells_handle = -1;
int
LevelSetMethodToolbox::s_D2_two_ghostcells_handle = -1;

int
LevelSetMethodToolbox::s_D1_three_ghostcells_handle = -1;
int
LevelSetMethodToolbox::s_D2_three_ghostcells_handle = -1;
int
LevelSetMethodToolbox::s_D3_three_ghostcells_handle = -1;


// parameters for computing unit normal vector
int
LevelSetMethodToolbox::s_compute_normal_grad_phi_handle = -1;
int
LevelSetMethodToolbox::s_compute_normal_grad_phi_plus_handle = -1;
int
LevelSetMethodToolbox::s_compute_normal_grad_phi_minus_handle = -1;


/****************************************************************
 *
 * Implementation of LevelSetMethodToolbox Methods
 *
 ****************************************************************/

/* computeUpwindSpatialDerivatives() */
void LevelSetMethodToolbox::computeUpwindSpatialDerivatives(
  boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  const int spatial_derivative_order,
  const int grad_phi_handle,
  const int phi_handle,
  const int upwind_function_handle,
  const int phi_component)
{

  initializeComputeSpatialDerivativesParameters(hierarchy);

  const int finest_level = hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeUpwindSpatialDerivatives(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // compute spatial derivatives for phi
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( grad_phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> upwind_function_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( upwind_function_handle ));

      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch->getPatchGeometry());

  int DIM = hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geom->getDx();
#else
      const double* dx_double = patch_geom->getDx();
      float *dx = new float[DIM];
      for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif

      hier::Box fillbox = grad_phi_data->getBox();
      const hier::IntVector grad_phi_fillbox_lower = fillbox.lower();
      const hier::IntVector grad_phi_fillbox_upper = fillbox.upper();

      hier::Box grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const hier::IntVector grad_phi_ghostbox_lower = grad_phi_ghostbox.lower();
      const hier::IntVector grad_phi_ghostbox_upper = grad_phi_ghostbox.upper();

      hier::Box phi_ghostbox = phi_data->getGhostBox();
      const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
      const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

      hier::Box upwind_fcn_ghostbox = upwind_function_data->getGhostBox();
      const hier::IntVector upwind_fcn_ghostbox_lower =
        upwind_fcn_ghostbox.lower();
      const hier::IntVector upwind_fcn_ghostbox_upper =
        upwind_fcn_ghostbox.upper();

      LSMLIB_REAL* grad_phi[LSM_DIM_MAX];
      LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
      LSMLIB_REAL* upwind_function[LSM_DIM_MAX];
      for (int dim = 0; dim < hierarchy->getDim().getValue(); dim++) {
        grad_phi[dim] = grad_phi_data->getPointer(dim);
        upwind_function[dim] = upwind_function_data->getPointer(dim);
      }

      switch (spatial_derivative_type) {
        case ENO: {
          switch (spatial_derivative_order) {
            case 1: {

              // prepare scratch PatchData
              patch->allocatePatchData( s_D1_one_ghostcell_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
               BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
               patch->getPatchData( s_D1_one_ghostcell_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();

              if ( hierarchy->getDim().getValue() == 3 ) {

                LSM3D_UPWIND_HJ_ENO1(
                  grad_phi[0], grad_phi[1], grad_phi[2],
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
                  upwind_function[0],
                  upwind_function[1],
                  upwind_function[2],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  &upwind_fcn_ghostbox_lower[1],
                  &upwind_fcn_ghostbox_upper[1],
                  &upwind_fcn_ghostbox_lower[2],
                  &upwind_fcn_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &grad_phi_fillbox_lower[2],
                  &grad_phi_fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_UPWIND_HJ_ENO1(
                  grad_phi[0], grad_phi[1],
                  &grad_phi_ghostbox_lower[0],
                  &grad_phi_ghostbox_upper[0],
                  &grad_phi_ghostbox_lower[1],
                  &grad_phi_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  upwind_function[0],
                  upwind_function[1],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  &upwind_fcn_ghostbox_lower[1],
                  &upwind_fcn_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_UPWIND_HJ_ENO1(
                  grad_phi[0],
                  &grad_phi_ghostbox_lower[0],
                  &grad_phi_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  upwind_function[0],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeUpwindSpatialDerivatives(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              } // end switch over dimensions

              // deallocate scratch PatchData
              patch->deallocatePatchData( s_D1_one_ghostcell_handle );

              break;
            }
            case 2: {

              // prepare scratch PatchData
              patch->allocatePatchData( s_D1_two_ghostcells_handle );
              patch->allocatePatchData( s_D2_two_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D1_two_ghostcells_handle));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D2_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D2_two_ghostcells_handle));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();
              hier::Box D2_ghostbox = D2_data->getGhostBox();
              const hier::IntVector D2_ghostbox_lower = D2_ghostbox.lower();
              const hier::IntVector D2_ghostbox_upper = D2_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();
              LSMLIB_REAL* D2 = D2_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_UPWIND_HJ_ENO2(
                  grad_phi[0], grad_phi[1], grad_phi[2],
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
                  upwind_function[0],
                  upwind_function[1],
                  upwind_function[2],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  &upwind_fcn_ghostbox_lower[1],
                  &upwind_fcn_ghostbox_upper[1],
                  &upwind_fcn_ghostbox_lower[2],
                  &upwind_fcn_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &D2_ghostbox_lower[2],
                  &D2_ghostbox_upper[2],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &grad_phi_fillbox_lower[2],
                  &grad_phi_fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_UPWIND_HJ_ENO2(
                  grad_phi[0], grad_phi[1],
                  &grad_phi_ghostbox_lower[0],
                  &grad_phi_ghostbox_upper[0],
                  &grad_phi_ghostbox_lower[1],
                  &grad_phi_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  upwind_function[0],
                  upwind_function[1],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  &upwind_fcn_ghostbox_lower[1],
                  &upwind_fcn_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_UPWIND_HJ_ENO2(
                  grad_phi[0],
                  &grad_phi_ghostbox_lower[0],
                  &grad_phi_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  upwind_function[0],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeUpwindSpatialDerivatives(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              } // end switch over dimensions

              // deallocate scratch PatchData
              patch->deallocatePatchData( s_D1_two_ghostcells_handle );
              patch->deallocatePatchData( s_D2_two_ghostcells_handle );

              break;
            }
            case 3: {

              // prepare scratch PatchData
              patch->allocatePatchData( s_D1_three_ghostcells_handle );
              patch->allocatePatchData( s_D2_three_ghostcells_handle );
              patch->allocatePatchData( s_D3_three_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D1_three_ghostcells_handle ));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D2_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D2_three_ghostcells_handle ));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D3_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D3_three_ghostcells_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();
              hier::Box D2_ghostbox = D2_data->getGhostBox();
              const hier::IntVector D2_ghostbox_lower = D2_ghostbox.lower();
              const hier::IntVector D2_ghostbox_upper = D2_ghostbox.upper();
              hier::Box D3_ghostbox = D3_data->getGhostBox();
              const hier::IntVector D3_ghostbox_lower = D3_ghostbox.lower();
              const hier::IntVector D3_ghostbox_upper = D3_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();
              LSMLIB_REAL* D2 = D2_data->getPointer();
              LSMLIB_REAL* D3 = D3_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_UPWIND_HJ_ENO3(
                  grad_phi[0], grad_phi[1], grad_phi[2],
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
                  upwind_function[0],
                  upwind_function[1],
                  upwind_function[2],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  &upwind_fcn_ghostbox_lower[1],
                  &upwind_fcn_ghostbox_upper[1],
                  &upwind_fcn_ghostbox_lower[2],
                  &upwind_fcn_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &D2_ghostbox_lower[2],
                  &D2_ghostbox_upper[2],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &D3_ghostbox_lower[1],
                  &D3_ghostbox_upper[1],
                  &D3_ghostbox_lower[2],
                  &D3_ghostbox_upper[2],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &grad_phi_fillbox_lower[2],
                  &grad_phi_fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_UPWIND_HJ_ENO3(
                  grad_phi[0], grad_phi[1],
                  &grad_phi_ghostbox_lower[0],
                  &grad_phi_ghostbox_upper[0],
                  &grad_phi_ghostbox_lower[1],
                  &grad_phi_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  upwind_function[0],
                  upwind_function[1],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  &upwind_fcn_ghostbox_lower[1],
                  &upwind_fcn_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &D3_ghostbox_lower[1],
                  &D3_ghostbox_upper[1],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_UPWIND_HJ_ENO3(
                  grad_phi[0],
                  &grad_phi_ghostbox_lower[0],
                  &grad_phi_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  upwind_function[0],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeUpwindSpatialDerivatives(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );
              }

              // deallocate scratch PatchData
              patch->deallocatePatchData( s_D1_three_ghostcells_handle );
              patch->deallocatePatchData( s_D2_three_ghostcells_handle );
              patch->deallocatePatchData( s_D3_three_ghostcells_handle );

              break;
            }
            default: {
              TBOX_ERROR(  "LevelSetMethodToolbox::"
                        << "computeUpwindSpatialDerivatives(): "
                        << "Unsupported order for ENO derivative.  "
                        << "Only ENO1, ENO2, and ENO3 supported."
                        << endl );
            }
          } // end switch on ENO spatial derivative order

          break;
        } // end case ENO

        case WENO: {
          switch (spatial_derivative_order) {
            case 5: {

              // prepare scratch PatchData
              patch->allocatePatchData( s_D1_three_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D1_three_ghostcells_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_UPWIND_HJ_WENO5(
                  grad_phi[0], grad_phi[1], grad_phi[2],
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
                  upwind_function[0],
                  upwind_function[1],
                  upwind_function[2],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  &upwind_fcn_ghostbox_lower[1],
                  &upwind_fcn_ghostbox_upper[1],
                  &upwind_fcn_ghostbox_lower[2],
                  &upwind_fcn_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &grad_phi_fillbox_lower[2],
                  &grad_phi_fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_UPWIND_HJ_WENO5(
                  grad_phi[0], grad_phi[1],
                  &grad_phi_ghostbox_lower[0],
                  &grad_phi_ghostbox_upper[0],
                  &grad_phi_ghostbox_lower[1],
                  &grad_phi_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  upwind_function[0],
                  upwind_function[1],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  &upwind_fcn_ghostbox_lower[1],
                  &upwind_fcn_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_UPWIND_HJ_WENO5(
                  grad_phi[0],
                  &grad_phi_ghostbox_lower[0],
                  &grad_phi_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  upwind_function[0],
                  &upwind_fcn_ghostbox_lower[0],
                  &upwind_fcn_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeUpwindSpatialDerivatives(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );
              }

              // deallocate scratch PatchData
              patch->deallocatePatchData( s_D1_three_ghostcells_handle );

              break;
            }
            default: {
              TBOX_ERROR(  "LevelSetMethodToolbox::"
                        << "computeUpwindSpatialDerivatives(): "
                        << "Unsupported order for WENO derivative.  "
                        << "Only WENO5 supported."
                        << endl );
            }

          } // end switch on WENO spatial derivative order

          break;
        } // end case WENO

        default: {
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeUpwindSpatialDerivatives(): "
                    << "Unsupported spatial derivative type.  "
                    << "Only ENO and WENO derivatives are supported."
                    << endl );
        }

      } // end switch on derivative type

#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computePlusAndMinusSpatialDerivatives() */
void LevelSetMethodToolbox::computePlusAndMinusSpatialDerivatives(
  boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  const int spatial_derivative_order,
  const int grad_phi_plus_handle,
  const int grad_phi_minus_handle,
  const int phi_handle,
  const int phi_component)
{
  // make sure that the scratch PatchData handles have been created
  initializeComputeSpatialDerivativesParameters(hierarchy);

  const int finest_level = hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computePlusAndMinusSpatialDerivatives(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // compute spatial derivatives for phi
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
           patch->getPatchGeometry());

      int DIM = hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geom->getDx();
#else
      const double* dx_double = patch_geom->getDx();
      float *dx = new float[DIM];
      for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif

      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_plus_data =
       BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
         patch->getPatchData( grad_phi_plus_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_minus_data =
       BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
         patch->getPatchData( grad_phi_minus_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
       BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
         patch->getPatchData( phi_handle ));

      hier::Box fillbox = grad_phi_plus_data->getBox();
      const hier::IntVector grad_phi_fillbox_lower = fillbox.lower();
      const hier::IntVector grad_phi_fillbox_upper = fillbox.upper();

      hier::Box grad_phi_plus_ghostbox = grad_phi_plus_data->getGhostBox();
      const hier::IntVector grad_phi_plus_ghostbox_lower =
        grad_phi_plus_ghostbox.lower();
      const hier::IntVector grad_phi_plus_ghostbox_upper =
        grad_phi_plus_ghostbox.upper();

      hier::Box grad_phi_minus_ghostbox = grad_phi_minus_data->getGhostBox();
      const hier::IntVector grad_phi_minus_ghostbox_lower =
        grad_phi_minus_ghostbox.lower();
      const hier::IntVector grad_phi_minus_ghostbox_upper =
        grad_phi_minus_ghostbox.upper();

      hier::Box phi_ghostbox = phi_data->getGhostBox();
      const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
      const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

      LSMLIB_REAL* grad_phi_plus[LSM_DIM_MAX];
      LSMLIB_REAL* grad_phi_minus[LSM_DIM_MAX];
      LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
      for (int dim = 0; dim < DIM; dim++) {
        grad_phi_plus[dim] = grad_phi_plus_data->getPointer(dim);
        grad_phi_minus[dim] = grad_phi_minus_data->getPointer(dim);
      }

      switch (spatial_derivative_type) {
        case ENO: {
          switch (spatial_derivative_order) {
            case 1: {

              // prepare scratch PatchData
              patch->allocatePatchData( s_D1_one_ghostcell_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
               BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
               patch->getPatchData( s_D1_one_ghostcell_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_HJ_ENO1(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &grad_phi_fillbox_lower[2],
                  &grad_phi_fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_ENO1(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_ENO1(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computePlusAndMinusSpatialDerivatives(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate scratch PatchData
              patch->deallocatePatchData( s_D1_one_ghostcell_handle );

              break;
            }

            case 2: {

              // prepare scratch PatchData
              patch->allocatePatchData( s_D1_two_ghostcells_handle );
              patch->allocatePatchData( s_D2_two_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
               BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
               patch->getPatchData( s_D1_two_ghostcells_handle));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D2_data =
               BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
               patch->getPatchData( s_D2_two_ghostcells_handle));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();
              hier::Box D2_ghostbox = D2_data->getGhostBox();
              const hier::IntVector D2_ghostbox_lower = D2_ghostbox.lower();
              const hier::IntVector D2_ghostbox_upper = D2_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();
              LSMLIB_REAL* D2 = D2_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_HJ_ENO2(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &D2_ghostbox_lower[2],
                  &D2_ghostbox_upper[2],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &grad_phi_fillbox_lower[2],
                  &grad_phi_fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_ENO2(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_ENO2(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computePlusAndMinusSpatialDerivatives(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate scratch PatchData
              patch->deallocatePatchData( s_D1_two_ghostcells_handle );
              patch->deallocatePatchData( s_D2_two_ghostcells_handle );

              break;
            }

            case 3: {

              // prepare scratch PatchData
              patch->allocatePatchData( s_D1_three_ghostcells_handle );
              patch->allocatePatchData( s_D2_three_ghostcells_handle );
              patch->allocatePatchData( s_D3_three_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                  patch->getPatchData( s_D1_three_ghostcells_handle ));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D2_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                  patch->getPatchData( s_D2_three_ghostcells_handle ));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D3_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                  patch->getPatchData( s_D3_three_ghostcells_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();
              hier::Box D2_ghostbox = D2_data->getGhostBox();
              const hier::IntVector D2_ghostbox_lower = D2_ghostbox.lower();
              const hier::IntVector D2_ghostbox_upper = D2_ghostbox.upper();
              hier::Box D3_ghostbox = D3_data->getGhostBox();
              const hier::IntVector D3_ghostbox_lower = D3_ghostbox.lower();
              const hier::IntVector D3_ghostbox_upper = D3_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();
              LSMLIB_REAL* D2 = D2_data->getPointer();
              LSMLIB_REAL* D3 = D3_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_HJ_ENO3(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &D2_ghostbox_lower[2],
                  &D2_ghostbox_upper[2],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &D3_ghostbox_lower[1],
                  &D3_ghostbox_upper[1],
                  &D3_ghostbox_lower[2],
                  &D3_ghostbox_upper[2],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &grad_phi_fillbox_lower[2],
                  &grad_phi_fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_ENO3(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &D3_ghostbox_lower[1],
                  &D3_ghostbox_upper[1],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_ENO3(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computePlusAndMinusSpatialDerivatives(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate scratch PatchData
              patch->deallocatePatchData( s_D1_three_ghostcells_handle );
              patch->deallocatePatchData( s_D2_three_ghostcells_handle );
              patch->deallocatePatchData( s_D3_three_ghostcells_handle );

              break;
            }
            default: {
              TBOX_ERROR(  "LevelSetMethodToolbox::"
                        << "computePlusAndMinusSpatialDerivatives(): "
                        << "Unsupported order for ENO derivative.  "
                        << "Only ENO1, ENO2, and ENO3 supported."
                        << endl );
            }
          } // end switch on ENO spatial derivative order

          break;
        } // end case ENO

        case WENO: {
          switch (spatial_derivative_order) {
            case 5: {

              // prepare scratch PatchData
              patch->allocatePatchData( s_D1_three_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
               BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
               patch->getPatchData( s_D1_three_ghostcells_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_HJ_WENO5(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &grad_phi_fillbox_lower[2],
                  &grad_phi_fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_WENO5(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &grad_phi_fillbox_lower[1],
                  &grad_phi_fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_WENO5(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &grad_phi_fillbox_lower[0],
                  &grad_phi_fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computePlusAndMinusSpatialDerivatives(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate scratch PatchData
              patch->deallocatePatchData( s_D1_three_ghostcells_handle );

              break;
            }
            default: {
              TBOX_ERROR(  "LevelSetMethodToolbox::"
                        << "computePlusAndMinusSpatialDerivatives(): "
                        << "Unsupported order for WENO derivative.  "
                        << "Only WENO5 supported."
                        << endl );
            }

          } // end switch on WENO spatial derivative order

          break;
        } // end case WENO

        default: {
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computePlusAndMinusSpatialDerivatives(): "
                    << "Unsupported spatial derivative type.  "
                    << "Only ENO and WENO derivatives are supported."
                    << endl );
        }

      } // end switch on derivative type
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computeCentralSpatialDerivatives() */
void LevelSetMethodToolbox::computeCentralSpatialDerivatives(
  boost::shared_ptr<hier::PatchHierarchy> hierarchy,
  const int spatial_derivative_order,
  const int grad_phi_handle,
  const int phi_handle,
  const int phi_component)
{

  // make sure that the scratch PatchData handles have been created
  initializeComputeSpatialDerivativesParameters(hierarchy);

  const int finest_level = hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
       // loop over patches
       boost::shared_ptr<hier::Patch> patch = *pi;
       if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeCentralSpatialDerivatives(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // compute spatial derivatives for phi
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( grad_phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( phi_handle ));

      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry());

      int DIM = hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geom->getDx();
#else
      const double* dx_double = patch_geom->getDx();
      float *dx = new float[DIM];
      for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif

      hier::Box fillbox = grad_phi_data->getBox();
      const hier::IntVector grad_phi_fillbox_lower = fillbox.lower();
      const hier::IntVector grad_phi_fillbox_upper = fillbox.upper();

      hier::Box grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const hier::IntVector grad_phi_ghostbox_lower =
        grad_phi_ghostbox.lower();
      const hier::IntVector grad_phi_ghostbox_upper =
        grad_phi_ghostbox.upper();

      hier::Box phi_ghostbox = phi_data->getGhostBox();
      const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
      const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

      LSMLIB_REAL* grad_phi[LSM_DIM_MAX];
      LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
      for (int dim = 0; dim < DIM; dim++) {
        grad_phi[dim] = grad_phi_data->getPointer(dim);
      }

     switch (spatial_derivative_order) {
       case 1:
       case 2: {

          if ( DIM == 3 ) {

            LSM3D_CENTRAL_GRAD_ORDER2(
              grad_phi[0], grad_phi[1], grad_phi[2],
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
              &grad_phi_fillbox_lower[0],
              &grad_phi_fillbox_upper[0],
              &grad_phi_fillbox_lower[1],
              &grad_phi_fillbox_upper[1],
              &grad_phi_fillbox_lower[2],
              &grad_phi_fillbox_upper[2],
              &dx[0], &dx[1], &dx[2]);

          } else if ( DIM == 2 ) {

            LSM2D_CENTRAL_GRAD_ORDER2(
              grad_phi[0], grad_phi[1],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              &grad_phi_ghostbox_lower[1],
              &grad_phi_ghostbox_upper[1],
              phi,
              &phi_ghostbox_lower[0],
              &phi_ghostbox_upper[0],
              &phi_ghostbox_lower[1],
              &phi_ghostbox_upper[1],
              &grad_phi_fillbox_lower[0],
              &grad_phi_fillbox_upper[0],
              &grad_phi_fillbox_lower[1],
              &grad_phi_fillbox_upper[1],
              &dx[0], &dx[1]);

          } else if ( DIM == 1 ) {

            LSM1D_CENTRAL_GRAD_ORDER2(
              grad_phi[0],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              phi,
              &phi_ghostbox_lower[0],
              &phi_ghostbox_upper[0],
              &grad_phi_fillbox_lower[0],
              &grad_phi_fillbox_upper[0],
              &dx[0]);

          } else {

            TBOX_ERROR(  "LevelSetMethodToolbox:: "
                      << "computeCentralSpatialDerivatives(): "
                      << "Invalid value of DIM.  "
                      << "Only DIM = 1, 2, and 3 are supported."
                      << endl );

          }

          break;
        }

        case 3:
        case 4: {

          if ( DIM == 3 ) {

            LSM3D_CENTRAL_GRAD_ORDER4(
              grad_phi[0], grad_phi[1], grad_phi[2],
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
              &grad_phi_fillbox_lower[0],
              &grad_phi_fillbox_upper[0],
              &grad_phi_fillbox_lower[1],
              &grad_phi_fillbox_upper[1],
              &grad_phi_fillbox_lower[2],
              &grad_phi_fillbox_upper[2],
              &dx[0], &dx[1], &dx[2]);

          } else if ( DIM == 2 ) {

            LSM2D_CENTRAL_GRAD_ORDER4(
              grad_phi[0], grad_phi[1],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              &grad_phi_ghostbox_lower[1],
              &grad_phi_ghostbox_upper[1],
              phi,
              &phi_ghostbox_lower[0],
              &phi_ghostbox_upper[0],
              &phi_ghostbox_lower[1],
              &phi_ghostbox_upper[1],
              &grad_phi_fillbox_lower[0],
              &grad_phi_fillbox_upper[0],
              &grad_phi_fillbox_lower[1],
              &grad_phi_fillbox_upper[1],
              &dx[0], &dx[1]);

          } else if ( DIM == 1 ) {

            LSM1D_CENTRAL_GRAD_ORDER4(
              grad_phi[0],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              phi,
              &phi_ghostbox_lower[0],
              &phi_ghostbox_upper[0],
              &grad_phi_fillbox_lower[0],
              &grad_phi_fillbox_upper[0],
              &dx[0]);

          } else {

            TBOX_ERROR(  "LevelSetMethodToolbox::"
                      << "computeCentralSpatialDerivatives(): "
                      << "Invalid value of DIM.  "
                      << "Only DIM = 1, 2, and 3 are supported."
                      << endl );

          }

          break;
        }

        default: {
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeCentralSpatialDerivatives(): "
                    << "Unsupported order for central derivative.  "
                    << "Only 1st, 2nd, 3rd, and 4th order supported."
                    << endl );
        }
      } // end switch on spatial derivative order
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
    } // end loop over Patches
  } // end loop over PatchLevels
}


/* TVDRK1Step() */
void LevelSetMethodToolbox::TVDRK1Step(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int u_next_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const LSMLIB_REAL dt,
  const int u_next_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK1Step(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_next_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( u_next_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_cur_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( u_cur_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> rhs_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( rhs_handle ));

      hier::Box u_next_ghostbox = u_next_data->getGhostBox();
      const hier::IntVector u_next_ghostbox_lower = u_next_ghostbox.lower();
      const hier::IntVector u_next_ghostbox_upper = u_next_ghostbox.upper();

      hier::Box u_cur_ghostbox = u_cur_data->getGhostBox();
      const hier::IntVector u_cur_ghostbox_lower = u_cur_ghostbox.lower();
      const hier::IntVector u_cur_ghostbox_upper = u_cur_ghostbox.upper();

      hier::Box rhs_ghostbox = rhs_data->getGhostBox();
      const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
      const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      hier::Box fillbox = rhs_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      LSMLIB_REAL* u_next = u_next_data->getPointer(u_next_component);
      LSMLIB_REAL* u_cur = u_cur_data->getPointer(u_cur_component);
      LSMLIB_REAL* rhs = rhs_data->getPointer(rhs_component);
      int DIM = patch_hierarchy->getDim().getValue();

      if (DIM == 3 ) {
        LSM3D_RK1_STEP(
          u_next,
          &u_next_ghostbox_lower[0],
          &u_next_ghostbox_upper[0],
          &u_next_ghostbox_lower[1],
          &u_next_ghostbox_upper[1],
          &u_next_ghostbox_lower[2],
          &u_next_ghostbox_upper[2],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          &u_cur_ghostbox_lower[2],
          &u_cur_ghostbox_upper[2],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2],
          &dt);

      } else if ( DIM == 2 ) {
        LSM2D_RK1_STEP(
          u_next,
          &u_next_ghostbox_lower[0],
          &u_next_ghostbox_upper[0],
          &u_next_ghostbox_lower[1],
          &u_next_ghostbox_upper[1],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &dt);

      } else if ( DIM == 1 ) {
        LSM1D_RK1_STEP(
          u_next,
          &u_next_ghostbox_lower[0],
          &u_next_ghostbox_upper[0],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &dt);

      } else {  // Unsupported dimension
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK1Step(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

    } // end loop over patches in level
  } // end loop over levels in hierarchy
}


/* TVDRK2Stage1() */
void LevelSetMethodToolbox::TVDRK2Stage1(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int u_stage1_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const LSMLIB_REAL dt,
  const int u_stage1_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK2Stage1(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_stage1_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( u_stage1_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_cur_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( u_cur_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> rhs_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( rhs_handle ));

      hier::Box u_stage1_ghostbox = u_stage1_data->getGhostBox();
      const hier::IntVector u_stage1_ghostbox_lower =
        u_stage1_ghostbox.lower();
      const hier::IntVector u_stage1_ghostbox_upper =
        u_stage1_ghostbox.upper();

      hier::Box u_cur_ghostbox = u_cur_data->getGhostBox();
      const hier::IntVector u_cur_ghostbox_lower =
        u_cur_ghostbox.lower();
      const hier::IntVector u_cur_ghostbox_upper =
        u_cur_ghostbox.upper();

      hier::Box rhs_ghostbox = rhs_data->getGhostBox();
      const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
      const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      hier::Box fillbox = u_stage1_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      LSMLIB_REAL* u_stage1 = u_stage1_data->getPointer(u_stage1_component);
      LSMLIB_REAL* u_cur = u_cur_data->getPointer(u_cur_component);
      LSMLIB_REAL* rhs = rhs_data->getPointer(rhs_component);
      int DIM = patch_hierarchy->getDim().getValue();

      if ( DIM == 3 ) {
        LSM3D_TVD_RK2_STAGE1(
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          &u_stage1_ghostbox_lower[1],
          &u_stage1_ghostbox_upper[1],
          &u_stage1_ghostbox_lower[2],
          &u_stage1_ghostbox_upper[2],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          &u_cur_ghostbox_lower[2],
          &u_cur_ghostbox_upper[2],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2],
          &dt);

      } else if ( DIM == 2 ) {
        LSM2D_TVD_RK2_STAGE1(
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          &u_stage1_ghostbox_lower[1],
          &u_stage1_ghostbox_upper[1],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &dt);

      } else if ( DIM == 1 ) {
        LSM1D_TVD_RK2_STAGE1(
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &dt);

      } else {  // Unsupported dimension
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK2Stage1(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

    } // end loop over patches in level
  } // end loop over levels in hierarchy
}


/* TVDRK2Stage2() */
void LevelSetMethodToolbox::TVDRK2Stage2(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int u_next_handle,
  const int u_stage1_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const LSMLIB_REAL dt,
  const int u_next_component,
  const int u_stage1_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK2Stage2(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_next_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( u_next_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_stage1_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( u_stage1_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_cur_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( u_cur_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> rhs_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( rhs_handle ));

      hier::Box u_next_ghostbox = u_next_data->getGhostBox();
      const hier::IntVector u_next_ghostbox_lower =
        u_next_ghostbox.lower();
      const hier::IntVector u_next_ghostbox_upper =
        u_next_ghostbox.upper();

      hier::Box u_stage1_ghostbox = u_stage1_data->getGhostBox();
      const hier::IntVector u_stage1_ghostbox_lower =
        u_stage1_ghostbox.lower();
      const hier::IntVector u_stage1_ghostbox_upper =
        u_stage1_ghostbox.upper();

      hier::Box u_cur_ghostbox = u_cur_data->getGhostBox();
      const hier::IntVector u_cur_ghostbox_lower =
        u_cur_ghostbox.lower();
      const hier::IntVector u_cur_ghostbox_upper =
        u_cur_ghostbox.upper();

      hier::Box rhs_ghostbox = rhs_data->getGhostBox();
      const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
      const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      hier::Box fillbox = u_next_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      LSMLIB_REAL* u_next = u_next_data->getPointer(u_next_component);
      LSMLIB_REAL* u_stage1 = u_stage1_data->getPointer(u_stage1_component);
      LSMLIB_REAL* u_cur = u_cur_data->getPointer(u_cur_component);
      LSMLIB_REAL* rhs = rhs_data->getPointer(rhs_component);
      int DIM = patch_hierarchy->getDim().getValue();

      if ( DIM == 3 ) {
        LSM3D_TVD_RK2_STAGE2(
          u_next,
          &u_next_ghostbox_lower[0],
          &u_next_ghostbox_upper[0],
          &u_next_ghostbox_lower[1],
          &u_next_ghostbox_upper[1],
          &u_next_ghostbox_lower[2],
          &u_next_ghostbox_upper[2],
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          &u_stage1_ghostbox_lower[1],
          &u_stage1_ghostbox_upper[1],
          &u_stage1_ghostbox_lower[2],
          &u_stage1_ghostbox_upper[2],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          &u_cur_ghostbox_lower[2],
          &u_cur_ghostbox_upper[2],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2],
          &dt);

      } else if ( DIM == 2 ) {
        LSM2D_TVD_RK2_STAGE2(
          u_next,
          &u_next_ghostbox_lower[0],
          &u_next_ghostbox_upper[0],
          &u_next_ghostbox_lower[1],
          &u_next_ghostbox_upper[1],
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          &u_stage1_ghostbox_lower[1],
          &u_stage1_ghostbox_upper[1],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &dt);

      } else if ( DIM == 1 ) {
        LSM1D_TVD_RK2_STAGE2(
          u_next,
          &u_next_ghostbox_lower[0],
          &u_next_ghostbox_upper[0],
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &dt);

      } else {  // Unsupported dimension
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK2Stage2(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

    } // end loop over patches in level
  } // end loop over levels in hierarchy

}


/* TVDRK3Stage1() */
void LevelSetMethodToolbox::TVDRK3Stage1(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int u_stage1_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const LSMLIB_REAL dt,
  const int u_stage1_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK3Stage1(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_stage1_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( u_stage1_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_cur_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( u_cur_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> rhs_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( rhs_handle ));

      hier::Box u_stage1_ghostbox = u_stage1_data->getGhostBox();
      const hier::IntVector u_stage1_ghostbox_lower =
        u_stage1_ghostbox.lower();
      const hier::IntVector u_stage1_ghostbox_upper =
        u_stage1_ghostbox.upper();

      hier::Box u_cur_ghostbox = u_cur_data->getGhostBox();
      const hier::IntVector u_cur_ghostbox_lower =
        u_cur_ghostbox.lower();
      const hier::IntVector u_cur_ghostbox_upper =
        u_cur_ghostbox.upper();

      hier::Box rhs_ghostbox = rhs_data->getGhostBox();
      const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
      const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      hier::Box fillbox = u_stage1_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      LSMLIB_REAL* u_stage1 = u_stage1_data->getPointer(u_stage1_component);
      LSMLIB_REAL* u_cur = u_cur_data->getPointer(u_cur_component);
      LSMLIB_REAL* rhs = rhs_data->getPointer(rhs_component);
      int DIM = patch_hierarchy->getDim().getValue();

      if ( DIM == 3 ) {
        LSM3D_TVD_RK3_STAGE1(
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          &u_stage1_ghostbox_lower[1],
          &u_stage1_ghostbox_upper[1],
          &u_stage1_ghostbox_lower[2],
          &u_stage1_ghostbox_upper[2],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          &u_cur_ghostbox_lower[2],
          &u_cur_ghostbox_upper[2],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2],
          &dt);

      } else if ( DIM == 2 ) {
        LSM2D_TVD_RK3_STAGE1(
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          &u_stage1_ghostbox_lower[1],
          &u_stage1_ghostbox_upper[1],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &dt);

      } else if ( DIM == 1 ) {
        LSM1D_TVD_RK3_STAGE1(
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &dt);

      } else {  // Unsupported dimension
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK3Stage1(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

    } // end loop over patches in level
  } // end loop over levels in hierarchy
}


/* TVDRK3Stage2() */
void LevelSetMethodToolbox::TVDRK3Stage2(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int u_stage2_handle,
  const int u_stage1_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const LSMLIB_REAL dt,
  const int u_stage2_component,
  const int u_stage1_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK3Stage2(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_stage2_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( u_stage2_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_stage1_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( u_stage1_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_cur_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( u_cur_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> rhs_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( rhs_handle ));

      hier::Box u_stage2_ghostbox = u_stage2_data->getGhostBox();
      const hier::IntVector u_stage2_ghostbox_lower =
        u_stage2_ghostbox.lower();
      const hier::IntVector u_stage2_ghostbox_upper =
        u_stage2_ghostbox.upper();

      hier::Box u_stage1_ghostbox = u_stage1_data->getGhostBox();
      const hier::IntVector u_stage1_ghostbox_lower =
        u_stage1_ghostbox.lower();
      const hier::IntVector u_stage1_ghostbox_upper =
        u_stage1_ghostbox.upper();

      hier::Box u_cur_ghostbox = u_cur_data->getGhostBox();
      const hier::IntVector u_cur_ghostbox_lower =
        u_cur_ghostbox.lower();
      const hier::IntVector u_cur_ghostbox_upper =
        u_cur_ghostbox.upper();

      hier::Box rhs_ghostbox = rhs_data->getGhostBox();
      const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
      const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      hier::Box fillbox = u_stage2_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      LSMLIB_REAL* u_stage2 = u_stage2_data->getPointer(u_stage2_component);
      LSMLIB_REAL* u_stage1 = u_stage1_data->getPointer(u_stage1_component);
      LSMLIB_REAL* u_cur = u_cur_data->getPointer(u_cur_component);
      LSMLIB_REAL* rhs = rhs_data->getPointer(rhs_component);
      int DIM = patch_hierarchy->getDim().getValue();

      if ( DIM == 3 ) {
        LSM3D_TVD_RK3_STAGE2(
          u_stage2,
          &u_stage2_ghostbox_lower[0],
          &u_stage2_ghostbox_upper[0],
          &u_stage2_ghostbox_lower[1],
          &u_stage2_ghostbox_upper[1],
          &u_stage2_ghostbox_lower[2],
          &u_stage2_ghostbox_upper[2],
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          &u_stage1_ghostbox_lower[1],
          &u_stage1_ghostbox_upper[1],
          &u_stage1_ghostbox_lower[2],
          &u_stage1_ghostbox_upper[2],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          &u_cur_ghostbox_lower[2],
          &u_cur_ghostbox_upper[2],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2],
          &dt);

      } else if ( DIM == 2 ) {
        LSM2D_TVD_RK3_STAGE2(
          u_stage2,
          &u_stage2_ghostbox_lower[0],
          &u_stage2_ghostbox_upper[0],
          &u_stage2_ghostbox_lower[1],
          &u_stage2_ghostbox_upper[1],
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          &u_stage1_ghostbox_lower[1],
          &u_stage1_ghostbox_upper[1],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &dt);

      } else if ( DIM == 1 ) {
        LSM1D_TVD_RK3_STAGE2(
          u_stage2,
          &u_stage2_ghostbox_lower[0],
          &u_stage2_ghostbox_upper[0],
          u_stage1,
          &u_stage1_ghostbox_lower[0],
          &u_stage1_ghostbox_upper[0],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &dt);

      } else {  // Unsupported dimension
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK3Stage2(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

    } // end loop over patches in level
  } // end loop over levels in hierarchy

}


/* TVDRK3Stage3() */
void LevelSetMethodToolbox::TVDRK3Stage3(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int u_next_handle,
  const int u_stage2_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const LSMLIB_REAL dt,
  const int u_next_component,
  const int u_stage2_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK3Stage3(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_next_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( u_next_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_stage2_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( u_stage2_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> u_cur_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( u_cur_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> rhs_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( rhs_handle ));

      hier::Box u_next_ghostbox = u_next_data->getGhostBox();
      const hier::IntVector u_next_ghostbox_lower =
        u_next_ghostbox.lower();
      const hier::IntVector u_next_ghostbox_upper =
        u_next_ghostbox.upper();

      hier::Box u_stage2_ghostbox = u_stage2_data->getGhostBox();
      const hier::IntVector u_stage2_ghostbox_lower =
        u_stage2_ghostbox.lower();
      const hier::IntVector u_stage2_ghostbox_upper =
        u_stage2_ghostbox.upper();

      hier::Box u_cur_ghostbox = u_cur_data->getGhostBox();
      const hier::IntVector u_cur_ghostbox_lower =
        u_cur_ghostbox.lower();
      const hier::IntVector u_cur_ghostbox_upper =
        u_cur_ghostbox.upper();

      hier::Box rhs_ghostbox = rhs_data->getGhostBox();
      const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
      const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      hier::Box fillbox = u_stage2_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      LSMLIB_REAL* u_next = u_next_data->getPointer(u_next_component);
      LSMLIB_REAL* u_stage2 = u_stage2_data->getPointer(u_stage2_component);
      LSMLIB_REAL* u_cur = u_cur_data->getPointer(u_cur_component);
      LSMLIB_REAL* rhs = rhs_data->getPointer(rhs_component);
      int DIM = patch_hierarchy->getDim().getValue();

      if ( DIM == 3 ) {
        LSM3D_TVD_RK3_STAGE3(
          u_next,
          &u_next_ghostbox_lower[0],
          &u_next_ghostbox_upper[0],
          &u_next_ghostbox_lower[1],
          &u_next_ghostbox_upper[1],
          &u_next_ghostbox_lower[2],
          &u_next_ghostbox_upper[2],
          u_stage2,
          &u_stage2_ghostbox_lower[0],
          &u_stage2_ghostbox_upper[0],
          &u_stage2_ghostbox_lower[1],
          &u_stage2_ghostbox_upper[1],
          &u_stage2_ghostbox_lower[2],
          &u_stage2_ghostbox_upper[2],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          &u_cur_ghostbox_lower[2],
          &u_cur_ghostbox_upper[2],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2],
          &dt);

      } else if ( DIM == 2 ) {
        LSM2D_TVD_RK3_STAGE3(
          u_next,
          &u_next_ghostbox_lower[0],
          &u_next_ghostbox_upper[0],
          &u_next_ghostbox_lower[1],
          &u_next_ghostbox_upper[1],
          u_stage2,
          &u_stage2_ghostbox_lower[0],
          &u_stage2_ghostbox_upper[0],
          &u_stage2_ghostbox_lower[1],
          &u_stage2_ghostbox_upper[1],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          &u_cur_ghostbox_lower[1],
          &u_cur_ghostbox_upper[1],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &dt);

      } else if ( DIM == 1 ) {
        LSM1D_TVD_RK3_STAGE3(
          u_next,
          &u_next_ghostbox_lower[0],
          &u_next_ghostbox_upper[0],
          u_stage2,
          &u_stage2_ghostbox_lower[0],
          &u_stage2_ghostbox_upper[0],
          u_cur,
          &u_cur_ghostbox_lower[0],
          &u_cur_ghostbox_upper[0],
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &dt);

      } else {  // Unsupported dimension
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "TVDRK3Stage3(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

    } // end loop over patches in level
  } // end loop over levels in hierarchy

}


/* computeUnitNormalVectorFromPhi() */
void LevelSetMethodToolbox::computeUnitNormalVectorFromPhi(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  const int spatial_derivative_order,
  const int normal_vector_handle,
  const int phi_handle,
  const UNIT_NORMAL_TYPE unit_normal_type,
  const int phi_component)
{

  // make sure that the scratch PatchData handles have been created
  if ( (s_compute_normal_grad_phi_handle < 0) ||
       (s_compute_normal_grad_phi_plus_handle < 0) ||
       (s_compute_normal_grad_phi_minus_handle < 0) ) {
    initializeComputeUnitNormalParameters(patch_hierarchy);
  }
  initializeComputeSpatialDerivativesParameters(patch_hierarchy);

  const int finest_level = patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                   << "computeUnitNormalVectorFromPhi(): "
                   << "Cannot find patch. Null patch pointer."
                   << endl );
      }

      // allocate scratch space for spatial derivatives
      patch->allocatePatchData(s_compute_normal_grad_phi_handle);
      patch->allocatePatchData(s_compute_normal_grad_phi_plus_handle);
      patch->allocatePatchData(s_compute_normal_grad_phi_minus_handle);

      // get PatchData
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> normal_vector_data=
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( normal_vector_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( s_compute_normal_grad_phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_plus_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( s_compute_normal_grad_phi_plus_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_minus_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( s_compute_normal_grad_phi_minus_handle ));

      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST <geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry());

  int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geom->getDx();
#else
      const double* dx_double = patch_geom->getDx();
      float *dx = new float[DIM];
     for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif

      hier::Box fillbox = normal_vector_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      hier::Box normal_vector_ghostbox = normal_vector_data->getGhostBox();
      const hier::IntVector normal_vector_ghostbox_lower =
        normal_vector_ghostbox.lower();
      const hier::IntVector normal_vector_ghostbox_upper =
        normal_vector_ghostbox.upper();

      hier::Box phi_ghostbox = phi_data->getGhostBox();
      const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
      const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

      hier::Box grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const hier::IntVector grad_phi_ghostbox_lower =
        grad_phi_ghostbox.lower();
      const hier::IntVector grad_phi_ghostbox_upper =
        grad_phi_ghostbox.upper();

      hier::Box grad_phi_plus_ghostbox = grad_phi_plus_data->getGhostBox();
      const hier::IntVector grad_phi_plus_ghostbox_lower =
        grad_phi_plus_ghostbox.lower();
      const hier::IntVector grad_phi_plus_ghostbox_upper =
        grad_phi_plus_ghostbox.upper();

      hier::Box grad_phi_minus_ghostbox = grad_phi_minus_data->getGhostBox();
      const hier::IntVector grad_phi_minus_ghostbox_lower =
        grad_phi_minus_ghostbox.lower();
      const hier::IntVector grad_phi_minus_ghostbox_upper =
        grad_phi_minus_ghostbox.upper();

      LSMLIB_REAL* normal_vector[LSM_DIM_MAX];
      LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
      LSMLIB_REAL* grad_phi[LSM_DIM_MAX];
      LSMLIB_REAL* grad_phi_plus[LSM_DIM_MAX];
      LSMLIB_REAL* grad_phi_minus[LSM_DIM_MAX];
      for (int dim = 0; dim < DIM; dim++) {
        normal_vector[dim] = normal_vector_data->getPointer(dim);
        grad_phi[dim] = grad_phi_data->getPointer(dim);
        grad_phi_plus[dim] = grad_phi_plus_data->getPointer(dim);
        grad_phi_minus[dim] = grad_phi_minus_data->getPointer(dim);
      }

      // compute plus and minus spatial derivatives
      switch (spatial_derivative_type) {
        case ENO: {
          switch (spatial_derivative_order) {
            case 1: {

              // prepare scratch PatchData for computing grad(phi)
              patch->allocatePatchData( s_D1_one_ghostcell_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                  patch->getPatchData( s_D1_one_ghostcell_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();

              if ( DIM == 3 ) {
                LSM3D_HJ_ENO1(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &fillbox_lower[2],
                  &fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_ENO1(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_ENO1(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeUnitNormalVectorFromPhi(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate PatchData for computing grad(phi)
              patch->deallocatePatchData( s_D1_one_ghostcell_handle );

              break;
            }

            case 2: {

              // prepare scratch PatchData for computing grad(phi)
              patch->allocatePatchData( s_D1_two_ghostcells_handle );
              patch->allocatePatchData( s_D2_two_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                  patch->getPatchData( s_D1_two_ghostcells_handle ));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D2_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                  patch->getPatchData( s_D2_two_ghostcells_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();
              hier::Box D2_ghostbox = D2_data->getGhostBox();
              const hier::IntVector D2_ghostbox_lower = D2_ghostbox.lower();
              const hier::IntVector D2_ghostbox_upper = D2_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();
              LSMLIB_REAL* D2 = D2_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_HJ_ENO2(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &D2_ghostbox_lower[2],
                  &D2_ghostbox_upper[2],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &fillbox_lower[2],
                  &fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_ENO2(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_ENO2(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeUnitNormalVectorFromPhi(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate PatchData for computing grad(phi)
              patch->deallocatePatchData( s_D1_two_ghostcells_handle );
              patch->deallocatePatchData( s_D2_two_ghostcells_handle );

              break;
            }
            case 3: {

              // prepare scratch PatchData for computing grad(phi)
              patch->allocatePatchData( s_D1_three_ghostcells_handle );
              patch->allocatePatchData( s_D2_three_ghostcells_handle );
              patch->allocatePatchData( s_D3_three_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D1_three_ghostcells_handle ));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D2_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D2_three_ghostcells_handle ));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D3_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D3_three_ghostcells_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();
              hier::Box D2_ghostbox = D2_data->getGhostBox();
              const hier::IntVector D2_ghostbox_lower = D2_ghostbox.lower();
              const hier::IntVector D2_ghostbox_upper = D2_ghostbox.upper();
              hier::Box D3_ghostbox = D3_data->getGhostBox();
              const hier::IntVector D3_ghostbox_lower = D3_ghostbox.lower();
              const hier::IntVector D3_ghostbox_upper = D3_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();
              LSMLIB_REAL* D2 = D2_data->getPointer();
              LSMLIB_REAL* D3 = D3_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_HJ_ENO3(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &D2_ghostbox_lower[2],
                  &D2_ghostbox_upper[2],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &D3_ghostbox_lower[1],
                  &D3_ghostbox_upper[1],
                  &D3_ghostbox_lower[2],
                  &D3_ghostbox_upper[2],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &fillbox_lower[2],
                  &fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_ENO3(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &D3_ghostbox_lower[1],
                  &D3_ghostbox_upper[1],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_ENO3(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeUnitNormalVectorFromPhi(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate PatchData for computing grad(phi)
              patch->deallocatePatchData( s_D1_three_ghostcells_handle );
              patch->deallocatePatchData( s_D2_three_ghostcells_handle );
              patch->deallocatePatchData( s_D3_three_ghostcells_handle );

              break;
            }
            default: {
              TBOX_ERROR(  "LevelSetMethodToolbox::"
                        << "computeUnitNormalVectorFromPhi(): "
                        << "Unsupported order for ENO derivative.  "
                        << "Only ENO1, ENO2, and ENO3 supported."
                        << endl );
            }
          } // end switch on ENO spatial derivative order

          break;
        } // end case ENO

        case WENO: {
          switch (spatial_derivative_order) {
            case 5: {

              // prepare scratch PatchData for computing grad(phi)
              patch->allocatePatchData( s_D1_three_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                  patch->getPatchData( s_D1_three_ghostcells_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_HJ_WENO5(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &fillbox_lower[2],
                  &fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_WENO5(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_WENO5(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeUnitNormalVectorFromPhi(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate PatchData for computing grad(phi)
              patch->deallocatePatchData( s_D1_three_ghostcells_handle );

              break;
            }
            default: {
              TBOX_ERROR(  "LevelSetMethodToolbox::"
                        << "computeUnitNormalVectorFromPhi(): "
                        << "Unsupported order for WENO derivative.  "
                        << "Only WENO5 supported."
                        << endl );
            }

          } // end switch on WENO spatial derivative order

          break;
        } // end case WENO

        default: {
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeUnitNormalVectorFromPhi(): "
                    << "Unsupported spatial derivative type.  "
                    << "Only ENO and WENO derivatives are supported."
                    << endl );
        }

      } // end computation of spatial derivatives

      // compute unit normal vector
      if ( DIM == 3 ) {
        // compute grad(phi) based on unit_normal_type
        switch (unit_normal_type) {
          case PHI_UPWIND: {
            LSM3D_PHI_UPWIND_GRAD_F(
              grad_phi[0], grad_phi[1], grad_phi[2],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              &grad_phi_ghostbox_lower[1],
              &grad_phi_ghostbox_upper[1],
              &grad_phi_ghostbox_lower[2],
              &grad_phi_ghostbox_upper[2],
              grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              &grad_phi_plus_ghostbox_lower[1],
              &grad_phi_plus_ghostbox_upper[1],
              &grad_phi_plus_ghostbox_lower[2],
              &grad_phi_plus_ghostbox_upper[2],
              grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              &grad_phi_minus_ghostbox_lower[1],
              &grad_phi_minus_ghostbox_upper[1],
              &grad_phi_minus_ghostbox_lower[2],
              &grad_phi_minus_ghostbox_upper[2],
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
              &fillbox_upper[2]);
            break;
          }
          case AVERAGE: {
            LSM3D_AVERAGE_GRAD_PHI(
              grad_phi[0], grad_phi[1], grad_phi[2],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              &grad_phi_ghostbox_lower[1],
              &grad_phi_ghostbox_upper[1],
              &grad_phi_ghostbox_lower[2],
              &grad_phi_ghostbox_upper[2],
              grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              &grad_phi_plus_ghostbox_lower[1],
              &grad_phi_plus_ghostbox_upper[1],
              &grad_phi_plus_ghostbox_lower[2],
              &grad_phi_plus_ghostbox_upper[2],
              grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              &grad_phi_minus_ghostbox_lower[1],
              &grad_phi_minus_ghostbox_upper[1],
              &grad_phi_minus_ghostbox_lower[2],
              &grad_phi_minus_ghostbox_upper[2],
              &fillbox_lower[0],
              &fillbox_upper[0],
              &fillbox_lower[1],
              &fillbox_upper[1],
              &fillbox_lower[2],
              &fillbox_upper[2]);
            break;
          }
          default: {
            TBOX_ERROR(  "LevelSetMethodToolbox::"
                      << "computeUnitNormalVectorFromPhi(): "
                      << "Unknown method for computing spatial derivatives "
                      << "used to compute the unit normal vector."
                      << "Only PHI_UPWIND and AVERAGE are supported."
                      << endl );
          }
        }

        // compute unit normal from grad(phi)
        LSM3D_COMPUTE_UNIT_NORMAL(
          normal_vector[0], normal_vector[1], normal_vector[2],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &normal_vector_ghostbox_lower[1],
          &normal_vector_ghostbox_upper[1],
          &normal_vector_ghostbox_lower[2],
          &normal_vector_ghostbox_upper[2],
          grad_phi[0], grad_phi[1], grad_phi[2],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          &grad_phi_ghostbox_lower[1],
          &grad_phi_ghostbox_upper[1],
          &grad_phi_ghostbox_lower[2],
          &grad_phi_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2]);

      } else if ( DIM == 2 ) {

        // compute grad(phi) based on unit_normal_type
        switch (unit_normal_type) {
          case PHI_UPWIND: {
            LSM2D_PHI_UPWIND_GRAD_F(
              grad_phi[0], grad_phi[1],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              &grad_phi_ghostbox_lower[1],
              &grad_phi_ghostbox_upper[1],
              grad_phi_plus[0], grad_phi_plus[1],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              &grad_phi_plus_ghostbox_lower[1],
              &grad_phi_plus_ghostbox_upper[1],
              grad_phi_minus[0], grad_phi_minus[1],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              &grad_phi_minus_ghostbox_lower[1],
              &grad_phi_minus_ghostbox_upper[1],
              phi,
              &phi_ghostbox_lower[0],
              &phi_ghostbox_upper[0],
              &phi_ghostbox_lower[1],
              &phi_ghostbox_upper[1],
              &fillbox_lower[0],
              &fillbox_upper[0],
              &fillbox_lower[1],
              &fillbox_upper[1]);
            break;
          }
          case AVERAGE: {
            LSM2D_AVERAGE_GRAD_PHI(
              grad_phi[0], grad_phi[1],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              &grad_phi_ghostbox_lower[1],
              &grad_phi_ghostbox_upper[1],
              grad_phi_plus[0], grad_phi_plus[1],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              &grad_phi_plus_ghostbox_lower[1],
              &grad_phi_plus_ghostbox_upper[1],
              grad_phi_minus[0], grad_phi_minus[1],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              &grad_phi_minus_ghostbox_lower[1],
              &grad_phi_minus_ghostbox_upper[1],
              &fillbox_lower[0],
              &fillbox_upper[0],
              &fillbox_lower[1],
              &fillbox_upper[1]);
            break;
          }
          default: {
            TBOX_ERROR(  "LevelSetMethodToolbox::"
                      << "computeUnitNormalVectorFromPhi(): "
                      << "Unknown method for computing spatial derivatives "
                      << "used to compute the unit normal vector."
                      << "Only PHI_UPWIND and AVERAGE are supported."
                      << endl );
          }
        }

        // compute unit normal from grad(phi)
        LSM2D_COMPUTE_UNIT_NORMAL(
          normal_vector[0], normal_vector[1],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &normal_vector_ghostbox_lower[1],
          &normal_vector_ghostbox_upper[1],
          grad_phi[0], grad_phi[1],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          &grad_phi_ghostbox_lower[1],
          &grad_phi_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1]);

      } else if ( DIM == 1 ) {

        // compute grad(phi) based on unit_normal_type
        switch (unit_normal_type) {
          case PHI_UPWIND: {
            LSM1D_PHI_UPWIND_GRAD_F(
              grad_phi[0],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              grad_phi_plus[0],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              grad_phi_minus[0],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              phi,
              &phi_ghostbox_lower[0],
              &phi_ghostbox_upper[0],
              &fillbox_lower[0],
              &fillbox_upper[0]);
            break;
          }
          case AVERAGE: {
            LSM1D_AVERAGE_GRAD_PHI(
              grad_phi[0],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              grad_phi_plus[0],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              grad_phi_minus[0],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              &fillbox_lower[0],
              &fillbox_upper[0]);
            break;
          }
          default: {
            TBOX_ERROR(  "LevelSetMethodToolbox::"
                      << "computeUnitNormalVectorFromPhi(): "
                      << "Unknown method for computing spatial derivatives "
                      << "used to compute the unit normal vector."
                      << "Only PHI_UPWIND and AVERAGE are supported."
                      << endl );
          }
        }

        // compute unit normal from grad(phi)
        LSM1D_COMPUTE_UNIT_NORMAL(
          normal_vector[0],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          grad_phi[0],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0]);

      } else {  // Unsupported dimension

        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeUnitNormalVectorFromPhi(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl );

      } // end switch on DIM

      // deallocate scratch space for spatial derivatives
      patch->deallocatePatchData(s_compute_normal_grad_phi_handle);
      patch->deallocatePatchData(s_compute_normal_grad_phi_plus_handle);
      patch->deallocatePatchData(s_compute_normal_grad_phi_minus_handle);

#ifndef LSMLIB_DOUBLE_PRECISION
 delete [] dx;
#endif
    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computeSignedUnitNormalVectorFromPhi() */
void LevelSetMethodToolbox::computeSignedUnitNormalVectorFromPhi(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  const int spatial_derivative_order,
  const int normal_vector_handle,
  const int phi_handle,
  const UNIT_NORMAL_TYPE unit_normal_type,
  const int phi_component)
{

  // make sure that the scratch PatchData handles have been created
  if ( (s_compute_normal_grad_phi_handle < 0) ||
       (s_compute_normal_grad_phi_plus_handle < 0) ||
       (s_compute_normal_grad_phi_minus_handle < 0) ) {
    initializeComputeUnitNormalParameters(patch_hierarchy);
  }
  initializeComputeSpatialDerivativesParameters(patch_hierarchy);

  const int finest_level = patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeSignedUnitNormalVectorFromPhi(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // allocate scratch space for spatial derivatives
      patch->allocatePatchData(s_compute_normal_grad_phi_handle);
      patch->allocatePatchData(s_compute_normal_grad_phi_plus_handle);
      patch->allocatePatchData(s_compute_normal_grad_phi_minus_handle);

      // get PatchData
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> normal_vector_data=
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( normal_vector_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( s_compute_normal_grad_phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_plus_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( s_compute_normal_grad_phi_plus_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_minus_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( s_compute_normal_grad_phi_minus_handle ));

      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry());

      int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geom->getDx();
#else
      const double* dx_double = patch_geom->getDx();
      float *dx = new float[DIM];
      for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif

      hier::Box fillbox = normal_vector_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      hier::Box normal_vector_ghostbox = normal_vector_data->getGhostBox();
      const hier::IntVector normal_vector_ghostbox_lower =
      normal_vector_ghostbox.lower();
      const hier::IntVector normal_vector_ghostbox_upper =
        normal_vector_ghostbox.upper();

      hier::Box phi_ghostbox = phi_data->getGhostBox();
      const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
      const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

      hier::Box grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const hier::IntVector grad_phi_ghostbox_lower =
        grad_phi_ghostbox.lower();
      const hier::IntVector grad_phi_ghostbox_upper =
        grad_phi_ghostbox.upper();

      hier::Box grad_phi_plus_ghostbox = grad_phi_plus_data->getGhostBox();
      const hier::IntVector grad_phi_plus_ghostbox_lower =
        grad_phi_plus_ghostbox.lower();
      const hier::IntVector grad_phi_plus_ghostbox_upper =
        grad_phi_plus_ghostbox.upper();

      hier::Box grad_phi_minus_ghostbox = grad_phi_minus_data->getGhostBox();
      const hier::IntVector grad_phi_minus_ghostbox_lower =
        grad_phi_minus_ghostbox.lower();
      const hier::IntVector grad_phi_minus_ghostbox_upper =
        grad_phi_minus_ghostbox.upper();

      LSMLIB_REAL* normal_vector[LSM_DIM_MAX];
      LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
      LSMLIB_REAL* grad_phi[LSM_DIM_MAX];
      LSMLIB_REAL* grad_phi_plus[LSM_DIM_MAX];
      LSMLIB_REAL* grad_phi_minus[LSM_DIM_MAX];
      for (int dim = 0; dim < DIM; dim++) {
        normal_vector[dim] = normal_vector_data->getPointer(dim);
        grad_phi[dim] = grad_phi_data->getPointer(dim);
        grad_phi_plus[dim] = grad_phi_plus_data->getPointer(dim);
        grad_phi_minus[dim] = grad_phi_minus_data->getPointer(dim);
      }

      // compute plus and minus spatial derivatives
      switch (spatial_derivative_type) {
        case ENO: {
          switch (spatial_derivative_order) {
            case 1: {

              // prepare scratch PatchData for computing grad(phi)
              patch->allocatePatchData( s_D1_one_ghostcell_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
               BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                 patch->getPatchData( s_D1_one_ghostcell_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();

              if ( DIM == 3 ) {
                LSM3D_HJ_ENO1(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &fillbox_lower[2],
                  &fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_ENO1(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_ENO1(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeSignedUnitNormalVectorFromPhi(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate PatchData for computing grad(phi)
              patch->deallocatePatchData( s_D1_one_ghostcell_handle );

              break;
            }

            case 2: {

              // prepare scratch PatchData for computing grad(phi)
              patch->allocatePatchData( s_D1_two_ghostcells_handle );
              patch->allocatePatchData( s_D2_two_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D1_two_ghostcells_handle ));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D2_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D2_two_ghostcells_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();
              hier::Box D2_ghostbox = D2_data->getGhostBox();
              const hier::IntVector D2_ghostbox_lower = D2_ghostbox.lower();
              const hier::IntVector D2_ghostbox_upper = D2_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();
              LSMLIB_REAL* D2 = D2_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_HJ_ENO2(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &D2_ghostbox_lower[2],
                  &D2_ghostbox_upper[2],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &fillbox_lower[2],
                  &fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_ENO2(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_ENO2(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeSignedUnitNormalVectorFromPhi(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate PatchData for computing grad(phi)
              patch->deallocatePatchData( s_D1_two_ghostcells_handle );
              patch->deallocatePatchData( s_D2_two_ghostcells_handle );

              break;
            }
            case 3: {

              // prepare scratch PatchData for computing grad(phi)
              patch->allocatePatchData( s_D1_three_ghostcells_handle );
              patch->allocatePatchData( s_D2_three_ghostcells_handle );
              patch->allocatePatchData( s_D3_three_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D1_three_ghostcells_handle ));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D2_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D2_three_ghostcells_handle ));
              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D3_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D3_three_ghostcells_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();
              hier::Box D2_ghostbox = D2_data->getGhostBox();
              const hier::IntVector D2_ghostbox_lower = D2_ghostbox.lower();
              const hier::IntVector D2_ghostbox_upper = D2_ghostbox.upper();
              hier::Box D3_ghostbox = D3_data->getGhostBox();
              const hier::IntVector D3_ghostbox_lower = D3_ghostbox.lower();
              const hier::IntVector D3_ghostbox_upper = D3_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();
              LSMLIB_REAL* D2 = D2_data->getPointer();
              LSMLIB_REAL* D3 = D3_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_HJ_ENO3(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  &D2_ghostbox_lower[2],
                  &D2_ghostbox_upper[2],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &D3_ghostbox_lower[1],
                  &D3_ghostbox_upper[1],
                  &D3_ghostbox_lower[2],
                  &D3_ghostbox_upper[2],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &fillbox_lower[2],
                  &fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_ENO3(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  &D2_ghostbox_lower[1],
                  &D2_ghostbox_upper[1],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &D3_ghostbox_lower[1],
                  &D3_ghostbox_upper[1],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_ENO3(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  D2,
                  &D2_ghostbox_lower[0],
                  &D2_ghostbox_upper[0],
                  D3,
                  &D3_ghostbox_lower[0],
                  &D3_ghostbox_upper[0],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeSignedUnitNormalVectorFromPhi(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate PatchData for computing grad(phi)
              patch->deallocatePatchData( s_D1_three_ghostcells_handle );
              patch->deallocatePatchData( s_D2_three_ghostcells_handle );
              patch->deallocatePatchData( s_D3_three_ghostcells_handle );

              break;
            }
            default: {
              TBOX_ERROR(  "LevelSetMethodToolbox::"
                        << "computeSignedUnitNormalVectorFromPhi(): "
                        << "Unsupported order for ENO derivative.  "
                        << "Only ENO1, ENO2, and ENO3 supported."
                        << endl );
            }
          } // end switch on ENO spatial derivative order

          break;
        } // end case ENO

        case WENO: {
          switch (spatial_derivative_order) {
            case 5: {

              // prepare scratch PatchData for computing grad(phi)
              patch->allocatePatchData( s_D1_three_ghostcells_handle );

              boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> D1_data =
                BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
                patch->getPatchData( s_D1_three_ghostcells_handle ));

              hier::Box D1_ghostbox = D1_data->getGhostBox();
              const hier::IntVector D1_ghostbox_lower = D1_ghostbox.lower();
              const hier::IntVector D1_ghostbox_upper = D1_ghostbox.upper();

              LSMLIB_REAL* D1 = D1_data->getPointer();

              if ( DIM == 3 ) {

                LSM3D_HJ_WENO5(
                  grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  &grad_phi_plus_ghostbox_lower[2],
                  &grad_phi_plus_ghostbox_upper[2],
                  grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  &grad_phi_minus_ghostbox_lower[2],
                  &grad_phi_minus_ghostbox_upper[2],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  &phi_ghostbox_lower[2],
                  &phi_ghostbox_upper[2],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &D1_ghostbox_lower[2],
                  &D1_ghostbox_upper[2],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &fillbox_lower[2],
                  &fillbox_upper[2],
                  &dx[0], &dx[1], &dx[2]);

              } else if ( DIM == 2 ) {

                LSM2D_HJ_WENO5(
                  grad_phi_plus[0], grad_phi_plus[1],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  &grad_phi_plus_ghostbox_lower[1],
                  &grad_phi_plus_ghostbox_upper[1],
                  grad_phi_minus[0], grad_phi_minus[1],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  &grad_phi_minus_ghostbox_lower[1],
                  &grad_phi_minus_ghostbox_upper[1],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  &phi_ghostbox_lower[1],
                  &phi_ghostbox_upper[1],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &D1_ghostbox_lower[1],
                  &D1_ghostbox_upper[1],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &fillbox_lower[1],
                  &fillbox_upper[1],
                  &dx[0], &dx[1]);

              } else if ( DIM == 1 ) {

                LSM1D_HJ_WENO5(
                  grad_phi_plus[0],
                  &grad_phi_plus_ghostbox_lower[0],
                  &grad_phi_plus_ghostbox_upper[0],
                  grad_phi_minus[0],
                  &grad_phi_minus_ghostbox_lower[0],
                  &grad_phi_minus_ghostbox_upper[0],
                  phi,
                  &phi_ghostbox_lower[0],
                  &phi_ghostbox_upper[0],
                  D1,
                  &D1_ghostbox_lower[0],
                  &D1_ghostbox_upper[0],
                  &fillbox_lower[0],
                  &fillbox_upper[0],
                  &dx[0]);

              } else {

                TBOX_ERROR(  "LevelSetMethodToolbox::"
                          << "computeSignedUnitNormalVectorFromPhi(): "
                          << "Invalid value of DIM.  "
                          << "Only DIM = 1, 2, and 3 are supported."
                          << endl );

              }

              // deallocate PatchData for computing grad(phi)
              patch->deallocatePatchData( s_D1_three_ghostcells_handle );

              break;
            }
            default: {
              TBOX_ERROR(  "LevelSetMethodToolbox::"
                        << "computeSignedUnitNormalVectorFromPhi(): "
                        << "Unsupported order for WENO derivative.  "
                        << "Only WENO5 supported."
                        << endl );
            }

          } // end switch on WENO spatial derivative order

          break;
        } // end case WENO

        default: {
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeSignedUnitNormalVectorFromPhi(): "
                    << "Unsupported spatial derivative type.  "
                    << "Only ENO and WENO derivatives are supported."
                    << endl );
        }

      } // end computation of spatial derivatives

      // compute unit normal vector
      if ( DIM == 3 ) {
        // compute grad(phi) based on unit_normal_type
        switch (unit_normal_type) {
          case PHI_UPWIND: {
            LSM3D_PHI_UPWIND_GRAD_F(
              grad_phi[0], grad_phi[1], grad_phi[2],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              &grad_phi_ghostbox_lower[1],
              &grad_phi_ghostbox_upper[1],
              &grad_phi_ghostbox_lower[2],
              &grad_phi_ghostbox_upper[2],
              grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              &grad_phi_plus_ghostbox_lower[1],
              &grad_phi_plus_ghostbox_upper[1],
              &grad_phi_plus_ghostbox_lower[2],
              &grad_phi_plus_ghostbox_upper[2],
              grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              &grad_phi_minus_ghostbox_lower[1],
              &grad_phi_minus_ghostbox_upper[1],
              &grad_phi_minus_ghostbox_lower[2],
              &grad_phi_minus_ghostbox_upper[2],
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
              &fillbox_upper[2]);
            break;
          }
          case AVERAGE: {
            LSM3D_AVERAGE_GRAD_PHI(
              grad_phi[0], grad_phi[1], grad_phi[2],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              &grad_phi_ghostbox_lower[1],
              &grad_phi_ghostbox_upper[1],
              &grad_phi_ghostbox_lower[2],
              &grad_phi_ghostbox_upper[2],
              grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              &grad_phi_plus_ghostbox_lower[1],
              &grad_phi_plus_ghostbox_upper[1],
              &grad_phi_plus_ghostbox_lower[2],
              &grad_phi_plus_ghostbox_upper[2],
              grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              &grad_phi_minus_ghostbox_lower[1],
              &grad_phi_minus_ghostbox_upper[1],
              &grad_phi_minus_ghostbox_lower[2],
              &grad_phi_minus_ghostbox_upper[2],
              &fillbox_lower[0],
              &fillbox_upper[0],
              &fillbox_lower[1],
              &fillbox_upper[1],
              &fillbox_lower[2],
              &fillbox_upper[2]);
            break;
          }
          default: {
            TBOX_ERROR(  "LevelSetMethodToolbox::"
                      << "computeSignedUnitNormalVectorFromPhi(): "
                      << "Unknown method for computing spatial derivatives "
                      << "used to compute the unit normal vector."
                      << "Only PHI_UPWIND and AVERAGE are supported."
                      << endl );
          }
        }

        // compute signed unit normal from grad(phi)
        LSM3D_COMPUTE_SIGNED_UNIT_NORMAL(
          normal_vector[0], normal_vector[1], normal_vector[2],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &normal_vector_ghostbox_lower[1],
          &normal_vector_ghostbox_upper[1],
          &normal_vector_ghostbox_lower[2],
          &normal_vector_ghostbox_upper[2],
          grad_phi[0], grad_phi[1], grad_phi[2],
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
          &dx[0], &dx[1], &dx[2]);

      } else if ( DIM == 2 ) {

        // compute grad(phi) based on unit_normal_type
        switch (unit_normal_type) {
          case PHI_UPWIND: {
            LSM2D_PHI_UPWIND_GRAD_F(
              grad_phi[0], grad_phi[1],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              &grad_phi_ghostbox_lower[1],
              &grad_phi_ghostbox_upper[1],
              grad_phi_plus[0], grad_phi_plus[1],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              &grad_phi_plus_ghostbox_lower[1],
              &grad_phi_plus_ghostbox_upper[1],
              grad_phi_minus[0], grad_phi_minus[1],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              &grad_phi_minus_ghostbox_lower[1],
              &grad_phi_minus_ghostbox_upper[1],
              phi,
              &phi_ghostbox_lower[0],
              &phi_ghostbox_upper[0],
              &phi_ghostbox_lower[1],
              &phi_ghostbox_upper[1],
              &fillbox_lower[0],
              &fillbox_upper[0],
              &fillbox_lower[1],
              &fillbox_upper[1]);
            break;
          }
          case AVERAGE: {
            LSM2D_AVERAGE_GRAD_PHI(
              grad_phi[0], grad_phi[1],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              &grad_phi_ghostbox_lower[1],
              &grad_phi_ghostbox_upper[1],
              grad_phi_plus[0], grad_phi_plus[1],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              &grad_phi_plus_ghostbox_lower[1],
              &grad_phi_plus_ghostbox_upper[1],
              grad_phi_minus[0], grad_phi_minus[1],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              &grad_phi_minus_ghostbox_lower[1],
              &grad_phi_minus_ghostbox_upper[1],
              &fillbox_lower[0],
              &fillbox_upper[0],
              &fillbox_lower[1],
              &fillbox_upper[1]);
            break;
          }
          default: {
            TBOX_ERROR(  "LevelSetMethodToolbox::"
                      << "computeSignedUnitNormalVectorFromPhi(): "
                      << "Unknown method for computing spatial derivatives "
                      << "used to compute the unit normal vector."
                      << "Only PHI_UPWIND and AVERAGE are supported."
                      << endl );
          }
        }

        // compute signed unit normal from grad(phi)
        LSM2D_COMPUTE_SIGNED_UNIT_NORMAL(
          normal_vector[0], normal_vector[1],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &normal_vector_ghostbox_lower[1],
          &normal_vector_ghostbox_upper[1],
          grad_phi[0], grad_phi[1],
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
          &dx[0], &dx[1]);

      } else if ( DIM == 1 ) {

        // compute grad(phi) based on unit_normal_type
        switch (unit_normal_type) {
          case PHI_UPWIND: {
            LSM1D_PHI_UPWIND_GRAD_F(
              grad_phi[0],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              grad_phi_plus[0],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              grad_phi_minus[0],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              phi,
              &phi_ghostbox_lower[0],
              &phi_ghostbox_upper[0],
              &fillbox_lower[0],
              &fillbox_upper[0]);
            break;
          }
          case AVERAGE: {
            LSM1D_AVERAGE_GRAD_PHI(
              grad_phi[0],
              &grad_phi_ghostbox_lower[0],
              &grad_phi_ghostbox_upper[0],
              grad_phi_plus[0],
              &grad_phi_plus_ghostbox_lower[0],
              &grad_phi_plus_ghostbox_upper[0],
              grad_phi_minus[0],
              &grad_phi_minus_ghostbox_lower[0],
              &grad_phi_minus_ghostbox_upper[0],
              &fillbox_lower[0],
              &fillbox_upper[0]);
            break;
          }
          default: {
            TBOX_ERROR(  "LevelSetMethodToolbox::"
                      << "computeSignedUnitNormalVectorFromPhi(): "
                      << "Unknown method for computing spatial derivatives "
                      << "used to compute the unit normal vector."
                      << "Only PHI_UPWIND and AVERAGE are supported."
                      << endl );
          }
        }

        // compute signed unit normal from grad(phi)
        LSM1D_COMPUTE_SIGNED_UNIT_NORMAL(
          normal_vector[0],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          grad_phi[0],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &dx[0]);

      } else {  // Unsupported dimension

        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeSignedUnitNormalVectorFromPhi(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl );

      } // end switch on DIM

      // deallocate scratch space for spatial derivatives
      patch->deallocatePatchData(s_compute_normal_grad_phi_handle);
      patch->deallocatePatchData(s_compute_normal_grad_phi_plus_handle);
      patch->deallocatePatchData(s_compute_normal_grad_phi_minus_handle);

#ifndef LSMLIB_DOUBLE_PRECISION
 delete [] dx;
#endif
    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computeUnitNormalVectorFromGradPhi() */
void LevelSetMethodToolbox::computeUnitNormalVectorFromGradPhi(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int normal_vector_handle,
  const int grad_phi_handle)
{
  const int finest_level = patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeUnitNormalVectorFromGradPhi(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // get PatchData
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> normal_vector_data=
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( normal_vector_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( grad_phi_handle ));

      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry());

      hier::Box fillbox = normal_vector_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      hier::Box normal_vector_ghostbox = normal_vector_data->getGhostBox();
      const hier::IntVector normal_vector_ghostbox_lower =
        normal_vector_ghostbox.lower();
      const hier::IntVector normal_vector_ghostbox_upper =
        normal_vector_ghostbox.upper();

      hier::Box grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const hier::IntVector grad_phi_ghostbox_lower =
        grad_phi_ghostbox.lower();
      const hier::IntVector grad_phi_ghostbox_upper =
        grad_phi_ghostbox.upper();

      LSMLIB_REAL* normal_vector[LSM_DIM_MAX];
      LSMLIB_REAL* grad_phi[LSM_DIM_MAX];
      for (int dim = 0; dim < patch_hierarchy->getDim().getValue(); dim++) {
        normal_vector[dim] = normal_vector_data->getPointer(dim);
        grad_phi[dim] = grad_phi_data->getPointer(dim);
      }

      // compute unit normal from grad(phi)
      if ( patch_hierarchy->getDim().getValue() == 3 ) {
        LSM3D_COMPUTE_UNIT_NORMAL(
          normal_vector[0], normal_vector[1], normal_vector[2],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &normal_vector_ghostbox_lower[1],
          &normal_vector_ghostbox_upper[1],
          &normal_vector_ghostbox_lower[2],
          &normal_vector_ghostbox_upper[2],
          grad_phi[0], grad_phi[1], grad_phi[2],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          &grad_phi_ghostbox_lower[1],
          &grad_phi_ghostbox_upper[1],
          &grad_phi_ghostbox_lower[2],
          &grad_phi_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2]);

      } else if ( patch_hierarchy->getDim().getValue()  == 2 ) {

        LSM2D_COMPUTE_UNIT_NORMAL(
          normal_vector[0], normal_vector[1],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &normal_vector_ghostbox_lower[1],
          &normal_vector_ghostbox_upper[1],
          grad_phi[0], grad_phi[1],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          &grad_phi_ghostbox_lower[1],
          &grad_phi_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1]);

      } else if ( patch_hierarchy->getDim().getValue()  == 1 ) {

        LSM1D_COMPUTE_UNIT_NORMAL(
          normal_vector[0],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          grad_phi[0],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0]);

      } else {  // Unsupported dimension

        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeUnitNormalVectorFromGradPhi(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl );

      } // end switch on DIM

    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computeSignedUnitNormalVectorFromGradPhi() */
void LevelSetMethodToolbox::computeSignedUnitNormalVectorFromGradPhi(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int normal_vector_handle,
  const int grad_phi_handle,
  const int phi_handle,
  const int phi_component)
{
  const int finest_level = patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeSignedUnitNormalVectorFromGradPhi(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // get PatchData
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> normal_vector_data=
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( normal_vector_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( grad_phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( phi_handle ));

      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST <geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry());
  int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geom->getDx();
#else
      const double* dx_double = patch_geom->getDx();
      float *dx = new float[DIM];
      for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif

      hier::Box fillbox = normal_vector_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      hier::Box normal_vector_ghostbox = normal_vector_data->getGhostBox();
      const hier::IntVector normal_vector_ghostbox_lower =
        normal_vector_ghostbox.lower();
      const hier::IntVector normal_vector_ghostbox_upper =
        normal_vector_ghostbox.upper();

      hier::Box grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const hier::IntVector grad_phi_ghostbox_lower =
        grad_phi_ghostbox.lower();
      const hier::IntVector grad_phi_ghostbox_upper =
        grad_phi_ghostbox.upper();

      hier::Box phi_ghostbox = phi_data->getGhostBox();
      const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
      const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

      LSMLIB_REAL* normal_vector[LSM_DIM_MAX];
      LSMLIB_REAL* grad_phi[LSM_DIM_MAX];
      LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
      for (int dim = 0; dim < DIM; dim++) {
        normal_vector[dim] = normal_vector_data->getPointer(dim);
        grad_phi[dim] = grad_phi_data->getPointer(dim);
      }

      // compute unit normal from grad(phi)
      if ( DIM == 3 ) {
        LSM3D_COMPUTE_SIGNED_UNIT_NORMAL(
          normal_vector[0], normal_vector[1], normal_vector[2],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &normal_vector_ghostbox_lower[1],
          &normal_vector_ghostbox_upper[1],
          &normal_vector_ghostbox_lower[2],
          &normal_vector_ghostbox_upper[2],
          grad_phi[0], grad_phi[1], grad_phi[2],
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
          &dx[0], &dx[1], &dx[2]);

      } else if ( DIM == 2 ) {

        LSM2D_COMPUTE_SIGNED_UNIT_NORMAL(
          normal_vector[0], normal_vector[1],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          &normal_vector_ghostbox_lower[1],
          &normal_vector_ghostbox_upper[1],
          grad_phi[0], grad_phi[1],
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
          &dx[0], &dx[1]);

      } else if ( DIM == 1 ) {

        LSM1D_COMPUTE_SIGNED_UNIT_NORMAL(
          normal_vector[0],
          &normal_vector_ghostbox_lower[0],
          &normal_vector_ghostbox_upper[0],
          grad_phi[0],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &dx[0]);

      } else {  // Unsupported dimension

        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeSignedUnitNormalVectorFromGradPhi(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl );

      } // end switch on DIM
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
    } // end loop over Patches
 } // end loop over PatchLevels
}


/* computeVolumeOfRegionDefinedByZeroLevelSet() */
LSMLIB_REAL LevelSetMethodToolbox::computeVolumeOfRegionDefinedByZeroLevelSet(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int phi_handle,
  const int control_volume_handle,
  const int region_indicator,
  const int phi_component,
  const int heaviside_width)
{
  LSMLIB_REAL volume = 0.0;

  // loop over PatchHierarchy and compute the integral on each Patch
  // by calling Fortran subroutines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  if (region_indicator > 0) { // integrate over region {x | phi(x) > 0}

    for ( int ln=0 ; ln < num_levels; ln++ ) {

      boost::shared_ptr<hier::PatchLevel> level =
          patch_hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
        // loop over patches
        boost::shared_ptr<hier::Patch> patch = *pi;
        if ( patch==NULL ) {
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeVolumeOfRegionDefinedByZeroLevelSet(): "
                    << "Cannot find patch. Null patch pointer."
                    << endl);
        }

        // get dx and epsilon
        boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
          BOOST_CAST <geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry());
  int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
        const double* dx = patch_geom->getDx();
#else
        const double* dx_double = patch_geom->getDx();
        float *dx = new float[DIM];
        for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif
        LSMLIB_REAL max_dx = dx[0];
        for (int k = 1; k < DIM; k++) {
          if (max_dx < dx[k]) max_dx = dx[k];
        }
        LSMLIB_REAL epsilon = heaviside_width*max_dx;

        // get pointers to data and index space ranges
        boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
          BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( phi_handle ));
        boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> control_volume_data =
          BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( control_volume_handle ));

        hier::Box phi_ghostbox = phi_data->getGhostBox();
        const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
        const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

        hier::Box control_volume_ghostbox = control_volume_data->getGhostBox();
        const hier::IntVector control_volume_ghostbox_lower =
          control_volume_ghostbox.lower();
        const hier::IntVector control_volume_ghostbox_upper =
          control_volume_ghostbox.upper();

        // interior box
        hier::Box interior_box = patch->getBox();
        const hier::IntVector interior_box_lower = interior_box.lower();
        const hier::IntVector interior_box_upper = interior_box.upper();

        LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
        LSMLIB_REAL* control_volume = control_volume_data->getPointer();
        LSMLIB_REAL volume_on_patch = 0.0;
        int control_volume_sgn = 1;

        if ( DIM == 3 ) {
          LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME(
            &volume_on_patch,
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            &phi_ghostbox_lower[1],
            &phi_ghostbox_upper[1],
            &phi_ghostbox_lower[2],
            &phi_ghostbox_upper[2],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_ghostbox_lower[1],
            &control_volume_ghostbox_upper[1],
            &control_volume_ghostbox_lower[2],
            &control_volume_ghostbox_upper[2],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &interior_box_lower[1],
            &interior_box_upper[1],
            &interior_box_lower[2],
            &interior_box_upper[2],
            &dx[0], &dx[1], &dx[2],
            &epsilon);

        } else if ( DIM == 2 ) {
          LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME(
            &volume_on_patch,
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            &phi_ghostbox_lower[1],
            &phi_ghostbox_upper[1],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_ghostbox_lower[1],
            &control_volume_ghostbox_upper[1],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &interior_box_lower[1],
            &interior_box_upper[1],
            &dx[0], &dx[1],
            &epsilon);

        } else if ( DIM == 1 ) {
          LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
            &volume_on_patch,
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &dx[0],
            &epsilon);

        } else {  // Unsupported dimension
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeVolumeOfRegionDefinedByZeroLevelSet(): "
                    << "Invalid value of DIM.  "
                    << "Only DIM = 1, 2, and 3 are supported."
                    << endl);
        }

        volume += volume_on_patch;
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
      } // end loop over patches in level
    } // end loop over levels in hierarchy

  } else { // integrate over region {x | phi(x) < 0}

    for ( int ln=0 ; ln < num_levels; ln++ ) {

      boost::shared_ptr<hier::PatchLevel> level =
          patch_hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
        // loop over patches
        boost::shared_ptr<hier::Patch> patch = *pi;
        if ( patch==NULL ) {
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeVolumeOfRegionDefinedByZeroLevelSet(): "
                    << "Cannot find patch. Null patch pointer."
                    << endl);
        }

        // get dx and epsilon
        boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
          BOOST_CAST <geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry());
  int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
        const double* dx = patch_geom->getDx();
#else
        const double* dx_double = patch_geom->getDx();
        float *dx = new float[DIM];
        for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif
        LSMLIB_REAL max_dx = dx[0];
        for (int k = 1; k < DIM; k++) {
          if (max_dx < dx[k]) max_dx = dx[k];
        }
        LSMLIB_REAL epsilon = heaviside_width*max_dx;

        // get pointers to data and index space ranges
        boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
          BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( phi_handle ));
        boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> control_volume_data =
          BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( control_volume_handle ));

        hier::Box phi_ghostbox = phi_data->getGhostBox();
        const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
        const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

        hier::Box control_volume_ghostbox = control_volume_data->getGhostBox();
        const hier::IntVector control_volume_ghostbox_lower =
          control_volume_ghostbox.lower();
        const hier::IntVector control_volume_ghostbox_upper =
          control_volume_ghostbox.upper();

        // interior box
        hier::Box interior_box = patch->getBox();
        const hier::IntVector interior_box_lower = interior_box.lower();
        const hier::IntVector interior_box_upper = interior_box.upper();

        LSMLIB_REAL* phi = phi_data->getPointer();
        LSMLIB_REAL* control_volume = control_volume_data->getPointer();
        LSMLIB_REAL volume_on_patch = 0.0;
        int control_volume_sgn = 1;

        if ( DIM == 3 ) {
          LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
            &volume_on_patch,
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            &phi_ghostbox_lower[1],
            &phi_ghostbox_upper[1],
            &phi_ghostbox_lower[2],
            &phi_ghostbox_upper[2],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_ghostbox_lower[1],
            &control_volume_ghostbox_upper[1],
            &control_volume_ghostbox_lower[2],
            &control_volume_ghostbox_upper[2],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &interior_box_lower[1],
            &interior_box_upper[1],
            &interior_box_lower[2],
            &interior_box_upper[2],
            &dx[0], &dx[1], &dx[2],
            &epsilon);

        } else if ( DIM == 2 ) {
          LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
            &volume_on_patch,
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            &phi_ghostbox_lower[1],
            &phi_ghostbox_upper[1],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_ghostbox_lower[1],
            &control_volume_ghostbox_upper[1],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &interior_box_lower[1],
            &interior_box_upper[1],
            &dx[0], &dx[1],
            &epsilon);

        } else if ( DIM == 1 ) {
          LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
            &volume_on_patch,
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &dx[0],
            &epsilon);

        } else {  // Unsupported dimension
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeVolumeOfRegionDefinedByZeroLevelSet(): "
                    << "Invalid value of DIM.  "
                    << "Only DIM = 1, 2, and 3 are supported."
                    << endl);
        }

        volume += volume_on_patch;

#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
        } // end loop over patches in level
    } // end loop over levels in hierarchy

  } // end if statement on (region_indicator > 0)

  int err = tbox::SAMRAI_MPI::getSAMRAIWorld().AllReduce(&volume,1, MPI_SUM);

  if(err!=0){
      tbox::perr << "Error AllReduce=" <<err<< endl;
  }
  return volume;


}


/* computeVolumeOfZeroLevelSet() */
LSMLIB_REAL LevelSetMethodToolbox::computeVolumeOfZeroLevelSet(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int phi_handle,
  const int grad_phi_handle,
  const int control_volume_handle,
  const int phi_component,
  const int delta_width)
{
  LSMLIB_REAL volume = 0.0;

  // loop over PatchHierarchy and compute the integral on each Patch
  // by calling Fortran subroutines
  const int num_levels = patch_hierarchy->getNumberOfLevels();

  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
       if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeVolumeOfZeroLevelSet(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get dx and epsilon
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST <geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry());
  int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geom->getDx();
#else
      const double* dx_double = patch_geom->getDx();
      float *dx = new float[DIM];
      for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif
      LSMLIB_REAL max_dx = dx[0];
      for (int k = 1; k < DIM; k++) {
        if (max_dx < dx[k]) max_dx = dx[k];
      }
      LSMLIB_REAL epsilon = delta_width*max_dx;

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( grad_phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> control_volume_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( control_volume_handle ));

      hier::Box phi_ghostbox = phi_data->getGhostBox();
      const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
      const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

      hier::Box grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const hier::IntVector grad_phi_ghostbox_lower =
        grad_phi_ghostbox.lower();
      const hier::IntVector grad_phi_ghostbox_upper =
        grad_phi_ghostbox.upper();

      hier::Box control_volume_ghostbox = control_volume_data->getGhostBox();
      const hier::IntVector control_volume_ghostbox_lower =
        control_volume_ghostbox.lower();
      const hier::IntVector control_volume_ghostbox_upper =
        control_volume_ghostbox.upper();

      // interior box
      hier::Box interior_box = patch->getBox();
      const hier::IntVector interior_box_lower = interior_box.lower();
      const hier::IntVector interior_box_upper = interior_box.upper();

      LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
      LSMLIB_REAL* control_volume = control_volume_data->getPointer();
      LSMLIB_REAL* grad_phi[LSM_DIM_MAX];
      for (int k=0; k<DIM; k++) {
        grad_phi[k] = grad_phi_data->getPointer(k);
      }
      LSMLIB_REAL volume_on_patch = 0.0;
      int control_volume_sgn = 1;

      if ( DIM == 3 ) {
        LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_CONTROL_VOLUME(
          &volume_on_patch,
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          &phi_ghostbox_lower[2],
          &phi_ghostbox_upper[2],
          grad_phi[0], grad_phi[1], grad_phi[2],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          &grad_phi_ghostbox_lower[1],
          &grad_phi_ghostbox_upper[1],
          &grad_phi_ghostbox_lower[2],
          &grad_phi_ghostbox_upper[2],
          control_volume,
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_ghostbox_lower[1],
          &control_volume_ghostbox_upper[1],
          &control_volume_ghostbox_lower[2],
          &control_volume_ghostbox_upper[2],
          &control_volume_sgn,
          &interior_box_lower[0],
          &interior_box_upper[0],
          &interior_box_lower[1],
          &interior_box_upper[1],
          &interior_box_lower[2],
          &interior_box_upper[2],
          &dx[0], &dx[1], &dx[2],
          &epsilon);

      } else if ( DIM == 2 ) {
        LSM2D_PERIMETER_ZERO_LEVEL_SET_CONTROL_VOLUME(
          &volume_on_patch,
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          grad_phi[0], grad_phi[1],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          &grad_phi_ghostbox_lower[1],
          &grad_phi_ghostbox_upper[1],
          control_volume,
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_ghostbox_lower[1],
          &control_volume_ghostbox_upper[1],
          &control_volume_sgn,
          &interior_box_lower[0],
          &interior_box_upper[0],
          &interior_box_lower[1],
          &interior_box_upper[1],
          &dx[0], &dx[1],
          &epsilon);

      } else if ( DIM == 1 ) {
        LSM1D_SIZE_ZERO_LEVEL_SET_CONTROL_VOLUME(
          &volume_on_patch,
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          grad_phi[0],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          control_volume,
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_sgn,
          &interior_box_lower[0],
          &interior_box_upper[0],
          &dx[0],
          &epsilon);

      } else {  // Unsupported dimension
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeVolumeOfZeroLevelSet(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

      volume += volume_on_patch;
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
    } // end loop over patches in level
  } // end loop over levels in hierarchy

  int err = tbox::SAMRAI_MPI::getSAMRAIWorld().AllReduce(&volume,1, MPI_SUM);

  if (err!=0) {
      tbox::perr << "Error AllReduce=" << err << endl;
  }
  return volume;
}


/* computeVolumeIntegral() */
LSMLIB_REAL LevelSetMethodToolbox::computeVolumeIntegral(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int F_handle,
  const int phi_handle,
  const int control_volume_handle,
  const int region_indicator,
  const int F_component,
  const int phi_component,
  const int heaviside_width)
{
  LSMLIB_REAL integral_F = 0.0;

  // loop over PatchHierarchy and compute the integral on each Patch
  // by calling Fortran subroutines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  if (region_indicator > 0) { // integrate over region {x | phi(x) > 0}

    for ( int ln=0 ; ln < num_levels; ln++ ) {

      boost::shared_ptr<hier::PatchLevel> level =
          patch_hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
        // loop over patches
        boost::shared_ptr<hier::Patch> patch = *pi;
        if ( patch==NULL ) {
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeVolumeIntegral(): "
                    << "Cannot find patch. Null patch pointer."
                    << endl);
        }

        // get dx and epsilon
        boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
          BOOST_CAST <geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry());
  int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geom->getDx();
#else
      const double* dx_double = patch_geom->getDx();
      float *dx = new float[DIM];
     for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif
        LSMLIB_REAL max_dx = dx[0];
        for (int k = 1; k < DIM; k++) {
          if (max_dx < dx[k]) max_dx = dx[k];
        }
        LSMLIB_REAL epsilon = heaviside_width*max_dx;

        // get pointers to data and index space ranges
        boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> F_data =
          BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( F_handle));
        boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
          BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( phi_handle ));
        boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> control_volume_data =
          BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( control_volume_handle ));

        hier::Box F_ghostbox = F_data->getGhostBox();
        const hier::IntVector F_ghostbox_lower = F_ghostbox.lower();
        const hier::IntVector F_ghostbox_upper = F_ghostbox.upper();

        hier::Box phi_ghostbox = phi_data->getGhostBox();
        const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
        const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

        hier::Box control_volume_ghostbox = control_volume_data->getGhostBox();
        const hier::IntVector control_volume_ghostbox_lower =
          control_volume_ghostbox.lower();
        const hier::IntVector control_volume_ghostbox_upper =
          control_volume_ghostbox.upper();

        // interior box
        hier::Box interior_box = patch->getBox();
        const hier::IntVector interior_box_lower = interior_box.lower();
        const hier::IntVector interior_box_upper = interior_box.upper();

        LSMLIB_REAL* F = F_data->getPointer(F_component);
        LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
        LSMLIB_REAL* control_volume = control_volume_data->getPointer();
        LSMLIB_REAL integral_F_on_patch = 0.0;
        int control_volume_sgn = 1;

        if ( DIM == 3 ) {
          LSM3D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME(
            &integral_F_on_patch,
            F,
            &F_ghostbox_lower[0],
            &F_ghostbox_upper[0],
            &F_ghostbox_lower[1],
            &F_ghostbox_upper[1],
            &F_ghostbox_lower[2],
            &F_ghostbox_upper[2],
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            &phi_ghostbox_lower[1],
            &phi_ghostbox_upper[1],
            &phi_ghostbox_lower[2],
            &phi_ghostbox_upper[2],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_ghostbox_lower[1],
            &control_volume_ghostbox_upper[1],
            &control_volume_ghostbox_lower[2],
            &control_volume_ghostbox_upper[2],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &interior_box_lower[1],
            &interior_box_upper[1],
            &interior_box_lower[2],
            &interior_box_upper[2],
            &dx[0], &dx[1], &dx[2],
            &epsilon);

        } else if ( DIM == 2 ) {
          LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME(
            &integral_F_on_patch,
            F,
            &F_ghostbox_lower[0],
            &F_ghostbox_upper[0],
            &F_ghostbox_lower[1],
            &F_ghostbox_upper[1],
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            &phi_ghostbox_lower[1],
            &phi_ghostbox_upper[1],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_ghostbox_lower[1],
            &control_volume_ghostbox_upper[1],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &interior_box_lower[1],
            &interior_box_upper[1],
            &dx[0], &dx[1],
            &epsilon);

        } else if ( DIM == 1 ) {
          LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME(
            &integral_F_on_patch,
            F,
            &F_ghostbox_lower[0],
            &F_ghostbox_upper[0],
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &dx[0],
            &epsilon);

        } else {  // Unsupported dimension
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeVolumeIntegral(): "
                    << "Invalid value of DIM.  "
                    << "Only DIM = 1, 2, and 3 are supported."
                    << endl);
        }

        integral_F += integral_F_on_patch;
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
      } // end loop over patches in level
    } // end loop over levels in hierarchy

  } else {

    for ( int ln=0 ; ln < num_levels; ln++ ) {

      boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
        // loop over patches
        boost::shared_ptr<hier::Patch> patch = *pi;
        if ( patch==NULL ) {
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeVolumeIntegral(): "
                    << "Cannot find patch. Null patch pointer."
                    << endl);
        }

        // get dx and epsilon
        boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
          BOOST_CAST <geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry());
  int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
        const double* dx = patch_geom->getDx();
#else
        const double* dx_double = patch_geom->getDx();
        float *dx = new float[DIM];
        for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif
        LSMLIB_REAL max_dx = dx[0];
        for (int k = 1; k < DIM; k++) {
          if (max_dx < dx[k]) max_dx = dx[k];
        }
        LSMLIB_REAL epsilon = heaviside_width*max_dx;

        // get pointers to data and index space ranges
        boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> F_data =
          BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( F_handle));
        boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
         BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
         patch->getPatchData( phi_handle ));
        boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> control_volume_data =
          BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
          patch->getPatchData( control_volume_handle ));

        hier::Box F_ghostbox = F_data->getGhostBox();
        const hier::IntVector F_ghostbox_lower = F_ghostbox.lower();
        const hier::IntVector F_ghostbox_upper = F_ghostbox.upper();

        hier::Box phi_ghostbox = phi_data->getGhostBox();
        const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
        const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

        hier::Box control_volume_ghostbox = control_volume_data->getGhostBox();
        const hier::IntVector control_volume_ghostbox_lower =
          control_volume_ghostbox.lower();
        const hier::IntVector control_volume_ghostbox_upper =
          control_volume_ghostbox.upper();

        // interior box
        hier::Box interior_box = patch->getBox();
        const hier::IntVector interior_box_lower = interior_box.lower();
        const hier::IntVector interior_box_upper = interior_box.upper();

        LSMLIB_REAL* F = F_data->getPointer();
        LSMLIB_REAL* phi = phi_data->getPointer();
        LSMLIB_REAL* control_volume = control_volume_data->getPointer();
        LSMLIB_REAL integral_F_on_patch = 0.0;
        int control_volume_sgn = 1;

        if ( DIM == 3 ) {
          LSM3D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
            &integral_F_on_patch,
            F,
            &F_ghostbox_lower[0],
            &F_ghostbox_upper[0],
            &F_ghostbox_lower[1],
            &F_ghostbox_upper[1],
            &F_ghostbox_lower[2],
            &F_ghostbox_upper[2],
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            &phi_ghostbox_lower[1],
            &phi_ghostbox_upper[1],
            &phi_ghostbox_lower[2],
            &phi_ghostbox_upper[2],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_ghostbox_lower[1],
            &control_volume_ghostbox_upper[1],
            &control_volume_ghostbox_lower[2],
            &control_volume_ghostbox_upper[2],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &interior_box_lower[1],
            &interior_box_upper[1],
            &interior_box_lower[2],
            &interior_box_upper[2],
            &dx[0], &dx[1], &dx[2],
            &epsilon);

        } else if ( DIM == 2 ) {
          LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
            &integral_F_on_patch,
            F,
            &F_ghostbox_lower[0],
            &F_ghostbox_upper[0],
            &F_ghostbox_lower[1],
            &F_ghostbox_upper[1],
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            &phi_ghostbox_lower[1],
            &phi_ghostbox_upper[1],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_ghostbox_lower[1],
            &control_volume_ghostbox_upper[1],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &interior_box_lower[1],
            &interior_box_upper[1],
            &dx[0], &dx[1],
            &epsilon);

        } else if ( DIM == 1 ) {
          LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME(
            &integral_F_on_patch,
            F,
            &F_ghostbox_lower[0],
            &F_ghostbox_upper[0],
            phi,
            &phi_ghostbox_lower[0],
            &phi_ghostbox_upper[0],
            control_volume,
            &control_volume_ghostbox_lower[0],
            &control_volume_ghostbox_upper[0],
            &control_volume_sgn,
            &interior_box_lower[0],
            &interior_box_upper[0],
            &dx[0],
            &epsilon);

        } else {  // Unsupported dimension
          TBOX_ERROR(  "LevelSetMethodToolbox::"
                    << "computeVolumeIntegral(): "
                    << "Invalid value of DIM.  "
                    << "Only DIM = 1, 2, and 3 are supported."
                    << endl);
        }

        integral_F += integral_F_on_patch;
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
      } // end loop over patches in level
    } // end loop over levels in hierarchy

  } // end if statement on (region_indicator > 0)

  int err = tbox::SAMRAI_MPI::getSAMRAIWorld().AllReduce(&integral_F,1, MPI_SUM);

  if(err!=0){
      tbox::perr << "Error AllReduce=" << err << endl;
  }
  return integral_F;
}


/* computeSurfaceIntegral() */
LSMLIB_REAL LevelSetMethodToolbox::computeSurfaceIntegral(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int F_handle,
  const int phi_handle,
  const int grad_phi_handle,
  const int control_volume_handle,
  const int F_component,
  const int phi_component,
  const int delta_width)
{
  LSMLIB_REAL integral_F = 0.0;

  // loop over PatchHierarchy and compute the integral on each Patch
  // by calling Fortran subroutines
  const int num_levels = patch_hierarchy->getNumberOfLevels();

  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeSurfaceIntegral(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get dx and epsilon
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        BOOST_CAST <geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry());
  int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
        const double* dx = patch_geom->getDx();
#else
        const double* dx_double = patch_geom->getDx();
        float *dx = new float[DIM];
        for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif
      LSMLIB_REAL max_dx = dx[0];
      for (int k = 1; k < DIM; k++) {
        if (max_dx < dx[k]) max_dx = dx[k];
      }
      LSMLIB_REAL epsilon = delta_width*max_dx;

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> F_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( F_handle));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( grad_phi_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> control_volume_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( control_volume_handle ));

      hier::Box F_ghostbox = F_data->getGhostBox();
      const hier::IntVector F_ghostbox_lower = F_ghostbox.lower();
      const hier::IntVector F_ghostbox_upper = F_ghostbox.upper();

      hier::Box phi_ghostbox = phi_data->getGhostBox();
      const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
      const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();

      hier::Box grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const hier::IntVector grad_phi_ghostbox_lower =
        grad_phi_ghostbox.lower();
      const hier::IntVector grad_phi_ghostbox_upper =
        grad_phi_ghostbox.upper();

      hier::Box control_volume_ghostbox = control_volume_data->getGhostBox();
      const hier::IntVector control_volume_ghostbox_lower =
        control_volume_ghostbox.lower();
      const hier::IntVector control_volume_ghostbox_upper =
        control_volume_ghostbox.upper();

      // interior box
      hier::Box interior_box = patch->getBox();
      const hier::IntVector interior_box_lower = interior_box.lower();
      const hier::IntVector interior_box_upper = interior_box.upper();

      LSMLIB_REAL* F = F_data->getPointer(F_component);
      LSMLIB_REAL* phi = phi_data->getPointer(phi_component);
      LSMLIB_REAL* control_volume = control_volume_data->getPointer();
      LSMLIB_REAL* grad_phi[LSM_DIM_MAX];
      for (int k=0; k<DIM; k++) {
        grad_phi[k] = grad_phi_data->getPointer(k);
      }
      LSMLIB_REAL integral_F_on_patch = 0.0;
      int control_volume_sgn = 1;

      if ( DIM == 3 ) {
        LSM3D_SURFACE_INTEGRAL_CONTROL_VOLUME(
          &integral_F_on_patch,
          F,
          &F_ghostbox_lower[0],
          &F_ghostbox_upper[0],
          &F_ghostbox_lower[1],
          &F_ghostbox_upper[1],
          &F_ghostbox_lower[2],
          &F_ghostbox_upper[2],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          &phi_ghostbox_lower[2],
          &phi_ghostbox_upper[2],
          grad_phi[0], grad_phi[1], grad_phi[2],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          &grad_phi_ghostbox_lower[1],
          &grad_phi_ghostbox_upper[1],
          &grad_phi_ghostbox_lower[2],
          &grad_phi_ghostbox_upper[2],
          control_volume,
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_ghostbox_lower[1],
          &control_volume_ghostbox_upper[1],
          &control_volume_ghostbox_lower[2],
          &control_volume_ghostbox_upper[2],
          &control_volume_sgn,
          &interior_box_lower[0],
          &interior_box_upper[0],
          &interior_box_lower[1],
          &interior_box_upper[1],
          &interior_box_lower[2],
          &interior_box_upper[2],
          &dx[0], &dx[1], &dx[2],
          &epsilon);

      } else if ( DIM == 2 ) {
        LSM2D_SURFACE_INTEGRAL_CONTROL_VOLUME(
          &integral_F_on_patch,
          F,
          &F_ghostbox_lower[0],
          &F_ghostbox_upper[0],
          &F_ghostbox_lower[1],
          &F_ghostbox_upper[1],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          &phi_ghostbox_lower[1],
          &phi_ghostbox_upper[1],
          grad_phi[0], grad_phi[1],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          &grad_phi_ghostbox_lower[1],
          &grad_phi_ghostbox_upper[1],
          control_volume,
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_ghostbox_lower[1],
          &control_volume_ghostbox_upper[1],
          &control_volume_sgn,
          &interior_box_lower[0],
          &interior_box_upper[0],
          &interior_box_lower[1],
          &interior_box_upper[1],
          &dx[0], &dx[1],
          &epsilon);

      } else if ( DIM == 1 ) {
        LSM1D_SURFACE_INTEGRAL_CONTROL_VOLUME(
          &integral_F_on_patch,
          F,
          &F_ghostbox_lower[0],
          &F_ghostbox_upper[0],
          phi,
          &phi_ghostbox_lower[0],
          &phi_ghostbox_upper[0],
          grad_phi[0],
          &grad_phi_ghostbox_lower[0],
          &grad_phi_ghostbox_upper[0],
          control_volume,
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_sgn,
          &interior_box_lower[0],
          &interior_box_upper[0],
          &dx[0],
          &epsilon);

      } else {  // Unsupported dimension
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeSurfaceIntegral(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

      integral_F += integral_F_on_patch;
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
    } // end loop over patches in level
  } // end loop over levels in hierarchy

  int err = tbox::SAMRAI_MPI::getSAMRAIWorld().AllReduce(&integral_F,1, MPI_SUM);

  if(err!=0){
      tbox::perr << "Error AllReduce=" << err << endl;
  }
  return integral_F;
}


/* computeStableAdvectionDt() */
LSMLIB_REAL LevelSetMethodToolbox::computeStableAdvectionDt(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int velocity_handle,
  const int control_volume_handle,
  const LSMLIB_REAL cfl_number)
{
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry =
    BOOST_CAST <geom::CartesianGridGeometry, hier::BaseGridGeometry>(
    patch_hierarchy->getGridGeometry());
  int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
  const double* dx_level0 = grid_geometry->getDx();
#else
  const double* dx_level0_double = grid_geometry->getDx();
  float *dx_level0 = new float[DIM];
  for (int i = 0; i < DIM; i++) dx_level0[i] = (float) dx_level0_double[i];
#endif
  LSMLIB_REAL max_advection_dt = LSMLIB_REAL_MAX;

  // loop over PatchHierarchy and compute the maximum stable
  // advection dt by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);
    const hier::IntVector ratio_to_coarsest = level->getRatioToCoarserLevel();

    LSMLIB_REAL dx[LSM_DIM_MAX];
    for (int dir = 0; dir < patch_hierarchy->getDim().getValue(); dir++) {
      dx[dir] = dx_level0[dir]/ratio_to_coarsest[dir];
    }
    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeStableAdvectionDt(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      LSMLIB_REAL max_advection_dt_on_patch = -1;  // bogus value overwritten
                                              // by Fortran subroutine

      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> vel_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( velocity_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> control_volume_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( control_volume_handle ));

      hier::Box vel_box = vel_data->getBox();
      const hier::IntVector vel_box_lower = vel_box.lower();
      const hier::IntVector vel_box_upper = vel_box.upper();

      hier::Box vel_ghostbox = vel_data->getGhostBox();
      const hier::IntVector vel_ghostbox_lower = vel_ghostbox.lower();
      const hier::IntVector vel_ghostbox_upper = vel_ghostbox.upper();

      hier::Box control_volume_ghostbox = control_volume_data->getGhostBox();
      const hier::IntVector control_volume_ghostbox_lower =
        control_volume_ghostbox.lower();
      const hier::IntVector control_volume_ghostbox_upper =
        control_volume_ghostbox.upper();

      int control_volume_sgn = 1;

      if ( DIM == 3 ){
        LSM3D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME(
          &max_advection_dt_on_patch,
          vel_data->getPointer(0),
          vel_data->getPointer(1),
          vel_data->getPointer(2),
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_ghostbox_lower[2],
          &vel_ghostbox_upper[2],
          control_volume_data->getPointer(),
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_ghostbox_lower[1],
          &control_volume_ghostbox_upper[1],
          &control_volume_ghostbox_lower[2],
          &control_volume_ghostbox_upper[2],
          &control_volume_sgn,
          &vel_box_lower[0],
          &vel_box_upper[0],
          &vel_box_lower[1],
          &vel_box_upper[1],
          &vel_box_lower[2],
          &vel_box_upper[2],
          &dx[0],
          &dx[1],
          &dx[2],
          &cfl_number);
      } else if ( DIM == 2 ) {
        LSM2D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME(
          &max_advection_dt_on_patch,
          vel_data->getPointer(0),
          vel_data->getPointer(1),
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          control_volume_data->getPointer(),
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_ghostbox_lower[1],
          &control_volume_ghostbox_upper[1],
          &control_volume_sgn,
          &vel_box_lower[0],
          &vel_box_upper[0],
          &vel_box_lower[1],
          &vel_box_upper[1],
          &dx[0],
          &dx[1],
          &cfl_number);
      } else if ( DIM == 1 ) {
        LSM1D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME(
          &max_advection_dt_on_patch,
          vel_data->getPointer(0),
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          control_volume_data->getPointer(),
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_sgn,
          &vel_box_lower[0],
          &vel_box_upper[0],
          &dx[0],
          &cfl_number);
      } else { // invalid DIM value
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeStableAdvectionDt(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl );
      } // end switch over dimension (DIM) of level set method calculation

      // update max_advection_dt
      if (max_advection_dt_on_patch < max_advection_dt)
        max_advection_dt = max_advection_dt_on_patch;

    } // end loop over patches in level
  } // end loop over levels in hierarchy
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx_level0;
#endif
  int err = tbox::SAMRAI_MPI::getSAMRAIWorld().AllReduce(&max_advection_dt,1, MPI_MAX);

  if(err!=0){
      tbox::perr << "Error AllReduce=" << err << endl;
  }
  return max_advection_dt;
}


/* computeStableNormalVelocityDt() */
LSMLIB_REAL LevelSetMethodToolbox::computeStableNormalVelocityDt(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int normal_velocity_handle,
  const int grad_phi_plus_handle,
  const int grad_phi_minus_handle,
  const int control_volume_handle,
  const LSMLIB_REAL cfl_number)
{
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry =
    BOOST_CAST <geom::CartesianGridGeometry, hier::BaseGridGeometry>(
    patch_hierarchy->getGridGeometry());
  int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
  const double* dx_level0 = grid_geometry->getDx();
#else
  const double* dx_level0_double = grid_geometry->getDx();
  float *dx_level0 = new float[DIM];
  for (int i = 0; i < DIM; i++) dx_level0[i] = (float) dx_level0_double[i];
#endif
  LSMLIB_REAL max_normal_vel_dt = LSMLIB_REAL_MAX;

  // loop over PatchHierarchy and compute the maximum stable
  // advection dt by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);
    const hier::IntVector ratio_to_coarsest = level->getRatioToCoarserLevel();

    LSMLIB_REAL dx[LSM_DIM_MAX];
    for (int dir = 0; dir < patch_hierarchy->getDim().getValue(); dir++) {
      dx[dir] = dx_level0[dir]/ratio_to_coarsest[dir];
    }
    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeStableNormalVelocityDt(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      LSMLIB_REAL max_normal_vel_dt_on_patch = -1;  // bogus value overwritten
                                               // by Fortran subroutine

      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> vel_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData(normal_velocity_handle));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_plus_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( grad_phi_plus_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> grad_phi_minus_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( grad_phi_minus_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> control_volume_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( control_volume_handle ));

      hier::Box vel_box = vel_data->getBox();
      const hier::IntVector vel_box_lower = vel_box.lower();
      const hier::IntVector vel_box_upper = vel_box.upper();

      hier::Box vel_ghostbox = vel_data->getGhostBox();
      const hier::IntVector vel_ghostbox_lower = vel_ghostbox.lower();
      const hier::IntVector vel_ghostbox_upper = vel_ghostbox.upper();

      hier::Box grad_phi_plus_ghostbox = grad_phi_plus_data->getGhostBox();
      const hier::IntVector grad_phi_plus_ghostbox_lower =
        grad_phi_plus_ghostbox.lower();
      const hier::IntVector grad_phi_plus_ghostbox_upper =
        grad_phi_plus_ghostbox.upper();
      hier::Box grad_phi_minus_ghostbox = grad_phi_minus_data->getGhostBox();
      const hier::IntVector grad_phi_minus_ghostbox_lower =
        grad_phi_minus_ghostbox.lower();
      const hier::IntVector grad_phi_minus_ghostbox_upper =
        grad_phi_minus_ghostbox.upper();

      hier::Box control_volume_ghostbox = control_volume_data->getGhostBox();
      const hier::IntVector control_volume_ghostbox_lower =
        control_volume_ghostbox.lower();
      const hier::IntVector control_volume_ghostbox_upper =
        control_volume_ghostbox.upper();

      int control_volume_sgn = 1;

      if ( DIM == 3 ){
        LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME(
          &max_normal_vel_dt_on_patch,
          vel_data->getPointer(),
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_ghostbox_lower[2],
          &vel_ghostbox_upper[2],
          grad_phi_plus_data->getPointer(0),
          grad_phi_plus_data->getPointer(1),
          grad_phi_plus_data->getPointer(2),
          &grad_phi_plus_ghostbox_lower[0],
          &grad_phi_plus_ghostbox_upper[0],
          &grad_phi_plus_ghostbox_lower[1],
          &grad_phi_plus_ghostbox_upper[1],
          &grad_phi_plus_ghostbox_lower[2],
          &grad_phi_plus_ghostbox_upper[2],
          grad_phi_minus_data->getPointer(0),
          grad_phi_minus_data->getPointer(1),
          grad_phi_minus_data->getPointer(2),
          &grad_phi_minus_ghostbox_lower[0],
          &grad_phi_minus_ghostbox_upper[0],
          &grad_phi_minus_ghostbox_lower[1],
          &grad_phi_minus_ghostbox_upper[1],
          &grad_phi_minus_ghostbox_lower[2],
          &grad_phi_minus_ghostbox_upper[2],
          control_volume_data->getPointer(),
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_ghostbox_lower[1],
          &control_volume_ghostbox_upper[1],
          &control_volume_ghostbox_lower[2],
          &control_volume_ghostbox_upper[2],
          &control_volume_sgn,
          &vel_box_lower[0],
          &vel_box_upper[0],
          &vel_box_lower[1],
          &vel_box_upper[1],
          &vel_box_lower[2],
          &vel_box_upper[2],
          &dx[0],
          &dx[1],
          &dx[2],
          &cfl_number);
      } else if ( DIM == 2 ) {
        LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME(
          &max_normal_vel_dt_on_patch,
          vel_data->getPointer(),
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          grad_phi_plus_data->getPointer(0),
          grad_phi_plus_data->getPointer(1),
          &grad_phi_plus_ghostbox_lower[0],
          &grad_phi_plus_ghostbox_upper[0],
          &grad_phi_plus_ghostbox_lower[1],
          &grad_phi_plus_ghostbox_upper[1],
          grad_phi_minus_data->getPointer(0),
          grad_phi_minus_data->getPointer(1),
          &grad_phi_minus_ghostbox_lower[0],
          &grad_phi_minus_ghostbox_upper[0],
          &grad_phi_minus_ghostbox_lower[1],
          &grad_phi_minus_ghostbox_upper[1],
          control_volume_data->getPointer(),
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_ghostbox_lower[1],
          &control_volume_ghostbox_upper[1],
          &control_volume_sgn,
          &vel_box_lower[0],
          &vel_box_upper[0],
          &vel_box_lower[1],
          &vel_box_upper[1],
          &dx[0],
          &dx[1],
          &cfl_number);
      } else if ( DIM == 1 ) {
        LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME(
          &max_normal_vel_dt_on_patch,
          vel_data->getPointer(),
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          grad_phi_plus_data->getPointer(0),
          &grad_phi_plus_ghostbox_lower[0],
          &grad_phi_plus_ghostbox_upper[0],
          grad_phi_minus_data->getPointer(0),
          &grad_phi_minus_ghostbox_lower[0],
          &grad_phi_minus_ghostbox_upper[0],
          control_volume_data->getPointer(),
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_sgn,
          &vel_box_lower[0],
          &vel_box_upper[0],
          &dx[0],
          &cfl_number);
      } else { // invalid DIM value
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeStableNormalVelocityDt(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl );
      } // end switch over dimension (DIM) of level set method calculation

      // update max_normal_vel_dt
      if (max_normal_vel_dt_on_patch < max_normal_vel_dt)
        max_normal_vel_dt = max_normal_vel_dt_on_patch;

    } // end loop over patches in level
  } // end loop over levels in hierarchy
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx_level0;
#endif
  int err = tbox::SAMRAI_MPI::getSAMRAIWorld().AllReduce(&max_normal_vel_dt,1, MPI_MAX);

  if(err!=0){
      tbox::perr << "Error AllReduce=" << err << endl;
  }
  return max_normal_vel_dt;
}


/* maxNormOfDifference() */
LSMLIB_REAL LevelSetMethodToolbox::maxNormOfDifference(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int field1_handle,
  const int field2_handle,
  const int control_volume_handle,
  const int field1_component,
  const int field2_component)
{
  LSMLIB_REAL max_norm_diff = 0;

  // loop over PatchHierarchy and compute the max norm of (field1-field2)
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberOfLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);
    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "maxNormOfDifference(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> field1_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( field1_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> field2_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( field2_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> control_volume_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( control_volume_handle ));

      hier::Box field1_ghostbox = field1_data->getGhostBox();
      const hier::IntVector field1_ghostbox_lower = field1_ghostbox.lower();
      const hier::IntVector field1_ghostbox_upper = field1_ghostbox.upper();

      hier::Box field2_ghostbox = field2_data->getGhostBox();
      const hier::IntVector field2_ghostbox_lower = field2_ghostbox.lower();
      const hier::IntVector field2_ghostbox_upper = field2_ghostbox.upper();

      hier::Box control_volume_ghostbox =
        control_volume_data->getGhostBox();
      const hier::IntVector control_volume_ghostbox_lower =
        control_volume_ghostbox.lower();
      const hier::IntVector control_volume_ghostbox_upper =
        control_volume_ghostbox.upper();

      // interior box
      hier::Box interior_box = field1_data->getBox();
      const hier::IntVector interior_box_lower = interior_box.lower();
      const hier::IntVector interior_box_upper = interior_box.upper();

      LSMLIB_REAL* field1 = field1_data->getPointer(field1_component);
      LSMLIB_REAL* field2 = field2_data->getPointer(field2_component);
      LSMLIB_REAL* control_volume = control_volume_data->getPointer();
      int control_volume_sgn = 1;

      LSMLIB_REAL max_norm_diff_on_patch = 0.0;

      if ( patch_hierarchy->getDim().getValue() == 3 ) {
        LSM3D_MAX_NORM_DIFF_CONTROL_VOLUME(
          &max_norm_diff_on_patch,
          field1,
          &field1_ghostbox_lower[0],
          &field1_ghostbox_upper[0],
          &field1_ghostbox_lower[1],
          &field1_ghostbox_upper[1],
          &field1_ghostbox_lower[2],
          &field1_ghostbox_upper[2],
          field2,
          &field2_ghostbox_lower[0],
          &field2_ghostbox_upper[0],
          &field2_ghostbox_lower[1],
          &field2_ghostbox_upper[1],
          &field2_ghostbox_lower[2],
          &field2_ghostbox_upper[2],
          control_volume,
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_ghostbox_lower[1],
          &control_volume_ghostbox_upper[1],
          &control_volume_ghostbox_lower[2],
          &control_volume_ghostbox_upper[2],
          &control_volume_sgn,
          &interior_box_lower[0],
          &interior_box_upper[0],
          &interior_box_lower[1],
          &interior_box_upper[1],
          &interior_box_lower[2],
          &interior_box_upper[2]);

      } else if ( patch_hierarchy->getDim().getValue() == 2 ) {
        LSM2D_MAX_NORM_DIFF_CONTROL_VOLUME(
          &max_norm_diff_on_patch,
          field1,
          &field1_ghostbox_lower[0],
          &field1_ghostbox_upper[0],
          &field1_ghostbox_lower[1],
          &field1_ghostbox_upper[1],
          field2,
          &field2_ghostbox_lower[0],
          &field2_ghostbox_upper[0],
          &field2_ghostbox_lower[1],
          &field2_ghostbox_upper[1],
          control_volume,
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_ghostbox_lower[1],
          &control_volume_ghostbox_upper[1],
          &control_volume_sgn,
          &interior_box_lower[0],
          &interior_box_upper[0],
          &interior_box_lower[1],
          &interior_box_upper[1]);

      } else if ( patch_hierarchy->getDim().getValue() == 1 ) {
        LSM1D_MAX_NORM_DIFF_CONTROL_VOLUME(
          &max_norm_diff_on_patch,
          field1,
          &field1_ghostbox_lower[0],
          &field1_ghostbox_upper[0],
          field2,
          &field2_ghostbox_lower[0],
          &field2_ghostbox_upper[0],
          control_volume,
          &control_volume_ghostbox_lower[0],
          &control_volume_ghostbox_upper[0],
          &control_volume_sgn,
          &interior_box_lower[0],
          &interior_box_upper[0]);

      } else {  // Unsupported dimension
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "maxNormOfDifference(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

      if (max_norm_diff < max_norm_diff_on_patch)
        max_norm_diff = max_norm_diff_on_patch;

    } // end loop over patches in level
  } // end loop over levels in hierarchy

  int err = tbox::SAMRAI_MPI::getSAMRAIWorld().AllReduce(&max_norm_diff,1, MPI_MAX);

  if(err!=0){
      tbox::perr << "Error AllReduce=" << err << endl;
  }
  return max_norm_diff;
}


/* computeControlVolumes() */
void LevelSetMethodToolbox::computeControlVolumes(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  int control_volume_handle)
{
  int finest_ln = patch_hierarchy->getFinestLevelNumber();
  int ln;
  for ( ln=finest_ln; ln >= 0; --ln ) {

    /*
     * On every level, first assign cell volume to vector weight.
     */

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry =
        BOOST_CAST <geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry());
 int DIM = patch_hierarchy->getDim().getValue();
#ifdef LSMLIB_DOUBLE_PRECISION
      const double* dx = patch_geometry->getDx();
#else
      const double* dx_double = patch_geometry->getDx();
      float *dx = new float[DIM];
      for (int i = 0; i < DIM; i++) dx[i] = (float) dx_double[i];
#endif
      LSMLIB_REAL cell_vol = dx[0];
      for (int i = 1; i<DIM; i++) {
        cell_vol *= dx[i];
      }

      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> cv_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData(control_volume_handle));
      if ( !cv_data ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeControlVolumes(): "
                  << "control_volume_handle must refer to a "
                  << "pdat_CellVariable"
                  << endl );
      }
      cv_data->fillAll(cell_vol);
#ifndef LSMLIB_DOUBLE_PRECISION
  delete [] dx;
#endif
    }

    /*
     * On all but the finest level, assign 0 to vector
     * weight to cells covered by finer cells.
     */

    if (ln < finest_ln) {

      /*
       * First get the boxes that describe index space of the next finer
       * level and coarsen them to describe corresponding index space
       * at this level.
       */

      boost::shared_ptr<hier::PatchLevel> next_finer_level
            = patch_hierarchy->getPatchLevel(ln+1);
      hier::BoxContainer coarsened_boxes = next_finer_level->getBoxes();
      hier::IntVector coarsen_ratio = next_finer_level->getRatioToCoarserLevel();
      coarsen_ratio /= level->getRatioToCoarserLevel();
      coarsened_boxes.coarsen(coarsen_ratio);

      /*
       * Then set vector weight to 0 wherever there is
       * a nonempty intersection with the next finer level.
       * Note that all assignments are local.
       */

      for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
        // loop over patches
        boost::shared_ptr<hier::Patch> patch = *pi;

        for (hier::BoxContainer::iterator ib = coarsened_boxes.begin();ib != coarsened_boxes.end(); ++ib) {

          hier::Box coarse_box = *ib;
          hier::Box intersection = coarse_box*(patch->getBox());
          if ( !intersection.empty() ) {
            boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> cv_data =
            BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
            patch->getPatchData(control_volume_handle));
            cv_data->fillAll(0.0, intersection);

          }  // assignment only in non-empty intersection
        }  // loop over coarsened boxes from finer level
      }  // loop over patches in level
    }  // all levels except finest
  }  // loop over levels
}


/* copySAMRAIData() */
void LevelSetMethodToolbox::copySAMRAIData(
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
  const int dst_handle,
  const int src_handle,
  const int dst_component,
  const int src_component)
{
  // loop over PatchHierarchy and copy data in Patch interiors
  const int num_levels = patch_hierarchy->getNumberOfLevels();

  for ( int ln=0 ; ln < num_levels; ln++ ) {

    boost::shared_ptr<hier::PatchLevel> level =
        patch_hierarchy->getPatchLevel(ln);

    for (hier::PatchLevel::Iterator pi(level->begin()); pi!=level->end(); pi++) {
      // loop over patches
      boost::shared_ptr<hier::Patch> patch = *pi;
      if ( patch==NULL ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "copySAMRAIData(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> dst_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( dst_handle ));
      boost::shared_ptr<pdat::CellData<LSMLIB_REAL>> src_data =
        BOOST_CAST<pdat::CellData<LSMLIB_REAL>, hier::PatchData>(
        patch->getPatchData( src_handle ));

      hier::Box dst_ghostbox = dst_data->getGhostBox();
      const hier::IntVector dst_ghostbox_lower = dst_ghostbox.lower();
      const hier::IntVector dst_ghostbox_upper = dst_ghostbox.upper();

      hier::Box src_ghostbox = src_data->getGhostBox();
      const hier::IntVector src_ghostbox_lower = src_ghostbox.lower();
      const hier::IntVector src_ghostbox_upper = src_ghostbox.upper();

      hier::Box fillbox = dst_data->getBox();
      const hier::IntVector fillbox_lower = fillbox.lower();
      const hier::IntVector fillbox_upper = fillbox.upper();

      LSMLIB_REAL* dst = dst_data->getPointer(dst_component);
      LSMLIB_REAL* src = src_data->getPointer(src_component);

      if (patch_hierarchy->getDim().getValue() == 3) {

        LSM3D_SAMRAI_UTILITIES_COPY_DATA(
          dst,
          &dst_ghostbox_lower[0],
          &dst_ghostbox_upper[0],
          &dst_ghostbox_lower[1],
          &dst_ghostbox_upper[1],
          &dst_ghostbox_lower[2],
          &dst_ghostbox_upper[2],
          src,
          &src_ghostbox_lower[0],
          &src_ghostbox_upper[0],
          &src_ghostbox_lower[1],
          &src_ghostbox_upper[1],
          &src_ghostbox_lower[2],
          &src_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2]);

      } else if (patch_hierarchy->getDim().getValue()  == 2) {

        LSM2D_SAMRAI_UTILITIES_COPY_DATA(
          dst,
          &dst_ghostbox_lower[0],
          &dst_ghostbox_upper[0],
          &dst_ghostbox_lower[1],
          &dst_ghostbox_upper[1],
          src,
          &src_ghostbox_lower[0],
          &src_ghostbox_upper[0],
          &src_ghostbox_lower[1],
          &src_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1]);

      } else if (patch_hierarchy->getDim().getValue() == 1) {

        LSM1D_SAMRAI_UTILITIES_COPY_DATA(
          dst,
          &dst_ghostbox_lower[0],
          &dst_ghostbox_upper[0],
          src,
          &src_ghostbox_lower[0],
          &src_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0]);

      } else {  // Unsupported dimension
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "copySAMRAIData(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      }

    }  // end loop over Patches
  }  // end loop over PatchLevels
}


/* initializeComputeSpatialDerivativesParameters() */
void LevelSetMethodToolbox::initializeComputeSpatialDerivativesParameters(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy)
{
  // do nothing if PatchData for grad(phi) has already been created
  if (s_D1_one_ghostcell_handle >= 0) {
    return;
  }

  // get pointer to VariableDatabase
  hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

  // create zero ghostcell width IntVector
  // TODO: SET CORRECT VALUES
  hier::IntVector one_ghostcell_width(hierarchy->getDim(), 0);
  hier::IntVector two_ghostcells_width(hierarchy->getDim(), 0);
  hier::IntVector three_ghostcells_width(hierarchy->getDim(), 0);

  // create Variables and VariableContext
  boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>> D1_variable;
  if (var_db->checkVariableExists("D1")){
    D1_variable = BOOST_CAST<pdat::CellVariable<LSMLIB_REAL>, hier::Variable>(
        var_db->getVariable("D1"));
   // VariableDatabaseD1_variable = var_db->getVariable(hierarchy->getDim(),"D1");
    D1_variable = boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>>(
        new pdat::CellVariable<LSMLIB_REAL>(hierarchy->getDim(),"D1", 1));
  }
  boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>> D2_variable;
  if (var_db->checkVariableExists("D2")) {
    D2_variable = BOOST_CAST<pdat::CellVariable<LSMLIB_REAL>, hier::Variable>(
        var_db->getVariable("D2"));
  } else {
    D2_variable = boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>>(
        new pdat::CellVariable<LSMLIB_REAL>(hierarchy->getDim(),"D2", 1));
  }
  boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>> D3_variable;
  if (var_db->checkVariableExists("D3")) {
    D3_variable = BOOST_CAST<pdat::CellVariable<LSMLIB_REAL>, hier::Variable>(
        var_db->getVariable("D3"));
  } else {
    D3_variable = boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>>(
        new pdat::CellVariable<LSMLIB_REAL>(hierarchy->getDim(),"D3", 1));
  }
  boost::shared_ptr<hier::VariableContext> one_ghostcell_context =
    var_db->getContext("LSM_TOOLBOX_SCRATCH::ONE_GHOSTCELL");
  boost::shared_ptr<hier::VariableContext> two_ghostcells_context =
    var_db->getContext("LSM_TOOLBOX_SCRATCH::TWO_GHOSTCELLS");
  boost::shared_ptr<hier::VariableContext> three_ghostcells_context =
    var_db->getContext("LSM_TOOLBOX_SCRATCH::THREE_GHOSTCELLS");

  // register (Variable,VariableContext) pairs with VariableDatabase
  s_D1_one_ghostcell_handle = var_db->registerVariableAndContext(
    D1_variable, one_ghostcell_context, one_ghostcell_width);

  s_D1_two_ghostcells_handle = var_db->registerVariableAndContext(
    D1_variable, two_ghostcells_context, two_ghostcells_width);
  s_D2_two_ghostcells_handle = var_db->registerVariableAndContext(
    D2_variable, two_ghostcells_context, two_ghostcells_width);

  s_D1_three_ghostcells_handle = var_db->registerVariableAndContext(
    D1_variable, three_ghostcells_context, three_ghostcells_width);
  s_D2_three_ghostcells_handle = var_db->registerVariableAndContext(
    D2_variable, three_ghostcells_context, three_ghostcells_width);
  s_D3_three_ghostcells_handle = var_db->registerVariableAndContext(
    D3_variable, three_ghostcells_context, three_ghostcells_width);

}


/* initializeComputeUnitNormalParameters() */
void LevelSetMethodToolbox::initializeComputeUnitNormalParameters(
    boost::shared_ptr<hier::PatchHierarchy> hierarchy)
{
  // do nothing if PatchData for grad(phi) has already been created
  // do nothing if PatchData for grad(phi) has already been created
  // do nothing if PatchData for grad(phi) has already been created
  if ( (s_compute_normal_grad_phi_handle >= 0) &&
       (s_compute_normal_grad_phi_plus_handle >= 0) &&
       (s_compute_normal_grad_phi_minus_handle >= 0) ) {
    return;
  }

  // get pointer to VariableDatabase
  hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

  // create zero ghostcell width IntVector
  hier::IntVector zero_ghostcell_width(hierarchy->getDim(),0);

  // create Variables and VariableContext
  boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>> grad_phi_variable;
  if (var_db->checkVariableExists("grad phi")) {
    grad_phi_variable =
        BOOST_CAST<pdat::CellVariable<LSMLIB_REAL>, hier::Variable>(
            var_db->getVariable("grad phi"));
  } else {
    grad_phi_variable = boost::shared_ptr<pdat::CellVariable<LSMLIB_REAL>>(
        new pdat::CellVariable<LSMLIB_REAL>(hierarchy->getDim(),"grad phi"));
  }
  boost::shared_ptr<hier::VariableContext> scratch_context =
    var_db->getContext("LSM_TOOLBOX_SCRATCH");
  boost::shared_ptr<hier::VariableContext> grad_phi_plus_context =
    var_db->getContext("LSM_TOOLBOX_GRAD_PHI_PLUS");
  boost::shared_ptr<hier::VariableContext> grad_phi_minus_context =
    var_db->getContext("LSM_TOOLBOX_GRAD_PHI_MINUS");

  // register (Variable,VariableContext) pairs with VariableDatabase
  s_compute_normal_grad_phi_handle =
    var_db->registerVariableAndContext(grad_phi_variable,
                                       scratch_context,
                                       zero_ghostcell_width);
  s_compute_normal_grad_phi_plus_handle =
    var_db->registerVariableAndContext(grad_phi_variable,
                                       grad_phi_plus_context,
                                       zero_ghostcell_width);
  s_compute_normal_grad_phi_minus_handle =
    var_db->registerVariableAndContext(grad_phi_variable,
                                       grad_phi_minus_context,
                                       zero_ghostcell_width);

}

} // end LSMLIB namespace

#endif
