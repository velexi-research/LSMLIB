/*
 * File:        LevelSetMethodToolbox.cc
 * Copyright:   (c) 2005-2006 Kevin T. Chu
 * Revision:    $Revision: 1.19 $
 * Modified:    $Date: 2006/11/02 17:23:54 $
 * Description: Implementation file for level set method toolbox class
 */

#ifndef included_LevelSetMethodToolbox_cc
#define included_LevelSetMethodToolbox_cc

// System Headers
#include <vector>
#include <float.h>

#include "LevelSetMethodToolbox.h" 
#include "LSMLIB_DefaultParameters.h"

// SAMRAI Headers
#include "Box.h"
#include "CartesianGridGeometry.h" 
#include "CartesianPatchGeometry.h"
#include "CellData.h" 
#include "CellVariable.h" 
#include "IntVector.h" 
#include "Patch.h" 
#include "PatchLevel.h" 
#include "VariableContext.h" 
#include "VariableDatabase.h" 
#include "tbox/Utilities.h"


// headers for level set method numerical kernels
extern "C" {
  #include "lsm_fast_marching_method.h"
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


// SAMRAI namespaces
using namespace std;
using namespace geom;
using namespace pdat;


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
template <int DIM> int 
LevelSetMethodToolbox<DIM>::s_D1_one_ghostcell_handle = -1;

template <int DIM> int 
LevelSetMethodToolbox<DIM>::s_D1_two_ghostcells_handle = -1;
template <int DIM> int 
LevelSetMethodToolbox<DIM>::s_D2_two_ghostcells_handle = -1;

template <int DIM> int 
LevelSetMethodToolbox<DIM>::s_D1_three_ghostcells_handle = -1;
template <int DIM> int 
LevelSetMethodToolbox<DIM>::s_D2_three_ghostcells_handle = -1;
template <int DIM> int 
LevelSetMethodToolbox<DIM>::s_D3_three_ghostcells_handle = -1;


// parameters for computing unit normal vector
template <int DIM> int 
LevelSetMethodToolbox<DIM>::s_compute_normal_grad_phi_handle = -1;
template <int DIM> int 
LevelSetMethodToolbox<DIM>::s_compute_normal_grad_phi_plus_handle = -1;
template <int DIM> int 
LevelSetMethodToolbox<DIM>::s_compute_normal_grad_phi_minus_handle = -1;


/****************************************************************
 *
 * Implementation of LevelSetMethodToolbox Methods
 *
 ****************************************************************/

/* computeUpwindSpatialDerivatives() */
template <int DIM> 
void LevelSetMethodToolbox<DIM>::computeUpwindSpatialDerivatives(
  Pointer< PatchHierarchy<DIM> > hierarchy,
  const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  const int spatial_derivative_order,
  const int grad_phi_handle,
  const int phi_handle,
  const int upwind_function_handle,
  const int phi_component)
{

  // make sure that the scratch PatchData handles have been created
  initializeComputeSpatialDerivativesParameters();

  const int finest_level = hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    Pointer< PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeUpwindSpatialDerivatives(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // compute spatial derivatives for phi
      Pointer< CellData<DIM,double> > grad_phi_data =
        patch->getPatchData( grad_phi_handle );
      Pointer< CellData<DIM,double> > phi_data =
        patch->getPatchData( phi_handle );
      Pointer< CellData<DIM,double> > upwind_function_data =
        patch->getPatchData( upwind_function_handle );
  
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
  
      Box<DIM> fillbox = grad_phi_data->getBox();
      const IntVector<DIM> grad_phi_fillbox_lower = fillbox.lower();
      const IntVector<DIM> grad_phi_fillbox_upper = fillbox.upper();

      Box<DIM> grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const IntVector<DIM> grad_phi_ghostbox_lower = grad_phi_ghostbox.lower();
      const IntVector<DIM> grad_phi_ghostbox_upper = grad_phi_ghostbox.upper();

      Box<DIM> phi_ghostbox = phi_data->getGhostBox();
      const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
      const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();

      Box<DIM> upwind_fcn_ghostbox = upwind_function_data->getGhostBox();
      const IntVector<DIM> upwind_fcn_ghostbox_lower = 
        upwind_fcn_ghostbox.lower();
      const IntVector<DIM> upwind_fcn_ghostbox_upper = 
        upwind_fcn_ghostbox.upper();

      double* grad_phi[LSM_DIM_MAX];
      double* phi = phi_data->getPointer(phi_component);
      double* upwind_function[LSM_DIM_MAX];
      for (int dim = 0; dim < DIM; dim++) {
        grad_phi[dim] = grad_phi_data->getPointer(dim);
        upwind_function[dim] = upwind_function_data->getPointer(dim);
      }

      switch (spatial_derivative_type) {
        case ENO: {
          switch (spatial_derivative_order) { 
            case 1: {

              // prepare scratch PatchData
              patch->allocatePatchData( s_D1_one_ghostcell_handle );

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_one_ghostcell_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();

              double* D1 = D1_data->getPointer();

              if ( DIM == 3 ) {

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_two_ghostcells_handle );
              Pointer< CellData<DIM,double> > D2_data =
                patch->getPatchData( s_D2_two_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();
              Box<DIM> D2_ghostbox = D2_data->getGhostBox();
              const IntVector<DIM> D2_ghostbox_lower = D2_ghostbox.lower();
              const IntVector<DIM> D2_ghostbox_upper = D2_ghostbox.upper();

              double* D1 = D1_data->getPointer();
              double* D2 = D2_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_three_ghostcells_handle );
              Pointer< CellData<DIM,double> > D2_data =
                patch->getPatchData( s_D2_three_ghostcells_handle );
              Pointer< CellData<DIM,double> > D3_data =
                patch->getPatchData( s_D3_three_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();
              Box<DIM> D2_ghostbox = D2_data->getGhostBox();
              const IntVector<DIM> D2_ghostbox_lower = D2_ghostbox.lower();
              const IntVector<DIM> D2_ghostbox_upper = D2_ghostbox.upper();
              Box<DIM> D3_ghostbox = D3_data->getGhostBox();
              const IntVector<DIM> D3_ghostbox_lower = D3_ghostbox.lower();
              const IntVector<DIM> D3_ghostbox_upper = D3_ghostbox.upper();

              double* D1 = D1_data->getPointer();
              double* D2 = D2_data->getPointer();
              double* D3 = D3_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_three_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();

              double* D1 = D1_data->getPointer();

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

    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computePlusAndMinusSpatialDerivatives() */
template <int DIM> 
void LevelSetMethodToolbox<DIM>::computePlusAndMinusSpatialDerivatives(
  Pointer< PatchHierarchy<DIM> > hierarchy,
  const SPATIAL_DERIVATIVE_TYPE spatial_derivative_type,
  const int spatial_derivative_order,
  const int grad_phi_plus_handle,
  const int grad_phi_minus_handle,
  const int phi_handle,
  const int phi_component)
{

  // make sure that the scratch PatchData handles have been created
  initializeComputeSpatialDerivativesParameters();

  const int finest_level = hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    Pointer< PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computePlusAndMinusSpatialDerivatives(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // compute spatial derivatives for phi
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
  
      Pointer< CellData<DIM,double> > grad_phi_plus_data =
        patch->getPatchData( grad_phi_plus_handle );
      Pointer< CellData<DIM,double> > grad_phi_minus_data =
        patch->getPatchData( grad_phi_minus_handle );
      Pointer< CellData<DIM,double> > phi_data =
        patch->getPatchData( phi_handle );
  
      Box<DIM> fillbox = grad_phi_plus_data->getBox();
      const IntVector<DIM> grad_phi_fillbox_lower = fillbox.lower();
      const IntVector<DIM> grad_phi_fillbox_upper = fillbox.upper();

      Box<DIM> grad_phi_plus_ghostbox = grad_phi_plus_data->getGhostBox();
      const IntVector<DIM> grad_phi_plus_ghostbox_lower = 
        grad_phi_plus_ghostbox.lower();
      const IntVector<DIM> grad_phi_plus_ghostbox_upper = 
        grad_phi_plus_ghostbox.upper();

      Box<DIM> grad_phi_minus_ghostbox = grad_phi_minus_data->getGhostBox();
      const IntVector<DIM> grad_phi_minus_ghostbox_lower = 
        grad_phi_minus_ghostbox.lower();
      const IntVector<DIM> grad_phi_minus_ghostbox_upper = 
        grad_phi_minus_ghostbox.upper();

      Box<DIM> phi_ghostbox = phi_data->getGhostBox();
      const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
      const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();

      double* grad_phi_plus[LSM_DIM_MAX];
      double* grad_phi_minus[LSM_DIM_MAX];
      double* phi = phi_data->getPointer(phi_component);
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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_one_ghostcell_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();

              double* D1 = D1_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_two_ghostcells_handle );
              Pointer< CellData<DIM,double> > D2_data =
                patch->getPatchData( s_D2_two_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();
              Box<DIM> D2_ghostbox = D2_data->getGhostBox();
              const IntVector<DIM> D2_ghostbox_lower = D2_ghostbox.lower();
              const IntVector<DIM> D2_ghostbox_upper = D2_ghostbox.upper();

              double* D1 = D1_data->getPointer();
              double* D2 = D2_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_three_ghostcells_handle );
              Pointer< CellData<DIM,double> > D2_data =
                patch->getPatchData( s_D2_three_ghostcells_handle );
              Pointer< CellData<DIM,double> > D3_data =
                patch->getPatchData( s_D3_three_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();
              Box<DIM> D2_ghostbox = D2_data->getGhostBox();
              const IntVector<DIM> D2_ghostbox_lower = D2_ghostbox.lower();
              const IntVector<DIM> D2_ghostbox_upper = D2_ghostbox.upper();
              Box<DIM> D3_ghostbox = D3_data->getGhostBox();
              const IntVector<DIM> D3_ghostbox_lower = D3_ghostbox.lower();
              const IntVector<DIM> D3_ghostbox_upper = D3_ghostbox.upper();

              double* D1 = D1_data->getPointer();
              double* D2 = D2_data->getPointer();
              double* D3 = D3_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_three_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();

              double* D1 = D1_data->getPointer();

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

    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computeCentralSpatialDerivatives() */
template <int DIM> 
void LevelSetMethodToolbox<DIM>::computeCentralSpatialDerivatives(
  Pointer< PatchHierarchy<DIM> > hierarchy,
  const int spatial_derivative_order,
  const int grad_phi_handle,
  const int phi_handle,
  const int phi_component)
{

  // make sure that the scratch PatchData handles have been created
  initializeComputeSpatialDerivativesParameters();

  const int finest_level = hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    Pointer< PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeCentralSpatialDerivatives(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // compute spatial derivatives for phi
      Pointer< CellData<DIM,double> > grad_phi_data =
        patch->getPatchData( grad_phi_handle );
      Pointer< CellData<DIM,double> > phi_data =
        patch->getPatchData( phi_handle );
  
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
  
      Box<DIM> fillbox = grad_phi_data->getBox();
      const IntVector<DIM> grad_phi_fillbox_lower = fillbox.lower();
      const IntVector<DIM> grad_phi_fillbox_upper = fillbox.upper();

      Box<DIM> grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const IntVector<DIM> grad_phi_ghostbox_lower = 
        grad_phi_ghostbox.lower();
      const IntVector<DIM> grad_phi_ghostbox_upper = 
        grad_phi_ghostbox.upper();

      Box<DIM> phi_ghostbox = phi_data->getGhostBox();
      const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
      const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();

      double* grad_phi[LSM_DIM_MAX];
      double* phi = phi_data->getPointer(phi_component);
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

    } // end loop over Patches
  } // end loop over PatchLevels
}


/* TVDRK1Step() */
template <int DIM> 
void LevelSetMethodToolbox<DIM>::TVDRK1Step(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int u_next_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const double dt,
  const int u_next_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "TVDRK1Step(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,double> > u_next_data =
        patch->getPatchData( u_next_handle );
      Pointer< CellData<DIM,double> > u_cur_data =
        patch->getPatchData( u_cur_handle );
      Pointer< CellData<DIM,double> > rhs_data =
        patch->getPatchData( rhs_handle );
  
      Box<DIM> u_next_ghostbox = u_next_data->getGhostBox();
      const IntVector<DIM> u_next_ghostbox_lower = u_next_ghostbox.lower();
      const IntVector<DIM> u_next_ghostbox_upper = u_next_ghostbox.upper();

      Box<DIM> u_cur_ghostbox = u_cur_data->getGhostBox();
      const IntVector<DIM> u_cur_ghostbox_lower = u_cur_ghostbox.lower();
      const IntVector<DIM> u_cur_ghostbox_upper = u_cur_ghostbox.upper();

      Box<DIM> rhs_ghostbox = rhs_data->getGhostBox();
      const IntVector<DIM> rhs_ghostbox_lower = rhs_ghostbox.lower();
      const IntVector<DIM> rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      Box<DIM> fillbox = rhs_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      double* u_next = u_next_data->getPointer(u_next_component);
      double* u_cur = u_cur_data->getPointer(u_cur_component);
      double* rhs = rhs_data->getPointer(rhs_component);

      if ( DIM == 3 ) {
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
template <int DIM> 
void LevelSetMethodToolbox<DIM>::TVDRK2Stage1(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int u_stage1_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const double dt,
  const int u_stage1_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "TVDRK2Stage1(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,double> > u_stage1_data =
        patch->getPatchData( u_stage1_handle );
      Pointer< CellData<DIM,double> > u_cur_data =
        patch->getPatchData( u_cur_handle );
      Pointer< CellData<DIM,double> > rhs_data =
        patch->getPatchData( rhs_handle );
  
      Box<DIM> u_stage1_ghostbox = u_stage1_data->getGhostBox();
      const IntVector<DIM> u_stage1_ghostbox_lower = 
        u_stage1_ghostbox.lower();
      const IntVector<DIM> u_stage1_ghostbox_upper = 
        u_stage1_ghostbox.upper();

      Box<DIM> u_cur_ghostbox = u_cur_data->getGhostBox();
      const IntVector<DIM> u_cur_ghostbox_lower = 
        u_cur_ghostbox.lower();
      const IntVector<DIM> u_cur_ghostbox_upper = 
        u_cur_ghostbox.upper();

      Box<DIM> rhs_ghostbox = rhs_data->getGhostBox();
      const IntVector<DIM> rhs_ghostbox_lower = rhs_ghostbox.lower();
      const IntVector<DIM> rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      Box<DIM> fillbox = u_stage1_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      double* u_stage1 = u_stage1_data->getPointer(u_stage1_component);
      double* u_cur = u_cur_data->getPointer(u_cur_component);
      double* rhs = rhs_data->getPointer(rhs_component);

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
template <int DIM> 
void LevelSetMethodToolbox<DIM>::TVDRK2Stage2(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int u_next_handle,
  const int u_stage1_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const double dt,
  const int u_next_component,
  const int u_stage1_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "TVDRK2Stage2(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,double> > u_next_data =
        patch->getPatchData( u_next_handle );
      Pointer< CellData<DIM,double> > u_stage1_data =
        patch->getPatchData( u_stage1_handle );
      Pointer< CellData<DIM,double> > u_cur_data =
        patch->getPatchData( u_cur_handle );
      Pointer< CellData<DIM,double> > rhs_data =
        patch->getPatchData( rhs_handle );
  
      Box<DIM> u_next_ghostbox = u_next_data->getGhostBox();
      const IntVector<DIM> u_next_ghostbox_lower = 
        u_next_ghostbox.lower();
      const IntVector<DIM> u_next_ghostbox_upper = 
        u_next_ghostbox.upper();

      Box<DIM> u_stage1_ghostbox = u_stage1_data->getGhostBox();
      const IntVector<DIM> u_stage1_ghostbox_lower = 
        u_stage1_ghostbox.lower();
      const IntVector<DIM> u_stage1_ghostbox_upper = 
        u_stage1_ghostbox.upper();

      Box<DIM> u_cur_ghostbox = u_cur_data->getGhostBox();
      const IntVector<DIM> u_cur_ghostbox_lower = 
        u_cur_ghostbox.lower();
      const IntVector<DIM> u_cur_ghostbox_upper = 
        u_cur_ghostbox.upper();

      Box<DIM> rhs_ghostbox = rhs_data->getGhostBox();
      const IntVector<DIM> rhs_ghostbox_lower = rhs_ghostbox.lower();
      const IntVector<DIM> rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      Box<DIM> fillbox = u_next_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      double* u_next = u_next_data->getPointer(u_next_component);
      double* u_stage1 = u_stage1_data->getPointer(u_stage1_component);
      double* u_cur = u_cur_data->getPointer(u_cur_component);
      double* rhs = rhs_data->getPointer(rhs_component);

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
template <int DIM> 
void LevelSetMethodToolbox<DIM>::TVDRK3Stage1(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int u_stage1_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const double dt,
  const int u_stage1_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "TVDRK3Stage1(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,double> > u_stage1_data =
        patch->getPatchData( u_stage1_handle );
      Pointer< CellData<DIM,double> > u_cur_data =
        patch->getPatchData( u_cur_handle );
      Pointer< CellData<DIM,double> > rhs_data =
        patch->getPatchData( rhs_handle );
  
      Box<DIM> u_stage1_ghostbox = u_stage1_data->getGhostBox();
      const IntVector<DIM> u_stage1_ghostbox_lower = 
        u_stage1_ghostbox.lower();
      const IntVector<DIM> u_stage1_ghostbox_upper = 
        u_stage1_ghostbox.upper();

      Box<DIM> u_cur_ghostbox = u_cur_data->getGhostBox();
      const IntVector<DIM> u_cur_ghostbox_lower = 
        u_cur_ghostbox.lower();
      const IntVector<DIM> u_cur_ghostbox_upper = 
        u_cur_ghostbox.upper();

      Box<DIM> rhs_ghostbox = rhs_data->getGhostBox();
      const IntVector<DIM> rhs_ghostbox_lower = rhs_ghostbox.lower();
      const IntVector<DIM> rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      Box<DIM> fillbox = u_stage1_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      double* u_stage1 = u_stage1_data->getPointer(u_stage1_component);
      double* u_cur = u_cur_data->getPointer(u_cur_component);
      double* rhs = rhs_data->getPointer(rhs_component);

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
template <int DIM> 
void LevelSetMethodToolbox<DIM>::TVDRK3Stage2(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int u_stage2_handle,
  const int u_stage1_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const double dt,
  const int u_stage2_component,
  const int u_stage1_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "TVDRK3Stage2(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,double> > u_stage2_data =
        patch->getPatchData( u_stage2_handle );
      Pointer< CellData<DIM,double> > u_stage1_data =
        patch->getPatchData( u_stage1_handle );
      Pointer< CellData<DIM,double> > u_cur_data =
        patch->getPatchData( u_cur_handle );
      Pointer< CellData<DIM,double> > rhs_data =
        patch->getPatchData( rhs_handle );
  
      Box<DIM> u_stage2_ghostbox = u_stage2_data->getGhostBox();
      const IntVector<DIM> u_stage2_ghostbox_lower = 
        u_stage2_ghostbox.lower();
      const IntVector<DIM> u_stage2_ghostbox_upper = 
        u_stage2_ghostbox.upper();

      Box<DIM> u_stage1_ghostbox = u_stage1_data->getGhostBox();
      const IntVector<DIM> u_stage1_ghostbox_lower = 
        u_stage1_ghostbox.lower();
      const IntVector<DIM> u_stage1_ghostbox_upper = 
        u_stage1_ghostbox.upper();

      Box<DIM> u_cur_ghostbox = u_cur_data->getGhostBox();
      const IntVector<DIM> u_cur_ghostbox_lower = 
        u_cur_ghostbox.lower();
      const IntVector<DIM> u_cur_ghostbox_upper = 
        u_cur_ghostbox.upper();

      Box<DIM> rhs_ghostbox = rhs_data->getGhostBox();
      const IntVector<DIM> rhs_ghostbox_lower = rhs_ghostbox.lower();
      const IntVector<DIM> rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      Box<DIM> fillbox = u_stage2_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      double* u_stage2 = u_stage2_data->getPointer(u_stage2_component);
      double* u_stage1 = u_stage1_data->getPointer(u_stage1_component);
      double* u_cur = u_cur_data->getPointer(u_cur_component);
      double* rhs = rhs_data->getPointer(rhs_component);

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
template <int DIM> 
void LevelSetMethodToolbox<DIM>::TVDRK3Stage3(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int u_next_handle,
  const int u_stage2_handle,
  const int u_cur_handle,
  const int rhs_handle,
  const double dt,
  const int u_next_component,
  const int u_stage2_component,
  const int u_cur_component,
  const int rhs_component)
{
  // loop over PatchHierarchy and take Runge-Kutta step
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "TVDRK3Stage3(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,double> > u_next_data =
        patch->getPatchData( u_next_handle );
      Pointer< CellData<DIM,double> > u_stage2_data =
        patch->getPatchData( u_stage2_handle );
      Pointer< CellData<DIM,double> > u_cur_data =
        patch->getPatchData( u_cur_handle );
      Pointer< CellData<DIM,double> > rhs_data =
        patch->getPatchData( rhs_handle );
  
      Box<DIM> u_next_ghostbox = u_next_data->getGhostBox();
      const IntVector<DIM> u_next_ghostbox_lower = 
        u_next_ghostbox.lower();
      const IntVector<DIM> u_next_ghostbox_upper = 
        u_next_ghostbox.upper();

      Box<DIM> u_stage2_ghostbox = u_stage2_data->getGhostBox();
      const IntVector<DIM> u_stage2_ghostbox_lower = 
        u_stage2_ghostbox.lower();
      const IntVector<DIM> u_stage2_ghostbox_upper = 
        u_stage2_ghostbox.upper();

      Box<DIM> u_cur_ghostbox = u_cur_data->getGhostBox();
      const IntVector<DIM> u_cur_ghostbox_lower = 
        u_cur_ghostbox.lower();
      const IntVector<DIM> u_cur_ghostbox_upper = 
        u_cur_ghostbox.upper();

      Box<DIM> rhs_ghostbox = rhs_data->getGhostBox();
      const IntVector<DIM> rhs_ghostbox_lower = rhs_ghostbox.lower();
      const IntVector<DIM> rhs_ghostbox_upper = rhs_ghostbox.upper();

      // fill box
      Box<DIM> fillbox = u_stage2_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      double* u_next = u_next_data->getPointer(u_next_component);
      double* u_stage2 = u_stage2_data->getPointer(u_stage2_component);
      double* u_cur = u_cur_data->getPointer(u_cur_component);
      double* rhs = rhs_data->getPointer(rhs_component);

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


/* computeDistanceFunctionUsingFMM() */
template <int DIM> 
void LevelSetMethodToolbox<DIM>::computeDistanceFunctionUsingFMM(
  Pointer< PatchHierarchy<DIM> > hierarchy,
  const int spatial_derivative_order,
  const int distance_function_handle,
  const int phi_handle,
  const int distance_function_component,
  const int phi_component)
{
  const int finest_level = hierarchy->getFinestLevelNumber();

  /*
   * error checking
   */

  // TEMP - only single level calculations currently supported
  if (finest_level > 0) {
    TBOX_ERROR(  "LevelSetMethodToolbox::"
              << "computeDistanceFunctionUsingFMM(): "
              << "Only single level calculations currently supported."
              << endl );
  }

  // TEMP - only first-order calculations currently supported
  if (spatial_derivative_order > 0) {
    TBOX_WARNING(  "LevelSetMethodToolbox::"
                << "computeDistanceFunctionUsingFMM(): "
                << "Only first-order calculations currently supported.  "
                << "Dropping to first-order calculations." 
                << endl );
  }

  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    Pointer< PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);

    // error checking...only single level calculations currently supported
    if (level->getNumberOfPatches() > 1) {
      TBOX_ERROR(  "LevelSetMethodToolbox::"
                << "computeDistanceFunctionUsingFMM(): "
                << "Only single patch calculations currently supported."
                << endl );
    }

    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeDistanceFunctionUsingFMM(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // get geometry information for patch
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
  
      // get PatchData for distance function and phi
      Pointer< CellData<DIM,double> > distance_function_data =
        patch->getPatchData( distance_function_handle );
      Pointer< CellData<DIM,double> > phi_data =
        patch->getPatchData( phi_handle );
  
      // get index space information for PatchData
      Box<DIM> ghostbox = distance_function_data->getGhostBox();
      IntVector<DIM> grid_dims_int_vector = 
        ghostbox.upper() - ghostbox.lower();
      int* grid_dims = grid_dims_int_vector;

      double* distance_function = 
        distance_function_data->getPointer(distance_function_component);
      double* phi = phi_data->getPointer(phi_component);

      // call computeExtensionFields*() from toolbox to carry out
      // computation
      if ( DIM == 3 ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeDistanceFunctionUsingFMM(): "
                  << "DIM = 3 is currently unsupported."
                  << endl );
        /*
        computeExtensionFields3d(
          distance_function,
          (double**) 0,
          phi,
          (double**) 0,
          0,
          grid_dims,
          (double*) dx);
        */
      } else if ( DIM == 2 ) {
        computeExtensionFields2d(
          distance_function,
          (double**) 0,
          phi,
          (double*) 0, // NULL mask field
          (double**) 0,
          0,
          spatial_derivative_order,  
          grid_dims,
          (double*) dx);
      } else if ( DIM == 1 ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeDistanceFunctionUsingFMM(): "
                  << "DIM = 1 is currently unsupported."
                  << endl );
      } else {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeDistanceFunctionUsingFMM(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 2 is currently supported."
                  << endl );

      } // end switch on DIM

    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computeExtensionFieldsUsingFMM() */
template <int DIM>
void LevelSetMethodToolbox<DIM>::computeExtensionFieldsUsingFMM(
  Pointer< PatchHierarchy<DIM> > hierarchy,
  const int spatial_derivative_order,
  const vector<int>& extension_field_handles,
  const int distance_function_handle,
  const vector<int>& source_field_handles,
  const int phi_handle,
  const int extension_field_component,
  const int distance_function_component,
  const int source_field_component,
  const int phi_component)
{
  const int finest_level = hierarchy->getFinestLevelNumber();

  /*
   * error checking
   */

  // TEMP - only single level calculations currently supported
  if (finest_level > 0) {
    TBOX_ERROR(  "LevelSetMethodToolbox::"
              << "computeExtensionFieldsUsingFMM(): "
              << "Only single level calculations currently supported."
              << endl );
  }

  // TEMP - only first-order calculations currently supported
  if (spatial_derivative_order > 0) {
    TBOX_WARNING(  "LevelSetMethodToolbox::"
                << "computeExtensionFieldsUsingFMM(): "
                << "Only first-order calculations currently supported.  "
                << "Dropping to first-order calculations." 
                << endl );
  }

  // allocate memory for extension fields and source fields
  int num_extension_fields = extension_field_handles.size();
  double** extension_fields = new double*[num_extension_fields];
  double** source_fields = new double*[num_extension_fields];

  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    Pointer< PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);

    // error checking...only single level calculations currently supported
    if (level->getNumberOfPatches() > 1) {
      TBOX_ERROR(  "LevelSetMethodToolbox::"
                << "computeExtensionFieldsUsingFMM(): "
                << "Only single patch calculations currently supported."
                << endl );
    }

    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeExtensionFieldsUsingFMM(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // get geometry information for patch
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
  
      // get PatchData for distance function and phi
      Pointer< CellData<DIM,double> > distance_function_data =
        patch->getPatchData( distance_function_handle );
      Pointer< CellData<DIM,double> > phi_data =
        patch->getPatchData( phi_handle );
      double* distance_function = 
        distance_function_data->getPointer(distance_function_component);
      double* phi = phi_data->getPointer(phi_component);
  
      // get PatchData for source fields and extension fields
      for (int k=0; k < num_extension_fields; k++) {
        Pointer< CellData<DIM,double> > extension_field_data = 
          patch->getPatchData( extension_field_handles[k] );
        Pointer< CellData<DIM,double> > source_field_data = 
          patch->getPatchData( source_field_handles[k] );

        extension_fields[k] = 
          extension_field_data->getPointer(extension_field_component);
        source_fields[k] = 
          source_field_data->getPointer(source_field_component);
      }

      // get index space information for PatchData
      Box<DIM> ghostbox = distance_function_data->getGhostBox();
      IntVector<DIM> grid_dims_int_vector = 
        ghostbox.upper() - ghostbox.lower();
      int* grid_dims = grid_dims_int_vector;

      // call computeExtensionFields*() from toolbox to carry out
      // computation
      if ( DIM == 3 ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeExtensionFieldsUsingFMM(): "
                  << "DIM = 3 is currently unsupported."
                  << endl );
        /*
        computeExtensionFields3d(
          distance_function,
          ext_field
          phi,
          (double*) 0,
          (double**) 0,
          0,
          grid_dims,
          (double*) dx);
        */
      } else if ( DIM == 2 ) {
        computeExtensionFields2d(
          distance_function,
          extension_fields,
          phi,
          (double*) 0,  // NULL mask field
          source_fields,
          num_extension_fields,
          spatial_derivative_order,
          grid_dims,
          (double*) dx);
      } else if ( DIM == 1 ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeExtensionFieldsUsingFMM(): "
                  << "DIM = 1 is currently unsupported."
                  << endl );
      } else {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeExtensionFieldsUsingFMM(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 2 is currently supported."
                  << endl );

      } // end switch on DIM

    } // end loop over Patches
  } // end loop over PatchLevels

  // free memory for extension fields and source fields
  delete [] extension_fields; 
  delete [] source_fields; 
}


/* computeUnitNormalVectorFromPhi() */
template <int DIM>
void LevelSetMethodToolbox<DIM>::computeUnitNormalVectorFromPhi(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
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
    initializeComputeUnitNormalParameters();
  }
  initializeComputeSpatialDerivativesParameters();

  const int finest_level = patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
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
      Pointer< CellData<DIM,double> > normal_vector_data=
        patch->getPatchData( normal_vector_handle );
      Pointer< CellData<DIM,double> > phi_data =
        patch->getPatchData( phi_handle );
      Pointer< CellData<DIM,double> > grad_phi_data =
        patch->getPatchData( s_compute_normal_grad_phi_handle );
      Pointer< CellData<DIM,double> > grad_phi_plus_data =
        patch->getPatchData( s_compute_normal_grad_phi_plus_handle );
      Pointer< CellData<DIM,double> > grad_phi_minus_data =
        patch->getPatchData( s_compute_normal_grad_phi_minus_handle );
  
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
  
      Box<DIM> fillbox = normal_vector_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      Box<DIM> normal_vector_ghostbox = normal_vector_data->getGhostBox();
      const IntVector<DIM> normal_vector_ghostbox_lower = 
        normal_vector_ghostbox.lower();
      const IntVector<DIM> normal_vector_ghostbox_upper = 
        normal_vector_ghostbox.upper();

      Box<DIM> phi_ghostbox = phi_data->getGhostBox();
      const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
      const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();

      Box<DIM> grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const IntVector<DIM> grad_phi_ghostbox_lower = 
        grad_phi_ghostbox.lower();
      const IntVector<DIM> grad_phi_ghostbox_upper = 
        grad_phi_ghostbox.upper();

      Box<DIM> grad_phi_plus_ghostbox = grad_phi_plus_data->getGhostBox();
      const IntVector<DIM> grad_phi_plus_ghostbox_lower = 
        grad_phi_plus_ghostbox.lower();
      const IntVector<DIM> grad_phi_plus_ghostbox_upper = 
        grad_phi_plus_ghostbox.upper();

      Box<DIM> grad_phi_minus_ghostbox = grad_phi_minus_data->getGhostBox();
      const IntVector<DIM> grad_phi_minus_ghostbox_lower = 
        grad_phi_minus_ghostbox.lower();
      const IntVector<DIM> grad_phi_minus_ghostbox_upper = 
        grad_phi_minus_ghostbox.upper();

      double* normal_vector[LSM_DIM_MAX];
      double* phi = phi_data->getPointer(phi_component);
      double* grad_phi[LSM_DIM_MAX];
      double* grad_phi_plus[LSM_DIM_MAX];
      double* grad_phi_minus[LSM_DIM_MAX];
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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_one_ghostcell_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();

              double* D1 = D1_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_two_ghostcells_handle );
              Pointer< CellData<DIM,double> > D2_data =
                patch->getPatchData( s_D2_two_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();
              Box<DIM> D2_ghostbox = D2_data->getGhostBox();
              const IntVector<DIM> D2_ghostbox_lower = D2_ghostbox.lower();
              const IntVector<DIM> D2_ghostbox_upper = D2_ghostbox.upper();

              double* D1 = D1_data->getPointer();
              double* D2 = D2_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_three_ghostcells_handle );
              Pointer< CellData<DIM,double> > D2_data =
                patch->getPatchData( s_D2_three_ghostcells_handle );
              Pointer< CellData<DIM,double> > D3_data =
                patch->getPatchData( s_D3_three_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();
              Box<DIM> D2_ghostbox = D2_data->getGhostBox();
              const IntVector<DIM> D2_ghostbox_lower = D2_ghostbox.lower();
              const IntVector<DIM> D2_ghostbox_upper = D2_ghostbox.upper();
              Box<DIM> D3_ghostbox = D3_data->getGhostBox();
              const IntVector<DIM> D3_ghostbox_lower = D3_ghostbox.lower();
              const IntVector<DIM> D3_ghostbox_upper = D3_ghostbox.upper();

              double* D1 = D1_data->getPointer();
              double* D2 = D2_data->getPointer();
              double* D3 = D3_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_three_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();

              double* D1 = D1_data->getPointer();

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

    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computeSignedUnitNormalVectorFromPhi() */
template <int DIM>
void LevelSetMethodToolbox<DIM>::computeSignedUnitNormalVectorFromPhi(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
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
    initializeComputeUnitNormalParameters();
  }
  initializeComputeSpatialDerivativesParameters();

  const int finest_level = patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
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
      Pointer< CellData<DIM,double> > normal_vector_data=
        patch->getPatchData( normal_vector_handle );
      Pointer< CellData<DIM,double> > phi_data =
        patch->getPatchData( phi_handle );
      Pointer< CellData<DIM,double> > grad_phi_data =
        patch->getPatchData( s_compute_normal_grad_phi_handle );
      Pointer< CellData<DIM,double> > grad_phi_plus_data =
        patch->getPatchData( s_compute_normal_grad_phi_plus_handle );
      Pointer< CellData<DIM,double> > grad_phi_minus_data =
        patch->getPatchData( s_compute_normal_grad_phi_minus_handle );
  
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
  
      Box<DIM> fillbox = normal_vector_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      Box<DIM> normal_vector_ghostbox = normal_vector_data->getGhostBox();
      const IntVector<DIM> normal_vector_ghostbox_lower = 
        normal_vector_ghostbox.lower();
      const IntVector<DIM> normal_vector_ghostbox_upper = 
        normal_vector_ghostbox.upper();

      Box<DIM> phi_ghostbox = phi_data->getGhostBox();
      const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
      const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();

      Box<DIM> grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const IntVector<DIM> grad_phi_ghostbox_lower = 
        grad_phi_ghostbox.lower();
      const IntVector<DIM> grad_phi_ghostbox_upper = 
        grad_phi_ghostbox.upper();

      Box<DIM> grad_phi_plus_ghostbox = grad_phi_plus_data->getGhostBox();
      const IntVector<DIM> grad_phi_plus_ghostbox_lower = 
        grad_phi_plus_ghostbox.lower();
      const IntVector<DIM> grad_phi_plus_ghostbox_upper = 
        grad_phi_plus_ghostbox.upper();

      Box<DIM> grad_phi_minus_ghostbox = grad_phi_minus_data->getGhostBox();
      const IntVector<DIM> grad_phi_minus_ghostbox_lower = 
        grad_phi_minus_ghostbox.lower();
      const IntVector<DIM> grad_phi_minus_ghostbox_upper = 
        grad_phi_minus_ghostbox.upper();

      double* normal_vector[LSM_DIM_MAX];
      double* phi = phi_data->getPointer(phi_component);
      double* grad_phi[LSM_DIM_MAX];
      double* grad_phi_plus[LSM_DIM_MAX];
      double* grad_phi_minus[LSM_DIM_MAX];
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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_one_ghostcell_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();

              double* D1 = D1_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_two_ghostcells_handle );
              Pointer< CellData<DIM,double> > D2_data =
                patch->getPatchData( s_D2_two_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();
              Box<DIM> D2_ghostbox = D2_data->getGhostBox();
              const IntVector<DIM> D2_ghostbox_lower = D2_ghostbox.lower();
              const IntVector<DIM> D2_ghostbox_upper = D2_ghostbox.upper();

              double* D1 = D1_data->getPointer();
              double* D2 = D2_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_three_ghostcells_handle );
              Pointer< CellData<DIM,double> > D2_data =
                patch->getPatchData( s_D2_three_ghostcells_handle );
              Pointer< CellData<DIM,double> > D3_data =
                patch->getPatchData( s_D3_three_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();
              Box<DIM> D2_ghostbox = D2_data->getGhostBox();
              const IntVector<DIM> D2_ghostbox_lower = D2_ghostbox.lower();
              const IntVector<DIM> D2_ghostbox_upper = D2_ghostbox.upper();
              Box<DIM> D3_ghostbox = D3_data->getGhostBox();
              const IntVector<DIM> D3_ghostbox_lower = D3_ghostbox.lower();
              const IntVector<DIM> D3_ghostbox_upper = D3_ghostbox.upper();

              double* D1 = D1_data->getPointer();
              double* D2 = D2_data->getPointer();
              double* D3 = D3_data->getPointer();

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

              Pointer< CellData<DIM,double> > D1_data =
                patch->getPatchData( s_D1_three_ghostcells_handle );

              Box<DIM> D1_ghostbox = D1_data->getGhostBox();
              const IntVector<DIM> D1_ghostbox_lower = D1_ghostbox.lower();
              const IntVector<DIM> D1_ghostbox_upper = D1_ghostbox.upper();

              double* D1 = D1_data->getPointer();

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

    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computeUnitNormalVectorFromGradPhi() */
template <int DIM>
void LevelSetMethodToolbox<DIM>::computeUnitNormalVectorFromGradPhi(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int normal_vector_handle,
  const int grad_phi_handle)
{
  const int finest_level = patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeUnitNormalVectorFromGradPhi(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // get PatchData
      Pointer< CellData<DIM,double> > normal_vector_data=
        patch->getPatchData( normal_vector_handle );
      Pointer< CellData<DIM,double> > grad_phi_data =
        patch->getPatchData( grad_phi_handle );
  
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
  
      Box<DIM> fillbox = normal_vector_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      Box<DIM> normal_vector_ghostbox = normal_vector_data->getGhostBox();
      const IntVector<DIM> normal_vector_ghostbox_lower = 
        normal_vector_ghostbox.lower();
      const IntVector<DIM> normal_vector_ghostbox_upper = 
        normal_vector_ghostbox.upper();

      Box<DIM> grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const IntVector<DIM> grad_phi_ghostbox_lower = 
        grad_phi_ghostbox.lower();
      const IntVector<DIM> grad_phi_ghostbox_upper = 
        grad_phi_ghostbox.upper();

      double* normal_vector[LSM_DIM_MAX];
      double* grad_phi[LSM_DIM_MAX];
      for (int dim = 0; dim < DIM; dim++) {
        normal_vector[dim] = normal_vector_data->getPointer(dim);
        grad_phi[dim] = grad_phi_data->getPointer(dim);
      }

      // compute unit normal from grad(phi)
      if ( DIM == 3 ) {
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
template <int DIM>
void LevelSetMethodToolbox<DIM>::computeSignedUnitNormalVectorFromGradPhi(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int normal_vector_handle,
  const int grad_phi_handle,
  const int phi_handle,
  const int phi_component)
{
  const int finest_level = patch_hierarchy->getFinestLevelNumber();
  for ( int ln=0 ; ln<=finest_level ; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeSignedUnitNormalVectorFromGradPhi(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      // get PatchData
      Pointer< CellData<DIM,double> > normal_vector_data=
        patch->getPatchData( normal_vector_handle );
      Pointer< CellData<DIM,double> > grad_phi_data =
        patch->getPatchData( grad_phi_handle );
      Pointer< CellData<DIM,double> > phi_data =
        patch->getPatchData( phi_handle );
  
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
  
      Box<DIM> fillbox = normal_vector_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      Box<DIM> normal_vector_ghostbox = normal_vector_data->getGhostBox();
      const IntVector<DIM> normal_vector_ghostbox_lower = 
        normal_vector_ghostbox.lower();
      const IntVector<DIM> normal_vector_ghostbox_upper = 
        normal_vector_ghostbox.upper();

      Box<DIM> grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const IntVector<DIM> grad_phi_ghostbox_lower = 
        grad_phi_ghostbox.lower();
      const IntVector<DIM> grad_phi_ghostbox_upper = 
        grad_phi_ghostbox.upper();

      Box<DIM> phi_ghostbox = phi_data->getGhostBox();
      const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
      const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();

      double* normal_vector[LSM_DIM_MAX];
      double* grad_phi[LSM_DIM_MAX];
      double* phi = phi_data->getPointer(phi_component);
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

    } // end loop over Patches
  } // end loop over PatchLevels
}


/* computeVolumeOfRegionDefinedByZeroLevelSet() */
template <int DIM> 
double LevelSetMethodToolbox<DIM>::computeVolumeOfRegionDefinedByZeroLevelSet(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int phi_handle,
  const int control_volume_handle,
  const int region_indicator,
  const int phi_component,
  const int heaviside_width)
{
  double volume = 0.0;

  // loop over PatchHierarchy and compute the integral on each Patch
  // by calling Fortran subroutines
  const int num_levels = patch_hierarchy->getNumberLevels();
  if (region_indicator > 0) { // integrate over region {x | phi(x) > 0}

    for ( int ln=0 ; ln < num_levels; ln++ ) {

      Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
      
      typename PatchLevel<DIM>::Iterator pi;
      for (pi.initialize(level); pi; pi++) { // loop over patches
        const int pn = *pi;
        Pointer< Patch<DIM> > patch = level->getPatch(pn);
        if ( patch.isNull() ) {
          TBOX_ERROR(  "LevelSetMethodToolbox::" 
                    << "computeVolumeOfRegionDefinedByZeroLevelSet(): "
                    << "Cannot find patch. Null patch pointer."
                    << endl);
        }

        // get dx and epsilon
        Pointer< CartesianPatchGeometry<DIM> > patch_geom =
          patch->getPatchGeometry();
        const double* dx = patch_geom->getDx();
        double max_dx = dx[0];
        for (int k = 1; k < DIM; k++) {
          if (max_dx < dx[k]) max_dx = dx[k];
        }
        double epsilon = heaviside_width*max_dx;

        // get pointers to data and index space ranges
        Pointer< CellData<DIM,double> > phi_data =
          patch->getPatchData( phi_handle );
        Pointer< CellData<DIM,double> > control_volume_data =
          patch->getPatchData( control_volume_handle );
    
        Box<DIM> phi_ghostbox = phi_data->getGhostBox();
        const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
        const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();
  
        Box<DIM> control_volume_ghostbox = control_volume_data->getGhostBox();
        const IntVector<DIM> control_volume_ghostbox_lower = 
          control_volume_ghostbox.lower();
        const IntVector<DIM> control_volume_ghostbox_upper = 
          control_volume_ghostbox.upper();
  
        // interior box
        Box<DIM> interior_box = patch->getBox();
        const IntVector<DIM> interior_box_lower = interior_box.lower();
        const IntVector<DIM> interior_box_upper = interior_box.upper();

        double* phi = phi_data->getPointer(phi_component);
        double* control_volume = control_volume_data->getPointer();
        double volume_on_patch = 0.0;
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

      } // end loop over patches in level
    } // end loop over levels in hierarchy

  } else { // integrate over region {x | phi(x) < 0}

    for ( int ln=0 ; ln < num_levels; ln++ ) {

      Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
      
      typename PatchLevel<DIM>::Iterator pi;
      for (pi.initialize(level); pi; pi++) { // loop over patches
        const int pn = *pi;
        Pointer< Patch<DIM> > patch = level->getPatch(pn);
        if ( patch.isNull() ) {
          TBOX_ERROR(  "LevelSetMethodToolbox::" 
                    << "computeVolumeOfRegionDefinedByZeroLevelSet(): "
                    << "Cannot find patch. Null patch pointer."
                    << endl);
        }

        // get dx and epsilon
        Pointer< CartesianPatchGeometry<DIM> > patch_geom =
          patch->getPatchGeometry();
        const double* dx = patch_geom->getDx();
        double max_dx = dx[0];
        for (int k = 1; k < DIM; k++) {
          if (max_dx < dx[k]) max_dx = dx[k];
        }
        double epsilon = heaviside_width*max_dx;

        // get pointers to data and index space ranges
        Pointer< CellData<DIM,double> > phi_data =
          patch->getPatchData( phi_handle );
        Pointer< CellData<DIM,double> > control_volume_data =
          patch->getPatchData( control_volume_handle );
  
        Box<DIM> phi_ghostbox = phi_data->getGhostBox();
        const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
        const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();
  
        Box<DIM> control_volume_ghostbox = control_volume_data->getGhostBox();
        const IntVector<DIM> control_volume_ghostbox_lower = 
          control_volume_ghostbox.lower();
        const IntVector<DIM> control_volume_ghostbox_upper = 
          control_volume_ghostbox.upper();
  
        // interior box
        Box<DIM> interior_box = patch->getBox();
        const IntVector<DIM> interior_box_lower = interior_box.lower();
        const IntVector<DIM> interior_box_upper = interior_box.upper();

        double* phi = phi_data->getPointer();
        double* control_volume = control_volume_data->getPointer();
        double volume_on_patch = 0.0;
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

      } // end loop over patches in level
    } // end loop over levels in hierarchy

  } // end if statement on (region_indicator > 0)

  return tbox::MPI::sumReduction(volume);
}


/* computeVolumeOfZeroLevelSet() */
template <int DIM> 
double LevelSetMethodToolbox<DIM>::computeVolumeOfZeroLevelSet(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int phi_handle,
  const int grad_phi_handle,
  const int control_volume_handle,
  const int phi_component,
  const int delta_width)
{
  double volume = 0.0;

  // loop over PatchHierarchy and compute the integral on each Patch
  // by calling Fortran subroutines
  const int num_levels = patch_hierarchy->getNumberLevels();

  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "computeVolumeOfZeroLevelSet(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get dx and epsilon
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
      double max_dx = dx[0];
      for (int k = 1; k < DIM; k++) {
        if (max_dx < dx[k]) max_dx = dx[k];
      }
      double epsilon = delta_width*max_dx;

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,double> > phi_data =
        patch->getPatchData( phi_handle );
      Pointer< CellData<DIM,double> > grad_phi_data =
        patch->getPatchData( grad_phi_handle );
      Pointer< CellData<DIM,double> > control_volume_data =
        patch->getPatchData( control_volume_handle );
  
      Box<DIM> phi_ghostbox = phi_data->getGhostBox();
      const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
      const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();

      Box<DIM> grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const IntVector<DIM> grad_phi_ghostbox_lower = 
        grad_phi_ghostbox.lower();
      const IntVector<DIM> grad_phi_ghostbox_upper = 
        grad_phi_ghostbox.upper();

      Box<DIM> control_volume_ghostbox = control_volume_data->getGhostBox();
      const IntVector<DIM> control_volume_ghostbox_lower = 
        control_volume_ghostbox.lower();
      const IntVector<DIM> control_volume_ghostbox_upper = 
        control_volume_ghostbox.upper();

      // interior box
      Box<DIM> interior_box = patch->getBox();
      const IntVector<DIM> interior_box_lower = interior_box.lower();
      const IntVector<DIM> interior_box_upper = interior_box.upper();

      double* phi = phi_data->getPointer(phi_component);
      double* control_volume = control_volume_data->getPointer();
      double* grad_phi[LSM_DIM_MAX];
      for (int k=0; k<DIM; k++) {
        grad_phi[k] = grad_phi_data->getPointer(k);
      }
      double volume_on_patch = 0.0;
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

    } // end loop over patches in level
  } // end loop over levels in hierarchy

  return tbox::MPI::sumReduction(volume);
}


/* computeVolumeIntegral() */
template <int DIM> 
double LevelSetMethodToolbox<DIM>::computeVolumeIntegral(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int F_handle,
  const int phi_handle,
  const int control_volume_handle,
  const int region_indicator,
  const int F_component,
  const int phi_component,
  const int heaviside_width)
{
  double integral_F = 0.0;

  // loop over PatchHierarchy and compute the integral on each Patch
  // by calling Fortran subroutines
  const int num_levels = patch_hierarchy->getNumberLevels();
  if (region_indicator > 0) { // integrate over region {x | phi(x) > 0}

    for ( int ln=0 ; ln < num_levels; ln++ ) {

      Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
      
      typename PatchLevel<DIM>::Iterator pi;
      for (pi.initialize(level); pi; pi++) { // loop over patches
        const int pn = *pi;
        Pointer< Patch<DIM> > patch = level->getPatch(pn);
        if ( patch.isNull() ) {
          TBOX_ERROR(  "LevelSetMethodToolbox::" 
                    << "computeVolumeIntegral(): "
                    << "Cannot find patch. Null patch pointer."
                    << endl);
        }

        // get dx and epsilon
        Pointer< CartesianPatchGeometry<DIM> > patch_geom =
          patch->getPatchGeometry();
        const double* dx = patch_geom->getDx();
        double max_dx = dx[0];
        for (int k = 1; k < DIM; k++) {
          if (max_dx < dx[k]) max_dx = dx[k];
        }
        double epsilon = heaviside_width*max_dx;

        // get pointers to data and index space ranges
        Pointer< CellData<DIM,double> > F_data =
          patch->getPatchData( F_handle);
        Pointer< CellData<DIM,double> > phi_data =
          patch->getPatchData( phi_handle );
        Pointer< CellData<DIM,double> > control_volume_data =
          patch->getPatchData( control_volume_handle );
    
        Box<DIM> F_ghostbox = F_data->getGhostBox();
        const IntVector<DIM> F_ghostbox_lower = F_ghostbox.lower();
        const IntVector<DIM> F_ghostbox_upper = F_ghostbox.upper();
  
        Box<DIM> phi_ghostbox = phi_data->getGhostBox();
        const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
        const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();
  
        Box<DIM> control_volume_ghostbox = control_volume_data->getGhostBox();
        const IntVector<DIM> control_volume_ghostbox_lower = 
          control_volume_ghostbox.lower();
        const IntVector<DIM> control_volume_ghostbox_upper = 
          control_volume_ghostbox.upper();
  
        // interior box
        Box<DIM> interior_box = patch->getBox();
        const IntVector<DIM> interior_box_lower = interior_box.lower();
        const IntVector<DIM> interior_box_upper = interior_box.upper();

        double* F = F_data->getPointer(F_component);
        double* phi = phi_data->getPointer(phi_component);
        double* control_volume = control_volume_data->getPointer();
        double integral_F_on_patch = 0.0;
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

      } // end loop over patches in level
    } // end loop over levels in hierarchy

  } else {

    for ( int ln=0 ; ln < num_levels; ln++ ) {

      Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
      
      typename PatchLevel<DIM>::Iterator pi;
      for (pi.initialize(level); pi; pi++) { // loop over patches
        const int pn = *pi;
        Pointer< Patch<DIM> > patch = level->getPatch(pn);
        if ( patch.isNull() ) {
          TBOX_ERROR(  "LevelSetMethodToolbox::" 
                    << "computeVolumeIntegral(): "
                    << "Cannot find patch. Null patch pointer."
                    << endl);
        }

        // get dx and epsilon
        Pointer< CartesianPatchGeometry<DIM> > patch_geom =
          patch->getPatchGeometry();
        const double* dx = patch_geom->getDx();
        double max_dx = dx[0];
        for (int k = 1; k < DIM; k++) {
          if (max_dx < dx[k]) max_dx = dx[k];
        }
        double epsilon = heaviside_width*max_dx;

        // get pointers to data and index space ranges
        Pointer< CellData<DIM,double> > F_data =
          patch->getPatchData( F_handle);
        Pointer< CellData<DIM,double> > phi_data =
          patch->getPatchData( phi_handle );
        Pointer< CellData<DIM,double> > control_volume_data =
          patch->getPatchData( control_volume_handle );
    
        Box<DIM> F_ghostbox = F_data->getGhostBox();
        const IntVector<DIM> F_ghostbox_lower = F_ghostbox.lower();
        const IntVector<DIM> F_ghostbox_upper = F_ghostbox.upper();
  
        Box<DIM> phi_ghostbox = phi_data->getGhostBox();
        const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
        const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();
  
        Box<DIM> control_volume_ghostbox = control_volume_data->getGhostBox();
        const IntVector<DIM> control_volume_ghostbox_lower = 
          control_volume_ghostbox.lower();
        const IntVector<DIM> control_volume_ghostbox_upper = 
          control_volume_ghostbox.upper();
  
        // interior box
        Box<DIM> interior_box = patch->getBox();
        const IntVector<DIM> interior_box_lower = interior_box.lower();
        const IntVector<DIM> interior_box_upper = interior_box.upper();

        double* F = F_data->getPointer();
        double* phi = phi_data->getPointer();
        double* control_volume = control_volume_data->getPointer();
        double integral_F_on_patch = 0.0;
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

      } // end loop over patches in level
    } // end loop over levels in hierarchy

  } // end if statement on (region_indicator > 0)

  return tbox::MPI::sumReduction(integral_F);
}


/* computeSurfaceIntegral() */
template <int DIM> 
double LevelSetMethodToolbox<DIM>::computeSurfaceIntegral(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int F_handle,
  const int phi_handle,
  const int grad_phi_handle,
  const int control_volume_handle,
  const int F_component,
  const int phi_component,
  const int delta_width)
{
  double integral_F = 0.0;

  // loop over PatchHierarchy and compute the integral on each Patch
  // by calling Fortran subroutines
  const int num_levels = patch_hierarchy->getNumberLevels();

  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "computeSurfaceIntegral(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get dx and epsilon
      Pointer< CartesianPatchGeometry<DIM> > patch_geom =
        patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
      double max_dx = dx[0];
      for (int k = 1; k < DIM; k++) {
        if (max_dx < dx[k]) max_dx = dx[k];
      }
      double epsilon = delta_width*max_dx;

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,double> > F_data =
        patch->getPatchData( F_handle);
      Pointer< CellData<DIM,double> > phi_data =
        patch->getPatchData( phi_handle );
      Pointer< CellData<DIM,double> > grad_phi_data =
        patch->getPatchData( grad_phi_handle );
      Pointer< CellData<DIM,double> > control_volume_data =
        patch->getPatchData( control_volume_handle );
  
      Box<DIM> F_ghostbox = F_data->getGhostBox();
      const IntVector<DIM> F_ghostbox_lower = F_ghostbox.lower();
      const IntVector<DIM> F_ghostbox_upper = F_ghostbox.upper();

      Box<DIM> phi_ghostbox = phi_data->getGhostBox();
      const IntVector<DIM> phi_ghostbox_lower = phi_ghostbox.lower();
      const IntVector<DIM> phi_ghostbox_upper = phi_ghostbox.upper();

      Box<DIM> grad_phi_ghostbox = grad_phi_data->getGhostBox();
      const IntVector<DIM> grad_phi_ghostbox_lower = 
        grad_phi_ghostbox.lower();
      const IntVector<DIM> grad_phi_ghostbox_upper = 
        grad_phi_ghostbox.upper();

      Box<DIM> control_volume_ghostbox = control_volume_data->getGhostBox();
      const IntVector<DIM> control_volume_ghostbox_lower = 
        control_volume_ghostbox.lower();
      const IntVector<DIM> control_volume_ghostbox_upper = 
        control_volume_ghostbox.upper();

      // interior box
      Box<DIM> interior_box = patch->getBox();
      const IntVector<DIM> interior_box_lower = interior_box.lower();
      const IntVector<DIM> interior_box_upper = interior_box.upper();

      double* F = F_data->getPointer(F_component);
      double* phi = phi_data->getPointer(phi_component);
      double* control_volume = control_volume_data->getPointer();
      double* grad_phi[LSM_DIM_MAX];
      for (int k=0; k<DIM; k++) {
        grad_phi[k] = grad_phi_data->getPointer(k);
      }
      double integral_F_on_patch = 0.0;
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

    } // end loop over patches in level
  } // end loop over levels in hierarchy

  return tbox::MPI::sumReduction(integral_F);
}


/* computeStableAdvectionDt() */
template <int DIM> 
double LevelSetMethodToolbox<DIM>::computeStableAdvectionDt(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int velocity_handle,
  const int control_volume_handle,
  const double cfl_number)
{
  Pointer< CartesianGridGeometry<DIM> > grid_geometry = 
    patch_hierarchy->getGridGeometry();
  const double* dx_level0 = grid_geometry->getDx();
  double max_advection_dt = DBL_MAX;

  // loop over PatchHierarchy and compute the maximum stable 
  // advection dt by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    const IntVector<DIM> ratio_to_coarsest = level->getRatio();
  
    double dx[LSM_DIM_MAX];
    for (int dir = 0; dir < DIM; dir++) {
      dx[dir] = dx_level0[dir]/ratio_to_coarsest[dir];
    }
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "computeStableAdvectionDt(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      double max_advection_dt_on_patch = -1;  // bogus value overwritten
                                              // by Fortran subroutine

      Pointer< CellData<DIM,double> > vel_data = 
        patch->getPatchData( velocity_handle );
      Pointer< CellData<DIM,double> > control_volume_data = 
        patch->getPatchData( control_volume_handle );

      Box<DIM> vel_box = vel_data->getBox();
      const IntVector<DIM> vel_box_lower = vel_box.lower();
      const IntVector<DIM> vel_box_upper = vel_box.upper();

      Box<DIM> vel_ghostbox = vel_data->getGhostBox();
      const IntVector<DIM> vel_ghostbox_lower = vel_ghostbox.lower();
      const IntVector<DIM> vel_ghostbox_upper = vel_ghostbox.upper();

      Box<DIM> control_volume_ghostbox = control_volume_data->getGhostBox();
      const IntVector<DIM> control_volume_ghostbox_lower =
        control_volume_ghostbox.lower();
      const IntVector<DIM> control_volume_ghostbox_upper =
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

  return tbox::MPI::minReduction(max_advection_dt);
}


/* computeStableNormalVelocityDt() */
template <int DIM> 
double LevelSetMethodToolbox<DIM>::computeStableNormalVelocityDt(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int normal_velocity_handle,
  const int grad_phi_plus_handle,
  const int grad_phi_minus_handle,
  const int control_volume_handle,
  const double cfl_number)
{
  Pointer< CartesianGridGeometry<DIM> > grid_geometry = 
    patch_hierarchy->getGridGeometry();
  const double* dx_level0 = grid_geometry->getDx();
  double max_normal_vel_dt = DBL_MAX;

  // loop over PatchHierarchy and compute the maximum stable 
  // advection dt by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    const IntVector<DIM> ratio_to_coarsest = level->getRatio();
  
    double dx[LSM_DIM_MAX];
    for (int dir = 0; dir < DIM; dir++) {
      dx[dir] = dx_level0[dir]/ratio_to_coarsest[dir];
    }
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "computeStableNormalVelocityDt(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      double max_normal_vel_dt_on_patch = -1;  // bogus value overwritten
                                               // by Fortran subroutine

      Pointer< CellData<DIM,double> > vel_data =
        patch->getPatchData(normal_velocity_handle);
      Pointer< CellData<DIM,double> > grad_phi_plus_data =
        patch->getPatchData( grad_phi_plus_handle );
      Pointer< CellData<DIM,double> > grad_phi_minus_data =
        patch->getPatchData( grad_phi_minus_handle );
      Pointer< CellData<DIM,double> > control_volume_data =
        patch->getPatchData( control_volume_handle );

      Box<DIM> vel_box = vel_data->getBox();
      const IntVector<DIM> vel_box_lower = vel_box.lower();
      const IntVector<DIM> vel_box_upper = vel_box.upper();

      Box<DIM> vel_ghostbox = vel_data->getGhostBox();
      const IntVector<DIM> vel_ghostbox_lower = vel_ghostbox.lower();
      const IntVector<DIM> vel_ghostbox_upper = vel_ghostbox.upper();

      Box<DIM> grad_phi_plus_ghostbox = grad_phi_plus_data->getGhostBox();
      const IntVector<DIM> grad_phi_plus_ghostbox_lower =
        grad_phi_plus_ghostbox.lower();
      const IntVector<DIM> grad_phi_plus_ghostbox_upper =
        grad_phi_plus_ghostbox.upper();
      Box<DIM> grad_phi_minus_ghostbox = grad_phi_minus_data->getGhostBox();
      const IntVector<DIM> grad_phi_minus_ghostbox_lower =
        grad_phi_minus_ghostbox.lower();
      const IntVector<DIM> grad_phi_minus_ghostbox_upper =
        grad_phi_minus_ghostbox.upper();

      Box<DIM> control_volume_ghostbox = control_volume_data->getGhostBox();
      const IntVector<DIM> control_volume_ghostbox_lower =
        control_volume_ghostbox.lower();
      const IntVector<DIM> control_volume_ghostbox_upper =
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

  return tbox::MPI::minReduction(max_normal_vel_dt);
}


/* maxNormOfDifference() */
template <int DIM> 
double LevelSetMethodToolbox<DIM>::maxNormOfDifference(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int field1_handle,
  const int field2_handle,
  const int control_volume_handle,
  const int field1_component,
  const int field2_component)
{
  double max_norm_diff = 0;

  // loop over PatchHierarchy and compute the max norm of (field1-field2)
  // by calling Fortran routines
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::" 
                  << "maxNormOfDifference(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,double> > field1_data =
        patch->getPatchData( field1_handle );
      Pointer< CellData<DIM,double> > field2_data =
        patch->getPatchData( field2_handle );
      Pointer< CellData<DIM,double> > control_volume_data =
        patch->getPatchData( control_volume_handle );
  
      Box<DIM> field1_ghostbox = field1_data->getGhostBox();
      const IntVector<DIM> field1_ghostbox_lower = field1_ghostbox.lower();
      const IntVector<DIM> field1_ghostbox_upper = field1_ghostbox.upper();

      Box<DIM> field2_ghostbox = field2_data->getGhostBox();
      const IntVector<DIM> field2_ghostbox_lower = field2_ghostbox.lower();
      const IntVector<DIM> field2_ghostbox_upper = field2_ghostbox.upper();

      Box<DIM> control_volume_ghostbox = 
        control_volume_data->getGhostBox();
      const IntVector<DIM> control_volume_ghostbox_lower = 
        control_volume_ghostbox.lower();
      const IntVector<DIM> control_volume_ghostbox_upper = 
        control_volume_ghostbox.upper();

      // interior box
      Box<DIM> interior_box = field1_data->getBox();
      const IntVector<DIM> interior_box_lower = interior_box.lower();
      const IntVector<DIM> interior_box_upper = interior_box.upper();

      double* field1 = field1_data->getPointer(field1_component);
      double* field2 = field2_data->getPointer(field2_component);
      double* control_volume = control_volume_data->getPointer();
      int control_volume_sgn = 1;

      double max_norm_diff_on_patch = 0.0;

      if ( DIM == 3 ) {
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

      } else if ( DIM == 2 ) {
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

      } else if ( DIM == 1 ) {
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

  return tbox::MPI::maxReduction(max_norm_diff);
}


/* computeControlVolumes() */
template <int DIM>
void LevelSetMethodToolbox<DIM>::computeControlVolumes(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  int control_volume_handle)
{
  int finest_ln = patch_hierarchy->getFinestLevelNumber();
  int ln;
  for ( ln=finest_ln; ln >= 0; --ln ) {

    /*
     * On every level, first assign cell volume to vector weight.
     */

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);

    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) {
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      Pointer< CartesianPatchGeometry<DIM> > patch_geometry =
        patch->getPatchGeometry();
      const double* dx = patch_geometry->getDx();
      double cell_vol = dx[0];
      for (int i = 1; i<DIM; i++) {
        cell_vol *= dx[i];
      } 

      Pointer< CellData<DIM,double> > cv_data =
        patch->getPatchData(control_volume_handle);
      if ( !cv_data ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "computeControlVolumes(): "
                  << "control_volume_handle must refer to a "
                  << "pdat_CellVariable"
                  << endl );
      }
      cv_data->fillAll(cell_vol);
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

      Pointer< PatchLevel<DIM> > next_finer_level
            = patch_hierarchy->getPatchLevel(ln+1);
      BoxArray<DIM> coarsened_boxes = next_finer_level->getBoxes();
      IntVector<DIM> coarsen_ratio = next_finer_level->getRatio();
      coarsen_ratio /= level->getRatio();
      coarsened_boxes.coarsen(coarsen_ratio);

      /*
       * Then set vector weight to 0 wherever there is
       * a nonempty intersection with the next finer level.
       * Note that all assignments are local.
       */

      typename PatchLevel<DIM>::Iterator pi;
      for (pi.initialize(level); pi; pi++) {
        const int pn = *pi;
        Pointer< Patch<DIM> > patch = level->getPatch(pn);
        for ( int i = 0; i < coarsened_boxes.getNumberOfBoxes(); i++ ) {

          Box<DIM> coarse_box = coarsened_boxes(i);
          Box<DIM> intersection = coarse_box*(patch->getBox());
          if ( !intersection.empty() ) {
            Pointer< CellData<DIM, double> > cv_data =
            patch->getPatchData(control_volume_handle);
            cv_data->fillAll(0.0, intersection);

          }  // assignment only in non-empty intersection
        }  // loop over coarsened boxes from finer level
      }  // loop over patches in level
    }  // all levels except finest
  }  // loop over levels
}


/* copySAMRAIData() */
template <int DIM> 
void LevelSetMethodToolbox<DIM>::copySAMRAIData(
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  const int dst_handle,
  const int src_handle,
  const int dst_component,
  const int src_component)
{
  // loop over PatchHierarchy and copy data in Patch interiors
  const int num_levels = patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  "LevelSetMethodToolbox::"
                  << "copySAMRAIData(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,double> > dst_data =
        patch->getPatchData( dst_handle );
      Pointer< CellData<DIM,double> > src_data =
        patch->getPatchData( src_handle );

      Box<DIM> dst_ghostbox = dst_data->getGhostBox();
      const IntVector<DIM> dst_ghostbox_lower = dst_ghostbox.lower();
      const IntVector<DIM> dst_ghostbox_upper = dst_ghostbox.upper();

      Box<DIM> src_ghostbox = src_data->getGhostBox();
      const IntVector<DIM> src_ghostbox_lower = src_ghostbox.lower();
      const IntVector<DIM> src_ghostbox_upper = src_ghostbox.upper();

      Box<DIM> fillbox = dst_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      double* dst = dst_data->getPointer(dst_component);
      double* src = src_data->getPointer(src_component);

      if (DIM == 3) {

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

      } else if (DIM == 2) {

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

      } else if (DIM == 1) {

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
template <int DIM>
void LevelSetMethodToolbox<DIM>::initializeComputeSpatialDerivativesParameters()
{
  // do nothing if PatchData for grad(phi) has already been created
  if (s_D1_one_ghostcell_handle >= 0) {
    return;
  }

  // get pointer to VariableDatabase
  VariableDatabase<DIM> *var_db = VariableDatabase<DIM>::getDatabase();

  // create zero ghostcell width IntVector
  IntVector<DIM> one_ghostcell_width(1);
  IntVector<DIM> two_ghostcells_width(2);
  IntVector<DIM> three_ghostcells_width(3);

  // create Variables and VariableContext 
  Pointer< CellVariable<DIM,double> > D1_variable;
  if (var_db->checkVariableExists("D1")) {
    D1_variable = var_db->getVariable("D1");
  } else {
    D1_variable = new CellVariable<DIM,double>("D1", 1);
  }
  Pointer< CellVariable<DIM,double> > D2_variable;
  if (var_db->checkVariableExists("D2")) {
    D2_variable = var_db->getVariable("D2");
  } else {
    D2_variable = new CellVariable<DIM,double>("D2", 1);
  }
  Pointer< CellVariable<DIM,double> > D3_variable;
  if (var_db->checkVariableExists("D3")) {
    D3_variable = var_db->getVariable("D3");
  } else {
    D3_variable = new CellVariable<DIM,double>("D3", 1);
  }
  Pointer<VariableContext> one_ghostcell_context = 
    var_db->getContext("LSM_TOOLBOX_SCRATCH::ONE_GHOSTCELL");
  Pointer<VariableContext> two_ghostcells_context = 
    var_db->getContext("LSM_TOOLBOX_SCRATCH::TWO_GHOSTCELLS");
  Pointer<VariableContext> three_ghostcells_context = 
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
template <int DIM>
void LevelSetMethodToolbox<DIM>::initializeComputeUnitNormalParameters()
{
  // do nothing if PatchData for grad(phi) has already been created
  if ( (s_compute_normal_grad_phi_handle >= 0) &&
       (s_compute_normal_grad_phi_plus_handle >= 0) &&
       (s_compute_normal_grad_phi_minus_handle >= 0) ) {
    return;
  }

  // get pointer to VariableDatabase
  VariableDatabase<DIM> *var_db = VariableDatabase<DIM>::getDatabase();

  // create zero ghostcell width IntVector
  IntVector<DIM> zero_ghostcell_width(0);

  // create Variables and VariableContext 
  Pointer< CellVariable<DIM,double> > grad_phi_variable;
  if (var_db->checkVariableExists("grad phi")) {
    grad_phi_variable = var_db->getVariable("grad phi");
  } else {
    grad_phi_variable = new CellVariable<DIM,double>("grad phi", DIM);
  }
  Pointer<VariableContext> scratch_context = 
    var_db->getContext("LSM_TOOLBOX_SCRATCH");
  Pointer<VariableContext> grad_phi_plus_context = 
    var_db->getContext("LSM_TOOLBOX_GRAD_PHI_PLUS");
  Pointer<VariableContext> grad_phi_minus_context = 
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
