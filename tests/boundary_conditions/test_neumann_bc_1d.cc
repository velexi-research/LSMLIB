/*
 * Unit tests for 1D homogeneous Neumann boundary condition functions.
 *
 * ---------------------------------------------------------------------
 * COPYRIGHT/LICENSE. This file is part of the LSMLIB package. It is
 * subject to the license terms in the LICENSE file found in the
 * top-level directory of this distribution. No part of the LSMLIB
 * package, including this file, may be copied, modified, propagated,
 * or distributed except according to the terms contained in the
 * LICENSE file.
 * ---------------------------------------------------------------------
 */

#include <math.h>                       // for fabs
#include <stddef.h>                     // for NULL
#include "gtest/gtest-message.h"        // for Message
#include "gtest/gtest-test-part.h"      // for TestPartResult
#include "gtest/gtest_pred_impl.h"      // for EXPECT_NEAR, SuiteApiResolver
#include "lsm_boundary_conditions1d.h"  // for LSM1D_HOMOGENEOUS_NEUMANN_ENO1
#include "lsm_spatial_derivatives1d.h"  // for LSM1D_HJ_ENO1, LSM1D_HJ_ENO2
#include "lsmlib_config.h"              // for LSMLIB_REAL

/*
 * Test fixtures
 */
class LSMBoundaryConditions1DTest : public ::testing::Test {
  protected:
    // Data members

    LSMLIB_REAL *phi;
    LSMLIB_REAL *phi_x_plus;
    LSMLIB_REAL *phi_x_minus;
    LSMLIB_REAL *D1 = NULL, *D2 = NULL, *D3 = NULL;
    int ghostcell_width = 3;
    int box_lower[1];
    int box_upper[1];
    int box_dims[1];
    int ghostbox_lower[1];
    int ghostbox_upper[1];
    int ghostbox_dims[1];
    LSMLIB_REAL dx;

    // Constructor
    LSMBoundaryConditions1DTest(){
        // set index space extents
        box_dims[0] = 25;
        box_lower[0] = 0;
        box_upper[0] = box_dims[0]-1;
        ghostbox_lower[0] = box_lower[0] - ghostcell_width;
        ghostbox_upper[0] = box_upper[0] + ghostcell_width;
        ghostbox_dims[0] = ghostbox_upper[0] - ghostbox_lower[0] + 1;

        // set grid spacing
        dx = 1.0/box_dims[0];

        // allocate memory for data
        phi = new LSMLIB_REAL[ghostbox_dims[0]];
        phi_x_plus = new LSMLIB_REAL[ghostbox_dims[0]];
        phi_x_minus = new LSMLIB_REAL[ghostbox_dims[0]];

        // set phi data in interior
        for (int i = 0; i < box_dims[0]; i++) {
            int idx = i+ghostcell_width;
            LSMLIB_REAL x = (i+0.5)*dx;
            phi[idx] = (x-0.25)*(x-0.25);
        }
    }

    // Destructor
    ~LSMBoundaryConditions1DTest(){
        // free heap memory
        delete[] phi;
        delete[] phi_x_plus;
        delete[] phi_x_minus;
        if (D1) {
            delete[] D1;
        }
        if (D2) {
            delete[] D2;
        }
        if (D3) {
            delete[] D3;
        }
    }
};


/*
 * Tests
 */
TEST_F(LSMBoundaryConditions1DTest, ENO1) {
    // allocate memory for divided differences
    D1 = new LSMLIB_REAL[ghostbox_dims[0]];

    // impose homogeneous Neumann boundary conditions
    int bdry_location_idx = 0;
    LSM1D_HOMOGENEOUS_NEUMANN_ENO1(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &box_lower[0],
        &box_upper[0],
        &bdry_location_idx);

    bdry_location_idx = 1;
    LSM1D_HOMOGENEOUS_NEUMANN_ENO1(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &box_lower[0],
        &box_upper[0],
        &bdry_location_idx);

    // compute spatial derivatives using ENO1
    LSM1D_HJ_ENO1(
        phi_x_plus,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        phi_x_minus,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        D1,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &box_lower[0],
        &box_upper[0],
        &dx);

    /*
     * check that derivatives normal to faces of computational domain are 0.
     */

    // x-lower
    LSMLIB_REAL err_x_lower = fabs(phi_x_minus[box_lower[0]+ghostcell_width]);
    EXPECT_NEAR(err_x_lower, 0, 1e-6);

    // x-upper
    LSMLIB_REAL err_x_upper = fabs(phi_x_plus[box_upper[0]+ghostcell_width]);
    EXPECT_NEAR(err_x_upper, 0, 1e-6);
}

TEST_F(LSMBoundaryConditions1DTest, ENO2) {
    // allocate memory for divided differences
    D1 = new LSMLIB_REAL[ghostbox_dims[0]];
    D2 = new LSMLIB_REAL[ghostbox_dims[0]];

    // impose homogeneous Neumann boundary conditions
    int bdry_location_idx = 0;
    LSM1D_HOMOGENEOUS_NEUMANN_ENO2(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &box_lower[0],
        &box_upper[0],
        &bdry_location_idx);

    bdry_location_idx = 1;
    LSM1D_HOMOGENEOUS_NEUMANN_ENO2(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &box_lower[0],
        &box_upper[0],
        &bdry_location_idx);

    // compute spatial derivatives using ENO2
    LSM1D_HJ_ENO2(
      phi_x_plus,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      phi_x_minus,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      phi,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      D1,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      D2,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      &box_lower[0],
      &box_upper[0],
      &dx);

    /*
     * check that derivatives normal to faces of computational domain are 0.
     */

    // x-lower
    LSMLIB_REAL err_x_lower = fabs(phi_x_minus[box_lower[0]+ghostcell_width]);
    EXPECT_NEAR(err_x_lower, 0, 1e-6);

    // x-upper
    LSMLIB_REAL err_x_upper = fabs(phi_x_plus[box_upper[0]+ghostcell_width]);
    EXPECT_NEAR(err_x_upper, 0, 1e-6);
}

TEST_F(LSMBoundaryConditions1DTest, ENO3) {
    // allocate memory for divided differences
    D1 = new LSMLIB_REAL[ghostbox_dims[0]];
    D2 = new LSMLIB_REAL[ghostbox_dims[0]];
    D3 = new LSMLIB_REAL[ghostbox_dims[0]];

    // impose homogeneous Neumann boundary conditions
    int bdry_location_idx = 0;
    LSM1D_HOMOGENEOUS_NEUMANN_ENO3(
      phi,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      &box_lower[0],
      &box_upper[0],
      &bdry_location_idx);

    bdry_location_idx = 1;
    LSM1D_HOMOGENEOUS_NEUMANN_ENO3(
      phi,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      &box_lower[0],
      &box_upper[0],
      &bdry_location_idx);

    // compute spatial derivatives using ENO3
    LSM1D_HJ_ENO3(
      phi_x_plus,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      phi_x_minus,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      phi,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      D1,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      D2,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      D3,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      &box_lower[0],
      &box_upper[0],
      &dx);

    /*
     * check that derivatives normal to faces of computational domain are 0.
     */

    // x-lower
    LSMLIB_REAL err_x_lower = fabs(phi_x_minus[box_lower[0]+ghostcell_width]);
    EXPECT_NEAR(err_x_lower, 0, 1e-6);

    // x-upper
    LSMLIB_REAL err_x_upper = fabs(phi_x_plus[box_upper[0]+ghostcell_width]);
    EXPECT_NEAR(err_x_upper, 0, 1e-6);
}

TEST_F(LSMBoundaryConditions1DTest, WENO5) {
    // allocate memory for divided differences
    D1 = new LSMLIB_REAL[ghostbox_dims[0]];

    // impose homogeneous Neumann boundary conditions
    int bdry_location_idx = 0;
    LSM1D_HOMOGENEOUS_NEUMANN_WENO5(
      phi,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      &box_lower[0],
      &box_upper[0],
      &bdry_location_idx);

    bdry_location_idx = 1;
    LSM1D_HOMOGENEOUS_NEUMANN_WENO5(
      phi,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      &box_lower[0],
      &box_upper[0],
      &bdry_location_idx);

    // compute spatial derivatives using ENO3
    LSM1D_HJ_WENO5(
      phi_x_plus,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      phi_x_minus,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      phi,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      D1,
      &ghostbox_lower[0],
      &ghostbox_upper[0],
      &box_lower[0],
      &box_upper[0],
      &dx);

    /*
     * check that derivatives normal to faces of computational domain are 0.
     */

    // x-lower
    LSMLIB_REAL err_x_lower = fabs(phi_x_minus[box_lower[0]+ghostcell_width]);
    EXPECT_NEAR(err_x_lower, 0, 1e-6);

    // x-upper
    LSMLIB_REAL err_x_upper = fabs(phi_x_plus[box_upper[0]+ghostcell_width]);
    EXPECT_NEAR(err_x_upper, 0, 1e-6);
}
