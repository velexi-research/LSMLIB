/*
 * Unit tests for 3D homogeneous Neumann boundary condition functions.
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

#include <gtest/gtest-message.h>        // for Message
#include <gtest/gtest-test-part.h>      // for TestPartResult
#include <gtest/gtest_pred_impl.h>      // for EXPECT_NEAR, SuiteApiResolver

#include "lsm_boundary_conditions3d.h"  // for LSM3D_HOMOGENEOUS_NEUMANN_ENO1
#include "lsm_spatial_derivatives3d.h"  // for LSM3D_HJ_ENO1, LSM3D_HJ_ENO2
#include "lsmlib_config.h"              // for LSMLIB_REAL

/*
 * Test fixtures
 */
class LSMBoundaryConditions3DTest : public ::testing::Test {
  protected:
    // Data members

    LSMLIB_REAL *phi;
    LSMLIB_REAL *phi_x_plus, *phi_y_plus, *phi_z_plus;
    LSMLIB_REAL *phi_x_minus, *phi_y_minus, *phi_z_minus;
    LSMLIB_REAL *D1 = NULL, *D2 = NULL, *D3 = NULL;
    int ghostcell_width = 3;
    int box_lower[3];
    int box_upper[3];
    int box_dims[3];
    int ghostbox_lower[3];
    int ghostbox_upper[3];
    int ghostbox_dims[3];
    LSMLIB_REAL dx, dy, dz;

    // Constructor
    LSMBoundaryConditions3DTest(){
        // set index space extents
        box_dims[0] = 25;
        box_dims[1] = 20;
        box_dims[2] = 30;
        box_lower[0] = 0;
        box_lower[1] = 0;
        box_lower[2] = 0;
        box_upper[0] = box_dims[0]-1;
        box_upper[1] = box_dims[1]-1;
        box_upper[2] = box_dims[2]-1;
        ghostbox_lower[0] = box_lower[0] - ghostcell_width;
        ghostbox_lower[1] = box_lower[1] - ghostcell_width;
        ghostbox_lower[2] = box_lower[2] - ghostcell_width;
        ghostbox_upper[0] = box_upper[0] + ghostcell_width;
        ghostbox_upper[1] = box_upper[1] + ghostcell_width;
        ghostbox_upper[2] = box_upper[2] + ghostcell_width;
        ghostbox_dims[0] = ghostbox_upper[0] - ghostbox_lower[0] + 1;
        ghostbox_dims[1] = ghostbox_upper[1] - ghostbox_lower[1] + 1;
        ghostbox_dims[2] = ghostbox_upper[2] - ghostbox_lower[2] + 1;

        // set grid spacing
        dx = 1.0/box_dims[0];
        dy = 1.0/box_dims[1];
        dz = 1.0/box_dims[2];

        // allocate memory for data
        phi = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                              ghostbox_dims[2]];
        phi_x_plus = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                                     ghostbox_dims[2]];
        phi_y_plus = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                                     ghostbox_dims[2]];
        phi_z_plus = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                                     ghostbox_dims[2]];
        phi_x_minus = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                                      ghostbox_dims[2]];
        phi_y_minus = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                                      ghostbox_dims[2]];
        phi_z_minus = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                                      ghostbox_dims[2]];

        // set phi data in interior
        for (int k = 0; k < box_dims[2]; k++) {
            for (int j = 0; j < box_dims[1]; j++) {
                for (int i = 0; i < box_dims[0]; i++) {
                    int idx = (i+ghostcell_width)
                            + (j+ghostcell_width)*ghostbox_dims[0]
                            + (k+ghostcell_width)*ghostbox_dims[0]*
                                                  ghostbox_dims[1];

                    LSMLIB_REAL x = (i+0.5)*dx;
                    LSMLIB_REAL y = (j+0.5)*dy;
                    LSMLIB_REAL z = (j+0.5)*dz;

                    phi[idx] = (x-0.25)*(x-0.25) + 4*(y-0.35)*(y-0.35)
                             + 0.25*(z-0.45)*(z-0.45) - 0.25;
                }
            }
        }
    }

    // Destructor
    ~LSMBoundaryConditions3DTest(){
        // free heap memory
        delete[] phi;
        delete[] phi_x_plus;
        delete[] phi_y_plus;
        delete[] phi_z_plus;
        delete[] phi_x_minus;
        delete[] phi_y_minus;
        delete[] phi_z_minus;
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
TEST_F(LSMBoundaryConditions3DTest, ENO1) {
    // variable declarations
    int i, j, k;
    int bdry_location_idx;

    // allocate memory for divided differences
    D1 = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                         ghostbox_dims[2]];

    // impose homogeneous Neumann boundary conditions
    bdry_location_idx = 0;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 1;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 2;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 3;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 4;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 5;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO1(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    // compute spatial derivatives using ENO1
    LSM3D_HJ_ENO1(
        phi_x_plus, phi_y_plus, phi_z_plus,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        phi_x_minus, phi_y_minus, phi_z_minus,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        D1,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &dx, &dy, &dz);

    /*
     * check that derivatives normal to faces of computational domain are 0.
     */

    // x-lower
    i = box_lower[0];
    for (k = 0; k < box_dims[2]; k++) {
        for (j = 0; j < box_dims[1]; j++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_x_lower = fabs(phi_x_minus[idx]);
            EXPECT_NEAR(err_x_lower, 0, 1e-6);
       }
   }

    // x-upper
    i = box_upper[0];
    for (k = 0; k < box_dims[2]; k++) {
        for (j = 0; j < box_dims[1]; j++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_x_upper = fabs(phi_x_plus[idx]);
            EXPECT_NEAR(err_x_upper, 0, 1e-6);
        }
    }

    // y-lower
    j = box_lower[1];
    for (k = 0; k < box_dims[2]; k++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_y_lower = fabs(phi_y_minus[idx]);
            EXPECT_NEAR(err_y_lower, 0, 1e-6);
        }
    }

    // y-upper
    j = box_upper[1];
    for (k = 0; k < box_dims[2]; k++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_y_upper = fabs(phi_y_plus[idx]);
            EXPECT_NEAR(err_y_upper, 0, 1e-6);
        }
    }

    // z-lower
    k = box_lower[2];
    for (j = 0; j < box_dims[1]; j++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_z_lower = fabs(phi_z_minus[idx]);
            EXPECT_NEAR(err_z_lower, 0, 1e-6);
       }
   }

    // z-upper
    k = box_upper[2];
    for (j = 0; j < box_dims[1]; j++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_z_upper = fabs(phi_z_plus[idx]);
            EXPECT_NEAR(err_z_upper, 0, 1e-6);
        }
    }
}

TEST_F(LSMBoundaryConditions3DTest, ENO2) {
    // variable declarations
    int i, j, k;
    int bdry_location_idx;

    // allocate memory for divided differences
    D1 = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                         ghostbox_dims[2]];
    D2 = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                         ghostbox_dims[2]];

    // impose homogeneous Neumann boundary conditions
    bdry_location_idx = 0;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO2(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 1;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO2(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 2;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO2(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 3;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO2(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 4;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO2(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 5;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO2(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    // compute spatial derivatives using ENO2
    LSM3D_HJ_ENO2(
        phi_x_plus, phi_y_plus, phi_z_plus,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        phi_x_minus, phi_y_minus, phi_z_minus,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        D1,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        D2,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &dx, &dy, &dz);

    /*
     * check that derivatives normal to faces of computational domain are 0.
     */

    // x-lower
    i = box_lower[0];
    for (k = 0; k < box_dims[2]; k++) {
        for (j = 0; j < box_dims[1]; j++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_x_lower = fabs(phi_x_minus[idx]);
            EXPECT_NEAR(err_x_lower, 0, 1e-6);
       }
   }

    // x-upper
    i = box_upper[0];
    for (k = 0; k < box_dims[2]; k++) {
        for (j = 0; j < box_dims[1]; j++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_x_upper = fabs(phi_x_plus[idx]);
            EXPECT_NEAR(err_x_upper, 0, 1e-6);
        }
    }

    // y-lower
    j = box_lower[1];
    for (k = 0; k < box_dims[2]; k++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_y_lower = fabs(phi_y_minus[idx]);
            EXPECT_NEAR(err_y_lower, 0, 1e-6);
        }
    }

    // y-upper
    j = box_upper[1];
    for (k = 0; k < box_dims[2]; k++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_y_upper = fabs(phi_y_plus[idx]);
            EXPECT_NEAR(err_y_upper, 0, 1e-6);
        }
    }

    // z-lower
    k = box_lower[2];
    for (j = 0; j < box_dims[1]; j++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_z_upper = fabs(phi_z_minus[idx]);
            EXPECT_NEAR(err_z_upper, 0, 1e-6);
       }
   }

    // z-upper
    k = box_upper[2];
    for (j = 0; j < box_dims[1]; j++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_z_upper = fabs(phi_z_plus[idx]);
            EXPECT_NEAR(err_z_upper, 0, 1e-6);
        }
    }
}

TEST_F(LSMBoundaryConditions3DTest, ENO3) {
    // variable declarations
    int i, j, k;
    int bdry_location_idx;

    // allocate memory for divided differences
    D1 = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                         ghostbox_dims[2]];
    D2 = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                         ghostbox_dims[2]];
    D3 = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                         ghostbox_dims[2]];

    // impose homogeneous Neumann boundary conditions
    bdry_location_idx = 0;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO3(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 1;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO3(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 2;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO3(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 3;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO3(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 4;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO3(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 5;
    LSM3D_HOMOGENEOUS_NEUMANN_ENO3(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    // compute spatial derivatives using ENO3
    LSM3D_HJ_ENO3(
        phi_x_plus, phi_y_plus, phi_z_plus,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        phi_x_minus, phi_y_minus, phi_z_minus,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        D1,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        D2,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        D3,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &dx, &dy, &dz);

    /*
     * check that derivatives normal to faces of computational domain are 0.
     */

    // x-lower
    i = box_lower[0];
    for (k = 0; k < box_dims[2]; k++) {
       for (j = 0; j < box_dims[1]; j++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_x_lower = fabs(phi_x_minus[idx]);
            EXPECT_NEAR(err_x_lower, 0, 1e-6);
      }
   }

    // x-upper
    i = box_upper[0];
    for (k = 0; k < box_dims[2]; k++) {
        for (j = 0; j < box_dims[1]; j++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_x_upper = fabs(phi_x_plus[idx]);
            EXPECT_NEAR(err_x_upper, 0, 1e-6);
        }
    }

    // y-lower
    j = box_lower[1];
    for (k = 0; k < box_dims[2]; k++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_y_lower = fabs(phi_y_minus[idx]);
            EXPECT_NEAR(err_y_lower, 0, 1e-6);
        }
    }

    // y-upper
    j = box_upper[1];
    for (k = 0; k < box_dims[2]; k++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_y_lower = fabs(phi_y_plus[idx]);
            EXPECT_NEAR(err_y_lower, 0, 1e-6);
        }
    }

    // z-lower
    k = box_lower[2];
    for (j = 0; j < box_dims[1]; j++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_z_lower = fabs(phi_z_minus[idx]);
            EXPECT_NEAR(err_z_lower, 0, 1e-6);
       }
   }

    // z-upper
    k = box_upper[2];
    for (j = 0; j < box_dims[1]; j++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_z_upper = fabs(phi_z_plus[idx]);
            EXPECT_NEAR(err_z_upper, 0, 1e-6);
        }
    }
}

TEST_F(LSMBoundaryConditions3DTest, WENO5) {
    // variable declarations
    int i, j, k;
    int bdry_location_idx;

    // allocate memory for divided differences
    D1 = new LSMLIB_REAL[ghostbox_dims[0] * ghostbox_dims[1] *
                         ghostbox_dims[2]];

    // impose homogeneous Neumann boundary conditions
    bdry_location_idx = 0;
    LSM3D_HOMOGENEOUS_NEUMANN_WENO5(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 1;
    LSM3D_HOMOGENEOUS_NEUMANN_WENO5(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 2;
    LSM3D_HOMOGENEOUS_NEUMANN_WENO5(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 3;
    LSM3D_HOMOGENEOUS_NEUMANN_WENO5(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 4;
    LSM3D_HOMOGENEOUS_NEUMANN_WENO5(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    bdry_location_idx = 5;
    LSM3D_HOMOGENEOUS_NEUMANN_WENO5(
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &bdry_location_idx);

    // compute spatial derivatives using WENO5
    LSM3D_HJ_WENO5(
        phi_x_plus, phi_y_plus, phi_z_plus,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        phi_x_minus, phi_y_minus, phi_z_minus,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        phi,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        D1,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &ghostbox_lower[2],
        &ghostbox_upper[2],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        &box_lower[2],
        &box_upper[2],
        &dx, &dy, &dz);

    /*
     * check that derivatives normal to faces of computational domain are 0.
     */

    // x-lower
    i = box_lower[0];
    for (k = 0; k < box_dims[2]; k++) {
        for (j = 0; j < box_dims[1]; j++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_x_lower = fabs(phi_x_minus[idx]);
            EXPECT_NEAR(err_x_lower, 0, 1e-6);
       }
   }

    // x-upper
    i = box_upper[0];
    for (k = 0; k < box_dims[2]; k++) {
        for (j = 0; j < box_dims[1]; j++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_x_upper = fabs(phi_x_plus[idx]);
            EXPECT_NEAR(err_x_upper, 0, 1e-6);
        }
    }

    // y-lower
    j = box_lower[1];
    for (k = 0; k < box_dims[2]; k++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_y_lower = fabs(phi_y_minus[idx]);
            EXPECT_NEAR(err_y_lower, 0, 1e-6);
        }
    }

    // y-upper
    j = box_upper[1];
    for (k = 0; k < box_dims[2]; k++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_y_upper = fabs(phi_y_plus[idx]);
            EXPECT_NEAR(err_y_upper, 0, 1e-6);
        }
    }

    // z-lower
    k = box_lower[2];
    for (j = 0; j < box_dims[1]; j++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_z_lower = fabs(phi_z_minus[idx]);
            EXPECT_NEAR(err_z_lower, 0, 1e-6);
       }
   }

    // z-upper
    k = box_upper[2];
    for (j = 0; j < box_dims[1]; j++) {
        for (i = 0; i < box_dims[0]; i++) {
            int idx = (i+ghostcell_width)
                    + (j+ghostcell_width)*ghostbox_dims[0]
                    + (k+ghostcell_width)*ghostbox_dims[0]*ghostbox_dims[1];

            LSMLIB_REAL err_z_upper = fabs(phi_z_plus[idx]);
            EXPECT_NEAR(err_z_upper, 0, 1e-6);
        }
    }
}
