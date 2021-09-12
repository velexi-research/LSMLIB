/*
 * Test program for delta and heaviside functions
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

#include "gtest/gtest-message.h"    // for Message
#include "gtest/gtest-test-part.h"  // for TestPartResult
#include "gtest/gtest_pred_impl.h"  // for SuiteApiResolver, EXPECT_NEAR

#include "lsmlib_config.h"
#include "lsm_calculus_toolbox.h"

/*
 * Test fixtures
 */
class LSMCalculusToolboxTest : public ::testing::Test {
    protected:
        // Data members
        LSMLIB_REAL x[20];
        LSMLIB_REAL d_phi;
        LSMLIB_REAL eps;


        // Constructor
        LSMCalculusToolboxTest(){
            // Set d_phi and eps (smoothing width)
            d_phi = 0.1;
            eps = 3.0 * d_phi;

            // Set grid values
            for (int i = 0; i < 20; i++) {
                x[i] = -1.0 + i * d_phi;
            }
    }
};

/*
 * Tests
 */
TEST_F(LSMCalculusToolboxTest, Heaviside) {
    LSMLIB_REAL h[20];
    LSMLIB_REAL h_expected[20] = {
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.028834442811219,
        0.195501109477885,
        0.500000000000000,
        0.804498890522115,
        0.971165557188782,
        1.000000000000000,
        1.000000000000000,
        1.000000000000000,
        1.000000000000000,
        1.000000000000000,
        1.000000000000000,
        1.000000000000000};

    // Check values of smoothed Heaviside function
    for (int i = 0; i < 20; i++) {
        EXPECT_NEAR(LSM_HEAVISIDE(x[i], eps), h_expected[i], 1e-6);
    }
}

TEST_F(LSMCalculusToolboxTest, DeltaFunction) {
    LSMLIB_REAL d[20];
    LSMLIB_REAL d_expected[20] = {
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.833333333333334,
        2.500000000000000,
        3.333333333333333,
        2.499999999999999,
        0.833333333333331,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000,
        0.000000000000000};

    // Check values of smoothed delta function
    for (int i = 0; i < 20; i++) {
        EXPECT_NEAR(LSM_DELTA_FUNCTION(x[i], eps), d_expected[i], 1e-6);
    }
}
