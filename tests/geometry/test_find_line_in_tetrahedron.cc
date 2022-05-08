/*
 * Test program for LSM3D_findLineInTetrahedron()
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

/*
 * This program tests that the LSM3D_findLineInTetrahedron() function
 * correctly finds the endpoints of the {phi=0,psi=0} line in the specified
 * tetrahedron.
 *
 * The linear functions for phi and psi were arbitrarily chosen by
 * selecting random alpha and beta vectors defined by:
 *
 *   phi = alpha_0 + alpha_1*x + alpha_2*y + alpha_3*z
 *   psi = beta_0 + beta_1*x + beta_2*y + beta_3*z
 *
 * with
 *
 *   alpha[0] = 0.5; alpha[1] = 1.0; alpha[2] = -0.25; alpha[3] = 1.0
 *
 *   beta[0] = 0.0; beta[1] = -2.0; beta[2] = 1.5; beta[3] = 0.5.
 *
 * For this choice of alpha and beta, the {phi=0,psi=0} line is
 *
 * (-0.75, -1.0, 0.0) + t (-1.625, -2.5, 1.0)
 *
 * The corners of the test tetrahedra were chosen to test the correctness
 * of the endpoints computed by the LSM3D_findLineInTetrahedron() function
 * for all possible cases of intersection of a line with a tetrahedron.
 */

#include <__algorithm/max.h>        // for max
#include <math.h>                   // for fabs
#include <stdio.h>                  // for printf
#include "gtest/gtest-message.h"    // for Message
#include "gtest/gtest-test-part.h"  // for TestPartResult
#include "gtest/gtest_pred_impl.h"  // for CmpHelperGT, SuiteApiResolver
#include "lsm_geometry3d.h"         // for LSM3D_findLineInTetrahedron
#include "lsmlib_config.h"          // for LSMLIB_ZERO_TOL, LSMLIB_REAL

/*
 * Helper function declarations
 */
void check_endpt_on_line(
    LSMLIB_REAL *endpt, LSMLIB_REAL *point_on_line, LSMLIB_REAL *true_line_dir);
void check_line_direction(
    LSMLIB_REAL *endpt1, LSMLIB_REAL *endpt2, LSMLIB_REAL *true_line_dir);

/*
 * Test fixtures
 */
class LSMGeometry3DTest : public ::testing::Test {
  protected:
    // --- Data members

    // Input variables
    LSMLIB_REAL alpha[4], beta[4];
    LSMLIB_REAL x1[3], x2[3], x3[3], x4[3];
    LSMLIB_REAL phi[4], psi[4];
    LSMLIB_REAL point_on_line[3], true_line_dir[3];

    // Output variables
    LSMLIB_REAL endpt1[3], endpt2[3];
};

/*
 * Tests
 */

// {phi=0,psi=0} line goes through single corner of tetrahedron
TEST_F(LSMGeometry3DTest, findLineInTetrahedron1) {
    printf("\nTesting case: {phi=0,psi=0} line goes through single corner\n");
    x1[0] = -0.75;
    x1[1] = -1.0;
    x1[2] = 0.0;
    x2[0] = -0.25;
    x2[1] = -1.0;
    x2[2] = 0.0;
    x3[0] = -0.75;
    x3[1] = -1.5;
    x3[2] = 0.0;
    x4[0] = -0.75;
    x4[1] = -1.0;
    x4[2] = -0.5;

    phi[0] = 0.0;
    phi[1] = 0.5;
    phi[2] = 0.125;
    phi[3] = -0.5;

    psi[0] = 0.0;
    psi[1] = -1.0;
    psi[2] = -0.75;
    psi[3] = -0.25;

    point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
    true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 3);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_NEAR(dist, 0, LSMLIB_ZERO_TOL);
}

// {phi=0,psi=0} line goes through two corners of tetrahedron
TEST_F(LSMGeometry3DTest, findLineInTetrahedron2) {
    printf("\nTesting case: {phi=0,psi=0} line goes through two corners\n");
    x1[0] = -0.75;
    x1[1] = -1.0;
    x1[2] = 0.0;
    x2[0] = -0.25;
    x2[1] = -1.0;
    x2[2] = 0.0;
    x3[0] = 0.0625;
    x3[1] = 0.25;
    x3[2] = -0.5;
    x4[0] = -0.75;
    x4[1] = -1.0;
    x4[2] = -0.5;

    phi[0] = 0.0;
    phi[1] = 0.5;
    phi[2] = 0.0;
    phi[3] = -0.5;

    psi[0] = 0.0;
    psi[1] = -1.0;
    psi[2] = 0.0;
    psi[3] = -0.25;

    point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
    true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 2);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line goes through one corner and non-adjacent edge of tetrahedron
TEST_F(LSMGeometry3DTest, findLineInTetrahedron3) {
    printf("\nTesting case: {phi=0,psi=0} line goes through one corner ");
    printf("and one non-adjacent edge\n");
    x1[0] = -0.75;
    x1[1] = -1.0;
    x1[2] = 0.0;
    x2[0] = 3.5;
    x2[1] = 4.0;
    x2[2] = -2.0;
    x3[0] = 1.5;
    x3[1] = 4.0;
    x3[2] = -2.0;
    x4[0] = -0.75;
    x4[1] = -1.0;
    x4[2] = 0.25;

    phi[0] = 0.0;
    phi[1] = 1.0;
    phi[2] = -1.0;
    phi[3] = 0.25;

    psi[0] = 0.0;
    psi[1] = -2.0;
    psi[2] = 2.0;
    psi[3] = 0.125;

    point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
    true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 3);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line goes through one corner and opposing face of tetrahedron
TEST_F(LSMGeometry3DTest, findLineInTetrahedron4) {
    printf("\nTesting case: {phi=0,psi=0} line goes through a corner ");
    printf("and the opposing face\n");
    x1[0] = -0.75;
    x1[1] = -1.0;
    x1[2] = 0.0;
    x2[0] = -0.25;
    x2[1] = -1.0;
    x2[2] = 0.0;
    x3[0] = -0.75;
    x3[1] = -0.5;
    x3[2] = 0.0;
    x4[0] = -0.75;
    x4[1] = -1.0;
    x4[2] = -0.5;

    phi[0] = 0.0;
    phi[1] = 0.5;
    phi[2] = -0.125;
    phi[3] = -0.5;

    psi[0] = 0.0;
    psi[1] = -1.0;
    psi[2] = 0.75;
    psi[3] = -0.25;

    point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
    true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 4);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line goes through a single point on a single edge
TEST_F(LSMGeometry3DTest, findLineInTetrahedron5) {
    printf("\nTesting case: {phi=0,psi=0} line goes through a single point ");
    printf("on a single edge\n");
    x1[0] = -1.0;
    x1[1] = -1.0;
    x1[2] = 0.0;
    x2[0] = -0.5;
    x2[1] = -1.0;
    x2[2] = 0.0;
    x3[0] = -0.75;
    x3[1] = -0.5;
    x3[2] = 0.0;
    x4[0] = -0.75;
    x4[1] = -1.0;
    x4[2] = 1.0;

    phi[0] = -0.25;
    phi[1] = 0.25;
    phi[2] = -0.125;
    phi[3] = 1.0;

    psi[0] = 0.5;
    psi[1] = -0.5;
    psi[2] = 0.75;
    psi[3] = 0.5;

    point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
    true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 2);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_NEAR(dist, 0, LSMLIB_ZERO_TOL);
}

// {phi=0,psi=0} line goes through two adjacent edges
TEST_F(LSMGeometry3DTest, findLineInTetrahedron6) {
    printf("\nTesting case: {phi=0,psi=0} line goes through two adjacent ");
    printf("edges\n");
    x1[0] = -3.0;
    x1[1] = -6.0;
    x1[2] = 2.0;
    x2[0] = 3.5;
    x2[1] = 4.0;
    x2[2] = -2.0;
    x3[0] = 1.5;
    x3[1] = 4.0;
    x3[2] = -2.0;
    x4[0] = -0.75;
    x4[1] = -1.0;
    x4[2] = -2.0;

    phi[0] = 1.0;
    phi[1] = 1.0;
    phi[2] = -1.0;
    phi[3] = -2.0;

    psi[0] = -2.0;
    psi[1] = -2.0;
    psi[2] = 2.0;
    psi[3] = -1.0;

    point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
    true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 2);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line goes through two non-adjacent edges
TEST_F(LSMGeometry3DTest, findLineInTetrahedron7) {
    printf("\nTesting case: {phi=0,psi=0} line goes through two non-adjacent ");
    printf("edges\n");
    x1[0] = -3.0;
    x1[1] = -6.0;
    x1[2] = 2.0;
    x2[0] = 1.5;
    x2[1] = 4.0;
    x2[2] = -2.0;
    x3[0] = 2.5;
    x3[1] = 3.0;
    x3[2] = -2.0;
    x4[0] = 2.5;
    x4[1] = 5.0;
    x4[2] = -2.0;

    phi[0] = 1.0;
    phi[1] = -1.0;
    phi[2] = 0.25;
    phi[3] = -0.25;

    psi[0] = -2.0;
    psi[1] = 2.0;
    psi[2] = -1.5;
    psi[3] = 1.5;

    point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
    true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 4);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line goes through an edge and a non-adjacent face
TEST_F(LSMGeometry3DTest, findLineInTetrahedron8) {
    printf("\nTesting case: {phi=0,psi=0} line goes through an edge and ");
    printf("a non-adjacent face\n");
    x1[0] = -3.0;
    x1[1] = -6.0;
    x1[2] = 2.0;
    x2[0] = 1.5;
    x2[1] = 4.0;
    x2[2] = -2.0;
    x3[0] = 2.5;
    x3[1] = 2.0;
    x3[2] = -3.0;
    x4[0] = 2.5;
    x4[1] = 5.0;
    x4[2] = -2.0;

    phi[0] = 1.0;
    phi[1] = -1.0;
    phi[2] = -0.5;
    phi[3] = -0.25;

    psi[0] = -2.0;
    psi[1] = 2.0;
    psi[2] = -3.5;
    psi[3] = 1.5;

    point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
    true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 3);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line goes through two faces
TEST_F(LSMGeometry3DTest, findLineInTetrahedron9) {
    printf("\nTesting case: {phi=0,psi=0} line goes through two faces\n");
    x1[0] = -3.0;
    x1[1] = -5.0;
    x1[2] = 2.0;
    x2[0] = 1.5;
    x2[1] = 4.0;
    x2[2] = -2.0;
    x3[0] = 2.5;
    x3[1] = 2.0;
    x3[2] = -3.0;
    x4[0] = 2.5;
    x4[1] = 5.0;
    x4[2] = -2.0;

    phi[0] = 0.75;
    phi[1] = -1.0;
    phi[2] = -0.5;
    phi[3] = -0.25;

    psi[0] = -0.5;
    psi[1] = 2.0;
    psi[2] = -3.5;
    psi[3] = 1.5;

    point_on_line[0] = -0.75; point_on_line[1] = -1.0; point_on_line[2] = 0.0;
    true_line_dir[0] = -1.625; true_line_dir[1] = -2.5; true_line_dir[2] = 1.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 2);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line parallel to x-axis
TEST_F(LSMGeometry3DTest, findLineInTetrahedron10) {
    printf("\nTesting case: {phi=0,psi=0} line parallel to x-axis\n");
    x1[0] = -3.0;
    x1[1] = -6.0;
    x1[2] = 2.0;
    x2[0] = 3.5;
    x2[1] = 4.0;
    x2[2] = -2.0;
    x3[0] = 1.5;
    x3[1] = 4.0;
    x3[2] = -2.0;
    x4[0] = -0.75;
    x4[1] = -1.0;
    x4[2] = -2.0;

    phi[0] = -3.25;
    phi[1] = 2.75;
    phi[2] = 2.75;
    phi[3] = -2.25;

    psi[0] = -7.75;
    psi[1] = 6.25;
    psi[2] = 6.25;
    psi[3] = 1.25;

    point_on_line[0] = 0.0; point_on_line[1] = -0.5; point_on_line[2] = -0.25;
    true_line_dir[0] = 1.0; true_line_dir[1] = 0.0; true_line_dir[2] = 0.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 2);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line parallel to y-axis
TEST_F(LSMGeometry3DTest, findLineInTetrahedron11) {
    printf("\nTesting case: {phi=0,psi=0} line parallel to y-axis\n");
    x1[0] = -3.0;
    x1[1] = -6.0;
    x1[2] = 2.0;
    x2[0] = 3.5;
    x2[1] = 4.0;
    x2[2] = -2.0;
    x3[0] = 1.5;
    x3[1] = 4.0;
    x3[2] = -2.0;
    x4[0] = -0.75;
    x4[1] = -1.0;
    x4[2] = -2.0;

    phi[0] = -0.25;
    phi[1] = 2.25;
    phi[2] = 0.25;
    phi[3] = -2.0;

    psi[0] = -4.75;
    psi[1] = 5.75;
    psi[2] = 3.75;
    psi[3] = 1.5;

    point_on_line[0] = -0.5; point_on_line[1] = 0.0; point_on_line[2] = -0.25;
    true_line_dir[0] = 0.0; true_line_dir[1] = 1.0; true_line_dir[2] = 0.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 2);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line parallel to z-axis
TEST_F(LSMGeometry3DTest, findLineInTetrahedron12) {
    printf("\nTesting case: {phi=0,psi=0} line parallel to z-axis\n");
    x1[0] = -3.0;
    x1[1] = -6.0;
    x1[2] = 2.0;
    x2[0] = 3.5;
    x2[1] = 4.0;
    x2[2] = -2.0;
    x3[0] = 1.5;
    x3[1] = 4.0;
    x3[2] = -2.0;
    x4[0] = -0.75;
    x4[1] = -1.0;
    x4[2] = -2.0;

    phi[0] = -9.75;
    phi[1] = 6.75;
    phi[2] = 4.75;
    phi[3] = -2.5;

    psi[0] = 2.75;
    psi[1] = -0.75;
    psi[2] = -2.75;
    psi[3] = 0.0;

    point_on_line[0] = 0.5; point_on_line[1] = 0.25; point_on_line[2] = 0.0;
    true_line_dir[0] = 0.0; true_line_dir[1] = 0.0; true_line_dir[2] = 1.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 2);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line in xy-plane with z = 0
TEST_F(LSMGeometry3DTest, findLineInTetrahedron13) {
    printf("\nTesting case: {phi=0,psi=0} line in xy-plane with z = 0\n");
    x1[0] = -0.25;
    x1[1] = -0.25;
    x1[2] = -0.5;
    x2[0] = 0.75;
    x2[1] = -0.25;
    x2[2] = -0.5;
    x3[0] = -0.25;
    x3[1] = 0.75;
    x3[2] = -0.5;
    x4[0] = -0.25;
    x4[1] = -0.25;
    x4[2] = 0.5;

    phi[0] = -1.5;
    phi[1] = -0.5;
    phi[2] = -0.5;
    phi[3] =  0.5;

    psi[0] = 0.0;
    psi[1] = 1.0;
    psi[2] = 1.0;
    psi[3] = -1.0;

    point_on_line[0] = 0.0; point_on_line[1] = 0.0; point_on_line[2] = 0.0;
    true_line_dir[0] = -3.0; true_line_dir[1] = 3.0; true_line_dir[2] = 0.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 2);

    // check that end points lie on {phi=0,psi=0} line by verifying that the cross product
    // of (endpt - point_on_line) and true_line_dir is 0
    check_endpt_on_line(endpt1, point_on_line, true_line_dir);
    check_endpt_on_line(endpt2, point_on_line, true_line_dir);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}

// {phi=0,psi=0} line in xy-plane with z = 0
TEST_F(LSMGeometry3DTest, findLineInTetrahedron14) {
    // Preparations
    x1[0] = -0.25;
    x1[1] = -0.25;
    x1[2] = -0.5;
    x2[0] = 0.75;
    x2[1] = -0.25;
    x2[2] = -0.5;
    x3[0] = -0.25;
    x3[1] = 0.75;
    x3[2] = -0.5;
    x4[0] = -0.25;
    x4[1] = -0.25;
    x4[2] = 0.5;

    phi[0] = -0.5;
    phi[1] =  0.5;
    phi[2] = -1.5;
    phi[3] =  0.5;

    psi[0] =  0.5;
    psi[1] =  1.5;
    psi[2] = -0.5;
    psi[3] = -0.5;

    point_on_line[0] = 0.0; point_on_line[1] = 0.0; point_on_line[2] = 0.0;
    true_line_dir[0] = 1.0; true_line_dir[1] = 1.0; true_line_dir[2] = 0.0;

    /*
     * Exercise functionality
     */

    int num_intersections = LSM3D_findLineInTetrahedron(
        endpt1, endpt2,
        x1, x2, x3, x4,
        phi, psi);

    /*
     * Check results
     */

    // check the number of intersections of the {phi=0,psi=0} line with the tetrahedron
    EXPECT_EQ(num_intersections, 3);

    // check that the two endpoints are not the same
    double dist = std::max({fabs(endpt1[0] - endpt2[0]),
                            fabs(endpt1[1] - endpt2[1]),
                            fabs(endpt1[2] - endpt2[2])});
    EXPECT_GT(dist, LSMLIB_ZERO_TOL);

    // check line direction
    check_line_direction(endpt1, endpt2, true_line_dir);
}


/*
 * Helper function definitions
 */
void check_endpt_on_line(
    LSMLIB_REAL *endpt, LSMLIB_REAL *point_on_line, LSMLIB_REAL *true_line_dir)
{
    double cross_product_l1_norm = fabs( (endpt[0] - point_on_line[0])*true_line_dir[1]
                                        -(endpt[1] - point_on_line[1])*true_line_dir[0] )
                                 + fabs( (endpt[2] - point_on_line[2])*true_line_dir[0]
                                        -(endpt[0] - point_on_line[0])*true_line_dir[2] )
                                 + fabs( (endpt[1] - point_on_line[1])*true_line_dir[2]
                                        -(endpt[2] - point_on_line[2])*true_line_dir[1] );

    EXPECT_NEAR(cross_product_l1_norm, 0, LSMLIB_ZERO_TOL);
}

void check_line_direction(
    LSMLIB_REAL *endpt1, LSMLIB_REAL *endpt2, LSMLIB_REAL *true_line_dir)
{
    double line_dir[3];

    line_dir[0] = endpt2[0] - endpt1[0];
    line_dir[1] = endpt2[1] - endpt1[1];
    line_dir[2] = endpt2[2] - endpt1[2];

    // Note: this implementation assumes that the inf-norm of (endpt1 and endpt2) is
    //       greater than LSMLIB_ZERO_TOL.
    if (fabs(line_dir[2]) > LSMLIB_ZERO_TOL) {
      line_dir[0] /= line_dir[2];
      line_dir[1] /= line_dir[2];
      line_dir[2] /= line_dir[2];
    } else if (fabs(line_dir[1]) > LSMLIB_ZERO_TOL) {
      line_dir[0] /= line_dir[1];
      line_dir[1] /= line_dir[1];
      line_dir[2] /= line_dir[1];
    } else {
      line_dir[0] /= line_dir[0];
      line_dir[1] /= line_dir[0];
      line_dir[2] /= line_dir[0];
    }

    // check that direction of the line joining endpt1 and endpt2 is the same as the
    // direction of the {phi=0,psi=0} line by verifying that the cross product
    // (endpt2 - endpt1) and true_line_dir is 0
    double cross_product_l1_norm =
          fabs(line_dir[0]*true_line_dir[1]-line_dir[1]*true_line_dir[0])
        + fabs(line_dir[0]*true_line_dir[2]-line_dir[2]*true_line_dir[0])
        + fabs(line_dir[1]*true_line_dir[2]-line_dir[2]*true_line_dir[1]);

    EXPECT_NEAR(cross_product_l1_norm, 0, LSMLIB_ZERO_TOL);
}
