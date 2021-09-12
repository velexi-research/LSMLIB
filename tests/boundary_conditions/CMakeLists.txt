# =============================================================================
# LSMLIB boundary conditions tests
# =============================================================================

# -----------------------------------------------------------------------------
# Build targets
# -----------------------------------------------------------------------------

# Add custom target for tests
set(TEST_PROGRAMS
    test_neumann_bc_ENO1_1d
    #    test_neumann_bc_ENO1_2d
    #test_neumann_bc_ENO1_3d
    #test_neumann_bc_ENO2_1d
    #test_neumann_bc_ENO2_2d
    #test_neumann_bc_ENO2_3d
    #test_neumann_bc_ENO3_1d
    #test_neumann_bc_ENO3_2d
    #test_neumann_bc_ENO3_3d
    #test_neumann_bc_WENO5_1d
    #test_neumann_bc_WENO5_2d
    #test_neumann_bc_WENO5_3d
)
add_custom_target(boundary-condition-tests DEPENDS ${TEST_PROGRAMS})

# Add build target for each test program
foreach(TEST_PROGRAM ${TEST_PROGRAMS})
    add_test_target(${TEST_PROGRAM} ${TEST_PROGRAM}.cc)
endforeach()

# -----------------------------------------------------------------------------
# Build preparations
# -----------------------------------------------------------------------------

# --- Set up test to run via GoogleTest

foreach(TEST_PROGRAM ${TEST_PROGRAMS})
    gtest_discover_tests(${TEST_PROGRAM})
endforeach()
