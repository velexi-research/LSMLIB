# =============================================================================
# LSMLIB toolbox Tests
# =============================================================================

# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------

# --- Build targets

# Add custom target for tests
set(TEST_PROGRAMS
    test_calculus_toolbox)
add_custom_target(toolbox-tests DEPENDS ${TEST_PROGRAMS})

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
