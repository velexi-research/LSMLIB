# =============================================================================
# LSMLIB fast marching method tests
# =============================================================================

# -----------------------------------------------------------------------------
# Test
# -----------------------------------------------------------------------------

# --- Targets

# Add custom target for tests
set(TEST_PROGRAMS
    test_FMM_Heap
    )
add_custom_target(fmm-tests DEPENDS ${TEST_PROGRAMS})

# Add build target for each test program
foreach(TEST_PROGRAM ${TEST_PROGRAMS})
    add_test_target(${TEST_PROGRAM} ${TEST_PROGRAM}.cc)
endforeach()

# --- GoogleTest configuration

# Set up tests to run via GoogleTest
foreach(TEST_PROGRAM ${TEST_PROGRAMS})
    gtest_discover_tests(${TEST_PROGRAM})
endforeach()
