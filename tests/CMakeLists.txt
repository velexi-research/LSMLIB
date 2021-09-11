# -----------------------------------------------------------------------------
# Testing customizations
# -----------------------------------------------------------------------------

# Build requirements: No customizations
# Parameters (Install): No customizations
# Targets (Install): No customizations

# -----------------------------------------------------------------------------
# Parameters (Build)
# -----------------------------------------------------------------------------

# Construct list of test programs
set(TEST_PROGRAMS
    #test_FMM_Heap
    test_calculus_toolbox)

# -----------------------------------------------------------------------------
# Targets (Build)
# -----------------------------------------------------------------------------

# Custom `tests` target to build test programs
add_custom_target(tests DEPENDS ${TEST_PROGRAMS})

# Set up build target for each test program
foreach(PROGRAM ${TEST_PROGRAMS})
    add_executable(
        ${PROGRAM}
        ${PROGRAM}.cc)
    add_dependencies(${PROGRAM} make-bin-dir lsm googletest)
    target_include_directories(${PROGRAM}
                               PUBLIC ${BINARY_INCLUDE_DIR}
                               PUBLIC ${GTEST_INCLUDE_DIRS})
    target_link_directories(${PROGRAM}
                            PUBLIC ${GTEST_LIB_DIRS})
    target_link_libraries(
        ${PROGRAM}
        lsm
        gtest
        gtest_main
    )
endforeach()

# -----------------------------------------------------------------------------
# Build preparations
# -----------------------------------------------------------------------------

# --- Set up test to run via GoogleTest

include(googletest)
foreach(PROGRAM ${TEST_PROGRAMS})
    gtest_discover_tests(${PROGRAM})
endforeach()
