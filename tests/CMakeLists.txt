# -----------------------------------------------------------------------------
# Testing customizations
# -----------------------------------------------------------------------------

# Build requirements: No customizations
# Parameters (Install): No customizations
# Targets (Install): No customizations

# -----------------------------------------------------------------------------
# Parameters (Build)
# -----------------------------------------------------------------------------

# Set header and library search paths
include_directories(${BINARY_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS})
link_directories(${GTEST_LIB_DIRS})

# Construct list of test programs
set(TEST_PROGRAMS
    test_calculus_toolbox)

# -----------------------------------------------------------------------------
# Targets (Build)
# -----------------------------------------------------------------------------

# Custom `tests` target to build test programs
add_custom_target(tests DEPENDS ${TEST_PROGRAMS})

# Set up build target for each test program
foreach(test_program ${TEST_PROGRAMS})
    add_executable(
        ${test_program}
        ${test_program}.cc
    )
    add_dependencies(${test_program} lsm googletest)
    target_link_libraries(
        ${test_program}
        lsm
        gtest
        gtest_main
    )
endforeach()

# -----------------------------------------------------------------------------
# Build preparations
# -----------------------------------------------------------------------------

# --- Create build directories

file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

# --- Set up test to run via GoogleTest

include(googletest)
foreach(test_program ${TEST_PROGRAMS})
    gtest_discover_tests(${test_program})
endforeach()
