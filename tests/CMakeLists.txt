# --------------------------------------------------------------
# Testing customizations
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Set header and library search paths
include_directories(${GTEST_INCLUDE_DIRS})
link_directories(${GTEST_LIB_DIRS})

# Construct list of test programs
set(TEST_PROGRAMS
    hello_test
    hello_test_2)

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Custom `tests` target to build test programs
add_custom_target(tests DEPENDS ${TEST_PROGRAMS})

# Set up build target for each test program
foreach(test_program ${TEST_PROGRAMS})
    add_executable(
        ${test_program} EXCLUDE_FROM_ALL
        ${test_program}.cc
    )
    add_dependencies(${test_program} lsm googletest)
    target_link_libraries(
        ${test_program}
        gtest
        gtest_main
    )
endforeach()

# --------------------------------------------------------------
# Set up test to run via GoogleTest
# --------------------------------------------------------------

include(googletest)
foreach(test_program ${TEST_PROGRAMS})
    gtest_discover_tests(${test_program})
endforeach()
