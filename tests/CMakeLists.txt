# =============================================================================
# LSMLIB Tests
# =============================================================================

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------

# LSMLIB
add_library(LSMLIB::lsm STATIC IMPORTED)
set_target_properties(LSMLIB::lsm PROPERTIES
  IMPORTED_LOCATION
    "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/lib/liblsm.a>"
)

# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------

# --- Build targets

# Construct list of test programs
set(TEST_PROGRAMS
    #test_FMM_Heap
    test_calculus_toolbox)

# Custom `tests` target to build test programs
add_custom_target(tests DEPENDS ${TEST_PROGRAMS})

# Set up build target for each test program
foreach(PROGRAM ${TEST_PROGRAMS})
    add_executable(
        ${PROGRAM}
        ${PROGRAM}.cc)
    add_dependencies(${PROGRAM} LSMLIB::lsm googletest)
    target_include_directories(${PROGRAM} PUBLIC ${GTEST_INCLUDE_DIRS})
    target_link_directories(${PROGRAM} PUBLIC ${GTEST_LIB_DIRS})
    target_link_libraries(
        ${PROGRAM}
        PRIVATE lsm
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
