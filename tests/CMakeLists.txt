# =============================================================================
# LSMLIB Tests
# =============================================================================

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------

# LSMLIB
add_library(LSMLIB::lsm STATIC IMPORTED)
set_property(TARGET LSMLIB::lsm PROPERTY
             IMPORTED_LOCATION "${LIBRARY_OUTPUT_PATH}/liblsm.a"
)

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

function(add_test_target TEST_PROGRAM SOURCE_FILES)
    add_executable(${TEST_PROGRAM} ${SOURCE_FILES})
    add_dependencies(${TEST_PROGRAM} LSMLIB::lsm)
    target_include_directories(${TEST_PROGRAM} PUBLIC ${GTEST_INCLUDE_DIRS})
    target_link_libraries(
        ${TEST_PROGRAM}
        PRIVATE lsm
        ${GTEST_LIBRARIES}
        ${GTEST_MAIN_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
    )
endfunction()

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

# Component tests
add_subdirectory(boundary_conditions)
add_subdirectory(fast_marching_method)
add_subdirectory(geometry)
add_subdirectory(toolbox)

# Custom `tests` target to build test programs
add_custom_target(tests DEPENDS
                  boundary-condition-tests
                  fmm-tests
                  geometry-tests
                  toolbox-tests)
