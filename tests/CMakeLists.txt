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

# GoogleTest
include(googletest)

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

function(add_test_target TEST_PROGRAM SOURCE_FILES)
    add_executable(${TEST_PROGRAM} ${SOURCE_FILES})
    add_dependencies(${TEST_PROGRAM} LSMLIB::lsm googletest)
    target_include_directories(${TEST_PROGRAM} PUBLIC ${GTEST_INCLUDE_DIRS})
    target_link_directories(${TEST_PROGRAM} PUBLIC ${GTEST_LIB_DIRS})
    target_link_libraries(
        ${TEST_PROGRAM}
        PRIVATE lsm
        gtest
        gtest_main
    )
endfunction()

# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------

# --- Sub-directories

add_subdirectory(toolbox)

# --- Build targets

# Custom `tests` target to build test programs
add_custom_target(tests DEPENDS toolbox-tests)
