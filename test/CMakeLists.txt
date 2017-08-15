# --------------------------------------------------------------
# Testing customizations
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Include search path
include_directories(${GTEST_INCLUDE_DIR})

# --------------------------------------------------------------
# lsm_parallel tests
# --------------------------------------------------------------

if (ENABLE_PARALLEL)

    # --------------------------------------------------------------
    # Build targets
    # --------------------------------------------------------------

    set(LSM_PARALLEL_TESTS ${CMAKE_PROJECT_NAME}-tests)

    # Source files containing tests
    file(GLOB TEST_SOURCE_FILES
         ${PROJECT_SOURCE_DIR}/test/parallel/*.cc)

    add_executable(${LSM_PARALLEL_TESTS} EXCLUDE_FROM_ALL ${TEST_SOURCE_FILES})
    add_dependencies(${LSM_PARALLEL_TESTS} lsm_parallel googletest)
    target_link_libraries(${LSM_PARALLEL_TESTS}
        lsm_parallel
        ${GTEST_LIBS_DIR}/libgtest.a
        ${GTEST_LIBS_DIR}/libgtest_main.a
        ${CMAKE_THREAD_LIBS_INIT})

    # --------------------------------------------------------------
    # Test setup
    # --------------------------------------------------------------
    add_test(${LSM_PARALLEL_TESTS}
             ${CMAKE_CURRENT_BINARY_DIR}/${LSM_PARALLEL_TESTS})
    add_dependencies(check ${LSM_PARALLEL_TESTS})

endif (ENABLE_PARALLEL)

