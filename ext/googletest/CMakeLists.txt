# --------------------------------------------------------------
# CMake initialization
# --------------------------------------------------------------

cmake_minimum_required(VERSION 3.5)
project(googletest C CXX)
include(ExternalProject)

# --------------------------------------------------------------
# Configure googletest build (as an external project)
# --------------------------------------------------------------

ExternalProject_Add(googletest
    URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip

    CMAKE_ARGS -Dgtest_force_shared_crt=ON
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
               -DBUILD_GTEST=ON

    PREFIX "${CMAKE_CURRENT_BINARY_DIR}"

    # Disable install step
    INSTALL_COMMAND ""
)

# --------------------------------------------------------------
# Define googletest variables required for main package build
# --------------------------------------------------------------

# Define GTEST_INCLUDE_DIRS and make it available in parent scope
ExternalProject_Get_Property(googletest source_dir)
set(GTEST_INCLUDE_DIRS ${source_dir}/googletest/include PARENT_SCOPE)
include_directories(${GTEST_INCLUDE_DIRS})

# Define GTEST_LIB_DIRS and make it available in parent scope
ExternalProject_Get_Property(googletest binary_dir)
set(GTEST_LIB_DIRS ${binary_dir}/lib PARENT_SCOPE)
