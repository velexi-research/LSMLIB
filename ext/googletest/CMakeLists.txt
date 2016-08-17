# --------------------------------------------------------------
# CMake initialization
# --------------------------------------------------------------

cmake_minimum_required(VERSION 3.0)
project(googletest_builder C CXX)
include(ExternalProject)

# --------------------------------------------------------------
# Configure googletest build (as an external project)
# --------------------------------------------------------------

ExternalProject_Add(googletest
    URL ${PROJECT_SOURCE_DIR}/googletest-release-1.7.0.tar.gz
    CMAKE_ARGS -Dgtest_force_shared_crt=ON
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
               # TODO: the following may need to be added for next release of
               #       googletest.
               -DBUILD_GTEST=ON

    PREFIX "${CMAKE_CURRENT_BINARY_DIR}"

    # Disable install step
    INSTALL_COMMAND ""
)

# --------------------------------------------------------------
# Define googletest variables required for main package build
# --------------------------------------------------------------

# Define GTEST_INCLUDE_DIR and make it available in parent scope
ExternalProject_Get_Property(googletest source_dir)
set(GTEST_INCLUDE_DIR ${source_dir}/include PARENT_SCOPE)

# Define GTEST_LIB_DIR and make it available in parent scope
ExternalProject_Get_Property(googletest binary_dir)
set(GTEST_LIBS_DIR ${binary_dir} PARENT_SCOPE)

# TODO: The above lines will need to be updated with next release of googletest.
#       The following lines should be closer to the new directory structure for
#       googletest.
#ExternalProject_Get_Property(googletest source_dir)
#set(GTEST_INCLUDE_DIR ${source_dir}/googletest/include PARENT_SCOPE)
#
#ExternalProject_Get_Property(googletest binary_dir)
#set(GTEST_LIBS_DIR ${binary_dir}/googlemock/gtest PARENT_SCOPE)
