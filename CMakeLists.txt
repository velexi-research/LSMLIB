# --------------------------------------------------------------
# CMake initialization
# --------------------------------------------------------------

cmake_minimum_required(VERSION 3.0)
project(lsmlib
        LANGUAGES C CXX Fortran
        VERSION 2.0.0)

# --------------------------------------------------------------
# Configuration options
# --------------------------------------------------------------

# Build shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Floating-point precision
option(SINGLE_PRECISION "Build single-precision libraries" OFF)
if (SINGLE_PRECISION)
    message("-- Setting floating-point precision to 'single'")
else (SINGLE_PRECISION)
    message("-- Setting floating-point precision to 'double'")
endif (SINGLE_PRECISION)

# Enable parallel library
option(ENABLE_PARALLEL "Enable parallel library" OFF)

# Boost
set(BOOST "/usr" CACHE FILEPATH "Location of Boost library")

# MPI
set(MPI "/usr" CACHE FILEPATH "Location of MPI library")

# SAMRAI
set(SAMRAI "/usr" CACHE FILEPATH "Location of SAMRAI library")

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Path variables
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
set(BINARY_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)

# Compiler flags
set(CMAKE_C_FLAGS "${CFLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS "${CXXFLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS "${FFLAGS} ${CMAKE_Fortran_FLAGS}")

# Search path for header files
include_directories(${BINARY_INCLUDE_DIR})

# CMake options
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# --------------------------------------------------------------
# Build requirements
# --------------------------------------------------------------

# --- Software with pre-defined CMake package modules

find_package(Threads REQUIRED)
find_package(Git)

# --- include-what-you-use

find_program(IWYU NAMES include-what-you-use iwyu)

if (IWYU)
    message("-- Found include-what-you-use: ${IWYU}")
    set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE ${IWYU} ${CMAKE_CXX_FLAGS})
else (IWYU)
    message("-- Could not find include-what-you-use.  "
            "Analysis of #includes will be skipped.")
endif (IWYU)

# Enable parallel library
if (ENABLE_PARALLEL)

    # --- BOOST

    find_file(BOOST_CONFIG_H config.hpp ${BOOST}/include/boost)
    find_library(LIBBOOST_PYTHON boost_python ${BOOST}/lib)

    # Verify that BOOST header files were found
    if (BOOST_CONFIG_H)
        message("-- Found Boost header files: TRUE (${BOOST_CONFIG_H})")
        set(BOOST_INCLUDE_DIR ${BOOST}/include)
    else (BOOST_CONFIG_H)
        message(FATAL_ERROR "Could NOT find Boost header files.")
    endif (BOOST_CONFIG_H)

    # Verify that BOOST library files were found
    if (LIBBOOST_PYTHON)
        message("-- Found Boost libraries: TRUE")
    else (LIBBOOST_PYTHON)
        message(FATAL_ERROR "Could NOT find Boost libraries.")
    endif (LIBBOOST_PYTHON)

    # --- MPI

    find_file(MPI_H mpi.h
              ${MPI}/include ${MPI}/include/mpi ${MPI}/include/openmpi)
    find_library(LIBMPI mpi ${MPI}/lib)

    # Verify that MPI header files were found
    if (MPI_H)
        message("-- Found MPI.h: TRUE (${MPI_H})")
        get_filename_component(MPI_INCLUDE_DIR ${MPI_H} DIRECTORY)
    else (MPI_H)
        message(FATAL_ERROR "Could NOT find MPI.h.")
    endif (MPI_H)

    # Verify that MPI library files were found
    if (LIBMPI)
        message("-- Found libmpi: TRUE (${LIBMPI})")
    else (LIBMPI)
        message(FATAL_ERROR "Could NOT find libmpi.")
    endif (LIBMPI)

    # --- SAMRAI

    find_file(SAMRAI_CONFIG_H SAMRAI/SAMRAI_config.h ${SAMRAI}/include)
    find_library(LIBSAMRAI_TBOX SAMRAI_tbox ${SAMRAI}/lib)

    # Verify that SAMRAI header files were found
    if (SAMRAI_CONFIG_H)
        message("-- Found SAMRAI header files: TRUE (${SAMRAI_CONFIG_H})")
        set(SAMRAI_INCLUDE_DIR ${SAMRAI}/include)
    else (SAMRAI_CONFIG_H)
        message(FATAL_ERROR "Could NOT find SAMRAI header files.")
    endif (SAMRAI_CONFIG_H)

    # Verify that SAMRAI library files were found
    if (LIBSAMRAI_TBOX)
        message("-- Found SAMRAI libraries: TRUE (${LIBSAMRAI_TBOX})")
    else (LIBSAMRAI_TBOX)
        message(FATAL_ERROR "Could NOT find SAMRAI libraries.")
    endif (LIBSAMRAI_TBOX)

endif (ENABLE_PARALLEL)

# --------------------------------------------------------------
# Process template files
# --------------------------------------------------------------

# Floating-point precision
if (SINGLE_PRECISION)

    set(LSMLIB_REAL "float")
    set(LSMLIB_ZERO_TOL 1.e-5)
    set(LSMLIB_REAL_MAX "FLT_MAX")
    set(LSMLIB_REAL_MIN "FLT_MIN")
    set(LSMLIB_REAL_EPSILON "FLT_EPSILON")

else (SINGLE_PRECISION)

    set(LSMLIB_DOUBLE_PRECISION ON)

    set(LSMLIB_REAL "double")
    set(LSMLIB_ZERO_TOL 1.e-11)
    set(LSMLIB_REAL_MAX "DBL_MAX")
    set(LSMLIB_REAL_MIN "DBL_MIN")
    set(LSMLIB_REAL_EPSILON "DBL_EPSILON")

    # Set Fortran compiler flags
    get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER}
                            NAME)
    message ("-- Fortran compiler: " ${CMAKE_Fortran_COMPILER})

    if (Fortran_COMPILER_NAME MATCHES "gfortran.*")

        # gfortran
        set(CMAKE_Fortran_FLAGS
            "-fdefault-real-8 -fdefault-double-8 ${CMAKE_Fortran_FLAG}")

    elseif (Fortran_COMPILER_NAME MATCHES "g77.*")

        # g77
        set(CMAKE_Fortran_FLAGS
            "-fdefault-real-8 -fdefault-double-8 ${CMAKE_Fortran_FLAG}")

    elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")

        # ifort
        set(CMAKE_Fortran_FLAGS "-r8 ${CMAKE_Fortran_FLAG}")

    else (Fortran_COMPILER_NAME MATCHES "gfortran.*")

        # Unknown fortran compiler
        message ("Trying GNU Fortran flags.")

        set(CMAKE_Fortran_FLAGS
            "-fdefault-real-8 -fdefault-double-8 ${CMAKE_Fortran_FLAG}")

    endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

endif (SINGLE_PRECISION)

# Generate configuration header files
configure_file(config/LSMLIB_config.h.in
               ${BINARY_INCLUDE_DIR}/LSMLIB_config.h)

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# lsm_toolbox
add_subdirectory(${PROJECT_SOURCE_DIR}/src/toolbox)

# lsm_serial
add_subdirectory(${PROJECT_SOURCE_DIR}/src/serial)

# lsm_parallel
if (ENABLE_PARALLEL)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/parallel)
endif (ENABLE_PARALLEL)

# examples
add_subdirectory(examples/serial EXCLUDE_FROM_ALL)

# --------------------------------------------------------------
# Testing
# --------------------------------------------------------------

# Enable CMake testing
enable_testing()

# Emulate 'make check'
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

# Add googletest directory to build
add_subdirectory(${EXT_PROJECTS_DIR}/googletest)

# Add test directory to build
add_subdirectory(${PROJECT_SOURCE_DIR}/test EXCLUDE_FROM_ALL)

# --------------------------------------------------------------
# Installation
# --------------------------------------------------------------

# Header files
#file(GLOB HEADER_FILES ${SRCDIR}/*.h)
#install(FILES ${HEADER_FILES} DESTINATION include)

# Library files
#install(TARGETS lsm_parallel DESTINATION lib)

# Executable files
# TODO

# --- Add `uninstall` target

# Generate cmake_uninstall.cmake file
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
)

# Define `uninstall` target
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND}
            -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)

# --------------------------------------------------------------
# Source code distribution
# --------------------------------------------------------------

if (GIT_FOUND)
    # Add target to build source code distribution
    set(ARCHIVE_NAME ${CMAKE_PROJECT_NAME}-${PROJECT_VERSION})
    set(ARCHIVE_TARBALL ${ARCHIVE_NAME}.tar.gz)
    add_custom_target(dist
        COMMAND git archive --prefix=${ARCHIVE_NAME}/ HEAD |
                gzip > ${CMAKE_BINARY_DIR}/${ARCHIVE_TARBALL}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    # Add target to verify source code distribution
    if (UNIX)
        file(GENERATE OUTPUT check-uninstall.sh
             CONTENT "
remaining=\"`find \${PWD}/${ARCHIVE_NAME}/build/_inst -type f | wc -l`\"
if test \"\${remaining}\" -ne 0; then
    echo \"** Stage directory not empty: \${remaining} file (s) remaining!\"
    exit 1
fi"
        )

        add_custom_target(distcheck
            COMMAND rm -rf ${ARCHIVE_NAME}
            COMMAND gzip -cd ${ARCHIVE_TARBALL} | tar xfz -
            COMMAND cd ${ARCHIVE_NAME} && cmake -H. -Bbuild
            COMMAND cd ${ARCHIVE_NAME}/build && make
            COMMAND cd ${ARCHIVE_NAME}/build && make test
            COMMAND cd ${ARCHIVE_NAME}/build && make DESTDIR=$${PWD}/_inst install
            COMMAND cd ${ARCHIVE_NAME}/build && make DESTDIR=$${PWD}/_inst uninstall
            COMMAND sh check-uninstall.sh
            COMMAND cd ${ARCHIVE_NAME}/build && make clean
            COMMAND rm -rf ${ARCHIVE_NAME}
            COMMAND echo "*** Package ${ARCHIVE_TARBALL} ready for distribution."
        )
        add_dependencies(distcheck dist)

    endif (UNIX)
endif (GIT_FOUND)

# --------------------------------------------------------------
# Build preparations
# --------------------------------------------------------------

# Create build output directories
# TODO: file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
file(MAKE_DIRECTORY ${LIBRARY_OUTPUT_PATH})
