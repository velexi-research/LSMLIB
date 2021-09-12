# =============================================================================
# Main CMakeLists.txt
# =============================================================================

# -----------------------------------------------------------------------------
# CMake initialization
# -----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.10)
project(LSMLIB
        LANGUAGES C CXX Fortran
        VERSION 2.0.0)

# -----------------------------------------------------------------------------
# CMake configuration
# -----------------------------------------------------------------------------

# CMake build type
set(CMAKE_BUILD_TYPE Release CACHE STRING "CMake build type")

# Compiler requirements
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Make cache variables for install destinations
include(GNUInstallDirs)

# -----------------------------------------------------------------------------
# External modules
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# Build options
# -----------------------------------------------------------------------------

# Build shared libraries
# TODO: fix shared library build
#option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Floating-point precision
option(USE_SINGLE_PRECISION "Build single-precision libraries" OFF)
if (USE_SINGLE_PRECISION)
    message("-- Setting floating-point precision to 'single'")
else (USE_SINGLE_PRECISION)
    message("-- Setting floating-point precision to 'double'")
endif (USE_SINGLE_PRECISION)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

# --- Build parameters

# Construct name to use for package-specific directories (e.g., header files)
execute_process(
    COMMAND "echo" ${CMAKE_PROJECT_NAME}
    COMMAND "tr" '[:upper:]' '[:lower:]'
    OUTPUT_VARIABLE get_project_name_lower_out
    RESULT_VARIABLE get_project_name_lower_result_var
    )
if (NOT "${get_project_name_lower_result_var}" STREQUAL 0)
    message(FATAL_ERROR "Error constructing directory name for package")
endif (NOT "${get_project_name_lower_result_var}" STREQUAL 0)
string(STRIP ${get_project_name_lower_out} PKG_DIR_NAME)

# Set path variables
set(BINARY_INCLUDEDIR include)
set(LIBRARY_OUTPUT_PATH lib)
set(EXECUTABLE_OUTPUT_PATH bin)
set(EXT_PROJECTS_DIR ext)

# Compiler flags
set(CMAKE_Fortran_FLAGS "${FFLAGS} ${CMAKE_Fortran_FLAGS}")
if (NOT USE_SINGLE_PRECISION)

    # Set Fortran compiler flags
    get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
    message("-- Fortran compiler: " ${CMAKE_Fortran_COMPILER})

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
        message("Trying GNU Fortran flags.")

        set(CMAKE_Fortran_FLAGS
            "-fdefault-real-8 -fdefault-double-8 ${CMAKE_Fortran_FLAG}")

    endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

endif (NOT USE_SINGLE_PRECISION)

# --- Install parameters

# Set directory (relative) for installed header files
set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}/${PKG_DIR_NAME})

# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------

# --- Build targets

# Libraries
add_subdirectory(src)

# Examples
add_subdirectory(examples EXCLUDE_FROM_ALL)

# --- Install targets

# Configuration header file
install(FILES ${CMAKE_BINARY_DIR}/${BINARY_INCLUDEDIR}/lsmlib_config.h
        DESTINATION ${INSTALL_INCLUDEDIR})

# --- `uninstall` target

# Generate cmake_uninstall.cmake file
configure_file(
    "${CMAKE_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
)

# Define `uninstall` target
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/cmake_uninstall.cmake
)

# -----------------------------------------------------------------------------
# Testing
# -----------------------------------------------------------------------------

# Enable CMake testing
enable_testing()

# GoogleTest
add_subdirectory(${EXT_PROJECTS_DIR}/googletest EXCLUDE_FROM_ALL)

# Unit tests
add_subdirectory(tests EXCLUDE_FROM_ALL)

# -----------------------------------------------------------------------------
# Packaging
# -----------------------------------------------------------------------------

# Create package version file
include(CMakePackageConfigHelpers)
set_property(TARGET lsm PROPERTY VERSION ${version})
set_property(TARGET lsm PROPERTY
    INTERFACE_LSMLIB_MAJOR_VERSION ${LSMLIB_VERSION_MAJOR})
set_property(TARGET lsm APPEND PROPERTY
    COMPATIBLE_INTERFACE_STRING LSMLIB_MAJOR_VERSION
)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/cmake/LSMLIBConfigVersion.cmake"
    VERSION "${version}"
    COMPATIBILITY AnyNewerVersion
)

# Create source code distribution
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

# -----------------------------------------------------------------------------
# Build preparations
# -----------------------------------------------------------------------------

# --- Process template files

# Floating-point precision
if (USE_SINGLE_PRECISION)

    set(LSMLIB_REAL "float")
    set(LSMLIB_ZERO_TOL 1.e-5)
    set(LSMLIB_REAL_MAX "FLT_MAX")
    set(LSMLIB_REAL_MIN "FLT_MIN")
    set(LSMLIB_REAL_EPSILON "FLT_EPSILON")

else (USE_SINGLE_PRECISION)

    set(LSMLIB_DOUBLE_PRECISION ON)

    set(LSMLIB_REAL "double")
    set(LSMLIB_ZERO_TOL 1.e-11)
    set(LSMLIB_REAL_MAX "DBL_MAX")
    set(LSMLIB_REAL_MIN "DBL_MIN")
    set(LSMLIB_REAL_EPSILON "DBL_EPSILON")

endif (USE_SINGLE_PRECISION)

# Generate configuration header files
configure_file(include/lsmlib_config.h.in include/lsmlib_config.h)
