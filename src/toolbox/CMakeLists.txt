# --------------------------------------------------------------
# Customizations for LSMLIB Toolbox
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Path variables
set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/toolbox)

# Search path for header files
include_directories(${SOURCE_DIR})

# --------------------------------------------------------------
# Process template files
# --------------------------------------------------------------

# --------------------------------------------------------------
# Process template files
# --------------------------------------------------------------

# Floating-point precision
if (SINGLE_PRECISION)

    set(lsmlib_zero_tol 1.e-5)
    set(tiny_nonzero_number 1.d-35)

else (SINGLE_PRECISION)

    set(lsmlib_zero_tol 1.e-11)
    set(tiny_nonzero_number 1.d-99)

endif (SINGLE_PRECISION)

# Generate fortran files
configure_file(lsm_calculus_toolbox2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_calculus_toolbox2d.f)
configure_file(lsm_calculus_toolbox2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_calculus_toolbox2d_local.f)
configure_file(lsm_calculus_toolbox3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_calculus_toolbox3d.f)
configure_file(lsm_curvature2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature2d.f)
configure_file(lsm_curvature2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature2d_local.f)
configure_file(lsm_curvature3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature3d.f)
configure_file(lsm_curvature3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature3d_local.f)
configure_file(lsm_geometry1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry1d.f)
configure_file(lsm_geometry2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry2d.f)
configure_file(lsm_geometry2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry2d_local.f)
configure_file(lsm_geometry3d_fort.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry3d_fort.f)
configure_file(lsm_level_set_evolution1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_level_set_evolution1d.f)
configure_file(lsm_level_set_evolution2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_level_set_evolution2d.f)
configure_file(lsm_level_set_evolution2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_level_set_evolution2d_local.f)
configure_file(lsm_level_set_evolution3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_level_set_evolution3d.f)
configure_file(lsm_level_set_evolution3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_level_set_evolution3d_local.f)
configure_file(lsm_reinitialization1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_reinitialization1d.f)
configure_file(lsm_reinitialization2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_reinitialization2d.f)
configure_file(lsm_reinitialization2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_reinitialization2d_local.f)
configure_file(lsm_reinitialization3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_reinitialization3d.f)
configure_file(lsm_reinitialization3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_reinitialization3d_local.f)
configure_file(lsm_spatial_derivatives1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_spatial_derivatives1d.f)
configure_file(lsm_spatial_derivatives2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_spatial_derivatives2d.f)
configure_file(lsm_spatial_derivatives2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_spatial_derivatives2d_local.f)
configure_file(lsm_spatial_derivatives3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_spatial_derivatives3d.f)
configure_file(lsm_spatial_derivatives3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_spatial_derivatives3d_local.f)
configure_file(lsm_utilities1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_utilities1d.f)
configure_file(lsm_utilities2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_utilities2d.f)
configure_file(lsm_utilities2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_utilities2d_local.f)
configure_file(lsm_utilities3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_utilities3d.f)
configure_file(lsm_utilities3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_utilities3d_local.f)

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Source files
file(GLOB SOURCE_FILES
     ${SOURCE_DIR}/*.c
     ${SOURCE_DIR}/*.f
     ${CMAKE_CURRENT_BINARY_DIR}/*.f)
add_library(lsm_toolbox ${SOURCE_FILES})

# --------------------------------------------------------------
# Installation
# --------------------------------------------------------------

# Header files
file(GLOB HEADER_FILES ${SRCDIR}/*.h)
install(FILES ${HEADER_FILES} DESTINATION include/lsmlib)

# Library files
install(TARGETS lsm_toolbox DESTINATION lib)

# --------------------------------------------------------------
# Build preparations
# --------------------------------------------------------------

# Create build output directories
file(MAKE_DIRECTORY ${LIBRARY_OUTPUT_PATH})
