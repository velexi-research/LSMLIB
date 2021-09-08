# -----------------------------------------------------------------------------
# Customizations for geometry components
# -----------------------------------------------------------------------------

# CMake initialization: No customizations
# Build - Requirements: No customizations
# Build - Parameters: No customizations
# Install - Parameters: No customizations

# -----------------------------------------------------------------------------
# Build - Targets
# -----------------------------------------------------------------------------

# Source files
file(GLOB LSM_GEOMETRY_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c
     ${CMAKE_CURRENT_SOURCE_DIR}/*.f
     ${CMAKE_CURRENT_BINARY_DIR}/*.f)
set(LSM_GEOMETRY_SOURCE_FILES
    ${LSM_GEOMETRY_SOURCE_FILES} PARENT_SCOPE)

# -----------------------------------------------------------------------------
# Build - Preparations
# -----------------------------------------------------------------------------

# --- Process template files

# Generate fortran files
configure_file(lsm_geometry1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry1d.f)
configure_file(lsm_geometry2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry2d.f)
configure_file(lsm_geometry2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry2d_local.f)
configure_file(lsm_geometry3d_fort.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry3d_fort.f)

configure_file(lsm_curvature2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature2d.f)
configure_file(lsm_curvature2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature2d_local.f)
configure_file(lsm_curvature3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature3d.f)
configure_file(lsm_curvature3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature3d_local.f)

# --- Copy header files to binary folder

file(COPY ${HEADER_FILES}
     DESTINATION ${BINARY_INCLUDE_DIR}/geometry)

# -----------------------------------------------------------------------------
# Install - Targets
# -----------------------------------------------------------------------------

# Header files
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
install(FILES ${HEADER_FILES}
        DESTINATION ${INSTALL_INCLUDE_DIR}/geometry)
