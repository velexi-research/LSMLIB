# -----------------------------------------------------------------------------
# Customizations for reinitialization components
# -----------------------------------------------------------------------------

# CMake initialization: No customizations
# Build - Requirements: No customizations
# Build - Parameters: No customizations
# Install - Parameters: No customizations

# -----------------------------------------------------------------------------
# Build - Targets
# -----------------------------------------------------------------------------

# Source files
file(GLOB LSM_REINITIALIZATION_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c
     ${CMAKE_CURRENT_SOURCE_DIR}/*.f
     ${CMAKE_CURRENT_BINARY_DIR}/*.f)
set(LSM_REINITIALIZATION_SOURCE_FILES
    ${LSM_REINITIALIZATION_SOURCE_FILES} PARENT_SCOPE)

# -----------------------------------------------------------------------------
# Build - Preparations
# -----------------------------------------------------------------------------

# --- Process template files

# Generate fortran files
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

# --- Copy header files to binary folder

file(COPY ${HEADER_FILES}
     DESTINATION ${BINARY_INCLUDE_DIR}/reinitialization)

# -----------------------------------------------------------------------------
# Install - Targets
# -----------------------------------------------------------------------------

# Header files
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
install(FILES ${HEADER_FILES}
        DESTINATION ${INSTALL_INCLUDE_DIR}/reinitialization)
