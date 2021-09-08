# --------------------------------------------------------------
# Customizations for reinitialization components
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Process template files
# --------------------------------------------------------------

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

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Source files
file(GLOB LSM_REINITIALIZATION_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c
     ${CMAKE_CURRENT_SOURCE_DIR}/*.f
     ${CMAKE_CURRENT_BINARY_DIR}/*.f)
set(LSM_REINITIALIZATION_SOURCE_FILES
    ${LSM_REINITIALIZATION_SOURCE_FILES} PARENT_SCOPE)

# --------------------------------------------------------------
# Installation
# --------------------------------------------------------------

# Header files
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
install(FILES ${HEADER_FILES}
        DESTINATION ${INSTALL_INCLUDE_DIR}/reinitialization)

# --------------------------------------------------------------
# Build preparations
# --------------------------------------------------------------

# None
