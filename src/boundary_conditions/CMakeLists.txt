# -----------------------------------------------------------------------------
# Customizations for boundary conditions components
# -----------------------------------------------------------------------------

# Build requirements: No customizations
# Parameters (Install): No customizations
# Targets (Build): No customizations
# Targets (Install): No customizations
# Build preparations: No customization

# -----------------------------------------------------------------------------
# Parameters (Build)
# -----------------------------------------------------------------------------

# Header files
file(GLOB LSM_BOUNDARY_CONDITIONS_HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
set(LSM_BOUNDARY_CONDITIONS_HEADER_FILES ${LSM_BOUNDARY_CONDITIONS_HEADER_FILES}
    PARENT_SCOPE)

# Source files
file(GLOB LSM_BOUNDARY_CONDITIONS_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c
     ${CMAKE_CURRENT_SOURCE_DIR}/*.f
     ${CMAKE_CURRENT_BINARY_DIR}/*.f)
set(LSM_BOUNDARY_CONDITIONS_SOURCE_FILES
    ${LSM_BOUNDARY_CONDITIONS_SOURCE_FILES} PARENT_SCOPE)
