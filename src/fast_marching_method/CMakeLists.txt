# -----------------------------------------------------------------------------
# Customizations for fast marching method components
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
file(GLOB LSM_FMM_HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
set(LSM_FMM_HEADER_FILES ${LSM_FMM_HEADER_FILES} PARENT_SCOPE)

# Source files
file(GLOB LSM_FMM_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
set(LSM_FMM_SOURCE_FILES ${LSM_FMM_SOURCE_FILES} PARENT_SCOPE)
