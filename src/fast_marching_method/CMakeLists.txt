# -----------------------------------------------------------------------------
# Customizations for fast marching method components
# -----------------------------------------------------------------------------

# Build requirements: No customizations
# Parameters (Install): No customizations

# -----------------------------------------------------------------------------
# Parameters (Build)
# -----------------------------------------------------------------------------

file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

# -----------------------------------------------------------------------------
# Targets (Build)
# -----------------------------------------------------------------------------

# Source files
file(GLOB LSM_FAST_MARCHING_METHOD_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
set(LSM_FAST_MARCHING_METHOD_SOURCE_FILES
    ${LSM_FAST_MARCHING_METHOD_SOURCE_FILES} PARENT_SCOPE)

# -----------------------------------------------------------------------------
# Targets (Install)
# -----------------------------------------------------------------------------

# Header files
install(FILES ${HEADER_FILES} DESTINATION ${INSTALL_HEADER_PATH})

# -----------------------------------------------------------------------------
# Build preparations
# -----------------------------------------------------------------------------

# Copy header files to binary folder
file(COPY ${HEADER_FILES} DESTINATION ${HEADER_OUTPUT_PATH})
