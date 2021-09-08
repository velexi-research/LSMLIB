# -----------------------------------------------------------------------------
# Customizations for field extension components
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
file(GLOB LSM_FIELD_EXTENSION_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c
     ${CMAKE_CURRENT_SOURCE_DIR}/*.f
     ${CMAKE_CURRENT_BINARY_DIR}/*.f)
set(LSM_FIELD_EXTENSION_SOURCE_FILES
    ${LSM_FIELD_EXTENSION_SOURCE_FILES} PARENT_SCOPE)

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
