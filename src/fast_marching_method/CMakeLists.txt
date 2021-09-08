# -----------------------------------------------------------------------------
# Customizations for fast marching method components
# -----------------------------------------------------------------------------

# CMake initialization: No customizations
# Build - Requirements: No customizations
# Build - Parameters: No customizations
# Install - Parameters: No customizations

# -----------------------------------------------------------------------------
# Build - Targets
# -----------------------------------------------------------------------------

# Source files
file(GLOB LSM_FAST_MARCHING_METHOD_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c
     ${CMAKE_CURRENT_SOURCE_DIR}/core/*.c)
set(LSM_FAST_MARCHING_METHOD_SOURCE_FILES
    ${LSM_FAST_MARCHING_METHOD_SOURCE_FILES} PARENT_SCOPE)

# -----------------------------------------------------------------------------
# Build - Preparations
# -----------------------------------------------------------------------------

# Copy header files to binary folder
file(COPY ${HEADER_FILES}
     DESTINATION ${BINARY_INCLUDE_DIR}/fast_marching_method)

# -----------------------------------------------------------------------------
# Install - Targets
# -----------------------------------------------------------------------------

# Header files
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
install(FILES ${HEADER_FILES}
        DESTINATION ${INSTALL_INCLUDE_DIR}/fast_marching_method)
