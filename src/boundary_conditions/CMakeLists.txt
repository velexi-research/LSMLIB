# --------------------------------------------------------------
# Customizations for boundary conditions components
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Source files
file(GLOB LSM_BOUNDARY_CONDITIONS_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c
     ${CMAKE_CURRENT_SOURCE_DIR}/*.f
     ${CMAKE_CURRENT_BINARY_DIR}/*.f)
set(LSM_BOUNDARY_CONDITIONS_SOURCE_FILES
    ${LSM_BOUNDARY_CONDITIONS_SOURCE_FILES} PARENT_SCOPE)

# --------------------------------------------------------------
# Installation
# --------------------------------------------------------------

# Header files
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
install(FILES ${HEADER_FILES}
        DESTINATION ${INSTALL_INCLUDE_DIR}/boundary_conditions)

# --------------------------------------------------------------
# Build preparations
# --------------------------------------------------------------

# None
