# --------------------------------------------------------------
# Customizations for utilities
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Source files
file(GLOB LSM_UTILS_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
set(LSM_UTILS_SOURCE_FILES
    ${LSM_UTILS_SOURCE_FILES} PARENT_SCOPE)

# --------------------------------------------------------------
# Installation
# --------------------------------------------------------------

# Header files
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
install(FILES ${HEADER_FILES}
        DESTINATION ${INSTALL_INCLUDE_DIR}/utils)

# --------------------------------------------------------------
# Build preparations
# --------------------------------------------------------------

# None
