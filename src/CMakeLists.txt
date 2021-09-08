# --------------------------------------------------------------
# Customizations for lsmlib library
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Template parameters
# --------------------------------------------------------------

# Floating-point precision
if (SINGLE_PRECISION)

    set(lsmlib_zero_tol 1.e-5)
    set(tiny_nonzero_number 1.d-35)

else (SINGLE_PRECISION)

    set(lsmlib_zero_tol 1.e-11)
    set(tiny_nonzero_number 1.d-99)

endif (SINGLE_PRECISION)

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions
                    ${CMAKE_CURRENT_SOURCE_DIR}/fast_marching_method
                    ${CMAKE_CURRENT_SOURCE_DIR}/fast_marching_method/core
                    ${CMAKE_CURRENT_SOURCE_DIR}/field_extension
                    ${CMAKE_CURRENT_SOURCE_DIR}/geometry
                    ${CMAKE_CURRENT_SOURCE_DIR}/localization
                    ${CMAKE_CURRENT_SOURCE_DIR}/reinitialization
                    ${CMAKE_CURRENT_SOURCE_DIR}/toolbox
                    ${CMAKE_CURRENT_SOURCE_DIR}/utils)


# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# boundary conditions
add_subdirectory(boundary_conditions)

# fast marching method
add_subdirectory(fast_marching_method)

# field extension
add_subdirectory(field_extension)

# geometry
add_subdirectory(geometry)

# localization
add_subdirectory(localization)

# reinitialization
add_subdirectory(reinitialization)

# toolbox
add_subdirectory(toolbox)

# utilities
add_subdirectory(utils)

add_library(lsm
            ${LSM_BOUNDARY_CONDITIONS_SOURCE_FILES}
            ${LSM_FAST_MARCHING_METHOD_SOURCE_FILES}
            ${LSM_FIELD_EXTENSION_SOURCE_FILES}
            ${LSM_GEOMETRY_SOURCE_FILES}
            ${LSM_LOCALIZATION_SOURCE_FILES}
            ${LSM_REINITIALIZATION_SOURCE_FILES}
            ${LSM_TOOLBOX_SOURCE_FILES}
            ${LSM_UTILS_SOURCE_FILES})

# --------------------------------------------------------------
# Installation
# --------------------------------------------------------------

# Library files
install(TARGETS lsm DESTINATION ${INSTALL_LIB_DIR})

# --------------------------------------------------------------
# Build preparations
# --------------------------------------------------------------

# Create build output directories
file(MAKE_DIRECTORY ${LIBRARY_OUTPUT_PATH})
