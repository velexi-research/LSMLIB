# -----------------------------------------------------------------------------
# Customizations for lsm library
# -----------------------------------------------------------------------------

# Build requirements: No customizations
# Parameters (Install): No customizations
# Build preparations: No customizations

# -----------------------------------------------------------------------------
# Parameters (Build)
# -----------------------------------------------------------------------------

# --- Template parameters

# Floating-point precision
if (SINGLE_PRECISION)

    set(lsmlib_zero_tol 1.e-5)
    set(tiny_nonzero_number 1.d-36)

else (SINGLE_PRECISION)

    set(lsmlib_zero_tol 1.e-11)
    set(tiny_nonzero_number 1.d-36)

endif (SINGLE_PRECISION)

# -----------------------------------------------------------------------------
# Targets (Build)
# -----------------------------------------------------------------------------

# --- Subdirectories

# boundary conditions
add_subdirectory(boundary_conditions)

# fast marching method
add_subdirectory(fast_marching_method)

# field extension
add_subdirectory(field_extension)

# geometry
add_subdirectory(geometry)

# reinitialization
add_subdirectory(reinitialization)

# toolbox
add_subdirectory(toolbox)

# utilities
add_subdirectory(utils)

# --- Header files and libraries

set(HEADER_FILES
    ${LSM_BOUNDARY_CONDITIONS_HEADER_FILES}
    ${LSM_FMM_HEADER_FILES}
    ${LSM_FIELD_EXTENSION_HEADER_FILES}
    ${LSM_GEOMETRY_HEADER_FILES}
    ${LSM_REINITIALIZATION_HEADER_FILES}
    ${LSM_TOOLBOX_HEADER_FILES}
    ${LSM_UTILS_HEADER_FILES})
add_custom_target(lsm-header-files
    COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HEADER_FILES} ${HEADER_OUTPUT_PATH})
add_dependencies(lsm-header-files make-header-dir)

add_library(lsm
            ${LSM_BOUNDARY_CONDITIONS_SOURCE_FILES}
            ${LSM_FMM_SOURCE_FILES}
            ${LSM_FIELD_EXTENSION_SOURCE_FILES}
            ${LSM_GEOMETRY_SOURCE_FILES}
            ${LSM_REINITIALIZATION_SOURCE_FILES}
            ${LSM_TOOLBOX_SOURCE_FILES}
            ${LSM_UTILS_SOURCE_FILES})
add_dependencies(lsm make-lib-dir lsm-header-files)

# -----------------------------------------------------------------------------
# Targets (Install)
# -----------------------------------------------------------------------------

# Header files
install(DIRECTORY ${HEADER_OUTPUT_PATH}/ DESTINATION ${INSTALL_HEADER_PATH})

# Library files
install(TARGETS lsm DESTINATION ${INSTALL_LIB_DIR})
