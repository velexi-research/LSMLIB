# =============================================================================
# LSMLIB::lsm library
# =============================================================================

# -----------------------------------------------------------------------------
# Template parameters
# -----------------------------------------------------------------------------

# Floating-point precision
if (USE_SINGLE_PRECISION)

    set(lsmlib_zero_tol 1.e-5)
    set(tiny_nonzero_number 1.d-36)

else (USE_SINGLE_PRECISION)

    set(lsmlib_zero_tol 1.e-11)
    set(tiny_nonzero_number 1.d-36)

endif (USE_SINGLE_PRECISION)

# -----------------------------------------------------------------------------
# Targets
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

# --- Build targets

add_library(lsm
            ${LSM_BOUNDARY_CONDITIONS_SOURCE_FILES}
            ${LSM_FMM_SOURCE_FILES}
            ${LSM_FIELD_EXTENSION_SOURCE_FILES}
            ${LSM_GEOMETRY_SOURCE_FILES}
            ${LSM_REINITIALIZATION_SOURCE_FILES}
            ${LSM_TOOLBOX_SOURCE_FILES}
            ${LSM_UTILS_SOURCE_FILES}
)
target_include_directories(lsm PUBLIC
    "$<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>"
    "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/toolbox>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/utils>"
)

# --- Install targets

# Library files
install(TARGETS lsm
        EXPORT LSMLIB
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${INSTALL_INCLUDEDIR}
)

# Header files
install(FILES
        ${LSM_TOOLBOX_HEADER_FILES}
        ${LSM_UTILS_HEADER_FILES}
        DESTINATION ${INSTALL_INCLUDEDIR}
)

# CMake files
install(EXPORT LSMLIB
        FILE LSMLIB.cmake
        NAMESPACE LSMLIB::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}
)

# -----------------------------------------------------------------------------
# Miscellaneous
# -----------------------------------------------------------------------------

# Create CMake script to support exporting from the build tree
export(EXPORT LSMLIB
       FILE "${CMAKE_BINARY_DIR}/cmake/LSMLIB.cmake"
       NAMESPACE LSMLIB::
)
