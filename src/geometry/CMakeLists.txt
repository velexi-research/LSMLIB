# =============================================================================
# geometry components
# =============================================================================

# -----------------------------------------------------------------------------
# Target parameters
# -----------------------------------------------------------------------------

# --- Build targets

# Source files
set(LSM_GEOMETRY_SOURCE_FILES)
foreach(FILE IN ITEMS
        lsm_geometry3d_c.c
       )
    list(APPEND LSM_GEOMETRY_SOURCE_FILES "geometry/${FILE}")
endforeach()

foreach(FILE IN ITEMS
        lsm_curvature2d.f
        lsm_curvature2d_local.f
        lsm_curvature3d.f
        lsm_curvature3d_local.f
        lsm_geometry1d.f
        lsm_geometry2d.f
        lsm_geometry2d_local.f
        lsm_geometry3d_fort.f
       )
        list(APPEND LSM_GEOMETRY_SOURCE_FILES
             "${CMAKE_CURRENT_BINARY_DIR}/${FILE}")
endforeach()
set(LSM_GEOMETRY_SOURCE_FILES ${LSM_GEOMETRY_SOURCE_FILES} PARENT_SCOPE)

# --- Install targets

# Header files
set(LSM_GEOMETRY_HEADER_FILES)
foreach(FILE IN ITEMS
        lsm_curvature2d.h
        lsm_curvature2d_local.h
        lsm_curvature3d.h
        lsm_curvature3d_local.h
        lsm_geometry1d.h
        lsm_geometry2d.h
        lsm_geometry2d_local.h
        lsm_geometry3d.h
       )
    list(APPEND LSM_GEOMETRY_HEADER_FILES "geometry/${FILE}")
endforeach()
set(LSM_GEOMETRY_HEADER_FILES ${LSM_GEOMETRY_HEADER_FILES} PARENT_SCOPE)

# -----------------------------------------------------------------------------
# Build preparations
# -----------------------------------------------------------------------------

# --- Process template files
#
#     Note: requires `lsmlib_zero_tol` and `tiny_nonzero_number` variables
#           to be set.

# Generate fortran files
configure_file(lsm_geometry1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry1d.f)
configure_file(lsm_geometry2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry2d.f)
configure_file(lsm_geometry2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry2d_local.f)
configure_file(lsm_geometry3d_fort.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_geometry3d_fort.f)

configure_file(lsm_curvature2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature2d.f)
configure_file(lsm_curvature2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature2d_local.f)
configure_file(lsm_curvature3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature3d.f)
configure_file(lsm_curvature3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_curvature3d_local.f)
