# =============================================================================
# reinitialization components
# =============================================================================

# -----------------------------------------------------------------------------
# Target parameters
# -----------------------------------------------------------------------------

# --- Build targets

# Source files
set(LSM_REINITIALIZATION_SOURCE_FILES)
foreach(FILE IN ITEMS
        lsm_reinitialization1d.f
        lsm_reinitialization2d.f
        lsm_reinitialization2d_local.f
        lsm_reinitialization3d.f
        lsm_reinitialization3d_local.f
       )
        list(APPEND LSM_REINITIALIZATION_SOURCE_FILES
             "${CMAKE_CURRENT_BINARY_DIR}/${FILE}")
endforeach()
set(LSM_REINITIALIZATION_SOURCE_FILES ${LSM_REINITIALIZATION_SOURCE_FILES}
    PARENT_SCOPE)

# --- Install targets

# Header files
set(LSM_REINITIALIZATION_HEADER_FILES)
foreach(FILE IN ITEMS
        lsm_reinitialization1d.h
        lsm_reinitialization2d.h
        lsm_reinitialization2d_local.h
        lsm_reinitialization3d.h
        lsm_reinitialization3d_local.h
       )
    list(APPEND LSM_REINITIALIZATION_HEADER_FILES "toolbox/${FILE}")
endforeach()
set(LSM_REINITIALIZATION_HEADER_FILES ${LSM_REINITIALIZATION_HEADER_FILES}
    PARENT_SCOPE)

# -----------------------------------------------------------------------------
# Build preparations
# -----------------------------------------------------------------------------

# --- Process template files
#
#     Note: requires `lsmlib_zero_tol` and `tiny_nonzero_number` variables
#           to be set.

# Generate fortran files
configure_file(lsm_reinitialization1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_reinitialization1d.f)
configure_file(lsm_reinitialization2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_reinitialization2d.f)
configure_file(lsm_reinitialization2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_reinitialization2d_local.f)
configure_file(lsm_reinitialization3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_reinitialization3d.f)
configure_file(lsm_reinitialization3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_reinitialization3d_local.f)
