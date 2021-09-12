# =============================================================================
# boundary conditions components
# =============================================================================

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

# --- Build parameters

# Source files
set(LSM_BOUNDARY_CONDITIONS_SOURCE_FILES)
foreach(FILE IN ITEMS
        lsm_boundary_conditions.c
        lsm_boundary_conditions1d.f
        lsm_boundary_conditions2d.f
        lsm_boundary_conditions3d.f
       )
    list(APPEND LSM_BOUNDARY_CONDITIONS_SOURCE_FILES
         "boundary_conditions/${FILE}")
endforeach()
set(LSM_BOUNDARY_CONDITIONS_SOURCE_FILES
    ${LSM_BOUNDARY_CONDITIONS_SOURCE_FILES} PARENT_SCOPE)

# --- Install parameters

# Header files
set(LSM_BOUNDARY_CONDITIONS_HEADER_FILES)
foreach(FILE IN ITEMS
        lsm_boundary_conditions.h
        lsm_boundary_conditions1d.h
        lsm_boundary_conditions2d.h
        lsm_boundary_conditions3d.h
       )
    list(APPEND LSM_BOUNDARY_CONDITIONS_HEADER_FILES
         "boundary_conditions/${FILE}")
endforeach()
set(LSM_BOUNDARY_CONDITIONS_HEADER_FILES
    ${LSM_BOUNDARY_CONDITIONS_HEADER_FILES} PARENT_SCOPE)
