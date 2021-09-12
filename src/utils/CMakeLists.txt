# =============================================================================
# utils components
# =============================================================================

# -----------------------------------------------------------------------------
# Target parameters
# -----------------------------------------------------------------------------

# --- Build targets

# Source files
set(LSM_UTILS_SOURCE_FILES)
foreach(FILE IN ITEMS
        lsm_data_arrays.c
        lsm_file.c
        lsm_grid.c
       )
    list(APPEND LSM_UTILS_SOURCE_FILES "utils/${FILE}")
endforeach()
set(LSM_UTILS_SOURCE_FILES ${LSM_UTILS_SOURCE_FILES} PARENT_SCOPE)

# --- Install targets

# Header files
set(LSM_UTILS_HEADER_FILES)
foreach(FILE IN ITEMS
        lsm_data_arrays.h
        lsm_file.h
        lsm_grid.h
        lsm_macros.h
       )
    list(APPEND LSM_UTILS_HEADER_FILES "utils/${FILE}")
endforeach()
set(LSM_UTILS_HEADER_FILES ${LSM_UTILS_HEADER_FILES} PARENT_SCOPE)
