# =============================================================================
# fast marching method components
# =============================================================================

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

# --- Build parameters

# Source files
set(LSM_FMM_SOURCE_FILES)
foreach(FILE IN ITEMS
        FMM_Core.c
        FMM_Heap.c
        lsm_FMM_eikonal2d.c
        lsm_FMM_eikonal3d.c
        lsm_FMM_field_extension2d.c
        lsm_FMM_field_extension3d.c
       )
    list(APPEND LSM_FMM_SOURCE_FILES "fast_marching_method/${FILE}")
endforeach()
set(LSM_FMM_SOURCE_FILES ${LSM_FMM_SOURCE_FILES} PARENT_SCOPE)

# --- Install parameters

# Header files
set(LSM_FMM_HEADER_FILES)
foreach(FILE IN ITEMS
        FMM_Callback_API.h
        FMM_Core.h
        FMM_Heap.h
        FMM_Macros.h
        lsm_FMM_eikonal.h
        lsm_FMM_field_extension.h
        lsm_fast_marching_method.h
       )
    list(APPEND LSM_FMM_HEADER_FILES "fast_marching_method/${FILE}")
endforeach()
set(LSM_FMM_HEADER_FILES ${LSM_FMM_HEADER_FILES} PARENT_SCOPE)
