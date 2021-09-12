# =============================================================================
# field extension components
# =============================================================================

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

# --- Build parameters

# Source files
set(LSM_FIELD_EXTENSION_SOURCE_FILES)
foreach(FILE IN ITEMS
        lsm_field_extension1d.f
        lsm_field_extension2d.f
        lsm_field_extension3d.f
       )
    list(APPEND LSM_FIELD_EXTENSION_SOURCE_FILES
         "field_extension/${FILE}")
endforeach()
set(LSM_FIELD_EXTENSION_SOURCE_FILES
    ${LSM_FIELD_EXTENSION_SOURCE_FILES} PARENT_SCOPE)

# --- Install parameters

# Header files
set(LSM_FIELD_EXTENSION_HEADER_FILES)
foreach(FILE IN ITEMS
        lsm_field_extension1d.h
        lsm_field_extension2d.h
        lsm_field_extension3d.h
       )
    list(APPEND LSM_FIELD_EXTENSION_HEADER_FILES
         "field_extension/${FILE}")
endforeach()
set(LSM_FIELD_EXTENSION_HEADER_FILES
    ${LSM_FIELD_EXTENSION_HEADER_FILES} PARENT_SCOPE)
