# =============================================================================
# LSMLIB Example: reinitialization example
# =============================================================================

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------

# LSMLIB
add_library(LSMLIB::lsm STATIC IMPORTED)
set_property(TARGET LSMLIB::lsm PROPERTY
             IMPORTED_LOCATION "${LIBRARY_OUTPUT_PATH}/liblsm.a"
)

# -----------------------------------------------------------------------------
# Programs
# -----------------------------------------------------------------------------

add_executable(example-reinitialization
               lsm_options.c
               lsm_options.h
               lsm_reinitialize.c
               lsm_reinitialization_medium2d.c
               lsm_reinitialization_medium2d.h
               )
add_dependencies(example-reinitialization LSMLIB::lsm)
target_link_libraries(example-reinitialization PRIVATE lsm)

# -----------------------------------------------------------------------------
# Custom Targets
# -----------------------------------------------------------------------------
add_custom_target(clean-reinit
    COMMAND rm -f data_init.gz data_final.gz grid.gz out_file
)
