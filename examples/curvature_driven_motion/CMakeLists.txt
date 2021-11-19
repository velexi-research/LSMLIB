# =============================================================================
# LSMLIB Example: curvature-driven motion
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

add_executable(
    example-curvature_driven_motion
    curvature_model.c
    curvature_model_top.c
    curvature_model3d.c
    curvature_model3d_local.c
    lsm_options.c
)
add_dependencies(example-curvature_driven_motion LSMLIB::lsm)
target_link_libraries(example-curvature_driven_motion PRIVATE lsm)

# -----------------------------------------------------------------------------
# Custom Targets
# -----------------------------------------------------------------------------
add_custom_target(clean-curvature
    COMMAND rm -f data_final.gz data_init.gz grid.gz mask.gz out_file
)
