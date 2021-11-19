# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------

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

set(PROGRAM curvature_driven_motion)
add_executable(
    ${PROGRAM}
    curvature_model.c
    curvature_model_top.c
    curvature_model3d.c
    curvature_model3d_local.c
    lsm_options.c
)
add_dependencies(${PROGRAM} LSMLIB::lsm)
target_link_libraries(${PROGRAM} PRIVATE lsm)

# -----------------------------------------------------------------------------
# Custom Targets
# -----------------------------------------------------------------------------
add_custom_target(spotless
    COMMAND rm -f data_final* data_init* grid* mask* out_file*
)
