# =============================================================================
# LSMLIB Example: grid management example
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

add_executable(example-grid_management grid_management.c)
add_dependencies(example-grid_management LSMLIB::lsm)
target_link_libraries(example-grid_management PRIVATE lsm)

# -----------------------------------------------------------------------------
# Custom Targets
# -----------------------------------------------------------------------------
add_custom_target(clean-grid
    COMMAND rm -f grid_2d_demo.*
    COMMAND rm -f grid_2d_demo.binary.*
    COMMAND rm -f grid_3d_demo.*
    COMMAND rm -f grid_3d_demo.binary.*
)
