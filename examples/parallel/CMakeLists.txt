# --------------------------------------------------------------
# Customizations for parallel examples
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

add_subdirectory(advection)

#field_extension        toolbox
# normal_velocity_motion vector_level_sets
#
#2d             3d
#
#toolbox
# normal_velocity_motion
#field_extension        orthogonalization

# --- All lsm_parallel examples

add_custom_target(examples-parallel
                  DEPENDS advection)
