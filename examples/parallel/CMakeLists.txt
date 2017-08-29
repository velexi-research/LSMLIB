# --------------------------------------------------------------
# Customizations for parallel examples
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# advection
add_subdirectory(advection/2d)
add_subdirectory(advection/3d)

# curvature-driven motion
add_subdirectory(curvature)

# field_extension/2d
# field_extension/3d

# toolbox/2d
# toolbox/3d

# normal_velocity_motion/2d
# normal_velocity_motion/3d

# vector_level_sets/2d

# orthogonalization

# --- All lsm_parallel examples

add_custom_target(examples-parallel
                  DEPENDS
                  advection-2d advection-3d)
