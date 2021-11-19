# --------------------------------------------------------------
# Customizations for serial examples
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# --- Individual examples

# curvature-driven motion
add_subdirectory(curvature_driven_motion)

# fast marching method example
add_subdirectory(fast_marching_method)

# grid management example
add_subdirectory(grid_management)

# reinitialization example
add_subdirectory(reinitialization)

# --- All lsm_serial examples

add_custom_target(examples
                  DEPENDS curvature_driven_motion fast_marching_method
                          grid_management reinitialization)
