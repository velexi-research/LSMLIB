# -----------------------------------------------------------------------------
# Customizations for toolbox components
# -----------------------------------------------------------------------------

# Build requirements: No customizations
# Parameters (Install): No customizations

# -----------------------------------------------------------------------------
# Parameters (Build)
# -----------------------------------------------------------------------------

file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

# -----------------------------------------------------------------------------
# Targets (Build)
# -----------------------------------------------------------------------------

# Source files
file(GLOB LSM_TOOLBOX_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c
     ${CMAKE_CURRENT_SOURCE_DIR}/*.f
     ${CMAKE_CURRENT_BINARY_DIR}/*.f)
set(LSM_TOOLBOX_SOURCE_FILES ${LSM_TOOLBOX_SOURCE_FILES} PARENT_SCOPE)

# -----------------------------------------------------------------------------
# Targets (Install)
# -----------------------------------------------------------------------------

# Header files
install(FILES ${HEADER_FILES} DESTINATION ${INSTALL_HEADER_PATH})

# -----------------------------------------------------------------------------
# Build preparations
# -----------------------------------------------------------------------------

# --- Process template files

configure_file(lsm_calculus_toolbox2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_calculus_toolbox2d.f)
configure_file(lsm_calculus_toolbox2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_calculus_toolbox2d_local.f)
configure_file(lsm_calculus_toolbox3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_calculus_toolbox3d.f)

configure_file(lsm_level_set_evolution1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_level_set_evolution1d.f)
configure_file(lsm_level_set_evolution2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_level_set_evolution2d.f)
configure_file(lsm_level_set_evolution2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_level_set_evolution2d_local.f)
configure_file(lsm_level_set_evolution3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_level_set_evolution3d.f)
configure_file(lsm_level_set_evolution3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_level_set_evolution3d_local.f)

configure_file(lsm_math_utils1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_math_utils1d.f)
configure_file(lsm_math_utils2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_math_utils2d.f)
configure_file(lsm_math_utils2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_math_utils2d_local.f)
configure_file(lsm_math_utils3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_math_utils3d.f)
configure_file(lsm_math_utils3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_math_utils3d_local.f)

configure_file(lsm_spatial_derivatives1d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_spatial_derivatives1d.f)
configure_file(lsm_spatial_derivatives2d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_spatial_derivatives2d.f)
configure_file(lsm_spatial_derivatives2d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_spatial_derivatives2d_local.f)
configure_file(lsm_spatial_derivatives3d.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_spatial_derivatives3d.f)
configure_file(lsm_spatial_derivatives3d_local.f.in
               ${CMAKE_CURRENT_BINARY_DIR}/lsm_spatial_derivatives3d_local.f)

# --- Copy header files to binary folder

file(COPY ${HEADER_FILES} DESTINATION ${HEADER_OUTPUT_PATH})
