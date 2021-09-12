# =============================================================================
# toolbox components
# =============================================================================

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

# --- Build parameters

# Source files
set(LSM_TOOLBOX_SOURCE_FILES)
foreach(FILE IN ITEMS
        lsm_initialization2d.c
        lsm_initialization3d.c
        lsm_calculus_toolbox.f
        lsm_localization2d.f
        lsm_localization3d.f
        lsm_tvd_runge_kutta1d.f
        lsm_tvd_runge_kutta2d.f
        lsm_tvd_runge_kutta2d_local.f
        lsm_tvd_runge_kutta3d.f
        lsm_tvd_runge_kutta3d_local.f
       )
    list(APPEND LSM_TOOLBOX_SOURCE_FILES "toolbox/${FILE}")
endforeach()

foreach(FILE IN ITEMS
        lsm_calculus_toolbox2d.f
        lsm_calculus_toolbox2d_local.f
        lsm_calculus_toolbox3d.f
        lsm_level_set_evolution1d.f
        lsm_level_set_evolution2d.f
        lsm_level_set_evolution2d_local.f
        lsm_level_set_evolution3d.f
        lsm_level_set_evolution3d_local.f
        lsm_math_utils1d.f
        lsm_math_utils2d.f
        lsm_math_utils2d_local.f
        lsm_math_utils3d.f
        lsm_math_utils3d_local.f
        lsm_spatial_derivatives1d.f
        lsm_spatial_derivatives2d.f
        lsm_spatial_derivatives2d_local.f
        lsm_spatial_derivatives3d.f
        lsm_spatial_derivatives3d_local.f
       )
        list(APPEND LSM_TOOLBOX_SOURCE_FILES
             "${CMAKE_CURRENT_BINARY_DIR}/${FILE}")
endforeach()
set(LSM_TOOLBOX_SOURCE_FILES ${LSM_TOOLBOX_SOURCE_FILES} PARENT_SCOPE)

# --- Install parameters

# Header files
set(LSM_TOOLBOX_HEADER_FILES)
foreach(FILE IN ITEMS
        lsm_calculus_toolbox.h
        lsm_calculus_toolbox2d.h
        lsm_calculus_toolbox2d_local.h
        lsm_calculus_toolbox3d.h
        lsm_initialization2d.h
        lsm_initialization3d.h
        lsm_level_set_evolution1d.h
        lsm_level_set_evolution2d.h
        lsm_level_set_evolution2d_local.h
        lsm_level_set_evolution3d.h
        lsm_level_set_evolution3d_local.h
        lsm_localization2d.h
        lsm_localization3d.h
        lsm_math_utils1d.h
        lsm_math_utils2d.h
        lsm_math_utils2d_local.h
        lsm_math_utils3d.h
        lsm_math_utils3d_local.h
        lsm_spatial_derivatives1d.h
        lsm_spatial_derivatives2d.h
        lsm_spatial_derivatives2d_local.h
        lsm_spatial_derivatives3d.h
        lsm_spatial_derivatives3d_local.h
        lsm_tvd_runge_kutta1d.h
        lsm_tvd_runge_kutta2d.h
        lsm_tvd_runge_kutta2d_local.h
        lsm_tvd_runge_kutta3d.h
        lsm_tvd_runge_kutta3d_local.h
       )
    list(APPEND LSM_TOOLBOX_HEADER_FILES "toolbox/${FILE}")
endforeach()
set(LSM_TOOLBOX_HEADER_FILES ${LSM_TOOLBOX_HEADER_FILES} PARENT_SCOPE)

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

# --- Process template files
#
#     Note: requires `lsmlib_zero_tol` and `tiny_nonzero_number` variables
#           to be set.

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
