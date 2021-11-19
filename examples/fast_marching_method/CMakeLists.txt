# =============================================================================
# LSMLIB Example: fast-marching method
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

# computeDistanceFunction2d
add_executable(example-computeDistanceFunction2d computeDistanceFunction2d.c)
add_dependencies(example-computeDistanceFunction2d LSMLIB::lsm)
target_link_libraries(example-computeDistanceFunction2d PRIVATE lsm)

# computeDistanceFunction3d
add_executable(example-computeDistanceFunction3d computeDistanceFunction3d.c)
add_dependencies(example-computeDistanceFunction3d LSMLIB::lsm)
target_link_libraries(example-computeDistanceFunction3d PRIVATE lsm)

# computeExtensionFields2d
add_executable(example-computeExtensionFields2d computeExtensionFields2d.c)
add_dependencies(example-computeExtensionFields2d LSMLIB::lsm)
target_link_libraries(example-computeExtensionFields2d PRIVATE lsm)

# solveEikonalEquation2d
add_executable(example-solveEikonalEquation2d solveEikonalEquation2d.c)
add_dependencies(example-solveEikonalEquation2d LSMLIB::lsm)
target_link_libraries(example-solveEikonalEquation2d PRIVATE lsm)

# -----------------------------------------------------------------------------
# Custom Targets
# -----------------------------------------------------------------------------
add_custom_target(clean-fmm
    COMMAND rm -f computeDistanceFunction2d.dat
    COMMAND rm -f computeDistanceFunction3d.dat
    COMMAND rm -f computeExtensionFields2d.dat
    COMMAND rm -f solveEikonalEquation2d.dat
)
