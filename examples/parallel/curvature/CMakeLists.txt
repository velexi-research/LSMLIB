# -----------------------------------------------------------
# Customizations for parallel curvature-driven motion example
# -----------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Source files
file(GLOB SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
     ${CMAKE_CURRENT_SOURCE_DIR}/fortran/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/fortran/*.f)

# Add executable and dependencies
add_executable(curvature-parallel ${SOURCE_FILES})
add_dependencies(curvature-parallel lsm_parallel samrai)

# Set search path for header files
target_include_directories(curvature-parallel
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/fortran
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/parallel
    PUBLIC ${Boost_INCLUDE_DIR}
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
    PUBLIC ${SAMRAI_INCLUDE_DIR})

# Set link libraries
target_link_libraries(curvature-parallel
    lsm_parallel
    lsm_toolbox
    ${SAMRAI_LIBS}
    ${HDF5_C_LIBRARY_hdf5}
    ${MPI_LIBRARY})

# Add files to copy directly to build directory
configure_file(curvature-2d-circle.input
               curvature-2d-circle.input COPYONLY)
configure_file(curvature-2d-starfish.input
               curvature-2d-starfish.input COPYONLY)
configure_file(curvature-3d-sphere.input
               curvature-3d-sphere.input COPYONLY)
configure_file(curvature-3d-bumpy-sphere.input
               curvature-3d-bumpy-sphere.input COPYONLY)
