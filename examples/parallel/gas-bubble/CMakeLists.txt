# -------------------------------------------------------------------
# Customizations for parallel gas bubble with surface-tension example
# -------------------------------------------------------------------

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
add_executable(gas-bubble-parallel ${SOURCE_FILES})
add_dependencies(gas-bubble-parallel lsm_parallel samrai)

# Set search path for header files
target_include_directories(gas-bubble-parallel
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/fortran
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/parallel
    PUBLIC ${Boost_INCLUDE_DIR}
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
    PUBLIC ${SAMRAI_INCLUDE_DIR})

# Set link libraries
target_link_libraries(gas-bubble-parallel
    lsm_parallel
    lsm_toolbox
    ${SAMRAI_LIBS}
    ${HDF5_C_LIBRARY_hdf5}
    ${MPI_LIBRARY})

# Add files to copy directly to build directory
configure_file(gas-bubble-2d-circle.input
               gas-bubble-2d-circle.input COPYONLY)
configure_file(gas-bubble-2d-starfish.input
               gas-bubble-2d-starfish.input COPYONLY)
configure_file(gas-bubble-3d-sphere.input
               gas-bubble-3d-sphere.input COPYONLY)
configure_file(gas-bubble-3d-bumpy-sphere.input
               gas-bubble-3d-bumpy-sphere.input COPYONLY)
