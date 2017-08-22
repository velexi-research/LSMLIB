# --------------------------------------------------------------
# Customizations for parallel advection-driven motion example
# --------------------------------------------------------------

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
add_executable(advection ${SOURCE_FILES})
add_dependencies(advection lsm_parallel samrai)

# Set search path for header files
target_include_directories(advection
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/fortran
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/parallel
    PUBLIC ${Boost_INCLUDE_DIR}
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
    PUBLIC ${SAMRAI_INCLUDE_DIR})

# Set link libraries
target_link_libraries(advection
    lsm_parallel
    lsm_toolbox
    ${SAMRAI_LIBS}
    ${HDF5_C_LIBRARY_hdf5}
    ${MPI_LIBRARY})

# Add files to copy directly to build directory
configure_file(advection-2d.input advection-2d.input COPYONLY)
