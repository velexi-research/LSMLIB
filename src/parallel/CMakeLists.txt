# --------------------------------------------------------------
# Customizations for lsm_parallel library
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Search path for header files
include_directories(${CMAKE_SOURCE_DIR}/src/toolbox)
include_directories(${CMAKE_SOURCE_DIR}/src/serial)
include_directories(${BOOST_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_DIR})
include_directories(${SAMRAI_INCLUDE_DIR})

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Source files
file(GLOB SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
     ${CMAKE_CURRENT_SOURCE_DIR}/*.f)

# Library
add_library(lsm_parallel ${SOURCE_FILES})

# --------------------------------------------------------------
# Installation
# --------------------------------------------------------------

# Header files
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
install(FILES ${HEADER_FILES} DESTINATION include/lsmlib)

# Library files
install(TARGETS lsm_parallel DESTINATION lib)

# --------------------------------------------------------------
# Build preparations
# --------------------------------------------------------------

# Create build output directories
file(MAKE_DIRECTORY ${LIBRARY_OUTPUT_PATH})
