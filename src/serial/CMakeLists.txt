# --------------------------------------------------------------
# Customizations for lsm_serial library
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Search path for header files
include_directories(${CMAKE_SOURCE_DIR}/src/toolbox)

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Source files
file(GLOB SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
list(REMOVE_ITEM SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/lsm_FMM_eikonal.c
     ${CMAKE_CURRENT_SOURCE_DIR}/lsm_FMM_field_extension.c)

# Library
add_library(lsm_serial ${SOURCE_FILES})

# --------------------------------------------------------------
# Installation
# --------------------------------------------------------------

# Header files
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
install(FILES ${HEADER_FILES} DESTINATION include/lsmlib)

# Library files
install(TARGETS lsm_serial DESTINATION lib)

# --------------------------------------------------------------
# Build preparations
# --------------------------------------------------------------

# Create build output directories
file(MAKE_DIRECTORY ${LIBRARY_OUTPUT_PATH})
