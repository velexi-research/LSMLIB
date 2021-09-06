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
target_link_libraries(lsm_serial lsm_toolbox)

# --------------------------------------------------------------
# Installation
# --------------------------------------------------------------

# Header files
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
install(FILES ${HEADER_FILES} DESTINATION ${INSTALL_INCLUDE_DIR})

# Library files
install(TARGETS lsm_serial DESTINATION ${INSTALL_LIB_DIR})

# --------------------------------------------------------------
# Build preparations
# --------------------------------------------------------------

# Create build output directories
file(MAKE_DIRECTORY ${LIBRARY_OUTPUT_PATH})
