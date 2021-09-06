# --------------------------------------------------------------
# Customizations for serial curvature-driven motion example
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Set header and library search paths
include_directories(${CMAKE_SOURCE_DIR}/src/toolbox
                    ${CMAKE_SOURCE_DIR}/src/serial)

# Source files
file(GLOB SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

# Add executable
add_executable(curvature_driven_motion ${SOURCE_FILES})
target_link_libraries(curvature_driven_motion lsm_serial lsm_toolbox)
