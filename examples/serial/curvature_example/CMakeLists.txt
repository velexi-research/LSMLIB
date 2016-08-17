# --------------------------------------------------------------
# Customizations for serial curvature examples
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
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

# Add executable
add_executable(curvature_example ${SOURCE_FILES})

# Search path for header files
target_include_directories(curvature_example
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/serial)

target_link_libraries(curvature_example lsm_serial lsm_toolbox m)
