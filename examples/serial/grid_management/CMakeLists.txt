# --------------------------------------------------------------
# Customizations for serial grid management example
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
add_executable(grid_management ${SOURCE_FILES})

# Search path for header files
target_include_directories(grid_management
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/serial)

target_link_libraries(grid_management lsm_serial lsm_toolbox)
