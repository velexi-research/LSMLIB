# --------------------------------------------------------------
# Customizations for serial reinitialization example
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
add_executable(reinitialization ${SOURCE_FILES})

# Search path for header files
target_include_directories(reinitialization
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/serial)

target_link_libraries(reinitialization lsm_serial lsm_toolbox)
