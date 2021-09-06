# --------------------------------------------------------------
# Customizations for serial fast marching method examples
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Add 'computeDistanceFunction2d' executable
add_executable(computeDistanceFunction2d computeDistanceFunction2d.c)
target_include_directories(computeDistanceFunction2d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/serial)
target_link_libraries(computeDistanceFunction2d lsm_serial lsm_toolbox)

# Add 'computeExtensionFields2d' executable
add_executable(computeExtensionFields2d computeExtensionFields2d.c)
target_include_directories(computeExtensionFields2d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/serial)
target_link_libraries(computeExtensionFields2d lsm_serial lsm_toolbox)

# Add 'solveEikonalEquation2d' executable
add_executable(solveEikonalEquation2d solveEikonalEquation2d.c)
target_include_directories(solveEikonalEquation2d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/serial)
target_link_libraries(solveEikonalEquation2d lsm_serial lsm_toolbox)

# Add 'computeDistanceFunction3d' executable
add_executable(computeDistanceFunction3d computeDistanceFunction3d.c)
target_include_directories(computeDistanceFunction3d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/serial)
target_link_libraries(computeDistanceFunction3d lsm_serial lsm_toolbox)
