# --------------------------------------------------------------
# Customizations for toolbox tests
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# --- boundary conditions

add_executable(neumann_bc_ENO1_1d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_ENO1_1d.c)
target_include_directories(neumann_bc_ENO1_1d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_ENO1_1d lsm_toolbox m)

add_executable(neumann_bc_ENO1_2d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_ENO1_2d.c)
target_include_directories(neumann_bc_ENO1_2d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_ENO1_2d lsm_toolbox m)

add_executable(neumann_bc_ENO1_3d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_ENO1_3d.c)
target_include_directories(neumann_bc_ENO1_3d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_ENO1_3d lsm_toolbox m)

add_executable(neumann_bc_ENO2_1d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_ENO2_1d.c)
target_include_directories(neumann_bc_ENO2_1d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_ENO2_1d lsm_toolbox m)

add_executable(neumann_bc_ENO2_2d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_ENO2_2d.c)
target_include_directories(neumann_bc_ENO2_2d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_ENO2_2d lsm_toolbox m)

add_executable(neumann_bc_ENO2_3d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_ENO2_3d.c)
target_include_directories(neumann_bc_ENO2_3d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_ENO2_3d lsm_toolbox m)

add_executable(neumann_bc_ENO3_1d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_ENO3_1d.c)
target_include_directories(neumann_bc_ENO3_1d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_ENO3_1d lsm_toolbox m)

add_executable(neumann_bc_ENO3_2d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_ENO3_2d.c)
target_include_directories(neumann_bc_ENO3_2d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_ENO3_2d lsm_toolbox m)

add_executable(neumann_bc_ENO3_3d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_ENO3_3d.c)
target_include_directories(neumann_bc_ENO3_3d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_ENO3_3d lsm_toolbox m)

add_executable(neumann_bc_WENO5_1d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_WENO5_1d.c)
target_include_directories(neumann_bc_WENO5_1d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_WENO5_1d lsm_toolbox m)

add_executable(neumann_bc_WENO5_2d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_WENO5_2d.c)
target_include_directories(neumann_bc_WENO5_2d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_WENO5_2d lsm_toolbox m)

add_executable(neumann_bc_WENO5_3d
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/neumann_bc_WENO5_3d.c)
target_include_directories(neumann_bc_WENO5_3d
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(neumann_bc_WENO5_3d lsm_toolbox m)

# --- calculus toolbox

add_executable(test_calculus_toolbox test_calculus_toolbox.c)
target_include_directories(test_calculus_toolbox
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(test_calculus_toolbox lsm_toolbox m)

# --- fast marching method

add_executable(test_FMM_Heap test_FMM_Heap.c)
target_include_directories(test_FMM_Heap
    PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
target_link_libraries(test_FMM_Heap lsm_toolbox)
