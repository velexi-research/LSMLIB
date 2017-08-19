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

set(BC_PROGRAMS
    neumann_bc_ENO1_1d
    neumann_bc_ENO1_2d
    neumann_bc_ENO1_3d
    neumann_bc_ENO2_1d
    neumann_bc_ENO2_2d
    neumann_bc_ENO2_3d
    neumann_bc_ENO3_1d
    neumann_bc_ENO3_2d
    neumann_bc_ENO3_3d
    neumann_bc_WENO5_1d
    neumann_bc_WENO5_2d
    neumann_bc_WENO5_3d)

foreach (PROG ${BC_PROGRAMS})
    add_executable(${PROG}
        ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/${PROG}.c)
    target_include_directories(${PROG}
        PUBLIC ${CMAKE_SOURCE_DIR}/src/toolbox)
    target_link_libraries(${PROG} lsm_toolbox)
endforeach (PROG)

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

# --- All lsm_toolbox examples

add_custom_target(examples-toolbox
                  DEPENDS ${BC_PROGRAMS} test_calculus_toolbox test_FMM_Heap)
