# --------------------------------------------------------------
# NOTES: TODO explain what is going on in this file.
# --------------------------------------------------------------

# --------------------------------------------------------------
# CMake initialization
# --------------------------------------------------------------

cmake_minimum_required(VERSION 3.0)
project(samrai)
include(ExternalProject)

# --------------------------------------------------------------
# Configure SAMRAI installation as an external project
# --------------------------------------------------------------

ExternalProject_Add(samrai
    URL ${CMAKE_CURRENT_SOURCE_DIR}/samrai-mock-ext-pkg.tar.gz

    # Disable install step
    INSTALL_COMMAND ""
)

# SAMRAI library list
set(SAMRAI_LIB_LIST
    libSAMRAI_algs.a
    libSAMRAI_appu.a
    libSAMRAI_geom.a
    libSAMRAI_hier.a
    libSAMRAI_math.a
    libSAMRAI_mesh.a
    libSAMRAI_pdat.a
    libSAMRAI_solv.a
    libSAMRAI_xfer.a
    libSAMRAI_tbox.a)

# --------------------------------------------------------------
# Verify SAMRAI files
# --------------------------------------------------------------

# SAMRAI library directory
set(SAMRAI_LIB_DIR ${SAMRAI}/lib)

# Verify SAMRAI header files
find_file(SAMRAI_CONFIG_H SAMRAI_config.h ${SAMRAI}/include/SAMRAI)
if (SAMRAI_CONFIG_H)
    message("-- Found SAMRAI header files: TRUE (${SAMRAI_CONFIG_H})")
    set(SAMRAI_INCLUDE_DIR ${SAMRAI}/include)
else (SAMRAI_CONFIG_H)
    message(FATAL_ERROR "Could NOT find SAMRAI header files.")
endif (SAMRAI_CONFIG_H)

# Verify SAMRAI library files
foreach (LIB ${SAMRAI_LIB_LIST})
    # Clear TMP_VAR variable
    unset(TMP_VAR CACHE)

    find_library(TMP_VAR ${LIB} ${SAMRAI_LIB_DIR})

    if (TMP_VAR)
        message("-- Found ${LIB}: TRUE (${TMP_VAR})")
    else (TMP_VAR)
        message(FATAL_ERROR "Could NOT find ${LIB} library.")
    endif (TMP_VAR)

    # Add library to SAMRAI_LIBS
    set(SAMRAI_LIBS ${SAMRAI_LIBS} ${TMP_VAR})
endforeach (LIB)

# Clear TMP_VAR variable
unset(TMP_VAR CACHE)

# --------------------------------------------------------------
# Define SAMRAI package variables
# --------------------------------------------------------------

set(SAMRAI_INCLUDE_DIR ${SAMRAI}/include PARENT_SCOPE)
set(SAMRAI_LIB_DIR ${SAMRAI}/lib PARENT_SCOPE)
set(SAMRAI_LIBS ${SAMRAI_LIBS} PARENT_SCOPE)
