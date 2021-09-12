# LSMLIBConfig.cmake

include(${CMAKE_CURRENT_LIST_DIR}/LSMLIBConfigVersion.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/LSMLIBTargets.cmake)

get_property(LSMLIB_INCLUDE_DIRS
             TARGET LSMLIB::lsm
             PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
get_property(LSMLIB_LIB_DIRS
             TARGET LSMLIB::lsm
             PROPERTY INTERFACE_LINK_DIRECTORIES)
