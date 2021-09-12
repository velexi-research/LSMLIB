# =============================================================================
# LSMLIBConfigExtras.cmake
# =============================================================================

# Set "exported" package variables
get_property(LSMLIB_INCLUDE_DIRS
             TARGET LSMLIB::lsm
             PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
get_property(LSMLIB_LIB_DIRS
             TARGET LSMLIB::lsm
             PROPERTY INTERFACE_LINK_DIRECTORIES)
