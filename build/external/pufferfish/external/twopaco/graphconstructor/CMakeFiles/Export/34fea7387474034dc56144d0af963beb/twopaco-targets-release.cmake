#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "twopaco" for configuration "RELEASE"
set_property(TARGET twopaco APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(twopaco PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libtwopaco.a"
  )

list(APPEND _cmake_import_check_targets twopaco )
list(APPEND _cmake_import_check_files_for_twopaco "${_IMPORT_PREFIX}/lib/libtwopaco.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
