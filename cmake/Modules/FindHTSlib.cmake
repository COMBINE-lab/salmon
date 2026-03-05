include(FindPackageHandleStandardArgs)

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_HTSLIB QUIET htslib)
endif()

find_path(HTSLIB_INCLUDE_DIR
  NAMES htslib/sam.h
  HINTS ${PC_HTSLIB_INCLUDE_DIRS})

find_library(HTSLIB_LIBRARY
  NAMES hts libhts
  HINTS ${PC_HTSLIB_LIBRARY_DIRS})

set(HTSLIB_VERSION ${PC_HTSLIB_VERSION})

find_package_handle_standard_args(HTSlib
  REQUIRED_VARS HTSLIB_INCLUDE_DIR HTSLIB_LIBRARY
  VERSION_VAR HTSLIB_VERSION)

if(HTSLIB_FOUND AND NOT TARGET HTSlib::HTSlib)
  add_library(HTSlib::HTSlib UNKNOWN IMPORTED)
  set_target_properties(HTSlib::HTSlib PROPERTIES
    IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}")
endif()
