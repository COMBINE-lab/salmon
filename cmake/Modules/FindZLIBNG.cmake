include(FindPackageHandleStandardArgs)

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_ZLIBNG QUIET zlib-ng zlib-ng-compat)
endif()

find_path(ZLIBNG_INCLUDE_DIR
  NAMES zlib-ng.h
  HINTS ${PC_ZLIBNG_INCLUDE_DIRS})

find_library(ZLIBNG_LIBRARY
  NAMES z zlib zlibstatic zlib-ng zlib-ng-compat
  HINTS ${PC_ZLIBNG_LIBRARY_DIRS})

set(ZLIBNG_VERSION ${PC_ZLIBNG_VERSION})

find_package_handle_standard_args(ZLIBNG
  REQUIRED_VARS ZLIBNG_INCLUDE_DIR ZLIBNG_LIBRARY
  VERSION_VAR ZLIBNG_VERSION)

if(ZLIBNG_FOUND AND NOT TARGET ZLIBNG::ZLIBNG)
  add_library(ZLIBNG::ZLIBNG UNKNOWN IMPORTED)
  set_target_properties(ZLIBNG::ZLIBNG PROPERTIES
    IMPORTED_LOCATION "${ZLIBNG_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${ZLIBNG_INCLUDE_DIR}")
endif()
