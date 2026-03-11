include(FindPackageHandleStandardArgs)

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_MIMALLOC QUIET mimalloc)
endif()

find_path(MIMALLOC_INCLUDE_DIR
  NAMES mimalloc.h
  HINTS ${PC_MIMALLOC_INCLUDE_DIRS})

find_library(MIMALLOC_LIBRARY
  NAMES mimalloc mimalloc-static
  HINTS ${PC_MIMALLOC_LIBRARY_DIRS})

set(MIMALLOC_VERSION ${PC_MIMALLOC_VERSION})

find_package_handle_standard_args(Mimalloc
  REQUIRED_VARS MIMALLOC_INCLUDE_DIR MIMALLOC_LIBRARY
  VERSION_VAR MIMALLOC_VERSION)

if(Mimalloc_FOUND AND NOT TARGET Mimalloc::Mimalloc)
  add_library(Mimalloc::Mimalloc UNKNOWN IMPORTED)
  set_target_properties(Mimalloc::Mimalloc PROPERTIES
    IMPORTED_LOCATION "${MIMALLOC_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MIMALLOC_INCLUDE_DIR}")
endif()
