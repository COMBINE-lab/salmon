###############################################################################
# Find Cereal
#
# This sets the following variables:
# CEREAL_FOUND - True if Cereal was found.
# CEREAL_INCLUDE_DIRS - Directories containing the Cereal include files.
# CEREAL_DEFINITIONS - Compiler flags for Cereal.

find_path(CEREAL_INCLUDE_DIR cereal
	HINTS "${CEREAL_ROOT}/include" "$ENV{CEREAL_ROOT}/include" "/usr/include" "$ENV{PROGRAMFILES}/cereal/include")

set(CEREAL_INCLUDE_DIRS ${CEREAL_INCLUDE_DIR})

if(CEREAL_INCLUDE_DIR)
  set(CEREAL_FOUND YES)
  set(CEREAL_VERSION_MAJOR 0)
  set(CEREAL_VERSION_MINOR 0)
  set(CEREAL_VERSION_PATCH 0)
  if(EXISTS "${CEREAL_INCLUDE_DIR}/cereal/version.hpp")
    # Read and parse cereal version header file for version number
    file(READ "${CEREAL_INCLUDE_DIR}/cereal/version.hpp"
        _CEREAL_HEADER_CONTENTS)
    string(REGEX REPLACE ".*#define CEREAL_VERSION_MAJOR ([0-9]+).*" "\\1"
        CEREAL_VERSION_MAJOR "${_CEREAL_HEADER_CONTENTS}")
    string(REGEX REPLACE ".*#define CEREAL_VERSION_MINOR ([0-9]+).*" "\\1"
        CEREAL_VERSION_MINOR "${_CEREAL_HEADER_CONTENTS}")
    string(REGEX REPLACE ".*#define CEREAL_VERSION_PATCH ([0-9]+).*" "\\1"
        CEREAL_VERSION_PATCH "${_CEREAL_HEADER_CONTENTS}")
    set(CEREAL_VERSION_STRING "${CEREAL_VERSION_MAJOR}.${CEREAL_VERSION_MINOR}.${CEREAL_VERSION_PATCH}")
  else()
    set(CEREAL_FOUND NO)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cereal 
  REQUIRED_VARS CEREAL_INCLUDE_DIR
  VERSION_VAR CEREAL_VERSION_STRING)

mark_as_advanced(CEREAL_INCLUDE_DIR)

if(CEREAL_FOUND)
  message(STATUS "cereal found (include: ${CEREAL_INCLUDE_DIRS})")
endif(CEREAL_FOUND)
