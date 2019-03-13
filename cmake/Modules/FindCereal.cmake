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

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cereal DEFAULT_MSG CEREAL_INCLUDE_DIR)

mark_as_advanced(CEREAL_INCLUDE_DIR)

if(CEREAL_FOUND)
  message(STATUS "Cereal found (include: ${CEREAL_INCLUDE_DIRS})")
endif(CEREAL_FOUND)
