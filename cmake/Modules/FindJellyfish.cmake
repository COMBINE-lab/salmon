###############################################################################
# Find Jellyfish 
#
# This sets the following variables:
# JELLYFISH_FOUND - True if Jellyfish was found.
# JELLYFISH_INCLUDE_DIRS - Directories containing the Jellyfish include files.
# JELLYFISH_DEFINITIONS - Compiler flags for Jellyfish.

find_path(JELLYFISH_INCLUDE_DIR jellyfish
	HINTS "${JELLYFISH_ROOT}/include" "$ENV{JELLYFISH_ROOT}/include" "/usr/include" "$ENV{PROGRAMFILES}/jellyfish/include")

set(JELLYFISH_INCLUDE_DIRS ${JELLYFISH_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
#message("Required Jellyfish version ${Jellyfish_FIND_VERSION}")
find_package_handle_standard_args(Jellyfish 
                                 DEFAULT_MSG 
                                 FOUND_VAR JELLYFISH_FOUND 
                                  REQUIRED_VARS JELLYFISH_INCLUDE_DIR 
                                  VERSION_VAR Jellyfish_FOUND_VERSION)

mark_as_advanced(JELLYFISH_INCLUDE_DIR)

if(JELLYFISH_FOUND)
    message(STATUS "Jellyfish found (include: ${JELLYFISH_INCLUDE_DIRS})")
endif(JELLYFISH_FOUND)
