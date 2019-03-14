
find_path(SUFFARRAY_INCLUDE_DIR divsufsort.h
	HINTS ${SUFFARRAY_ROOT} ENV SUFFARRAY_ROOT
  PATH_SUFFIXES include)

find_library(SUFFARRAY_LIBRARY NAMES divsufsort libdivsufsort
  HINTS ${SUFFARRAY_ROOT} ENV SUFFARRAY_ROOT PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libdivsufsort DEFAULT_MSG
                                  SUFFARRAY_LIBRARY 
                                  SUFFARRAY_INCLUDE_DIR)

mark_as_advanced(SUFFARRAY_INCLUDE_DIR SUFFARRAY_LIBRARY)
