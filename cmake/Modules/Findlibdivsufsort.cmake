
find_path(SUFFARRAY_INCLUDE_DIR divsufsort64.h
	HINTS ${SUFFARRAY_ROOT} ENV SUFFARRAY_ROOT
  PATH_SUFFIXES include)

find_library(SUFFARRAY_LIBRARY NAMES divsufsort divsufsort64 libdivsufsort libdivsufsort64
  HINTS ${SUFFARRAY_ROOT} ENV SUFFARRAY_ROOT PATH_SUFFIXES lib lib64)

find_library(SUFFARRAY_LIBRARY64 NAMES divsufsort64 libdivsufsort64
  HINTS ${SUFFARRAY_ROOT} ENV SUFFARRAY_ROOT PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libdivsufsort DEFAULT_MSG
                                  SUFFARRAY_LIBRARY 
                                  SUFFARRAY_LIBRARY64
                                  SUFFARRAY_INCLUDE_DIR)

mark_as_advanced(SUFFARRAY_INCLUDE_DIR SUFFARRAY_LIBRARY)
