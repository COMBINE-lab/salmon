
find_path(GFF_INCLUDE_DIR gff.h
  HINTS ${GFF_ROOT} ENV GFF_ROOT
  PATH_SUFFIXES include)

find_library(GFF_LIBRARY NAMES gff libgff
  HINTS ${GFF_ROOT} ENV GFF_ROOT PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libgff DEFAULT_MSG 
                                  GFF_INCLUDE_DIR GFF_LIBRARY)

mark_as_advanced(GFF_INCLUDE_DIR GFF_LIBRARY)
