# From: https://raw.githubusercontent.com/STEllAR-GROUP/hpx/master/cmake/FindJemalloc.cmake
# Copyright (c)      2014 Thomas Heller
# Copyright (c) 2007-2012 Hartmut Kaiser
# Copyright (c) 2010-2011 Matt Anderson
# Copyright (c) 2011      Bryce Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

find_package(PkgConfig)
pkg_check_modules(PC_JEMALLOC QUIET libjemalloc)

find_path(JEMALLOC_INCLUDE_DIR jemalloc/jemalloc.h
  HINTS
    ${JEMALLOC_ROOT} ENV JEMALLOC_ROOT
    ${PC_JEMALLOC_MINIMAL_INCLUDEDIR}
    ${PC_JEMALLOC_MINIMAL_INCLUDE_DIRS}
    ${PC_JEMALLOC_INCLUDEDIR}
    ${PC_JEMALLOC_INCLUDE_DIRS}
  PATH_SUFFIXES include)

find_library(JEMALLOC_LIBRARY NAMES jemalloc libjemalloc
  HINTS
    ${JEMALLOC_ROOT} ENV JEMALLOC_ROOT
    ${PC_JEMALLOC_MINIMAL_LIBDIR}
    ${PC_JEMALLOC_MINIMAL_LIBRARY_DIRS}
    ${PC_JEMALLOC_LIBDIR}
    ${PC_JEMALLOC_LIBRARY_DIRS}
  PATH_SUFFIXES lib lib64)

if(JEMALLOC_INCLUDE_DIR)
  set(_version_regex "^#define[ \t]+JEMALLOC_VERSION[ \t]+\"([^\"]+)\".*")
  file(STRINGS "${JEMALLOC_INCLUDE_DIR}/jemalloc/jemalloc.h"
    JEMALLOC_VERSION REGEX "${_version_regex}")
  string(REGEX REPLACE "${_version_regex}" "\\1"
    JEMALLOC_VERSION "${JEMALLOC_VERSION}")
  unset(_version_regex)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set JEMALLOC_FOUND to TRUE
# if all listed variables are TRUE and the requested version matches.
find_package_handle_standard_args(Jemalloc REQUIRED_VARS
                                  JEMALLOC_LIBRARY JEMALLOC_INCLUDE_DIR
                                  VERSION_VAR JEMALLOC_VERSION)


if(JEMALLOC_FOUND)
  set(JEMALLOC_LIBRARIES    ${JEMALLOC_LIBRARY})
  set(JEMALLOC_INCLUDE_DIRS ${JEMALLOC_INCLUDE_DIR})
endif()

mark_as_advanced(JEMALLOC_INCLUDE_DIR JEMALLOC_LIBRARY)
