# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/zzare/salmon/external/jemalloc-5.2.1"
  "/home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-build"
  "/home/zzare/salmon/external/install"
  "/home/zzare/salmon/build/libjemalloc-prefix/tmp"
  "/home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp"
  "/home/zzare/salmon/external"
  "/home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp${cfgdir}") # cfgdir has leading slash
endif()
