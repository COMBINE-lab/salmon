# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/zzare/salmon/external/oneTBB-2021.5.0"
  "/home/zzare/salmon/build/libtbb-prefix/src/libtbb-build"
  "/home/zzare/salmon/external/install"
  "/home/zzare/salmon/build/libtbb-prefix/tmp"
  "/home/zzare/salmon/build/libtbb-prefix/src/libtbb-stamp"
  "/home/zzare/salmon/external"
  "/home/zzare/salmon/build/libtbb-prefix/src/libtbb-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/zzare/salmon/build/libtbb-prefix/src/libtbb-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/zzare/salmon/build/libtbb-prefix/src/libtbb-stamp${cfgdir}") # cfgdir has leading slash
endif()
