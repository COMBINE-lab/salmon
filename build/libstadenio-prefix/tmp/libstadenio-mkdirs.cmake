# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/zzare/salmon/external/staden-io_lib"
  "/home/zzare/salmon/build/libstadenio-prefix/src/libstadenio-build"
  "/home/zzare/salmon/external/install"
  "/home/zzare/salmon/build/libstadenio-prefix/tmp"
  "/home/zzare/salmon/build/libstadenio-prefix/src/libstadenio-stamp"
  "/home/zzare/salmon/external"
  "/home/zzare/salmon/build/libstadenio-prefix/src/libstadenio-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/zzare/salmon/build/libstadenio-prefix/src/libstadenio-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/zzare/salmon/build/libstadenio-prefix/src/libstadenio-stamp${cfgdir}") # cfgdir has leading slash
endif()
