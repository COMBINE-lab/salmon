# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/zzare/salmon/external/libgff-2.0.0"
  "/home/zzare/salmon/external/libgff-2.0.0/build"
  "/home/zzare/salmon/external/install"
  "/home/zzare/salmon/build/libgff-prefix/tmp"
  "/home/zzare/salmon/build/libgff-prefix/src/libgff-stamp"
  "/home/zzare/salmon/external"
  "/home/zzare/salmon/build/libgff-prefix/src/libgff-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/zzare/salmon/build/libgff-prefix/src/libgff-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/zzare/salmon/build/libgff-prefix/src/libgff-stamp${cfgdir}") # cfgdir has leading slash
endif()
