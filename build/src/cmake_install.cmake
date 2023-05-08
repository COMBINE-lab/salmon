# Install script for directory: /home/zzare/salmon/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/zzare/salmon")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES
    "/home/zzare/salmon/external/install/lib/libtbb.so"
    "/home/zzare/salmon/external/install/lib/libtbb.so.12"
    "/home/zzare/salmon/external/install/lib/libtbb.so.12.5"
    "/home/zzare/salmon/external/install/lib/libtbbbind_2_0.so"
    "/home/zzare/salmon/external/install/lib/libtbbbind_2_0.so.3"
    "/home/zzare/salmon/external/install/lib/libtbbbind_2_0.so.3.5"
    "/home/zzare/salmon/external/install/lib/libtbbmalloc.so"
    "/home/zzare/salmon/external/install/lib/libtbbmalloc.so.2"
    "/home/zzare/salmon/external/install/lib/libtbbmalloc.so.2.5"
    "/home/zzare/salmon/external/install/lib/libtbbmalloc_proxy.so"
    "/home/zzare/salmon/external/install/lib/libtbbmalloc_proxy.so.2"
    "/home/zzare/salmon/external/install/lib/libtbbmalloc_proxy.so.2.5"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/zzare/salmon/build/src/salmon")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/salmon" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/salmon")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/salmon")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/zzare/salmon/build/src/libsalmon_core.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  
    execute_process(COMMAND "/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake"
                            -DCMAKE_SYSTEM_NAME=Linux
                            -DCMAKE_INSTALL_PREFIX=/home/zzare/salmon
                            -P "/home/zzare/salmon/cmake/PostInstall.cmake")
    
endif()

