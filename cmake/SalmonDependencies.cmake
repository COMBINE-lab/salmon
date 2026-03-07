include_guard(GLOBAL)

include(ExternalProject)
include(FetchContent)

set(SALMON_DEPS_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}/_deps/local"
    CACHE PATH "Install prefix for dependency fallbacks")

set(SALMON_PUFFERFISH_GIT_REPOSITORY
    "https://github.com/COMBINE-lab/pufferfish.git"
    CACHE STRING "Git repository used when fetching pufferfish")
set(SALMON_PUFFERFISH_GIT_TAG
    "5ec513ba07fce97636867b3ed2c4351706a7e1f8"
    CACHE STRING "Immutable git commit used when fetching pufferfish")
set(SALMON_PUFFERFISH_SOURCE_DIR
    ""
    CACHE PATH "Optional local pufferfish source checkout to use instead of fetching")
set(SALMON_FQFEEDER_GIT_REPOSITORY
    "https://github.com/rob-p/FQFeeder.git"
    CACHE STRING "Git repository used when fetching FQFeeder")
set(SALMON_FQFEEDER_GIT_TAG
    "08eeed92a50902d670e5c064b4e3015ebb8bd1cc"
    CACHE STRING "Immutable git commit used when fetching FQFeeder")
set(SALMON_FQFEEDER_SOURCE_DIR
    ""
    CACHE PATH "Optional local FQFeeder source checkout to use instead of fetching")

if(DEFINED CUSTOM_BOOST_PATH)
  set(CMAKE_INCLUDE_PATH ${CUSTOM_BOOST_PATH} ${CMAKE_INCLUDE_PATH})
  set(CMAKE_LIBRARY_PATH ${CUSTOM_BOOST_PATH}/lib ${CMAKE_LIBRARY_PATH})
endif()

if(USE_SHARED_LIBS)
  set(Boost_USE_STATIC_LIBS OFF)
else()
  set(Boost_USE_STATIC_LIBS ON)
endif()

set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_DEBUG_RUNTIME OFF)

function(salmon_pick_existing_target out_var)
  foreach(candidate IN LISTS ARGN)
    if(TARGET "${candidate}")
      set(${out_var} "${candidate}" PARENT_SCOPE)
      return()
    endif()
  endforeach()
  set(${out_var} "" PARENT_SCOPE)
endfunction()

set(SALMON_ZLIB_INCLUDE_DIRS "")
set(SALMON_ZLIB_LIBRARIES "")
find_package(ZLIBNG QUIET)
if(ZLIBNG_FOUND)
  message(STATUS "Found zlib-ng compatibility backend: ${ZLIBNG_LIBRARY}")
  set(SALMON_ZLIB_INCLUDE_DIRS ${ZLIBNG_INCLUDE_DIR})
  set(SALMON_ZLIB_LIBRARIES ZLIBNG::ZLIBNG)
  set(ZLIB_INCLUDE_DIR ${ZLIBNG_INCLUDE_DIR})
  set(ZLIB_LIBRARY ZLIBNG::ZLIBNG)
elseif(SALMON_FETCH_MISSING_DEPS)
  message(STATUS "zlib-ng not found; fetching pinned zlib-ng release in compatibility mode")
  set(ZLIB_COMPAT ON CACHE BOOL "" FORCE)
  set(ZLIB_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
  set(WITH_GTEST OFF CACHE BOOL "" FORCE)
  set(WITH_FUZZERS OFF CACHE BOOL "" FORCE)
  set(WITH_BENCHMARKS OFF CACHE BOOL "" FORCE)
  set(WITH_BENCHMARK_APPS OFF CACHE BOOL "" FORCE)
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(salmon_zlibng
    URL https://github.com/zlib-ng/zlib-ng/archive/refs/tags/2.2.5.tar.gz
    URL_HASH SHA256=5b3b022489f3ced82384f06db1e13ba148cbce38c7941e424d6cb414416acd18
  )
  FetchContent_MakeAvailable(salmon_zlibng)
  salmon_pick_existing_target(SALMON_ZLIB_TARGET
    zlib
    zlibstatic
    ZLIB::ZLIB
    ZLIBNG::ZLIBNG)
  if(NOT SALMON_ZLIB_TARGET)
    message(FATAL_ERROR "Fetched zlib-ng did not expose an expected CMake target.")
  endif()
  set(FETCHED_ZLIBNG TRUE)
  set(SALMON_ZLIB_INCLUDE_DIRS
      ${salmon_zlibng_SOURCE_DIR}
      ${salmon_zlibng_BINARY_DIR})
  set(SALMON_ZLIB_LIBRARIES ${SALMON_ZLIB_TARGET})
  set(ZLIB_INCLUDE_DIR ${SALMON_ZLIB_INCLUDE_DIRS})
  set(ZLIB_LIBRARY ${SALMON_ZLIB_TARGET})
else()
  message(FATAL_ERROR "zlib-ng is required. Install a zlib-ng compatibility package or enable SALMON_FETCH_MISSING_DEPS.")
endif()

if("${CMAKE_INCLUDE_PATH}" STREQUAL "")
  set(EXTRA_CMAKE_INCLUDE_FLAGS "")
else()
  set(EXTRA_CMAKE_INCLUDE_FLAGS "-I${CMAKE_INCLUDE_PATH}")
endif()

if("${CMAKE_LIBRARY_PATH}" STREQUAL "")
  set(EXTRA_CMAKE_LIBRARY_FLAGS "")
else()
  set(EXTRA_CMAKE_LIBRARY_FLAGS "-L${CMAKE_LIBRARY_PATH}")
endif()

find_package(Iconv REQUIRED)
if(NOT Iconv_IS_BUILT_IN)
  set(ICONV_LIB Iconv::Iconv)
endif()

find_package(LibLZMA)
if(NOT LIBLZMA_FOUND)
  message(STATUS "Will attempt to fetch and build liblzma")
  externalproject_add(liblzma
    PREFIX ${CMAKE_BINARY_DIR}/_deps/liblzma
    URL https://tukaani.org/xz/xz-5.2.2.tar.gz
    URL_HASH SHA256=73df4d5d34f0468bd57d09f2d8af363e95ed6cc3a4a86129d2f2c366259902a2
    SOURCE_SUBDIR .
    INSTALL_DIR ${SALMON_DEPS_INSTALL_PREFIX}
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS=${EXTRA_CMAKE_INCLUDE_FLAGS} CPPFLAGS=${EXTRA_CMAKE_INCLUDE_FLAGS} LDFLAGS=${EXTRA_CMAKE_LIBRARY_FLAGS}
    BUILD_COMMAND make ${QUIET_MAKE}
    INSTALL_COMMAND make ${QUIET_MAKE} install
  )
  set(LIBLZMA_LIBRARIES ${SALMON_DEPS_INSTALL_PREFIX}/lib/liblzma.a)
  set(FETCHED_LIBLZMA TRUE)
else()
  message(STATUS "Found liblzma library: ${LIBLZMA_LIBRARIES}")
endif()

find_package(BZip2)
if(NOT BZIP2_FOUND)
  message(STATUS "Will attempt to fetch and build libbz2")
  externalproject_add(libbz2
    PREFIX ${CMAKE_BINARY_DIR}/_deps/libbz2
    URL https://sourceware.org/pub/bzip2/bzip2-1.0.6.tar.gz
    URL_HASH SHA256=a2848f34fcd5d6cf47def00461fcb528a0484d8edef8208d6d2e2909dc61d9cd
    SOURCE_SUBDIR .
    INSTALL_DIR ${SALMON_DEPS_INSTALL_PREFIX}
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make ${QUIET_MAKE} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    INSTALL_COMMAND make ${QUIET_MAKE} install PREFIX=<INSTALL_DIR>
  )
  set(BZIP2_LIBRARIES ${SALMON_DEPS_INSTALL_PREFIX}/lib/libbz2.a)
  set(FETCHED_LIBBZ2 TRUE)
else()
  message(STATUS "Found libbz2 library: ${BZIP2_LIBRARIES}")
endif()

if(FETCH_BOOST OR BOOST_RECONFIGURE OR BOOST_WILL_RECONFIGURE)
  message(FATAL_ERROR
    "Legacy Boost knobs (FETCH_BOOST / BOOST_RECONFIGURE / BOOST_WILL_RECONFIGURE) are no longer supported. "
    "Install Boost through your package manager and reconfigure.")
endif()

find_package(Boost 1.59.0 REQUIRED COMPONENTS iostreams system filesystem timer chrono program_options)
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries: ${Boost_LIBRARIES}")

set(EXTERNAL_LIBRARY_PATH $CMAKE_CURRENT_SOURCE_DIR/lib)

find_package(cereal "1.3.2" QUIET)
if(NOT CEREAL_FOUND)
  message(STATUS "Build system will fetch and build the cereal serialization library")
  set(BUILD_DOC OFF CACHE BOOL "" FORCE)
  set(BUILD_SANDBOX OFF CACHE BOOL "" FORCE)
  set(SKIP_PERFORMANCE_COMPARISON ON CACHE BOOL "" FORCE)
  set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(CEREAL_INSTALL OFF CACHE BOOL "" FORCE)
  set(JUST_INSTALL_CEREAL ON CACHE BOOL "" FORCE)
  FetchContent_Declare(salmon_cereal
    URL https://github.com/USCiLab/cereal/archive/refs/tags/v1.3.2.tar.gz
    URL_HASH SHA256=16a7ad9b31ba5880dac55d62b5d6f243c3ebc8d46a3514149e56b5e7ea81f85f
  )
  FetchContent_MakeAvailable(salmon_cereal)
  set(CEREAL_INCLUDE_DIRS ${salmon_cereal_SOURCE_DIR}/include)
  set(FETCHED_CEREAL TRUE)
elseif(CEREAL_INCLUDE_DIRS AND NOT TARGET cereal::cereal)
  add_library(cereal INTERFACE)
  target_include_directories(cereal INTERFACE ${CEREAL_INCLUDE_DIRS})
  add_library(cereal::cereal ALIAS cereal)
endif()

if(TBB_RECONFIGURE OR TBB_WILL_RECONFIGURE)
  message(FATAL_ERROR
    "Legacy TBB knobs (TBB_RECONFIGURE / TBB_WILL_RECONFIGURE) are no longer supported. "
    "Install oneTBB through your package manager and reconfigure.")
endif()

find_package(TBB 2021.4 REQUIRED COMPONENTS tbb)
if(TBB_VERSION VERSION_LESS 2021.4)
  message(FATAL_ERROR "Found TBB version ${TBB_VERSION}, but Salmon requires >= 2021.4.")
endif()
set(TBB_TARGET_EXISTED TRUE)
get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Found suitable TBB version: ${TBB_VERSION}")


set(PUFFERFISH_EMBEDDED ON CACHE BOOL "" FORCE)
set(BUILD_PUFF_FOR_SALMON ON CACHE BOOL "" FORCE)
if(SALMON_PUFFERFISH_SOURCE_DIR)
  FetchContent_Declare(salmon_pufferfish
    SOURCE_DIR "${SALMON_PUFFERFISH_SOURCE_DIR}"
  )
else()
  FetchContent_Declare(salmon_pufferfish
    GIT_REPOSITORY ${SALMON_PUFFERFISH_GIT_REPOSITORY}
    GIT_TAG ${SALMON_PUFFERFISH_GIT_TAG}
    GIT_SHALLOW FALSE
  )
endif()
FetchContent_MakeAvailable(salmon_pufferfish)
set(SALMON_PUFFERFISH_SOURCE_DIR "${salmon_pufferfish_SOURCE_DIR}" CACHE INTERNAL "" FORCE)
set(SALMON_PUFFERFISH_BINARY_DIR "${salmon_pufferfish_BINARY_DIR}" CACHE INTERNAL "" FORCE)

if(SALMON_FQFEEDER_SOURCE_DIR)
  FetchContent_Declare(salmon_fqfeeder
    SOURCE_DIR "${SALMON_FQFEEDER_SOURCE_DIR}"
  )
else()
  FetchContent_Declare(salmon_fqfeeder
    GIT_REPOSITORY ${SALMON_FQFEEDER_GIT_REPOSITORY}
    GIT_TAG ${SALMON_FQFEEDER_GIT_TAG}
    GIT_SHALLOW FALSE
  )
endif()
FetchContent_GetProperties(salmon_fqfeeder)
if(NOT salmon_fqfeeder_POPULATED)
  if(POLICY CMP0169)
    cmake_policy(PUSH)
    cmake_policy(SET CMP0169 OLD)
  endif()
  FetchContent_Populate(salmon_fqfeeder)
  if(POLICY CMP0169)
    cmake_policy(POP)
  endif()
endif()
set(SALMON_FQFEEDER_SOURCE_DIR "${salmon_fqfeeder_SOURCE_DIR}" CACHE INTERNAL "" FORCE)

find_package(libgff 2.0.1 HINTS ${LIB_GFF_PATH} ${GFF_ROOT})
if(libgff_FOUND)
  if(GFF_INCLUDE_DIR)
    set(LIB_GFF_INCLUDE_DIR ${GFF_INCLUDE_DIR})
  endif()
  if(GFF_LIBRARY)
    get_filename_component(LIB_GFF_LIBRARY_DIR ${GFF_LIBRARY} DIRECTORY)
    if(NOT TARGET gff)
      add_library(gff UNKNOWN IMPORTED)
      set_target_properties(gff PROPERTIES
        IMPORTED_LOCATION "${GFF_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${LIB_GFF_INCLUDE_DIR}")
    endif()
  endif()
  message(STATUS "libgff ver. ${LIB_GFF_VERSION} found.")
  message(STATUS "    include: ${LIB_GFF_INCLUDE_DIR}")
  message(STATUS "    lib    : ${LIB_GFF_LIBRARY_DIR}")
endif()

if(NOT libgff_FOUND)
  message(STATUS "Build system will compile libgff")
  FetchContent_Declare(salmon_libgff
    URL https://github.com/COMBINE-lab/libgff/archive/v2.0.1.tar.gz
    URL_HASH SHA256=96d2bda64aaf9cf7b6c1a42205e408b0ef2a353ba42dad560db215e7ec105e2e
  )
  FetchContent_MakeAvailable(salmon_libgff)
  set(FETCHED_GFF TRUE)
  set(LIB_GFF_PATH ${salmon_libgff_SOURCE_DIR})
  set(LIB_GFF_INCLUDE_DIR ${salmon_libgff_SOURCE_DIR}/include)
  set(LIB_GFF_LIBRARY_DIR ${salmon_libgff_BINARY_DIR})
endif()

find_package(CURL)
find_package(HTSlib QUIET)
if(HTSlib_FOUND)
  message(STATUS "Found htslib: ${HTSLIB_LIBRARY}")
  set(HTSLIB_INCLUDE_DIR ${HTSLIB_INCLUDE_DIR})
  set(HTSLIB_LIBRARIES HTSlib::HTSlib)
  find_library(LIBDEFLATE_LIBRARY NAMES deflate libdeflate)
  if(LIBDEFLATE_LIBRARY)
    list(APPEND HTSLIB_LIBRARIES ${LIBDEFLATE_LIBRARY})
  endif()
elseif(SALMON_FETCH_MISSING_DEPS)
  message(STATUS "htslib not found; fetching pinned htslib release")
  get_filename_component(_salmon_zlib_link_dir "${salmon_zlibng_BINARY_DIR}" ABSOLUTE)
  set(_salmon_htslib_cppflags "")
  foreach(_zinc IN LISTS SALMON_ZLIB_INCLUDE_DIRS)
    string(APPEND _salmon_htslib_cppflags " -I${_zinc}")
  endforeach()
  set(_salmon_htslib_ldflags "-L${_salmon_zlib_link_dir}")
  externalproject_add(libhtslib
    PREFIX ${CMAKE_BINARY_DIR}/_deps/libhtslib
    URL https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2
    URL_HASH SHA256=6250c1df297db477516e60ac8df45ed75a652d1f25b0f37f12f5b17269eafde9
    SOURCE_SUBDIR .
    INSTALL_DIR ${SALMON_DEPS_INSTALL_PREFIX}
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR> --disable-libcurl CPPFLAGS=${_salmon_htslib_cppflags} LDFLAGS=${_salmon_htslib_ldflags} CC=${CMAKE_C_COMPILER}
    BUILD_COMMAND make ${QUIET_MAKE}
    INSTALL_COMMAND make ${QUIET_MAKE} install
  )
  if(FETCHED_ZLIBNG AND TARGET ${SALMON_ZLIB_TARGET})
    add_dependencies(libhtslib ${SALMON_ZLIB_TARGET})
  endif()
  if(FETCHED_LIBBZ2)
    add_dependencies(libhtslib libbz2)
  endif()
  if(FETCHED_LIBLZMA)
    add_dependencies(libhtslib liblzma)
  endif()
  set(FETCHED_HTSLIB TRUE)
  set(HTSLIB_INCLUDE_DIR ${SALMON_DEPS_INSTALL_PREFIX}/include)
  if(NOT TARGET HTSlib::HTSlib)
    add_library(HTSlib::HTSlib UNKNOWN IMPORTED)
    set_target_properties(HTSlib::HTSlib PROPERTIES
      IMPORTED_LOCATION "${SALMON_DEPS_INSTALL_PREFIX}/lib/libhts.a"
      INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}")
  endif()
  set(HTSLIB_LIBRARIES HTSlib::HTSlib)
else()
  message(FATAL_ERROR "htslib is required. Install htslib or enable SALMON_FETCH_MISSING_DEPS.")
endif()

if(ASAN_BUILD OR SALMON_USE_MIMALLOC STREQUAL "OFF")
  set(FAST_MALLOC_LIB "")
  set(FAST_MALLOC_OBJECT "")
  set(HAVE_FAST_MALLOC TRUE)
else()
  set(FAST_MALLOC_LIB "")
  set(FAST_MALLOC_OBJECT "")
  set(HAVE_FAST_MALLOC FALSE)
  find_package(Mimalloc QUIET)
  if(Mimalloc_FOUND)
    message(STATUS "Found mimalloc: ${MIMALLOC_LIBRARY}")
    set(FAST_MALLOC_LIB Mimalloc::Mimalloc)
    set(HAVE_FAST_MALLOC TRUE)
  endif()
endif()

if(NOT HAVE_FAST_MALLOC AND SALMON_FETCH_MISSING_DEPS)
  message(STATUS "mimalloc not found; fetching pinned mimalloc release")
  string(TOUPPER "${SALMON_MIMALLOC_OVERRIDE}" SALMON_MIMALLOC_OVERRIDE_UPPER)
  string(TOUPPER "${SALMON_MIMALLOC_OSX_ZONE}" SALMON_MIMALLOC_OSX_ZONE_UPPER)
  string(TOUPPER "${SALMON_MIMALLOC_OSX_INTERPOSE}" SALMON_MIMALLOC_OSX_INTERPOSE_UPPER)
  set(MI_OVERRIDE OFF CACHE BOOL "" FORCE)
  if(SALMON_MIMALLOC_OVERRIDE_UPPER STREQUAL "ON")
    set(MI_OVERRIDE ON CACHE BOOL "" FORCE)
  elseif(NOT SALMON_MIMALLOC_OVERRIDE_UPPER STREQUAL "OFF")
    message(FATAL_ERROR "SALMON_MIMALLOC_OVERRIDE must be ON or OFF")
  endif()
  if(APPLE)
    set(MI_OSX_ZONE OFF CACHE BOOL "" FORCE)
    if(SALMON_MIMALLOC_OSX_ZONE_UPPER STREQUAL "ON")
      set(MI_OSX_ZONE ON CACHE BOOL "" FORCE)
    elseif(NOT SALMON_MIMALLOC_OSX_ZONE_UPPER STREQUAL "OFF")
      message(FATAL_ERROR "SALMON_MIMALLOC_OSX_ZONE must be ON or OFF")
    endif()
    set(MI_OSX_INTERPOSE OFF CACHE BOOL "" FORCE)
    if(SALMON_MIMALLOC_OSX_INTERPOSE_UPPER STREQUAL "ON")
      set(MI_OSX_INTERPOSE ON CACHE BOOL "" FORCE)
    elseif(NOT SALMON_MIMALLOC_OSX_INTERPOSE_UPPER STREQUAL "OFF")
      message(FATAL_ERROR "SALMON_MIMALLOC_OSX_INTERPOSE must be ON or OFF")
    endif()
  endif()
  message(STATUS "mimalloc fetched config: MI_OVERRIDE=${MI_OVERRIDE} MI_OSX_ZONE=${MI_OSX_ZONE} MI_OSX_INTERPOSE=${MI_OSX_INTERPOSE}")
  set(MI_BUILD_SHARED OFF CACHE BOOL "" FORCE)
  set(MI_BUILD_STATIC ON CACHE BOOL "" FORCE)
  set(MI_BUILD_OBJECT ON CACHE BOOL "" FORCE)
  set(MI_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(salmon_mimalloc
    URL https://github.com/microsoft/mimalloc/archive/refs/tags/v3.2.8.tar.gz
    URL_HASH SHA256=68163666575518c213a6593850099adce3863b340ca2751103dbd1f253664e05
  )
  FetchContent_MakeAvailable(salmon_mimalloc)
  set(FETCHED_MIMALLOC TRUE)
  if(TARGET mimalloc-static)
    set(FAST_MALLOC_LIB mimalloc-static)
  endif()
  if(TARGET mimalloc-obj-target)
    set(FAST_MALLOC_OBJECT "${salmon_mimalloc_BINARY_DIR}/mimalloc${CMAKE_C_OUTPUT_EXTENSION}")
    set(FAST_MALLOC_LIB "")
    message(STATUS "Using mimalloc static override object: ${FAST_MALLOC_OBJECT}")
  endif()
  set(HAVE_FAST_MALLOC TRUE)
elseif(NOT HAVE_FAST_MALLOC AND SALMON_USE_MIMALLOC STREQUAL "REQUIRED")
  message(FATAL_ERROR "mimalloc is required but was not found. Install mimalloc or enable SALMON_FETCH_MISSING_DEPS.")
endif()
