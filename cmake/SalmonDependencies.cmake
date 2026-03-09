include_guard(GLOBAL)

include(ExternalProject)
include(FetchContent)
include(CheckCXXSourceCompiles)

set(SALMON_DEPS_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}/_deps/local"
    CACHE PATH "Install prefix for dependency fallbacks")

set(SALMON_PUFFERFISH_GIT_REPOSITORY
    "https://github.com/COMBINE-lab/pufferfish.git"
    CACHE STRING "Git repository used when fetching pufferfish")
set(SALMON_PUFFERFISH_GIT_TAG
    "c5007eb02c22e4c783d778f3370dcd33ecb0df25"
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

set(SALMON_DEPENDENCY_PREFIX
    ""
    CACHE PATH "Optional dependency prefix used to seed CMAKE_PREFIX_PATH for find_package()")
if(SALMON_DEPENDENCY_PREFIX)
  list(PREPEND CMAKE_PREFIX_PATH "${SALMON_DEPENDENCY_PREFIX}")
  message(STATUS "Using SALMON_DEPENDENCY_PREFIX for package discovery: ${SALMON_DEPENDENCY_PREFIX}")
endif()

if(DEFINED CUSTOM_BOOST_PATH)
  message(DEPRECATION "CUSTOM_BOOST_PATH is deprecated; use SALMON_DEPENDENCY_PREFIX or CMAKE_PREFIX_PATH instead.")
  list(PREPEND CMAKE_PREFIX_PATH "${CUSTOM_BOOST_PATH}")
endif()

if(DEFINED USE_SHARED_LIBS)
  message(DEPRECATION "USE_SHARED_LIBS is deprecated; use SALMON_BOOST_USE_STATIC_LIBS instead.")
endif()

if(NOT DEFINED SALMON_BOOST_USE_STATIC_LIBS)
  if(DEFINED USE_SHARED_LIBS)
    if(USE_SHARED_LIBS)
      set(_salmon_boost_static_default OFF)
    else()
      set(_salmon_boost_static_default ON)
    endif()
  else()
    # Prefer static Boost by default.
    set(_salmon_boost_static_default ON)
  endif()
  set(SALMON_BOOST_USE_STATIC_LIBS
      ${_salmon_boost_static_default}
      CACHE BOOL "Prefer static Boost component libraries")
endif()
set(Boost_USE_STATIC_LIBS ${SALMON_BOOST_USE_STATIC_LIBS})
set(FETCHED_BOOST FALSE)

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
if(SALMON_USE_SYSTEM_DEPS)
  find_package(ZLIBNG QUIET)
endif()
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

find_package(Iconv REQUIRED)
if(NOT Iconv_IS_BUILT_IN)
  set(ICONV_LIB Iconv::Iconv)
endif()

find_package(LibLZMA)
if(NOT LIBLZMA_FOUND)
  if(SALMON_FETCH_MISSING_DEPS)
    message(STATUS "Will attempt to fetch and build liblzma")
    externalproject_add(liblzma
      PREFIX ${CMAKE_BINARY_DIR}/_deps/liblzma
      URL https://tukaani.org/xz/xz-5.2.2.tar.gz
      URL_HASH SHA256=73df4d5d34f0468bd57d09f2d8af363e95ed6cc3a4a86129d2f2c366259902a2
      SOURCE_SUBDIR .
      INSTALL_DIR ${SALMON_DEPS_INSTALL_PREFIX}
      BUILD_IN_SOURCE TRUE
      CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
      BUILD_COMMAND make
      INSTALL_COMMAND make install
    )
    set(LIBLZMA_LIBRARIES ${SALMON_DEPS_INSTALL_PREFIX}/lib/liblzma.a)
    set(FETCHED_LIBLZMA TRUE)
  else()
    message(FATAL_ERROR "liblzma is required by htslib but was not found. Install liblzma or enable SALMON_FETCH_MISSING_DEPS.")
  endif()
else()
  message(STATUS "Found liblzma library: ${LIBLZMA_LIBRARIES}")
endif()

find_package(BZip2)
if(NOT BZIP2_FOUND)
  if(SALMON_FETCH_MISSING_DEPS)
    message(STATUS "Will attempt to fetch and build libbz2")
    externalproject_add(libbz2
      PREFIX ${CMAKE_BINARY_DIR}/_deps/libbz2
      URL https://sourceware.org/pub/bzip2/bzip2-1.0.6.tar.gz
      URL_HASH SHA256=a2848f34fcd5d6cf47def00461fcb528a0484d8edef8208d6d2e2909dc61d9cd
      SOURCE_SUBDIR .
      INSTALL_DIR ${SALMON_DEPS_INSTALL_PREFIX}
      BUILD_IN_SOURCE TRUE
      CONFIGURE_COMMAND ""
      BUILD_COMMAND make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
      INSTALL_COMMAND make install PREFIX=<INSTALL_DIR>
    )
    set(BZIP2_LIBRARIES ${SALMON_DEPS_INSTALL_PREFIX}/lib/libbz2.a)
    set(FETCHED_LIBBZ2 TRUE)
  else()
    message(FATAL_ERROR "bzip2 is required by htslib but was not found. Install bzip2 or enable SALMON_FETCH_MISSING_DEPS.")
  endif()
else()
  message(STATUS "Found libbz2 library: ${BZIP2_LIBRARIES}")
endif()

set(_salmon_boost_components system filesystem timer chrono program_options)
if(SALMON_USE_SYSTEM_DEPS)
  find_package(Boost 1.59.0 QUIET COMPONENTS ${_salmon_boost_components})
endif()
if(Boost_FOUND)
  message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
  message(STATUS "Boost libraries: ${Boost_LIBRARIES}")
elseif(SALMON_FETCH_MISSING_DEPS)
  message(STATUS "Static Boost components not found; fetching pinned Boost release")
  set(_salmon_boost_prefix "${SALMON_DEPS_INSTALL_PREFIX}")
  set(_salmon_boost_url "https://archives.boost.io/release/1.84.0/source/boost_1_84_0.tar.bz2")
  set(_salmon_boost_libs "system,filesystem,timer,chrono,program_options,atomic")
  externalproject_add(libboost
    PREFIX ${CMAKE_BINARY_DIR}/_deps/libboost
    URL ${_salmon_boost_url}
    SOURCE_SUBDIR .
    INSTALL_DIR ${_salmon_boost_prefix}
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ./bootstrap.sh --prefix=<INSTALL_DIR> --with-libraries=${_salmon_boost_libs}
    BUILD_COMMAND ./b2 install --prefix=<INSTALL_DIR> --layout=system link=static threading=multi runtime-link=shared variant=release cxxstd=17
    INSTALL_COMMAND ""
  )
  set(FETCHED_BOOST TRUE)
  set(Boost_INCLUDE_DIRS "${_salmon_boost_prefix}/include")
  set(Boost_LIBRARIES
      "${_salmon_boost_prefix}/lib/libboost_system.a"
      "${_salmon_boost_prefix}/lib/libboost_filesystem.a"
      "${_salmon_boost_prefix}/lib/libboost_timer.a"
      "${_salmon_boost_prefix}/lib/libboost_chrono.a"
      "${_salmon_boost_prefix}/lib/libboost_program_options.a"
      "${_salmon_boost_prefix}/lib/libboost_atomic.a")
  message(STATUS "Using fetched static Boost from ${_salmon_boost_prefix}")
else()
  message(FATAL_ERROR "Boost (static=${SALMON_BOOST_USE_STATIC_LIBS}) is required. Install static Boost components or enable SALMON_FETCH_MISSING_DEPS.")
endif()

if(SALMON_USE_SYSTEM_DEPS)
  find_package(cereal "1.3.2" QUIET)
endif()
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

if(SALMON_ENABLE_TESTS)
  if(SALMON_USE_SYSTEM_DEPS)
    find_package(Catch2 3 QUIET CONFIG)
  endif()
  if(NOT Catch2_FOUND)
    if(SALMON_FETCH_MISSING_DEPS)
      message(STATUS "Catch2 not found; fetching pinned Catch2 release")
      set(CATCH_BUILD_TESTING OFF CACHE BOOL "" FORCE)
      set(CATCH_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
      set(CATCH_INSTALL_EXTRAS OFF CACHE BOOL "" FORCE)
      FetchContent_Declare(salmon_catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.12.0
        GIT_SHALLOW FALSE
      )
      FetchContent_MakeAvailable(salmon_catch2)
    else()
      message(FATAL_ERROR "Catch2 >= 3 is required for tests. Install Catch2 or enable SALMON_FETCH_MISSING_DEPS.")
    endif()
  endif()

  if(TARGET Catch2::Catch2WithMain)
    set(SALMON_CATCH2_TARGET Catch2::Catch2WithMain CACHE INTERNAL "" FORCE)
  else()
    message(FATAL_ERROR "Catch2::Catch2WithMain target was not found after dependency resolution.")
  endif()
endif()

set(FETCHED_TBB FALSE)
if(SALMON_USE_SYSTEM_DEPS)
  # Prefer package config mode to avoid brittle module-mode discovery.
  find_package(TBB 2021.4 QUIET CONFIG COMPONENTS tbb)
  if(NOT TBB_FOUND)
    find_package(TBB 2021.4 QUIET COMPONENTS tbb)
  endif()
endif()

if(TBB_FOUND AND TARGET TBB::tbb)
  if(TBB_VERSION VERSION_LESS 2021.4)
    message(FATAL_ERROR "Found TBB version ${TBB_VERSION}, but Salmon requires >= 2021.4.")
  endif()
  get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Found suitable TBB version: ${TBB_VERSION}")
elseif(SALMON_FETCH_MISSING_DEPS)
  message(STATUS "TBB >= 2021.4 not found; fetching pinned oneTBB release")
  set(TBB_TEST OFF CACHE BOOL "" FORCE)
  set(TBB_STRICT OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(salmon_tbb
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
    GIT_TAG v2022.3.0
    GIT_SHALLOW FALSE
  )
  FetchContent_MakeAvailable(salmon_tbb)
  if(TARGET tbb AND NOT TARGET TBB::tbb)
    add_library(TBB::tbb ALIAS tbb)
  endif()
  if(NOT TARGET TBB::tbb)
    message(FATAL_ERROR "Fetched oneTBB but TBB::tbb target is unavailable.")
  endif()
  set(FETCHED_TBB TRUE)
  get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Using fetched oneTBB target")
else()
  message(FATAL_ERROR "TBB >= 2021.4 is required. Install oneTBB or enable SALMON_FETCH_MISSING_DEPS.")
endif()


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

if(SALMON_USE_SYSTEM_DEPS)
  find_package(libgff 2.0.1 HINTS ${LIB_GFF_PATH} ${GFF_ROOT})
endif()
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
if(SALMON_USE_SYSTEM_DEPS)
  find_package(HTSlib QUIET)
endif()
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
    CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR> --disable-libcurl --disable-ref-cache CPPFLAGS=${_salmon_htslib_cppflags} LDFLAGS=${_salmon_htslib_ldflags} CC=${CMAKE_C_COMPILER}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
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
  if(SALMON_USE_SYSTEM_DEPS)
    find_package(Mimalloc QUIET)
  endif()
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

# Enforce a modern spdlog/fmt floor for runtime-format safety and consistent
# behavior across Salmon and embedded pufferfish headers.
unset(SALMON_SPDLOG_FMT_OK CACHE)
set(_SALMON_PREV_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES}")
set(CMAKE_REQUIRED_INCLUDES
    "${SALMON_PUFFERFISH_SOURCE_DIR}/include;${GAT_SOURCE_DIR}/include")
check_cxx_source_compiles(
  "
  #include <spdlog/version.h>
  #include <spdlog/fmt/fmt.h>
  #if !defined(SPDLOG_VERSION) || (SPDLOG_VERSION < 11700)
  #error spdlog too old
  #endif
  #if !defined(FMT_VERSION) || (FMT_VERSION < 120000)
  #error fmt too old
  #endif
  int main() { return 0; }
  "
  SALMON_SPDLOG_FMT_OK)
set(CMAKE_REQUIRED_INCLUDES "${_SALMON_PREV_REQUIRED_INCLUDES}")
unset(_SALMON_PREV_REQUIRED_INCLUDES)
if(NOT SALMON_SPDLOG_FMT_OK)
  message(FATAL_ERROR "Active spdlog/fmt headers do not satisfy required minimum versions (spdlog>=1.17, fmt>=12).")
endif()

# Keep phmap aligned across Salmon and pufferfish and enforce a modern floor.
unset(SALMON_PHMAP_OK CACHE)
set(_SALMON_PREV_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES}")
set(CMAKE_REQUIRED_INCLUDES
    "${SALMON_PUFFERFISH_SOURCE_DIR}/include;${GAT_SOURCE_DIR}/include")
check_cxx_source_compiles(
  "
  #include <parallel_hashmap/phmap_config.h>
  #if !defined(PHMAP_VERSION_MAJOR) || (PHMAP_VERSION_MAJOR < 2)
  #error phmap too old
  #endif
  int main() { return 0; }
  "
  SALMON_PHMAP_OK)
set(CMAKE_REQUIRED_INCLUDES "${_SALMON_PREV_REQUIRED_INCLUDES}")
unset(_SALMON_PREV_REQUIRED_INCLUDES)
if(NOT SALMON_PHMAP_OK)
  message(FATAL_ERROR "Active parallel_hashmap headers do not satisfy required minimum version (major>=2).")
endif()

# Enforce minimum nlohmann/json version for serialized metadata compatibility.
unset(SALMON_JSON_OK CACHE)
set(_SALMON_PREV_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES}")
set(CMAKE_REQUIRED_INCLUDES
    "${SALMON_PUFFERFISH_SOURCE_DIR}/include;${GAT_SOURCE_DIR}/include")
check_cxx_source_compiles(
  "
  #include <salmon/vendor/json.hpp>
  #if !defined(NLOHMANN_JSON_VERSION_MAJOR) || !defined(NLOHMANN_JSON_VERSION_MINOR)
  #error nlohmann/json version macros missing
  #endif
  #if (NLOHMANN_JSON_VERSION_MAJOR < 3) || ((NLOHMANN_JSON_VERSION_MAJOR == 3) && (NLOHMANN_JSON_VERSION_MINOR < 11))
  #error nlohmann/json too old
  #endif
  int main() { return 0; }
  "
  SALMON_JSON_OK)
set(CMAKE_REQUIRED_INCLUDES "${_SALMON_PREV_REQUIRED_INCLUDES}")
unset(_SALMON_PREV_REQUIRED_INCLUDES)
if(NOT SALMON_JSON_OK)
  message(FATAL_ERROR "Active nlohmann/json headers do not satisfy required minimum version (>=3.11).")
endif()

# Enforce Eigen floor for numerical stability/performance features used in
# Salmon and pufferfish.
unset(SALMON_EIGEN_OK CACHE)
set(_SALMON_PREV_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES}")
set(CMAKE_REQUIRED_INCLUDES
    "${GAT_SOURCE_DIR}/include/eigen3;${SALMON_PUFFERFISH_SOURCE_DIR}/include")
check_cxx_source_compiles(
  "
  #include <Eigen/Core>
  #if !defined(EIGEN_WORLD_VERSION) || !defined(EIGEN_MAJOR_VERSION) || !defined(EIGEN_MINOR_VERSION)
  #error Eigen version macros missing
  #endif
  #if (EIGEN_WORLD_VERSION < 3) || ((EIGEN_WORLD_VERSION == 3) && (EIGEN_MAJOR_VERSION < 4))
  #error Eigen too old
  #endif
  int main() { return 0; }
  "
  SALMON_EIGEN_OK)
set(CMAKE_REQUIRED_INCLUDES "${_SALMON_PREV_REQUIRED_INCLUDES}")
unset(_SALMON_PREV_REQUIRED_INCLUDES)
if(NOT SALMON_EIGEN_OK)
  message(FATAL_ERROR "Active Eigen headers do not satisfy required minimum version (>=3.4).")
endif()
