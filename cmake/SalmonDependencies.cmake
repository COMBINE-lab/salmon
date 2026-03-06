include_guard(GLOBAL)

include(ExternalProject)
include(FetchContent)

set(SALMON_PUFFERFISH_GIT_REPOSITORY
    "https://github.com/COMBINE-lab/pufferfish.git"
    CACHE STRING "Git repository used when fetching pufferfish")
set(SALMON_PUFFERFISH_GIT_TAG
    "1b15f5c6ac7cea73d30dea9664daaa0907044f34"
    CACHE STRING "Immutable git commit used when fetching pufferfish")
set(SALMON_PUFFERFISH_SOURCE_DIR
    ""
    CACHE PATH "Optional local pufferfish source checkout to use instead of fetching")

if(DEFINED CUSTOM_BOOST_PATH)
  set(CMAKE_INCLUDE_PATH ${CUSTOM_BOOST_PATH} ${CMAKE_INCLUDE_PATH})
  set(CMAKE_LIBRARY_PATH ${CUSTOM_BOOST_PATH}/lib ${CMAKE_LIBRARY_PATH})
endif()

if(CONDA_BUILD)
  set(Boost_USE_STATIC_LIBS OFF)
elseif(USE_SHARED_LIBS)
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
  if(NOT TARGET salmon_stage_zlibng)
    add_custom_target(salmon_stage_zlibng
      COMMAND ${CMAKE_COMMAND} --install ${salmon_zlibng_BINARY_DIR} --prefix ${CMAKE_CURRENT_SOURCE_DIR}/external/install
      COMMENT "Staging zlib-ng for external dependency consumers")
    add_dependencies(salmon_stage_zlibng ${SALMON_ZLIB_TARGET})
  endif()
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
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L http://tukaani.org/xz/xz-5.2.2.tar.gz -o xz-5.2.2.tar.gz &&
      ${SHASUM} 73df4d5d34f0468bd57d09f2d8af363e95ed6cc3a4a86129d2f2c366259902a2 xz-5.2.2.tar.gz &&
      tar -xzvf xz-5.2.2.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/xz-5.2.2
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/xz-5.2.2/configure --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS=${EXTRA_CMAKE_INCLUDE_FLAGS} CPPFLAGS=${EXTRA_CMAKE_INCLUDE_FLAGS} LDFLAGS=${EXTRA_CMAKE_LIBRARY_FLAGS}
    BUILD_COMMAND make ${QUIET_MAKE}
    INSTALL_COMMAND make ${QUIET_MAKE} install
  )
  set(LIBLZMA_LIBRARIES ${GAT_SOURCE_DIR}/external/install/lib/liblzma.a)
  set(LIBSTADEN_LDFLAGS "-L${GAT_SOURCE_DIR}/external/install/lib")
  set(LIBSTADEN_CFLAGS "-I${GAT_SOURCE_DIR}/external/install/include")
  set(FETCHED_LIBLZMA TRUE)
else()
  message(STATUS "Found liblzma library: ${LIBLZMA_LIBRARIES}")
endif()

find_package(BZip2)
if(NOT BZIP2_FOUND)
  message(STATUS "Will attempt to fetch and build libbz2")
  externalproject_add(libbz2
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://sourceware.org/pub/bzip2/bzip2-1.0.6.tar.gz -o bzip2-1.0.6.tar.gz &&
      ${SHASUM} a2848f34fcd5d6cf47def00461fcb528a0484d8edef8208d6d2e2909dc61d9cd bzip2-1.0.6.tar.gz &&
      tar -xzvf  bzip2-1.0.6.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bzip2-1.0.6
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make ${QUIET_MAKE} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    INSTALL_COMMAND make ${QUIET_MAKE} install PREFIX=<INSTALL_DIR>
  )
  set(BZIP2_LIBRARIES ${GAT_SOURCE_DIR}/external/install/lib/libbz2.a)
  set(LIBSTADEN_LDFLAGS "-L${GAT_SOURCE_DIR}/external/install/lib -I${GAT_SOURCE_DIR}/external/install/include")
  set(LIBSTADEN_CFLAGS "-I${GAT_SOURCE_DIR}/external/install/include")
  set(FETCHED_LIBBZ2 TRUE)
else()
  message(STATUS "Found libbz2 library: ${BZIP2_LIBRARIES}")
endif()

set(Boost_ADDITIONAL_VERSIONS
    "1.59.0" "1.60.0" "1.61.0" "1.62.0" "1.63.0" "1.64.0" "1.65.0" "1.66.0"
    "1.67.0" "1.68.0" "1.69.0" "1.70.0" "1.71.0" "1.72.0" "1.73.0" "1.74.0"
    "1.75.0" "1.76.0" "1.77.0" "1.78.0" "1.79.0" "1.80.0" "1.81.0" "1.82.0"
    "1.83.0" "1.84.0")
if(NOT BOOST_RECONFIGURE)
  find_package(Boost 1.59.0 COMPONENTS iostreams system filesystem timer chrono program_options)
  message(STATUS "BOOST_INCLUDEDIR = ${BOOST_INCLUDEDIR}")
  message(STATUS "BOOST_LIBRARYDIR = ${BOOST_LIBRARYDIR}")
  message(STATUS "Boost_FOUND = ${Boost_FOUND}")
endif()

if(BOOST_RECONFIGURE)
  message(STATUS "Executing Boost Reconfiguration")
  unset(Boost_FOUND CACHE)
  unset(Boost_INCLUDE_DIR CACHE)
  unset(Boost_INCLUDE_DIRS CACHE)
  unset(Boost_LIBRARY_DIRS CACHE)
  unset(Boost_LIBRARIES CACHE)
  unset(BOOST_ROOT CACHE)
  unset(CMAKE_PREFIX_PATH CACHE)
  unset(Boost::diagnostic_definitions CACHE)
  unset(Boost::disable_autolinking CACHE)
  unset(Boost::dynamic_linking CACHE)
  set(BOOST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
  set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
  set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
  set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
  find_package(Boost 1.59.0 COMPONENTS iostreams system filesystem timer chrono program_options locale REQUIRED)
  set(FETCH_BOOST FALSE)
endif()

if((NOT Boost_FOUND) AND (NOT FETCH_BOOST))
  message(FATAL_ERROR
    "Salmon cannot be compiled without Boost.\n"
    "It is recommended to visit http://www.boost.org/ and install Boost according to those instructions.\n"
    "This build system can also download and install a local version of boost for you (this takes a lot of time).\n"
    "To fetch and build boost locally, call cmake with -DFETCH_BOOST=TRUE")
elseif(FETCH_BOOST)
  if(NOT DEFINED BOOST_BUILD_THREADS)
    set(BOOST_BUILD_THREADS 2)
  endif()

  set(BOOST_LIB_SUBSET --with-iostreams --with-atomic --with-chrono --with-container --with-date_time --with-exception
      --with-filesystem --with-graph --with-graph_parallel --with-math
      --with-program_options --with-system --with-locale --with-timer)
  set(BOOST_WILL_RECONFIGURE TRUE)
  set(FETCH_BOOST FALSE)
  set(BOOST_FETCHED_VERSION "1_72_0")
  message(STATUS "Build system will fetch and build Boost")
  externalproject_add(libboost
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.gz/download -o boost_1_72_0.tar.gz &&
      ${SHASUM} c66e88d5786f2ca4dbebb14e06b566fb642a1a6947ad8cc9091f9f445134143f boost_${BOOST_FETCHED_VERSION}.tar.gz &&
      tar xzf boost_${BOOST_FETCHED_VERSION}.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_${BOOST_FETCHED_VERSION}
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_${BOOST_FETCHED_VERSION}/bootstrap.sh ${BOOST_CONFIGURE_TOOLSET} ${BOOST_BUILD_LIBS} --prefix=<INSTALL_DIR>
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_${BOOST_FETCHED_VERSION}/tools/build/src/user-config.jam
      PRE_BUILD
      COMMAND echo "using gcc : ${CC_VERSION} : ${CMAKE_CXX_COMPILER} ;"
    )
    BUILD_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_${BOOST_FETCHED_VERSION}/b2 -d0 -j${BOOST_BUILD_THREADS} ${BOOST_LIB_SUBSET} toolset=${BOOST_TOOLSET} ${BOOST_EXTRA_FLAGS} cxxflags=${BOOST_CXX_FLAGS} link=static install
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
  )

  externalproject_add_step(libboost makedir
    COMMAND mkdir -p <SOURCE_DIR>/build
    COMMENT "Make build directory"
    DEPENDEES download
    DEPENDERS configure)

  set(RECONFIG_FLAGS ${RECONFIG_FLAGS} -DBOOST_WILL_RECONFIGURE=FALSE -DBOOST_RECONFIGURE=TRUE -DFETCH_BOOST=FALSE)
  externalproject_add_step(libboost reconfigure
    COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
    DEPENDEES install
  )
  set(FETCHED_BOOST TRUE)
endif()

if(BOOST_WILL_RECONFIGURE)
  message(STATUS "Setting temporary Boost paths")
  set(Boost_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
  set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
  set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
  set(Boost_FOUND TRUE)
endif()

message(STATUS "BOOST ROOT = ${BOOST_ROOT}")
message(STATUS "BOOST INCLUDE DIR = ${Boost_INCLUDE_DIR}")
message(STATUS "BOOST INCLUDE DIRS = ${Boost_INCLUDE_DIRS}")
message(STATUS "BOOST LIB DIR = ${Boost_LIBRARY_DIRS}")
message(STATUS "BOOST LIBRARIES = ${Boost_LIBRARIES}")

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

find_package(TBB 2021.4
  HINTS ${TBB_ROOT_SEARCH}
  COMPONENTS tbb)

if(TBB_FOUND)
  if(TBB_VERSION VERSION_GREATER_EQUAL 2021.4)
    message(STATUS "FOUND SUITABLE TBB VERSION : ${TBB_VERSION}")
    set(TBB_TARGET_EXISTED TRUE)
    get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
  else()
    set(TBB_TARGET_EXISTED FALSE)
  endif()
else()
  set(TBB_TARGET_EXISTED FALSE)
endif()

if(NOT TBB_TARGET_EXISTED)
  set(TBB_WILL_RECONFIGURE TRUE)
  if(CLANG)
    set(TBB_COMPILER "clang")
  else()
    set(TBB_COMPILER "gcc")
  endif()

  message(STATUS "Build system will fetch and build Intel Threading Building Blocks")
  set(TBB_SOURCE_DIR ${GAT_SOURCE_DIR}/external/oneTBB-2021.11.0)
  set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
  set(TBB_CXXFLAGS "${TBB_CXXFLAGS} ${CXXSTDFLAG} ${SCHAR_FLAG}")

  ExternalProject_Add(libtbb
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2021.11.0.tar.gz -o v2021.11.0.tar.gz &&
      ${SHASUM} 782ce0cab62df9ea125cdea253a50534862b563f1d85d4cda7ad4e77550ac363 v2021.11.0.tar.gz &&
      tar -xzvf v2021.11.0.tar.gz
    SOURCE_DIR ${TBB_SOURCE_DIR}
    INSTALL_DIR ${TBB_INSTALL_DIR}
    PATCH_COMMAND "${TBB_PATCH_STEP}"
    CMAKE_ARGS -DCMAKE_CXX_FLAGS=${TBB_CXXFLAGS} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DTBB_TEST=OFF -DTBB_EXAMPLES=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    BUILD_IN_SOURCE TRUE
  )

  set(RECONFIG_FLAGS ${RECONFIG_FLAGS} -DTBB_WILL_RECONFIGURE=FALSE -DTBB_RECONFIGURE=TRUE)
  ExternalProject_Add_Step(libtbb reconfigure
    COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
    DEPENDEES install
  )
  set(FETCHED_TBB TRUE)
  set(TBB_ROOT_SEARCH ${CMAKE_SOURCE_DIR}/external/install)

  if(FETCHED_BOOST)
    add_dependencies(libtbb libboost)
  endif()
endif()

if(TBB_WILL_RECONFIGURE)
  set(TBB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
  set(TBB_INCLUDE_DIRS ${TBB_INSTALL_DIR}/include)
  set(TBB_INCLUDE_DIR ${TBB_INSTALL_DIR}/include)
  set(TBB_LIBRARY_DIRS ${TBB_INSTALL_DIR}/lib)
  set(TBB_LIBRARY ${TBB_INSTALL_DIR}/lib)
  set(TBB_LIB_DIR ${TBB_INSTALL_DIR}/lib)
  set(TBB_LIBRARIES
      ${TBB_INSTALL_DIR}/lib/libtbb.${SHARED_LIB_EXTENSION})
  message(STATUS "TBB_INCLUDE_DIRS = ${TBB_INCLUDE_DIRS}")
  message(STATUS "TBB_LIBRARY_DIRS = ${TBB_LIBRARY_DIRS}")
endif()

if(TBB_RECONFIGURE)
  unset(TBB_FOUND CACHE)
  unset(TBB_INSTALL_DIR CACHE)
  unset(CMAKE_PREFIX_PATH CACHE)
  unset(TBB_INCLUDE_DIRS CACHE)
  unset(TBB_INCLUDE_DIR CACHE)
  unset(TBB_LIBRARY_DIRS CACHE)
  unset(TBB_LIBRARY CACHE)
  unset(TBB_LIBRARIES CACHE)
  set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
  set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
  set(TBB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
  set(TBB_INCLUDE_DIRS ${TBB_INSTALL_DIR}/include)
  set(TBB_INCLUDE_DIR ${TBB_INSTALL_DIR}/include)
  set(TBB_LIBRARY_DIRS ${TBB_INSTALL_DIR}/lib)
  set(TBB_LIBRARY ${TBB_INSTALL_DIR}/lib)
  set(TBB_LIB_DIR ${TBB_INSTALL_DIR}/lib)
  message(STATUS "TBB_INSTALL_DIR = ${TBB_INSTALL_DIR}")
  find_package(TBB 2021.4
    HINTS ${TBB_ROOT_SEARCH}
    COMPONENTS tbb)
  message(STATUS "[in TBB_RECONFIGURE] TBB_LIBRARIES = ${TBB_LIBRARIES}")
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
set(FETCHED_PUFFERFISH TRUE CACHE BOOL "Has pufferfish been fetched?" FORCE)
set(SALMON_PUFFERFISH_SOURCE_DIR "${salmon_pufferfish_SOURCE_DIR}" CACHE INTERNAL "" FORCE)
set(SALMON_PUFFERFISH_BINARY_DIR "${salmon_pufferfish_BINARY_DIR}" CACHE INTERNAL "" FORCE)

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
  externalproject_add(libhtslib
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 -o htslib-1.22.tar.bz2 &&
      ${SHASUM} 6250c1df297db477516e60ac8df45ed75a652d1f25b0f37f12f5b17269eafde9 htslib-1.22.tar.bz2 &&
      tar -xjf htslib-1.22.tar.bz2
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/htslib-1.22
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR> --disable-libcurl CPPFLAGS=-I${CMAKE_CURRENT_SOURCE_DIR}/external/install/include LDFLAGS=-L${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib CC=${CMAKE_C_COMPILER}
    BUILD_COMMAND make ${QUIET_MAKE}
    INSTALL_COMMAND make ${QUIET_MAKE} install
  )
  if(FETCHED_ZLIBNG)
    add_dependencies(libhtslib salmon_stage_zlibng)
  endif()
  if(FETCHED_LIBBZ2)
    add_dependencies(libhtslib libbz2)
  endif()
  if(FETCHED_LIBLZMA)
    add_dependencies(libhtslib liblzma)
  endif()
  set(FETCHED_HTSLIB TRUE)
  set(HTSLIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
  if(NOT TARGET HTSlib::HTSlib)
    add_library(HTSlib::HTSlib UNKNOWN IMPORTED)
    set_target_properties(HTSlib::HTSlib PROPERTIES
      IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib/libhts.a"
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
