cmake_minimum_required(VERSION 3.15)

if(DEFINED ENV{CC})
  set(CC $ENV{CC})
else()
  set(CC gcc)
endif()
message("CC: ${CC}")

set(CC_VERSION "")
if(${CC} MATCHES ^gcc-)
  string(REGEX REPLACE "gcc-" "" CC_VERSION ${CC})
endif()
message("CC version: ${CC_VERSION}")

enable_testing()

project(Salmon)

# detect host architecture 
if(NOT DEFINED USE_ARM)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*|ARM64.*)")
    message("Detected 64-bit ARM host. Setting USE_ARM to true.")
    set(USE_ARM TRUE)
    # set char to be signed
    add_compile_options(-fsigned-char)
    set(SCHAR_FLAG "-fsigned-char")
  else()
    message("Detected non-ARM host. Setting USE_ARM to false.")
    set(USE_ARM FALSE)
    set(SCHAR_FLAG "")
  endif()
endif()

option(USE_SHARED_LIBS "Use shared instead of static libraries" OFF)

# auto-populate version:
# from https://stackoverflow.com/questions/47066115/cmake-get-version-from-multi-line-text-file
file(READ "current_version.txt" ver)

string(REGEX MATCH "VERSION_MAJOR ([0-9]*)" _ ${ver})
set(ver_major ${CMAKE_MATCH_1})

string(REGEX MATCH "VERSION_MINOR ([0-9]*)" _ ${ver})
set(ver_minor ${CMAKE_MATCH_1})

string(REGEX MATCH "VERSION_PATCH ([0-9]*)" _ ${ver})
set(ver_patch ${CMAKE_MATCH_1})

set(CPACK_PACKAGE_VERSION_MAJOR ${ver_major})
set(CPACK_PACKAGE_VERSION_MINOR ${ver_minor})
set(CPACK_PACKAGE_VERSION_PATCH ${ver_patch})

set(CPACK_PACKAGE_VERSION "${ver_major}.${ver_minor}.${ver_patch}")
message("version: ${CPACK_PACKAGE_VERSION}")

set(PROJECT_VERSION ${CPACK_PACKAGE_VERSION})
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "University of Maryland")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Salmon - Wicked-fast RNA-seq isoform quantification using selective alignment")
set(CPACK_PACKAGE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-Source")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  #set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
  #  "Debug" "Release")
endif()

message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

## Set the standard required compile flags
set(REMOVE_WARNING_FLAGS "-Wno-unused-function;-Wno-unused-local-typedefs")
set(TGT_COMPILE_FLAGS "${SCHAR_FLAG};-ftree-vectorize;-funroll-loops;-fPIC;-fomit-frame-pointer;-O3;-DNDEBUG;-DSTX_NO_STD_STRING_VIEW;-D__STDC_FORMAT_MACROS")
set(TGT_WARN_FLAGS "-Wall;-Wno-unknown-pragmas;-Wno-reorder;-Wno-unused-variable;-Wreturn-type;-Werror=return-type;${REMOVE_WARNING_FLAGS}")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

###
# Sanitizers BEGIN
###
if (ASAN_BUILD)
  list(APPEND TGT_COMPILE_FLAGS "-fsanitize=address")
#list(APPEND TGT_COMPILE_FLAGS "-fsanitize=undefined")
#set(CMAKE_LINK_FLAGS "-fsanitize=address")
#list(APPEND CMAKE_LINK_FLAGS "-fsanitize=undefined")
  set(ASAN_LIB "asan")
else()
  set(ASAN_LIB "")
endif()
###
# Sanitizers END
###

if(APPLE)
set(WARNING_IGNORE_FLAGS "-Wno-deprecated-register")
list(APPEND TGT_WARN_FLAGS -Wno-deprecated-register)
else()
set(WARNING_IGNORE_FLAGS "")
endif()

## Prefer static to dynamic libraries
if(NOT USE_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  set(MALLOC_STATIC_BUILD_FLAG "--enable-static")
endif()

include(CheckIPOSupported)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(ICU_LIBS "")
set(ICU_INC_DIRS "")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CXXSTDFLAG "-std=c++14")
set(GCCVERSION "5.2")

if(CONDA_BUILD)
  message("Building with CONDA_BUILD flag")
  if(APPLE)
    # Do we require all these components?  Any others?
    find_package(ICU COMPONENTS data i18n io uc REQUIRED)
    if(ICU_FOUND)
      message(STATUS "ICU_INCLUDE_DIRS = ${ICU_INCLUDE_DIRS}")
      message(STATUS "ICU_LIBRARIES = ${ICU_LIBRARIES}")
    endif()
    set(ICU_INC_DIRS ${ICU_INCLUDE_DIRS})
    set(ICU_LIBS ${ICU_LIBRARIES})
  endif()
endif()

set(BOOST_CXX_FLAGS "${WARNING_IGNORE_FLAGS} ${CXXSTDFLAG}")
if(FETCH_BOOST)
  set(BOOST_CXX_FLAGS "${BOOST_CXX_FLAGS} -I${CMAKE_CURRENT_SOURCE_DIR}/external/install/include -L${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib")
endif()

##
# OSX is strange (some might say, stupid in this regard).  Deal with it's quirkines here.
##
if(APPLE)
  # To allow ourselves to build a dynamic library, we have to tell the compiler
  # that, yes, the symbols will be around at runtime.
  list(APPEND TGT_COMPILE_FLAGS "-undefined dynamic_lookup;-Wno-unused-command-line-argument")
  # set(LIBSALMON_LINKER_FLAGS "-all_load")
  # In order to "think different", we also have to use non-standard suffixes
  # for our shared libraries
  set(SHARED_LIB_EXTENSION "dylib")
else()
  # We're in sane linux world
  set(SHARED_LIB_EXTENSION "so")
  set(LIBSALMON_LINKER_FLAGS "")
endif()

set( BOOST_EXTRA_FLAGS "--layout=tagged" )
## this get's set differently below if we
## are on clang & apple
set(NON_APPLECLANG_LIBS gomp)

if(UNIX AND NOT APPLE)
  set(LIBRT rt)
endif()

set(PTHREAD_LIB)

##
# Let us check the sha sum of our pacakges if we have the right tools
##
set(SHASUM ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_shasum.sh)

##
# Compiler-specific C++11/14 activation.
# http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake
##
##
# First take care of what to do if we have gcc
##
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    # If we're on OSX
    if(APPLE AND NOT (GCC_VERSION VERSION_GREATER ${GCCVERSION} OR GCC_VERSION VERSION_EQUAL ${GCCVERSION}))
      message(FATAL_ERROR "When building under OSX, ${PROJECT_NAME} requires "
                            "either clang or g++ >= ${GCCVERSION}")
    elseif(NOT (GCC_VERSION VERSION_GREATER ${GCCVERSION} OR GCC_VERSION VERSION_EQUAL ${GCCVERSION}))
      message(FATAL_ERROR "${PROJECT_NAME} requires g++ ${GCCVERSION} or greater.")
    endif()
    
    if(GCC_VERSION VERSION_GREATER_EQUAL "7.1")
      list(APPEND TGT_WARN_FLAGS "-Wno-int-in-bool-context")
    endif()
    
    if(GCC_VERSION VERSION_GREATER_EQUAL "9.1")
      list(APPEND TGT_WARN_FLAGS "-Wno-deprecated-copy")
    endif()

    set(GCC TRUE)

    # Put complete static linking on hold for the time-being
    # If we're not on OSX, make an attempt to compile everything statically
    #if(NOT APPLE)
    #set(CMAKE_EXE_LINK_FLAGS "-static")
    set(PTHREAD_LIB "pthread")
    #endif()

    # If we're on Linux (i.e. not OSX) and we're using
    # gcc, then set the -static-libstdc++ flag
    if(NOT APPLE)
        list(APPEND TGT_COMPILE_FLAGS -static-libstdc++)
    endif()

    set(WARNING_IGNORE_FLAGS "${WARNING_IGNORE_FLAGS} -Wno-unused-local-typedefs")
    set(BOOST_TOOLSET "${CC}")
    set(BOOST_CONFIGURE_TOOLSET "--with-toolset=gcc")
    set(BCXX_FLAGS "${CXXSTDFLAG} ${SCHAR_FLAG}")
    set(BOOST_EXTRA_FLAGS toolset=${CC} cxxflags=${BCXX_FLAGS})
# Tentatively, we support clang now
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CLANG TRUE)
    # If we have libc++, then try and use it
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-stdlib=libc++ HAVE_LIBCPP)
    if(HAVE_LIBCPP)
      message("It appears that you're compiling with clang and that libc++ is available, so I'll use that")
      list(APPEND TGT_COMPILE_FLAGS -stdlib=libc++)
      set(BOOST_TOOLSET "clang")
      set(BOOST_CONFIGURE_TOOLSET "--with-toolset=clang")
      set(BCXX_FLAGS "-stdlib=libc++ -DBOOST_HAS_INT128 ${SCHAR_FLAG}")
      set(BOOST_EXTRA_FLAGS toolset=clang cxxflags=${BCXX_FLAGS} linkflags="-stdlib=libc++")
    # Otherwise, use libstdc++ (and make it static)
    else()
      list(APPEND TGT_COMPILE_FLAGS -static-libstdc++)
    endif()
    # There's currently a bug with clang-3.4 & Boost 1.55 -- this hack fixes it
    # but we should do something better (does this break things if CPU doesn't
    # have 128-bit support)?
    list(APPEND TGT_COMPILE_FLAGS -DBOOST_HAS_INT128)

    if(APPLE)
      set(NON_APPLECLANG_LIBS "")
    else()
      set(PTHREAD_LIB "pthread")
    endif()
else()
  message(FATAL_ERROR "Your C++ compiler does not support C++14.")
endif()

if(DO_QUIET_MAKE)
  set(QUIET_MAKE "--silent")
else()
  set(QUIET_MAKE "")
endif()

## TODO: Figure out how to detect this automatically
# If the "assembler" is too old, tell TBB not to compile
# with -mrtm
if(NO_RTM)
  set(TBB_CXXFLAGS "-mno-rtm")
endif()

include(ExternalProject)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("Making Debug build")
elseif(CMAKE_BUILD_TYPE MATCHES Release)
  message("Making Release build")
else()
  message("Making Default build type")
endif()

##
# Record this top-level path
##
set(GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Have CMake tell us what it's doing
# set(CMAKE_VERBOSE_MAKEFILE true)

###
# check if numeric_limits<__int128_t> is defined
###
try_compile(HAVE_INT128_NUMERIC_LIMITS ${CMAKE_BINARY_DIR} 
            SOURCES ${GAT_SOURCE_DIR}/tests/compile_tests/int128_numeric_limits.cpp
            CXX_STANDARD 14
            CXX_STANDARD_REQUIRED ON 
            )
if(HAVE_INT128_NUMERIC_LIMITS)
  message("setting -DHAVE_NUMERIC_LIMITS128")
  list(APPEND TGT_COMPILE_FLAGS "-DHAVE_NUMERIC_LIMITS128")
else()
  message("not setting -DHAVE_NUMERIC_LIMITS128")
endif()

###
#
#  Grab pufferfish source --- DURING CONFIGURE TIME!
#
####
if(NOT FETCHED_PUFFERFISH)
  exec_program(${CMAKE_CURRENT_SOURCE_DIR}/scripts/fetchPufferfish.sh RETURN_VALUE FETCH_PF_SCRIPT_RET)
  message(STATUS "fetch PUFFERFISH exit code ${FETCH_PF_SCRIPT_RET}")
  if(NOT (FETCH_PF_SCRIPT_RET EQUAL 0))
    message(FATAL_ERROR "Could not fetch pufferfish source [fetchPufferfish.sh returned exit code ${FETCH_PF_SCRIPT_RET}].")
  endif()
  set(FETCHED_PUFFERFISH TRUE CACHE BOOL "Has pufferfish been fetched?" FORCE)
endif()

##
# Super-secret override
##
if( DEFINED CUSTOM_BOOST_PATH )
  set(CMAKE_INCLUDE_PATH ${CUSTOM_BOOST_PATH} ${CMAKE_INCLUDE_PATH})
  set(CMAKE_LIBRARY_PATH ${CUSTOM_BOOST_PATH}/lib ${CMAKE_LIBRARY_PATH})
endif()

##
# We want static, multithreaded boost libraries
##
if(CONDA_BUILD)
  set(Boost_USE_STATIC_LIBS OFF)
elseif(USE_SHARED_LIBS)   # CI failed when using an OR statement above...
  set(Boost_USE_STATIC_LIBS OFF)
else()
  set(Boost_USE_STATIC_LIBS ON)
endif()

set(Boost_USE_MULTITHREADED ON)
#set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_USE_DEBUG_RUNTIME OFF)

find_package(ZLIB)
if(NOT ZLIB_FOUND)
  message(FATAL_ERROR "zlib must be installed before configuration & building can proceed")
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
  message("Will attempt to fetch and build liblzma")
  message("=======================================")
externalproject_add(liblzma
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
  ##
  DOWNLOAD_COMMAND curl -k -L http://tukaani.org/xz/xz-5.2.2.tar.gz -o xz-5.2.2.tar.gz &&
    ${SHASUM} 73df4d5d34f0468bd57d09f2d8af363e95ed6cc3a4a86129d2f2c366259902a2 xz-5.2.2.tar.gz &&
    tar -xzvf xz-5.2.2.tar.gz
  #URL http://tukaani.org/xz/xz-5.2.2.tar.gz
  #URL_HASH SHA1=14663612422ab61386673be78fbb2556f50a1f08
  ##
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/xz-5.2.2
  INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
  BUILD_IN_SOURCE TRUE
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/xz-5.2.2/configure --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS=${EXTRA_CMAKE_INCLUDE_FLAGS} CPPFLAGS=${EXTRA_CMAKE_INCLUDE_FLAGS} LDFLAGS=${EXTRA_CMAKE_LIBRARY_FLAGS}
  BUILD_COMMAND make ${QUIET_MAKE}
  INSTALL_COMMAND make ${QUIET_MAKE} install
)

# Tell cmake that the external project generated a library so we can
# add dependencies here instead of later
set(LIBLZMA_LIBRARIES ${GAT_SOURCE_DIR}/external/install/lib/liblzma.a)
set(LIBSTADEN_LDFLAGS "-L${GAT_SOURCE_DIR}/external/install/lib")
set(LIBSTADEN_CFLAGS "-I${GAT_SOURCE_DIR}/external/install/include")
set(FETCHED_LIBLZMA TRUE)
else()
  message("Found liblzma library: ${LIBLZMA_LIBRARIES}")
  message("===========================================")
endif()

find_package(BZip2)
if(NOT BZIP2_FOUND)
  message("Will attempt to fetch and build libbz2")
  message("=======================================")
externalproject_add(libbz2
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
  DOWNLOAD_COMMAND curl -k -L https://sourceware.org/pub/bzip2/bzip2-1.0.6.tar.gz -o bzip2-1.0.6.tar.gz &&
    ${SHASUM} a2848f34fcd5d6cf47def00461fcb528a0484d8edef8208d6d2e2909dc61d9cd bzip2-1.0.6.tar.gz &&
    tar -xzvf  bzip2-1.0.6.tar.gz
  #URL http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
  #URL_HASH SHA1=3f89f861209ce81a6bab1fd1998c0ef311712002
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bzip2-1.0.6
  INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
  BUILD_IN_SOURCE TRUE
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make ${QUIET_MAKE} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
  INSTALL_COMMAND make ${QUIET_MAKE} install PREFIX=<INSTALL_DIR>
)
# Tell cmake that the external project generated a library so we can
# add dependencies here instead of later
set(BZIP2_LIBRARIES ${GAT_SOURCE_DIR}/external/install/lib/libbz2.a)
set(LIBSTADEN_LDFLAGS "-L${GAT_SOURCE_DIR}/external/install/lib -I${GAT_SOURCE_DIR}/external/install/include")
set(LIBSTADEN_CFLAGS "-I${GAT_SOURCE_DIR}/external/install/include")
set(FETCHED_LIBBZ2 TRUE)
else()
  message("Found libbz2 library: ${BZIP2_LIBRARIES}")
  message("===========================================")
endif()

##
# Set the latest version and look for what we need
##
set(Boost_ADDITIONAL_VERSIONS "1.59.0" "1.60.0" "1.61.0" "1.62.0" "1.63.0" "1.64.0" "1.65.0" "1.66.0" "1.67.0" "1.68.0" "1.69.0" "1.70.0" "1.71.0" "1.72.0" "1.73.0" "1.74.0"  "1.75.0" "1.76.0" "1.77.0" "1.78.0")
if (NOT BOOST_RECONFIGURE)
find_package(Boost 1.59.0 COMPONENTS iostreams system filesystem timer chrono program_options)
message("BOOST_INCLUDEDIR = ${BOOST_INCLUDEDIR}")
message("BOOST_LIBRARYDIR = ${BOOST_LIBRARYDIR}")
message("Boost_FOUND = ${Boost_FOUND}")
endif()

include(ExternalProject)

##
#  If we had to fetch Boost, the reconfigure step will re-run cmake.  The second configuration
#  pass is executed with the BOOST_RECONFIGURE flag set. This should allow our newly
#  installed Boost to be found by CMake.
##
if(BOOST_RECONFIGURE)
  message("Executing Boost Reconfiguration")
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



##
# Either inform the user of how to obtain Boost, or, if they passed in the FETCH_BOOST
# option, go and grab it for them.
##
if((NOT Boost_FOUND) AND (NOT FETCH_BOOST))
  message(FATAL_ERROR
    "Salmon cannot be compiled without Boost.\n"
    "It is recommended to visit http://www.boost.org/ and install Boost according to those instructions.\n"
    "This build system can also download and install a local version of boost for you (this takes a lot of time).\n"
    "To fetch and build boost locally, call cmake with -DFETCH_BOOST=TRUE"
  )
elseif(FETCH_BOOST)
  if(NOT DEFINED BOOST_BUILD_THREADS)
    set(BOOST_BUILD_THREADS 2)
  endif()

  ## Let the rest of the build process know we're going to be fetching boost
  set(BOOST_LIB_SUBSET --with-iostreams --with-atomic --with-chrono --with-container --with-date_time --with-exception
    --with-filesystem --with-graph --with-graph_parallel --with-math
    --with-program_options --with-system  --with-locale
    --with-timer)
  set(BOOST_WILL_RECONFIGURE TRUE)
  set(FETCH_BOOST FALSE)
  set(BOOST_FETCHED_VERSION "1_72_0")
  message("Build system will fetch and build Boost")
  message("==================================================================")
  externalproject_add(libboost
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.gz/download -o boost_1_72_0.tar.gz &&
      #${SHASUM} 96b34f7468f26a141f6020efb813f1a2f3dfb9797ecf76a7d7cbd843cc95f5bd boost_1_71_0.tar.gz &&
      ${SHASUM} c66e88d5786f2ca4dbebb14e06b566fb642a1a6947ad8cc9091f9f445134143f boost_${BOOST_FETCHED_VERSION}.tar.gz &&
      tar xzf boost_${BOOST_FETCHED_VERSION}.tar.gz 
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_${BOOST_FETCHED_VERSION}
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    #PATCH_COMMAND patch -p2 < ${CMAKE_CURRENT_SOURCE_DIR}/external/boost156.patch
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

    ##
    # After we've installed boost,
    ##
    set(RECONFIG_FLAGS ${RECONFIG_FLAGS} -DBOOST_WILL_RECONFIGURE=FALSE -DBOOST_RECONFIGURE=TRUE -DFETCH_BOOST=FALSE)
    externalproject_add_step(libboost reconfigure
      COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
      DEPENDEES install
    )
    set(FETCHED_BOOST TRUE)
endif()

##
# If we're fetching boost and we need to have dummy paths for these variables
# so that CMake won't complain
##
if(BOOST_WILL_RECONFIGURE)
  message("Setting Temporary Boost paths")
  set(Boost_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
  set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
  set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
  set(Boost_FOUND TRUE)
endif()

message("BOOST ROOT = ${BOOST_ROOT}")
message("BOOST INCLUDE DIR = ${Boost_INCLUDE_DIR}")
message("BOOST INCLUDE DIRS = ${Boost_INCLUDE_DIRS}")
message("BOOST LIB DIR = ${Boost_LIBRARY_DIRS}")
message("BOOST LIBRARIES = ${Boost_LIBRARIES}")

set(EXTERNAL_LIBRARY_PATH $CMAKE_CURRENT_SOURCE_DIR/lib)

#find_package(libdivsufsort)
#if(NOT LIBDIVSUFSORT_FOUND)
#  message("Build system will build libdivsufsort")
#  message("==================================================================")
#  include(ExternalProject)
#  externalproject_add(libdivsufsort
#    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
#    URL ${CMAKE_CURRENT_SOURCE_DIR}/external/libdivsufsort.zip
#    # Note: This zip comes from the fetched rapmap.zip, whose SHA we check
#    # so we souldn't need to check this one separately.
#    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libdivsufsort-master
#    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
#    #UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
#    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libdivsufsort-master/build
#    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DBUILD_DIVSUFSORT64=TRUE -DUSE_OPENMP=TRUE -DBUILD_SHARED_LIBS=FALSE
#  )
#  externalproject_add_step(libdivsufsort makedir
#    COMMAND mkdir -p <SOURCE_DIR>/build
#    COMMENT "Make build directory"
#    DEPENDEES download
#    DEPENDERS configure)
#
#  set(SUFFARRAY_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
#  set(FETCHED_LIBDIVSUFSORT TRUE)
#else()
#  message("SUFFARRAY_LIB = ${SUFFARRAY_LIBRARY}")
#  set(SUFFARRAY_LIB ${SUFFARRAY_LIBRARY})
#  message("SUFFARRAY_LIB64 = ${SUFFARRAY_LIBRARY64}")
#  set(SUFFARRAY_LIB64 ${SUFFARRAY_LIBRARY64})
#  set(SUFFARRAY_INCLUDE_DIRS ${SUFFARRAY_INCLUDE_DIR})
#endif()


find_package(Cereal)
if (NOT CEREAL_FOUND)
  message("Build system will fetch and build the Cereal serialization library")
  message("==================================================================")
  include(ExternalProject)
  externalproject_add(libcereal
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/USCiLab/cereal/archive/v1.3.0.tar.gz -o cereal-v1.3.0.tar.gz &&
      ${SHASUM} 329ea3e3130b026c03a4acc50e168e7daff4e6e661bc6a7dfec0d77b570851d5 cereal-v1.3.0.tar.gz &&
      tar -xzvf cereal-v1.3.0.tar.gz

    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cereal-1.3.0
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    #UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cereal-1.3.0/build
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND sh -c "mkdir -p <INSTALL_DIR>/include && cp -r <SOURCE_DIR>/include/cereal <INSTALL_DIR>/include"
  )
  externalproject_add_step(libcereal makedir
    COMMAND mkdir -p <SOURCE_DIR>/build
    COMMENT "Make build directory"
    DEPENDEES download
    DEPENDERS configure)

  set(FETCHED_CEREAL TRUE)
endif()

## Try and find TBB first
find_package(TBB 2021.4
  HINTS ${TBB_ROOT_SEARCH}
  COMPONENTS tbb tbbmalloc tbbmalloc_proxy)

if (${TBB_FOUND})
  if (${TBB_VERSION} VERSION_GREATER_EQUAL 2021.4)
  message("FOUND SUITABLE TBB VERSION : ${TBB_VERSION}")
  set(TBB_TARGET_EXISTED TRUE)
  else()
    set(TBB_TARGET_EXISTED FALSE)
  endif()
else()
  set(TBB_TARGET_EXISTED FALSE)
endif()

##
#
# Fetch and build Intel's Threading Building Blocks library.
#
##
if(NOT ${TBB_TARGET_EXISTED})

set(TBB_WILL_RECONFIGURE TRUE)
# Set the appropriate compiler
if(CLANG)
  set(TBB_COMPILER "clang")
else()
  set(TBB_COMPILER "gcc")
endif()

message("Build system will fetch and build Intel Threading Building Blocks")
message("==================================================================")
# These are useful for the custom install step we'll do later
set(TBB_SOURCE_DIR ${GAT_SOURCE_DIR}/external/oneTBB-2021.5.0)
set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)

if("${TBB_COMPILER}" STREQUAL "gcc")
  ## Don't know why it's a problem yet, but if we're using
  ## GCC, get rid of the DO_ITT_NOTIFY flag
  # set(TBB_CXXFLAGS "${TBB_CXXFLAGS} -UDO_ITT_NOTIFY")
endif()

set(TBB_CXXFLAGS "${TBB_CXXFLAGS} ${CXXSTDFLAG} ${SCHAR_FLAG}")

ExternalProject_Add(libtbb
DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
DOWNLOAD_COMMAND curl -k -L https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2021.5.0.tar.gz -o v2021.5.tar.gz &&
${SHASUM} e5b57537c741400cf6134b428fc1689a649d7d38d9bb9c1b6d64f092ea28178a v2021.5.tar.gz &&
tar -xzvf v2021.5.tar.gz
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

  if(${FETCHED_BOOST})
    add_dependencies(libtbb libboost)
  endif()

endif() # end of fetch tbb

##
# If we're fetching tbb, we need to have dummy paths for these variables
# so that CMake won't complain
##
if(TBB_WILL_RECONFIGURE)
    set(TBB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    set(TBB_INCLUDE_DIRS ${TBB_INSTALL_DIR}/include)
    set(TBB_INCLUDE_DIR ${TBB_INSTALL_DIR}/include)
    set(TBB_LIBRARY_DIRS ${TBB_INSTALL_DIR}/lib)
    set(TBB_LIBRARY ${TBB_INSTALL_DIR}/lib)
    set(TBB_LIB_DIR ${TBB_INSTALL_DIR}/lib)
    #set(TBB_LIBRARIES tbb tbbmalloc)
    set(TBB_LIBRARIES ${TBB_INSTALL_DIR}/lib/libtbb.${SHARED_LIB_EXTENSION}
    			${TBB_INSTALL_DIR}/lib/libtbbmalloc.${SHARED_LIB_EXTENSION}
          ${TBB_INSTALL_DIR}/lib/libtbbmalloc_proxy.${SHARED_LIB_EXTENSION}
	   )
    message("TBB_INCLUDE_DIRS = ${TBB_INCLUDE_DIRS}")
    message("TBB_LIBRARY_DIRS = ${TBB_LIBRARY_DIRS}")
endif()

##
#  Similar to the Boost trick above, the libtbb reconfigure should force this code
#  to be run on the second configuration pass, where it should appropriately set the
#  TBB_INSTALL_DIR variable.
##
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
  message("TBB_INSTALL_DIR = ${TBB_INSTALL_DIR}")
  find_package(TBB 2021.4
  HINTS ${TBB_ROOT_SEARCH}
  COMPONENTS tbb tbbmalloc tbbmalloc_proxy)
  message("[in TBB_RECONFIGURE] TBB_LIBRARIES = ${TBB_LIBRARIES}")
endif()

#message("TBB_LIBRARIES = ${TBB_LIBRARIES}")
#message("TBB_FOUND        ${TBB_FOUND}       ")
#message("TBB_INSTALL_DIR  ${TBB_INSTALL_DIR}")
#message("TBB_INCLUDE_DIRS ${TBB_INCLUDE_DIRS}")
#message("TBB_INCLUDE_DIR  ${TBB_INCLUDE_DIR} ")
#message("TBB_LIBRARY_DIRS ${TBB_LIBRARY_DIRS}")
#message("TBB_LIBRARIES    ${TBB_LIBRARIES}   ")

find_package(libgff 2.0.0 
HINTS ${LIB_GFF_PATH} ${GFF_ROOT}
)
if(libgff_FOUND)
  message(STATUS "libgff ver. ${LIB_GFF_VERSION} found.")
  message(STATUS "    include: ${LIB_GFF_INCLUDE_DIR}")
  message(STATUS "    lib    : ${LIB_GFF_LIBRARY_DIR}")
endif()

if(NOT libgff_FOUND)
  message("Build system will compile libgff")
  message("==================================================================")
  externalproject_add(libgff
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/COMBINE-lab/libgff/archive/v2.0.0.tar.gz -o libgff.tgz &&
      ${SHASUM} 7656b19459a7ca7d2fd0fcec4f2e0fd0deec1b4f39c703a114e8f4c22d82a99c libgff.tgz &&
      tar -xzvf libgff.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff-2.0.0
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff-2.0.0/build
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} 
  )
  externalproject_add_step(libgff makedir
    COMMAND mkdir -p <SOURCE_DIR>/build
    COMMENT "Make build directory"
    DEPENDEES download
    DEPENDERS configure)
  set(FETCHED_GFF TRUE)
  set(LIB_GFF_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
endif()

# Because of the way that Apple has changed SIP
# in el capitan, some headers may be in a new location
if(APPLE)
  set(STADEN_INC "-I/usr/local/include")
  set(STADEN_LIB "-L/usr/local/lib")
endif()

if(CONDA_BUILD)
  set(LZFLAG "-lz")
else()
  set(LZFLAG "")
endif()

find_package(CURL)

if (FETCH_STADEN)
  set(LIBSTADEN_FOUND FALSE)
else ()
  find_package(libstadenio 1.14.15)
endif()

if (NOT LIBSTADENIO_FOUND)
  message("Build system will compile Staden IOLib")
  message("==================================================================")
  externalproject_add(libstadenio
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/jkbonfield/io_lib/releases/download/io_lib-1-14-15/io_lib-1.14.15.tar.gz -o staden-io_lib-v1.14.15.tar.gz &&
      ${SHASUM} 20814c4365e1e2fe6630fb11d0df370dec4c5688af3871de7f1cb0129671401e staden-io_lib-v1.14.15.tar.gz &&
      mkdir -p staden-io_lib-1.14.15 &&
      tar -xzf staden-io_lib-v1.14.15.tar.gz --strip-components=1 -C staden-io_lib-1.14.15 &&
      rm -fr staden-io_lib &&
      mv -f staden-io_lib-1.14.15 staden-io_lib
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/staden-io_lib
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ./configure --enable-shared=no --without-libcurl --prefix=<INSTALL_DIR> LDFLAGS=${LIBSTADEN_LDFLAGS} CFLAGS=${LIBSTADEN_CFLAGS} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_COMMAND make ${QUIET_MAKE} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS+=${STADEN_INC} CFLAGS+=${STADEN_LIB} LDFLAGS+=${EXTRA_CMAKE_LIBRARY_FLAGS} CFLAGS+=${EXTRA_CMAKE_INCLUDE_FLAGS} CFLAGS+=${LZFLAG} CFLAGS+=${SCHAR_FLAG}


    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
  )
  if(NOT LIBLZMA_FOUND)
  	ExternalProject_Add_StepDependencies(libstadenio build liblzma)
  endif()

  set(FETCHED_STADEN TRUE)
  set(STADEN_LIBRARIES "${GAT_SOURCE_DIR}/external/install/lib/libhtscodecs.a;${GAT_SOURCE_DIR}/external/install/lib/libstaden-read.a")
endif()

if (ASAN_BUILD)
  set(FAST_MALLOC_LIB "")
  set(HAVE_FAST_MALLOC TRUE)
else()
  set(FAST_MALLOC_LIB "")
  set(HAVE_FAST_MALLOC FALSE)

  # See if we have Jemalloc
  find_package(Jemalloc)
  if(Jemalloc_FOUND)
    ##
    # Don't be so stringent about the version yet
    ##
    #if (NOT (${JEMALLOC_VERSION} VERSION_LESS 5.2.1))
      message("Found Jemalloc library --- using this memory allocator")
      set(FAST_MALLOC_LIB ${JEMALLOC_LIBRARIES})
      set(HAVE_FAST_MALLOC TRUE)
    #else()
    #  message("Fond Jemalloc version ${JEMALLOC_VERSION}, but require >= 5.2.1. Downloading newer version")
    #endif()
  endif()
endif()

if(CONDA_BUILD)
  set(JEMALLOC_FLAGS "CC=${CMAKE_C_COMPILER} CFLAGS=\\\"-fPIC ${SCHAR_FLAG}\\\" CPPFLAGS=\\\"-fPIC ${SCHAR_FLAG}\\\"")
else()
  set(JEMALLOC_FLAGS "CC=${CMAKE_C_COMPILER} CFLAGS=${SCHAR_FLAG} CPPFLAGS=${SCHAR_FLAG}")
endif()

if(NOT HAVE_FAST_MALLOC)
  message("Build system will fetch and use JEMalloc")
  message("==================================================================")
  externalproject_add(libjemalloc
      DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
      DOWNLOAD_COMMAND curl -k -L https://github.com/jemalloc/jemalloc/archive/5.2.1.tar.gz -o jemalloc-5.2.1.tar.gz &&
        ${SHASUM} ed51b0b37098af4ca6ed31c22324635263f8ad6471889e0592a9c0dba9136aea jemalloc-5.2.1.tar.gz &&
        tar -xzf jemalloc-5.2.1.tar.gz

      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/jemalloc-5.2.1
      BUILD_IN_SOURCE TRUE
      INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
      CONFIGURE_COMMAND sh -c "${JEMALLOC_FLAGS} ./autogen.sh --disable-debug ${MALLOC_STATIC_BUILD_FLAG} --prefix=<INSTALL_DIR>"
      #CONFIGURE_COMMAND sh -c "${JEMALLOC_FLAGS} ./autogen.sh ${MALLOC_STATIC_BUILD_FLAG} --prefix=<INSTALL_DIR>"
      INSTALL_COMMAND cp -r lib <INSTALL_DIR>/ && cp -r include <INSTALL_DIR>/
      )

  set(FAST_MALLOC_LIB ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib/libjemalloc.a)
  set(HAVE_FAST_MALLOC TRUE)
  set(FETCHED_JEMALLOC TRUE)
  if(FETCHED_LIBBZ2)
    add_dependencies(libjemalloc libbz2)
  endif()

  if(FETCHED_LIBLZMA)
    add_dependencies(libjemalloc liblzma)
  endif()
endif()

###
#
# Done building external dependencies.
#
###

set(CPACK_SOURCE_IGNORE_FILES
"/src/PCA.cpp"
"/src/PCAUtils.cpp"
"/build/"
"/scripts/AggregateToGeneLevel.py"
"/scripts/ExpressionTools.py"
"/scripts/GenerateExpressionFiles.sh"
"/scripts/ParseSoftFile.py"
"/scripts/PlotCorrelation.py"
"/scripts/junk"
"/scripts/sfstrace.log"
"/scripts/SFPipeline.py"
"/bin/"
"/lib/"
"/sample_data/"
"PublishREADMEToWebsite.sh"
"/external/"
"/src/obsolete/"
"/include/obsolete/"
"WebsiteHeader.txt"
"/experimental_configs/"
".git/")

message("CPACK_SOURCE_IGNORE_FILES = ${CPACK_SOURCE_IGNORE_FILES}")

# we will use this property later
define_property(TARGET PROPERTY COMPACT_VECTOR_DIR INHERITED 
   BRIEF_DOCS "the path to the directory containing the compact_vector include tree"
   FULL_DOCS "the path to the directory containing the compact_vector include tree")

# Recurse into pufferfish source directory 
# and build the library 
set(BUILD_PUFF_FOR_SALMON TRUE)
add_subdirectory(external/pufferfish)

# make sure we know the path to compact_vector
get_property(COMPACT_VECTOR_INCLUDE_PATH TARGET graphdump PROPERTY COMPACT_VECTOR_DIR)
message("fetched path for compact_vector as [${COMPACT_VECTOR_INCLUDE_PATH}]")

# and then the main salmon source directory
add_subdirectory(src)

#add_dependencies(salmon RapMap)
# build a CPack driven installer package
include(CPack)

set(ARCHIVE_NAME ${CMAKE_PROJECT_NAME}-${PROJECT_VERSION})
add_custom_target(dist
  COMMAND git archive --prefix=${ARCHIVE_NAME}/ HEAD
      | gzip > ${CMAKE_BINARY_DIR}/${ARCHIVE_NAME}.tar.gz
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
