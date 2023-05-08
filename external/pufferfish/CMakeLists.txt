
cmake_minimum_required(VERSION 3.9)
#cmake_policy(SET CMP0069 NEW)

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

if(NOT DEFINED USE_ARM)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*|ARM64.*)")
    message("Detected 64-bit ARM host. Setting USE_ARM to true.")
    set(USE_ARM TRUE)
  else()
    message("Detected non-ARM host. Setting USE_ARM to false.")
    set(USE_ARM FALSE)
  endif()
endif()

project(pufferfish LANGUAGES CXX C)
include(ExternalProject)


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
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Pufferfish — a tool for indexing and alignment based on the compacted colored de Bruijn Graph")
set(CPACK_PACKAGE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-Source")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

if (NOT CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE "Release")
endif()

# We require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -static-libgcc -static-libstdc++")
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
set(CXXSTDFLAG "-std=c++14")

set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Threads REQUIRED)

if (BUILD_PUFF_FOR_SALMON) 
    message("Building basic pufferfish components for salmon")
else()
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "^(Clang|AppleClang)$")
  	set(CLANG TRUE)
  endif()
endif()

if (NOT GAT_SOURCE_DIR)
    message("Top-level source directory variable not set externally; setting it to ${CMAKE_CURRENT_SOURCE_DIR}")
    set(GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif()

set(KSW_FLAGS "-DHAVE_KALLOC")
#set(PF_CPP_FLAGS "-msse4;-W;-Wall;-Wextra;-Wpointer-arith;-Wunused;-Wwrite-strings;-Wno-unknown-pragmas;-Wno-unused-function")
set(PF_CPP_FLAGS "-W;-Wall;-Wextra;-Wpointer-arith;-Wunused;-Wwrite-strings;-Wno-unknown-pragmas;-Wno-unused-function")

if (USE_ARM)
    list(APPEND PF_CPP_FLAGS "-fsigned-char")
    add_compile_options(-fsigned-char)
endif()

if (BUILD_PUFF_FOR_SALMON)
    #list(APPEND PF_CPP_FLAGS "-Wno-unused-parameter") 
    #list(APPEND PF_CPP_FLAGS "-Wno-unused-variable") 
    list(APPEND PF_CPP_FLAGS "-Wno-reorder") 
    #list(APPEND PF_CPP_FLAGS "-g")
    #list(APPEND PF_CPP_FLAGS "-Wno-sign-compare") 
    #list(APPEND PF_CPP_FLAGS "-Wno-sequence-point")
   list(APPEND PF_CPP_FLAGS "-DPUFFERFISH_SALMON_SUPPORT")
endif()

###
# Sanitizers BEGIN
###
if (ASAN_BUILD)
  list(APPEND PF_CPP_FLAGS "-fsanitize=address")
  list(APPEND PF_CPP_FLAGS "-fsanitize=memory")
#list(APPEND PF_CPP_FLAGS "-fsanitize=undefined")
#set(CMAKE_LINK_FLAGS "-fsanitize=address")
#list(APPEND CMAKE_LINK_FLAGS "-fsanitize=undefined")
# list(APPEND PF_CPP_FLAGS "-g")
  set(ASAN_LIB "asan")
else()
  set(ASAN_LIB "")
endif()
###
# Sanitizers END
###


set(WARN_ALL_THINGS "-fdiagnostics-color=always;-Wall;-Wcast-align;-Wcast-qual;-Wconversion;-Wctor-dtor-privacy;-Wdisabled-optimization;-Wdouble-promotion;-Wextra;-Wformat=2;-Winit-self;-Wlogical-op;-Wmissing-declarations;-Wmissing-include-dirs;-Wno-sign-conversion;-Wnoexcept;-Wold-style-cast;-Woverloaded-virtual;-Wpedantic;-Wredundant-decls;-Wshadow;-Wstrict-aliasing=1;-Wstrict-null-sentinel;-Wstrict-overflow=5;-Wswitch-default;-Wundef;-Wno-unknown-pragmas;-Wuseless-cast;-Wno-unused-parameter")

#set(WARN_ALL_THINGS "-fdiagnostics-color=always -Wall -Wcast-align -Wcast-qual -Wconversion -Wctor-dtor-privacy -Wdisabled-optimization -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Wextra -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wno-sign-conversion -Wnoexcept -Wnull-dereference -Wold-style-cast -Woverloaded-virtual -Wpedantic -Wredundant-decls -Wrestrict -Wshadow -Wstrict-aliasing=1 -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unknown-pragmas -Wuseless-cast") 



set(OPT_FLAGS "-DHAVE_SIMDE;-D__STDC_FORMAT_MACROS;-DSTX_NO_STD_STRING_VIEW;-O3;-fPIC;-DNDEBUG;-funroll-loops;-ftree-vectorize;-fno-strict-aliasing")
set(DEBUG_FLAGS "-DHAVE_SIMDE;-D__STDC_FORMAT_MACROS;-DSTX_NO_STD_STRING_VIEW;-pg;-g;-gstabs")

##
# OSX is strange (some might say, stupid in this regard).  Deal with it's quirkines here.
##
if(APPLE)
    # To allow ourselves to build a dynamic library, we have to tell the compiler
    # that, yes, the symbols will be around at runtime.
    list(APPEND TGT_COMPILE_FLAGS "-undefined dynamic_lookup;-Wno-unused-command-line-argument")
    set(LIBSALMON_LINKER_FLAGS "-all_load")
    # In order to "think different", we also have to use non-standard suffixes
    # for our shared libraries
    set(SHARED_LIB_EXTENSION "dylib")

	message("SHARED_LIB_EXT = ${SHARED_LIB_EXTENSION}")

else()
    # We're in sane linux world
    set(SHARED_LIB_EXTENSION "so")
    set(LIBSALMON_LINKER_FLAGS "")
#    if (NOT BUILD_PUFF_FOR_SALMON)
#        set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
#        list(APPEND PF_CPP_FLAGS -static-libstdc++)
#        message("PF_CPP_FLAGS = ${PF_CPP_FLAGS}")
#    endif()
endif()

###
# check if numeric_limits<__int128_t> is defined
###
if(HAVE_INT128_NUMERIC_LIMITS)

  message("setting -DHAVE_NUMERIC_LIMITS128")
  list(APPEND PF_CPP_FLAGS "-DHAVE_NUMERIC_LIMITS128")
else()
  try_compile(HAVE_INT128_NUMERIC_LIMITS ${CMAKE_BINARY_DIR} 
              SOURCES ${GAT_SOURCE_DIR}/tests/compile_tests/int128_numeric_limits.cpp
              CXX_STANDARD 14
              CXX_STANDARD_REQUIRED ON
              )
  if(HAVE_INT128_NUMERIC_LIMITS)
    message("setting -DHAVE_NUMERIC_LIMITS128")
    list(APPEND PF_CPP_FLAGS "-DHAVE_NUMERIC_LIMITS128")
  else()
    message("not setting -DHAVE_NUMERIC_LIMITS128")
  endif()
endif()

##
# Let us check the sha sum of our packages if we have the right tools
##
set(SHASUM ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_shasum.sh)


find_package(Jemalloc)
if(JEMALLOC_FOUND)
  include_directories(SYSTEM ${JEMALLOC_INCLUDE_DIRS})
endif()

if(DO_QUIET_MAKE)
  set(QUIET_MAKE "--silent")
else()
  set(QUIET_MAKE "")
endif()


if (NOT BUILD_PUFF_FOR_SALMON)
## Try and find TBB first
find_package(TBB 2021.4 
  HINTS ${TBB_ROOT_SEARCH}
  COMPONENTS tbb tbbmalloc tbbmalloc_proxy)


if (${TBB_FOUND})
  if (${TBB_VERSION} VERSION_GREATER_EQUAL 2021.4)
    message("FOUND SUITABLE TBB VERSION : ${TBB_VERSION}")
    set(TBB_TARGET_EXISTED TRUE)
    get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
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
    set(TBB_INSTALL_DIR ${GAT_SOURCE_DIR}/external/install)
    message("fetching oneAPI::TBB : TBB_SOURCE_DIR = ${TBB_SOURCE_DIR}, TBB_INSTALL_DIR = ${TBB_INSTALL_DIR}")

    if(NO_RTM)
  	set(TBB_CXXFLAGS "-mno-rtm")
    endif()

    if("${TBB_COMPILER}" STREQUAL "gcc")
        ## Don't know why it's a problem yet, but if we're using
        ## GCC, get rid of the DO_ITT_NOTIFY flag
        # set(TBB_CXXFLAGS "${TBB_CXXFLAGS} -UDO_ITT_NOTIFY")
    endif()

    set(TBB_CXXFLAGS "${TBB_CXXFLAGS} ${CXXSTDFLAG}")

    ExternalProject_Add(libtbb
            DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
            DOWNLOAD_COMMAND curl -k -L https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2021.5.0.tar.gz -o v2021.5.tar.gz &&
            ${SHASUM} e5b57537c741400cf6134b428fc1689a649d7d38d9bb9c1b6d64f092ea28178a v2021.5.tar.gz &&
            tar -xzvf v2021.5.tar.gz
            SOURCE_DIR ${TBB_SOURCE_DIR}
            INSTALL_DIR ${TBB_INSTALL_DIR}
            CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_FLAGS=${TBB_CXXFLAGS} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DTBB_TEST=OFF -DTBB_EXAMPLES=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=14
            BUILD_IN_SOURCE TRUE
            )

    set(RECONFIG_FLAGS ${RECONFIG_FLAGS} -DTBB_WILL_RECONFIGURE=FALSE -DTBB_RECONFIGURE=TRUE)
    ExternalProject_Add_Step(libtbb reconfigure
    COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
    DEPENDEES install 
    )
    set(FETCHED_TBB TRUE)
    set(TBB_ROOT_SEARCH ${CMAKE_SOURCE_DIR}/external/install)    
endif()

find_package(LibLZMA)
if(NOT LIBLZMA_FOUND)

  message("Will attempt to fetch and build liblzma")

  message("=======================================")
  externalproject_add(liblzma
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
  DOWNLOAD_COMMAND wget http://tukaani.org/xz/xz-5.2.2.tar.gz -o xz-5.2.2.tar.gz &&
    ${SHASUM} 73df4d5d34f0468bd57d09f2d8af363e95ed6cc3a4a86129d2f2c366259902a2 xz-5.2.2.tar.gz &&
    tar -xzvf xz-5.2.2.tar.gz
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



  #endif()

##
# If we're fetching tbb, we need to have dummy paths for these variables
# so that CMake won't complain
##
if(TBB_WILL_RECONFIGURE)
    set(TBB_INCLUDE_DIRS ${TBB_INSTALL_DIR}/include)
    set(TBB_INCLUDE_DIR ${TBB_INSTALL_DIR}/include)
    set(TBB_LIBRARY_DIRS ${TBB_INSTALL_DIR}/lib)
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

#  The libtbb reconfigure should force this code
#  to be run on the second configuration pass, where it should appropriately set the
#  TBB_INSTALL_DIR variable.
##
if(TBB_RECONFIGURE)
    unset(TBB_FOUND CACHE)
    unset(TBB_INSTALL_DIR CACHE)
    unset(CMAKE_PREFIX_PATH CACHE)

    set(CMAKE_PREFIX_PATH ${GAT_SOURCE_DIR}/external/install)
    set(TBB_INSTALL_DIR ${GAT_SOURCE_DIR}/external/install)

    message("TBB_INSTALL_DIR = ${TBB_INSTALL_DIR}")
  
    find_package(TBB 2021.4 
      HINTS ${TBB_ROOT_SEARCH}
      COMPONENTS tbb tbbmalloc tbbmalloc_proxy
    )
    
    set(TBB_LIBRARIES 
		${TBB_INSTALL_DIR}/lib/libtbbmalloc.${SHARED_LIB_EXTENSION}
		${TBB_INSTALL_DIR}/lib/libtbbmalloc_proxy.${SHARED_LIB_EXTENSION}
    ${TBB_INSTALL_DIR}/lib/libtbb.${SHARED_LIB_EXTENSION}
	   )
endif()
message("TBB_LIBRARIES = ${TBB_LIBRARIES}")

ExternalProject_Add(libseqlib
GIT_REPOSITORY https://github.com/COMBINE-lab/SeqLib.git
GIT_TAG        master
UPDATE_COMMAND ""
UPDATE_DISCONNECTED 1
BUILD_IN_SOURCE TRUE
DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/seqlib
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/seqlib
INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
CONFIGURE_COMMAND ./configure 
BUILD_COMMAND     make CXXFLAGS='-std=c++14' ${LIBSTADEN_LDFLAGS} ${LIBSTADEN_CFLAGS}
INSTALL_COMMAND   mkdir -p <INSTALL_DIR>/lib && mkdir -p <INSTALL_DIR>/include && cp src/libseqlib.a <INSTALL_DIR>/lib && 
                  cp htslib/libhts.a <INSTALL_DIR>/lib &&
                  cp -r SeqLib <INSTALL_DIR>/include &&
                  cp -r json <INSTALL_DIR>/include &&
                  cp -r htslib <INSTALL_DIR>/include
)

if (FETCHED_TBB)
add_dependencies(libseqlib libtbb)
endif()

if(NOT LIBLZMA_FOUND)
  ExternalProject_Add_StepDependencies(libseqlib build liblzma)
endif()


if(NOT BZIP2_FOUND)
  ExternalProject_Add_StepDependencies(libseqlib build bz2)
endif()

ExternalProject_Add(libsetcover
        GIT_REPOSITORY https://github.com/martin-steinegger/setcover.git
        GIT_TAG        master
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        UPDATE_DISCONNECTED 1
        BUILD_IN_SOURCE TRUE
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/setcover
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/setcover
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        BUILD_COMMAND     make -j8
        INSTALL_COMMAND     ar rcs libSetCover.a SetCover.o &&
        mkdir -p <INSTALL_DIR>/lib && mkdir -p <INSTALL_DIR>/include &&
        cp libSetCover.a <INSTALL_DIR>/lib &&
        cp SetCover.h <INSTALL_DIR>/include &&
        cp SetElement.h <INSTALL_DIR>/include &&

        cp LinearMultiArray.h <INSTALL_DIR>/include
        )
if(FETCHED_TBB)
add_dependencies(libsetcover libtbb)
endif()

endif() # end of NOT BUILD_PUFF_FOR_SALMON

set(LIB_FLAGS "-L ${CMAKE_CURRENT_SOURCE_DIR}/external/lib")
set(PUFF_RELEASE_FLAGS "${DEFINE_FLAGS};${KSW_FLAGS};${PF_CPP_FLAGS};${OPT_FLAGS}")
set(PUFF_DEBUG_FLAGS "${DEFINE_FLAGS};${KSW_FLAGS};${PF_CPP_FLAGS};${DEBUG_FLAGS}")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/ntcard)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/twopaco)

if (NOT BUILD_PUFF_FOR_SALMON) 
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/asyncplusplus)

  if (NOT TBB_TARGET_EXISTED)
    if (FETCHED_TBB)
    add_dependencies(Async++ libtbb)
  endif()
  endif()

endif()

include(CPack)

include_directories(include)
link_directories(lib)
add_subdirectory(src)
