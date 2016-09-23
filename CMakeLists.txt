cmake_minimum_required (VERSION 2.8)

enable_testing()

project (Salmon)

set(CPACK_PACKAGE_VERSION "0.7.3")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "7")
set(CPACK_PACKAGE_VERSION_PATCH "3")
set(PROJECT_VERSION ${CPACK_PACKAGE_VERSION})
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "Stony Brook University")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Salmon - Wicked-fast RNA-seq isoform quantification using lightweight mapping")
set(CPACK_PACKAGE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-Source")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set (WARNING_IGNORE_FLAGS "-Wno-deprecated-register")
set (BOOST_CXX_FLAGS "-Wno-deprecated-register -std=c++11")
## Prefer static to dynamic libraries
SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

## Set the standard required compile flags
# Nov 18th --- removed -DHAVE_CONFIG_H
set (CMAKE_CXX_FLAGS "-pthread -ftree-vectorize -funroll-loops -fPIC -fomit-frame-pointer -Ofast -DRAPMAP_SALMON_SUPPORT -DHAVE_ANSI_TERM -DHAVE_SSTREAM -Wall -Wno-unknown-pragmas -Wno-reorder -Wno-unused-variable -std=c++11 -Wreturn-type -Werror=return-type")

##
# OSX is strange (some might say, stupid in this regard).  Deal with it's quirkines here.
##
if (APPLE)
    # To allow ourselves to build a dynamic library, we have to tell the compiler
    # that, yes, the symbols will be around at runtime.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
    set (LIBSALMON_LINKER_FLAGS "-all_load")
    # In order to "think different", we also have to use non-standard suffixes
    # for our shared libraries
    set(SHARED_LIB_EXTENSION "dylib")
else()
    # We're in sane linux world
   set (SHARED_LIB_EXTENSION "so")
   set (LIBSALMON_LINKER_FLAGS "")
endif()

set( BOOST_EXTRA_FLAGS "--layout=tagged" )
## this get's set differently below if we
## are on clang & apple
set (NON_APPLECLANG_LIBS gomp rt)
set (PTHREAD_LIB)

##
# Compiler-specific C++11 activation.
# http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake
##
##
# First take care of what to do if we have gcc
##
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    # If we're on OSX
    if (APPLE AND NOT (GCC_VERSION VERSION_GREATER 4.8.2 OR GCC_VERSION VERSION_EQUAL 4.8.2))
        message(FATAL_ERROR "When building under OSX, ${PROJECT_NAME} requires "
                            "either clang or g++ >= 4.8.2")
    elseif (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()

    set (GCC TRUE)

    # Put complete static linking on hold for the time-being
    # If we're not on OSX, make an attempt to compile everything statically
    #if (NOT APPLE)
    #set (CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
    #set (CMAKE_EXE_LINK_FLAGS "-static")
    set (PTHREAD_LIB "pthread")
    #endif()

    # If we're on Linux (i.e. not OSX) and we're using
    # gcc, then set the -static-libstdc++ flag
    if (NOT APPLE)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")
    endif()

    set (WARNING_IGNORE_FLAGS "${WARNING_IGNORE_FLAGS} -Wno-unused-local-typedefs")
    set (BOOST_TOOLSET "gcc")
    set (BOOST_CONFIGURE_TOOLSET "--with-toolset=gcc")
	set (BCXX_FLAGS "-std=c++11")
    set (BOOST_EXTRA_FLAGS toolset=gcc cxxflags=${BCXX_FLAGS})
# Tentatively, we support clang now
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CLANG TRUE)
    # If we have libc++, then try and use it
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-stdlib=libc++ HAVE_LIBCPP)
    if (HAVE_LIBCPP)
        message ("It appears that you're compiling with clang and that libc++ is available, so I'll use that")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
	    set (BOOST_TOOLSET "clang")
        set (BOOST_CONFIGURE_TOOLSET "--with-toolset=clang")
	    set (BCXX_FLAGS "-stdlib=libc++ -DBOOST_HAS_INT128")
	    set (BOOST_EXTRA_FLAGS toolset=clang cxxflags=${BCXX_FLAGS} linkflags="-stdlib=libc++")
        set (JELLYFISH_CXX_FLAGS "-stdlib=libc++")
    # Otherwise, use libstdc++ (and make it static)
    else()
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")
    endif()
    # There's currently a bug with clang-3.4 & Boost 1.55 -- this hack fixes it
    # but we should do something better (does this break things if CPU doesn't
    # have 128-bit support)?
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_HAS_INT128")

    if (APPLE)
        set (NON_APPLECLANG_LIBS "")
    else()
        set (PTHREAD_LIB "pthread")
    endif()
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

## TODO: Figure out how to detect this automatically
# If the "assembler" is too old, tell TBB not to compile
# with -mrtm
if (NO_RTM)
    set (TBB_CXXFLAGS "-mno-rtm")
endif()

include(ExternalProject)

##
#  Update the CXX flags according to the system and compiler
##
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_IGNORE_FLAGS}")

if (CMAKE_BUILD_TYPE MATCHES Debug)
    message ("Making Debug build")
    set (CMAKE_CXX_FLAGS_DEBUG "-g ${CMAKE_CXX_FLAGS}")
elseif (CMAKE_BUILD_TYPE MATCHES Release)
    message ("Making Release build")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS}")
else ()
    message ("Making Default build type")
endif ()

##
# Record this top-level path
##
set (GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Have CMake tell us what it's doing
set (CMAKE_VERBOSE_MAKEFILE true)

###
#
#  Grab RapMap sources for quasi-mapping code --- DURING CONFIGURE TIME!
#
####
if(NOT FETCHED_RAPMAP)
  exec_program(${CMAKE_CURRENT_SOURCE_DIR}/scripts/fetchRapMap.sh)
  set(FETCHED_RAPMAP TRUE CACHE BOOL "Has RapMap been fetched?" FORCE)
endif()

##
# Super-secret override
##
if ( DEFINED CUSTOM_BOOST_PATH )
	set (CMAKE_INCLUDE_PATH ${CUSTOM_BOOST_PATH} ${CMAKE_INCLUDE_PATH})
    set (CMAKE_LIBRARY_PATH ${CUSTOM_BOOST_PATH}/lib ${CMAKE_LIBRARY_PATH})
endif ( DEFINED CUSTOM_BOOST_PATH )

##
# We want static, multithreaded boost libraries
##
if(CONDA_BUILD)
  set (Boost_USE_STATIC_LIBS OFF)
else ()
  set (Boost_USE_STATIC_LIBS ON)
endif(CONDA_BUILD)

set (Boost_USE_MULTITHREADED ON)
#set (Boost_USE_STATIC_RUNTIME OFF)

find_package (ZLIB)
if (NOT ZLIB_FOUND)
	message (FATAL_ERROR "zlib must be installed before configuration & building can proceed")
endif()

find_package (LibLZMA)
if (NOT LIBLZMA_FOUND)
    message ("Will attempt to fetch and build liblzma")
    message ("=======================================")
ExternalProject_Add(liblzma
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://tukaani.org/xz/xz-5.2.2.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/xz-5.2.2
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/xz-5.2.2/configure --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)
# Tell cmake that the external project generated a library so we can
# add dependencies here instead of later
set (LIBLZMA_LIBRARIES ${GAT_SOURCE_DIR}/external/install/lib/liblzma.a)
set (LIBSTADEN_LDFLAGS "-L${GAT_SOURCE_DIR}/external/install/lib")
set (LIBSTADEN_CFLAGS "-I${GAT_SOURCE_DIR}/external/install/include")
else()
    message("Found liblzma library: ${LIBLZMA_LIBRARIES}")
    message("===========================================")
endif()

find_package (BZip2)
if (NOT BZIP2_FOUND)
    message ("Will attempt to fetch and build libbz2")
    message ("=======================================")
ExternalProject_Add(libbz2
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bzip2-1.0.6
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    INSTALL_COMMAND make install PREFIX=<INSTALL_DIR>
)
# Tell cmake that the external project generated a library so we can
# add dependencies here instead of later
set (BZIP2_LIBRARIES ${GAT_SOURCE_DIR}/external/install/lib/libbz2.a)
set (LIBSTADEN_LDFLAGS "-L${GAT_SOURCE_DIR}/external/install/lib -I${GAT_SOURCE_DIR}/external/install/include")
set (LIBSTADEN_CFLAGS "-I${GAT_SOURCE_DIR}/external/install/include")
else()
    message("Found libbz2 library: ${BZIP2_LIBRARIES}")
    message("===========================================")
endif()

##
# Set the latest version and look for what we need
##
set(Boost_ADDITIONAL_VERSIONS "1.53" "1.53.0" "1.54" "1.55" "1.56" "1.57.0" "1.58" "1.59" "1.60" "1.61")
find_package(Boost 1.53.0 COMPONENTS iostreams filesystem system thread timer chrono program_options serialization)
message("BOOST_INCLUDEDIR = ${BOOST_INCLUDEDIR}")
message("BOOST_LIBRARYDIR = ${BOOST_LIBRARYDIR}")
message("Boost_FOUND = ${Boost_FOUND}")
include(ExternalProject)

##
#  If we had to fetch Boost, the reconfigure step will re-run cmake.  The second configuration
#  pass is executed with the BOOST_RECONFIGURE flag set. This should allow our newly
#  installed Boost to be found by CMake.
##
if (BOOST_RECONFIGURE)
    message("Executing Boost Reconfiguration")
    unset(Boost_FOUND CACHE)
    unset(Boost_INCLUDE_DIR CACHE)
    unset(Boost_INCLUDE_DIRS CACHE)
    unset(Boost_LIBRARY_DIRS CACHE)
    unset(Boost_LIBRARIES CACHE)
    unset(BOOST_ROOT CACHE)
    unset(CMAKE_PREFIX_PATH CACHE)

    set(BOOST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
    find_package(Boost 1.53.0 COMPONENTS iostreams filesystem system thread timer chrono program_options REQUIRED)
    set(FETCH_BOOST FALSE)
endif()



##
# Either inform the user of how to obtain Boost, or, if they passed in the FETCH_BOOST
# option, go and grab it for them.
##
if ((NOT Boost_FOUND) AND (NOT FETCH_BOOST))
	message(FATAL_ERROR
        "Salmon cannot be compiled without Boost.\n"
        "It is recommended to visit http://www.boost.org/ and install Boost according to those instructions.\n"
        "This build system can also download and install a local version of boost for you (this takes a lot of time).\n"
        "To fetch and build boost locally, call cmake with -DFETCH_BOOST=TRUE"
    )
elseif(FETCH_BOOST)
    ## Let the rest of the build process know we're going to be fetching boost
    set (BOOST_LIB_SUBSET --with-iostreams --with-atomic --with-chrono --with-container --with-date_time --with-exception
                          --with-filesystem --with-graph --with-graph_parallel --with-math
                          --with-program_options --with-system --with-thread
                          --with-timer)
    set (BOOST_WILL_RECONFIGURE TRUE)
    set (FETCH_BOOST FALSE)
    message("Build system will fetch and build Boost")
    message("==================================================================")
    ExternalProject_Add(libboost
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
        DOWNLOAD_COMMAND curl -k -L  http://downloads.sourceforge.net/project/boost/boost/1.59.0/boost_1_59_0.tar.gz -o boost_1_59_0.tar.gz &&
                tar xzf boost_1_59_0.tar.gz
        #URL http://downloads.sourceforge.net/project/boost/boost/1.57.0/boost_1_57_0.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_59_0
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        #PATCH_COMMAND patch -p2 < ${CMAKE_CURRENT_SOURCE_DIR}/external/boost156.patch
	CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_59_0/bootstrap.sh ${BOOST_CONFIGURE_TOOLSET} ${BOOST_BUILD_LIBS} --prefix=<INSTALL_DIR>
	BUILD_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_59_0/b2 -d0 -j2 ${BOOST_LIB_SUBSET} toolset=${BOOST_TOOLSET} ${BOOST_EXTRA_FLAGS} cxxflags=${BOOST_CXX_FLAGS} link=static install
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ""
    )

    ##
    # After we've installed boost,
    ##
    SET( RECONFIG_FLAGS ${RECONFIG_FLAGS} -DBOOST_WILL_RECONFIGURE=FALSE -DBOOST_RECONFIGURE=TRUE -DFETCH_BOOST=FALSE)
    ExternalProject_Add_Step(libboost reconfigure
        COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
        DEPENDEES install
    )
endif()

##
# If we're fetching boost and we need to have dummy paths for these variables
# so that CMake won't complain
##
if (BOOST_WILL_RECONFIGURE)
    message("Setting Temporary Boost paths")
    set(Boost_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
    set(Boost_FOUND TRUE)
endif()


message("BOOST INCLUDE DIR = ${Boost_INCLUDE_DIR}")
message("BOOST INCLUDE DIRS = ${Boost_INCLUDE_DIRS}")
message("BOOST LIB DIR = ${Boost_LIBRARY_DIRS}")
message("BOOST LIBRARIES = ${Boost_LIBRARIES}")

set(EXTERNAL_LIBRARY_PATH $CMAKE_CURRENT_SOURCE_DIR/lib)

message("Build system will build libdivsufsort")
message("==================================================================")
include(ExternalProject)
ExternalProject_Add(libdivsufsort
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL ${CMAKE_CURRENT_SOURCE_DIR}/external/libdivsufsort.zip
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libdivsufsort-master
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libdivsufsort-master/build
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DBUILD_DIVSUFSORT64=TRUE -DUSE_OPENMP=TRUE -DBUILD_SHARED_LIBS=FALSE
)
set(SUFFARRAY_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)

message("Build system will fetch and build the Cereal serialization library")
message("==================================================================")
include(ExternalProject)
ExternalProject_Add(libcereal
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/USCiLab/cereal/archive/v1.1.2.tar.gz -o cereal-v1.1.2.tar.gz &&
		tar -xzvf cereal-v1.1.2.tar.gz
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cereal-1.1.2
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cereal-1.1.2/build
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND sh -c "mkdir -p <INSTALL_DIR>/include && cp -r <SOURCE_DIR>/include/cereal <INSTALL_DIR>/include"
)

message("Build system will fetch and build BWA (for Salmon)")
message("==================================================================")
include(ExternalProject)
ExternalProject_Add(libbwa
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/COMBINE-lab/bwa/archive/0.7.12.3.tar.gz -o bwa-master.tar.gz &&
                     mkdir -p bwa-master &&
                     tar -xzvf bwa-master.tar.gz --strip-components=1 -C bwa-master
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bwa-master
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ""
    BUILD_COMMAND sh -c "make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}"
    INSTALL_COMMAND sh -c "mkdir -p <INSTALL_DIR>/lib && mkdir -p <INSTALL_DIR>/include/bwa && cp libbwa.a <INSTALL_DIR>/lib && cp *.h <INSTALL_DIR>/include/bwa && cp is.c bwtindex.c bwt_gen.c QSufSort.c ${CMAKE_CURRENT_SOURCE_DIR}/src/"
    BUILD_IN_SOURCE TRUE
)

find_package(Jellyfish 2.2.6)

if (NOT JELLYFISH_FOUND)
message("Build system will fetch and build Jellyfish")
message("==================================================================")
ExternalProject_Add(libjellyfish
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz -o jellyfish-2.2.6.tgz &&
    	rm -fr jellyfish-2.2.6 &&
     	tar -xzvf jellyfish-2.2.6.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-2.2.6
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-2.2.6/configure --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${JELLYFISH_CXX_FLAGS}
    BUILD_COMMAND ${MAKE} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${JELLYFISH_CXX_FLAGS}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
)
endif()

find_package(TBB)

##
#
# Fetch and build Intel's Threading Building Blocks library.
#
##
if(NOT TBB_FOUND)

set(TBB_WILL_RECONFIGURE TRUE)
# Set the appropriate compiler
if (CLANG)
    set(TBB_COMPILER "clang")
else()
    set(TBB_COMPILER "gcc")
endif()

message("Build system will fetch and build Intel Threading Building Blocks")
message("==================================================================")
# These are useful for the custom install step we'll do later
set(TBB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb44_20160526oss)
set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)

if ("${TBB_COMPILER}" STREQUAL "gcc")
    ## Don't know why it's a problem yet, but if we're using
    ## GCC, get rid of the DO_ITT_NOTIFY flag
    set(TBB_CXXFLAGS "${TBB_CXXFLAGS} -UDO_ITT_NOTIFY")
endif()

ExternalProject_Add(libtbb
	DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL https://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb44_20160526oss_src_0.tgz
    DOWNLOAD_COMMAND curl -k -L https://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb44_20160526oss_src_0.tgz -o tbb_20160526oss_src.tgz &&
                    tar -xzvf tbb_20160526oss_src.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb44_20160526oss
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    PATCH_COMMAND "${TBB_PATCH_STEP}"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make CXXFLAGS=${TBB_CXXFLAGS} lambdas=1 compiler=${TBB_COMPILER} cfg=release tbb_build_prefix=LIBS
    INSTALL_COMMAND sh -c "cp ${TBB_SOURCE_DIR}/build/LIBS_release/*.${SHARED_LIB_EXTENSION}* ${TBB_INSTALL_DIR}/lib && cp -r ${TBB_SOURCE_DIR}/include/* ${TBB_INSTALL_DIR}/include"
    BUILD_IN_SOURCE 1
)

SET( RECONFIG_FLAGS ${RECONFIG_FLAGS} -DTBB_WILL_RECONFIGURE=FALSE -DTBB_RECONFIGURE=TRUE)
ExternalProject_Add_Step(libtbb reconfigure
        COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
        DEPENDEES install
)
endif()

##
# If we're fetching tbb, we need to have dummy paths for these variables
# so that CMake won't complain
##
if(TBB_WILL_RECONFIGURE)
    set(TBB_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(TBB_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
    #set(TBB_LIBRARIES tbb tbbmalloc)
    set(TBB_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib/libtbb.${SHARED_LIB_EXTENSION}
                      ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib/libtbbmalloc.${SHARED_LIB_EXTENSION}
    )
endif()

##
#  Similar to the Boost trick above, the libtbb reconfigure should force this code
#  to be run on the second configuration pass, where it should appropriately set the
#  TBB_INSTALL_DIR variable.
##
if (TBB_RECONFIGURE)
    unset(TBB_FOUND CACHE)
    unset(TBB_INSTALL_DIR CACHE)
    unset(CMAKE_PREFIX_PATH CACHE)
    set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    message("TBB_INSTALL_DIR = ${TBB_INSTALL_DIR}")
    find_package(TBB)
endif()


message("TBB_LIBRARIES = ${TBB_LIBRARIES}")

message("Build system will compile libgff")
message("==================================================================")
ExternalProject_Add(libgff
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/Kingsford-Group/libgff/archive/v1.0.tar.gz -o libgff.tgz &&
    tar -xzvf libgff.tgz &&
    rm -fr libgff &&
    mv libgff-1.0 libgff
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff
    UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff/build
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_SOURCE_DIR}/external/install
)

# Because of the way that Apple has changed SIP
# in el capitan, some headers may be in a new location
if (APPLE)
    set(STADEN_INC "-I/usr/local/include")
    set(STADEN_LIB "-L/usr/local/lib")
endif()

if (CONDA_BUILD)
  set(LZFLAG "-lz")
else ()
  set(LZFLAG "")
endif (CONDA_BUILD)

message("Build system will compile Staden IOLib")
message("==================================================================")
ExternalProject_Add(libstadenio
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/COMBINE-lab/staden-io_lib/archive/v1.14.8.tar.gz -o staden-io_lib-v1.14.8.tar.gz &&
                     mkdir -p staden-io_lib-1.14.8 &&
                     tar -xzf staden-io_lib-v1.14.8.tar.gz --strip-components=1 -C staden-io_lib-1.14.8 &&
                     rm -fr staden-io_lib &&
                     mv -f staden-io_lib-1.14.8 staden-io_lib
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/staden-io_lib
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ./configure --enable-shared=no --without-libcurl --prefix=<INSTALL_DIR> LDFLAGS=${LIBSTADEN_LDFLAGS} CFLAGS=${LIBSTADEN_CFLAGS} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_COMMAND make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS+=${STADEN_INC} CFLAGS+=${STADEN_LIB} CFLAGS+=${LZFLAG}

    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
)

message("Build system will fetch SPDLOG")
message("==================================================================")
ExternalProject_Add(libspdlog
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/COMBINE-lab/spdlog/archive/v1.12.tar.gz -o spdlog-v1.12.tar.gz &&
                     tar -xzf spdlog-v1.12.tar.gz &&
                     rm -fr spdlog &&
                     mv -f  spdlog-1.12 spdlog
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/spdlog
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND cp -r <SOURCE_DIR>/include/spdlog <INSTALL_DIR>/include
)

set (FAST_MALLOC_LIB "")
set (HAVE_FAST_MALLOC FALSE)

# See if we have Jemalloc
find_package(Jemalloc)
if (Jemalloc_FOUND)
    message("Found Jemalloc library --- using this memory allocator")
    set (FAST_MALLOC_LIB ${JEMALLOC_LIBRARIES})
    set (HAVE_FAST_MALLOC TRUE)
endif()

if (NOT HAVE_FAST_MALLOC)
    # See if we have Tcmalloc
    find_package(Tcmalloc)
    if (Tcmalloc_FOUND)
        message("Fount TCMalloc library --- using this memory allocator")
        set (TCMALLOC_LIB ${Tcmalloc_LIBRARIES})
        set (FAST_MALLOC_LIB ${TCMALLOC_LIB})
        set (HAVE_FAST_MALLOC TRUE)
    endif()
endif()

if(CONDA_BUILD)
  set (JEMALLOC_FLAGS "CC=${CMAKE_C_COMPILER} CFLAGS=-fPIC CPPFLAGS=-fPIC")
else ()
  set (JEMALLOC_FLAGS "CC=${CMAKE_C_COMPILER}")
endif()

if (NOT HAVE_FAST_MALLOC)
    message("Build system will fetch and use JEMalloc")
    message("==================================================================")
    ExternalProject_Add(libjemalloc
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
        DOWNLOAD_COMMAND curl -k -L https://github.com/COMBINE-lab/jemalloc/archive/4.0.4.tar.gz -o jemalloc-4.0.4.tar.gz &&
        tar -xzf jemalloc-4.0.4.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/jemalloc-4.0.4
        BUILD_IN_SOURCE TRUE
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        CONFIGURE_COMMAND sh -c "${JEMALLOC_FLAGS} ./autogen.sh --prefix=<INSTALL_DIR>"
        INSTALL_COMMAND cp -r lib <INSTALL_DIR>/ && cp -r include <INSTALL_DIR>/
        )

    set (FAST_MALLOC_LIB ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib/libjemalloc.a)
    set (HAVE_FAST_MALLOC TRUE)
endif ()


##
## This depenency is for RapMap
##
#message("Build system will fetch and build SparseHash")
#message("==================================================================")
#ExternalProject_Add(libsparsehash
#    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
#    DOWNLOAD_COMMAND curl -k -L https://github.com/COMBINE-lab/sparsehash/archive/sparsehash-2.0.2.tar.gz -o sparsehash-2.0.2.tar.gz &&
#        tar -xzf sparsehash-2.0.2.tar.gz
#    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/sparsehash-sparsehash-2.0.2
#    BUILD_IN_SOURCE TRUE
#    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
#    CONFIGURE_COMMAND sh -c "CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ./configure --prefix=<INSTALL_DIR>"
#    INSTALL_COMMAND make install
#)

###
#
# Done building external dependencies.
#
###

set (CPACK_SOURCE_IGNORE_FILES
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

# Recurse into Salmon source directory
add_subdirectory ( src )
#add_dependencies(salmon RapMap)

# build a CPack driven installer package
include (CPack)

set(ARCHIVE_NAME ${CMAKE_PROJECT_NAME}-${PROJECT_VERSION})
add_custom_target(dist
    COMMAND git archive --prefix=${ARCHIVE_NAME}/ HEAD
        | gzip > ${CMAKE_BINARY_DIR}/${ARCHIVE_NAME}.tar.gz
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
