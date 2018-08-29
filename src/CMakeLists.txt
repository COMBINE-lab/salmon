include_directories(
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/include/eigen3
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/cereal/include
${GAT_SOURCE_DIR}/external/install/include
#${GAT_SOURCE_DIR}/external/install/include/bwa
${ZLIB_INCLUDE_DIR}
${TBB_INCLUDE_DIRS}
${Boost_INCLUDE_DIRS}
${GAT_SOURCE_DIR}/external/install/include/rapmap
${GAT_SOURCE_DIR}/external/install/include/rapmap/digestpp
${ICU_INC_DIRS}
)

set ( SALMON_MAIN_SRCS
${GAT_SOURCE_DIR}/external/install/src/rapmap/RapMapFileSystem.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/RapMapSAIndexer.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/RapMapSAIndex.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/RapMapUtils.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/HitManager.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/rank9b.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/bit_array.c
EMUtils.cpp
CollapsedEMOptimizer.cpp
CollapsedCellOptimizer.cpp
CollapsedGibbsSampler.cpp
Salmon.cpp
BuildSalmonIndex.cpp
Dedup.cpp
Alevin.cpp
SalmonAlevin.cpp
WhiteList.cpp
SalmonQuantify.cpp
FragmentLengthDistribution.cpp
FragmentStartPositionDistribution.cpp
SequenceBiasModel.cpp
GZipWriter.cpp
SalmonQuantMerge.cpp
ProgramOptionsGenerator.cpp
)

set (SALMON_ALIGN_SRCS
FASTAParser.cpp
AlignmentModel.cpp
FragmentLengthDistribution.cpp
SalmonQuantifyAlignments.cpp
)

set (ALEVIN_LIB_SRCS
  edlib.cpp
  AlevinUtils.cpp
)

set (SALMON_LIB_SRCS
${GAT_SOURCE_DIR}/src/jellyfish/mer_dna.cc
xxhash.c
TranscriptGroup.cpp
EffectiveLengthStats.cpp
LibraryFormat.cpp
GenomicFeature.cpp
VersionChecker.cpp
SBModel.cpp
FastxParser.cpp
StadenUtils.cpp
SalmonUtils.cpp
DistributionUtils.cpp
SalmonExceptions.cpp
SalmonStringUtils.cpp
SimplePosBias.cpp
SGSmooth.cpp
${GAT_SOURCE_DIR}/src/metro/metrohash64.cpp
)

set (KSW2PP_BASIC_LIB_SRCS
ksw2pp/kalloc.c
ksw2pp/ksw2_extd.c
ksw2pp/ksw2_extz.c
ksw2pp/ksw2_gg.c
ksw2pp/ksw2_gg2.c
ksw2pp/ksw2_gg2_sse.c
ksw2pp/KSW2Aligner.cpp
)

set (KSW2PP_ADVANCED_LIB_SRCS
ksw2pp/ksw2_extd2_sse.c
ksw2pp/ksw2_extf2_sse.c
ksw2pp/ksw2_extz2_sse.c
)

# check if we know how to do IPO
check_ipo_supported(RESULT HAS_IPOHAS_IPO)

add_library(ksw2pp_sse2 OBJECT ${KSW2PP_ADVANCED_LIB_SRCS})
add_library(ksw2pp_sse4 OBJECT ${KSW2PP_ADVANCED_LIB_SRCS})
add_library(ksw2pp_basic OBJECT ${KSW2PP_BASIC_LIB_SRCS})

set_target_properties(ksw2pp_sse2 PROPERTIES COMPILE_FLAGS "-O3 -msse2 -mno-sse4.1")
set_target_properties(ksw2pp_sse2 PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;KSW_SSE2_ONLY;HAVE_KALLOC")
set_target_properties(ksw2pp_sse4 PROPERTIES COMPILE_FLAGS "-O3 -msse4.1")
set_target_properties(ksw2pp_sse4 PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")
set_target_properties(ksw2pp_basic PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")

set_target_properties(ksw2pp_basic PROPERTIES INCLUDE_DIRECTORIES ${GAT_SOURCE_DIR}/include)
set_target_properties(ksw2pp_sse4 PROPERTIES INCLUDE_DIRECTORIES ${GAT_SOURCE_DIR}/include)

# Build the ksw2pp library
add_library(ksw2pp STATIC $<TARGET_OBJECTS:ksw2pp_sse2> $<TARGET_OBJECTS:ksw2pp_sse4> $<TARGET_OBJECTS:ksw2pp_basic>)
set_target_properties(ksw2pp PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")
if(NOT CONDA_BUILD AND HAS_IPO)
  set_property(TARGET ksw2pp PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

set ( UNIT_TESTS_SRCS
    ${GAT_SOURCE_DIR}/tests/UnitTests.cpp
    FragmentLengthDistribution.cpp
    ${GAT_SOURCE_DIR}/external/install/src/rapmap/rank9b.cpp
    ${GAT_SOURCE_DIR}/external/install/src/rapmap/bit_array.c
)


link_directories(
${GAT_SOURCE_DIR}/lib
${GAT_SOURCE_DIR}/external/install/lib
${Boost_LIBRARY_DIRS}
${TBB_LIBRARY_DIRS}
${LAPACK_LIBRARY_DIR}
${BLAS_LIBRARY_DIR}
)

message("TBB_LIBRARIES = ${TBB_LIBRARIES}")

# Set the RPATH
if (NOT APPLE)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:$ORIGIN/:$ORIGIN/../../external/install/lib")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
else()
  # use, i.e. do not skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

  # the RPATH to be used when installing
  set(CMAKE_INSTALL_RPATH "")

  # don't add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
endif()


set (TGT_RELEASE_FLAGS "${TGT_COMPILE_FLAGS};${TGT_WARN_FLAGS}")
set (TGT_DEBUG_FLAGS "-g;${TGT_COMPILE_FLAGS};${TGT_WARN_FLAGS}")

# Build the Salmon library
add_library(salmon_core STATIC ${SALMON_LIB_SRCS} )
target_compile_definitions(salmon_core PUBLIC
RAPMAP_SALMON_SUPPORT=1
HAVE_ANSI_TERM=1
HAVE_SSTREAM=1
)
target_compile_options(salmon_core PUBLIC "$<$<CONFIG:DEBUG>:${TGT_DEBUG_FLAGS}>")
target_compile_options(salmon_core PUBLIC "$<$<CONFIG:RELEASE>:${TGT_RELEASE_FLAGS}>")
if(NOT CONDA_BUILD AND HAS_IPO)
  set_property(TARGET salmon_core PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

# Build the Alevin library
add_library(alevin_core STATIC ${ALEVIN_LIB_SRCS})
target_compile_definitions(alevin_core PUBLIC
RAPMAP_SALMON_SUPPORT=1
HAVE_ANSI_TERM=1
HAVE_SSTREAM=1
span_FEATURE_MAKE_SPAN_TO_STD=14
)
target_compile_options(alevin_core PUBLIC "$<$<CONFIG:DEBUG>:${TGT_DEBUG_FLAGS}>")
target_compile_options(alevin_core PUBLIC "$<$<CONFIG:RELEASE>:${TGT_RELEASE_FLAGS}>")

if(NOT CONDA_BUILD AND HAS_IPO)
  set_property(TARGET alevin_core PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

# Build the salmon executable
add_executable(salmon ${SALMON_MAIN_SRCS} ${SALMON_ALIGN_SRCS})

if(NOT CONDA_BUILD AND HAS_IPO)
  set_property(TARGET salmon PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

add_executable(unitTests ${UNIT_TESTS_SRCS})

#add_executable(salmon-read ${SALMON_READ_SRCS})
#set_target_properties(salmon-read PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_LIBPTHREAD -D_PBGZF_USE -fopenmp"
#    LINK_FLAGS "-DHAVE_LIBPTHREAD -D_PBGZF_USE -fopenmp")

#set_target_properties(salmon_core salmon PROPERTIES LINK_SEARCH_END_STATIC TRUE)

# our suffix array construction libraries
set (SUFFARRAY_LIB ${GAT_SOURCE_DIR}/external/install/lib/libdivsufsort.a)
set (SUFFARRAY64_LIB ${GAT_SOURCE_DIR}/external/install/lib/libdivsufsort64.a)


# Link the executable
target_link_libraries(salmon
    Threads::Threads
    salmon_core
    alevin_core
    gff
    ${Boost_LIBRARIES}
    ${ICU_LIBS}
    ${GAT_SOURCE_DIR}/external/install/lib/libstaden-read.a
    ${ZLIB_LIBRARY}
    ${SUFFARRAY_LIB}
    ${SUFFARRAY64_LIB}
    #${GAT_SOURCE_DIR}/external/install/lib/libbwa.a
    m
    ${LIBLZMA_LIBRARIES}
    ${BZIP2_LIBRARIES}
    ${TBB_LIBRARIES}
    ${LIBSALMON_LINKER_FLAGS}
    ${NON_APPLECLANG_LIBS}
    ${FAST_MALLOC_LIB}
    ${LIBRT}
    ksw2pp
    ${CMAKE_DL_LIBS}
)

# Link the executable
target_link_libraries(unitTests
    Threads::Threads
    salmon_core
    alevin_core
    gff
    ${Boost_LIBRARIES}
    ${ICU_LIBS}
    ${GAT_SOURCE_DIR}/external/install/lib/libstaden-read.a
    ${ZLIB_LIBRARY}
    ${SUFFARRAY_LIB}
    ${SUFFARRAY64_LIB}
    #${GAT_SOURCE_DIR}/external/install/lib/libbwa.a
    m
    ${LIBLZMA_LIBRARIES}
    ${BZIP2_LIBRARIES}
    ${TBB_LIBRARIES}
    ${LIBSALMON_LINKER_FLAGS}
    ${NON_APPLECLANG_LIBS}
    ${LIBRT}
    ${CMAKE_DL_LIBS}
    )

add_dependencies(salmon unitTests)

##
# External dependencies of salmon_core and salmon
##
if (${FETCHED_JEMALLOC})
   add_dependencies(salmon_core libjemalloc)
   add_dependencies(salmon libjemalloc)
endif()

if (${FETCHED_BOOST})
   add_dependencies(alevin_core libboost)
   add_dependencies(salmon_core libboost)
   add_dependencies(salmon libboost)
endif()

if (${FETCHED_TBB})
   add_dependencies(salmon_core libtbb)
   add_dependencies(alevin_core libtbb)
   add_dependencies(salmon libtbb)
endif()

add_dependencies(salmon_core libcereal)
add_dependencies(salmon libcereal)
add_dependencies(alevin_core libstadenio)
add_dependencies(salmon_core libstadenio)
add_dependencies(salmon libstadenio)
#add_dependencies(salmon_core libbwa)
#add_dependencies(salmon libbwa)
add_dependencies(salmon_core libgff)
add_dependencies(salmon libgff)
add_dependencies(salmon_core libdivsufsort)
add_dependencies(salmon libdivsufsort)

### No need for this, I think
##  This ensures that the salmon executable should work with or without `make install`
###
## Grumble grumble . . . OSX
#if (APPLE)
#    # only attempt install_name_tool for tbb if we installed it
#    if (${TBB_LIBRARY_DIRS} MATCHES ${GAT_SOURCE_DIR}/external/install/lib)
#        add_custom_command(TARGET salmon
#            POST_BUILD
#            COMMAND install_name_tool -change libtbb.dylib @rpath/libtbb.dylib ${GAT_SOURCE_DIR}/build/src/salmon
#            COMMAND install_name_tool -change libtbbmalloc.dylib @rpath/libtbbmalloc.dylib ${GAT_SOURCE_DIR}/build/src/salmon
#            COMMAND install_name_tool -change libtbbmalloc_proxy.dylib @rpath/libtbbmalloc_proxy.dylib ${GAT_SOURCE_DIR}/build/src/salmon
#            COMMAND install_name_tool -add_rpath  ${GAT_SOURCE_DIR}/external/install/lib ${GAT_SOURCE_DIR}/build/src/salmon
#            )
#        add_custom_command(TARGET unitTests
#            POST_BUILD
#            COMMAND install_name_tool -change libtbb.dylib @rpath/libtbb.dylib ${GAT_SOURCE_DIR}/build/src/unitTests
#            COMMAND install_name_tool -change libtbbmalloc.dylib @rpath/libtbbmalloc.dylib ${GAT_SOURCE_DIR}/build/src/unitTests
#            COMMAND install_name_tool -change libtbbmalloc_proxy.dylib @rpath/libtbbmalloc_proxy.dylib ${GAT_SOURCE_DIR}/build/src/unitTests
#            COMMAND install_name_tool -add_rpath  ${GAT_SOURCE_DIR}/external/install/lib ${GAT_SOURCE_DIR}/build/src/unitTests
#            )
#    endif()
#else()
#    # related to complete static linking --- on hold    
#    set (BOOST_THREAD_LIBRARY)
#endif()

#if (APPLE)
#	add_custom_command(TARGET salmon
#		POST_BUILD
#		COMMAND install_name_tool -add_rpath ${GAT_SOURCE_DIR}/external/install/lib salmon
#	COMMAND install_name_tool -add_rpath @executable_path/../lib salmon
#		)
#endif()

##### ======================================

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
    "${GAT_SOURCE_DIR}" CACHE PATH "Default install prefix" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(INSTALL_LIB_DIR lib )
set(INSTALL_BIN_DIR bin )
set(INSTALL_INCLUDE_DIR include )

install(DIRECTORY
        ${GAT_SOURCE_DIR}/external/install/lib/
        DESTINATION ${INSTALL_LIB_DIR}
	    FILES_MATCHING PATTERN "libtbb*.${SHARED_LIB_EXTENSION}*"
    )

# install(FILES ${Boost_LIBRARIES}
# 	           DESTINATION ${INSTALL_LIB_DIR})

install(TARGETS salmon salmon_core
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
        )

add_custom_command(TARGET unitTests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:unitTests> ${GAT_SOURCE_DIR}/tests/$<TARGET_FILE_NAME:unitTests>
    COMMENT "Copying unitTests"
)

set(POST_INSTALL_SCRIPT ${GAT_SOURCE_DIR}/cmake/PostInstall.cmake)

install(
    CODE
    "
    execute_process(COMMAND \"${CMAKE_COMMAND}\"
                            -DCMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}
                            -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                            -P \"${POST_INSTALL_SCRIPT}\")
    "
)

include(InstallRequiredSystemLibraries)
add_test( NAME unit_tests COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${CMAKE_INSTALL_PREFIX} -P ${GAT_SOURCE_DIR}/cmake/UnitTests.cmake )
add_test( NAME salmon_read_test_quasi COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/TestSalmonQuasi.cmake )

# Remove this test since we are removing support for the FMD index. 
# add_test( NAME salmon_read_test_fmd COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/TestSalmonFMD.cmake )
