include_directories(
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/include/eigen3
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/cereal/include
${GAT_SOURCE_DIR}/external/install/include
${GAT_SOURCE_DIR}/external/install/include/jellyfish-2.2.6
${GAT_SOURCE_DIR}/external/install/include/bwa
${ZLIB_INCLUDE_DIR}
${TBB_INCLUDE_DIRS}
${Boost_INCLUDE_DIRS}
${GAT_SOURCE_DIR}/external/install/include/rapmap
)

if (JELLYFISH_FOUND)
    include_directories(${JELLYFISH_INCLUDE_DIR})
else()
    include_directories(${GAT_SOURCE_DIR}/external/install/include/jellyfish-2.2.6)
endif()

set ( SALMON_MAIN_SRCS
QSufSort.c
is.c
bwt_gen.c
bwtindex.c
xxhash.c
${GAT_SOURCE_DIR}/external/install/src/rapmap/RapMapFileSystem.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/RapMapSAIndexer.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/RapMapSAIndex.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/RapMapUtils.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/HitManager.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/rank9b.cpp
${GAT_SOURCE_DIR}/external/install/src/rapmap/bit_array.c
CollapsedEMOptimizer.cpp
CollapsedGibbsSampler.cpp
Salmon.cpp
BuildSalmonIndex.cpp
SalmonQuantify.cpp
FragmentLengthDistribution.cpp
FragmentStartPositionDistribution.cpp
SequenceBiasModel.cpp
TranscriptGroup.cpp
GZipWriter.cpp
#${GAT_SOURCE_DIR}/external/install/src/rapmap/sais.c
)

set (SALMON_ALIGN_SRCS
FASTAParser.cpp
ErrorModel.cpp
AlignmentModel.cpp
FragmentLengthDistribution.cpp
SalmonQuantifyAlignments.cpp
)

set (SALMON_LIB_SRCS
BWAUtils.cpp
LibraryFormat.cpp
GenomicFeature.cpp
VersionChecker.cpp
SBModel.cpp
FastxParser.cpp
StadenUtils.cpp
SalmonUtils.cpp
DistributionUtils.cpp
SalmonStringUtils.cpp
SimplePosBias.cpp
SGSmooth.cpp
)

set ( UNIT_TESTS_SRCS
    ${GAT_SOURCE_DIR}/tests/UnitTests.cpp
    FragmentLengthDistribution.cpp
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
message("Boost_LIBRARIES = ${Boost_LIBRARIES}")

# Set the RPATH
if (NOT APPLE)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:$ORIGIN/:$ORIGIN/../../external/install/lib")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
else()
  # use, i.e. don't skip the full RPATH for the build tree
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

# Build the Salmon library
add_library(salmon_core STATIC ${SALMON_LIB_SRCS} )


# Build the salmon executable
add_executable(salmon ${SALMON_MAIN_SRCS} ${SALMON_ALIGN_SRCS})

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
    salmon_core
    gff
    ${PTHREAD_LIB}
    ${Boost_LIBRARIES}
    ${GAT_SOURCE_DIR}/external/install/lib/libstaden-read.a
    ${ZLIB_LIBRARY}
    ${SUFFARRAY_LIB}
    ${SUFFARRAY64_LIB}
    ${GAT_SOURCE_DIR}/external/install/lib/libjellyfish-2.0.a
    ${GAT_SOURCE_DIR}/external/install/lib/libbwa.a
    m
    ${LIBLZMA_LIBRARIES}
    ${BZIP2_LIBRARIES}
    ${TBB_LIBRARIES}
    ${LIBSALMON_LINKER_FLAGS}
    ${NON_APPLECLANG_LIBS}
    ${FAST_MALLOC_LIB}
    ${LIBRT}
)

##
# External dependencies of salmon_core and salmon
##
if (${FETCHED_JEMALLOC})
   add_dependencies(salmon_core libjemalloc)
   add_dependencies(salmon libjemalloc)
endif()

if (${FETCHED_BOOST})
   add_dependencies(salmon_core libboost)
   add_dependencies(salmon libboost)
endif()

if (${FETCHED_JELLYFISH})
   add_dependencies(salmon_core libjellyfish)
   add_dependencies(salmon libjellyfish)
endif()

if (${FETCHED_TBB})
   add_dependencies(salmon_core libtbb)
   add_dependencies(salmon libtbb)
endif()

add_dependencies(salmon_core libcereal)
add_dependencies(salmon libcereal)
add_dependencies(salmon_core libstadenio)
add_dependencies(salmon libstadenio)
add_dependencies(salmon_core libspdlog)
add_dependencies(salmon libspdlog)
add_dependencies(salmon_core libbwa)
add_dependencies(salmon libbwa)
add_dependencies(salmon_core libgff)
add_dependencies(salmon libgff)
add_dependencies(salmon_core libdivsufsort)
add_dependencies(salmon libdivsufsort)

# Link the executable
target_link_libraries(unitTests
    salmon_core
    gff
    ${PTHREAD_LIB}
    ${Boost_LIBRARIES}
    ${GAT_SOURCE_DIR}/external/install/lib/libstaden-read.a
    ${ZLIB_LIBRARY}
    ${SUFFARRAY_LIB}
    ${SUFFARRAY64_LIB}
    ${GAT_SOURCE_DIR}/external/install/lib/libjellyfish-2.0.a
    ${GAT_SOURCE_DIR}/external/install/lib/libbwa.a
    m
    ${LIBLZMA_LIBRARIES}
    ${BZIP2_LIBRARIES}
    ${TBB_LIBRARIES}
    ${LIBSALMON_LINKER_FLAGS}
    ${NON_APPLECLANG_LIBS}
    ${FAST_MALLOC_LIB}
    ${LIBRT}
    )

add_dependencies(salmon unitTests)

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
add_test( NAME salmon_read_test_fmd COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/TestSalmonFMD.cmake )
add_test( NAME salmon_read_test_quasi COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/TestSalmonQuasi.cmake )

####
#
# Deprecated or currently unused
#
####


# # use, i.e. don't skip the full RPATH for the build tree
# SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# # when building, don't use the install RPATH already
# # (but later on when installing)
# SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# MESSAGE( "LINK DIRECTORIES : ${CMAKE_LIBRARY_PATH}")
# SET(CMAKE_INSTALL_RPATH "${CMAKE_LIBRARY_PATH}")

# # add the automatically determined parts of the RPATH
# # which point to directories outside the build tree to the install RPATH
# SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# # the RPATH to be used when installing, but only if it's not a system directory
# LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
# IF("${isSystemDir}" STREQUAL "-1")
#    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# ENDIF("${isSystemDir}" STREQUAL "-1")


# set( BUILD_LUT_SRCS
# BuildLUT.cpp
# IndexedCounter.cpp
# PerfectHashIndexer.cpp
# )

# add_executable( BuildLUT ${BUILD_LUT_SRCS} )
# target_link_libraries( BuildLUT
# 	${Boost_LIBRARIES} ${ZLIB_LIBRARY}
# 	cmph # perfect hashing library
# 	jellyfish-2.0 pthread
#     gomp lib_activeobject lib_g2logger m tbb)

# INSTALL(TARGETS BuildLUT DESTINATION ${PROJECT_SOURCE_DIR}/bin COMPONENT comp_buildlut)

# set ( BUILD_INDEX_SRCS Indexer.cpp )

# set( ANALYZE_TRANSCRIPT_GRAPH
# AnalyzeTranscriptGraph.cpp
# IndexedCounter.cpp
# PerfectHashIndexer.cpp
# )

# set ( PROBES_TO_GENES_SRC MapProbesToGenes.cpp affymetrix_utils.cpp )
# set( COMPUTE_TGRAPH_SRCS
# ComputeTranscriptGraph.cpp
# IndexedCounter.cpp
# PerfectHashIndexer.cpp
# )

# add_executable( ComputeTranscriptGraph ${COMPUTE_TGRAPH_SRCS} )
# target_link_libraries( ComputeTranscriptGraph
# 	${Boost_LIBRARIES} ${ZLIB_LIBRARY}
# 	cmph # perfect hashing library
# 	jellyfish pthread
#     gomp lib_activeobject lib_g2logger blas m tbb nnls)


# add_executable( AnalyzeTranscriptGraph ${ANALYZE_TRANSCRIPT_GRAPH} )
# target_link_libraries( AnalyzeTranscriptGraph
# 	${Boost_LIBRARIES} ${ZLIB_LIBRARY}
# 	cmph # perfect hashing library
# 	jellyfish pthread
#     gomp lib_activeobject lib_g2logger blas m tbb)


#add_executable( BuildIndex ${BUILD_INDEX_SRCS} )
#target_link_libraries( BuildIndex ${Boost_LIBRARIES} jellyfish pthread m tbb)


#add_executable( TestTsnnls ${TEST_SRCS} )
#target_link_libraries( TestTsnnls ${Boost_LIBRARIES} ${ZLIB_LIBRARY} blas lapack gomp lib_activeobject lib_g2logger tsnnls )
