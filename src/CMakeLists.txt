set ( SALMON_MAIN_SRCS
QSufSort.c
is.c
bwt_gen.c
bwtindex.c
CollapsedEMOptimizer.cpp
Salmon.cpp
BuildSalmonIndex.cpp
SalmonQuantify.cpp
FragmentLengthDistribution.cpp 
FragmentStartPositionDistribution.cpp
SequenceBiasModel.cpp
SalmonUtils.cpp
StadenUtils.cpp
TranscriptGroup.cpp
)

set (SALMON_ALIGN_SRCS
FASTAParser.cpp 
ErrorModel.cpp 
AlignmentModel.cpp
FragmentLengthDistribution.cpp 
SalmonQuantifyAlignments.cpp 
)

set (SALMON_LIB_SRCS
LibraryFormat.cpp
GenomicFeature.cpp
VersionChecker.cpp
SalmonUtils.cpp
SalmonStringUtils.cpp
ComputeBiasFeatures.cpp
PerformBiasCorrection.cpp
cokus.cpp
merge_files.cc
format.cc
)

include_directories( 
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/include/eigen3
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/cereal/include
${GAT_SOURCE_DIR}/external/install/include
${GAT_SOURCE_DIR}/external/install/include/jellyfish-2.1.3
${GAT_SOURCE_DIR}/external/install/include/bwa
${ZLIB_INCLUDE_DIR}
${TBB_INCLUDE_DIRS}
${Boost_INCLUDE_DIRS}
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
if (APPLE)
    ## This DOES NOT do what I / any one sane, expects.  Setting the 
    ## linker path on OSX is messed up.  Just tell the user to use 
    ## DYLD_FALLBACK_LIBRARY_PATH for now
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:$ORIGIN/:$ORIGIN/../../external/install/lib")
endif()

set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# Build the Salmon library
add_library(salmon_core STATIC ${SALMON_LIB_SRCS} )

# Build the salmon executable
add_executable(salmon ${SALMON_MAIN_SRCS} ${SALMON_ALIGN_SRCS})

#add_executable(salmon-read ${SALMON_READ_SRCS})
#set_target_properties(salmon-read PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_LIBPTHREAD -D_PBGZF_USE -fopenmp"
#    LINK_FLAGS "-DHAVE_LIBPTHREAD -D_PBGZF_USE -fopenmp")

#set_target_properties(salmon_core salmon PROPERTIES LINK_SEARCH_END_STATIC TRUE)

# Grumble grumble . . . OSX
if (APPLE)
    # only attempt install_name_tool for tbb if we installed it
    if (${TBB_LIBRARY_DIRS} MATCHES ${GAT_SOURCE_DIR}/external/install/lib)
        add_custom_command(TARGET salmon
            PRE_LINK
            COMMAND install_name_tool -id @rpath/libtbb.dylib ${TBB_LIBRARY_DIRS}/libtbb.dylib
            COMMAND install_name_tool -id @rpath/libtbbmalloc.dylib ${TBB_LIBRARY_DIRS}/libtbbmalloc.dylib
        )
    endif()

else()
    # related to complete static linking --- on hold	
    set (BOOST_THREAD_LIBRARY)
endif()

# Link the executable
target_link_libraries(salmon
    salmon_core 
    gff
    ${PTHREAD_LIB}
    ${Boost_LIBRARIES} 
    ${GAT_SOURCE_DIR}/external/install/lib/libstaden-read.a
    ${ZLIB_LIBRARY} 
    ${GAT_SOURCE_DIR}/external/install/lib/libjellyfish-2.0.a 
    ${GAT_SOURCE_DIR}/external/install/lib/libbwa.a
    m
    ${LIBLZMA_LIBRARIES}
    ${BZIP2_LIBRARIES}
    ${TBB_LIBRARIES}
    ${LIBSALMON_LINKER_FLAGS}
    ${NON_APPLECLANG_LIBS}
    ${FAST_MALLOC_LIB}
)

add_dependencies(salmon libbwa)

##
#  This ensures that the salmon executable should work with or without `make install`
##
if (APPLE)
	add_custom_command(TARGET salmon
		POST_BUILD
		COMMAND install_name_tool -add_rpath ${GAT_SOURCE_DIR}/external/install/lib salmon 
		)
endif()

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
if (APPLE)
install(DIRECTORY 
        ${GAT_SOURCE_DIR}/external/install/lib/ 
        DESTINATION ${INSTALL_LIB_DIR}
	    FILES_MATCHING PATTERN "libcmph*.${SHARED_LIB_EXTENSION}*"
       )
endif()


# install(FILES ${Boost_LIBRARIES}
# 	           DESTINATION ${INSTALL_LIB_DIR})

install(TARGETS salmon salmon_core
                RUNTIME DESTINATION bin 
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
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
add_test( NAME salmon_read_test COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/TestSalmon.cmake )

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

