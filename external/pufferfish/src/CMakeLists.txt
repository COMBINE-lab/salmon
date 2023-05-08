# check if we know how to do IPO
include(CheckIPOSupported)
check_ipo_supported(RESULT HAS_IPO)

link_directories(
		${GAT_SOURCE_DIR}/external/install/lib
    ${GAT_SOURCE_DIR}/external/install/lib64
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

add_library(ksw2pp_sse2 OBJECT ${KSW2PP_ADVANCED_LIB_SRCS})
add_library(ksw2pp_sse4 OBJECT ${KSW2PP_ADVANCED_LIB_SRCS})
add_library(ksw2pp_basic OBJECT ${KSW2PP_BASIC_LIB_SRCS})

if (USE_ARM)
    set_target_properties(ksw2pp_sse2 PROPERTIES COMPILE_FLAGS "-O3")
    set_target_properties(ksw2pp_sse2 PROPERTIES COMPILE_DEFINITIONS "__SSE2__;KSW_CPU_DISPATCH;KSW_USE_ARM;KSW_SSE2_ONLY;HAVE_KALLOC")
    set_target_properties(ksw2pp_sse4 PROPERTIES COMPILE_FLAGS "-O3")
    set_target_properties(ksw2pp_sse4 PROPERTIES COMPILE_DEFINITIONS "__SSE2__;__SSE4_1__;KSW_CPU_DISPATCH;KSW_USE_ARM;HAVE_KALLOC")
    set_target_properties(ksw2pp_basic PROPERTIES COMPILE_DEFINITIONS "__SSE2__;KSW_CPU_DISPATCH;KSW_USE_ARM;HAVE_KALLOC")
else()
    set_target_properties(ksw2pp_sse2 PROPERTIES COMPILE_FLAGS "-O3 -msse -msse2 -mno-sse4.1")
    set_target_properties(ksw2pp_sse2 PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;KSW_SSE2_ONLY;HAVE_KALLOC")
    set_target_properties(ksw2pp_sse4 PROPERTIES COMPILE_FLAGS "-O3 -msse -msse2 -msse3 -mssse3 -msse4 -msse4.1")
    set_target_properties(ksw2pp_sse4 PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")
    set_target_properties(ksw2pp_basic PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")
endif()

set_target_properties(ksw2pp_basic PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include)
set_target_properties(ksw2pp_sse4 PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include) 
set_target_properties(ksw2pp_sse2 PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include) 

# Build the ksw2pp library
add_library(ksw2pp STATIC $<TARGET_OBJECTS:ksw2pp_sse2> $<TARGET_OBJECTS:ksw2pp_sse4> $<TARGET_OBJECTS:ksw2pp_basic>)
set_target_properties(ksw2pp PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")
if(HAS_IPO)
  #set_property(TARGET ksw2pp PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

#set_target_properties(ksw2pp PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${GAT_SOURCE_DIR}/external/install/lib)

set(pufferfish_lib_srcs 
    PufferfishIndexer.cpp 
    PufferfishBaseIndex.cpp
    PufferfishIndex.cpp 
    PufferfishSparseIndex.cpp 
    PufferfishLossyIndex.cpp
    edlib.cpp
    Util.cpp
		rank9sel.cpp
    rank9b.cpp
    PufferfishValidate.cpp
	  PufferfishStats.cpp
    PufferfishTestLookup.cpp 
    PufferfishExamine.cpp
    PufferfishKmerQuery.cpp
    FastxParser.cpp 
	  PufferfishBinaryGFAReader.cpp
    PufferFS.cpp
    xxhash.c
    FixFasta.cpp
    MemCollector.cpp
	  MemChainer.cpp
		PuffAligner.cpp
	  PufferfishAligner.cpp
	  RefSeqConstructor.cpp
	  metro/metrohash64.cpp
)
if (BUILD_PUFF_FOR_SALMON)
 list(APPEND pufferfish_lib_srcs "${GAT_SOURCE_DIR}/src/LibraryFormat.cpp")
endif()

message("add library puffer : ${pufferfish_lib_srcs}")
add_library(puffer STATIC ${pufferfish_lib_srcs})
target_compile_options(puffer PUBLIC $<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>)
target_compile_options(puffer PUBLIC $<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>)
if (BUILD_PUFF_FOR_SALMON)
   target_include_directories(puffer PUBLIC "${GAT_SOURCE_DIR}/include")
endif()
if (FETCHED_TBB)
add_dependencies(puffer libtbb)
else()
add_dependencies(puffer TBB::tbb)
add_dependencies(puffer TBB::tbbmalloc)
endif()

if (NOT BUILD_PUFF_FOR_SALMON)
message("Building extra components")

# The actual executable need not be built for salmon
# Also, the `-static-libgcc` flag is a problem under OSX/clang.
add_executable(pufferfish Pufferfish.cpp)
target_compile_options(pufferfish PUBLIC $<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>)
target_compile_options(pufferfish PUBLIC $<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>)

if (FETCHED_TBB)
add_dependencies(pufferfish libtbb)
else()
add_dependencies(pufferfish TBB::tbb)
add_dependencies(pufferfish TBB::tbbmalloc)
endif()



#if (BUILD_PUFF_FOR_SALMON)
#   target_include_directories(pufferfish PUBLIC "${GAT_SOURCE_DIR}/include")
#endif()


add_dependencies(pufferfish puffer)
add_dependencies(pufferfish ksw2pp) 
add_dependencies(pufferfish twopaco)
add_dependencies(pufferfish graphdump)
add_dependencies(pufferfish ntcard)

if(TBB_RECONFIGURE OR TBB_TARGET_EXISTED)
  target_link_libraries(pufferfish 
    ksw2pp
    puffer 
    z 
    twopaco 
    graphdump 
    ntcard 
    ${LIBRT} 
    ${ASAN_LIB} 
    ${JEMALLOC_LIBRARIES} 
    ksw2pp
    puffer 
    ${CMAKE_DL_LIBS} 
    TBB::tbb
    TBB::tbbmalloc
    Threads::Threads 
  )
endif()

add_executable(bcalm_pufferize BCALMPufferizer.cpp FastxParser.cpp)
target_compile_options(bcalm_pufferize PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(bcalm_pufferize PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

#add_executable(edgedensity EdgeDensity.cpp)
#target_compile_options(edgedensity PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
#target_compile_options(edgedensity PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
#add_dependencies(edgedensity puffer)
#

#add_executable(edgedensity2 EdgeDensityCustom.cpp)
#target_compile_options(edgedensity2 PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
#target_compile_options(edgedensity2 PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
#add_dependencies(edgedensity2 puffer)

add_executable(filtersam testing/SAMFilter.cpp)
target_include_directories(filtersam PRIVATE
    ${GAT_SOURCE_DIR}/external/install/include
    ${GAT_SOURCE_DIR}/external/install/include/htslib
    )
target_compile_options(filtersam PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(filtersam PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
add_dependencies(filtersam libseqlib)
#add_dependencies(filtersam puffer)


#add_executable(krakmap KrakMap.cpp cedar/Taxa.cpp)
#target_compile_options(krakmap PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
#target_compile_options(krakmap PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

add_executable(cedar cedar/TargetGroup.cpp cedar/Cedar.cpp cedar/Taxa.cpp)

target_include_directories(cedar PRIVATE
    ${GAT_SOURCE_DIR}/external/install/include
    ${GAT_SOURCE_DIR}/external/install/include/htslib
)
target_compile_options(cedar PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(cedar PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
target_link_libraries(cedar Async++ Threads::Threads seqlib hts puffer bz2 lzma z SetCover ${ASAN_LIB})
add_dependencies(cedar libseqlib)
add_dependencies(cedar libsetcover)
#add_dependencies(ceder Async++)
#add_dependencies(ceder seqlib)

#get_target_property(dirs filtersam INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()

#add_executable(readKrk readKrk.cpp taxa.cpp)


#[[
#add_executable(rank_test rank_test.cpp rank9b.cpp rank9sel.cpp)
#target_compile_options(rank_test PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
#target_compile_options(rank_test PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
#target_link_libraries(rank_test)
#]]

add_executable(getLineage constructLineage.cpp)
target_compile_options(getLineage PUBLIC $<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>)
target_compile_options(getLineage PUBLIC $<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>)


target_link_libraries(bcalm_pufferize Threads::Threads z ${ASAN_LIB})
#target_link_libraries(fixFasta Threads::Threads z)


#if (${TBB_TARGET_EXISTED})
#target_link_libraries(edgedensity Threads::Threads z puffer ${ASAN_LIB})
#target_link_libraries(edgedensity2 Threads::Threads z puffer ${ASAN_LIB})
target_link_libraries(filtersam Threads::Threads seqlib hts puffer bz2 lzma z ${ASAN_LIB})
add_dependencies(filtersam libseqlib)
#endif()


endif()
