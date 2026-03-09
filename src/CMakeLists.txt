if(${TBB_FOUND})
  get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
endif()

set (SALMON_ENTRY_SRCS
cli/Salmon.cpp
)

set (SALMON_CLI_SRCS
cli/VersionChecker.cpp
commands/AlevinRemoved.cpp
commands/SalmonQuantMerge.cpp
)

set (SALMON_INDEXING_SRCS
index/BuildSalmonIndex.cpp
)

set (SALMON_OPTIONS_SRCS
cli/ProgramOptionsGenerator.cpp
)

set (SALMON_OUTPUT_SRCS
output/GZipWriter.cpp
)

set (SALMON_MODEL_SRCS
model/FragmentLengthDistribution.cpp
model/FragmentStartPositionDistribution.cpp
model/EffectiveLengthStats.cpp
model/TranscriptGroup.cpp
model/SBModel.cpp
model/SimplePosBias.cpp
model/SGSmooth.cpp
model/SequenceBiasModel.cpp
)

set (SALMON_ALIGNMENT_SRCS
alignment/AlignmentModel.cpp
alignment/ONTAlignmentModel.cpp
alignment/AlignmentCommon.cpp
alignment/SalmonQuantifyAlignments.cpp
alignment/pipeline/AlignmentPipelineStages.cpp
alignment/BAMUtils.cpp
)

set (SALMON_INFERENCE_SRCS
inference/EMUtils.cpp
inference/CollapsedEMOptimizer.cpp
inference/CollapsedGibbsSampler.cpp
)

set (SALMON_QUANT_SRCS
quant/SalmonQuantify.cpp
quant/pipeline/MappingPipelineStages.cpp
)

set (SALMON_IO_SRCS
io/FASTAParser.cpp
io/AlignmentIO.cpp
${SALMON_FQFEEDER_SOURCE_DIR}/src/FastxParser.cpp
)
set_source_files_properties(${SALMON_FQFEEDER_SOURCE_DIR}/src/FastxParser.cpp
  PROPERTIES COMPILE_DEFINITIONS "fastx_parser=salmon_fqfeeder")

set (SALMON_UTIL_SRCS
util/InputFileUtils.cpp
util/LibraryTypeUtils.cpp
util/QuantOptionsUtils.cpp
util/SalmonUtils.cpp
util/DistributionUtils.cpp
util/SalmonExceptions.cpp
util/SalmonStringUtils.cpp
util/TemplateInstantiations.cpp
)

set (SALMON_SUPPORT_SRCS
${GAT_SOURCE_DIR}/src/jellyfish/mer_dna.cc
backtrace.cc
xxhash.c
edlib.cpp
model/GenomicFeature.cpp
)

# check if we know how to do IPO
check_ipo_supported(RESULT HAS_IPO)

if(HAS_IPO AND SALMON_ENABLE_LTO)
  set_property(TARGET ksw2pp PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

if(SALMON_ENABLE_TESTS)
  set (UNIT_TESTS_INDIVIDUAL_SRCS
    ${GAT_SOURCE_DIR}/tests/UnitTests.cpp
    ${GAT_SOURCE_DIR}/src/model/FragmentLengthDistribution.cpp
    ${GAT_SOURCE_DIR}/tests/GCSampleTests.cpp
    ${GAT_SOURCE_DIR}/tests/LibraryTypeTests.cpp
  )
endif()

# Set the RPATH
if (NOT APPLE)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:$ORIGIN/")
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

add_library(salmon_project_options INTERFACE)
target_compile_features(salmon_project_options INTERFACE cxx_std_17)
target_compile_definitions(salmon_project_options INTERFACE
  RAPMAP_SALMON_SUPPORT=1
  PUFFERFISH_SALMON_SUPPORT=1
  HAVE_ANSI_TERM=1
  HAVE_SSTREAM=1
  STX_NO_STD_STRING_VIEW=1
  span_FEATURE_MAKE_SPAN_TO_STD=14
)
target_compile_options(salmon_project_options INTERFACE
  "$<$<CONFIG:DEBUG>:${TGT_DEBUG_FLAGS}>"
  "$<$<CONFIG:RELEASE>:${TGT_RELEASE_FLAGS}>"
)

add_library(salmon_repo_headers INTERFACE)
target_include_directories(salmon_repo_headers INTERFACE
  ${GAT_SOURCE_DIR}/include
)

add_library(salmon_pufferfish_preferred_headers INTERFACE)
target_include_directories(salmon_pufferfish_preferred_headers INTERFACE
  ${SALMON_PUFFERFISH_SOURCE_DIR}/include
)

add_library(salmon_fqfeeder_headers INTERFACE)
target_include_directories(salmon_fqfeeder_headers INTERFACE
  ${SALMON_FQFEEDER_SOURCE_DIR}/include
)

add_library(salmon_eigen_headers INTERFACE)
target_include_directories(salmon_eigen_headers INTERFACE
  ${GAT_SOURCE_DIR}/include/eigen3
)

add_library(salmon_vendor_headers INTERFACE)
target_include_directories(salmon_vendor_headers INTERFACE
  ${GAT_SOURCE_DIR}/external
  ${GAT_SOURCE_DIR}/external/cereal/include
  ${SALMON_PUFFERFISH_SOURCE_DIR}/include
  ${COMPACT_VECTOR_INCLUDE_PATH}
)

add_library(salmon_system_headers INTERFACE)
target_include_directories(salmon_system_headers SYSTEM INTERFACE
  ${SALMON_ZLIB_INCLUDE_DIRS}
  ${HTSLIB_INCLUDE_DIR}
  ${TBB_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${CEREAL_INCLUDE_DIRS}
  ${LIB_GFF_INCLUDE_DIR}
  ${ICU_INC_DIRS}
)

add_library(salmon_common_link_deps INTERFACE)
set(SALMON_COMMON_LINK_DIRS
  ${Boost_LIBRARY_DIRS}
  ${TBB_LIBRARY_DIRS}
  ${LAPACK_LIBRARY_DIR}
  ${BLAS_LIBRARY_DIR}
  ${LIB_GFF_LIBRARY_DIR}
)
set(SALMON_EXISTING_LINK_DIRS "")
foreach(salmon_link_dir IN LISTS SALMON_COMMON_LINK_DIRS)
  if(salmon_link_dir AND EXISTS "${salmon_link_dir}")
    list(APPEND SALMON_EXISTING_LINK_DIRS "${salmon_link_dir}")
  endif()
endforeach()
if(SALMON_EXISTING_LINK_DIRS)
  target_link_directories(salmon_common_link_deps INTERFACE
    ${SALMON_EXISTING_LINK_DIRS}
  )
endif()

function(salmon_configure_target target_name)
  target_link_libraries(${target_name}
    PRIVATE
      salmon_project_options
      salmon_fqfeeder_headers
      salmon_pufferfish_preferred_headers
      salmon_repo_headers
  )
  if (USE_ARM)
    target_compile_definitions(${target_name} PUBLIC KSW_USE_ARM=1)
  endif()
  if(HAS_IPO AND SALMON_ENABLE_LTO)
    set_property(TARGET ${target_name} PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
  endif()
endfunction()

function(salmon_attach_external_headers target_name)
  target_link_libraries(${target_name}
    PRIVATE
      salmon_vendor_headers
      salmon_system_headers
  )
endfunction()

function(salmon_attach_eigen_headers target_name)
  target_link_libraries(${target_name}
    PRIVATE
      salmon_eigen_headers
  )
endfunction()

function(salmon_add_module_library target_name)
  set(options USE_EIGEN)
  set(one_value_args)
  set(multi_value_args SOURCES LINK_LIBS)
  cmake_parse_arguments(SALMON_MODULE "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  add_library(${target_name} STATIC ${SALMON_MODULE_SOURCES})
  salmon_configure_target(${target_name})
  salmon_attach_external_headers(${target_name})
  if(SALMON_MODULE_USE_EIGEN)
    salmon_attach_eigen_headers(${target_name})
  endif()
  if(SALMON_MODULE_LINK_LIBS)
    target_link_libraries(${target_name}
      PRIVATE
        ${SALMON_MODULE_LINK_LIBS}
    )
  endif()
endfunction()

salmon_add_module_library(salmon_io
  USE_EIGEN
  SOURCES ${SALMON_IO_SRCS}
  LINK_LIBS salmon_util
)
target_include_directories(salmon_io BEFORE PRIVATE
  ${SALMON_FQFEEDER_SOURCE_DIR}/include
)

salmon_add_module_library(salmon_model
  USE_EIGEN
  SOURCES ${SALMON_MODEL_SRCS}
  LINK_LIBS salmon_util
)

salmon_add_module_library(salmon_util
  USE_EIGEN
  SOURCES ${SALMON_UTIL_SRCS}
)

salmon_add_module_library(salmon_output
  USE_EIGEN
  SOURCES ${SALMON_OUTPUT_SRCS}
  LINK_LIBS salmon_util salmon_model
)

salmon_add_module_library(salmon_alignment
  USE_EIGEN
  SOURCES ${SALMON_ALIGNMENT_SRCS}
  LINK_LIBS salmon_io salmon_model salmon_util salmon_output salmon_inference salmon_options
)

salmon_add_module_library(salmon_inference
  USE_EIGEN
  SOURCES ${SALMON_INFERENCE_SRCS}
  LINK_LIBS salmon_model salmon_util
)

salmon_add_module_library(salmon_cli
  USE_EIGEN
  SOURCES ${SALMON_CLI_SRCS}
)

salmon_add_module_library(salmon_indexing
  USE_EIGEN
  SOURCES ${SALMON_INDEXING_SRCS}
  LINK_LIBS salmon_util
)

salmon_add_module_library(salmon_options
  SOURCES ${SALMON_OPTIONS_SRCS}
)

salmon_add_module_library(salmon_quant
  USE_EIGEN
  SOURCES ${SALMON_QUANT_SRCS}
  LINK_LIBS salmon_io salmon_model salmon_util salmon_output salmon_inference salmon_options
)

salmon_add_module_library(salmon_core
  SOURCES ${SALMON_SUPPORT_SRCS}
  LINK_LIBS salmon_io salmon_model salmon_util
)

add_library(salmon_support ALIAS salmon_core)

# Build the salmon executable
add_executable(salmon ${SALMON_ENTRY_SRCS})
salmon_configure_target(salmon)
salmon_attach_external_headers(salmon)
salmon_attach_eigen_headers(salmon)
if(FAST_MALLOC_OBJECT)
  # Prefer mimalloc's static override object on the final executable link line.
  target_link_options(salmon PRIVATE "${FAST_MALLOC_OBJECT}")
  if(TARGET mimalloc-obj-target)
    add_dependencies(salmon mimalloc-obj-target)
  endif()
endif()
set(_salmon_module_libs
  salmon_cli
  salmon_indexing
  salmon_quant
  salmon_alignment
  salmon_core
  # Be explicit for static-link toolchains where transitive module deps from
  # intermediate static libraries are not reliably pulled in.
  salmon_io
  salmon_model
  salmon_util
  salmon_output
  salmon_inference
  salmon_options
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND UNIX AND NOT APPLE)
  target_link_libraries(salmon
    PRIVATE
      -Wl,--start-group
      ${_salmon_module_libs}
      -Wl,--end-group
  )
else()
  target_link_libraries(salmon
    PRIVATE
      ${_salmon_module_libs}
  )
endif()

if(SALMON_ENABLE_TESTS)
  add_executable(unitTests ${UNIT_TESTS_INDIVIDUAL_SRCS})
  salmon_configure_target(unitTests)
  salmon_attach_external_headers(unitTests)
  salmon_attach_eigen_headers(unitTests)
  add_custom_command(TARGET unitTests POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:unitTests> ${GAT_SOURCE_DIR}/tests/$<TARGET_FILE_NAME:unitTests>
      COMMENT "Copying unitTests"
  )
endif()

# Some pufferfish static dependencies are not always propagated to downstream
# consumers under all toolchain/static-link combinations.
set(_salmon_puffer_aux_libs)
foreach(_puff_aux IN ITEMS ksw2pp twopaco graphdump ntcard)
  if(TARGET ${_puff_aux})
    list(APPEND _salmon_puffer_aux_libs ${_puff_aux})
  endif()
endforeach()

target_link_libraries(salmon_common_link_deps INTERFACE
    Threads::Threads
    gff
    ${Boost_LIBRARIES}
    ${ICU_LIBS}
    ${CURL_LIBRARIES}
    ${SALMON_ZLIB_LIBRARIES}
    m
    ${HTSLIB_LIBRARIES}
    ${LIBLZMA_LIBRARIES}
    ${BZIP2_LIBRARIES}
    ${LIBSALMON_LINKER_FLAGS}
    ${NON_APPLECLANG_LIBS}
    pufferfish::puffer
    ${_salmon_puffer_aux_libs}
    ${ASAN_LIB}
    ${FAST_MALLOC_LIB}
    TBB::tbb
    ${SALMON_TBB_MALLOC_LIB}
    ${LIBRT}
    ${CMAKE_DL_LIBS}
)

# GNU ld with static archives can fail to resolve some intra-/inter-archive
# symbols in a single pass for fetched fallback dependencies.
if(UNIX AND NOT APPLE)
  set(_salmon_rescan_libs
    gff
    ${Boost_LIBRARIES}
    ${HTSLIB_LIBRARIES}
    pufferfish::puffer
    ${_salmon_puffer_aux_libs}
  )

  # Prefer CMake's native link-group abstraction when supported.
  target_link_libraries(salmon_common_link_deps INTERFACE
    "$<LINK_GROUP:RESCAN,${_salmon_rescan_libs}>"
  )

  # Fallback for toolchains/CMake variants where LINK_GROUP does not materialize
  # in the final link command: inject explicit GNU ld group markers.
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(salmon_common_link_deps INTERFACE
      -Wl,--start-group
      ${_salmon_rescan_libs}
      -Wl,--end-group
    )
  endif()
endif()

# Link the executable
target_link_libraries(salmon
    PRIVATE
      salmon_common_link_deps
)

add_dependencies(salmon salmon_core)

if(SALMON_ENABLE_TESTS)
target_link_libraries(unitTests
    PRIVATE
      ${SALMON_CATCH2_TARGET}
      salmon_common_link_deps
      salmon_core
)
endif()

if(NOT Iconv_IS_BUILT_IN)
  if(SALMON_ENABLE_TESTS)
    target_link_libraries(unitTests PRIVATE Iconv::Iconv)
  endif()
  target_link_libraries(salmon PRIVATE Iconv::Iconv)
endif()

if(SALMON_ENABLE_TESTS)
  add_dependencies(salmon unitTests)
  add_dependencies(unitTests salmon_core)
endif()

if(APPLE)
  find_program(SALMON_OTOOL_EXECUTABLE otool)
  if(SALMON_OTOOL_EXECUTABLE)
    add_custom_command(TARGET salmon POST_BUILD
      COMMAND ${CMAKE_COMMAND}
              -DSALMON_BIN=$<TARGET_FILE:salmon>
              -DSALMON_OTOOL=${SALMON_OTOOL_EXECUTABLE}
              -P ${GAT_SOURCE_DIR}/cmake/CheckNoSystemZlib.cmake
      COMMENT "Verifying salmon does not link against /usr/lib/libz.1.dylib"
      VERBATIM)
  endif()
endif()

if(UNIX)
  set(_salmon_enforce_static_tbb OFF)
  if(FETCHED_TBB AND NOT APPLE)
    set(_salmon_enforce_static_tbb ON)
  endif()
  add_custom_command(TARGET salmon POST_BUILD
    COMMAND ${CMAKE_COMMAND}
            -DSALMON_BIN=$<TARGET_FILE:salmon>
            -DSALMON_OTOOL=${SALMON_OTOOL_EXECUTABLE}
            -DSALMON_ENFORCE_STATIC_TBB=${_salmon_enforce_static_tbb}
            -P ${GAT_SOURCE_DIR}/cmake/CheckTBBLinkage.cmake
    COMMENT "Auditing TBB linkage mode for salmon"
    VERBATIM)
endif()

##
# External dependencies of salmon_core and salmon
##
set(SALMON_EXTERNAL_DEP_TARGETS
  salmon_io
  salmon_model
  salmon_util
  salmon_output
  salmon_alignment
  salmon_inference
  salmon_cli
  salmon_indexing
  salmon_options
  salmon_quant
  salmon_core
  salmon
)
if(SALMON_ENABLE_TESTS)
  list(APPEND SALMON_EXTERNAL_DEP_TARGETS unitTests)
endif()

function(salmon_add_external_dependency dependency_target)
  foreach(salmon_external_target IN LISTS SALMON_EXTERNAL_DEP_TARGETS)
    add_dependencies(${salmon_external_target} ${dependency_target})
  endforeach()
endfunction()

if (${FETCHED_HTSLIB})
  salmon_add_external_dependency(libhtslib)
endif()

if (${FETCHED_BOOST})
  salmon_add_external_dependency(libboost)
endif()

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
    "${GAT_SOURCE_DIR}" CACHE PATH "Default install prefix" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(INSTALL_LIB_DIR lib )
set(INSTALL_BIN_DIR bin )
set(INSTALL_INCLUDE_DIR include )

if(TARGET TBB::tbb)
#set(TBB_SOURCE_DIR $<TARGET_FILE:TBB::tbb>)
#add_custom_target(genexdebug COMMAND ${CMAKE_COMMAND} -E echo "$<TARGET_LINKER_FILE:TBB::tbb>")
get_target_property(TBB_LIB_INSTALL_NAME TBB::tbb IMPORTED_LOCATION_RELEASE)
if(NOT TBB_LIB_INSTALL_NAME)
  get_target_property(TBB_LIB_INSTALL_NAME TBB::tbb IMPORTED_LOCATION)
endif()
if(TBB_LIB_INSTALL_NAME)
get_filename_component(TBB_LIB_INSTALL_DIR ${TBB_LIB_INSTALL_NAME} DIRECTORY)
message("TBB_LIB_INSTALL_DIR = ${TBB_LIB_INSTALL_DIR}")
file(GLOB TBB_FILES ${TBB_LIB_INSTALL_DIR}/libtbb*.${SHARED_LIB_EXTENSION}*)
message("TBBGLOBS = ${TBB_FILES}")

install(FILES 
  ${TBB_LIB_INSTALL_NAME}
  DESTINATION ${INSTALL_LIB_DIR}
) 
endif()
endif()

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
