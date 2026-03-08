cmake_minimum_required(VERSION 3.24)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()

project(Salmon LANGUAGES C CXX)
find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND CMAKE_MODULE_PATH "${GAT_SOURCE_DIR}/cmake/Modules")

include("${GAT_SOURCE_DIR}/cmake/SalmonVersion.cmake")
include("${GAT_SOURCE_DIR}/cmake/SalmonProjectOptions.cmake")
include("${GAT_SOURCE_DIR}/cmake/SalmonDependencies.cmake")

set(CPACK_SOURCE_IGNORE_FILES
  "/src/util/PCA.cpp"
  "/src/util/PCAUtils.cpp"
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

message(STATUS "CPACK_SOURCE_IGNORE_FILES = ${CPACK_SOURCE_IGNORE_FILES}")

define_property(TARGET PROPERTY COMPACT_VECTOR_DIR INHERITED
  BRIEF_DOCS "the path to the directory containing the compact_vector include tree"
  FULL_DOCS "the path to the directory containing the compact_vector include tree")

get_property(COMPACT_VECTOR_INCLUDE_PATH TARGET graphdump PROPERTY COMPACT_VECTOR_DIR)
message(STATUS "fetched path for compact_vector as [${COMPACT_VECTOR_INCLUDE_PATH}]")

add_subdirectory(src)

add_custom_target(include_policy_check
  COMMAND ${CMAKE_COMMAND} -E env
          PYTHONPATH=
          ${Python3_EXECUTABLE} ${GAT_SOURCE_DIR}/scripts/check_include_policy.py
  WORKING_DIRECTORY ${GAT_SOURCE_DIR}
  COMMENT "Checking include policy (forbidden top-level vendor + legacy internal headers)")

add_custom_target(refactor_policy_check
  COMMAND ${CMAKE_COMMAND} -E env
          PYTHONPATH=
          ${Python3_EXECUTABLE} ${GAT_SOURCE_DIR}/scripts/check_refactor_policy.py
  WORKING_DIRECTORY ${GAT_SOURCE_DIR}
  COMMENT "Checking refactor policy (legacy knobs + forbidden raw symbols)")

if(SALMON_ENABLE_BENCHMARKS)
  add_custom_target(benchmark_smoke
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR} -R salmon_benchmark_smoke --output-on-failure
    COMMENT "Running Salmon smoke benchmark test")
endif()

include(CPack)

set(ARCHIVE_NAME ${CMAKE_PROJECT_NAME}-${PROJECT_VERSION})
add_custom_target(dist
  COMMAND git archive --prefix=${ARCHIVE_NAME}/ HEAD | gzip > ${CMAKE_BINARY_DIR}/${ARCHIVE_NAME}.tar.gz
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

include("${GAT_SOURCE_DIR}/cmake/SalmonInstall.cmake")
