cmake_minimum_required(VERSION 3.12)
project(gff)

message(STATUS "Install prefix : ${CMAKE_INSTALL_PREFIX}")

set(INCLUDE_INSTALL_DIR include/)
set(LIB_INSTALL_DIR lib/)

include(CMakePackageConfigHelpers)

set(ver_major 2)
set(ver_minor 0)
set(ver_patch 0)

set(CPACK_PACKAGE_VERSION_MAJOR ${ver_major})
set(CPACK_PACKAGE_VERSION_MINOR ${ver_minor})
set(CPACK_PACKAGE_VERSION_PATCH ${ver_patch})

set(CPACK_PACKAGE_VERSION "${ver_major}.${ver_minor}.${ver_patch}")
message("version: ${CPACK_PACKAGE_VERSION}")

set(PROJECT_VERSION ${CPACK_PACKAGE_VERSION})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -Wall -D_FILE_OFFSET_BITS=64" )

set(GFFLib_SRCS
    ${PROJECT_SOURCE_DIR}/codons.cpp
    ${PROJECT_SOURCE_DIR}/GArgs.cpp
    ${PROJECT_SOURCE_DIR}/GBase.cpp
    ${PROJECT_SOURCE_DIR}/gdna.cpp
    ${PROJECT_SOURCE_DIR}/GFaSeqGet.cpp
    ${PROJECT_SOURCE_DIR}/GFastaIndex.cpp
    ${PROJECT_SOURCE_DIR}/gff.cpp
    #${PROJECT_SOURCE_DIR}/gff_utils.cpp
    ${PROJECT_SOURCE_DIR}/GStr.cpp)

include_directories(${PROJECT_INCLUDE_DIR})

add_library(${PROJECT_NAME} STATIC ${GFFLib_SRCS})
add_executable(TestGFFParse ${PROJECT_SOURCE_DIR}/TestGFFParse.cpp)
target_link_libraries(TestGFFParse ${PROJECT_NAME})



configure_package_config_file(libgffConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/libgffConfig.cmake
  INSTALL_DESTINATION ${LIB_INSTALL_DIR}/libgff/cmake
  PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/libgffConfigVersion.cmake
  VERSION ${ver_major}.${ver_minor}.${ver_patch}
  COMPATIBILITY SameMajorVersion )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgffConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/libgffConfigVersion.cmake
        DESTINATION ${LIB_INSTALL_DIR}/libgff/cmake )

install(FILES ${CMAKE_BINARY_DIR}/lib${PROJECT_NAME}.a DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY ${PROJECT_INCLUDE_DIR} DESTINATION ${CMAKE_INSTALL_PREFIX})
