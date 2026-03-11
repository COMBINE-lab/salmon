include_guard(GLOBAL)

file(READ "${GAT_SOURCE_DIR}/current_version.txt" ver)

string(REGEX MATCH "VERSION_MAJOR ([0-9]*)" _ "${ver}")
set(ver_major "${CMAKE_MATCH_1}")

string(REGEX MATCH "VERSION_MINOR ([0-9]*)" _ "${ver}")
set(ver_minor "${CMAKE_MATCH_1}")

string(REGEX MATCH "VERSION_PATCH ([0-9]*)" _ "${ver}")
set(ver_patch "${CMAKE_MATCH_1}")

set(CPACK_PACKAGE_VERSION_MAJOR "${ver_major}")
set(CPACK_PACKAGE_VERSION_MINOR "${ver_minor}")
set(CPACK_PACKAGE_VERSION_PATCH "${ver_patch}")
set(CPACK_PACKAGE_VERSION
    "${ver_major}.${ver_minor}.${ver_patch}")
set(PROJECT_VERSION "${CPACK_PACKAGE_VERSION}")

message(STATUS "Salmon version: ${PROJECT_VERSION}")

set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "University of Maryland")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "Salmon - Wicked-fast RNA-seq isoform quantification using selective alignment")
set(CPACK_PACKAGE_NAME
    "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-Source")
