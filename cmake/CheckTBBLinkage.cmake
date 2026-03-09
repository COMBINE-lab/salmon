if(NOT DEFINED SALMON_BIN)
  message(FATAL_ERROR "SALMON_BIN is required")
endif()

if(NOT EXISTS "${SALMON_BIN}")
  message(FATAL_ERROR "SALMON_BIN does not exist: ${SALMON_BIN}")
endif()

set(_salmon_tbb_needs "")

if(APPLE)
  if(NOT DEFINED SALMON_OTOOL)
    find_program(SALMON_OTOOL otool)
  endif()
  if(NOT SALMON_OTOOL)
    message(WARNING "Skipping TBB linkage audit: otool not found")
    return()
  endif()
  execute_process(
    COMMAND "${SALMON_OTOOL}" -L "${SALMON_BIN}"
    OUTPUT_VARIABLE _salmon_link_info
    ERROR_VARIABLE _salmon_link_err
    RESULT_VARIABLE _salmon_link_rc
  )
  if(NOT _salmon_link_rc EQUAL 0)
    message(FATAL_ERROR "Failed to inspect linkage with otool: ${_salmon_link_err}")
  endif()
elseif(UNIX)
  find_program(SALMON_READELF readelf)
  if(SALMON_READELF)
    execute_process(
      COMMAND "${SALMON_READELF}" -d "${SALMON_BIN}"
      OUTPUT_VARIABLE _salmon_link_info
      ERROR_VARIABLE _salmon_link_err
      RESULT_VARIABLE _salmon_link_rc
    )
    if(NOT _salmon_link_rc EQUAL 0)
      message(FATAL_ERROR "Failed to inspect linkage with readelf: ${_salmon_link_err}")
    endif()
  else()
    find_program(SALMON_LDD ldd)
    if(NOT SALMON_LDD)
      message(WARNING "Skipping TBB linkage audit: neither readelf nor ldd found")
      return()
    endif()
    execute_process(
      COMMAND "${SALMON_LDD}" "${SALMON_BIN}"
      OUTPUT_VARIABLE _salmon_link_info
      ERROR_VARIABLE _salmon_link_err
      RESULT_VARIABLE _salmon_link_rc
    )
    if(NOT _salmon_link_rc EQUAL 0)
      message(FATAL_ERROR "Failed to inspect linkage with ldd: ${_salmon_link_err}")
    endif()
  endif()
else()
  message(STATUS "Skipping TBB linkage audit: unsupported platform")
  return()
endif()

string(REGEX MATCHALL "[^\n]+" _salmon_lines "${_salmon_link_info}")
foreach(_line IN LISTS _salmon_lines)
  if(_line MATCHES "libtbb|libtbbmalloc|tbb\\.so|tbbmalloc\\.so|tbb\\.dylib|tbbmalloc\\.dylib")
    string(APPEND _salmon_tbb_needs "${_line}\n")
  endif()
endforeach()

if(_salmon_tbb_needs STREQUAL "")
  message(STATUS "TBB linkage audit: no dynamic TBB dependency detected in salmon")
  return()
endif()

message(STATUS "TBB linkage audit: dynamic TBB dependency detected")
string(REGEX REPLACE "\n$" "" _salmon_tbb_needs "${_salmon_tbb_needs}")
message(STATUS "TBB linkage entries:\n${_salmon_tbb_needs}")

if(DEFINED SALMON_ENFORCE_STATIC_TBB AND SALMON_ENFORCE_STATIC_TBB)
  message(FATAL_ERROR "Dynamic TBB dependency detected while static TBB was required")
endif()
