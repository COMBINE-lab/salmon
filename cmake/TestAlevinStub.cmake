if(NOT DEFINED SALMON_EXECUTABLE)
  message(FATAL_ERROR "SALMON_EXECUTABLE is required")
endif()

execute_process(
  COMMAND ${SALMON_EXECUTABLE} alevin
  RESULT_VARIABLE ALV_RESULT
  OUTPUT_VARIABLE ALV_STDOUT
  ERROR_VARIABLE ALV_STDERR
)

set(ALV_OUTPUT "${ALV_STDOUT}${ALV_STDERR}")

if(ALV_RESULT EQUAL 0)
  message(FATAL_ERROR "Expected `salmon alevin` to fail with non-zero exit status")
endif()

if(ALV_OUTPUT MATCHES "alevin-fry" AND ALV_OUTPUT MATCHES "v1\\.10\\.2")
  message("salmon alevin stub behavior is correct")
else()
  message(FATAL_ERROR
    "Stub output missing required migration breadcrumbs.\nOutput:\n${ALV_OUTPUT}")
endif()
