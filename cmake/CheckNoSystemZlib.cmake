if(NOT DEFINED SALMON_BIN OR NOT DEFINED SALMON_OTOOL)
  message(FATAL_ERROR "CheckNoSystemZlib.cmake requires SALMON_BIN and SALMON_OTOOL")
endif()

execute_process(
  COMMAND "${SALMON_OTOOL}" -L "${SALMON_BIN}"
  OUTPUT_VARIABLE salmon_otool_out
  RESULT_VARIABLE salmon_otool_rc
)

if(NOT salmon_otool_rc EQUAL 0)
  message(FATAL_ERROR "Failed to inspect linked libraries for ${SALMON_BIN}")
endif()

string(FIND "${salmon_otool_out}" "/usr/lib/libz.1.dylib" salmon_has_system_zlib)
if(NOT salmon_has_system_zlib EQUAL -1)
  message(FATAL_ERROR
    "salmon links against /usr/lib/libz.1.dylib; zlib-ng compatibility backend is required.")
endif()
