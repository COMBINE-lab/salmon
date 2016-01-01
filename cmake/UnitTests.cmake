set(TEST_COMMAND ${TOPLEVEL_DIR}/tests/unitTests)
execute_process(COMMAND ${TEST_COMMAND}
	            WORKING_DIRECTORY ${TOPLEVEL_DIR}/tests
                RESULT_VARIABLE UNIT_TEST_RESULT
                )
if (UNIT_TEST_RESULT)
    message(FATAL_ERROR "Error running ${UNIT_TEST_RESULT}")
endif()


