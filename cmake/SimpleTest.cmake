execute_process(COMMAND tar xzvf sample_data.tgz
                WORKING_DIRECTORY ${TOPLEVEL_DIR}
                RESULT_VARIABLE TAR_RESULT
               )

if (TAR_RESULT)
    message(FATAL_ERROR "Error untarring sample_data.tgz")
endif()

set(INDEX_CMD ${TOPLEVEL_DIR}/build/src/sailfish --no-version-check index -t transcripts.fasta -k 20 -o sample_index --force)
execute_process(COMMAND ${INDEX_CMD}
                WORKING_DIRECTORY ${TOPLEVEL_DIR}/sample_data
                RESULT_VARIABLE INDEX_RESULT
                )

if (INDEX_RESULT)
    message(FATAL_ERROR "Error running ${INDEX_COMMAND}")
endif()

set(QUANT_COMMAND ${TOPLEVEL_DIR}/build/src/sailfish --no-version-check quant -i sample_index --noBiasCorrect -l IU -1 reads_1.fastq -2 reads_2.fastq -o sample_quant)
execute_process(COMMAND ${QUANT_COMMAND}
	            WORKING_DIRECTORY ${TOPLEVEL_DIR}/sample_data
	            RESULT_VARIABLE QUANT_RESULT
                )
if (QUANT_RESULT)
    message(FATAL_ERROR "Error running ${QUANT_RESULT}")
endif()

if (EXISTS ${TOPLEVEL_DIR}/sample_data/sample_quant/quant.sf)
	message("Sailfish ran successfully")
else()
	message(FATAL_ERROR "Sailfish failed to produce output")
endif()


