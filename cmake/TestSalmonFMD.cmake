execute_process(COMMAND tar xzvf sample_data.tgz
                WORKING_DIRECTORY ${TOPLEVEL_DIR}
                RESULT_VARIABLE TAR_RESULT
               )

if (TAR_RESULT)
    message(FATAL_ERROR "Error untarring sample_data.tgz")
endif()

set(SALMON_FMD_INDEX_CMD ${CMAKE_BINARY_DIR}/salmon index -t transcripts.fasta -i sample_salmon_fmd_index --type fmd)
execute_process(COMMAND ${SALMON_FMD_INDEX_CMD}
                WORKING_DIRECTORY ${TOPLEVEL_DIR}/sample_data
                RESULT_VARIABLE SALMON_FMD_INDEX_RESULT
                )

if (SALMON_FMD_INDEX_RESULT)
    message(FATAL_ERROR "Error running ${SALMON_FMD_INDEX_COMMAND}")
endif()

set(SALMON_QUANT_COMMAND ${CMAKE_BINARY_DIR}/salmon quant -i sample_salmon_fmd_index -l IU -1 reads_1.fastq -2 reads_2.fastq -o sample_salmon_fmd_quant)
execute_process(COMMAND ${SALMON_QUANT_COMMAND}
	            WORKING_DIRECTORY ${TOPLEVEL_DIR}/sample_data
                RESULT_VARIABLE SALMON_FMD_QUANT_RESULT
                )
if (SALMON_FMD_QUANT_RESULT)
    message(FATAL_ERROR "Error running ${SALMON_FMD_QUANT_RESULT}")
endif()

if (EXISTS ${TOPLEVEL_DIR}/sample_data/sample_salmon_fmd_quant/quant.sf)
    message("Salmon (read) ran successfully")
else()
    message(FATAL_ERROR "Salmon (read --- fmd-index) failed to produce output")
endif()
