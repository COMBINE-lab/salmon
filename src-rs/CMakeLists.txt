set(CARGO_CMD cargo build --release)
set(TARGET_DIR "release")

set(ALEVINX_SO "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_DIR}/libalvtools.so")

add_custom_target(alevinx ALL
    COMMENT "Compiling alevinx rust module"
    COMMAND CARGO_TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR} ${CARGO_CMD} 
    COMMAND cp ${ALEVINX_SO} ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(alevinx PROPERTIES LOCATION ${CMAKE_CURRENT_BINARY_DIR})
