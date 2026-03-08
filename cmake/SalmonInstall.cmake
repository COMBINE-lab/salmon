include_guard(GLOBAL)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${GAT_SOURCE_DIR}" CACHE PATH "Default install prefix" FORCE)
endif()

set(INSTALL_LIB_DIR lib)
set(INSTALL_BIN_DIR bin)
set(INSTALL_INCLUDE_DIR include)

install(TARGETS salmon salmon_core
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

set(POST_INSTALL_SCRIPT ${GAT_SOURCE_DIR}/cmake/PostInstall.cmake)
install(
  CODE
  "
  execute_process(COMMAND \"${CMAKE_COMMAND}\"
                          -DCMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}
                          -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                          -P \"${POST_INSTALL_SCRIPT}\")
  ")

include(InstallRequiredSystemLibraries)
add_test(
  NAME unit_tests
  COMMAND ${CMAKE_COMMAND}
          -DTOPLEVEL_DIR=${GAT_SOURCE_DIR}
          -DUNIT_TEST_EXECUTABLE=$<TARGET_FILE:unitTests>
          -P ${GAT_SOURCE_DIR}/cmake/UnitTests.cmake)
add_test(
  NAME salmon_read_test_quasi
  COMMAND ${CMAKE_COMMAND}
          -DTOPLEVEL_DIR=${GAT_SOURCE_DIR}
          -DSALMON_EXECUTABLE=$<TARGET_FILE:salmon>
          -P ${GAT_SOURCE_DIR}/cmake/TestSalmonQuasi.cmake)
add_test(
  NAME salmon_alevin_stub
  COMMAND ${CMAKE_COMMAND}
          -DSALMON_EXECUTABLE=$<TARGET_FILE:salmon>
          -P ${GAT_SOURCE_DIR}/cmake/TestAlevinStub.cmake)

if(SALMON_ENABLE_BENCHMARKS)
  add_test(
    NAME salmon_benchmark_smoke
    COMMAND ${CMAKE_COMMAND}
            -DTOPLEVEL_DIR=${GAT_SOURCE_DIR}
            -DSALMON_EXECUTABLE=$<TARGET_FILE:salmon>
            -DPYTHON_EXECUTABLE=${Python3_EXECUTABLE}
            -DBENCHMARK_OUTPUT=${CMAKE_BINARY_DIR}/benchmark-smoke.json
            -P ${GAT_SOURCE_DIR}/cmake/TestBenchmarkSmoke.cmake)
endif()
