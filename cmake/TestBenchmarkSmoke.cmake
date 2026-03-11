if(NOT DEFINED SALMON_EXECUTABLE)
  message(FATAL_ERROR "SALMON_EXECUTABLE is required")
endif()

if(NOT DEFINED TOPLEVEL_DIR)
  message(FATAL_ERROR "TOPLEVEL_DIR is required")
endif()

if(NOT DEFINED PYTHON_EXECUTABLE)
  message(FATAL_ERROR "PYTHON_EXECUTABLE is required")
endif()

if(NOT DEFINED BENCHMARK_OUTPUT)
  set(BENCHMARK_OUTPUT "${TOPLEVEL_DIR}/tests/benchmarks/smoke-benchmark.json")
endif()

execute_process(
  COMMAND ${PYTHON_EXECUTABLE}
          ${TOPLEVEL_DIR}/tests/benchmarks/smoke_benchmark.py
          --salmon ${SALMON_EXECUTABLE}
          --repo-root ${TOPLEVEL_DIR}
          --output ${BENCHMARK_OUTPUT}
  RESULT_VARIABLE BENCH_RESULT
)

if(BENCH_RESULT)
  message(FATAL_ERROR "Smoke benchmark failed with exit code ${BENCH_RESULT}")
endif()
