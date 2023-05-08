# CMake generated Testfile for 
# Source directory: /home/zzare/salmon/src
# Build directory: /home/zzare/salmon/build/src
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(unit_tests "/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake" "-DTOPLEVEL_DIR=/home/zzare/salmon" "-P" "/home/zzare/salmon/cmake/UnitTests.cmake")
set_tests_properties(unit_tests PROPERTIES  _BACKTRACE_TRIPLES "/home/zzare/salmon/src/CMakeLists.txt;459;add_test;/home/zzare/salmon/src/CMakeLists.txt;0;")
add_test(salmon_read_test_quasi "/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake" "-DTOPLEVEL_DIR=/home/zzare/salmon" "-P" "/home/zzare/salmon/cmake/TestSalmonQuasi.cmake")
set_tests_properties(salmon_read_test_quasi PROPERTIES  _BACKTRACE_TRIPLES "/home/zzare/salmon/src/CMakeLists.txt;460;add_test;/home/zzare/salmon/src/CMakeLists.txt;0;")
