# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake

# The command to remove a file.
RM = /home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zzare/salmon

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zzare/salmon/build

# Include any dependencies generated for this target.
include external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/compiler_depend.make

# Include the progress variables for this target.
include external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/progress.make

# Include the compile flags for this target's objects.
include external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/flags.make

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/flags.make
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.o: /home/zzare/salmon/external/pufferfish/src/ksw2pp/kalloc.c
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.o"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.o -MF CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.o.d -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.o -c /home/zzare/salmon/external/pufferfish/src/ksw2pp/kalloc.c

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.i"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/zzare/salmon/external/pufferfish/src/ksw2pp/kalloc.c > CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.i

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.s"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/zzare/salmon/external/pufferfish/src/ksw2pp/kalloc.c -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.s

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/flags.make
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.o: /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_extd.c
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.o"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.o -MF CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.o.d -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.o -c /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_extd.c

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.i"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_extd.c > CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.i

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.s"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_extd.c -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.s

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/flags.make
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.o: /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_extz.c
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.o"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.o -MF CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.o.d -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.o -c /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_extz.c

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.i"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_extz.c > CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.i

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.s"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_extz.c -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.s

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/flags.make
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.o: /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg.c
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.o"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.o -MF CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.o.d -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.o -c /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg.c

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.i"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg.c > CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.i

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.s"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg.c -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.s

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/flags.make
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.o: /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg2.c
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.o"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.o -MF CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.o.d -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.o -c /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg2.c

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.i"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg2.c > CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.i

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.s"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg2.c -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.s

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/flags.make
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.o: /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg2_sse.c
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.o"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.o -MF CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.o.d -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.o -c /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg2_sse.c

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.i"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg2_sse.c > CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.i

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.s"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/zzare/salmon/external/pufferfish/src/ksw2pp/ksw2_gg2_sse.c -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.s

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/flags.make
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.o: /home/zzare/salmon/external/pufferfish/src/ksw2pp/KSW2Aligner.cpp
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.o: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.o"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.o -MF CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.o.d -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.o -c /home/zzare/salmon/external/pufferfish/src/ksw2pp/KSW2Aligner.cpp

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.i"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zzare/salmon/external/pufferfish/src/ksw2pp/KSW2Aligner.cpp > CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.i

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.s"
	cd /home/zzare/salmon/build/external/pufferfish/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zzare/salmon/external/pufferfish/src/ksw2pp/KSW2Aligner.cpp -o CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.s

ksw2pp_basic: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.o
ksw2pp_basic: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.o
ksw2pp_basic: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.o
ksw2pp_basic: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.o
ksw2pp_basic: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.o
ksw2pp_basic: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.o
ksw2pp_basic: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.o
ksw2pp_basic: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/build.make
.PHONY : ksw2pp_basic

# Rule to build all files generated by this target.
external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/build: ksw2pp_basic
.PHONY : external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/build

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/clean:
	cd /home/zzare/salmon/build/external/pufferfish/src && $(CMAKE_COMMAND) -P CMakeFiles/ksw2pp_basic.dir/cmake_clean.cmake
.PHONY : external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/clean

external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/depend:
	cd /home/zzare/salmon/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zzare/salmon /home/zzare/salmon/external/pufferfish/src /home/zzare/salmon/build /home/zzare/salmon/build/external/pufferfish/src /home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/depend

