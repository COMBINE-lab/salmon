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
include external/pufferfish/src/CMakeFiles/ksw2pp.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include external/pufferfish/src/CMakeFiles/ksw2pp.dir/compiler_depend.make

# Include the progress variables for this target.
include external/pufferfish/src/CMakeFiles/ksw2pp.dir/progress.make

# Include the compile flags for this target's objects.
include external/pufferfish/src/CMakeFiles/ksw2pp.dir/flags.make

# Object files for target ksw2pp
ksw2pp_OBJECTS =

# External object files for target ksw2pp
ksw2pp_EXTERNAL_OBJECTS = \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_sse2.dir/ksw2pp/ksw2_extd2_sse.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_sse2.dir/ksw2pp/ksw2_extf2_sse.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_sse2.dir/ksw2pp/ksw2_extz2_sse.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_sse4.dir/ksw2pp/ksw2_extd2_sse.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_sse4.dir/ksw2pp/ksw2_extf2_sse.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_sse4.dir/ksw2pp/ksw2_extz2_sse.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.o" \
"/home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.o"

external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_sse2.dir/ksw2pp/ksw2_extd2_sse.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_sse2.dir/ksw2pp/ksw2_extf2_sse.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_sse2.dir/ksw2pp/ksw2_extz2_sse.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_sse4.dir/ksw2pp/ksw2_extd2_sse.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_sse4.dir/ksw2pp/ksw2_extf2_sse.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_sse4.dir/ksw2pp/ksw2_extz2_sse.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/kalloc.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extd.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_extz.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/ksw2_gg2_sse.c.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp_basic.dir/ksw2pp/KSW2Aligner.cpp.o
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp.dir/build.make
external/pufferfish/src/libksw2pp.a: external/pufferfish/src/CMakeFiles/ksw2pp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX static library libksw2pp.a"
	cd /home/zzare/salmon/build/external/pufferfish/src && $(CMAKE_COMMAND) -P CMakeFiles/ksw2pp.dir/cmake_clean_target.cmake
	cd /home/zzare/salmon/build/external/pufferfish/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ksw2pp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/pufferfish/src/CMakeFiles/ksw2pp.dir/build: external/pufferfish/src/libksw2pp.a
.PHONY : external/pufferfish/src/CMakeFiles/ksw2pp.dir/build

external/pufferfish/src/CMakeFiles/ksw2pp.dir/clean:
	cd /home/zzare/salmon/build/external/pufferfish/src && $(CMAKE_COMMAND) -P CMakeFiles/ksw2pp.dir/cmake_clean.cmake
.PHONY : external/pufferfish/src/CMakeFiles/ksw2pp.dir/clean

external/pufferfish/src/CMakeFiles/ksw2pp.dir/depend:
	cd /home/zzare/salmon/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zzare/salmon /home/zzare/salmon/external/pufferfish/src /home/zzare/salmon/build /home/zzare/salmon/build/external/pufferfish/src /home/zzare/salmon/build/external/pufferfish/src/CMakeFiles/ksw2pp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/pufferfish/src/CMakeFiles/ksw2pp.dir/depend

