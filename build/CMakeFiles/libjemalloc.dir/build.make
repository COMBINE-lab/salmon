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

# Utility rule file for libjemalloc.

# Include any custom commands dependencies for this target.
include CMakeFiles/libjemalloc.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/libjemalloc.dir/progress.make

CMakeFiles/libjemalloc: CMakeFiles/libjemalloc-complete

CMakeFiles/libjemalloc-complete: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-install
CMakeFiles/libjemalloc-complete: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-mkdir
CMakeFiles/libjemalloc-complete: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-download
CMakeFiles/libjemalloc-complete: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-update
CMakeFiles/libjemalloc-complete: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-patch
CMakeFiles/libjemalloc-complete: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-configure
CMakeFiles/libjemalloc-complete: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-build
CMakeFiles/libjemalloc-complete: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-install
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Completed 'libjemalloc'"
	/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E make_directory /home/zzare/salmon/build/CMakeFiles
	/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E touch /home/zzare/salmon/build/CMakeFiles/libjemalloc-complete
	/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E touch /home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-done

libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-build: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-configure
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Performing build step for 'libjemalloc'"
	cd /home/zzare/salmon/external/jemalloc-5.2.1 && $(MAKE)
	cd /home/zzare/salmon/external/jemalloc-5.2.1 && /home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E touch /home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-build

libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-configure: libjemalloc-prefix/tmp/libjemalloc-cfgcmd.txt
libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-configure: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-patch
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Performing configure step for 'libjemalloc'"
	cd /home/zzare/salmon/external/jemalloc-5.2.1 && sh -c "CC=/usr/bin/gcc CFLAGS= CPPFLAGS= ./autogen.sh --disable-debug --enable-static --prefix=/home/zzare/salmon/external/install"
	cd /home/zzare/salmon/external/jemalloc-5.2.1 && /home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E touch /home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-configure

libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-download: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-custominfo.txt
libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-download: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-mkdir
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Performing download step for 'libjemalloc'"
	cd /home/zzare/salmon/external && curl -k -L https://github.com/jemalloc/jemalloc/archive/5.2.1.tar.gz -o jemalloc-5.2.1.tar.gz && /home/zzare/salmon/scripts/check_shasum.sh ed51b0b37098af4ca6ed31c22324635263f8ad6471889e0592a9c0dba9136aea jemalloc-5.2.1.tar.gz && tar -xzf jemalloc-5.2.1.tar.gz
	cd /home/zzare/salmon/external && /home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E touch /home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-download

libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-install: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-build
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Performing install step for 'libjemalloc'"
	cd /home/zzare/salmon/external/jemalloc-5.2.1 && cp -r lib /home/zzare/salmon/external/install/ && cp -r include /home/zzare/salmon/external/install/
	cd /home/zzare/salmon/external/jemalloc-5.2.1 && /home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E touch /home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-install

libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-mkdir:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Creating directories for 'libjemalloc'"
	/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -Dcfgdir= -P /home/zzare/salmon/build/libjemalloc-prefix/tmp/libjemalloc-mkdirs.cmake
	/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E touch /home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-mkdir

libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-patch: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-update
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "No patch step for 'libjemalloc'"
	/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E echo_append
	/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E touch /home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-patch

libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-update: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/zzare/salmon/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "No update step for 'libjemalloc'"
	/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E echo_append
	/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/bin/cmake -E touch /home/zzare/salmon/build/libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-update

libjemalloc: CMakeFiles/libjemalloc
libjemalloc: CMakeFiles/libjemalloc-complete
libjemalloc: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-build
libjemalloc: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-configure
libjemalloc: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-download
libjemalloc: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-install
libjemalloc: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-mkdir
libjemalloc: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-patch
libjemalloc: libjemalloc-prefix/src/libjemalloc-stamp/libjemalloc-update
libjemalloc: CMakeFiles/libjemalloc.dir/build.make
.PHONY : libjemalloc

# Rule to build all files generated by this target.
CMakeFiles/libjemalloc.dir/build: libjemalloc
.PHONY : CMakeFiles/libjemalloc.dir/build

CMakeFiles/libjemalloc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/libjemalloc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/libjemalloc.dir/clean

CMakeFiles/libjemalloc.dir/depend:
	cd /home/zzare/salmon/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zzare/salmon /home/zzare/salmon /home/zzare/salmon/build /home/zzare/salmon/build /home/zzare/salmon/build/CMakeFiles/libjemalloc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/libjemalloc.dir/depend

