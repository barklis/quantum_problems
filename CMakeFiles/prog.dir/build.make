# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems

# Include any dependencies generated for this target.
include CMakeFiles/prog.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/prog.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/prog.dir/flags.make

CMakeFiles/prog.dir/polar_system.c.o: CMakeFiles/prog.dir/flags.make
CMakeFiles/prog.dir/polar_system.c.o: polar_system.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/prog.dir/polar_system.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/prog.dir/polar_system.c.o   -c /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems/polar_system.c

CMakeFiles/prog.dir/polar_system.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/prog.dir/polar_system.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems/polar_system.c > CMakeFiles/prog.dir/polar_system.c.i

CMakeFiles/prog.dir/polar_system.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/prog.dir/polar_system.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems/polar_system.c -o CMakeFiles/prog.dir/polar_system.c.s

CMakeFiles/prog.dir/polar_system.c.o.requires:

.PHONY : CMakeFiles/prog.dir/polar_system.c.o.requires

CMakeFiles/prog.dir/polar_system.c.o.provides: CMakeFiles/prog.dir/polar_system.c.o.requires
	$(MAKE) -f CMakeFiles/prog.dir/build.make CMakeFiles/prog.dir/polar_system.c.o.provides.build
.PHONY : CMakeFiles/prog.dir/polar_system.c.o.provides

CMakeFiles/prog.dir/polar_system.c.o.provides.build: CMakeFiles/prog.dir/polar_system.c.o


# Object files for target prog
prog_OBJECTS = \
"CMakeFiles/prog.dir/polar_system.c.o"

# External object files for target prog
prog_EXTERNAL_OBJECTS =

prog: CMakeFiles/prog.dir/polar_system.c.o
prog: CMakeFiles/prog.dir/build.make
prog: /usr/lib/x86_64-linux-gnu/libgsl.so
prog: /usr/lib/x86_64-linux-gnu/libgslcblas.so
prog: CMakeFiles/prog.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable prog"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/prog.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/prog.dir/build: prog

.PHONY : CMakeFiles/prog.dir/build

CMakeFiles/prog.dir/requires: CMakeFiles/prog.dir/polar_system.c.o.requires

.PHONY : CMakeFiles/prog.dir/requires

CMakeFiles/prog.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/prog.dir/cmake_clean.cmake
.PHONY : CMakeFiles/prog.dir/clean

CMakeFiles/prog.dir/depend:
	cd /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems /home/barklis/Dokumenty/IT/Programowanie/C/Nowe/Kwanty/quantum_problems/CMakeFiles/prog.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/prog.dir/depend

