# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

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
CMAKE_COMMAND = /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Downloads/cmake-3.9.0-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Downloads/cmake-3.9.0-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox/bin

# Include any dependencies generated for this target.
include CMakeFiles/ballpivoting.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ballpivoting.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ballpivoting.dir/flags.make

CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o: CMakeFiles/ballpivoting.dir/flags.make
CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o: ../src/srtools/BPAexe.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o -c /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox/src/srtools/BPAexe.cpp

CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox/src/srtools/BPAexe.cpp > CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.i

CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox/src/srtools/BPAexe.cpp -o CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.s

CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o.requires:

.PHONY : CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o.requires

CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o.provides: CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o.requires
	$(MAKE) -f CMakeFiles/ballpivoting.dir/build.make CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o.provides.build
.PHONY : CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o.provides

CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o.provides.build: CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o


# Object files for target ballpivoting
ballpivoting_OBJECTS = \
"CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o"

# External object files for target ballpivoting
ballpivoting_EXTERNAL_OBJECTS =

ballpivoting: CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o
ballpivoting: CMakeFiles/ballpivoting.dir/build.make
ballpivoting: CMakeFiles/ballpivoting.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ballpivoting"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ballpivoting.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ballpivoting.dir/build: ballpivoting

.PHONY : CMakeFiles/ballpivoting.dir/build

CMakeFiles/ballpivoting.dir/requires: CMakeFiles/ballpivoting.dir/src/srtools/BPAexe.cpp.o.requires

.PHONY : CMakeFiles/ballpivoting.dir/requires

CMakeFiles/ballpivoting.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ballpivoting.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ballpivoting.dir/clean

CMakeFiles/ballpivoting.dir/depend:
	cd /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox/bin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox/bin /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox/bin /nfs/punlinuxvault2/fcdata13/visitor/lgiaccar/Surface-Reconstruction-Toolbox/bin/CMakeFiles/ballpivoting.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ballpivoting.dir/depend
