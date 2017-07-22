cmake_minimum_required(VERSION 2.8.9)

file(GLOB SOURCES "./src/*.cpp")
include_directories("./src/")
project (Surface-Reconstruction-Tools)
add_executable(robustcrust "./src/srtools/RCexe.cpp")
add_executable(scbmesher "./src/srtools/SCBexe.cpp")
add_executable(ballpivoting "./src/srtools/BPAexe.cpp")
add_executable(quickterrain "./src/srtools/QuickTerrain_exe.cpp")



