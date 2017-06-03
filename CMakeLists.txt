project(bclest CXX)
cmake_minimum_required(VERSION 2.8.0)

# vscode cmake helper
include(CMakeToolsHelpers OPTIONAL)

# for additional custom cmake modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# custom compiler flags
include(set_compiler_flags)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

# bclest executable
add_subdirectory(bclest)

# unit testing
enable_testing()
add_subdirectory(tests)