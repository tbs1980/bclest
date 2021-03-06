project(bclest CXX)
cmake_minimum_required(VERSION 2.8.0)

# assume that we are in debug mode if build type not specified
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# vscode cmake helper
include(CMakeToolsHelpers OPTIONAL)

# for additional custom cmake modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# coverage report
if(CMAKE_COMPILER_IS_GNUCXX)
  include(CodeCoverage)
  setup_target_for_coverage(${PROJECT_NAME}_coverage test coverage)
endif()

# custom compiler flags
include(set_compiler_flags)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

# bclest executable
add_subdirectory(bclest)

# unit testing
include(CTest)
enable_testing()
add_subdirectory(tests)