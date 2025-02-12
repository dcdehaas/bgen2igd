cmake_minimum_required(VERSION 3.10)
project(bgen2igd)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds not allowed. Please make a build directory and run CMake from there.\n")
endif()

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/bgen/)

include_directories(${CMAKE_CURRENT_LIST_DIR} "bgen/include/")

add_executable(bgen2igd bgen2igd.cpp)
target_link_libraries(bgen2igd bgen)