# Required version of CMake
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# Project name
project(Tinned
        VERSION 0.1.0
        LANGUAGES CXX)

# Options for making the library
option(BUILD_SHARED_LIBS "Build shared library." OFF)
option(BUILD_TESTING "Build tests." ON)

# Find SymEngine, SymEngine_DIR can be set by users in advance
set(CMAKE_PREFIX_PATH ${SymEngine_DIR} ${CMAKE_PREFIX_PATH})
find_package(SymEngine REQUIRED CONFIG
             PATH_SUFFIXES lib/cmake/symengine cmake/symengine CMake/)
message("SymEngine_DIR: " ${SymEngine_DIR})
message("SymEngine version: " ${SymEngine_VERSION})

# Use the same build mode and C++ flags of SymEngine
set(CMAKE_BUILD_TYPE ${SYMENGINE_BUILD_TYPE})
set(CMAKE_CXX_FLAGS_RELEASE ${SYMENGINE_CXX_FLAGS_RELEASE})
set(CMAKE_CXX_FLAGS_DEBUG ${SYMENGINE_CXX_FLAGS_DEBUG})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYMENGINE_CXX_FLAGS}")

# Settings from host programs. Host programs can include this file by setting
# LIB_TINNED_PATH
if("${LIB_TINNED_PATH}" STREQUAL "")
    set(LIB_TINNED_PATH ${PROJECT_SOURCE_DIR})
endif()
if(NOT "${PARENT_DEFINITIONS}" STREQUAL "")
    foreach(_definition ${PARENT_DEFINITIONS})
        add_definitions(${_definition})
    endforeach()
endif()
if(NOT "${PARENT_INCLUDE_DIR}" STREQUAL "")
    include_directories(${PARENT_INCLUDE_DIR})
endif()
#if(NOT "${PARENT_MODULE_DIR}" STREQUAL "")
#    set(CMAKE_Fortran_MODULE_DIRECTORY ${PARENT_MODULE_DIR})
#endif()

# Include Tinned and SymEngine header file directories
set(TINNED_INCLUDE_DIRS ${LIB_TINNED_PATH}/include)
include_directories(${TINNED_INCLUDE_DIRS}
                    ${SYMENGINE_INCLUDE_DIRS})

# Tinned library source codes
add_subdirectory(src)

# Testing only in the main project, and not if it is included through
# add_subdirectory(). Learnt from Modern CMake at
# https://gitlab.com/CLIUtils/modern-cmake
if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) AND BUILD_TESTING)
    include(CTest)
    add_subdirectory(tests)
endif()

# Summary
message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message("CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message("CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message("LIB_TINNED_PATH: ${LIB_TINNED_PATH}")
message("TINNED_INCLUDE_DIRS: ${TINNED_INCLUDE_DIRS}")
message("BUILD_SHARED_LIBS: ${BUILD_SHARED_LIBS}")
message("BUILD_TESTING: ${BUILD_TESTING}")
