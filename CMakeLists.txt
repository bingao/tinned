# Required version of CMake
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# Project name
project(Tinned
        VERSION 0.1.0
        LANGUAGES CXX)

# Host program could include this file by setting LIB_TINNED_PATH
if("${LIB_TINNED_PATH}" STREQUAL "")
    set(LIB_TINNED_PATH ${PROJECT_SOURCE_DIR})
endif()

# Sets the build mode: if the user specifies -DCMAKE_BUILD_TYPE on the command
# line, take their definition and dump it in the cache along with proper
# documentation, otherwise set CMAKE_BUILD_TYPE prior to calling project(),
# from http://www.cmake.org/pipermail/cmake/2008-September/023808.html
if(DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
else()
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif()

# Tinned needs C++11 standard
set(CMAKE_CXX_STANDARD 11)

# Options for making the library
option(BUILD_SHARED_LIBS "Build shared library." ON)

# Find SymEngine
set(CMAKE_PREFIX_PATH ${SymEngine_DIR} ${CMAKE_PREFIX_PATH})
find_package(SymEngine REQUIRED CONFIG
             PATH_SUFFIXES lib/cmake/symengine cmake/symengine CMake/)
message("SymEngine_DIR: " ${SymEngine_DIR})
message("SymEngine Version: " ${SymEngine_VERSION})

#set(CMAKE_BUILD_TYPE ${SYMENGINE_BUILD_TYPE})
set(CMAKE_CXX_FLAGS_RELEASE ${SYMENGINE_CXX_FLAGS_RELEASE})
set(CMAKE_CXX_FLAGS_DEBUG ${SYMENGINE_CXX_FLAGS_DEBUG})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYMENGINE_CXX_FLAGS}")

# Include directory of header files
include_directories(${LIB_TINNED_PATH}/include
                    ${SYMENGINE_INCLUDE_DIRS})

add_library(tinned
            ${LIB_TINNED_PATH}/src/Perturbation.cpp
            ${LIB_TINNED_PATH}/src/NonElecFunction.cpp
            ${LIB_TINNED_PATH}/src/ElectronicState.cpp
            ${LIB_TINNED_PATH}/src/OneElecDensity.cpp
            ${LIB_TINNED_PATH}/src/OneElecOperator.cpp
            ${LIB_TINNED_PATH}/src/TwoElecOperator.cpp
            ${LIB_TINNED_PATH}/src/ExchCorrEnergy.cpp
            ${LIB_TINNED_PATH}/src/ExchCorrPotential.cpp
            ${LIB_TINNED_PATH}/src/TemporumOperator.cpp
            ${LIB_TINNED_PATH}/src/TemporumOverlap.cpp
            ${LIB_TINNED_PATH}/src/KeepVisitor.cpp
            ${LIB_TINNED_PATH}/src/RemoveVisitor.cpp
            #${LIB_TINNED_PATH}/src/ReplaceVisitor.cpp
            ${LIB_TINNED_PATH}/src/StringifyVisitor.cpp)

#target_link_libraries()
