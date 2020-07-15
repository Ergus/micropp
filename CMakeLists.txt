
cmake_minimum_required(VERSION 3.0)

project(MicroPP VERSION 0.1 LANGUAGES C CXX Fortran)

if (NANOS6)
  find_program(MCXX mcxx)
  find_program(MCC mcc)
  find_program(MFC mfc)

  if (MCXX AND MCC AND MFC)

    message(STATUS "Building for nanos6")
    set(CMAKE_CXX_COMPILER ${MCXX} CACHE INTERNAL "" FORCE)
    set(CMAKE_Fortran_COMPILER ${MFC} CACHE INTERNAL "" FORCE)
    set(CMAKE_C_COMPILER ${MCC} CACHE INTERNAL "" FORCE)

    add_definitions("--ompss-2 -k -DNANOS6")
    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} --ompss-2 -k -DNANOS6")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --ompss-2 -k -DNANOS6")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} --ompss-2 -k -DNANOS6")

    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --ompss-2")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --ompss-2")
  else()
    message(FATAL_ERROR "Mercurium not found")
  endif()
endif()

project(MicroPP VERSION 0.1 LANGUAGES C CXX Fortran)

# Use C++11
if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
  endif ()
else ()
  set(CMAKE_CXX_STANDARD 11)
endif ()

option(TIMER "Enable time benchmarks instrumentation" OFF)
option(CGDEBUG "Enable CG algorithm debug" OFF)
option(FILTER "Enable Filtering on Average Stress and Ctan" OFF)
option(OPENMP "Enables OpenMP paralelization" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE)

set_property(CACHE CMAKE_BUILD_TYPE
  PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")

if(OPENMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

if (TIMER)
  add_definitions(-DTIMER)
endif()

if (CGDEBUG)
  add_definitions(-DCGDEBUG)
endif()

if (FILTER)
  add_definitions(-DFILTER)
endif()

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR})

# wildcard all the sources in src
file(GLOB SOURCESLIB src/*.c src/*.cpp src/*.f95)
#file(GLOB SOURCESLIB src/*.cpp)

add_library(micropp ${SOURCESLIB})

enable_testing()
add_subdirectory(test)

