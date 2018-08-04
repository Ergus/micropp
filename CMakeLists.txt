
cmake_minimum_required(VERSION 3.0)


# find_program(MCXX mcxx)
# find_program(MCC mcc)
# find_program(MFC mfc)

# SET(CMAKE_CXX_COMPILER ${MCXX} CACHE INTERNAL "" FORCE)
# SET(CMAKE_Fortran_COMPILER ${MFC} CACHE INTERNAL "" FORCE)
# SET(CMAKE_C_COMPILER ${MCC} CACHE INTERNAL "" FORCE)

project(MicroPP VERSION 0.1
  LANGUAGES C CXX Fortran)


# add_definitions(--ompss-2 -DNANOS6)
# set(CMAKE_EXE_LINKER_FLAGS --ompss-2)


# Use C++11
if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

option(TIMER "Enable time benchmarks instrumentation" OFF)

if (TIMER)
  add_definitions(-DTIMER)
endif()

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR})

# wildcard all the sources in src
#file(GLOB SOURCESLIB src/*.cpp src/*.f95)
file(GLOB SOURCESLIB src/*.cpp)

# Library
add_library(micropp ${SOURCESLIB})

# Enable auto create tests
enable_testing()
add_subdirectory(test)
