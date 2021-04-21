cmake_minimum_required(VERSION 3.10)

# set the project name
project(M2C VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# use mpicc and mpicxx to compile
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# compiler flags (turn on most of the warnings, except unused variables and parameters
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-reorder -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# folders for executable and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# find packages 
if(NOT DEFINED ENV{PETSC_DIR})
  message(FATAL_ERROR "ERROR: Could not find Petsc...")
endif()

# add include and link directories 
include_directories("${CMAKE_SOURCE_DIR}")
link_directories(${CMAKE_SOURCE_DIR}/parser)

# compile the parser
add_custom_target(
   extern_lib
   COMMAND $(MAKE)
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/parser
)

# define the macro that allows 1D Riemann solver to output solution
add_compile_definitions(PRINT_RIEMANN_SOLUTION=1)

# add the executable
add_executable(m2c 
Main.cpp
IoData.cpp
ExactRiemannSolverBase.cpp
Utils.cpp)

# link to libraries
target_link_libraries(m2c petsc mpi parser)
add_dependencies(m2c extern_lib)
