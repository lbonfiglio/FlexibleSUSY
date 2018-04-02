cmake_minimum_required(VERSION 3.1)

project(FlexibleSUSY VERSION 2.0.1 LANGUAGES CXX Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FlexibleSUSYColors)
include(FlexibleSUSYBuildOptions)
include(FlexibleSUSYOptionalPackages)

if(ENABLE_META)
  find_package(Mathematica 7.0 REQUIRED)
endif()

if(ENABLE_COMPILE)
  # required packages
  find_package(Boost 1.37.0 COMPONENTS unit_test_framework REQUIRED)
  find_package(Eigen3 3.1 REQUIRED NO_MODULE)
  find_package(GSL REQUIRED)
  # optional packages
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_optional_package(LAPACK ENABLE_LAPACK)
  find_optional_package(Threads ENABLE_THREADS)
endif()

if("${CMAKE_EXECUTABLE_SUFFIX}" STREQUAL "")
  set(CMAKE_EXECUTABLE_SUFFIX .x)
endif()

if(ENABLE_COMPILE)
  add_subdirectory(addons)
  add_subdirectory(config)
  add_subdirectory(fflite)
  add_subdirectory(model_specific)
  add_subdirectory(slhaea)
  add_subdirectory(src)
endif()

if(ENABLE_TEST AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)
  enable_testing()
  add_subdirectory(test)
endif()

add_subdirectory(doc)
