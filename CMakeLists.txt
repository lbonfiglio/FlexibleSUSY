cmake_minimum_required(VERSION 3.0)

project(FlexibleSUSY VERSION 2.0.1 LANGUAGES CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FlexibleSUSYBuildOptions)

if(ENABLE_META)
  find_package(Mathematica 7.0 REQUIRED)
endif()

if(ENABLE_COMPILE)
  find_package(Boost 1.37.0 COMPONENTS unit_test_framework REQUIRED)
  find_package(Eigen3 3.1 REQUIRED NO_MODULE)
  find_package(GSL REQUIRED)
  if(ENABLE_LAPACK)
    find_package(LAPACK REQUIRED)
  endif()
endif()

if(ENABLE_COMPILE)
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
