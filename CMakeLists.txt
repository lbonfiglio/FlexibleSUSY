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
  find_package(Boost 1.37.0 REQUIRED)
  find_package(Eigen3 3.1 REQUIRED NO_MODULE)
endif()

add_subdirectory(config)
add_subdirectory(doc)

if(ENABLE_COMPILE)
  add_subdirectory(fflite)
  add_subdirectory(slhaea)
  add_subdirectory(src)
endif()
