# Distributed under the GPLv3 License.
# Author: Alexander Voigt
#
# FindLoopTools
# -------------
#
# Finds the LoopTools library [http://www.feynarts.de/looptools/].
#
# This module reads the following variables:
#
# LoopTools_BUILD_DIR     - LoopTools build/ directory
#
# This module defines the following variables:
#
# LoopTools_FOUND         - set if LoopTools has been found
# LoopTools_INCLUDE_DIR   - LoopTools include directory
# LoopTools_LIBRARY       - LoopTools library
#
# and defines the following imported targets:
#
# LoopTools::LoopTools

# search first in ${LoopTools_BUILD_DIR}
find_path(LoopTools_INCLUDE_DIR
  NAMES looptools.h
  PATHS
    ${LoopTools_BUILD_DIR}
  PATH_SUFFIXES
    build
    build-quad
    build-quad-alpha
    build-quad-gfortran
    build-quad-gfortran10
    build-quad-ifort
    build-xlf
  NO_DEFAULT_PATH
  NO_CMAKE_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  NO_CMAKE_FIND_ROOT_PATH
)

find_path(LoopTools_INCLUDE_DIR
  NAMES looptools.h
)

# search first in ${LoopTools_BUILD_DIR}
find_library(LoopTools_LIBRARY
  NAMES ooptools
  PATHS
    ${LoopTools_BUILD_DIR}
  PATH_SUFFIXES
    build
    build-quad
    build-quad-alpha
    build-quad-gfortran
    build-quad-gfortran10
    build-quad-ifort
    build-xlf
  NO_DEFAULT_PATH
  NO_CMAKE_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  NO_CMAKE_FIND_ROOT_PATH
)

find_library(LoopTools_LIBRARY
  NAMES ooptools
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(LoopTools
  FOUND_VAR LoopTools_FOUND
  REQUIRED_VARS
    LoopTools_LIBRARY
    LoopTools_INCLUDE_DIR
)

if(LoopTools_FOUND AND NOT TARGET LoopTools::LoopTools)
  add_library(LoopTools::LoopTools UNKNOWN IMPORTED)
  set_target_properties(LoopTools::LoopTools PROPERTIES
    IMPORTED_LOCATION "${LoopTools_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${LoopTools_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(
  LoopTools_INCLUDE_DIR
  LoopTools_LIBRARY
)
