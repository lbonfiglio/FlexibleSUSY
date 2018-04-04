# Distributed under the GPLv3 License.
# Author: Alexander Voigt
#
# FindSARAH
# ---------
#
# Finds the SARAH package.
# [http://sarah.hepforge.org]
#
# This module defines the following variables:
#
# SARAH_FOUND         - set if SARAH has been found
# SARAH_VERSION       - SARAH version

# set(Mathematica_DEBUG ON)

Mathematica_EXECUTE(
  CODE "Needs[\"SARAH\`\"]; If[!ValueQ[SA\`Version] || !StringQ[$sarahDir], WriteString[{\"stderr\"}, \"?.?.?\"], WriteString[{\"stderr\"}, SA\`Version]];"
  OUTPUT_VARIABLE SARAH_OUTPUT
  ERROR_VARIABLE SARAH_VERSION
  KERNEL_OPTIONS "-noprompt"
  TIMEOUT 30
)

if(${SARAH_VERSION} MATCHES "[0-9]+\.[0-9]+\.[0-9]+")
  set(SARAH_FOUND TRUE)
else()
  set(SARAH_FOUND FALSE)
  if(Mathematica_DEBUG)
    message(STATUS "The script to determine the SARAH version produced the following output:")
    message(STATUS "${SARAH_OUTPUT}")
  endif()
endif()

if(NOT SARAH_FOUND)
  message(FATAL_ERROR "${Red}Could not determine SARAH version.  "
    "Please make sure SARAH is installed on your system and can be loaded via Needs[\"SARAH`\"].  "
    "The install-sarah script can be used to install SARAH on unix-like systems:\n"
    "   ./install-sarah\n"
    "Please see `./install-sarah --help' for more information.  "
    "You can also re-run cmake with -DMathematica_DEBUG=ON to show some debug output.${ColourReset}")
endif()

if(SARAH_FIND_VERSION)
    if(SARAH_FIND_VERSION_EXACT AND NOT ${SARAH_VERSION} VERSION_EQUAL ${SARAH_FIND_VERSION})
      message(FATAL_ERROR "SARAH version ${SARAH_VERSION} found, "
        "but exact version ${SARAH_FIND_VERSION} is required.")
    elseif(${SARAH_VERSION} VERSION_LESS ${SARAH_FIND_VERSION})
      message(FATAL_ERROR "SARAH version ${SARAH_VERSION} found, "
        "but at least version ${SARAH_FIND_VERSION} is required.")
    endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(SARAH
  FOUND_VAR SARAH_FOUND
  REQUIRED_VARS
    SARAH_VERSION
)

mark_as_advanced(SARAH_OUTPUT)
