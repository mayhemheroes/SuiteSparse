#-------------------------------------------------------------------------------
# SuiteSparse/SuiteSparse_config/cmake_modules/FindAMD.cmake
#-------------------------------------------------------------------------------

# Copyright (c) 2022, Timothy A. Davis.  All Rights Reserved.
# SPDX-License-Identifier: BSD-3-clause

# Finds the AMD include file and compiled library and sets:

# AMD_INCLUDE_DIR - where to find amd.h
# AMD_LIBRARY     - compiled AMD library
# AMD_LIBRARIES   - libraries when using AMD
# AMD_FOUND       - true if AMD found

# set ``AMD_ROOT`` to an AMD installation root to
# tell this module where to look.

# To use this file in your application, copy this file into MyApp/cmake_modules
# where MyApp is your application and add the following to your
# MyApp/CMakeLists.txt file:
#
#   set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")
#
# or, assuming MyApp and SuiteSparse sit side-by-side in a common folder, you
# can leave this file in place and use this command (revise as needed):
#
#   set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
#       "${CMAKE_SOURCE_DIR}/../SuiteSparse/SuiteSparse_config/cmake_modules")

#-------------------------------------------------------------------------------

# include files for AMD
find_path ( AMD_INCLUDE_DIR
    NAMES amd.h
    HINTS ${CMAKE_SOURCE_DIR}/..
    HINTS ${CMAKE_SOURCE_DIR}/../SuiteSparse/AMD
    HINTS ${CMAKE_SOURCE_DIR}/../AMD
    PATH_SUFFIXES include Include
)

# compiled libraries AMD
find_library ( AMD_LIBRARY
    NAMES amd
    HINTS ${CMAKE_SOURCE_DIR}/..
    HINTS ${CMAKE_SOURCE_DIR}/../SuiteSparse/AMD
    HINTS ${CMAKE_SOURCE_DIR}/../AMD
    PATH_SUFFIXES lib build
)

# get version of the library
foreach (_VERSION MAIN_VERSION SUB_VERSION SUBSUB_VERSION)
  file (STRINGS ${AMD_INCLUDE_DIR}/amd.h _VERSION_LINE REGEX "define[ ]+AMD_${_VERSION}")
  if (_VERSION_LINE)
    string (REGEX REPLACE ".*define[ ]+AMD_${_VERSION}[ ]+([0-9]*).*" "\\1" _AMD_${_VERSION} "${_VERSION_LINE}")
  endif ()
  unset (_VERSION_LINE)
endforeach ()
set (AMD_VERSION "${_AMD_MAIN_VERSION}.${_AMD_SUB_VERSION}.${_AMD_SUBSUB_VERSION}")
set (AMD_LIBRARIES ${AMD_LIBRARY})

include (FindPackageHandleStandardArgs)

find_package_handle_standard_args ( AMD
    REQUIRED_VARS AMD_LIBRARIES AMD_INCLUDE_DIR
    VERSION_VAR AMD_VERSION
)

mark_as_advanced (
    AMD_INCLUDE_DIR
    AMD_LIBRARY
    AMD_LIBRARIES
)

if ( AMD_FOUND )
    message ( STATUS "AMD include dir: ${AMD_INCLUDE_DIR}")
    message ( STATUS "AMD library:     ${AMD_LIBRARY}")
    message ( STATUS "AMD version:     ${AMD_VERSION}" )
else ( )
    message ( STATUS "AMD not found" )
endif ( )

