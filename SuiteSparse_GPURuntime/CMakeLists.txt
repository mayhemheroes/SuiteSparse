#-------------------------------------------------------------------------------
# SuiteSparse/SuiteSparse_GPURuntime/CMakeLists.txt
#-------------------------------------------------------------------------------

# SuiteSparse_GPURuntime, Copyright (c) 2013-2022, Timothy A Davis,
# Sencer Nuri Yeralan, and Sanjay Ranka.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0+

#-------------------------------------------------------------------------------
# get the version
#-------------------------------------------------------------------------------

cmake_minimum_required ( VERSION 3.19 )

set ( SUITESPARSE_GPURUNTIME_DATE "Dec 9, 2022" )
set ( SUITESPARSE_GPURUNTIME_VERSION_MAJOR 2 )
set ( SUITESPARSE_GPURUNTIME_VERSION_MINOR 0 )
set ( SUITESPARSE_GPURUNTIME_VERSION_SUB   2 )

message ( STATUS "Building SUITESPARSE_GPURUNTIME version: v"
    ${SUITESPARSE_GPURUNTIME_VERSION_MAJOR}.
    ${SUITESPARSE_GPURUNTIME_VERSION_MINOR}.
    ${SUITESPARSE_GPURUNTIME_VERSION_SUB} " (" ${SUITESPARSE_GPURUNTIME_DATE} ")" )

#-------------------------------------------------------------------------------
# SuiteSparse policies
#-------------------------------------------------------------------------------

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake_modules
    ${CMAKE_SOURCE_DIR}/../SuiteSparse_config/cmake_modules )

option ( ENABLE_CUDA "Enable CUDA acceleration" on )

include ( SuiteSparsePolicy )

#-------------------------------------------------------------------------------
# define the project
#-------------------------------------------------------------------------------

if ( SUITESPARSE_CUDA )
    project ( suitesparse_gpuruntime 
        VERSION "${SUITESPARSE_GPURUNTIME_VERSION_MAJOR}.${SUITESPARSE_GPURUNTIME_VERSION_MINOR}.${SUITESPARSE_GPURUNTIME_VERSION_SUB}"
        LANGUAGES C CXX CUDA )
else ( )
    # CUDA not available
    project ( suitesparse_gpuruntime 
        VERSION "${SUITESPARSE_GPURUNTIME_VERSION_MAJOR}.${SUITESPARSE_GPURUNTIME_VERSION_MINOR}.${SUITESPARSE_GPURUNTIME_VERSION_SUB}"
        LANGUAGES C CXX )
    message ( STATUS "CUDA not enabled:  library will be empty" )
endif ( )

#-------------------------------------------------------------------------------
# find library dependencies
#-------------------------------------------------------------------------------

find_package ( SuiteSparse_config 6.0.2 REQUIRED )

#-------------------------------------------------------------------------------
# configure files
#-------------------------------------------------------------------------------

configure_file ( "Config/SuiteSparse_GPURuntime.hpp.in" "${PROJECT_SOURCE_DIR}/Include/SuiteSparse_GPURuntime.hpp")

#-------------------------------------------------------------------------------

if ( SUITESPARSE_CUDA )
    set ( CMAKE_CUDA_FLAGS "-cudart=static -lineinfo -DSUITESPARSE_CUDA" )
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSUITESPARSE_CUDA" )
    message ( STATUS "C++ flags for CUDA:  ${CMAKE_CXX_FLAGS}" )
    message ( STATUS "nvcc flags for CUDA: ${CMAKE_CUDA_FLAGS}" )
endif ( )

set ( CMAKE_CXX_STANDARD 17 )
set ( CMAKE_C_STANDARD 11 )

file ( GLOB SUITESPARSE_GPURUNTIME_SOURCES "Source/*.cpp" )

set ( SUITESPARSE_GPURUNTIME_INCLUDES Include 
    ${SUITESPARSE_CONFIG_INCLUDE_DIR} )

include_directories ( ${SUITESPARSE_GPURUNTIME_INCLUDES} ${CUDAToolkit_INCLUDE_DIRS} )

#-------------------------------------------------------------------------------
# dynamic suitesparse_gpuruntime library properties
#-------------------------------------------------------------------------------

add_library ( suitesparse_gpuruntime SHARED ${SUITESPARSE_GPURUNTIME_SOURCES} )

set_target_properties ( suitesparse_gpuruntime PROPERTIES
    VERSION ${SUITESPARSE_GPURUNTIME_VERSION_MAJOR}.${SUITESPARSE_GPURUNTIME_VERSION_MINOR}.${SUITESPARSE_GPURUNTIME_VERSION_SUB}
    CXX_STANDARD_REQUIRED 17
    SOVERSION ${SUITESPARSE_GPURUNTIME_VERSION_MAJOR}
    PUBLIC_HEADER "Include/SuiteSparse_GPURuntime.hpp" )

target_include_directories ( suitesparse_gpuruntime PUBLIC
        ${CUDAToolkit_INCLUDE_DIRS}
        ${SUITESPARSE_GPURUNTIME_INCLUDES} )

if ( SUITESPARSE_CUDA )
    set_target_properties ( suitesparse_gpuruntime PROPERTIES POSITION_INDEPENDENT_CODE ON )
    set_target_properties ( suitesparse_gpuruntime PROPERTIES CUDA_SEPARABLE_COMPILATION ON )
    target_link_libraries ( suitesparse_gpuruntime CUDA::nvrtc CUDA::cudart_static
        CUDA::nvToolsExt CUDA::cublas )
endif ( )

#-------------------------------------------------------------------------------
# static suitesparse_gpuruntime library properties
#-------------------------------------------------------------------------------

if ( NOT NSTATIC )
add_library ( suitesparse_gpuruntime_static STATIC ${SUITESPARSE_GPURUNTIME_SOURCES} )

set_target_properties ( suitesparse_gpuruntime_static PROPERTIES
    VERSION ${SUITESPARSE_GPURUNTIME_VERSION_MAJOR}.${SUITESPARSE_GPURUNTIME_VERSION_MINOR}.${SUITESPARSE_GPURUNTIME_VERSION_SUB}
    CXX_STANDARD_REQUIRED 17
    SOVERSION ${SUITESPARSE_GPURUNTIME_VERSION_MAJOR} )

target_include_directories ( suitesparse_gpuruntime_static PUBLIC
        ${CUDAToolkit_INCLUDE_DIRS}
        ${SUITESPARSE_GPURUNTIME_INCLUDES} )

if ( SUITESPARSE_CUDA )
    set_target_properties ( suitesparse_gpuruntime_static PROPERTIES CUDA_SEPARABLE_COMPILATION on )
    set_target_properties ( suitesparse_gpuruntime_static PROPERTIES POSITION_INDEPENDENT_CODE on )
    target_link_libraries ( suitesparse_gpuruntime_static CUDA::nvrtc CUDA::cudart_static
        CUDA::nvToolsExt CUDA::cublas )
endif ( )
endif ( )

#-------------------------------------------------------------------------------
# SuiteSparse_GPURuntime installation location
#-------------------------------------------------------------------------------

if ( GLOBAL_INSTALL )
    # install in /usr/local/lib and /usr/local/include.
    # requires "sudo make install"
    message ( STATUS "Installation will be system-wide (requires 'sudo make install')" )
    install ( TARGETS suitesparse_gpuruntime
        LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )
    install ( FILES
        ${CMAKE_SOURCE_DIR}/cmake_modules/FindSuiteSparse_GPURuntime.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SuiteSparse 
        COMPONENT Development )
    if ( NOT NSTATIC )
    install ( TARGETS suitesparse_gpuruntime_static
        ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR} )
    endif ( )
endif ( )

if ( INSIDE_SUITESPARSE )
    # also install in SuiteSparse/lib and SuiteSparse/include;
    # does not require "sudo make install", just "make install"
    message ( STATUS "Installation in ../../lib and ../../include," )
    message ( STATUS "  with 'make local ; make install'. No 'sudo' required." )
    install ( TARGETS suitesparse_gpuruntime
        LIBRARY       DESTINATION ${SUITESPARSE_LIBDIR}
        PUBLIC_HEADER DESTINATION ${SUITESPARSE_INCLUDEDIR} )
    install ( FILES
        ${CMAKE_SOURCE_DIR}/cmake_modules/FindSuiteSparse_GPURuntime.cmake
        DESTINATION ${SUITESPARSE_LIBDIR}/cmake/SuiteSparse 
        COMPONENT Development )
    if ( NOT NSTATIC )
    install ( TARGETS suitesparse_gpuruntime_static
        ARCHIVE       DESTINATION ${SUITESPARSE_LIBDIR} )
    endif ( )
endif ( )

#-------------------------------------------------------------------------------
# report status
#-------------------------------------------------------------------------------

include ( SuiteSparseReport )

