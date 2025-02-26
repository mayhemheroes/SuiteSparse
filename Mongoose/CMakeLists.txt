#-------------------------------------------------------------------------------
# Mongoose/CMakeLists.txt: instructions for building Mongoose with cmake
#-------------------------------------------------------------------------------
#
# Mongoose, Copyright (c) 2018-2022, All Rights Reserved.
#   Nuri Yeralan, Microsoft Research
#   Scott Kolodziej, Texas A&M University
#   Tim Davis, Texas A&M University
#   William Hager, University of Florida.
#
# SPDX-License-Identifier: GPL-3.0-only
#
# http://suitesparse.com   See Mongoose/Doc/License.txt for license.
#
# To compile the Mongoose library and demo programs just do "make"
# in this directory.  The Mongoose/Makefile will then run cmake
# in the build directory.
#
# To install Mongoose in /usr/local/lib and /usr/local/include,
# and also in ../lib, ../include, and ../bin:
#
#   sudo make install
#
# To install only in ../lib and ../include (no sudo required):
#
#   make local ; make install
#
# To remove all compiled files and libraries (except installed ones):
#
#   make distclean

cmake_minimum_required ( VERSION 3.19 )

#-------------------------------------------------------------------------------
# SuiteSparse policies
#-------------------------------------------------------------------------------

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake_modules
    ${CMAKE_SOURCE_DIR}/../SuiteSparse_config/cmake_modules )

# always build the static libraries (needed by the mongoose program)
set ( NSTATIC false )

include ( SuiteSparsePolicy )

#-------------------------------------------------------------------------------

set(Mongoose_DATE "Dec 9, 2022")
set(Mongoose_NUMERIC_DATE "2022-12-09")
set(Mongoose_VERSION_MAJOR 3)
set(Mongoose_VERSION_MINOR 0)
set(Mongoose_VERSION_PATCH 2)

project(Mongoose
        VERSION "${Mongoose_VERSION_MAJOR}.${Mongoose_VERSION_MINOR}.${Mongoose_VERSION_PATCH}"
        LANGUAGES CXX C)

# Mongoose installation location
include ( GNUInstallDirs )

# Configure files to pass CMake information
configure_file (
        "Version/Mongoose_Version.hpp.in"
        "${PROJECT_SOURCE_DIR}/Include/Mongoose_Version.hpp"
)
configure_file (
        "Version/title-info.tex.in"
        "${PROJECT_SOURCE_DIR}/Doc/title-info.tex"
)
configure_file (
        "Version/codemeta.json.in"
        "${PROJECT_SOURCE_DIR}/codemeta.json"
)

find_package ( SuiteSparse_config 6.0.2 REQUIRED )

include_directories("${PROJECT_BINARY_DIR}")

set(MONGOOSE_FILES
        Include/Mongoose_BoundaryHeap.hpp
        Include/Mongoose_Coarsening.hpp
        Include/Mongoose_CSparse.hpp
        Include/Mongoose_CutCost.hpp
        Include/Mongoose_Debug.hpp
        Include/Mongoose_EdgeCutOptions.hpp
        Include/Mongoose_EdgeCutProblem.hpp
        Include/Mongoose_EdgeCut.hpp
        Include/Mongoose_Graph.hpp
        Include/Mongoose_GuessCut.hpp
        Include/Mongoose_ImproveFM.hpp
        Include/Mongoose_ImproveQP.hpp
        Include/Mongoose_Internal.hpp
        Include/Mongoose_IO.hpp
        Include/Mongoose_Logger.hpp
        Include/Mongoose_Matching.hpp
        Include/Mongoose_Random.hpp
        Include/Mongoose_Refinement.hpp
        Include/Mongoose_Sanitize.hpp
        Include/Mongoose_Version.hpp
        Include/Mongoose_Waterdance.hpp
        Source/Mongoose_BoundaryHeap.cpp
        Source/Mongoose_Coarsening.cpp
        Source/Mongoose_CSparse.cpp
        Source/Mongoose_Debug.cpp
        Source/Mongoose_EdgeCut.cpp
        Source/Mongoose_Graph.cpp
        Source/Mongoose_GuessCut.cpp
        Source/Mongoose_ImproveFM.cpp
        Source/Mongoose_ImproveQP.cpp
        Source/Mongoose_IO.cpp
        Source/Mongoose_Logger.cpp
        Source/Mongoose_Matching.cpp
        Source/Mongoose_EdgeCutOptions.cpp
        Source/Mongoose_EdgeCutProblem.cpp
        Source/Mongoose_EdgeCut.cpp
        Source/Mongoose_Random.cpp
        Source/Mongoose_Refinement.cpp
        Source/Mongoose_Sanitize.cpp
        Source/Mongoose_Version.cpp
        Source/Mongoose_Waterdance.cpp
        )

set(NAPHEAP_FILES
        Include/Mongoose_QPBoundary.hpp
        Include/Mongoose_QPDelta.hpp
        Include/Mongoose_QPGradProj.hpp
        Include/Mongoose_QPLinks.hpp
        Include/Mongoose_QPMaxHeap.hpp
        Include/Mongoose_QPMinHeap.hpp
        Include/Mongoose_QPNapDown.hpp
        Include/Mongoose_QPNapsack.hpp
        Include/Mongoose_QPNapUp.hpp
        Source/Mongoose_QPBoundary.cpp
        Source/Mongoose_QPDelta.cpp
        Source/Mongoose_QPGradProj.cpp
        Source/Mongoose_QPLinks.cpp
        Source/Mongoose_QPMaxHeap.cpp
        Source/Mongoose_QPMinHeap.cpp
        Source/Mongoose_QPNapDown.cpp
        Source/Mongoose_QPNapsack.cpp
        Source/Mongoose_QPNapUp.cpp
        )

set(MMIO_FILES
        External/mmio/Source/mmio.c
        External/mmio/Include/mmio.h)

set(MONGOOSE_LIB_FILES
        ${MONGOOSE_FILES}
        ${NAPHEAP_FILES}
        ${MMIO_FILES})

set(EXE_FILES
        Executable/mongoose.cpp)

set(DEMO_FILES
        Demo/demo.cpp
        Include/Mongoose.hpp)

# Specify the Include directory for *.hpp files
include_directories(Include
        External/mmio/Include
        )

# Require C++11... Someday.
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories ( ${SUITESPARSE_CONFIG_INCLUDE_DIR} )

# set the output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Build the Mongoose library
add_library(mongoose_lib ${MONGOOSE_LIB_FILES})
set_property(TARGET mongoose_lib PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(mongoose_lib
        PROPERTIES OUTPUT_NAME mongoose)

target_link_libraries(mongoose_lib ${SUITESPARSE_CONFIG_LIBRARY})

if (UNIX AND NOT APPLE)
    target_link_libraries(mongoose_lib rt)
endif ()

# Build the Mongoose library for dynamic linking
add_library(mongoose_dylib SHARED ${MONGOOSE_LIB_FILES})
set_property(TARGET mongoose_dylib PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(mongoose_dylib
        PROPERTIES OUTPUT_NAME mongoose)

target_link_libraries(mongoose_dylib ${SUITESPARSE_CONFIG_LIBRARY})

# if(CMAKE_MAJOR_VERSION GREATER 2)
#     set_target_properties(mongoose_dylib PROPERTIES VERSION ${PROJECT_VERSION})
#     set_target_properties(mongoose_dylib PROPERTIES SOVERSION ${Mongoose_VERSION_MAJOR})
# else()
    set_target_properties(mongoose_dylib PROPERTIES VERSION 
        ${Mongoose_VERSION_MAJOR}.${Mongoose_VERSION_MINOR}.${Mongoose_VERSION_PATCH})
    set_target_properties(mongoose_dylib PROPERTIES SOVERSION
        ${Mongoose_VERSION_MAJOR})
# endif()

set_target_properties(mongoose_dylib PROPERTIES PUBLIC_HEADER Include/Mongoose.hpp)
target_include_directories(mongoose_dylib PRIVATE .)

if (UNIX AND NOT APPLE)
    target_link_libraries(mongoose_dylib rt)
endif ()

#-------------------------------------------------------------------------------

# Build the Mongoose library and all tests with coverage flags
add_library(mongoose_lib_dbg ${MONGOOSE_LIB_FILES})
set_target_properties(mongoose_lib_dbg
        PROPERTIES OUTPUT_NAME mongoose_dbg)

target_link_libraries(mongoose_lib_dbg ${SUITESPARSE_CONFIG_LIBRARY})

if (UNIX AND NOT APPLE)
    target_link_libraries(mongoose_lib_dbg rt)
endif ()

# Build the Mongoose executable
add_executable(mongoose_exe ${EXE_FILES})
set_target_properties(mongoose_exe
        PROPERTIES OUTPUT_NAME mongoose)
target_link_libraries(mongoose_exe mongoose_lib)

# Build the Demo executable
add_executable(demo_exe ${DEMO_FILES})
set_target_properties(demo_exe
        PROPERTIES OUTPUT_NAME demo)
target_link_libraries(demo_exe mongoose_lib)

# Coverage and Unit Testing Setup
enable_testing()
set(TESTING_OUTPUT_PATH ${CMAKE_BINARY_DIR}/tests)

# I/O Tests
add_executable(mongoose_test_io
        Tests/Mongoose_Test_IO.cpp
        Tests/Mongoose_Test_IO_exe.cpp)
target_link_libraries(mongoose_test_io mongoose_lib_dbg)
set_target_properties(mongoose_test_io PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TESTING_OUTPUT_PATH})

add_test(IO_Test ./runTests -min 1 -max 15 -t io -k)

# Edge Separator Tests
add_executable(mongoose_test_edgesep
        Tests/Mongoose_Test_EdgeSeparator.cpp
        Tests/Mongoose_Test_EdgeSeparator_exe.cpp)
target_link_libraries(mongoose_test_edgesep mongoose_lib_dbg)
set_target_properties(mongoose_test_edgesep PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TESTING_OUTPUT_PATH})

add_test(Edge_Separator_Test ./runTests -min 1 -max 15 -t edgesep)
add_test(Edge_Separator_Test_2 ./runTests -t edgesep -i 21 39 191 1557 1562 353 2468 1470 1380 505 182 201 2331 760 1389 2401 2420 242 250 1530 1533 360 1437)
add_test(Weighted_Edge_Separator_Test ./runTests -t edgesep -i 2624)
add_test(Target_Split_Test ./runTests -min 1 -max 15 -t edgesep -s 0.3)

# Memory Tests
add_executable(mongoose_test_memory
        Tests/Mongoose_Test_Memory.cpp
        Tests/Mongoose_Test_Memory_exe.cpp)
target_link_libraries(mongoose_test_memory mongoose_lib_dbg)
set_target_properties(mongoose_test_memory PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TESTING_OUTPUT_PATH})

add_test(Memory_Test ./runTests -min 1 -max 15 -t memory)

# Performance Test
add_executable(mongoose_test_performance
        Tests/Mongoose_Test_Performance.cpp
        Tests/Mongoose_Test_Performance_exe.cpp)
target_link_libraries(mongoose_test_performance mongoose_lib)
set_target_properties(mongoose_test_performance PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TESTING_OUTPUT_PATH})
add_test(Performance_Test ./runTests -min 1 -max 15 -t performance -p)
add_test(Performance_Test_2 ./runTests -t performance -i 21 39 1557 1562 353 2468 1470 1380 505 182 201 2331 760 1389 2401 2420 242 250 1530 1533 -p)

# Reference Test
add_executable(mongoose_test_reference
        Tests/Mongoose_Test_Reference.cpp
        Tests/Mongoose_Test_Reference_exe.cpp)
target_link_libraries(mongoose_test_reference mongoose_lib)
set_target_properties(mongoose_test_reference PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TESTING_OUTPUT_PATH})

# Unit Tests
add_executable(mongoose_unit_test_io
        Tests/Mongoose_UnitTest_IO_exe.cpp)
target_link_libraries(mongoose_unit_test_io mongoose_lib_dbg)
set_target_properties(mongoose_unit_test_io PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TESTING_OUTPUT_PATH})
add_test(Unit_Test_IO ./tests/mongoose_unit_test_io)

add_executable(mongoose_unit_test_graph
        Tests/Mongoose_UnitTest_Graph_exe.cpp)
target_link_libraries(mongoose_unit_test_graph mongoose_lib_dbg)
set_target_properties(mongoose_unit_test_graph PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TESTING_OUTPUT_PATH})
add_test(Unit_Test_Graph ./tests/mongoose_unit_test_graph)

add_executable(mongoose_unit_test_edgesep
        Tests/Mongoose_UnitTest_EdgeSep_exe.cpp)
target_link_libraries(mongoose_unit_test_edgesep mongoose_lib_dbg)
set_target_properties(mongoose_unit_test_edgesep PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TESTING_OUTPUT_PATH})
add_test(Unit_Test_EdgeSep ./tests/mongoose_unit_test_edgesep)

message(STATUS "CMAKE_CXX_COMPILER: " ${BoldBlue} ${CMAKE_CXX_COMPILER_ID} ${ColourReset})
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    # using Clang
    SET(CMAKE_CXX_FLAGS "-O3 -fwrapv")
    # Debug flags for Clang
    SET(CMAKE_CXX_FLAGS_DEBUG "--coverage -g -fwrapv")
    SET(CMAKE_C_FLAGS_DEBUG "--coverage -g")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "--coverage -g")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    SET(CMAKE_CXX_FLAGS "-O3 -fwrapv")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-but-set-variable -Wno-unused-variable" )
    # Debug flags for GCC
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
        message(WARNING "${BoldRed} WARNING:${ColourReset} Your compiler does not support the C++11 ISO standard. Future versions of Mongoose will require a compiler with C++11 support. We recommend you upgrade to at least GCC 4.6.")
    endif ()
    SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -fwrapv -W -Wshadow -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
    SET(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall -fwrapv -W -fprofile-arcs -ftest-coverage")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "-fprofile-arcs -ftest-coverage")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    SET(CMAKE_CXX_FLAGS "-O3 -no-prec-div -xHOST -ipo -fwrapv")
    # Debug flags for Intel
    SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -fwrapv")
    SET(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall")
#elseif ( MSVC )
    # using Visual Studio C++
endif ()

set_target_properties(mongoose_lib_dbg PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
set_target_properties(mongoose_lib_dbg PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")

# Add debug compile/linker flags
set_target_properties(mongoose_test_io PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
set_target_properties(mongoose_test_io PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
set_target_properties(mongoose_test_memory PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
set_target_properties(mongoose_test_memory PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
set_target_properties(mongoose_test_edgesep PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
set_target_properties(mongoose_test_edgesep PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
set_target_properties(mongoose_unit_test_io PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
set_target_properties(mongoose_unit_test_io PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
set_target_properties(mongoose_unit_test_graph PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
set_target_properties(mongoose_unit_test_graph PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
set_target_properties(mongoose_unit_test_edgesep PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
set_target_properties(mongoose_unit_test_edgesep PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")

set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1) # Necessary for gcov - prevents file.cpp.gcda instead of file.gcda

# Copy over runTest.py to build folder for ease of use
file(COPY Tests/runTests DESTINATION ${CMAKE_BINARY_DIR})

add_custom_command(TARGET mongoose_lib
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:mongoose_lib> ${PROJECT_SOURCE_DIR}/Lib
                   COMMENT "Copying libmongoose (static) to root Lib directory"
        )

add_custom_command(TARGET mongoose_dylib
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:mongoose_dylib> ${PROJECT_SOURCE_DIR}/Lib
        COMMENT "Copying libmongoose (dynamic) to root Lib directory"
        )

add_custom_target(purge
        COMMAND rm -rf ${CMAKE_BINARY_DIR}/*
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Lib/libmongoose.*
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Doc/title-info.tex
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Doc/title-info.tex.aux
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Doc/Mongoose_UserGuide.idx
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Doc/Mongoose_UserGuide.log
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Doc/Mongoose_UserGuide.out
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Doc/Mongoose_UserGuide.aux
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Doc/Mongoose_UserGuide.toc
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Doc/Mongoose_UserGuide.bbl
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Doc/Mongoose_UserGuide.blg
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Matrix/*.tar.gz
        COMMAND rm -f ${PROJECT_SOURCE_DIR}/Matrix/*.csv
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )

add_custom_target(userguide
        COMMAND make
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Doc
        )

#-------------------------------------------------------------------------------
# Mongoose installation location
#-------------------------------------------------------------------------------

if ( GLOBAL_INSTALL )
    # install in /usr/local/lib, /usr/local/include, and /usr/local/bin.
    # requires "sudo make install"
    message ( STATUS "Installation will be system-wide (requires 'sudo make install')." )
    message ( STATUS "If you do not have 'sudo' privilege, do 'make local' instead." )
    install ( TARGETS mongoose_dylib
        LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )
    install ( TARGETS mongoose_lib
        ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR} )
    install ( TARGETS mongoose_exe
        RUNTIME       DESTINATION ${CMAKE_INSTALL_BINDIR} )
    install ( FILES ${CMAKE_SOURCE_DIR}/cmake_modules/FindMongoose.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SuiteSparse 
        COMPONENT Development )
endif ( )

if ( INSIDE_SUITESPARSE )
    # also install in SuiteSparse/lib and SuiteSparse/include;
    # does not require "sudo make install", just "make install"
    message ( STATUS "Installation in ../lib and ../include," )
    message ( STATUS "  with 'make local ; make install'. No 'sudo' required." )
    install ( TARGETS mongoose_dylib
        LIBRARY       DESTINATION ${SUITESPARSE_LIBDIR}
        PUBLIC_HEADER DESTINATION ${SUITESPARSE_INCLUDEDIR} )
    install ( TARGETS mongoose_lib
        ARCHIVE       DESTINATION ${SUITESPARSE_LIBDIR} )
    install ( TARGETS mongoose_exe
        RUNTIME       DESTINATION ${SUITESPARSE_BINDIR} )
    install ( FILES ${CMAKE_SOURCE_DIR}/cmake_modules/FindMongoose.cmake
        DESTINATION ${SUITESPARSE_LIBDIR}/cmake/SuiteSparse 
        COMPONENT Development )
endif ( )

#-------------------------------------------------------------------------------
# report status
#-------------------------------------------------------------------------------

include ( SuiteSparseReport )

