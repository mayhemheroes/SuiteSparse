#-------------------------------------------------------------------------------
# SuiteSparse/SPEX/Makefile
#-------------------------------------------------------------------------------

# SPEX_Left_LU: (c) 2019-2022, Chris Lourenco (US Naval Academy), Jinhao Chen,
# Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

# A simple Makefile for SPEX, which relies on cmake to do the
# actual build.  All the work is done in cmake so this Makefile is just for
# convenience.

# To compile with an alternate compiler:
#
#       make CC=gcc CXX=g++
#
# To compile/install for system-wide usage:
#
#       make
#       sudo make install
#
# To compile/install for local usage (SuiteSparse/lib and SuiteSparse/include):
#
#       make local
#       make install
#
# To clean up the files:
#
#       make clean

JOBS ?= 8

default: library

# default is to install only in /usr/local
library:
	( cd build && cmake $(CMAKE_OPTIONS) .. && cmake --build . -j${JOBS} )

# install only in SuiteSparse/lib and SuiteSparse/include
local:
	( cd build && cmake $(CMAKE_OPTIONS) -DGLOBAL_INSTALL=0 -DLOCAL_INSTALL=1 .. && cmake --build . -j${JOBS} )

# install only in /usr/local (default)
global:
	( cd build && cmake $(CMAKE_OPTIONS) -DGLOBAL_INSTALL=1 -DLOCAL_INSTALL=0 .. && cmake --build . -j${JOBS} )

# install in SuiteSparse/lib both and SuiteSparse/include and /usr/local
both:
	( cd build && cmake $(CMAKE_OPTIONS) -DGLOBAL_INSTALL=1 -DLOCAL_INSTALL=1 .. && cmake --build . -j${JOBS} )

debug:
	( cd build && cmake $(CMAKE_OPTIONS) -DCMAKE_BUILD_TYPE=Debug .. && cmake --build . -j${JOBS} )

all: library

demos: library
	( cd build && cmake $(CMAKE_OPTIONS) -DDEMO=1 .. && cmake --build . -j${JOBS} )
	( cd SPEX_Left_LU/Demo && ../../build/example  )
	( cd SPEX_Left_LU/Demo && ../../build/example2 )
	( cd SPEX_Left_LU/Demo && ../../build/spexlu_demo )

cov:
	( cd SPEX_Left_LU/Tcov && $(MAKE) )

# just compile after running cmake; do not run cmake again
remake:
	( cd build && cmake --build . -j${JOBS} )

# just run cmake to set things up
setup:
	( cd build && cmake $(CMAKE_OPTIONS) .. )

install:
	( cd build && cmake --install . )

# remove any installed libraries and #include files
uninstall:
	- xargs rm < build/install_manifest.txt

# remove all files not in the distribution
clean:
	- $(RM) -rf build/* Config/*.tmp MATLAB/*.o MATLAB/*.mex*
	( cd SPEX_Left_LU/Tcov && $(MAKE) clean )

purge: clean

distclean: clean

docs:
	( cd Doc && $(MAKE) )
