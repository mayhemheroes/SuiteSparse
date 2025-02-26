#-------------------------------------------------------------------------------
# SuiteSparse/UMFPACK/Tcov/Makefile
#-------------------------------------------------------------------------------

# Copyright (c) 1995-2022, Timothy A. Davis.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0+

#-------------------------------------------------------------------------------

# statement coverage test

# choose the test options:  also see ../Source/umf_config.h
#
# NBLAS: do not use the BLAS; use plain C code instead
# NRECIPROCAL: do not multiply by the reciprocal when scaling
# TEST_FOR_INTEGER_OVERFLOW: set the max mark to 3n, to test clearing mark
# NO_DIVIDE_BY_ZERO: if true, never divide by zero
#
# OPTION = 
# OPTION = -DNBLAS -DNRECIPROCAL
  OPTION = -DNBLAS -DNRECIPROCAL -DTEST_FOR_INTEGER_OVERFLOW
# OPTION = -DNBLAS -DNRECIPROCAL -DTEST_FOR_INTEGER_OVERFLOW -DNO_DIVIDE_BY_ZERO

SUITESPARSE ?= $(realpath $(CURDIR)/../..)

CF = $(OPTION) -g --coverage -fprofile-abs-path -std=c11 -DTESTING -fopenmp \
    -DBLAS32
I = -I../../include -I../Include -I../../UMFPACK/Source -I../../AMD/Include

BLAS = -lblas
CLIB = -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd
FLIB = $(LAPACK) $(BLAS)

# LIBS = $(SUITESPARSE)/SuiteSparse_config/build/libsuitesparseconfig.a \
#     $(SUITESPARSE)/CHOLMOD/build/libcholmod.a \
#     $(SUITESPARSE)/AMD/build/libamd.a \
#     $(SUITESPARSE)/CAMD/build/libcamd.a \
#     $(SUITESPARSE)/COLAMD/build/libcolamd.a \
#     $(SUITESPARSE)/CCOLAMD/build/libccolamd.a \
#     -lblas -lm

LIBS = -L../../lib $(FLIB) $(CLIB) \
    -Wl,--rpath=$(SUITESPARSE)/lib -lm -lrt

SRC = ../Source2/*.c

SRC2 = $(notdir $(wildcard $(SRC)))

OBJ = $(SRC2:.c=.o)

default: $(OBJ) ut.c
	gcc ut_di.c $(CF) $(I) $(OBJ) $(LIBS) -o ut_di
	gcc ut_dl.c $(CF) $(I) $(OBJ) $(LIBS) -o ut_dl
	gcc ut_zi.c $(CF) $(I) $(OBJ) $(LIBS) -o ut_zi
	gcc ut_zl.c $(CF) $(I) $(OBJ) $(LIBS) -o ut_zl
	./ut_di > ut_di.out
	./ut_dl > ut_dl.out
	./ut_zi > ut_zi.out
	./ut_zl > ut_zl.out
	- ./gcovs umf*.c
	- ./cover umf*v > cover.out
	- echo -n "statments not yet tested: "
	- grep -c "#####" cover.out

.PRECIOUS: %.c

%.o: ../Source2/%.c
	- ln -s $< $(basename $@).c
	gcc $(CF) $(I) -c $(basename $@).c -o $(basename $@).o

purge:
	$(RM) umf*.c *.o ut_?? *.gcno *.gcda *.gcov debug.amd *.umf A500 *.out

