CC=gcc
CPP=g++

# flags for naive/eigen where the core LU routines need to be compiled
CFLAGSEIGEN=-Wall -O3 -DNDEBUG 
CFLAGSNAIVE =-Wall -O3 -funroll-loops #-ffast-math -ftree-vectorizer-verbose=0
ARCH=-mfpmath=sse -march=native

# flags for blas/atlas/goto where the core LU routines are given by a library
CFLAGS=-O2 -Wall

PROGS=lu_goto lu_atlas lu_blas lu_naive lu_eigen lu_goto_pt lu_atlas_pt lu_eigen_pt 

all:		${PROGS}


lu_goto:	lu_blas.o
		${CC} lu_blas.o -lgoto2 -lgfortran -lm -o lu_goto

lu_goto_pt:	lu_blas.o
		${CC} lu_blas.o -lgoto2_pt -lgfortran -lm -o lu_goto_pt

lu_atlas:	lu_blas.o
		${CC} lu_blas.o -llapack -lcblas -lf77blas -latlas -lgfortran -lm -o lu_atlas

lu_atlas_pt:	lu_blas.o
		${CC} lu_blas.o ${LATLAS} -llapack -lptcblas -lptf77blas -latlas -lgfortran -lm -lpthread -o lu_atlas_pt

lu_blas:	lu_blas.o
		${CC} lu_blas.o -llapack -lblas -lm -o lu_blas



lu_naive:	lu_naive.c
		${CC} ${CFLAGSNAIVE} ${ARCH} lu_naive.c -lm -o lu_naive

lu_eigen:	lu_eigen.cpp
		${CPP} ${CFLAGSEIGEN} ${ARCH} lu_eigen.cpp -o lu_eigen

lu_eigen_pt:	lu_eigen.cpp
		${CPP} ${CFLAGSEIGEN} ${ARCH} -fopenmp  lu_eigen.cpp -o lu_eigen_pt



# ------------------------------------------
%.o: %.c
		${CC} ${CFLAGS} -c $< -o $@
%.o: %.cpp
		${CPP} ${CFLAGS} -c $< -o $@

clean:
		rm -f ${PROGS} a.out *.o
