/*
 * ---------------------------------------------------------------------
 * Copyright (C) 2012 Tino Kluge (tino.kluge@hrz.tu-chemnitz.de)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "helper.h"


int main(int argc, char** argv) {
   // values needed for dgesv (see also man 1 dgesv)
   int n;		/* order of matrix A, i.e. number of equations */
   int lda;		/* leading dimension of A, should be n */
   int ldb;		/* leading dimension of b, should be n */
   int nrhs = 1;	/* number of right hand sides */
   double **A;		/* matrix A, want to solve Ax=b */
   double *b;
   double *x;
   int *ipiv;		/* permutation vector for LU decomposition */
   int info;		/* result of dgesv, success=0 */

   // values needed for dgemv
   char	trans='N';	/* is A transposed: 'N' no, 'T' transposed, 'C' conj */
   double   alpha=1.0;
   double   beta=0.0;
   int	incx=1,incy=1;

   // other values
   int i,j;
   double   norm=0.0;      // norm of the error
   double   t[5];          // timer: 0 alloc, 1 init, 2 A*x, 3 solve, 4 total
   int      csvout;

   t[4]=stoptime();

   // read input parameters
   if(argc<2) {
      printf("usage: %s <order of matrix> [csvout]\n", argv[0]);
      exit(EXIT_FAILURE);
   }
   if(argc==2){ csvout=0; } else { csvout=1; }
   n=atoi(argv[1]);
   if(n<=0) n=1;
   lda=n;
   ldb=n;

   // reserver memory
   t[0]=stoptime();
   A=dmatrix(n,n);
   b=dvector(n);
   x=dvector(n);
   ipiv=ivector(n);
   t[0]=stoptime()-t[0];


   // randomly initialise the matrix
   t[1]=stoptime();
   for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
         A[j][i]=randf();
      }
      x[i]=randf();
   }
   t[1]=stoptime()-t[1];

   /* b = 1.0*A*x + 0.0*b */
   t[2]=stoptime();
   dgemv_(&trans,&n,&n,&alpha,&A[0][0],&lda,&x[0],&incx,&beta,&b[0],&incy);
   t[2]=stoptime()-t[2];


   // solve the linear system
   // since fortran stores matrices in different order than C, we
   // provide the transposed matrix of A (see indexing erlier)
   if(!csvout) { printf("solving ... "); fflush(stdout); }
   t[3]=stoptime();
   dgesv_(&n, &nrhs, &A[0][0], &lda, &ipiv[0], &b[0], &ldb, &info);
   t[3]=stoptime()-t[3];
   if(!csvout) { printf("finished\n"); fflush(stdout); }
   // -------------------------------------------------------------

   // check for success
   if(info == 0) {
      // check the numerical error
      norm=0.0;
      for(i = 0; i < n; i++) {
         norm+=(x[i]-b[i])*(x[i]-b[i]);
      }
      norm=sqrt(norm);
      if(!csvout) { printf("The norm of the error is %e\n",norm); }
   } else {
      // write an error message
      fprintf(stderr,"dgesv returned error %i\n", info) ;
   }
   t[4]=stoptime()-t[4];

   print_timing(argv[0],t,n,norm,csvout);

   return info;
}
