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
#include <time.h>
#include "helper.h"

// randomly initialise the matrix
void rand_matrix(int n, double** A) {
   int	i,j;
   for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
         A[i][j]=randf();
      }
   }
}

void rand_vector(int n, double* x) {
   int	i;
   for(i=0; i<n; i++) {
      x[i]=randf();
   }
}

/* b=alpha*x+y */
void d_axpy(int n, double* b, double alpha, double* x, double* y) {
   int	i;
   for(i=0; i<n; i++) {
      b[i]=alpha*x[i]+y[i];
   }
}

/* y=A*x */
void d_Ax(int n, int m, double* y, double** A, double* x) {
   int	i,j;
   double	sum;

   for(i=0; i<n; i++) {
      sum=0;
      for(j=0; j<m; j++) {
         sum+=A[i][j]*x[j];
      }
      y[i]=sum;
   }
}

/* C=A*B */
void d_AB(int n, int m, double** C, double** A, double** B) {
   int	i,j,k;
   double	sum;

   for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
         sum=0;
         for(k=0; k<m; k++) {
            sum+=A[i][k]*B[k][j];
         }
         C[i][j]=sum;
      }
   }
}

/* A=L*U, the result is stored in A, where the diagonal belongs to U (diag(L)=1) */
/* A=L*U  <--> L^(-1)*A=U, so we do a Gauss elimination to obtain U */
void d_lu(int n, double** A) {
   int	i,j,k;
   double	factor;

   for(i=0; i<n; i++) {
      // use line i to eliminate values in lines i+1,...,n
      for(k=i+1; k<n; k++) {
         factor=-A[k][i]/A[i][i];
         // write into U
         for(j=i+1; j<n; j++) {
            A[k][j]+=factor*A[i][j];
         }
         // write into L
         A[k][i]=-factor;
      }
   }
}

/* solves Lx=b, assuming diag(L)=1 */
void sweep_forward(int n, double *x, double** A, double* b) {
   int	i,j;
   double	sum;

   for(i=0; i<n; i++) {
      sum=0;
      for(j=0; j<i; j++) {
         sum+=A[i][j]*x[j];
      }
      x[i]=b[i]-sum;
   }
}

/* solves Ux=b, assuming diag(L)!=1 */
void sweep_backward(int n, double *x, double** A, double* b) {
   int	i,j;
   double	sum;

   for(i=n-1; i>=0; i--) {
      sum=0;
      for(j=i+1; j<n; j++) {
         sum+=A[i][j]*x[j];
      }
      x[i]=(b[i]-sum)/A[i][i];
   }
}

/* solve equation system Ax=b */
void solve(int n, double* x, double** A, double* b) {
   double* y;
   // A=LU
   d_lu(n,A);
   // LUx=b <--> Ly=b and Ux=b
   y=dvector(n);
   sweep_forward (n,y,A,b);
   sweep_backward(n,x,A,y);
}

double norm2(int n, double *x) {
   int	i;
   double	sum=0;
   for(i=0; i<n; i++) {
      sum+=x[i]*x[i];
   }
   return sqrt(sum);
}

int main(int argc, char** argv) {
   // values needed for dgesv (see also man 1 dgesv)
   int n;			/* order of matrix A, i.e. number of equations */
   double **A;		/* matrix A, want to solve Ax=b */
   double *b;
   double *x, *x_sol, *x_err;

   // other values
   double norm=0.0;	/* norm of the error */
   double t[5];	        // timer: 0 alloc, 1 init, 2 A*x, 3 solve, 4 total
   int   csvout;

   t[4]=stoptime();

   // read input parameters
   if(argc<2) {
      printf("usage: %s <matrix order> [csvout]\n", argv[0]);
      exit(EXIT_FAILURE);
   }
   if(argc==2){ csvout=0; } else { csvout=1; }
   n=atoi(argv[1]);
   if(n<=0) n=1;

   // reserver memory
   t[0]=stoptime();
   A=dmatrix(n,n);
   b=dvector(n);
   x=dvector(n);
   x_sol=dvector(n);
   x_err=dvector(n);
   t[0]=stoptime()-t[0];

   // init random matrix
   t[1]=stoptime();
   srand(time(NULL));
   rand_matrix(n,A);
   rand_vector(n,x_sol);
   t[1]=stoptime()-t[1];

   /* b = Ax */
   t[2]=stoptime();
   d_Ax(n,n,b,A,x_sol);
   t[2]=stoptime()-t[2];


   // solve the linear system
   if(!csvout) { printf("solving ... "); fflush(stdout); }
   t[3]=stoptime();
   solve(n,x,A,b);
   t[3]=stoptime()-t[3];
   if(!csvout) { printf("finished\n"); fflush(stdout); }
   // ------------------------------------------------------------------------

   // check for success
   // check the numerical error
   // b=-x+b
   d_axpy(n,x_err,-1.0,x_sol, x);
   norm=norm2(n,x_err);
   if(!csvout) { printf("The norm of the error is %e\n",norm); }

   t[4]=stoptime()-t[4];


   print_timing(argv[0],t,n,norm,csvout);

   return EXIT_SUCCESS;
}

