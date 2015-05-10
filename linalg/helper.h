#include <sys/time.h>


double stoptime(void) {
 struct timeval	t;
 gettimeofday(&t,NULL); 
 return (double) t.tv_sec + t.tv_usec/1000000.0;
}

/* dynamically allocate memory for a matrix (n rows, m colums) */
double **dmatrix(int n, int m){
 double **A;
 int i;

 /* allocate memory for pointers to each row */
 A= (double **) calloc(n,sizeof(double*));

 /* allocate a contingent block for the matrix and let A[0] point to it */
 A[0]= (double*) calloc(n*m,sizeof(double));
 if(A[0]==NULL){
   fprintf(stderr,"couldn't reserve memory for an %i x %i matrix\n",n,m);
   exit(EXIT_FAILURE);
 }

 /* let A[i] point to the correct memory position */
 for(i=1;i<n;i++){
    A[i]=A[i-1]+m;
 }

 /* return pointer to the matrix */
 return A;
}

/* dynamically allocate memory for a vector (n rows) */
double *dvector(int n){
 double *v;

 /* allocate a contingent block for the vector and let v point to it */
 v = (double*) calloc(n,sizeof(double));
 if(v==NULL){
   fprintf(stderr,"couldn't reserve memory for an %i vector\n",n);
   exit(EXIT_FAILURE);
 }
 /* return pointer to the matrix */
 return v;
}

/* dynamically allocate memory for a vector (n rows) */
int *ivector(int n){
 int *v;

 /* allocate a contingent block for the vector and let v point to it */
 v = (int*) calloc(n,sizeof(int));
 if(v==NULL){
   fprintf(stderr,"couldn't reserve memory for an %i vector\n",n);
   exit(EXIT_FAILURE);
 }
 /* return pointer to the matrix */
 return v;
}

extern void dgemv_(char *trans, int *m, int *n,
                   double *alpha, double *a, int *lda, double *x, int *incx,
                   double *beta,  double *y, int *incy );

extern void dgesv_(int *n, int *nrhs, double *a, int *lda,
                   int *ipiv, double *b, int *ldb, int *info);

/* (0,1) equally distributed random variable */
inline double randf(void) {
   return (double) (0.5+rand())/(1.0+RAND_MAX);
}

// number of floating point operations
double ops_solve(int n) {
   // LU: sum_1^{n-1} [(n-i) ( 1div, (n-i) mul, (n-i) add )] 
   //     div: sum_1^{n-1} i   = 1/2 ( (n-1)^2 + (n-1) )
   //     add: sum_1^{n-1} i^2 = 1/6 ( 2(n-1)^3 + 3(n-1)^2 + (n-1) )
   //     mul: same as add
   // pivoting: ...
   // L-solve: 
   // U-solve: 
   return 2.0/3.0*n*n*n;         // only correct for n^3
}
double ops_matvec(int n){
   // add: n^2
   // mul: n*(n-1);
   return 2.0*n*n - n;
}

void print_timing(char* prog, double t[5], int n, double norm, int csvout) {
   // time statistics
   if(!csvout){
      printf("\nreal time information:\n");
      printf("memory alloc\t matrix init\t A*x\t\t solving\t total\n");
      printf("%f\t %f\t %f\t %f\t %f\n",t[0],t[1],t[2],t[3],t[4]);
      printf("\nmflops: %f\n\n",ops_solve(n)/t[3]/1e6);
      printf("\n");
   } else {
      printf("%s, %i, %.3e, %.3e, %.4e, %.4e, %.3e\n", prog, n, 
         ops_solve(n)/t[3], ops_matvec(n)/t[2], t[3], t[2], norm);
   }
}

