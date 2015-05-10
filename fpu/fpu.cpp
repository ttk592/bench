/*
 * fpu - a simple performance test for several floating point operations
 *
 * execute: ./fpu <MHz of cpu> <mill ops>, e.g. ./fpu 2667 100
 *
 * note: inspect the generated assembly code fpu.s to see whether the
 *    c++ compiler has optimised the code appropriately
 *
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
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cassert>
#include "stoptime.h"

#if __GNUC__ <= 2
   #define COMPAT_GCC2
#endif

#ifndef COMPAT_GCC2
   #include <fenv.h>
   #define ASM_COMMENT(X)  asm volatile ("#" X)
#else
   #define ASM_COMMENT(X)  asm ("# userdefined label")
#endif


#define PRINT_FLOPS(KFLOPS) \
{                                                              \
   if(kflops<1e2){                                             \
      printf("%7.3f kflops, ", kflops);                        \
   } else if(kflops<1e5) {                                     \
      printf("%7.3f Mflops, ", kflops/1e3);                    \
   } else {                                                    \
      printf("%7.3f Gflops, ", kflops/1e6);                    \
   }                                                           \
   printf("%.3f flops/cycle\n",kflops/1e3/mhz);                \
}


#define RUN_TEST(ROUTINE, VAR, OPS) \
{                                                              \
   double t=stoptime();                                        \
   double res=ROUTINE(VAR,OPS);                                \
   t=stoptime()-t;                                             \
   printf("%s:\t%.3fs, res=%10.3e,\t",                         \
         # ROUTINE,t,res);                                     \
   double kflops=(double)OPS/t/1e3;                            \
   PRINT_FLOPS(kflops);                                        \
}

#define RUN_TEST2(ROUTINE, VAR1, VAR2, OPS) \
{                                                              \
   double t=stoptime();                                        \
   double res=ROUTINE(VAR1, VAR2, OPS);                        \
   t=stoptime()-t;                                             \
   printf("%s:\t%.3fs, res=%10.3e,\t",                         \
         # ROUTINE,t,res);                                     \
   double kflops=(double)OPS/t/1e3;                            \
   PRINT_FLOPS(kflops);                                        \
}


double rand01() {
   return random()/(RAND_MAX+1.0);
}

// direct assembly input (circumvents poor optimisation of old gcc compilers)
#ifndef ASM_INTEL
double addx87(double a, size_t ops) {
   double   x=0.0;

   // this only compiles if at&t syntax is enabled
   // gcc compiler options -masm=intel and -masm=att do not set a ifdef flag
   // if masm=intel is activated we could bracket the asm code with
   //    .att_syntax prefix
   //    ...
   //    .intel_syntax no_prefix
   // however, this doesn't work, as asm() will replace %0, %1, etc with
   // intel syntax

   asm ( "# start addx87\n\t"
         "movl %2, %%ebx\n\t"        // ebx = ops
         "xorl %%eax, %%eax\n\t"     // eax = 0
         "fldl %1\n\t"              // st(1) = a
         "fldz\n"                   // st(0) = 0.0
         "1:\n\t"                   // local label
         "fadd    %%st(1)\n\t"
         "fadd    %%st(1)\n\t"
         "fadd    %%st(1)\n\t"
         "fadd    %%st(1)\n\t"
         "fadd    %%st(1)\n\t"
         "fadd    %%st(1)\n\t"
         "fadd    %%st(1)\n\t"
         "fadd    %%st(1)\n\t"
         "addl    $8, %%eax\n\t"
         "cmpl    %%eax, %%ebx\n\t" 
         "jg      1b\n\t"          // goto local label 1 (backwards)
         "fstpl %0\n\t"              // x=st(0), pop
         "fstp %%st(0)\n\t"          // pop (store to itself and pop)
         "# end addx87\n"
         : "=m" (x) : "m" (a), "m" (ops) : "%eax", "%ebx");

   return x-a*ops;
}
#endif // ASM_INTEL



// calculates x +  a + a + ... + a (for many (Num) randomly intialised x's)
// compiler won't simplify to x+n*a due to rounding issues
// (unless compiler flags -ffast-math / -funsafe-math-optimizations are used)
template <size_t Num> 
double add(double a, size_t ops) {
   assert(ops>0);
   double   X[Num];
   double   expected=0.0;
   for(size_t i=0; i<Num; i++){ X[i]=rand01(); expected+=X[i]; }

   // main add loop
   const size_t loops = ops/Num;
   ASM_COMMENT("start add");
   for(size_t k=0; k<loops; k++) {
      // compiler must fully unroll below loop, or code won't be efficient
      for(size_t i=0; i<Num; i++){ X[i]+=a; }
   }
   ASM_COMMENT("end add");

   expected += a*Num*loops;
   double actual=0.0;
   for(size_t i=0; i<Num; i++){ actual+=X[i]; }
   return actual-expected;
}
// specialisation of add<1>, just to check compiler optimisation
double add1(double a, size_t ops){
   assert(ops>0);
   double   x0=rand01();
   double   expected=x0;
   const size_t loops = ops;

   ASM_COMMENT("start add1");
   for(size_t k=0; k<loops; k++) {
      x0+=a;
   }
   ASM_COMMENT("end add1");

   expected += a*loops;
   double actual=x0;
   return actual-expected;
}
// specialisation of add<2>, just to check compiler optimisation
double add2(double a, size_t ops){
   assert(ops>0);
   double   x0=rand01();
   double   x1=rand01();
   double   expected=x0+x1;
   const size_t loops = ops/2;

   ASM_COMMENT("start add2");
   for(size_t k=0; k<loops; k++) {
      x0+=a; x1+=a;
   }
   ASM_COMMENT("end add2");

   expected += 2*a*loops;
   double actual=x0+x1;
   return actual-expected;
}


// calculates x * a * a * ... * a (for many (Num) randomly intialised x's)
// compiler won't simplify
template <size_t Num> 
double mul(double a, size_t ops) {
   assert(ops>0);
   double   X[Num];
   double   expected=0.0;
   for(size_t i=0; i<Num; i++){ X[i]=rand01(); expected+=X[i]; }

   // main add loop
   const size_t loops = ops/Num;
   ASM_COMMENT("start mul");
   for(size_t k=0; k<loops; k++) {
      for(size_t i=0; i<Num; i++){ X[i]*=a; }
   }
   ASM_COMMENT("end mul");

   expected *= pow(a, (double)loops);
   double actual=0.0;
   for(size_t i=0; i<Num; i++){ actual+=X[i]; }
   return actual-expected;
}
// specialisation of mul<1>, just to check compiler optimisation
double mul1(double a, size_t ops) {
   assert(ops>0);
   double   x0=rand01();
   double   expected=x0;
   const size_t loops = ops;

   ASM_COMMENT("start mul1");
   for(size_t k=0; k<loops; k++) {
      x0*=a;
   }
   ASM_COMMENT("end mul1");

   expected *= pow(a, (double)loops);
   double actual=x0;
   return actual-expected;
}


// calculates x / a / a / ... / a (for many (Num) randomly inititialised x's)
// compiler doesn't simplify to x * (1/a) * (1/a) * ... due to rounding issues?
template <size_t Num> 
double div(double a, size_t ops) {
   assert(ops>0);
   double   X[Num];
   double   expected=0.0;
   for(size_t i=0; i<Num; i++){ X[i]=rand01(); expected+=X[i]; }

   // main add loop
   const size_t loops = ops/Num;
   ASM_COMMENT("start div");
   for(size_t k=0; k<loops; k++) {
      // compiler most likely fully unrolls this loop
      for(size_t i=0; i<Num; i++){ X[i]/=a; }
   }
   ASM_COMMENT("end div");

   expected *= pow(a, -((double)loops));
   double actual=0.0;
   for(size_t i=0; i<Num; i++){ actual+=X[i]; }
   return actual-expected;
}
// specialisation of div<1>, just to check compiler optimisation
double div1(double a, size_t ops) {
   assert(ops>0);
   double   x0=rand01();
   double   expected=x0;
   const size_t loops = ops;

   ASM_COMMENT("start div1");
   for(size_t k=0; k<loops; k++) {
      x0/=a;
   }
   ASM_COMMENT("end div1");

   expected *= pow(a, -((double)loops));
   double actual=x0;
   return actual-expected;
}

// calculates sqrt(sqrt(sqrt(...(a)...))), 100 (innerloop) times,
// outer loop periodically rescales result so that we don't converge to 1.0
double sqroot(double a, size_t ops) {
   assert(ops>0);
   double x=1.0;
   const size_t innerloop=100;
   const size_t loops=ops/innerloop;

   // main add loop
   for(size_t k=0; k<loops; k++){
      x*=a;
      ASM_COMMENT("start sqrt");
      for(size_t i=0; i<innerloop; i++) {
         x=sqrt(x);
      }
      ASM_COMMENT("end sqrt");
   }
   return x-1.0;
}

// calculates (...((a^b)^b)^ ... ^b), 100 (innerloop) times,
// outer loop periodically rescales result so that we don't converge or explode
double power(double a, double b, size_t ops) {
   assert(ops>0);
   double x=1.0;
   const size_t innerloop=100;
   const size_t loops=ops/innerloop;

   // main add loop
   for(size_t k=0; k<loops; k++){
      x*=a;
      ASM_COMMENT("start pow");
      for(size_t i=0; i<innerloop; i++) {
         x=pow(x, b);
      }
      ASM_COMMENT("end pow");
   }
   return x-1.0;
}

// calculates iteration x_{i+1} = x_i+exp(-x_i), x_0=a
// additions are not counted as one operation (as much faster than exp)
double expon(double a, size_t ops) {
   assert(ops>0);
   double x=a;

   ASM_COMMENT("start exp");
   for(size_t k=0; k<ops; k++){
      x+=exp(-x);
   }
   ASM_COMMENT("end exp");
   return x;
}

double logs(double a, size_t ops) {
   assert(ops>0);
   double x=a;

   ASM_COMMENT("start log");
   for(size_t k=0; k<ops; k++){
      x+=log(x);
   }
   ASM_COMMENT("end log");
   return x;
}

double logexp(double a, size_t ops) {
   assert(ops>0);
   const size_t loops=ops/2;
   double x=a;

   ASM_COMMENT("start logexp");
   for(size_t k=0; k<loops; k++){
      x=exp(x)+1.1e-5;
      x=log(x);
   }
   ASM_COMMENT("end logexp");
   return x/a;
}

double sinus(double a, size_t ops) {
   assert(ops>0);
   double x=a;

   ASM_COMMENT("start sin");
   for(size_t k=0; k<ops; k++){
      x=sin(x)+0.31;
   }
   ASM_COMMENT("end sin");
   return x;
}

double erfunc(double a, size_t ops) {
   assert(ops>0);
   double x=a;

   ASM_COMMENT("start erf");
   for(size_t k=0; k<ops; k++){
      x=erf(x)+0.31;
   }
   ASM_COMMENT("end erf");
   return x;
}

template <size_t Num> 
double addmul(double scale, size_t ops){
   const size_t loops=ops/Num/2;
   double Fac[Num];
   double Prod[Num];
   double Sum[Num];

   for(size_t i=0; i<Num; i++){
      Fac[i]  = 1.0+scale*rand01();
      Prod[i] = 1.0;
      Sum[i]  = 0.0;
   }

   ASM_COMMENT("start addmul");
   for(size_t k=0; k<loops; k++) {
      for(size_t i=0; i<Num; i++){
         Prod[i] *= Fac[i];
         Sum[i]  += Prod[i];
      }
   }
   ASM_COMMENT("end addmul");

   double actual=0.0;
   for(size_t i=0; i<Num; i++){ actual +=   Sum[i]; }
   return actual;
}


// calculates (x + a + ... + a)  + (y * m * ... * m) - (x*n*a) - (y*m^n)
template <size_t Num> 
double addmul2(double add, double mul, size_t ops){
   //add=mul;
   const size_t loops=ops/Num/2;
   double X[Num];
   double Y[Num];

   double sum1=0, sum2=0;
   for(size_t i=0; i<Num; i++){
      X[i]=rand01(); Y[i]=rand01();
      sum1+=X[i]; sum2+=Y[i];
   }
   double expected = sum1+add*loops*Num + sum2*pow(mul,loops);

   ASM_COMMENT("start addmul2");
   for(size_t k=0; k<loops; k++) {
      for(size_t i=0; i<Num; i++){
         X[i]+=add;
         Y[i]*=mul;
      }
   }
   ASM_COMMENT("end addmul2");

   double actual=0.0;
   for(size_t i=0; i<Num; i++){ actual+=(X[i]+Y[i]); }
   return actual-expected;
}


// if sse2 instructions are enabled (e.g. -msse2)
#ifdef __SSE2__
#include <emmintrin.h>  // SSE2
#define ALIGN( var ) var __attribute__ ( (__aligned__ ( 32 ) ) )
#define _mm_extract_pd(R, I) (*((double*)(&R)+I))     // dodgy I think
// use _MM_EXTRACT_FLOAT for floats
// http://stackoverflow.com/questions/3130169/
// MSVC compiler gives access via x.m128d_f64[0] to individual elements:
// x.m128_f32[0]
/*
union my_m128d {
    __m128d v;
    double  e[2];
};
typedef union __declspec(align(16)) __m128d { 
    double m128d_f64[2];
} __m128d;
*/

// calculates sum_{i=0}^n x^i, for many (2*Num) different x's
// x's are randomly initialised as 1+scale*rand01()
// this calculation is very similar to scalar vec*vec product as
// we basically calculate:
//    prod*=fac; sum+=prod;
template <size_t Num> 
double addmulsse(double scale, size_t ops){
   const size_t loops=ops/Num/4;
   __m128d Fac[Num];
   __m128d Prod[Num];
   __m128d Sum[Num];

   for(size_t i=0; i<Num; i++){
      Fac[i]  = _mm_set_pd(1.0+scale*rand01(), 1.0+scale*rand01());
      Prod[i] = _mm_set_pd(1.0, 1.0);
      Sum[i]  = _mm_set_pd(0.0, 0.0);
   }

   ASM_COMMENT("start addmulsse");
   for(size_t k=0; k<loops; k++) {
      for(size_t i=0; i<Num; i++){
         Prod[i] *= Fac[i];
         Sum[i]  += Prod[i];
      }
   }
   ASM_COMMENT("end addmulsse");

   double actual=0.0;
   for(size_t i=0; i<Num; i++){
      actual +=   _mm_extract_pd(Sum[i],0) + _mm_extract_pd(Sum[i],1);
   }
   return actual;
}


template <size_t Num> 
double addmul2sse(__m128d add, __m128d mul, size_t ops){
   add=mul;
   const size_t loops=ops/Num/4;
   __m128d X[Num];
   __m128d Y[Num];

   double   sum1=0, sum2=0;
   for(size_t i=0; i<Num; i++){
      X[i]=_mm_set_pd(rand01(), rand01());
      Y[i]=_mm_set_pd(rand01(), rand01());
      sum1 += _mm_extract_pd(X[i],0) + _mm_extract_pd(X[i],1);
      sum2 += _mm_extract_pd(Y[i],0) + _mm_extract_pd(Y[i],1);
   }
   double expected = sum1 +
         (_mm_extract_pd(add,0)+_mm_extract_pd(add,0))*loops*Num +
          sum2*(  pow(_mm_extract_pd(mul,0),loops) *
                  pow(_mm_extract_pd(mul,1),loops) );

   ASM_COMMENT("start addmul2sse");
   for(size_t k=0; k<loops; k++) {
      for(size_t i=0; i<Num; i++){
         X[i] += add;
         Y[i] *= mul;
//         X[i] = _mm_add_pd(X[i], add);
//         Y[i] = _mm_mul_pd(Y[i], mul);
      }
   }
   ASM_COMMENT("end addmul2sse");

   double actual=0.0;
   for(size_t i=0; i<Num; i++){
      actual +=   _mm_extract_pd(X[i],0) + _mm_extract_pd(X[i],1) +
                  _mm_extract_pd(Y[i],0) + _mm_extract_pd(Y[i],1) ;
   }
   return actual-expected;
}


#endif   // SSE2


// if avx instructions are enabled (e.g. -mavx)
#ifdef __AVX__
#include <immintrin.h>     // AVX intrinsics

template <size_t Num> 
double addmulavx(double scale, size_t ops){
   const size_t loops=ops/Num/8;
   __m256d Fac[Num];
   __m256d Prod[Num];
   __m256d Sum[Num];

   for(size_t i=0; i<Num; i++){
      Fac[i]  = _mm256_set_pd(1.0+scale*rand01(), 1.0+scale*rand01(),
                                1.0+scale*rand01(), 1.0+scale*rand01());
      Prod[i] = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
      Sum[i]  = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
   }

   ASM_COMMENT("start addmulavx");
   for(size_t k=0; k<loops; k++) {
      for(size_t i=0; i<Num; i++){
         Prod[i] *= Fac[i];
         Sum[i]  += Prod[i];
      }
   }
   ASM_COMMENT("end addmulavx");

   double actual=0.0;
   for(size_t i=0; i<Num; i++){
      double ALIGN( sumd[4] );
      _mm256_store_pd(sumd, Sum[i]);
      actual +=   sumd[0] + sumd[1] + sumd[2] + sumd[3];
   }
   return actual;
}


#endif   // AVX

// define this global because of simple print macro
double   mhz;

// ======================================================================
//         Name:  main
//  Description:  main function, executed at runtime
// ======================================================================
int main(int argc, char** argv) {
#ifndef COMPAT_GCC2
   // we don't want performance degradation due to nan's and inf's
   // so if anything occurs we exit
   feenableexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
#endif

   if(argc!=3) {
      printf("usage: %s <cpu MHz> <mill ops>\n", argv[0]);
      exit(EXIT_FAILURE);
   }
   mhz=atof(argv[1]);
   int      n=atoi(argv[2])*1000000;
   if(n<=0) n=10000;

   double a,b;


   printf("Note, if the compiler has properly optimised this code then\n");
   printf("addx87, add1 and add<1> should perform equally fast. If not,\n");
   printf("gcc doesn't efficiently optimise and this test is unreliable.\n\n");

   // addition
   a=1.5217;
   printf("calculating: x<k> + a + a + ... + a - (x<k>+n*a), with a=%f\n",a);
#ifndef ASM_INTEL
   RUN_TEST(addx87,a,n);
#endif

   RUN_TEST(add1,a,n);
   RUN_TEST(add2,a,n);

   RUN_TEST(add<1>,a,n);
   RUN_TEST(add<2>,a,n);
   RUN_TEST(add<3>,a,n);
   RUN_TEST(add<4>,a,n);
   RUN_TEST(add<5>,a,n);
   RUN_TEST(add<6>,a,n);
   RUN_TEST(add<7>,a,n);
   RUN_TEST(add<8>,a,n);
   RUN_TEST(add<9>,a,n);
   printf("\n");

   // multiplication
   a=1.0+1.7e-8;
   printf("calculating: x<k> * a * a * ... * a - (x<k>*a^n), with a=%f\n",a);
   RUN_TEST(mul1,a,n);
   RUN_TEST(mul<1>,a,n);
   RUN_TEST(mul<2>,a,n);
   RUN_TEST(mul<3>,a,n);
   RUN_TEST(mul<4>,a,n);
   RUN_TEST(mul<5>,a,n);
   RUN_TEST(mul<6>,a,n);
   RUN_TEST(mul<7>,a,n);
   RUN_TEST(mul<8>,a,n);
   RUN_TEST(mul<9>,a,n);
   printf("\n");

   // division
   a=1.0-1.7e-8;
   printf("calculating: x<k> / a / a / ... / a - (x<k>/a^n), with a=%f\n",a);
   RUN_TEST(div1,a,n/10);
   RUN_TEST(div<1>,a,n/10);
   RUN_TEST(div<2>,a,n/10);
   RUN_TEST(div<3>,a,n/10);
   printf("\n");

   // sqrt
   a=1.53e156;
   RUN_TEST(sqroot,a,n/10);
   // exp
   a=1.1;
   RUN_TEST(expon,a,n/40);
   // log
   a=1.1;
   RUN_TEST(logs,a,n/40);
   // combined log/exp
   a=1.1;
   RUN_TEST(logexp,a,n/40);
   // sinus
   a=0.1;
   RUN_TEST(sinus,a,n/40);
   // erf
   a=1.1;
   RUN_TEST(erfunc,a,n/40);
   // pow
   a=1.53e156; b=0.86;
   RUN_TEST2(power,a,b,n/100);

   printf("\n");




   // combined addition and multiplication
   a=1.31e-8;
   printf("calculating sum_i x^i, with <k> different x's randomly ");
   printf("of order 1+%e)\n",a);
   printf("basically:  prod*=x; sum+=prod;\n");
   RUN_TEST(addmul<1>,a,n);
   RUN_TEST(addmul<2>,a,n);
   RUN_TEST(addmul<3>,a,n);
   RUN_TEST(addmul<4>,a,n);
   RUN_TEST(addmul<5>,a,n);
   RUN_TEST(addmul<6>,a,n);
   RUN_TEST(addmul<7>,a,n);
   RUN_TEST(addmul<8>,a,n);
   RUN_TEST(addmul<9>,a,n);

   a=4.12;
   b=1.0+1.31e-8;
   printf("calculating x<k>+a+...+a + y<k>*b*...*b - exact\n");  
   printf("basically:  x+=a; y*=b;\n");
   RUN_TEST2(addmul2<1>,a,b,n);
   RUN_TEST2(addmul2<2>,a,b,n);
   RUN_TEST2(addmul2<3>,a,b,n);
   RUN_TEST2(addmul2<4>,a,b,n);
   RUN_TEST2(addmul2<5>,a,b,n);
   RUN_TEST2(addmul2<6>,a,b,n);
   RUN_TEST2(addmul2<7>,a,b,n);
   RUN_TEST2(addmul2<8>,a,b,n);
   RUN_TEST2(addmul2<9>,a,b,n);
   printf("\n");

#ifdef __SSE2__
   // combined addition and multiplication sse
   a=1e-10;
   printf("sse version of addmul\n");
   RUN_TEST(addmulsse<1>,a,n);
   RUN_TEST(addmulsse<2>,a,n);
   RUN_TEST(addmulsse<3>,a,n);
   RUN_TEST(addmulsse<4>,a,n);
   RUN_TEST(addmulsse<5>,a,n);
   RUN_TEST(addmulsse<6>,a,n);
   RUN_TEST(addmulsse<7>,a,n);
   RUN_TEST(addmulsse<8>,a,n);
   RUN_TEST(addmulsse<9>,a,n);

   __m128d ma=_mm_set_pd(1.1,3.31);
   __m128d mb=_mm_set_pd(1.0+1.77e-10,1.0+3.59e-9);
   printf("sse version of addmul2, but with a=b hard coded\n");
   printf("basically:  x+=a; y*=a;\n");
   RUN_TEST2(addmul2sse<1>,ma,mb,n);
   RUN_TEST2(addmul2sse<2>,ma,mb,n);
   RUN_TEST2(addmul2sse<3>,ma,mb,n);
   RUN_TEST2(addmul2sse<4>,ma,mb,n);
   RUN_TEST2(addmul2sse<5>,ma,mb,n);
   RUN_TEST2(addmul2sse<6>,ma,mb,n);
   RUN_TEST2(addmul2sse<7>,ma,mb,n);
   RUN_TEST2(addmul2sse<8>,ma,mb,n);
   RUN_TEST2(addmul2sse<9>,ma,mb,n);

#endif   // SSE2

#ifdef __AVX__
   // combined addition and multiplication sse
   a=1e-10;
   printf("avx version of addmul\n");
   RUN_TEST(addmulavx<1>,a,n);
   RUN_TEST(addmulavx<2>,a,n);
   RUN_TEST(addmulavx<3>,a,n);
   RUN_TEST(addmulavx<4>,a,n);
   RUN_TEST(addmulavx<5>,a,n);
   RUN_TEST(addmulavx<6>,a,n);
   RUN_TEST(addmulavx<7>,a,n);
   RUN_TEST(addmulavx<8>,a,n);
   RUN_TEST(addmulavx<9>,a,n);
#endif // AVX

   return EXIT_SUCCESS;
}

