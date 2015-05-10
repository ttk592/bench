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

#include <cstdio>
#include <cstdlib>
#include <Eigen/Core>
#include <Eigen/LU>
#include "helper.h"


// ======================================================================
//         Name:  main
//  Description:  main function, executed at runtime
// ======================================================================
int main(int argc, char** argv) {
   if(argc<2) {
      printf("usage: %s <order of matrix> [csvout]\n", argv[0]);
      exit(EXIT_FAILURE);
   }
   int   n=atoi(argv[1]);
   if(n<=0) n=1;
   bool  csvout=false;
   if(argc!=2) csvout=true;


   double t[5];          // timer: 0 alloc, 1 init, 2 A*x, 3 solve, 4 total

   t[4]=stoptime();

   // reserve memory
   t[0]=stoptime();
   Eigen::MatrixXd   A(n,n);
   Eigen::VectorXd   b(n);
   Eigen::VectorXd   x(n);
   Eigen::VectorXd   x_sol(n);
   Eigen::VectorXd   x_err(n);
   t[0]=stoptime()-t[0];

   // init random matrix
   t[1]=stoptime();
   A.setRandom();
   x_sol.setRandom();
   t[1]=stoptime()-t[1];

   // b = Ax
   t[2]=stoptime();
   b=A*x_sol;
   t[2]=stoptime()-t[2];

   // solve the linear system
   if(!csvout) { printf("solving ... "); fflush(stdout); }
   t[3]=stoptime();
   //x=A.fullPivLu().solve(b);
   x=A.lu().solve(b);
   t[3]=stoptime()-t[3];
   if(!csvout) { printf("finished\n"); fflush(stdout); }

   // check the numerical error
   x_err=x-x_sol;
   if(!csvout) { printf("The norm of the error is %e\n",x_err.norm()); }
   t[4]=stoptime()-t[4];

   // time statistics
   print_timing(argv[0],t,n,x_err.norm(),csvout);

   return EXIT_SUCCESS;
}


