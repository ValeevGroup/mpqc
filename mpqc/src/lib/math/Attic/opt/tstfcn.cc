static char rcsid[] = "$Header$";

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#define WANT_MATH
#define WANT_STREAM
#include "nlp.h"
static Real square(Real x) { return x*x; }

void quad(int mode, int n, ColumnVector x, double& fx, ColumnVector& g, 
	  SymmetricMatrix& H)
/* Simple quadratic to test out optimization routines */
/* f = x1**2 + x2**2 + ...  + xn**2 */

{
  fx = 0;
  for (int i=1; i<=n; ++i) {
    fx += x(i)*x(i);
    g(i) = 2.0*x(i);
  }
  if (mode == 2) {
    H = 0.0;
    for(i=1; i<=n; ++i) H(i,i) = 2.0;
  } 
}

void ds_ex9(int mode, int n, ColumnVector x, double& fx, ColumnVector& g, 
	     SymmetricMatrix& H)
// Dennis & Schnabel Ex 71. p.164, n= 2

{
  double x1, x2, f1, f2, f1sq;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = x1 - 2.0;
  f2 = x2 + 1.0;
  f1sq = f1*f1;
  fx = f1sq*f1sq + f1sq*x2*x2 + f2*f2;

  g(1) = 4.0*f1sq*f1 + 2.0*f1*x2*x2;
  g(2) = 2.0*f1sq*x2 + 2.0*f2;

}


void ds_ex71(int mode, int n, ColumnVector x, double& fx, ColumnVector& g, 
	     SymmetricMatrix& H)
// Dennis & Schnabel Ex 71. p.164, n= 2

{
  double x1, x2, x1sq, x2sq;

  if (n != 2) return;

  x1 = x(1)-1.e6;
  x2 = x(2)-1.e-6;
  x1sq = x1*x1;
  x2sq = x2*x2;

  fx = x1sq*x1sq + x1sq + x2sq*x2sq + x2sq;
  g(1) = 4.0*x1sq*x1 + 2.0*x1;
  g(2) = 4.0*x2sq*x2 + 2.0*x2;

}


void rosen(int mode, int n, ColumnVector x, double& fx, ColumnVector& g, 
	   SymmetricMatrix& H)
// Rosenbrock's function, n= 2

{ 
  double f1, f2, x1, x2;
  
  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = (x2 - x1 * x1);
  f2 = 1. - x1;
  
  if (mode == 0 || mode == 3) {
    fx  = 100.* f1*f1 + f2*f2;
  }
  if (mode == 1 || mode == 3) {
    g(1) = -400.*f1*x1 - 2.*f2; 
    g(2) = 200.*f1;
  }
  
  if (mode == 2 || mode == 3) {
    f1 = (x2 - 3.0*x1*x1);
    
    H(1,1) = -400.0*f1 + 2.0;
    H(2,1) = -400.0*x1;
    H(2,2) = 200.0;
  }
}
void erosen(int mode, int n, ColumnVector x, double& fx, ColumnVector& g, 
	    SymmetricMatrix& H)
{
//  erosen
//  Extended rosenbrock function for any even value of n. 
//                          2 
//     f    (x) = 10(x   - x    ) 
//      2i-1          2i    2i-1 
//
//     f  (x) = 1 - x 
//      2i           2i-1 
//                                           2         2 
//     f(x) = sum from i = 1 to n/2 (f    (x)  + f  (x) ). 
//                                    2i-1        2i 
//
//     x  = (epsilon ) where epsilon    = -1.2, epsilon  = 1 
//      0           i               2i-1               2i 
//
//     f(x ) = 0 at x  = (1,...,1) 
//        *          * 
//
//     parameters      
//
//     input 
//        n             dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  static int i;
  static double f1, f2, x1, x2;
  

  fx = 0.;

  if (n % 2 != 0) {
    cerr << "erosen: n must be even, n= \n" << n;
    return;
  }

  for (i = 1; i <= n/2; ++i) {
    x1 = x(2*i-1);
    x2 = x(2*i);
    
    f1 = (x2 - x1 * x1);
    f2 = 1. - x1;
    
    if (mode == 0 || mode == 3){
      fx += 100.*square(f1) + square(f2);
    }
    if (mode == 1 || mode == 3) {
      g(2*i-1) = -400.*f1*x1 - 2.*f2; 
      g(2*i)   = 200.*f1;
    }
  }
  
}
void epowell(int mode, int n, ColumnVector x, double& fx, ColumnVector& g, 
	  SymmetricMatrix& H)
// Not implemented yet
{
  fx = 0;
}

void trig(int mode, int n, ColumnVector x, double& fx, ColumnVector& g, 
	  SymmetricMatrix& H)

{
  fx = 0;
  if (mode == 0 || mode == 3) {
    for (int i=1; i<=n; ++i) {
      double sum = 0.0;
      for (int j=1; j<=n; ++j) {
	sum += cos(x(j)) + double(i)*(1.0-cos(x(i)) - sin(x(i)));
      }
      fx -= sum;
    }  
    fx += n;
  }

  if (mode == 2) {
    H = 0.0;
    for(int i=1; i<=n; ++i) H(i,i) = 2.0;
  } 
}

void helical(int mode, int n, ColumnVector x, double& fx, ColumnVector& g, 
	  SymmetricMatrix& H)

{
  double x1, x2, x3;
  double f1, f2, f3;

  fx = 0.0;
  if (n != 3) return;

  x1 = x(1); x2 = x(2); x3 = x(3);
  f1 = 10.0*(x3 - 10.0*atan2(x1,x2)/(2.0*M_PI));
  f2 = 10.0*(sqrt(square(x1)+square(x2)) - 1.0);
  f3 = x3;
  if (mode == 0 || mode == 3) {
    fx = square(f1) + square(f2) + square(f3);
  }

  if (mode == 1 || mode == 3) {
    for (int i=1; i<=n; ++i) {
      g(i) = 2.0*x(i);
    }
  }
  if (mode == 2 || mode == 3) {
    H = 0.0;
    for(int i=1; i<=n; ++i) {H(i,i) = 1.0;}
  } 
}

void wood(int mode, int n, ColumnVector x, double& fx, ColumnVector& g, 
	  SymmetricMatrix& H)
{
  double x1, x2, x3, x4;
  double f1, f2, f3, f4, f5, f6;

  if (n != 4) return;

  x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);

  f1 = 1.0 - x1;
  f2 = 1.0 - x2;
  f3 = 1.0 - x3;
  f4 = 1.0 - x4;
  f5 = x1*x1 - x2;
  f6 = x3*x3 - x4;

  if (mode == 0 || mode == 3) {
    fx = 100.0*square(f5) + square(f1) + 90.0*square(f6) + square(f3) + 
    10.1*(square(f2)+square(f4)) + 19.8*f2*f4;
  }

  if (mode == 1 || mode == 3) {
    g(1) =  400.*f5*x1 - 2.0*f1;
    g(2) = -200.*f5 - 20.2*f2 - 19.8*f4;
    g(3) =  360.*f6*x3 - 2.0*f3;
    g(4) = -180.*f6 - 20.2*f4 - 19.8*f2 ;
  }
}
