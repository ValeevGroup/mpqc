//$$ hholder.cxx                       Householder triangularisation

// Copyright (C) 1991,2,3: R B Davies

#define WANT_MATH

#include "include.h"

#include "newmatap.h"


/********************** householder triangularisation *********************/

inline Real square(Real x) { return x*x; }

void HHDecompose(Matrix& X, LowerTriangularMatrix& L)
{
   Tracer et("HHDecompose(1)");
   int n = X.Ncols(); int s = X.Nrows(); L.ReDimension(s);
   Real* xi = X.Store(); int k;
   for (int i=0; i<s; i++)
   {
      Real sum = 0.0;
      Real* xi0=xi; k=n; while(k--) { sum += square(*xi++); }
      sum = sqrt(sum);
      L.element(i,i) = sum;
      if (sum==0.0) Throw(SingularException(L));
      Real* xj0=xi0; k=n; while(k--) { *xj0++ /= sum; }
      for (int j=i+1; j<s; j++)
      {
         sum=0.0;
	 xi=xi0; Real* xj=xj0; k=n; while(k--) { sum += *xi++ * *xj++; }
	 xi=xi0; k=n; while(k--) { *xj0++ -= sum * *xi++; }
	 L.element(j,i) = sum;
      }
   }
}

void HHDecompose(const Matrix& X, Matrix& Y, Matrix& M)
{
   Tracer et("HHDecompose(1)");
   int n = X.Ncols(); int s = X.Nrows(); int t = Y.Nrows();
   if (Y.Ncols() != n)
   { Throw(ProgramException("Unequal row lengths",X,Y)); }
   M.ReDimension(t,s);
   Real* xi = X.Store(); int k;
   for (int i=0; i<s; i++)
   {
      Real* xj0 = Y.Store(); Real* xi0 = xi;
      for (int j=0; j<t; j++)
      {
         Real sum=0.0;
         xi=xi0; Real* xj=xj0; k=n; while(k--) { sum += *xi++ * *xj++; }
	 xi=xi0; k=n; while(k--) { *xj0++ -= sum * *xi++; }
	 M.element(j,i) = sum;
      }
   }
}

