//$$ cholesky.cxx                     cholesky decomposition

// Copyright (C) 1991,2,3: R B Davies

#define WANT_MATH

#include "include.h"

#include "newmat.h"


/********* Cholesky decomposition of a positive definite matrix *************/

// Suppose S is symmetrix and positive definite. Then there exists a unique
// lower triangular matrix L such that L L.t() = S;

static Real square(Real x) { return x*x; }

ReturnMatrix Cholesky(const SymmetricMatrix& S)
{
   Tracer trace("Cholesky");
   int nr = S.Nrows();
   LowerTriangularMatrix T(nr);
   Real* s = S.Store(); Real* t = T.Store(); Real* ti = t;
   for (int i=0; i<nr; i++)
   {
      Real* tj = t; Real sum; int k;
      for (int j=0; j<i; j++)
      {
	 Real* tk = ti; sum = 0.0; k = j;
	 while (k--) { sum += *tj++ * *tk++; }
	 *tk = (*s++ - sum) / *tj++;
      }
      sum = 0.0; k = i;
      while (k--) { sum += square(*ti++); }
      Real d = *s++ - sum;
      if (d<=0.0) Throw(NPDException(S));
      *ti++ = sqrt(d);
   }

// Mike Colvin made this change to get example to work with gcc 2.3.1
//   T.Release(); return (ReturnMatrix)T;
   T.Release(); return T;
}

ReturnMatrix Cholesky(const SymmetricBandMatrix& S)
{
   Tracer trace("Band-Cholesky");
   int nr = S.Nrows(); int m = S.lower;
   LowerBandMatrix T(nr,m);
   Real* s = S.Store(); Real* t = T.Store(); Real* ti = t;

   for (int i=0; i<nr; i++)
   {
      Real* tj = t; Real sum; int l;
      if (i<m) { l = m-i; s += l; ti += l; l = i; }
      else { t += (m+1); l = m; }

      for (int j=0; j<l; j++)
      {
	 Real* tk = ti; sum = 0.0; int k = j; tj += (m-j);
	 while (k--) { sum += *tj++ * *tk++; }
	 *tk = (*s++ - sum) / *tj++;
      }
      sum = 0.0;
      while (l--) { sum += square(*ti++); }
      Real d = *s++ - sum;
      if (d<=0.0)  Throw(NPDException(S));
      *ti++ = sqrt(d);
   }

// Mike Colvin made this change to get example to work with gcc 2.3.1
//   T.Release(); return (ReturnMatrix)T;
   T.Release(); return T;
}

