static char rcsid[] = "$Header$";
//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include <stdio.h>
#include "newmat.h"
#include "precisio.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

//------------------------------------------------------------------------
//
// Perturbed  Cholesky decomposition 
//
//------------------------------------------------------------------------

static Real square(Real x) { return x*x; }
ReturnMatrix PertChol(SymmetricMatrix&, Real, Real&);

ReturnMatrix MCholesky(SymmetricMatrix& S)
{
   Tracer trace("MCholesky");
   int nr = S.Nrows();
   LowerTriangularMatrix L(nr);
   Real mcheps = FloatingPointPrecision::Epsilon();

   Real* s = S.Store(); Real* l = L.Store(); Real* li = l;

   Real maxadd = 0.0;

   int i, j;

   Real sqrteps = sqrt(mcheps);
   Real maxdiag = 0.0;
   Real mindiag = 1.0e10;
   Real maxoff  = 0.0;
   for (i=1; i<=nr; ++i) {
     maxdiag = max(maxdiag,S(i,i));
     mindiag = min(mindiag,S(i,i));
     for (j=i; j<=i; ++j) {
       maxoff = max(maxoff,S(i,j));
     }
   }

   Real maxposdiag = max(0.0,maxdiag);
   Real mu;

   if (mindiag <= sqrteps*maxposdiag) {
     mu = 2.0*(maxposdiag-mindiag)*sqrteps - mindiag;
     maxdiag = maxdiag + mu;
   }
   else mu = 0.0;

   if (maxoff*(1.0 + 2.0*sqrteps) > maxdiag) {
     mu = mu + (maxoff-maxdiag) + 2.0*sqrteps*maxoff;
     maxdiag = maxoff * (1.0 + 2.0*sqrteps);
   }

   if (maxdiag == 0.0) {
     mu = 1.0;
     maxdiag = 1.0;
   }
   if (mu > 0.0) {
     for (i=1; i<=nr; ++i) S(i,i) = S(i,i) + mu;
   }

   Real maxoffl = sqrt(max(maxdiag,(maxoff/nr)));

   L = PertChol(S,maxoffl,maxadd);
   
   
   if (maxadd > 0.0) {
     printf("MCholesky: Hessian not positive definite\n");
     printf("MCholesky: maxadd = %12.4e\n", maxadd);
     Real maxev = S(1,1);
     Real minev = S(1,1);
     for (i=1; i<=nr; ++i) {
       Real offrow = 0.0;
       for(j=1; j<=i-1; ++j) offrow += fabs(S(j,i));
       for(j=i+1; j<=nr; ++j) offrow += fabs(S(i,j));
       maxev = max(maxev,(S(i,i)+offrow));
       minev = min(minev,(S(i,i)-offrow));
      }
     Real sdd = (maxev -minev) * sqrteps - minev;
     sdd = max(sdd,0.0);
     mu = min(maxadd,sdd);
     for (i=1; i<=nr; ++i) S(i,i) = S(i,i) + mu;
     
     L = PertChol(S,0.0,maxadd);
   }
       
   L.Release(); return L;
 }



ReturnMatrix PertChol(SymmetricMatrix& S, Real maxoffl, Real& maxadd)
{
  Tracer trace("PertChol");
  int nr = S.Nrows();
  LowerTriangularMatrix L(nr);
  Real mcheps = FloatingPointPrecision::Epsilon();
  
  Real* s = S.Store(); Real* l = L.Store(); Real* li = l;
  Real sum;
  Real minl2  = 0.0;
  
  int j, k;
  Real minl = pow(mcheps,.25)*maxoffl;
  
  if (maxoffl == 0.0) {
    Real maxdiag = 0.0;
    for (int i=1;i<=nr;++i) maxdiag = max(maxdiag,fabs(S(i,i)));
    maxoffl = sqrt(maxdiag);
    minl2 = sqrt(mcheps)*maxoffl;
  }
  maxadd = 0.0;
  
  for (j=1; j<=nr; j++) {
    sum = 0.0;
    for (int i=1; i<=j-1; ++i) {
      sum += square(L(j,i)); 
    }
    Real ljj = S(j,j) - sum;
    
    Real  minljj = 0.0;
    
    for (i=j+1; i<=nr; ++i) {
      sum = 0.0;
      for(k=1; k<=j-1; ++k) {
	sum += L(i,k) * L(j,k);
      }
      L(i,j) = S(j,i) - sum;
      minljj = max(fabs(L(i,j)),minljj);
    }
    minljj = max((minljj/maxoffl),minl);
    
    if (ljj > square(minljj)) { // Normal Cholesky
      L(j,j) = sqrt(ljj);
    }
    else {//    Modify ljj since it is too small
      if (minljj < minl2) minljj = minl2;
      maxadd = max(maxadd,(square(minljj)-ljj));
      L(j,j) = minljj;
    }
    for (i=j+1; i<=nr; ++i){
      L(i,j) = L(i,j) / L(j,j);
    }
    
  }

  L.Release();
  return L;
}


