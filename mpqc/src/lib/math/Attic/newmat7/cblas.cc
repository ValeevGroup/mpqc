static char rcsid[] = "$Header$";

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#define WANT_STREAM
#include "newmat.h"
#include "cblas.h"

#define REPORT {}
//#define REPORT { static ExeCounter ExeCount(__LINE__,8); ExeCount++; }


Real ColumnVector::Dot(ColumnVector& y)
{
  REPORT                                         // not accessed
  int n1 = storage; int n2 = y.storage;
  if (n1 != n2) {
    cerr << "ColumnVector::Dot: Vectors not of the same size\n";
    return 0.0;
  }
  if (n1<=0) return 0.0;
  
  Real* el1=store; Real* el2=y.store;
  Real sum = ddot(n1,el1,1,el2,1);
  return sum;
}



