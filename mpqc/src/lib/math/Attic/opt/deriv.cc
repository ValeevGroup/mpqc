static char rcsid[] = "$Header$";
//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#define WANT_MATH
#define WANT_STREAM
#include <stdio.h>
#include <string.h>
#include "opt.h"
#include <math/newmat7/cblas.h>
#include <math/newmat7/precisio.h>

int OptCGLike::CheckDeriv() // Check the analytic gradient with FD gradient
{return 0;}

int OptCG::CheckDeriv() // Check the analytic gradient with FD gradient
{return 0;}

int OptNewtonLike::CheckDeriv() // Check the analytic gradient with FD gradient
{return 0;}

int OptNewton::CheckDeriv() // Check the analytic gradient with FD gradient
{
  Real mcheps = FloatingPointPrecision::Epsilon();

  int i;
  int retcode = TRUE;
  int n = dim;

  Real maxerr; Real third = 0.33333;
  ColumnVector grad(n), fd_grad(n), error(n);
  SymmetricMatrix Hess(n), FDHess(n), ErrH(n);


  fd_grad = nlp->FDGrad(sx);    // Evaluate gradient using finite differences
  grad = nlp->GetGrad(); // Now use the analytical functions

  double gnorm = grad.NormInfinity();
  double eta   = pow(mcheps,third)*max(1,gnorm);

  printf("Check_Deriv: Checking gradients versus finite-differences\n");
  printf("    i    gradient     fd grad       error\n");
  for (i=1; i<=n; i++) {
    error(i) = fabs(grad(i)-fd_grad(i));
    printf("%5d %12.4e %12.4e %12.4e\n",
	   i, grad(i), fd_grad(i), error(i));
  }
  maxerr = error.NormInfinity();		
  printf("maxerror = %12.4e, tolerance =  %12.4e\n",
	 maxerr, eta);
  if (maxerr > eta) retcode = FALSE;

  printf("\nCheck_Deriv: Checking Hessian versus finite-differences\n");
  FDHess = nlp->FDHessian(sx); 
  Hess   = nlp->GetHess();
  ErrH   = Hess - FDHess;
  Print(ErrH);
  maxerr = ErrH.NormInfinity();
  printf("maxerror = %12.4e, tolerance =  %12.4e\n",
	 maxerr, eta);
  if (maxerr > eta) retcode = FALSE;

  return retcode;
}

int OptQNewton::CheckDeriv() // Check the analytic gradient with FD gradient
{return 0;}

// 
// int OptFDNewton::CheckDeriv() // Check the analytic gradient with FD gradient
// {return 0;}
//

//----------------------------------------------------------------------------
// Evaluate the Hessian using finite differences
// Assume that analytical gradients are available 
//----------------------------------------------------------------------------

SymmetricMatrix NLP1::FDHessian(ColumnVector sx) 
{
  Tracer trace("Optimize::FDHessian");
  Real mcheps = FloatingPointPrecision::Epsilon();
  Real sqrteps = sqrt(mcheps);

  int i;
  double hi;
  double xtmp;

  int nr = GetDim();

  ColumnVector gx(nr), grad(nr), xc(nr);
  Matrix Htmp(nr,nr);
  SymmetricMatrix H(nr);
   
  xc = GetXc();
  gx = GetGrad();
  
  for (i=1; i<=nr; i++) {

    hi = sqrteps*max(fabs(xc(i)),sx(i));
    copysign(hi,xc(i));
    xtmp = xc(i);
    xc(i) = xtmp + hi;
    grad = EvalG(xc);
    Htmp.Column(i) << (grad - gx) / hi;
    xc(i) = xtmp;
  }
  
  H << (Htmp.t() + Htmp)/2.0;
  return H;
}

//----------------------------------------------------------------------------
// Evaluate the Hessian using finite differences
// No analytical gradients available so use function values
//----------------------------------------------------------------------------

SymmetricMatrix NLP0::FD2Hessian(ColumnVector sx) 
{
  Tracer trace("Optimize::FD2Hessian");
  Real mcheps = FloatingPointPrecision::Epsilon();

  int i;
  int nr = GetDim();

  Real eta = pow(mcheps,0.333333);
  double xtmpi, xtmpj;
  double fii, fij, fx;
  
  ColumnVector fhi(nr), step(nr);
  SymmetricMatrix H(nr);

  xc = GetXc();
  fx = GetF();
  
  for (i=1; i<=nr; i++) {
    step(i) = eta*max(fabs(xc(i)),sx(i)); copysign(step(i),xc(i));
    xtmpi = xc(i);
    xc(i) = xtmpi + step(i);
    fhi(i) = EvalF(xc);
    xc(i) = xtmpi;
  }
  
  for (i=1; i<=nr; i++) {
    xtmpi = xc(i);
    xc(i) = xc(i) + step(i)*2.0;
    fii = EvalF(xc); 
    H(i,i) = ((fx - fhi(i)) + (fii - fhi(i))) / (step(i)*step(i));
    xc(i) = xtmpi + step(i);
    for (int j=i+1; j<=nr; ++j) {
      xtmpj = xc(j);
      xc(j) = xc(j) + step(j);
      fij = EvalF(xc);
      H(i,j) = ((fx - fhi(i)) + (fij - fhi(j))) / (step(i)*step(j));
      xc(j) = xtmpj;
    }
    xc(i) = xtmpi;
  } 
  return H;
}
ColumnVector NLP0::FDGrad(ColumnVector sx) // Compute gradient using finite differences
{
  Real mcheps = FloatingPointPrecision::Epsilon();
  double sqrteps = sqrt(mcheps);

  int i, n;
  double xtmp, fplus, fx;
  double hi;
  
  n = dim;
  ColumnVector grad(n);

  fx = GetF();

  for (i=1; i<=n; i++) {
    hi = sqrteps*max(fabs(xc(i)),sx(i));
    copysign(hi,xc(i));
    xtmp = xc(i);
    xc(i) = xtmp + hi;
    fplus = EvalF(xc);
    grad(i) = (fplus - fx) / hi;
    xc(i) = xtmp;
  }
  return grad;
}

