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

//------------------------------------------------------------------------
//
//   Finite-Difference Newton Method member functions
//
//------------------------------------------------------------------------

void OptFDNewton::init(HessianUpdate*u)
{
  strcpy(method,"Finite-Difference Newton");
  if (!u) {
      update = new BFGSupdate();
      must_delete_update = 1;
    }
  else {
      udpate = u;
      must_delete_update = 0;
   }
}

OptFDNewton::~OptFDNewton()
{
  if (must_delete_update) delete update;
}

OptFDNewton::OptFDNewton()
{
  init(0);
}

OptFDNewton::OptFDNewton(NLP1* p):
  OptNewtonLike(p->GetDim()),
  nlp(p)
{
  init(0);
}

OptFDNewton::OptFDNewton(NLP1* p, TOLS* t):
  OptNewtonLike(p->GetDim(),t),
  nlp(p)
{
  init(0);
}

OptFDNewton::OptFDNewton(HessianUpdate*u,NewtonStep*s):
  OptNewtonLike(s)
{
  init(u);
}

OptFDNewton::OptFDNewton(NLP1* p,HessianUpdate*u,NewtonStep*s):
  OptNewtonLike(p->GetDim(),s),
  nlp(p)
{
  init(u);
}

OptFDNewton::OptFDNewton(NLP1* p, TOLS* t,HessianUpdate*u,NewtonStep*s):
  OptNewtonLike(p->GetDim(),t,s),
  nlp(p)
{
  init(u);
}

void OptFDNewton::DefProblem(NLP1 *p) 
{
  nlp =  p; 
  dim = nlp->GetDim();
  ColumnVector xtmp(dim), sxtmp(dim), sfxtmp(dim), gtmp(dim);
  SymmetricMatrix Htmp(dim);

  sxtmp = 1.0;
  xtmp  = 1.0;

  sx    = sxtmp;
  sfx   = sfxtmp;

  xprev = xtmp;
  fprev = 1.0e30;
  step_length = 1.0;
  gprev = gtmp;
  Hessian = Htmp;

  ret_code = iter_taken = fcn_evals = grad_evals = 0;
}

void OptFDNewton::optimize()
//---------------------------------------------------------------------------- 
// Finite-Difference Newton Method
// 
// Given a nonlinear operator nlp find the minimizer using a
// Finite-Difference Newton method
//
// Uses
//   mcsrch, mcstep  More and Thuente line search routine
//
//
// Based on the newmat Matrix Class Package
//
// Juan C. Meza
// Sandia National Laboratories
// meza@sandia.llnl.gov
//
//---------------------------------------------------------------------------- 
{

  int k;
  int maxfev, maxiter, nfev;
  int convgd = 0;

  double fvalue;
  double gnorm;
  double ftol = 1.e-4;
  double stp_length;
  double zero = 0.;
  
// Allocate local vectors 

  int n = dim;
  ColumnVector sk(n);
  SymmetricMatrix Hk(n);
  LowerTriangularMatrix L(n);

// Initialize iteration
// Evaluate Function, Gradient, and Hessian

  nlp->Eval();
  xprev   = nlp->GetXc();
  fvalue  = fprev = nlp->GetF();
  gprev   = nlp->GetGrad();  gnorm = Norm2(gprev);
  update->Initial(Hessian,fvalue,xprev,gprev);
  Hk = Hessian;
#if DEBUG
  nlp->PrintState("newton: Initial Guess");
#endif

  printf("\n\t\t%s Method",method);
  printf("\nIter\t  ||step||\t F(x) \t   ||grad|| \n");
  printf(" %5d %12.4e %12.4e %12.4e\n",0,0.0,fvalue,gnorm);

// Check for convergence. Need to take into account that this is the
// zeroth iteration
//  convgd = objfcn.Check_Convg(tol,stat);
//  if (convgd > 0) {
//    stat.ret_code = convgd;
//    return;
//  }
  
  maxiter = tol->GetMaxIter();
  for (k=1; k <= maxiter; k++) {

    iter_taken = k;

//  Solve for the Newton direction
//  H * step = -grad;

    stepper->step(sk,gprev,Hk);

//  Find an appropriate steplength
//  mcsrch will do a line search in the direction sk from the current
//  point and update all the relevant information.
//  the step length taken will be returned in step_length

    maxfev     = tol->GetMaxFeval();
    stp_length = 1.0;
    if (mcsrch(nlp, sk, stp_length, nfev, maxfev, ftol) < 0) {
      SetMesg("Error in line search mcsrch");
      ret_code = 9;
      return;
    }


//  At this point xc, fvalue and grad should have the new information
//  corresponding to the point that mcsrch found
//  Compute some quantities that will be used later 

    fvalue = nlp->GetF();
    gnorm = Norm2(nlp->GetGrad());
    printf(" %5d %12.4e %12.4e %12.4e\n",
	   k,stp_length,fvalue,gnorm);
    
//  Test for Convergence

    convgd = CheckConvg();
    if (convgd > 0) {
      ret_code = convgd;
      return;
    }

// Update all the variables

    step_length = stp_length;

    update->Update(Hessian,fvalue,nlp->GetXc(),xprev,nlp->GetGrad(),gprev);
    Hk = Hessian;

    xprev = nlp->GetXc();
    fprev = nlp->GetF();
    gprev = nlp->GetGrad();

  }

  SetMesg("Maximum number of iterations in newton");
  ret_code = 4;
}

int OptFDNewton::CheckConvg() // Check convergence
{
  int    n;
  double stol, ftol, rftol, gtol, rgtol;
  double xnorm, snorm, gnorm;
  ColumnVector xc, grad;

  double step_tol, fvalue;

  n  = nlp->GetDim();
  xc = nlp->GetXc();
  fvalue = nlp->GetF();

  xnorm =  Norm2(xc);
  
// Test 1. step tolerance 

  ColumnVector step(n);
  step = xc - xprev;
  step_tol = tol->GetStepTol();
  snorm = Norm2(step);
  stol  = step_tol*max(1.0,xnorm);
  if (snorm  <= stol) {
    strcpy(mesg,"Step tolerance test passed");
    printf("CheckConvg: snorm = %12.4e, stol = %12.4e\n", snorm, stol);
    return 1;
  }
  
// Test 2. function tolerance
  Real deltaf = fprev - fvalue;
  ftol = tol->GetFtol();
  rftol = ftol*max(1.0,fabs(fvalue));
  if (deltaf <= rftol) {
    strcpy(mesg,"Function tolerance test passed");
    printf("CheckConvg: deltaf = %12.4e, ftol = %12.4e\n", deltaf, ftol);
    return 2;
  }
  

// Test 3. gradient tolerance 

  grad = nlp->GetGrad();
  gtol = tol->GetGtol();
  rgtol = gtol*max(1.0,fabs(fvalue));
  gnorm = Norm2(grad);
  if (gnorm <= rgtol) {
    strcpy(mesg,"Gradient tolerance test passed");
    printf("CheckConvg: gnorm = %12.4e, gtol = %12.4e\n", gnorm, gtol);
    return 3;
  }
  

// Test 4. absolute gradient tolerance 

  if (gnorm <= gtol) {
    strcpy(mesg,"Gradient tolerance test passed");
    printf("CheckConvg: gnorm = %12.4e, gtol = %12.4e\n", gnorm, gtol);
    return 4;
  }
  
  // Nothing to report 

  return 0;

}

//---------------------------------------------------------------------------- 
//
// Update Hessian using a Finite-Difference approximation
//
//---------------------------------------------------------------------------- 
SymmetricMatrix OptFDNewton::UpdateH(SymmetricMatrix& Hk, int k) 
{
  Tracer trace("FDNewton::UpdateH");

#ifdef DEBUG
  fprintf(stderr, "Inside OptFDNewton::UpdateH"\n);
  Print(Hk);
#endif
  Real mcheps = FloatingPointPrecision::Epsilon();
  Real sqrteps = sqrt(mcheps);

  int i;
  double hi;
  double xtmp;

  int nr = nlp->GetDim();

  ColumnVector gx(nr), grad(nr), xc;
  Matrix Htmp(nr,nr);
   
  xc = nlp->GetXc();
  gx = nlp->GetGrad();
  
  for (i=1; i<=nr; i++) {
    hi = sqrteps*max(fabs(xc(i)),sx(i));
    copysign(hi,xc(i));
    xtmp = xc(i);
    xc(i) = xtmp + hi;
    nlp->SetX(xc);
    grad = nlp->EvalG(); 
    Htmp.Column(i) << (grad - gx) / hi;
    xc(i) = xtmp;
  }
  
  nlp->SetX(xc);
  Hk << (Htmp.t() + Htmp)/2.0;
//   Print(Hessian);
  Hessian = Hk;
  return Hk;
}

