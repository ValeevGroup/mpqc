static char rcsid[] = "$Header$";

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#define WANT_STREAM
#define WANT_MATH
#include <stdio.h>
#include <string.h>
#include "opt.h"
#include <math/newmat7/cblas.h>

//------------------------------------------------------------------------
//
//   Newton Method member functions
//
//------------------------------------------------------------------------

OptNewton::OptNewton()
{strcpy(method,"Newton");}

OptNewton::OptNewton(NLP2* p):
  OptNewtonLike(p->GetDim()),
  nlp(p)
{strcpy(method,"Newton");}

OptNewton::OptNewton(NLP2* p, TOLS* t):
  OptNewtonLike(p->GetDim(),t),
  nlp(p)
{strcpy(method,"Newton");}

OptNewton::OptNewton(NewtonStep*s):
  OptNewtonLike(s)
{strcpy(method,"Newton");}

OptNewton::OptNewton(NLP2* p,NewtonStep*s):
  OptNewtonLike(p->GetDim(),s),
  nlp(p)
{strcpy(method,"Newton");}

OptNewton::OptNewton(NLP2* p, TOLS* t,NewtonStep*s):
  OptNewtonLike(p->GetDim(),t,s),
  nlp(p)
{strcpy(method,"Newton");}

void OptNewton::optimize()
//---------------------------------------------------------------------------- 
// Newton Method
// 
// Given a nonlinear operator nlp find the minimizer using a
// Newton method
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

// Initialize iteration
// Evaluate Function, Gradient, and Hessian

  nlp->Eval();
  xprev  = nlp->GetXc();
  fvalue = fprev = nlp->GetF();
  gprev  = nlp->GetGrad();  gnorm = Norm2(gprev);
  Hk     = nlp->EvalH();

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
    if (mcsrch(nlp.pointer(), sk, stp_length, nfev, maxfev, ftol) < 0) {
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

    Hk = nlp->EvalH();

    xprev = nlp->GetXc();
    fprev = nlp->GetF();
    gprev = nlp->GetGrad();

  }

  SetMesg("Maximum number of iterations in newton");
  ret_code = 4;
}

int OptNewton::CheckConvg() // Check convergence
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

