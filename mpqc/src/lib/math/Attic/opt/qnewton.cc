static char rcsid[] = "$Header$";

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

//
// I took out scaling because it was fully implemented and I wasn't sure
// exactly how it should be implemented when I rearranged the code
//

#define WANT_STREAM
#define WANT_MATH
#include <stdio.h>
#include <string.h>
#include "opt.h"
#include <math/newmat7/cblas.h>
#include <math/newmat7/precisio.h>
//------------------------------------------------------------------------
//
//   Quasi-Newton Method member functions
//
//------------------------------------------------------------------------
void OptQNewton::init(RefHessianUpdate u)
{
  strcpy(method,"Quasi-Newton");
  if (u.null()) {
      update = new BFGSupdate();
    }
  else {
      update = u;
   }
}

OptQNewton::~OptQNewton()
{
}

OptQNewton::OptQNewton()
{
  init(0);
}

OptQNewton::OptQNewton(NLP1* p):
  OptNewtonLike(p->GetDim()),
  nlp(p)
{
  init(0);
}

OptQNewton::OptQNewton(NLP1* p, TOLS* t):
  OptNewtonLike(p->GetDim(), t),
  nlp(p)
{
  init(0);
}

OptQNewton::OptQNewton(HessianUpdate*u,NewtonStep*s):
  OptNewtonLike(s)
{
  init(u);
}

OptQNewton::OptQNewton(NLP1* p,HessianUpdate*u,NewtonStep*s):
  OptNewtonLike(p->GetDim(),s),
  nlp(p)
{
  init(u);
}

OptQNewton::OptQNewton(NLP1* p, TOLS* t,HessianUpdate*u,NewtonStep*s):
  OptNewtonLike(p->GetDim(), t, s),
  nlp(p)
{
  init(u);
}

void OptQNewton::DefProblem(NLP1 *p) 
{
  nlp =  p; 

  xprev = 0.0;
  fprev = 1.0e30;
  gprev = 0.0;
  Hessian = 0.0;

  ret_code = iter_taken = fcn_evals = grad_evals = 0;
}

void OptQNewton::optimize()
//---------------------------------------------------------------------------- 
// 
// Given a nonlinear operator nlp find the minimizer using a
// Quasi-Newton method
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
  update->Initial(Hessian,*nlp);
  Hk = Hessian;

#if DEBUG
  nlp->PrintState("newton: Initial Guess");
#endif

  printf("\n\t\t%s Method", method);
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

    update->Update(Hessian,*nlp);
    Hk = Hessian;

    xprev = nlp->GetXc();
    fprev = nlp->GetF();
    gprev = nlp->GetGrad();

  }

  SetMesg("Maximum number of iterations in newton");
  ret_code = 4;
}

int OptQNewton::CheckConvg() // Check convergence
{
  OptProgress op(nlp->GetF(),fprev,nlp->GetXc(),xprev,nlp->GetGrad());

  int converged = tol->Converged(op);

  strcpy(mesg,op.GetStatus());

  return converged;
}
