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
#include "opt.h"
#include <math/newmat7/cblas.h>

void OptCG::optimize()
// Nonlinear Conjugate Gradient
// 
// Given a nonlinear operator objfcn find the minimizer using a
// nonlinear conjugate gradient method
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
     
{

  int nlcg_iter;
  int maxfev, nfev;
  int convgd = 0;

  double gnorm, ggprev,beta, gg;
  double ftol = 1.e-1;
  double stp;
  double zero = 0.;
  
// Allocate local vectors 

  int n = dim;
  int maxiter;
  double fvalue;
  ColumnVector search(n), grad(n);

// Evaluate Function

  if (nlp.nonnull()) nlp->Eval();
  else {
      fprintf(stderr,"OptCG::optimize: nlp is null\n");
      abort();
    }

#if DEBUG
  objfcn.Print_State("nlcg: Initial Guess");
#endif


// Initialize iteration

  xprev  = nlp->GetXc();
  fprev  = fvalue = nlp->GetF();
  grad   = nlp->GetGrad();
  gprev  = grad;

  search = -grad;
  gg     = Dot(gprev,gprev), ggprev = gg, gnorm = sqrt(gg);

  printf("\n\t\tNonlinear CG");
  printf("\nIter\t   beta \t F(x) \t   ||grad|| \t ||step||\n");
  printf(" %5d %12.4e %12.4e %12.4e \n",0,0.0,fvalue,gnorm);

  maxiter = tol->GetMaxIter();
  for (nlcg_iter=1; nlcg_iter <= maxiter; nlcg_iter++) {

    iter_taken = nlcg_iter;

    maxfev = tol->GetMaxFeval();
    stp    = 1.0;

//  mcsrch will do a line search in the direction search from the current
//  point and update all the relevant information.
//  the step taken will be returned in stp

    if (mcsrch(nlp, search, stp, nfev, maxfev, ftol) < 0) {
      const char *msg = "Error in line search mcsrch";
      printf(msg);
      SetMesg(msg);
      ret_code = 9;
      return;
    }
    
    fvalue = nlp->GetF();
    grad = nlp->GetGrad();
    gg    = Dot(grad,grad), gnorm = sqrt(gg);
    beta  = max(zero,gg/ggprev);

    printf(" %5d %12.4e %12.4e %12.4e %12.4e\n",
	   nlcg_iter,beta,fvalue,gnorm,stp);
    
    convgd = CheckConvg();
    if (convgd > 0) {
        printf("Converged.\n");
      ret_code = convgd;
      return;
    }

// Update search direction and norms 
  
    xprev = nlp->GetXc();
    fprev = fvalue;
    gprev = grad;

    search = -grad + search*beta;
    ggprev = gg;

  }

  const char *msg = "Maximum number of iterations in nlcg";
  printf(msg);
  SetMesg(msg);
  ret_code = 4;
}
