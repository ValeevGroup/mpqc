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
#include <util/keyval/keyval.h>

#define CLASSNAME OptCG
#define PARENTS public OptCGLike
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
OptCG::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { OptCGLike::_castdown(cd) };
  return do_castdowns(casts,cd);
}
OptCG::OptCG(StateIn&s):
  SavableState(s,class_desc_),
  OptCGLike(s)
{
  nlp = NLP1::restore_state(s);
}
void
OptCG::save_data_state(StateOut&s)
{
  OptCGLike::save_data_state(s);
  nlp->save_state(s);
}
OptCG::OptCG(KeyVal&kv):
  OptCGLike(kv)
{
  //nlp = kv.describedclassvalue("nlp");
  RefDescribedClass nlpdc = kv.describedclassvalue("nlp");
  printf("nlpdc.nonnull() = %d\n",nlpdc.nonnull());
  if (nlpdc.nonnull()) {
      printf("classname = \"%s\"\n",nlpdc->class_desc()->name());
    }
  nlp = nlpdc;
  dim = nlp->GetDim();
  printf("nlp.nonnull() = %d\n",nlp.nonnull());
}

OptCG::OptCG(NLP1* p): OptCGLike(p->GetDim()),nlp(p)
{
  printf("OptCG: Constructing an OptCG class, dim = %d\n",dim);
  strcpy(method,"CG");
}

OptCG::OptCG(NLP1* p, TOLS* t): OptCGLike(p->GetDim(),t),nlp(p)
{
  printf("OptCG: Constructing an OptCG class, dim = %d\n",dim);
  strcpy(method,"CG");
}

void OptCG::PrintStatus() // Set Message
{
  printf("%s\n",mesg);

  printf("Dimension of the problem = %d\n",dim);
  printf("Optimization method      = %s\n",method);

  printf("Return code = %d\n",ret_code);
  printf("Number of iterations taken    = %d\n",iter_taken);
  printf("Number of function evalutions = %d\n",fcn_evals);
  printf("Number of gradient evalutions = %d\n",grad_evals);

/*   Print(Hessian); */
}


int OptCG::CheckConvg() // Check convergence
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
    printf("snorm = %12.4e, stol = %12.4e\n", snorm, stol);
    return 1;
  }
  
// Test 2. function tolerance
  Real deltaf = fprev - fvalue;
  ftol = tol->GetFtol();
  rftol = ftol*max(1.0,fabs(fvalue));
  if (deltaf <= rftol) {
    strcpy(mesg,"Function tolerance test passed");
    printf("deltaf = %12.4e, ftol = %12.4e\n", deltaf, ftol);
    return 2;
  }
  

// Test 3. gradient tolerance 

  grad = nlp->GetGrad();
  gtol = tol->GetGtol();
  rgtol = gtol*max(1.0,fabs(fvalue));
  gnorm = Norm2(grad);
  if (gnorm <= rgtol) {
    strcpy(mesg,"Gradient tolerance test passed");
    printf("gnorm = %12.4e, gtol = %12.4e\n", gnorm, gtol);
    return 3;
  }
  

// Test 4. absolute gradient tolerance 

  if (gnorm <= gtol) {
    strcpy(mesg,"Gradient tolerance test passed");
    printf("gnorm = %12.4e, gtol = %12.4e\n", gnorm, gtol);
    return 4;
  }
  
  // Nothing to report 

  return 0;

}

