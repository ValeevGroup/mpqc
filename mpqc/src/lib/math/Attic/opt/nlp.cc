static char rcsid[] = "$Header$";

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include <math.h>

#define WANT_MATH
#define WANT_STREAM
#include <stdio.h>
#include <float.h>
#include <string.h>
#include "nlp.h"
#include <math/newmat7/cblas.h>
#include <math/nihmatrix/nihmatrix.h>
#include <math/nihmatrix/newmat.h>
#include <util/keyval/keyval.h>

////////////////////////////////////////////////////////////////////////

#define CLASSNAME NLP0
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
NLP0::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}
void
NLP0::save_data_state(StateOut&s)
{
  s.put(dim);
  s.put(deriv_level);
  s.put(fvalue.computed());
  s.put(fvalue.compute());
  if (fvalue.computed()) {
      s.put(fvalue.result_noupdate());
    }
}
NLP0::NLP0(KeyVal&kv):
  fvalue(this)
{
  // set up the initial value for x
  DVector initx(PrefixKeyVal("initialx",kv));
  SetDim(initx.dim());
  ColumnVector x;
  Convert(initx,x);
  SetX(x);
}
NLP0::NLP0(StateIn&s):
  SavableState(s,class_desc_),
  fvalue(this)
{
  s.get(dim);
  s.get(deriv_level);
  s.get(fvalue.computed());
  s.get(fvalue.compute());
  if (fvalue.computed()) {
      s.get(fvalue.result_noupdate());
    }
  nfevals = 0;
}  
DescribedClass_REF_def(NLP0);
#define CLASSNAME NLP1
#define PARENTS public NLP0
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
NLP1::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { NLP0::_castdown(cd) };
  return do_castdowns(casts,cd);
}
void
NLP1::save_data_state(StateOut&s)
{
  NLP0::save_data_state(s);
  s.put(grad.computed());
  s.put(grad.compute());
  if (grad.computed()) {
      DVector d(GetDim());
      ColumnVector* v = &(grad.result_noupdate());
      for (int i=0; i<GetDim(); i++) {
          d[i] = v->element(i);
        }
      d.save_object_state(s);
    }
}
NLP1::NLP1(KeyVal&kv):
  NLP0(kv),
  grad(this)
{
}
NLP1::NLP1(StateIn&s):
  SavableState(s,class_desc_),
  NLP0(s),
  grad(this)
{
  s.get(grad.computed());
  s.get(grad.compute());
  if (grad.computed()) {
      DVector d(s);
      ColumnVector* v = &(grad.result_noupdate());
      v->ReDimension(GetDim());
      for (int i=0; i<GetDim(); i++) {
          v->element(i) = d(i);
        }
    }
  ngevals = 0;
}
DescribedClass_REF_def(NLP1);
#define CLASSNAME NLP2
#define PARENTS public NLP1
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
NLP2::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { NLP1::_castdown(cd) };
  return do_castdowns(casts,cd);
}
void
NLP2::save_data_state(StateOut&s)
{
  NLP1::save_data_state(s);
  s.put(Hessian.computed());
  s.put(Hessian.compute());
  if (Hessian.computed()) {
      DMatrix d(GetDim(),GetDim());
      SymmetricMatrix* v = &(Hessian.result_noupdate());
      for (int i=0; i<GetDim(); i++) {
          for (int j=0; j<=i; j++) {
              d(i,j) = v->element(i,j);
            }
        }
      d.save_object_state(s);
    }
}
NLP2::NLP2(KeyVal&kv):
  NLP1(kv),
  Hessian(this)
{
}
NLP2::NLP2(StateIn&s):
  SavableState(s,class_desc_),
  NLP1(s),
  Hessian(this)
{
  s.get(Hessian.computed());
  s.get(Hessian.compute());
  if (Hessian.computed()) {
      DMatrix d(s);
      SymmetricMatrix* v = &(Hessian.result_noupdate());
      v->ReDimension(GetDim());
      for (int i=0; i<GetDim(); i++) {
          for (int j=0; j<=i; j++) {
              v->element(i,j) = d(i,j);
            }
        }
    }
  nhevals = 0;
}
DescribedClass_REF_def(NLP2);

////////////////////////////////////////////////////////////////////////

NLP0::NLP0():
  dim(0),
  deriv_level(0),
  fvalue(this)
{
}

NLP0::NLP0(int ndim):
  dim(ndim),
  deriv_level(0),
  xc(ndim),
  fvalue(this)
{
  xc = 0;
}

double
NLP0::value()
{
  return fvalue;
}

NLP0::~NLP0() {}               // Destructor

void NLP0::SetDim(int ndim)
{
  dim = ndim;
  xc.ReDimension(ndim);
}

void NLP0::SetDerivLevel(int k) 
{ 
// Set the derivative level
// 0 means only function values are available
// 1 means function and 1st derivatives are available
// 2 means function, 1st and 2nd derivatives are available

  if ( k < 0 || k > 3) {
    printf("NLP0::SetDerivLevel. Derivative level out of bounds.\n");
    printf("NLP0::SetDerivLevel. Setting Derivative level = 0.\n");
    deriv_level = 0;
  }
  else deriv_level = k;
};

double
NLP0::EvalF() { return fvalue; }

//-------------------------------------------------------------------------
// EvalF() Evaluate f at xc
//-------------------------------------------------------------------------

void
NLF0::compute() // Evaluate function
{
  if (fvalue.compute()) {
      fcn(dim, xc, fvalue.result_noupdate());
      nfevals++;
    }
}
double NLF0::EvalF(ColumnVector x) // Evaluate function at x
{
  double fx;
  fcn(dim, x, fx);
  nfevals++;
  return fx;
}
//
//
//
void
NLF1::compute() // Evaluate function
{
  if (fvalue.compute()) {
      ColumnVector gtmp(dim);
      fcn(0, dim, xc, fvalue.result_noupdate(), gtmp);
      fvalue.computed() = 1;
      nfevals++;
    }
  if (grad.compute()) {
      double fx;
      fcn(1, dim, xc, fx, grad.result_noupdate());
      grad.computed() = 1;
      ngevals++;
    }
}
double NLF1::EvalF(ColumnVector x) // Evaluate function at x
{
  double fx;
  ColumnVector gtmp(dim);
  fcn(0, dim, x, fx, gtmp);
  nfevals++;
  return fx;
}
ColumnVector NLF1::EvalG(ColumnVector x) // Evaluate the gradient at x
{
  double fx;
  ColumnVector gx(dim);

  fcn(1, dim, x, fx, gx);
  ngevals++;
  return gx;
}

//
//
//
void
NLF2::compute() // Evaluate function
{
  if (fvalue.compute()) {
      ColumnVector gtmp(dim);
      SymmetricMatrix Htmp(dim);

      fcn(0, dim, xc, fvalue.result_noupdate(), gtmp, Htmp);
      nfevals++;
    }
  if (grad.compute()) {
      double fx;
      SymmetricMatrix Htmp(dim);

      fcn(1, dim, xc, fx, grad.result_noupdate(), Htmp);
      ngevals++;
    }
  if (Hessian.compute()) {
      double fx;
      ColumnVector gtmp(dim);

      fcn(2, dim, xc, fx, gtmp, Hessian.result_noupdate());
      nhevals++;
    }
}
double NLF2::EvalF(ColumnVector x) // Evaluate function at x
{
  double fx;
  ColumnVector gtmp(dim);
  SymmetricMatrix Htmp(dim);

  fcn(0, dim, x, fx, gtmp, Htmp);
  nfevals++;
  return fx;
}
ColumnVector NLF2::EvalG(ColumnVector x) // Evaluate the gradient at x
{
  double fx;
  ColumnVector gx(dim);
  SymmetricMatrix Htmp(dim);

  fcn(1, dim, x, fx, gx, Htmp);
  ngevals++;
  return gx;
}

//-------------------------------------------------------------------------
// Output Routines
//-------------------------------------------------------------------------

void NLP0::PrintState(char * s) 
{ // Print out current state: x current, gradient and function value
  printf("\nNLP0: %s\n",s);
  printf("\n    i\t    xc \t\t grad \n");
  for (int i=1; i<=dim; i++) printf("%6d %12.4e\n", 
				   i, xc(i));
  printf("f(x) = %12.4e \n\n", fvalue.result());
};

NLP1::NLP1():
  grad(this)
{
  deriv_level=1;
}

NLP1::NLP1(int ndim):
  NLP0(ndim),
  grad(this)
{
  grad.result_noupdate().ReDimension(ndim);
  deriv_level = 1;
}

NLP1::~NLP1()
{
}

void NLP1::SetDim(int ndim)
{
  NLP0::SetDim(ndim);
  grad.result_noupdate().ReDimension(ndim);
}

void NLP1::PrintState(char * s) 
{ // Print out current state: x current, gradient and function value
  printf("\nNLP1: %s\n",s);
  printf("\n    i\t    xc \t\t grad \n");
  for (int i=1; i<=dim; i++) printf("%6d %12.4e %12.4e\n", 
				   i, xc(i),grad.result()(i));
  printf("f(x) = %12.4e \n\n", fvalue.result());
};

NLP2::NLP2():
  Hessian(this)
{
  deriv_level = 2;
}

NLP2::NLP2(int ndim):
  NLP1(ndim),
  Hessian(this)
{
  Hessian.result_noupdate().ReDimension(ndim);
  deriv_level = 2;
}

NLP2::~NLP2()
{
}

void NLP2::SetDim(int ndim)
{
  NLP1::SetDim(ndim);
  Hessian.result_noupdate().ReDimension(ndim);
}

void NLP2::PrintState(char * s) 
{ // Print out current state: x current, gradient and function value
  printf("\nNLP2: %s\n",s);
  printf("\n    i\t    xc \t\t grad \n");
  for (int i=1; i<=dim; i++) printf("%6d %12.4e %12.4e\n", 
				   i, xc(i),grad.result()(i));
  printf("f(x) = %12.4e \n\n", fvalue.result());
};

//-------------------------------------------------------------------------
// Routines for OptProgress
//-------------------------------------------------------------------------

OptProgress::OptProgress(double f,double oldf,
              const ColumnVector&x,const ColumnVector&oldx,
              const ColumnVector&g)
{
  ColumnVector s = x - oldx;
  _adeltaf = fabs(f-oldf);
  _agnorm = Norm2(g);
  _asnorm = Norm2(s);
  _xnorm = Norm2(x);
  _fvalue = f;
  ResetStatus();
}

void OptProgress::ResetStatus() { SetStatus(""); }
void OptProgress::SetStatus(const char*m) { _status = m; }
const char* OptProgress::GetStatus() { return _status; }
double OptProgress::XNorm() { return _xnorm; }
double OptProgress::ASNorm() { return _asnorm; }
double OptProgress::RSNorm() { return _asnorm/max(1.0,_xnorm); }
double OptProgress::ADeltaF() { return _adeltaf; }
double OptProgress::RDeltaF() { return _adeltaf/max(1.0,fabs(_fvalue)); }
double OptProgress::AGNorm() { return _agnorm; }
double OptProgress::RGNorm() { return _agnorm/max(1.0,fabs(_fvalue)); }

//-------------------------------------------------------------------------
// Routines for Tolerance setting 
//-------------------------------------------------------------------------

DescribedClass_REF_def(TOLS);

#define CLASSNAME TOLS
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
TOLS::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}
TOLS::TOLS()
{
  SetDefaultTol();
}
TOLS::TOLS(KeyVal&kv)
{
  SetDefaultTol();
#define getdoubleTOL(name) \
  { \
  double tmp = kv.doublevalue( # name ); \
  if (kv.error() == KeyVal::OK) name = tmp; \
  }
#define getintTOL(name) \
  { \
  int tmp = kv.intvalue( # name ); \
  if (kv.error() == KeyVal::OK) name = tmp; \
  }
#define getbooleanTOL(name) \
  { \
  int tmp = kv.booleanvalue( # name ); \
  if (kv.error() == KeyVal::OK) name = tmp; \
  }
  getdoubleTOL(max_step);
  getintTOL(max_iter);
  getintTOL(max_feval);
  getdoubleTOL(step_tol);
  getdoubleTOL(fcn_tol);
  getdoubleTOL(grad_tol);
  getbooleanTOL(use_ftol);
  getbooleanTOL(use_gtol);
  getbooleanTOL(use_agtol);
  getbooleanTOL(use_steptol);
  getbooleanTOL(or_tol);
#undef getdoubleTOL
#undef getintTOL
#undef getbooleanTOL
}
TOLS::TOLS(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(mcheps);
  s.get(max_step);
  s.get(max_iter);
  s.get(max_feval);
  s.get(step_tol);
  s.get(fcn_tol);
  s.get(grad_tol);
  s.get(use_ftol);
  s.get(use_gtol);
  s.get(use_agtol);
  s.get(use_steptol);
  s.get(or_tol);
}
void
TOLS::save_data_state(StateOut&s)
{
  s.put(mcheps);
  s.put(max_step);
  s.put(max_iter);
  s.put(max_feval);
  s.put(step_tol);
  s.put(fcn_tol);
  s.put(grad_tol);
  s.put(use_ftol);
  s.put(use_gtol);
  s.put(use_agtol);
  s.put(use_steptol);
  s.put(or_tol);
}

void TOLS::SetDefaultTol()  // Set the default tolerances
{
  mcheps     = DBL_EPSILON;
  max_step   = 100.;
  max_iter   = 25;
  max_feval  = 250;
  step_tol   = mcheps;
  fcn_tol    = sqrt(mcheps);
  grad_tol   = pow(mcheps,1.0/3.0);
  use_ftol   = 1;
  use_gtol   = 1;
  use_agtol  = 1;
  use_steptol= 1;
  or_tol     = 1;
}

void TOLS::SetFtol(double x) {fcn_tol  = x;} // Set the function tolerance
void TOLS::SetGtol(double x) {grad_tol = x;} // Set the gradient tolerance
void TOLS::SetStepTol(double x) {step_tol = x;} // Set the step tolerance
void TOLS::SetMaxIter(int k) {max_iter = k;} // Set the max number of iters
void TOLS::SetMaxFeval(int k) {max_feval = k;} // Set the max number of iters


void TOLS::SetUseFtol(int a) { use_ftol = a; }
void TOLS::SetUseGtol(int a) { use_gtol = a; }
void TOLS::SetUseAGtol(int a) { use_agtol = a; }
void TOLS::SetUseStepTol(int a) { use_steptol = a; }
int  TOLS::UseFtol() { return use_ftol; }
int  TOLS::UseGtol() { return use_gtol; }
int  TOLS::UseAGtol() { return use_agtol; }
int  TOLS::UseStepTol() { return use_steptol; }
void TOLS::SetOrTol(int a) { or_tol = a; }
void TOLS::SetAndTol(int a) { or_tol = !a; }
int  TOLS::OrTol() { return or_tol; }
int  TOLS::AndTol() { return !or_tol; }


void TOLS::PrintTol()  // Set the default tolerances
{
  printf("Tolerances\n");
  printf("machine epsilon    = %e\n",mcheps);
  printf("maximum step       = %e\n",max_step);
  printf("maximum iter       = %d\n",max_iter);
  printf("maximum fcn eval   = %d\n",max_feval);
  printf("step tolerance     = %e\n",step_tol);
  printf("function tolerance = %e\n",fcn_tol);
  printf("gradient tolerance = %e\n",grad_tol);

}

int TOLS::Converged(OptProgress& op)
{
  op.ResetStatus();

// Test 1. step tolerance 
  int StepTolMet;
  if (UseStepTol()) {
    if (op.RSNorm() <= GetStepTol()) {
      op.SetStatus("step tolerance test passed");
      printf("CheckConvg: snorm = %12.4e, stol = %12.4e\n", op.RSNorm(), GetStepTol());
      if (OrTol()) return 1;
      StepTolMet = 1;
    }
    else StepTolMet = 0;
  }
  else StepTolMet = 1;
  
// Test 2. function tolerance
  int FTolMet;
  if (UseFtol()) {
    if (op.RDeltaF() <= GetFtol()) {
      op.SetStatus("function tolerance test passed");
      printf("CheckConvg: deltaf = %12.4e, ftol = %12.4e\n", op.RDeltaF(), GetFtol());
      if (OrTol()) return 2;
      FTolMet = 1;
    }
    else FTolMet = 0;
  }
  else FTolMet = 1;

// Test 3. gradient tolerance 
  int GTolMet;
  if (UseGtol()) {
    if (op.RGNorm() <= GetGtol()) {
      op.SetStatus("gradient tolerance test passed");
      printf("CheckConvg: gnorm = %12.4e, gtol = %12.4e\n", op.RGNorm(), GetGtol());
      if (OrTol()) return 3;
      GTolMet = 1;
    }
    else GTolMet = 0;
  }
  else GTolMet = 1;

// Test 4. absolute gradient tolerance 
  int AGTolMet;
  if (UseAGtol()) {
    if (op.AGNorm() <= GetGtol()) {
      op.SetStatus("absolute gradient tolerance test passed");
      printf("CheckConvg: gnorm = %12.4e, gtol = %12.4e\n", op.AGNorm(), GetGtol());
      if (OrTol()) return 4;
      AGTolMet = 1;
    }
    else AGTolMet = 0;
  }
  else AGTolMet = 1;
  
  if (StepTolMet && FTolMet && GTolMet && AGTolMet) {
    op.SetStatus("all tolerance tests passed");
    return 1;
  }

  return 0;
}


// added by clj
double NLP0::EvalF(ColumnVector x)
{
  ColumnVector xtmp = xc;
  xc = x;
  return EvalF();
}

ColumnVector NLP1::EvalG(ColumnVector x)
{
  ColumnVector xtmp = xc;
  xc = x;
  return EvalG();
}

// these are inlines that were removed from nlp.h
  void NLP0::Eval() {double tmp = EvalF();}   
  void NLP0::SetX(int i, double x)  {xc(i) = x; obsolete(); }
  void NLP0::SetX(ColumnVector& x) {xc = x; obsolete(); }
  int NLP0::GetDim() {return dim;}
  int NLP0::GetFevals() {return nfevals;}
  ColumnVector NLP0::GetXc() {return xc;}
  const ColumnVector& NLP0::GetX() const { return xc; }
  double NLP0::GetF() {return fvalue;}

  void NLP1::Eval() {double tmp = EvalF(); ColumnVector gtmp = EvalG();} 
  ColumnVector NLP1::GetGrad() {return grad;}
  int NLP1::GetGevals() {return ngevals;}
const ColumnVector&
NLP1::gradient()
{
  return grad;
}
ColumnVector
NLP1::EvalG()
{
  return grad;
}

  void NLP2::Eval() {double tmp = EvalF(); ColumnVector gtmp = EvalG(); 
		       SymmetricMatrix Htmp = EvalH();}
  SymmetricMatrix NLP2::GetHess() {return Hessian;}
  int NLP2::GetHevals() {return nhevals;}
const SymmetricMatrix&
NLP2::hessian()
{
  return Hessian;
}
SymmetricMatrix
NLP2::EvalH()
{
  return Hessian;
}

  NLF0::NLF0() {}  // Constructor
  NLF0::~NLF0() {}               // Destructor
  NLF0::NLF0(int ndim):NLP0(ndim){}
  NLF0::NLF0(int ndim, USERFCN0 f):NLP0(ndim){fcn = f;}
  void NLF0::SetF(USERFCN0 f)    {fcn = f;} 

  NLF1::NLF1() {}  // Constructor
  NLF1::NLF1(int ndim): NLP1(ndim){}
  NLF1::NLF1(int ndim, USERFCN1 f):NLP1(ndim){fcn = f;}
  NLF1::~NLF1() {}                     // Destructor
  void NLF1::SetF(USERFCN1 f)    {fcn = f;} 

  NLF2::NLF2() {}                      // Constructor
  NLF2::NLF2(int ndim): NLP2(ndim){}
  NLF2::NLF2(int ndim, USERFCN2 f): NLP2(ndim) {fcn = f;}
  NLF2::~NLF2() {}                     // Destructor
  void NLF2::SetF(USERFCN2 f)    {fcn = f;} 

  double  TOLS::GetStepTol() {return step_tol;}
  double  TOLS::GetFtol() {return fcn_tol;}
  double  TOLS::GetGtol() {return grad_tol;}
  int     TOLS::GetMaxIter() {return max_iter;}
  int     TOLS::GetMaxFeval() {return max_feval;}
