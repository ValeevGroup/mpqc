#ifndef _math_opt_opt_h
#define _math_opt_opt_h

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <math/opt/nlp.h>
#include <math/newmat7/newmatap.h>
#include <math/opt/tstfcn.h>

//  Some useful typdefs

//typedef void (*USERFCN)(int, int, ColumnVector, double&, 
//			ColumnVector&, SymmetricMatrix&);

//------------------------------------------------------------------------
//
// Base Optimization Class
//
//------------------------------------------------------------------------

class Optimize: virtual public SavableState {
#   define CLASSNAME Optimize
#   include <util/state/stated.h>
#   include <util/class/classd.h>
protected:
  int  dim;
  RefTOLS tol;
  
  ColumnVector sx;       // Diagonal Scaling Matrix for x
  ColumnVector sfx;      // Diagonal Scaling Matrix for f
  
  ColumnVector xprev;
  double       fprev;
  double       step_length;

  char       method[80];  // What method is being used
  char       mesg[80];    // Optional message
  int        ret_code;    // Return code from optimization class
  int        iter_taken;      // number of iterations taken 
  int        fcn_evals;         // number of function evaluations taken
  
public:
  FILE *outfile;

  Optimize();
  Optimize(int n);
  Optimize(TOLS* t);
  Optimize(int n, TOLS* t);
  virtual ~Optimize();
    
    Optimize(StateIn&);
    Optimize(KeyVal&);
    void save_data_state(StateOut&);

// Problem definition
  
  virtual void SetTols(TOLS *t);

// Set various properties
  
  void SetOutput(FILE *fp);
  void SetTol(TOLS *t);
  void SetMesg(const char *s);     // Set message 
  void SetMethod(const char *s);   // Set method of choice

  ColumnVector Getsx();
  ColumnVector Getsfx();

// Finally routines that each method will have to define for themselves
  
  virtual void  optimize()    = 0;
  virtual int   CheckConvg()  = 0;
  virtual void  PrintStatus() = 0;

};
DescribedClass_REF_dec(Optimize);

//------------------------------------------------------------------------
// Derived Classes from Optimize 
//------------------------------------------------------------------------

class OptDirect: public Optimize {
#   define CLASSNAME OptDirect
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 public:
  OptDirect();
  virtual ~OptDirect();
};

//----------------------------------------------------------------------
// Parallel Direct Search Method
//----------------------------------------------------------------------

class OptPDS: public OptDirect {
#   define CLASSNAME OptPDS
#   include <util/state/stated.h>
#   include <util/class/classd.h>
protected:
  RefNLP0 nlp;
 public:
  OptPDS();
  virtual ~OptPDS();
};

//----------------------------------------------------------------------
// CG-Like Methods
//----------------------------------------------------------------------

class OptCGLike: public Optimize {
#   define CLASSNAME OptCGLike
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 protected:
  ColumnVector gprev;
  int grad_evals;
 public:
  OptCGLike();
    OptCGLike(StateIn&);
    OptCGLike(KeyVal&);
    void save_data_state(StateOut&);
  OptCGLike(int n);
  OptCGLike(int n, TOLS* t);
  virtual ~OptCGLike();
  virtual int CheckConvg();
  virtual int CheckDeriv();
};

class OptCG: public OptCGLike {
#   define CLASSNAME OptCG
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 protected:
  RefNLP1 nlp;
  //ColumnVector gprev;
 public:
  OptCG();
    OptCG(KeyVal&);
    OptCG(StateIn&);
    void save_data_state(StateOut&);
  OptCG(NLP1* p);
  OptCG(NLP1* p, TOLS* t);
  virtual ~OptCG();

// These are defined elsewhere

  void optimize();
  int CheckConvg();
  int CheckDeriv();
  void PrintStatus();

  void DefProblem(NLP1 *p);
  NLP1* GetProblem();

};

//----------------------------------------------------------------------
// Helper class for doing the newton update
//----------------------------------------------------------------------

class NewtonStep: public virtual SavableState
{
#   define CLASSNAME NewtonStep
#   include <util/state/stated.h>
#   include <util/class/classd.h>    
  public:
    NewtonStep();
    virtual ~NewtonStep();
    NewtonStep(StateIn&);
    NewtonStep(KeyVal&);
    void save_data_state(StateOut&);
    virtual void step(ColumnVector&x,
                      ColumnVector&grad,
                      SymmetricMatrix&hess) = 0;
};
DescribedClass_REF_dec(NewtonStep);

class MCholeskyNewtonStep: public NewtonStep
{
#   define CLASSNAME MCholeskyNewtonStep
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>    
 public:
    MCholeskyNewtonStep();
    MCholeskyNewtonStep(StateIn&);
    MCholeskyNewtonStep(KeyVal&);
    void save_data_state(StateOut&);
    ~MCholeskyNewtonStep();
  void step(ColumnVector&x,ColumnVector&grad,SymmetricMatrix&hess);
};

class GeneralizedNewtonStep: public NewtonStep
{
#   define CLASSNAME GeneralizedNewtonStep
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>    
 private:
  double tol;
  int expected_dim;
 public:
  GeneralizedNewtonStep(double tol);
  GeneralizedNewtonStep(int tol);
  GeneralizedNewtonStep(StateIn&);
  GeneralizedNewtonStep(KeyVal&);
  void save_data_state(StateOut&);
  void step(ColumnVector&x,ColumnVector&grad,SymmetricMatrix&hess);
};

//----------------------------------------------------------------------
// Newton-Like Methods
//----------------------------------------------------------------------

class OptNewtonLike: public Optimize {
#   define CLASSNAME OptNewtonLike
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 protected:
  ColumnVector gprev;
  SymmetricMatrix Hessian;
  int grad_evals;
  RefNewtonStep stepper;
 public:
  OptNewtonLike();
  OptNewtonLike(int n);
  OptNewtonLike(int n, TOLS* t);
  OptNewtonLike(NewtonStep*s);
  OptNewtonLike(int n,NewtonStep*s);
  OptNewtonLike(int n, TOLS* t,NewtonStep*s);
  virtual ~OptNewtonLike();
    OptNewtonLike(StateIn&);
    OptNewtonLike(KeyVal&);
    void save_data_state(StateOut&);

  virtual int CheckDeriv();
  virtual int CheckConvg();

//  virtual int GetFevals() = 0;
//  virtual int GetGevals() = 0;
  void PrintStatus();
};

//----------------------------------------------------------------------
// Full Newton Method
//----------------------------------------------------------------------

class OptNewton: public OptNewtonLike {
#   define CLASSNAME OptNewton
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 protected:
  RefNLP2 nlp;
 public:
  ~OptNewton();
  OptNewton();
  OptNewton(NLP2* p);
  OptNewton(NLP2* p, TOLS* t);
  OptNewton(NewtonStep*s);
  OptNewton(NLP2* p,NewtonStep*s);
  OptNewton(NLP2* p, TOLS* t,NewtonStep*s);
    OptNewton(StateIn&);
    OptNewton(KeyVal&);
    void save_data_state(StateOut&);

  void DefProblem(NLP2 *p);
  NLP2* GetProblem();
//  int GetFevals() {fcn_evals = nlp->GetFevals(); return fcn_evals;}
//  int GetGevals() {grad_evals = nlp->GetGevals();return grad_evals;}

// These are defined elsewhere

  void optimize();
  int CheckDeriv();
  int CheckConvg();
};

//----------------------------------------------------------------------
// Hessian Update Methods for Newton Methods like Quasi and Finite Diff
//----------------------------------------------------------------------

class HessianUpdate: virtual public SavableState {
#   define CLASSNAME HessianUpdate
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    SymmetricMatrix*guessH;
  public:
    HessianUpdate();
    HessianUpdate(StateIn&);
    HessianUpdate(KeyVal&);
    void save_data_state(StateOut&);
    HessianUpdate(SymmetricMatrix&guess);
    virtual ~HessianUpdate();
    // this routine should modify H
    virtual void Initial(SymmetricMatrix&H,NLP1&nlp);
    virtual void Update(SymmetricMatrix&H,NLP1&nlp) = 0;
};
DescribedClass_REF_dec(HessianUpdate);

class BFGSupdate: public HessianUpdate {
#   define CLASSNAME BFGSupdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    ColumnVector xprev;
    ColumnVector gprev;
  public:
    BFGSupdate();
    BFGSupdate(SymmetricMatrix&guess);
    BFGSupdate(StateIn&);
    BFGSupdate(KeyVal&);
    void save_data_state(StateOut&);
    ~BFGSupdate();
    void Initial(SymmetricMatrix&H,NLP1&nlp);
    void Update(SymmetricMatrix&H,NLP1&nlp);
};

class NoUpdate: public HessianUpdate {
#   define CLASSNAME NoUpdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    NoUpdate();
    NoUpdate(SymmetricMatrix&guess);
    NoUpdate(StateIn&);
    NoUpdate(KeyVal&);
    void save_data_state(StateOut&);
    ~NoUpdate();
    void Update(SymmetricMatrix&H,NLP1&nlp);
};

class FDupdate: public HessianUpdate {
// #   define CLASSNAME FDupdate
// #   define HAVE_CTOR
// #   define HAVE_KEYVAL_CTOR
// #   define HAVE_STATEIN_CTOR
// #   include <util/state/stated.h>
// #   include <util/class/classd.h>
  public:
    FDupdate();
    FDupdate(StateIn&);
    FDupdate(KeyVal&);
    FDupdate(SymmetricMatrix&guess);
    ~FDupdate();
    void Initial(SymmetricMatrix&H,NLP1&nlp);
    void Update(SymmetricMatrix&H,NLP1&nlp);
};

//----------------------------------------------------------------------
// Quasi-Newton Method
//----------------------------------------------------------------------

class OptQNewton: public OptNewtonLike {
#   define CLASSNAME OptQNewton
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 private:
  void init(RefHessianUpdate);
 protected:
  RefNLP1 nlp;
  RefHessianUpdate update;
 public:
  OptQNewton();
  OptQNewton(NLP1* p);
  OptQNewton(NLP1* p, TOLS* t);
  OptQNewton(HessianUpdate*,NewtonStep*s);
  OptQNewton(NLP1* p,HessianUpdate*,NewtonStep*s);
  OptQNewton(NLP1* p, TOLS* t,HessianUpdate*,NewtonStep*s);
  ~OptQNewton();
    OptQNewton(StateIn&);
    OptQNewton(KeyVal&);
    void save_data_state(StateOut&);

  void DefProblem(NLP1 *p);
  NLP1* GetProblem();

// These are defined elsewhere

  void optimize();
  SymmetricMatrix UpdateH(SymmetricMatrix& H, int k);
  int CheckDeriv();
  int CheckConvg();
};

//----------------------------------------------------------------------
// Finite-Difference Newton Method
//----------------------------------------------------------------------
// 
// class OptFDNewton: public OptNewtonLike {
//  protected:
//   NLP1 *nlp;
//   HessianUpdate*update;
//  public:
// 
//   OptFDNewton();
//   OptFDNewton(NLP1* p);
//   OptFDNewton(NLP1* p, TOLS* t);
//   OptFDNewton(HessianUpdate*,NewtonStep*s);
//   OptFDNewton(NLP1* p,HessianUpdate*,NewtonStep*s);
//   OptFDNewton(NLP1* p, TOLS* t,HessianUpate*,NewtonStep*s);
//   ~OptFDNewton();
// 
//   void DefProblem(NLP1 *p);
//   NLP1* GetProblem() { return nlp; }
// 
// // These are defined elsewhere
// 
//   void optimize();
//   int CheckDeriv();
//   int CheckConvg();
// };

  
//-------------------------------------------------------------------------
// Various optimization methods and related support routines
//-------------------------------------------------------------------------

int mcstep(double&, double&, double&, double&, double&, 
	   double&, double&, double&, double&, long int&, 
	   double&, double&, int&);

int mcsrch(NLP1*, ColumnVector, double&, int&, int,
	   double ftol = 1.e-4, double xtol = 2.2e-16, 
	   double gtol = 0.9, double stpmax = 1.e3, 
	   double stpmin = 1.e-9);

ReturnMatrix PertChol(SymmetricMatrix&, Real, Real&);
ReturnMatrix MCholesky(SymmetricMatrix&);

#endif
