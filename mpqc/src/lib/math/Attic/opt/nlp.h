#ifndef _math_opt_nlp_h
#define _math_opt_nlp_h
//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include <util/state/state.h>
#include <util/misc/compute.h>
#include <math/newmat7/newmatap.h>
#include <math/opt/tstfcn.h>

//  Some useful typdefs

typedef void (*USERFCN0)(int, ColumnVector, double&);

typedef void (*USERFCN1)(int, int, ColumnVector, double&, ColumnVector&);

typedef void (*USERFCN2)(int, int, ColumnVector, double&, 
			ColumnVector&, SymmetricMatrix&);

class Optimize;
//------------------------------------------------------------------------
//
// Base Class for NonLinear Programming Problem
// For NLP0 the only assumption on the objective function is
// that it be continuous
// Note that NLP0-2 are abstract data types
//
//------------------------------------------------------------------------

Result_dec(ColumnVector);
Result_dec(SymmetricMatrix);

class NLP0: virtual public SavableState, public Compute {
#   define CLASSNAME NLP0
#   include <util/state/stated.h>
#   include <util/class/classd.h>
protected:
  int          dim;          // Dimension of the problem
  int          deriv_level;  // Number of derivatives available
  ColumnVector xc;           // Current point
  Resultdouble       fvalue; // Function value of objective function at xc
  double  fcn_accrcy;        // Accuracy with which the function can be evaluated
  int     nfevals;
public:
  NLP0();
  virtual ~NLP0();
  NLP0(int ndim);
    double value();

    NLP0(StateIn&);
    NLP0(KeyVal&);
    void save_data_state(StateOut&);

  virtual void SetDim(int dim);

  virtual double EvalF();
  virtual double EvalF(ColumnVector x);

// Evaluate everything (just f in this case)
  virtual void Eval();

// Functions for setting various properties of this NLP problem
// and for evaluating the objective function

  virtual void SetX(int i, double x);
  virtual void SetX(ColumnVector& x);

  void SetDerivLevel(int k);
  int GetDim();
  int GetFevals();
  ColumnVector GetXc();
  const ColumnVector& GetX() const;
  double GetF();

// These are defined elsewhere

  virtual ColumnVector FDGrad(ColumnVector);        // Finite-difference gradient
  virtual SymmetricMatrix FD2Hessian(ColumnVector); // Finite-difference Hessian

  virtual void PrintState(char *); // Print state

// Friends
  
  friend  Optimize;
};
DescribedClass_REF_dec(NLP0);

//------------------------------------------------------------------------
// Derived Classes from NonLinear Programming Problem
// NLP0: NLP  + No derivative info
// NLP1: NLP0 + First derivatives
// NLP2: NLP1 + Second derivatives
//------------------------------------------------------------------------

class NLP1: public NLP0 {
#   define CLASSNAME NLP1
#   include <util/state/stated.h>
#   include <util/class/classd.h>
protected:
  ResultColumnVector grad;     // Gradient of objective function at xc
  int     ngevals;

public:
  NLP1();
  NLP1(int ndim);
  virtual ~NLP1();
    const ColumnVector& gradient();

    NLP1(StateIn&);
    NLP1(KeyVal&);
    void save_data_state(StateOut&);

  void SetDim(int dim);

// Make this an abstract data type
  virtual ColumnVector EvalG();
  virtual ColumnVector EvalG(ColumnVector x);

  virtual void Eval();

// Evaluate grad at x

  ColumnVector GetGrad();
  void    PrintState(char *);    // Print state
  int GetGevals();
  virtual SymmetricMatrix FDHessian(ColumnVector);  // Finite-difference Hessian

};
DescribedClass_REF_dec(NLP1);

class NLP2: public NLP1 {
#   define CLASSNAME NLP2
#   include <util/state/stated.h>
#   include <util/class/classd.h>
protected:
  // Hessian of objective function (or approximation)
  ResultSymmetricMatrix Hessian;
  int     nhevals;

public:
  NLP2();
  NLP2(int ndim);
  virtual ~NLP2();

    NLP2(StateIn&);
    NLP2(KeyVal&);
    void save_data_state(StateOut&);

  void SetDim(int dim);
  
// Make this an abstract data type
  virtual SymmetricMatrix EvalH();
  virtual void Eval();

// Functions for evaluating objective function and derivatives

    const SymmetricMatrix& hessian();
  SymmetricMatrix GetHess();
  void    PrintState(char *);    // Print state
  int GetHevals();

};
DescribedClass_REF_dec(NLP2);

//
//
//  Derived from NLP's
//
class NLF0: public NLP0 {
protected:
  USERFCN0 fcn;
  
 public:
  NLF0();
  virtual ~NLF0();
  NLF0(int ndim);
  NLF0(int ndim, USERFCN0 f);

    void compute();
  double EvalF(ColumnVector x); // Evaluate f and return the value
  void SetF(USERFCN0 f);
};

class NLF1: public NLP1 {
protected:
  USERFCN1 fcn;
  
public:
  NLF1();
  NLF1(int ndim);
  NLF1(int ndim, USERFCN1 f);
  virtual ~NLF1();

    void compute();
  double EvalF(ColumnVector x);         // Evaluate f and return the value
  ColumnVector EvalG(ColumnVector x);   // Evaluate grad
  void SetF(USERFCN1 f);

};

class NLF2: public NLP2 {
protected:
  USERFCN2 fcn;

public:
  NLF2();
  NLF2(int ndim);
  NLF2(int ndim, USERFCN2 f);
  virtual ~NLF2();

    void compute();
  double EvalF(ColumnVector x);         // Evaluate f and return the value
  ColumnVector EvalG(ColumnVector x);   // Evaluate grad
  void SetF(USERFCN2 f);
};

//------------------------------------------------------------------------
//
// Class for Monitoring Optimization Progress
//
//------------------------------------------------------------------------

class OptProgress {
 private:
  double _asnorm;
  double _adeltaf;
  double _agnorm;
  double _fvalue;
  double _xnorm;
  const char* _status;
 public:
  OptProgress(double f,double oldf,
              const ColumnVector&x,const ColumnVector&oldx,
              const ColumnVector&g);
  double XNorm();
  double ASNorm(); // absolute SNorm
  double RSNorm(); // relative SNorm (to xnorm)
  double ADeltaF(); // absolute Delta F
  double RDeltaF(); // relative Delta F (to function value)
  double AGNorm(); // absolute G norm
  double RGNorm(); // relative G norm (to function value)
  void ResetStatus();
  void SetStatus(const char*);
  const char* GetStatus();
};

//------------------------------------------------------------------------
//
// Base Class for Tolerances
//
//------------------------------------------------------------------------

class TOLS: virtual public SavableState
{
#   define CLASSNAME TOLS
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
  double  mcheps;         // Machine epsilon
  double  max_step;       // Maximum step allowed
  double  step_tol;       // Step tolerance used for convergence
  double  fcn_tol;        // Function tolerance used for convergence
  double  grad_tol;       // Gradient tolerance used for convergence
  int     max_iter;       // Maximum number of iterations allowed
  int     max_feval;      // Maximum number of function evaluations allowed

  int use_ftol;
  int use_gtol;
  int use_agtol;
  int use_steptol;

  int or_tol;
public:
    TOLS();
    TOLS(StateIn&);
    TOLS(KeyVal&);

    void save_data_state(StateOut&);

  void    SetDefaultTol();
  void    SetStepTol(double x);
  void    SetFtol(double x);
  void    SetGtol(double x);
  void    SetMaxIter(int k);
  void    SetMaxFeval(int k);

  void    SetUseFtol(int);
  void    SetUseGtol(int);
  void    SetUseAGtol(int);
  void    SetUseStepTol(int);

  int     UseFtol();
  int     UseGtol();
  int     UseAGtol();
  int     UseStepTol();

  void    SetOrTol(int);
  void    SetAndTol(int);
  int     OrTol();
  int     AndTol();

  double  GetStepTol();
  double  GetFtol();
  double  GetGtol();
  int     GetMaxIter();
  int     GetMaxFeval();

  void    PrintTol();

  int Converged(OptProgress&);

  friend  Optimize;
};
DescribedClass_REF_dec(TOLS);

//------------------------------------------------------------------------
// Print various quantities in newmat classes
//------------------------------------------------------------------------

// void Print(const Matrix& X);
// void Print(const UpperTriangularMatrix& X);
// void Print(const DiagonalMatrix& X);
// void Print(const SymmetricMatrix& X);
// void Print(const LowerTriangularMatrix& X);

#endif
