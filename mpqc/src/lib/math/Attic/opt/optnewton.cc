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
#include <math/nihmatrix/nihmatrix.h>
#include <math/nihmatrix/newmat.h>

//------------------------------------------------------------------------
//
//   NewtonStep and derivatives' member functions
//
//------------------------------------------------------------------------

DescribedClass_REF_def(NewtonStep);

#define CLASSNAME NewtonStep
#define PARENTS public virtual SavableState
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
NewtonStep::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

#define CLASSNAME MCholeskyNewtonStep
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public NewtonStep
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
MCholeskyNewtonStep::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { NewtonStep::_castdown(cd) };
  return do_castdowns(casts,cd);
}

#define CLASSNAME GeneralizedNewtonStep
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public NewtonStep
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
GeneralizedNewtonStep::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { NewtonStep::_castdown(cd) };
  return do_castdowns(casts,cd);
}


NewtonStep::NewtonStep()
{
}

NewtonStep::NewtonStep(StateIn&s):
  SavableState(s,class_desc_)
{
}

NewtonStep::NewtonStep(KeyVal&s)
{
}

NewtonStep::~NewtonStep()
{
}

void
NewtonStep::save_data_state(StateOut&s)
{
}

MCholeskyNewtonStep::MCholeskyNewtonStep()
{
}

MCholeskyNewtonStep::MCholeskyNewtonStep(StateIn&s):
  SavableState(s,class_desc_)
{
}

MCholeskyNewtonStep::MCholeskyNewtonStep(KeyVal&s):
  NewtonStep(s)
{
}

MCholeskyNewtonStep::~MCholeskyNewtonStep()
{
}

void
MCholeskyNewtonStep::save_data_state(StateOut&s)
{
}

void MCholeskyNewtonStep::step(ColumnVector&x,
			       ColumnVector&grad,
			       SymmetricMatrix&hess)
{
  LowerTriangularMatrix L(hess.Nrows());
  L = MCholesky(hess);
  x = -(L.t().i()*(L.i()*grad));
}

  /* Do a quick sort of the double data in item, by the integer   *
   * indices in index.  (data in item remains unchanged)          */
static void
iqs(double*item,int*index,int left,int right)
{
  register int i,j;
  double x;
  int y;
 
  i=left; j=right;
  x=item[index[(left+right)/2]];
 
  do {
    while(item[index[i]]<x && i<right) i++;
    while(x<item[index[j]] && j>left) j--;
 
    if (i<=j) {
      if (item[index[i]] != item[index[j]]) {
        y=index[i];
        index[i]=index[j];
        index[j]=y;
        }
      i++; j--;
      }
    } while(i<=j);
       
  if (left<j) iqs (item,index,left,j);
  if (i<right) iqs (item,index,i,right);
  }
void
static iquicksort (double*item,int*index,int n)
{
  int i;
  if (n<=0) return;
  for (i=0; i<n; i++) {
    index[i] = i;
    }
  iqs (item,index,0,n-1);
  }

GeneralizedNewtonStep::GeneralizedNewtonStep(KeyVal&kv)
{
  if (kv.exists("tol")) {
      tol = kv.doublevalue("tol");
      expected_dim = 0;
    }
  else {
      expected_dim = kv.intvalue("dim");
    }
}

GeneralizedNewtonStep::GeneralizedNewtonStep(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(expected_dim);
  if (!expected_dim) s.get(tol);
}

void
GeneralizedNewtonStep::save_data_state(StateOut&s)
{
  s.put(expected_dim);
  if (!expected_dim) s.put(tol);
}

GeneralizedNewtonStep::GeneralizedNewtonStep(double t):
expected_dim(0),
tol(t)
{
}

GeneralizedNewtonStep::GeneralizedNewtonStep(int dim):
expected_dim(dim)
{
}

void GeneralizedNewtonStep::step(ColumnVector&x,
				 ColumnVector&grad,
				 SymmetricMatrix&hess)
{
  int i;
  int dim = hess.Nrows();

  Matrix V(dim,dim);
  DiagonalMatrix D(dim);
  EigenValues(hess,D,V);

  // remove eigenvectors based on tolerances
  if (expected_dim==0) {
      for (i=1; i<=dim; i++) {
	  if (fabs(D(i))>tol) D(i) = 1.0/D(i);
	  else D(i) = 0.0;
	}
    }
  // remove eigenvectors based on the expected dimension
  else {
      // sort the absolute values of the eigenvalues
      double* dabs = new double[dim];
      int* index = new int[dim];
      for (i=0; i<dim; i++) {
	  dabs[i] = fabs(D(i+1));
	  index[i] = i;
	}
      iquicksort(dabs,index,dim);
      for (i=0; i<expected_dim; i++) {
	  D(index[i]+1) = 1.0/D(index[i]+1);
	}
      for (; i<dim; i++) {
	  D(index[i]+1) = 0.0;
	}
   }

  // newmat cannot handle a symmetric matrix here
  Matrix hessi = V * D * V.t();

  x = - hessi * grad;

}

//------------------------------------------------------------------------
//
//   Newton-Like Base Class functions
//
//------------------------------------------------------------------------

#define CLASSNAME OptNewtonLike
#define PARENTS public Optimize
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
OptNewtonLike::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { Optimize::_castdown(cd) };
  return do_castdowns(casts,cd);
}
OptNewtonLike::OptNewtonLike(StateIn&s):
  SavableState(s,class_desc_),
  Optimize(s)
{
  DVector d(s);
  Convert(d,gprev);
  DMatrix h(s);
  Convert(h,Hessian);
  stepper = NewtonStep::restore_state(s);
}
void
OptNewtonLike::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  DVector v;
  Convert(gprev,v);
  v.save_object_state(s);
  DMatrix h;
  Convert(Hessian,h);
  h.save_object_state(s);
  stepper->save_state(s);
}
OptNewtonLike::OptNewtonLike(KeyVal&kv):
  Optimize(kv)
{
  stepper = kv.describedclassvalue("step");
}

OptNewtonLike::OptNewtonLike():
  stepper(new MCholeskyNewtonStep)
{
}

OptNewtonLike::OptNewtonLike(int n):
  Optimize(n),
  gprev(n),
  Hessian(n),
  grad_evals(0),
  stepper(new MCholeskyNewtonStep)
{
}

OptNewtonLike::OptNewtonLike(int n, TOLS* t):
  Optimize(n,t),
  gprev(n),
  Hessian(n),
  grad_evals(0),
  stepper(new MCholeskyNewtonStep)
{
}

OptNewtonLike::OptNewtonLike(NewtonStep*s):
stepper(s)
{
}

OptNewtonLike::OptNewtonLike(int n,NewtonStep*s):
  Optimize(n),
  gprev(n),
  Hessian(n),
  grad_evals(0),
  stepper(s)
{
}

OptNewtonLike::OptNewtonLike(int n, TOLS* t,NewtonStep*s):
  Optimize(n,t),
  gprev(n),
  Hessian(n),
  grad_evals(0),
  stepper(s)
{
}

OptNewtonLike::~OptNewtonLike()
{
}

void OptNewtonLike::PrintStatus() // Set Message
{
  printf("%s\n",mesg);

  printf("Dimension of the problem = %d\n",dim);
  printf("Optimization method      = %s\n",method);

  printf("Return code = %d\n",ret_code);
  printf("Number of iterations taken    = %d\n",iter_taken);
//  printf("Number of function evalutions = %d\n",GetFevals());
//  printf("Number of gradient evalutions = %d\n",GetGevals());

  Print(Hessian);
}

