
#include "nlp.h"
#include "opt.h"
#include <math/nihmatrix/nihmatrix.h>
#include <math/nihmatrix/newmat.h>
#include <util/keyval/keyval.h>

DescribedClass_REF_def(Optimize);

#define CLASSNAME Optimize
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Optimize::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}
void
Optimize::save_data_state(StateOut&s)
{
  s.put(dim);
  tol->save_state(s);
  DVector d;
  Convert(sx,d);
  d.save_object_state(s);
  Convert(sfx,d);
  d.save_object_state(s);
  Convert(xprev,d);
  d.save_object_state(s);
  s.put(fprev);
  s.put(step_length);
  s.put_array_char(method,80);
  s.put_array_char(mesg,80);
  s.put(ret_code);
  s.put(iter_taken);
}
Optimize::Optimize(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(dim);
  tol = TOLS::restore_state(s);
  DVector d(s);
  Convert(d,sx);
  DVector d2(s);
  Convert(d2,sfx);
  DVector d3(s);
  Convert(d3,xprev);
  s.get(fprev);
  s.get(step_length);
  s.get_array_char(method,80);
  s.get_array_char(mesg,80);
  s.get(ret_code);
  s.get(iter_taken);
  outfile = 0;
}
Optimize::Optimize(KeyVal&kv):
  dim(0),outfile(0)
{
  tol = kv.describedclassvalue("tols");
}

  Optimize::Optimize():dim(0),outfile(0),tol(0) {}
  Optimize::Optimize(int n):dim(n),outfile(0),tol(0),sx(n),sfx(n),xprev(n),fcn_evals(0)
    {sx  = 1.0; sfx = 1.0; xprev = 0.0;}
  Optimize::Optimize(TOLS* t): dim(0), outfile(0), tol(t) {}
  Optimize::Optimize(int n, TOLS* t):dim(n),outfile(0),sx(n),sfx(n),
  xprev(n),fcn_evals(0),tol(t) {sx  = 1.0; sfx = 1.0; xprev = 0.0;}
  Optimize::~Optimize() {}
  void Optimize::SetTols(TOLS *t)    { tol =  t; }
  void Optimize::SetOutput(FILE *fp) { outfile = fp; }
  void Optimize::SetTol(TOLS *t) { tol =  t; }
  void Optimize::SetMesg(const char *s)   {strcpy(mesg,s);}     // Set message 
  void Optimize::SetMethod(const char *s) {strcpy(method,s);}   // Set method of choice

  ColumnVector Optimize::Getsx()  {return sx;}
  ColumnVector Optimize::Getsfx() {return sfx;}

  OptDirect::OptDirect(){}
  OptDirect::~OptDirect(){}

  OptPDS::OptPDS(){}
  OptPDS::~OptPDS(){}

#define CLASSNAME OptNewton
#define PARENTS public OptNewtonLike
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
OptNewton::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { OptNewtonLike::_castdown(cd) };
  return do_castdowns(casts,cd);
}
OptNewton::OptNewton(KeyVal&kv):
  OptNewtonLike(kv)
{
  nlp = kv.describedclassvalue("nlp");
  dim = nlp->GetDim();
}
OptNewton::OptNewton(StateIn&s):
  SavableState(s,class_desc_),
  OptNewtonLike(s)
{
  nlp = NLP2::restore_state(s);
}
void
OptNewton::save_data_state(StateOut&s)
{
  OptNewtonLike::save_data_state(s);
  nlp->save_state(s);
}

#define CLASSNAME OptQNewton
#define PARENTS public OptNewtonLike
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
OptQNewton::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { OptNewtonLike::_castdown(cd) };
  return do_castdowns(casts,cd);
}
OptQNewton::OptQNewton(KeyVal&kv):
  OptNewtonLike(kv)
{
  nlp = kv.describedclassvalue("nlp");
  RefHessianUpdate u = kv.describedclassvalue("update");
  dim = nlp->GetDim();
  init(u);
}
OptQNewton::OptQNewton(StateIn&s):
  SavableState(s,class_desc_),
  OptNewtonLike(s)
{
  nlp = NLP1::restore_state(s);
  RefHessianUpdate u = HessianUpdate::restore_state(s);
  init(u);
}
void
OptQNewton::save_data_state(StateOut&s)
{
  OptNewtonLike::save_data_state(s);
  nlp->save_state(s);
  update->save_state(s);
}

#define CLASSNAME OptDirect
#define PARENTS public Optimize
// #define HAVE_CTOR
// #define HAVE_KEYVAL_CTOR
// #define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
OptDirect::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { Optimize::_castdown(cd) };
  return do_castdowns(casts,cd);
}

#define CLASSNAME OptPDS
#define PARENTS public OptDirect
// #define HAVE_CTOR
// #define HAVE_KEYVAL_CTOR
// #define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
OptPDS::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { OptDirect::_castdown(cd) };
  return do_castdowns(casts,cd);
}

#define CLASSNAME OptCGLike
#define PARENTS public Optimize
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
OptCGLike::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { Optimize::_castdown(cd) };
  return do_castdowns(casts,cd);
}
OptCGLike::OptCGLike(StateIn&s):
  SavableState(s,class_desc_),
  Optimize(s)
{
  DVector d(s);
  Convert(d,gprev);
  grad_evals = 0;
}
void
OptCGLike::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  DVector d;
  Convert(gprev,d);
  d.save_object_state(s);
}
OptCGLike::OptCGLike(KeyVal&kv):
  Optimize(kv)
{
}
  OptCGLike::OptCGLike(){}
  OptCGLike::OptCGLike(int n): Optimize(n),gprev(n),grad_evals(0){}
  OptCGLike::OptCGLike(int n, TOLS* t): Optimize(n,t),gprev(n),grad_evals(0){}
  OptCGLike::~OptCGLike(){}
  int OptCGLike::CheckConvg() {return 0;}

  OptCG::OptCG(){strcpy(method,"Nonlinear CG");}
  OptCG::~OptCG(){}
  void OptCG::DefProblem(NLP1 *p) {nlp =  p; dim = nlp->GetDim();}
  NLP1* OptCG::GetProblem() { return nlp; }

  int OptNewtonLike::CheckConvg() {return 0;}

  void OptNewton::DefProblem(NLP2 *p) {nlp =  p; dim = nlp->GetDim();}
  NLP2* OptNewton::GetProblem() { return nlp; }

  NLP1* OptQNewton::GetProblem() { return nlp; }

OptNewton::~OptNewton(){}
