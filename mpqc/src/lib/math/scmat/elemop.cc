
#ifdef __GNUC__
#pragma implementation
#endif

#include <iostream.h>
#include <stdlib.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/elemop.h>

#ifdef PARAGON
extern "C" {
    double drand48();
}
#endif

/////////////////////////////////////////////////////////////////////////////
// SCElementOp member functions

SavableState_REF_def(SCElementOp);

#define CLASSNAME SCElementOp
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCElementOp::SCElementOp()
{
}

void *
SCElementOp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCElementOp::~SCElementOp()
{
}

int
SCElementOp::has_collect()
{
  return 0;
}

int
SCElementOp::has_side_effects()
{
  return 0;
}

void
SCElementOp::collect(RefSCElementOp&)
{
}

void
SCElementOp::process(SCMatrixBlock* a)
{
  a->process(this);
}

// If specializations of SCElementOp do not handle a particle
// block type, then these functions will be called and will
// set up an appropiate block iterator which specializations
// of SCElementOp must handle since it is pure virtual.

void
SCElementOp::process(SCMatrixRectBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixRectBlockIter(a);
  SCMatrixBlockIter&r=*i;
  process(r);
  // this causes a SCMatrixRectBlock::operator int() to be
  // called with this = 0x0 using gcc 2.5.6
  // process(*i,b);
  delete i;
}
void
SCElementOp::process(SCMatrixLTriBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixLTriBlockIter(a);
  process(*i);
  delete i;
}
void
SCElementOp::process(SCMatrixDiagBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixDiagBlockIter(a);
  process(*i);
  delete i;
}
void
SCElementOp::process(SCVectorSimpleBlock* a)
{
  SCMatrixBlockIter*i = new SCVectorSimpleBlockIter(a);
  process(*i);
  delete i;
}

/////////////////////////////////////////////////////////////////////////////
// SCElementOp2 member functions

SavableState_REF_def(SCElementOp2);

#define CLASSNAME SCElementOp2
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCElementOp2::SCElementOp2()
{
}

void *
SCElementOp2::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCElementOp2::~SCElementOp2()
{
}

int
SCElementOp2::has_collect()
{
  return 0;
}

int
SCElementOp2::has_side_effects()
{
  return 0;
}

int
SCElementOp2::has_side_effects_in_arg()
{
  return 0;
}

void
SCElementOp2::collect(RefSCElementOp2&)
{
}

void
SCElementOp2::process(SCMatrixBlock* a, SCMatrixBlock* b)
{
  a->process(this, b);
}

// If specializations of SCElementOp2 do not handle a particle
// block type, then these functions will be called and will
// set up an appropiate block iterator which specializations
// of SCElementOp2 must handle since it is pure virtual.

void
SCElementOp2::process(SCMatrixRectBlock* a,SCMatrixRectBlock* b)
{
  SCMatrixBlockIter*i = new SCMatrixRectBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixRectBlockIter(b);
  process(*i,*j);
  // this causes a SCMatrixRectBlock::operator int() to be
  // called with this = 0x0 using gcc 2.5.6
  // process(*i,b);
  delete i;
  delete j;
}
void
SCElementOp2::process(SCMatrixLTriBlock* a,SCMatrixLTriBlock* b)
{
  SCMatrixBlockIter*i = new SCMatrixLTriBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixLTriBlockIter(b);
  process(*i,*j);
  delete i;
  delete j;
}
void
SCElementOp2::process(SCMatrixDiagBlock* a,SCMatrixDiagBlock* b)
{
  SCMatrixBlockIter*i = new SCMatrixDiagBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixDiagBlockIter(b);
  process(*i,*j);
  delete i;
  delete j;
}
void
SCElementOp2::process(SCVectorSimpleBlock* a,SCVectorSimpleBlock* b)
{
  SCMatrixBlockIter*i = new SCVectorSimpleBlockIter(a);
  SCMatrixBlockIter*j = new SCVectorSimpleBlockIter(b);
  process(*i,*j);
  delete i;
  delete j;
}

/////////////////////////////////////////////////////////////////////////////
// SCElementOp3 member functions

SavableState_REF_def(SCElementOp3);

#define CLASSNAME SCElementOp3
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCElementOp3::SCElementOp3()
{
}

void *
SCElementOp3::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCElementOp3::~SCElementOp3()
{
}

int
SCElementOp3::has_collect()
{
  return 0;
}

int
SCElementOp3::has_side_effects()
{
  return 0;
}

int
SCElementOp3::has_side_effects_in_arg1()
{
  return 0;
}

int
SCElementOp3::has_side_effects_in_arg2()
{
  return 0;
}

void
SCElementOp3::collect(RefSCElementOp3&)
{
}

void
SCElementOp3::process(SCMatrixBlock* a,
                      SCMatrixBlock* b,
                      SCMatrixBlock* c)
{
  a->process(this, b, c);
}

// If specializations of SCElementOp3 do not handle a particle
// block type, then these functions will be called and will
// set up an appropiate block iterator which specializations
// of SCElementOp3 must handle since it is pure virtual.

void
SCElementOp3::process(SCMatrixRectBlock* a,
                      SCMatrixRectBlock* b,
                      SCMatrixRectBlock* c)
{
  SCMatrixBlockIter*i = new SCMatrixRectBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixRectBlockIter(b);
  SCMatrixBlockIter*k = new SCMatrixRectBlockIter(c);
  process(*i,*j,*k);
  delete i;
  delete j;
  delete k;
}
void
SCElementOp3::process(SCMatrixLTriBlock* a,
                      SCMatrixLTriBlock* b,
                      SCMatrixLTriBlock* c)
{
  SCMatrixBlockIter*i = new SCMatrixLTriBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixLTriBlockIter(b);
  SCMatrixBlockIter*k = new SCMatrixLTriBlockIter(c);
  process(*i,*j,*k);
  delete i;
  delete j;
  delete k;
}
void
SCElementOp3::process(SCMatrixDiagBlock* a,
                      SCMatrixDiagBlock* b,
                      SCMatrixDiagBlock* c)
{
  SCMatrixBlockIter*i = new SCMatrixDiagBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixDiagBlockIter(b);
  SCMatrixBlockIter*k = new SCMatrixDiagBlockIter(c);
  process(*i,*j,*k);
  delete i;
  delete j;
  delete k;
}
void
SCElementOp3::process(SCVectorSimpleBlock* a,
                      SCVectorSimpleBlock* b,
                      SCVectorSimpleBlock* c)
{
  SCMatrixBlockIter*i = new SCVectorSimpleBlockIter(a);
  SCMatrixBlockIter*j = new SCVectorSimpleBlockIter(b);
  SCMatrixBlockIter*k = new SCVectorSimpleBlockIter(c);
  process(*i,*j,*k);
  delete i;
  delete j;
  delete k;
}

/////////////////////////////////////////////////////////////////////////
// SCElementScale members

#define CLASSNAME SCElementScale
#define PARENTS   public SCElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementScale::SCElementScale(double a):scale(a) {}
SCElementScale::SCElementScale(StateIn&s):
  SCElementOp(s)
{
  s.get(scale);
}
void
SCElementScale::save_data_state(StateOut&s)
{
  s.put(scale);
}
void *
SCElementScale::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}
SCElementScale::~SCElementScale() {}
void
SCElementScale::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.set(scale*i.get());
    }
}

int
SCElementScale::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementScalarProduct members

#define CLASSNAME SCElementScalarProduct
#define PARENTS   public SCElementOp2
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SCElementScalarProduct::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp2::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCElementScalarProduct::SCElementScalarProduct():
  product(0.0)
{
}

SCElementScalarProduct::SCElementScalarProduct(StateIn&s):
  SCElementOp2(s)
{
  s.get(product);
}

void
SCElementScalarProduct::save_data_state(StateOut&s)
{
  s.put(product);
}

SCElementScalarProduct::~SCElementScalarProduct()
{
}

void
SCElementScalarProduct::process(SCMatrixBlockIter&i,
                                SCMatrixBlockIter&j)
{
  for (i.reset(),j.reset(); i; ++i,++j) {
      product += i.get()*j.get();
    }
}

int
SCElementScalarProduct::has_collect()
{
  return 1;
}

void
SCElementScalarProduct::collect(RefSCElementOp&op)
{
  RefSCElementScalarProduct ma(op);
  product += ma->product;
}

double
SCElementScalarProduct::result()
{
  return product;
}

/////////////////////////////////////////////////////////////////////////
// SCDestructiveElementProduct members

#define CLASSNAME SCDestructiveElementProduct
#define PARENTS   public SCElementOp2
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCDestructiveElementProduct::SCDestructiveElementProduct() {}
SCDestructiveElementProduct::SCDestructiveElementProduct(StateIn&s):
  SCElementOp2(s)
{
}
void
SCDestructiveElementProduct::save_data_state(StateOut&s)
{
}
void *
SCDestructiveElementProduct::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp2::_castdown(cd);
  return do_castdowns(casts,cd);
}
SCDestructiveElementProduct::~SCDestructiveElementProduct() {}
void
SCDestructiveElementProduct::process(SCMatrixBlockIter&i,
                                     SCMatrixBlockIter&j)
{
  for (i.reset(),j.reset(); i; ++i,++j) {
      i.set(i.get()*j.get());
    }
}

int
SCDestructiveElementProduct::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementInvert members

#define CLASSNAME SCElementInvert
#define PARENTS   public SCElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementInvert::SCElementInvert() {}
SCElementInvert::SCElementInvert(double a) {}
SCElementInvert::SCElementInvert(StateIn&s):
  SCElementOp(s)
{
}
void
SCElementInvert::save_data_state(StateOut&s)
{
}
void *
SCElementInvert::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}
SCElementInvert::~SCElementInvert() {}
void
SCElementInvert::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.set(1.0/i.get());
    }
}

int
SCElementInvert::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementSquareRoot members

#define CLASSNAME SCElementSquareRoot
#define PARENTS   public SCElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementSquareRoot::SCElementSquareRoot() {}
SCElementSquareRoot::SCElementSquareRoot(double a) {}
SCElementSquareRoot::SCElementSquareRoot(StateIn&s):
  SCElementOp(s)
{
}
void
SCElementSquareRoot::save_data_state(StateOut&s)
{
}
void *
SCElementSquareRoot::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}
SCElementSquareRoot::~SCElementSquareRoot() {}
void
SCElementSquareRoot::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.set(sqrt(i.get()));
    }
}

int
SCElementSquareRoot::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementMaxAbs members

SavableState_REF_def(SCElementMaxAbs);
#define CLASSNAME SCElementMaxAbs
#define PARENTS   public SCElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCElementMaxAbs::SCElementMaxAbs():r(0.0) {}
SCElementMaxAbs::SCElementMaxAbs(StateIn&s):
  SCElementOp(s)
{
  s.get(r);
}
void
SCElementMaxAbs::save_data_state(StateOut&s)
{
  s.put(r);
}
void *
SCElementMaxAbs::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}
SCElementMaxAbs::~SCElementMaxAbs() {}
void
SCElementMaxAbs::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      if (fabs(i.get()) > r) r = fabs(i.get());
    }
}
double
SCElementMaxAbs::result()
{
  return r;
}
int
SCElementMaxAbs::has_collect()
{
  return 1;
}
void
SCElementMaxAbs::collect(RefSCElementOp&op)
{
  RefSCElementMaxAbs ma(op);
  if (ma->r > r) r = ma->r;
}

/////////////////////////////////////////////////////////////////////////
// SCElementAssign members

#define CLASSNAME SCElementAssign
#define PARENTS   public SCElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementAssign::SCElementAssign(double a):assign(a) {}
SCElementAssign::SCElementAssign(StateIn&s):
  SCElementOp(s)
{
  s.get(assign);
}
void
SCElementAssign::save_data_state(StateOut&s)
{
  s.put(assign);
}
void *
SCElementAssign::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}
SCElementAssign::~SCElementAssign() {}
void
SCElementAssign::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.set(assign);
    }
}

int
SCElementAssign::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementRandomize members

#define CLASSNAME SCElementRandomize
#define PARENTS   public SCElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementRandomize::SCElementRandomize() {}
SCElementRandomize::SCElementRandomize(StateIn&s):
  SCElementOp(s)
{
}
void
SCElementRandomize::save_data_state(StateOut&s)
{
}
void *
SCElementRandomize::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}
SCElementRandomize::~SCElementRandomize() {}
void
SCElementRandomize::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.set(drand48()*(drand48()<0.5?1.0:-1.0));
    }
}

int
SCElementRandomize::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementShiftDiagonal members

#define CLASSNAME SCElementShiftDiagonal
#define PARENTS   public SCElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementShiftDiagonal::SCElementShiftDiagonal(double a):shift_diagonal(a) {}
SCElementShiftDiagonal::SCElementShiftDiagonal(StateIn&s):
  SCElementOp(s)
{
  s.get(shift_diagonal);
}
void
SCElementShiftDiagonal::save_data_state(StateOut&s)
{
  s.put(shift_diagonal);
}
void *
SCElementShiftDiagonal::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}
SCElementShiftDiagonal::~SCElementShiftDiagonal() {}
void
SCElementShiftDiagonal::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      if (i.i() == i.j()) i.set(shift_diagonal+i.get());
    }
}

int
SCElementShiftDiagonal::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementScaleDiagonal members

#define CLASSNAME SCElementScaleDiagonal
#define PARENTS   public SCElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementScaleDiagonal::SCElementScaleDiagonal(double a):scale_diagonal(a) {}
SCElementScaleDiagonal::SCElementScaleDiagonal(StateIn&s):
  SCElementOp(s)
{
  s.get(scale_diagonal);
}
void
SCElementScaleDiagonal::save_data_state(StateOut&s)
{
  s.put(scale_diagonal);
}
void *
SCElementScaleDiagonal::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}
SCElementScaleDiagonal::~SCElementScaleDiagonal() {}
void
SCElementScaleDiagonal::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      if (i.i() == i.j()) i.set(scale_diagonal*i.get());
    }
}

int
SCElementScaleDiagonal::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementDot members

#define CLASSNAME SCElementDot
#define PARENTS   public SCElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SCElementDot::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCElementDot::SCElementDot(double**a, double**b, int n):
  avects(a),
  bvects(b),
  length(n)
{
}

SCElementDot::SCElementDot(StateIn&s)
{
  fprintf(stderr,"SCElementDot does not permit StateIn CTOR\n");
  abort();
}

void
SCElementDot::save_data_state(StateOut&s)
{
  fprintf(stderr,"SCElementDot does not permit save_data_state\n");
  abort();
}

int
SCElementDot::has_side_effects()
{
  return 1;
}

void
SCElementDot::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      double tmp = i.get();
      double* a = avects[i.i()];
      double* b = bvects[i.j()];
      for (int j = length; j; j--, a++, b++) {
          tmp += *a * *b;
        }
      i.set(tmp);
    }
}
