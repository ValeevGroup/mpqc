
#include <math.h>
#include <math/scmat/abstract.h>
#include <math/scmat/blkiter.h>

/////////////////////////////////////////////////////////////////////////
// These member are used by the abstract SCMatrix classes.
/////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////
// SCDimension members

#define CLASSNAME SCDimension
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCDimension::SCDimension()
{
}

void *
SCDimension::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SCDimension::~SCDimension()
{
}

/////////////////////////////////////////////////////////////////////////
// SCElementScale members

#define CLASSNAME SCElementScale
#define PARENTS       virtual public SCDiagElementOp, \
                      virtual public SCSymmElementOp, \
                      virtual public SCRectElementOp, \
                      virtual public SCVectorElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementScale::SCElementScale(double a):scale(a) {}
SCElementScale::SCElementScale(StateIn&s):
  SavableState(s,class_desc_)
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
  void* casts[] =  { SCDiagElementOp::_castdown(cd),
                     SCSymmElementOp::_castdown(cd),
                     SCRectElementOp::_castdown(cd),
                     SCVectorElementOp::_castdown(cd)
                   };
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

/////////////////////////////////////////////////////////////////////////
// SCElementInvert members

#define CLASSNAME SCElementInvert
#define PARENTS       virtual public SCDiagElementOp, \
                      virtual public SCSymmElementOp, \
                      virtual public SCRectElementOp, \
                      virtual public SCVectorElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementInvert::SCElementInvert(double a) {}
SCElementInvert::SCElementInvert(StateIn&s):
  SavableState(s,class_desc_)
{
}
void
SCElementInvert::save_data_state(StateOut&s)
{
}
void *
SCElementInvert::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCDiagElementOp::_castdown(cd),
                     SCSymmElementOp::_castdown(cd),
                     SCRectElementOp::_castdown(cd),
                     SCVectorElementOp::_castdown(cd)
                   };
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

/////////////////////////////////////////////////////////////////////////
// SCElementSquareRoot members

#define CLASSNAME SCElementSquareRoot
#define PARENTS       virtual public SCDiagElementOp, \
                      virtual public SCSymmElementOp, \
                      virtual public SCRectElementOp, \
                      virtual public SCVectorElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementSquareRoot::SCElementSquareRoot(double a) {}
SCElementSquareRoot::SCElementSquareRoot(StateIn&s):
  SavableState(s,class_desc_)
{
}
void
SCElementSquareRoot::save_data_state(StateOut&s)
{
}
void *
SCElementSquareRoot::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCDiagElementOp::_castdown(cd),
                     SCSymmElementOp::_castdown(cd),
                     SCRectElementOp::_castdown(cd),
                     SCVectorElementOp::_castdown(cd)
                   };
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

/////////////////////////////////////////////////////////////////////////
// SCElementMaxAbs members

SavableState_REF_def(SCElementMaxAbs);
#define CLASSNAME SCElementMaxAbs
#define PARENTS       virtual public SCDiagElementOp, \
                      virtual public SCSymmElementOp, \
                      virtual public SCRectElementOp, \
                      virtual public SCVectorElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCElementMaxAbs::SCElementMaxAbs():r(0.0) {}
SCElementMaxAbs::SCElementMaxAbs(StateIn&s):
  SavableState(s,class_desc_)
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
  void* casts[] =  { SCDiagElementOp::_castdown(cd),
                     SCSymmElementOp::_castdown(cd),
                     SCRectElementOp::_castdown(cd),
                     SCVectorElementOp::_castdown(cd)
                   };
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
#define PARENTS       virtual public SCDiagElementOp, \
                      virtual public SCSymmElementOp, \
                      virtual public SCRectElementOp, \
                      virtual public SCVectorElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementAssign::SCElementAssign(double a):assign(a) {}
SCElementAssign::SCElementAssign(StateIn&s):
  SavableState(s,class_desc_)
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
  void* casts[] =  { SCDiagElementOp::_castdown(cd),
                     SCSymmElementOp::_castdown(cd),
                     SCRectElementOp::_castdown(cd),
                     SCVectorElementOp::_castdown(cd)
                   };
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

/////////////////////////////////////////////////////////////////////////
// SCElementShiftDiagonal members

#define CLASSNAME SCElementShiftDiagonal
#define PARENTS       virtual public SCDiagElementOp, \
                      virtual public SCSymmElementOp, \
                      virtual public SCRectElementOp, \
                      virtual public SCVectorElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementShiftDiagonal::SCElementShiftDiagonal(double a):shift_diagonal(a) {}
SCElementShiftDiagonal::SCElementShiftDiagonal(StateIn&s):
  SavableState(s,class_desc_)
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
  void* casts[] =  { SCDiagElementOp::_castdown(cd),
                     SCSymmElementOp::_castdown(cd),
                     SCRectElementOp::_castdown(cd),
                     SCVectorElementOp::_castdown(cd)
                   };
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

/////////////////////////////////////////////////////////////////////////
// SCMatrix members

#define CLASSNAME SCMatrix
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCMatrix::SCMatrix()
{
}

SCMatrix::~SCMatrix()
{
}

SCMatrix::SCMatrix(StateIn&s):
  SavableState(s,class_desc_)
{
}

void
SCMatrix::save_data_state(StateOut&s)
{
}

void *
SCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

int
SCMatrix::nrow()
{
  return rowdim()->n();
}

int
SCMatrix::ncol()
{
  return coldim()->n();
}

double
SCMatrix::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCRectElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
SCMatrix::assign(double a)
{
  RefSCRectElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
SCMatrix::scale(double a)
{
  RefSCRectElementOp op = new SCElementScale(a);
  this->element_op(op);
}

void
SCMatrix::shift_diagonal(double a)
{
  RefSCRectElementOp op = new SCElementShiftDiagonal(a);
  this->element_op(op);
}

void
SCMatrix::unit()
{
  this->assign(0.0);
  this->shift_diagonal(1.0);
}

void
SCMatrix::assign(SCMatrix*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
SCMatrix::assign(const double*a)
{
  int n = nrow();
  const double **v = new double*[n];
  for (int i=0; i<n; i++) {
      v[i] = &a[i*n];
    }
  assign(v);
  delete[] v;
}

void
SCMatrix::assign(const double**a)
{
  int nr = nrow();
  int nc = ncol();
  for (int i=0; i<nr; i++) {
      for (int j=0; j<nc; i++) {
          set_element(i,j,a[i][j]);
        }
    }
}

void
SCMatrix::convert(double*a)
{
  int n = nrow();
  double **v = new double*[n];
  for (int i=0; i<n; i++) {
      v[i] = &a[i*n];
    }
  convert(v);
  delete[] v;
}

void
SCMatrix::convert(double**a)
{
  int nr = nrow();
  int nc = ncol();
  for (int i=0; i<nr; i++) {
      for (int j=0; j<nc; i++) {
          a[i][j] = get_element(i,j);
        }
    }
}

void
SCMatrix::accumulate_product(SymmSCMatrix*a,SCMatrix*b)
{
  SCMatrix *t = b->copy();
  t->transpose_this();
  SCMatrix *t2 = this->copy();
  t2->transpose_this();
  t2->accumulate_product(t,a);
  delete t;
  t2->transpose_this();
  assign(t2);
  delete t2;
}

void
SCMatrix::accumulate_product(DiagSCMatrix*a,SCMatrix*b)
{
  SCMatrix *t = b->copy();
  t->transpose_this();
  SCMatrix *t2 = this->copy();
  t2->accumulate_product(t,a);
  delete t;
  t2->transpose_this();
  assign(t2);
  delete t2;
}

void
SCMatrix::print(ostream&o)
{
  print(0, o, 10);
}

SCMatrix*
SCMatrix::clone()
{
  return rowdim()->create_matrix(coldim());
}

SCMatrix*
SCMatrix::copy()
{
  SCMatrix* result = clone();
  result->assign(this);
  return result;
}

/////////////////////////////////////////////////////////////////////////
// SymmSCMatrix member functions

#define CLASSNAME SymmSCMatrix
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SymmSCMatrix::SymmSCMatrix()
{
}

SymmSCMatrix::SymmSCMatrix(StateIn&s):
  SavableState(s,class_desc_)
{
}

void
SymmSCMatrix::save_data_state(StateOut&s)
{
}

void *
SymmSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

int
SymmSCMatrix::n()
{
  return dim()->n();
}

double
SymmSCMatrix::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCSymmElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
SymmSCMatrix::assign(double a)
{
  RefSCSymmElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
SymmSCMatrix::assign(const double*a)
{
  int nr = n();
  const double **v = new double*[nr];
  int ioff= 0;
  for (int i=0; i<nr; i++) {
      v[i] = &a[ioff];
      ioff += i;
    }
  assign(v);
  delete[] v;
}

void
SymmSCMatrix::assign(const double**a)
{
  int nr = n();
  for (int i=0; i<nr; i++) {
      for (int j=0; j<=i; i++) {
          set_element(i,j,a[i][j]);
        }
    }
}

void
SymmSCMatrix::convert(double*a)
{
  int nr = n();
  double **v = new double*[nr];
  int ioff= 0;
  for (int i=0; i<nr; i++) {
      v[i] = &a[ioff];
      ioff += i;
    }
  convert(v);
  delete[] v;
}

void
SymmSCMatrix::convert(double**a)
{
  int nr = n();
  for (int i=0; i<nr; i++) {
      for (int j=0; j<=i; i++) {
          a[i][j] = get_element(i,j);
        }
    }
}

void
SymmSCMatrix::scale(double a)
{
  RefSCSymmElementOp op = new SCElementScale(a);
  this->element_op(op);
}

void
SymmSCMatrix::shift_diagonal(double a)
{
  RefSCSymmElementOp op = new SCElementShiftDiagonal(a);
  this->element_op(op);
}

void
SymmSCMatrix::unit()
{
  this->assign(0.0);
  this->shift_diagonal(1.0);
}

void
SymmSCMatrix::assign(SymmSCMatrix*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
SymmSCMatrix::print(ostream&o)
{
  print(0, o, 10);
}

SymmSCMatrix*
SymmSCMatrix::clone()
{
  return dim()->create_symmmatrix();
}

SymmSCMatrix*
SymmSCMatrix::copy()
{
  SymmSCMatrix* result = clone();
  result->assign(this);
  return result;
}

/////////////////////////////////////////////////////////////////////////
// DiagSCMatrix member functions

#define CLASSNAME DiagSCMatrix
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

DiagSCMatrix::DiagSCMatrix()
{
}

DiagSCMatrix::DiagSCMatrix(StateIn&s):
  SavableState(s,class_desc_)
{
}

void
DiagSCMatrix::save_data_state(StateOut&s)
{
}

void *
DiagSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

int
DiagSCMatrix::n()
{
  return dim()->n();
}

double
DiagSCMatrix::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCDiagElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
DiagSCMatrix::assign(double a)
{
  RefSCDiagElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
DiagSCMatrix::assign(const double*a)
{
  int nr = n();
  for (int i=0; i<nr; i++) {
      set_element(i,a[i]);
    }
}

void
DiagSCMatrix::convert(double*a)
{
  int nr = n();
  for (int i=0; i<nr; i++) {
      a[i] = get_element(i);
    }
}

void
DiagSCMatrix::scale(double a)
{
  RefSCDiagElementOp op = new SCElementScale(a);
  this->element_op(op);
}

void
DiagSCMatrix::assign(DiagSCMatrix*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
DiagSCMatrix::print(ostream&o)
{
  print(0, o, 10);
}

DiagSCMatrix*
DiagSCMatrix::clone()
{
  return dim()->create_diagmatrix();
}

DiagSCMatrix*
DiagSCMatrix::copy()
{
  DiagSCMatrix* result = clone();
  result->assign(this);
  return result;
}

/////////////////////////////////////////////////////////////////////////
// These member are used by the abstract SCVector classes.
/////////////////////////////////////////////////////////////////////////

#define CLASSNAME SCVector
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCVector::SCVector()
{
}

SCVector::~SCVector()
{
}

SCVector::SCVector(StateIn&s):
  SavableState(s,class_desc_)
{
}

void
SCVector::save_data_state(StateOut&s)
{
}

void *
SCVector::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

int
SCVector::n()
{
  return dim()->n();
}

double
SCVector::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCRectElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
SCVector::assign(double a)
{
  RefSCVectorElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
SCVector::assign(const double*a)
{
  int nr = n();
  for (int i=0; i<nr; i++) {
      set_element(i,a[i]);
    }
}

void
SCVector::convert(double*a)
{
  int nr = n();
  for (int i=0; i<nr; i++) {
      a[i] = get_element(i);
    }
}

void
SCVector::scale(double a)
{
  RefSCVectorElementOp op = new SCElementScale(a);
  this->element_op(op);
}

void
SCVector::assign(SCVector*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
SCVector::print(ostream&o)
{
  print(0, o, 10);
}

void
SCVector::normalize()
{
  double norm = scalar_product(this);
  norm = 1.0/sqrt(norm);
  scale(norm);
}

SCVector*
SCVector::clone()
{
  return dim()->create_vector();
}

SCVector*
SCVector::copy()
{
  SCVector* result = clone();
  result->assign(this);
  return result;
}
