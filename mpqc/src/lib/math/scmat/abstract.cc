
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blkiter.h>

/////////////////////////////////////////////////////////////////////////
// These member are used by the abstract SCMatrix classes.
/////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////
// SCDimension members

#define CLASSNAME SCDimension
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCDimension::SCDimension(const char* name)
{
  if (name) name_ = strcpy(new char[strlen(name)+1], name);
  else name_ = 0;
}

SCDimension::SCDimension(StateIn&s):
  SavableState(s)
{
  s.getstring(name_);
}

void
SCDimension::save_data_state(StateOut&s)
{
  s.putstring(name_);
}

void *
SCDimension::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCMatrix*
SCDimension::create_matrix(const RefSCDimension&d)
{
  return create_matrix(d.pointer());
}

SCDimension::~SCDimension()
{
  if (name_) delete[] name_;
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

/////////////////////////////////////////////////////////////////////////
// SCMatrix members

#define CLASSNAME SCMatrix
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCMatrix::SCMatrix()
{
}

SCMatrix::~SCMatrix()
{
}

SCMatrix::SCMatrix(StateIn&s):
  SavableState(s)
{
}

void
SCMatrix::save_data_state(StateOut&s)
{
}

void *
SCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

int
SCMatrix::nrow()
{
  int result = rowdim()->n();
  return result;
}

int
SCMatrix::ncol()
{
  int result = coldim()->n();
  return result;
}

double
SCMatrix::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
SCMatrix::assign(double a)
{
  RefSCElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
SCMatrix::scale(double a)
{
  RefSCElementOp op = new SCElementScale(a);
  this->element_op(op);
}

void
SCMatrix::shift_diagonal(double a)
{
  RefSCElementOp op = new SCElementShiftDiagonal(a);
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
  int i;
  int nr = nrow();
  int nc = ncol();
  // some compilers need the following cast
  const double **v = (const double**) new double*[nr];
  for (i=0; i<nr; i++) {
      v[i] = &a[i*nc];
    }
  assign(v);
  delete[] v;
}

void
SCMatrix::assign(const double**a)
{
  int i;
  int j;
  int nr;
  int nc;
  nr = nrow();
  nc = ncol();
  for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
          set_element(i,j,a[i][j]);
        }
    }
}

void
SCMatrix::convert(double*a)
{
  int i;
  int nr = nrow();
  int nc = ncol();
  double **v = new double*[nr];
  for (i=0; i<nr; i++) {
      v[i] = &a[i*nc];
    }
  convert(v);
  delete[] v;
}

void
SCMatrix::convert(double**a)
{
  int i, j;
  int nr, nc;
  nr = nrow();
  nc = ncol();
  for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
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
  return rowdim()->create_matrix(coldim().pointer());
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
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SymmSCMatrix::SymmSCMatrix()
{
}

SymmSCMatrix::SymmSCMatrix(StateIn&s):
  SavableState(s)
{
}

void
SymmSCMatrix::save_data_state(StateOut&s)
{
}

void *
SymmSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

int
SymmSCMatrix::n()
{
  int result = dim()->n();
  return result;
}

double
SymmSCMatrix::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
SymmSCMatrix::assign(double a)
{
  RefSCElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
SymmSCMatrix::assign(const double*a)
{
  int i;
  int nr = n();
  // some compilers need the following cast
  const double **v = (const double **) new double*[nr];
  int ioff= 0;
  for (i=0; i<nr; i++) {
      ioff += i;
      v[i] = &a[ioff];
//    ioff += i;
    }
  assign(v);
  delete[] v;
}

void
SymmSCMatrix::assign(const double**a)
{
  int i;
  int j;
  int nr;
  nr = n();
  for (i=0; i<nr; i++) {
      for (j=0; j<=i; j++) {
          set_element(i,j,a[i][j]);
        }
    }
}

void
SymmSCMatrix::convert(double*a)
{
  int i;
  int nr = n();
  double **v = new double*[nr];
  int ioff= 0;
  for (i=0; i<nr; i++) {
      ioff += i;
      v[i] = &a[ioff];
//    ioff += i;
    }
  convert(v);
  delete[] v;
}

void
SymmSCMatrix::convert(double**a)
{
  int i;
  int j;
  int nr;
  nr = n();
  for (i=0; i<nr; i++) {
      for (j=0; j<=i; j++) {
          a[i][j] = get_element(i,j);
        }
    }
}

void
SymmSCMatrix::scale(double a)
{
  RefSCElementOp op = new SCElementScale(a);
  this->element_op(op);
}

void
SymmSCMatrix::shift_diagonal(double a)
{
  RefSCElementOp op = new SCElementShiftDiagonal(a);
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
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

DiagSCMatrix::DiagSCMatrix()
{
}

DiagSCMatrix::DiagSCMatrix(StateIn&s):
  SavableState(s)
{
}

void
DiagSCMatrix::save_data_state(StateOut&s)
{
}

void *
DiagSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

int
DiagSCMatrix::n()
{
  int result = dim()->n();
  return result;
}

double
DiagSCMatrix::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
DiagSCMatrix::assign(double a)
{
  RefSCElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
DiagSCMatrix::assign(const double*a)
{
  int i;
  int nr = n();
  for (i=0; i<nr; i++) {
      set_element(i,a[i]);
    }
}

void
DiagSCMatrix::convert(double*a)
{
  int i;
  int nr = n();
  for (i=0; i<nr; i++) {
      a[i] = get_element(i);
    }
}

void
DiagSCMatrix::scale(double a)
{
  RefSCElementOp op = new SCElementScale(a);
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
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCVector::SCVector()
{
}

SCVector::~SCVector()
{
}

SCVector::SCVector(StateIn&s):
  SavableState(s)
{
}

void
SCVector::save_data_state(StateOut&s)
{
}

void *
SCVector::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

int
SCVector::n()
{
  int result = dim()->n();
  return result;
}

double
SCVector::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
SCVector::assign(double a)
{
  RefSCElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
SCVector::assign(const double*a)
{
  int i;
  int nr = n();
  for (i=0; i<nr; i++) {
      set_element(i,a[i]);
    }
}

void
SCVector::convert(double*a)
{
  int i;
  int nr = n();
  for (i=0; i<nr; i++) {
      a[i] = get_element(i);
    }
}

void
SCVector::scale(double a)
{
  RefSCElementOp op = new SCElementScale(a);
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
  double norm = sqrt(scalar_product(this));
  if (norm > 1.e-20) norm = 1.0/norm;
  else {
      fprintf(stderr,"SCVector::normalize: tried to normalize tiny vector\n");
      abort();
    }
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
