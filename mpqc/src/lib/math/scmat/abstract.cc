
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////
// These member are used by the abstract SCMatrix classes.
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
// SCMatrixKit members

#define CLASSNAME SCMatrixKit
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
SCMatrixKit::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCMatrixKit::SCMatrixKit()
{
}

SCMatrixKit::SCMatrixKit(const RefKeyVal&)
{
}

SCMatrixKit::SCMatrixKit(StateIn&s):
  SavableState(s)
{
}

void
SCMatrixKit::save_data_state(StateOut&s)
{
}

SCMatrix*
SCMatrixKit::matrix(const RefSCDimension&d1, const RefSCDimension&d2)
{
  return d1->create_matrix(d2.pointer());
}

SymmSCMatrix*
SCMatrixKit::symmmatrix(const RefSCDimension&d)
{
  return d->create_symmmatrix();
}

DiagSCMatrix*
SCMatrixKit::diagmatrix(const RefSCDimension&d)
{
  return d->create_diagmatrix();
}

SCVector*
SCMatrixKit::vector(const RefSCDimension&d)
{
  return d->create_vector();
}

SCMatrixKit::~SCMatrixKit()
{
}

SCMatrix*
SCMatrixKit::restore_matrix(StateIn& s,
                            const RefSCDimension& d1,
                            const RefSCDimension& d2)
{
  SCMatrix *r = matrix(d1,d2);
  r->restore(s);
  return r;
}

SymmSCMatrix*
SCMatrixKit::restore_symmmatrix(StateIn& s, const RefSCDimension& d)
{
  SymmSCMatrix *r = symmmatrix(d);
  r->restore(s);
  return r;
}

DiagSCMatrix*
SCMatrixKit::restore_diagmatrix(StateIn& s, const RefSCDimension& d)
{
  DiagSCMatrix *r = diagmatrix(d);
  r->restore(s);
  return r;
}

SCVector*
SCMatrixKit::restore_vector(StateIn& s, const RefSCDimension& d)
{
  SCVector *r = vector(d);
  r->restore(s);
  return r;
}

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

void
SCMatrix::save(StateOut&s)
{
  int nr = nrow();
  int nc = ncol();
  s.put(nr);
  s.put(nc);
  for (int i=0; i<nr; i++) {
      for (int j=0; j<nc; j++) {
          s.put(get_element(i,j));
        }
    }
}

void
SCMatrix::restore(StateIn& s)
{
  int nrt, nr = nrow();
  int nct, nc = ncol();
  s.get(nrt);
  s.get(nct);
  if (nrt != nr || nct != nc) {
      cerr << "SCMatrix::restore(): bad dimensions" << endl;
      abort();
    }
  for (int i=0; i<nr; i++) {
      for (int j=0; j<nc; j++) {
          double tmp;
          s.get(tmp);
          set_element(i,j, tmp);
        }
    }
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
SCMatrix::randomize()
{
  RefSCElementOp op = new SCElementRandomize();
  this->element_op(op);
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
SCMatrix::scale_diagonal(double a)
{
  RefSCElementOp op = new SCElementScaleDiagonal(a);
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

void
SCMatrix::svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V)
{
  fprintf(stderr,"%s: SVD not implemented\n", class_name());
  abort();
}

void
SCMatrix::accumulate_product(SCMatrix*a,SymmSCMatrix*b)
{
  SCMatrix *brect = b->dim()->create_matrix(b->dim());
  brect->assign(0.0);
  brect->accumulate(b);
  accumulate_product(a,brect);
  delete brect;
}

void
SCMatrix::accumulate_product(SCMatrix*a,DiagSCMatrix*b)
{
  SCMatrix *brect = b->dim()->create_matrix(b->dim());
  brect->assign(0.0);
  brect->accumulate(b);
  accumulate_product(a,brect);
  delete brect;
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

void
SymmSCMatrix::save(StateOut&s)
{
  int nr = n();
  s.put(nr);
  for (int i=0; i<nr; i++) {
      for (int j=0; j<=i; j++) {
          s.put(get_element(i,j));
        }
    }
}

void
SymmSCMatrix::restore(StateIn& s)
{
  int nrt, nr = n();
  s.get(nrt);
  if (nrt != nr) {
      cerr << "SymmSCMatrix::restore(): bad dimension" << endl;
      abort();
    }
  for (int i=0; i<nr; i++) {
      for (int j=0; j<=i; j++) {
          double tmp;
          s.get(tmp);
          set_element(i,j, tmp);
        }
    }
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
SymmSCMatrix::randomize()
{
  RefSCElementOp op = new SCElementRandomize();
  this->element_op(op);
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
SymmSCMatrix::scale_diagonal(double a)
{
  RefSCElementOp op = new SCElementScaleDiagonal(a);
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

void
SymmSCMatrix::print(const char* title, ostream& out, int i)
{
  RefSCMatrix m = dim()->create_matrix(dim());
  m->assign(0.0);
  m->accumulate(this);
  m->print(title, out, i);
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

void
SymmSCMatrix::accumulate_symmetric_product(SCMatrix *a)
{
  RefSCMatrix at = a->copy();
  at->transpose_this();
  RefSCMatrix m = a->rowdim()->create_matrix(a->rowdim());
  m->assign(0.0);
  m->accumulate_product(a, at.pointer());
  scale(2.0);
  accumulate_symmetric_sum(m.pointer());
  scale(0.5);
}

void
SymmSCMatrix::accumulate_transform(SCMatrix *a, SymmSCMatrix *b)
{
  RefSCMatrix m = a->rowdim()->create_matrix(a->rowdim());
  RefSCMatrix brect = b->dim()->create_matrix(b->dim());
  brect->assign(0.0);
  brect->accumulate(b);

  RefSCMatrix tmp = a->clone();
  tmp->assign(0.0);
  tmp->accumulate_product(a,brect);
  brect = 0;

  RefSCMatrix at = a->copy();
  at->assign(0.0);
  at->transpose_this();

  RefSCMatrix res = dim()->create_matrix(dim());
  res->assign(0.0);
  res->accumulate_product(tmp.pointer(), at.pointer());
  tmp = 0;
  at = 0;

  scale(2.0);
  accumulate_symmetric_sum(res.pointer());
  scale(0.5);
}

void
SymmSCMatrix::accumulate_transform(SCMatrix *a, DiagSCMatrix *b)
{
  RefSCMatrix m = a->rowdim()->create_matrix(a->rowdim());
  RefSCMatrix brect = b->dim()->create_matrix(b->dim());
  brect->assign(0.0);
  brect->accumulate(b);

  RefSCMatrix tmp = a->clone();
  tmp->assign(0.0);
  tmp->accumulate_product(a,brect);
  brect = 0;

  RefSCMatrix at = a->copy();
  at->transpose_this();

  RefSCMatrix res = dim()->create_matrix(dim());
  res->assign(0.0);
  res->accumulate_product(tmp.pointer(), at.pointer());
  tmp = 0;
  at = 0;

  scale(2.0);
  accumulate_symmetric_sum(res.pointer());
  scale(0.5);
}

void
SymmSCMatrix::accumulate_transform(SymmSCMatrix *a, SymmSCMatrix *b)
{
  RefSCMatrix m = a->dim()->create_matrix(a->dim());
  m->assign(0.0);
  m->accumulate(a);
  accumulate_transform(m.pointer(),b);
}

void
SymmSCMatrix::accumulate_symmetric_outer_product(SCVector *v)
{
  RefSCMatrix m = dim()->create_matrix(dim());
  m->assign(0.0);
  m->accumulate_outer_product(v,v);

  scale(2.0);
  accumulate_symmetric_sum(m.pointer());
  scale(0.5);
}

double
SymmSCMatrix::scalar_product(SCVector *v)
{
  RefSCVector v2 = dim()->create_vector();
  v2->assign(0.0);
  v2->accumulate_product(this,v);
  return v2->scalar_product(v);
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

void
DiagSCMatrix::save(StateOut&s)
{
  int nr = n();
  s.put(nr);
  for (int i=0; i<nr; i++) {
      s.put(get_element(i));
    }
}

void
DiagSCMatrix::restore(StateIn& s)
{
  int nrt, nr = n();
  s.get(nrt);
  if (nrt != nr) {
      cerr << "DiagSCMatrix::restore(): bad dimension" << endl;
      abort();
    }
  for (int i=0; i<nr; i++) {
      double tmp;
      s.get(tmp);
      set_element(i, tmp);
    }
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
DiagSCMatrix::randomize()
{
  RefSCElementOp op = new SCElementRandomize();
  this->element_op(op);
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

void
DiagSCMatrix::print(const char* title, ostream& out, int i)
{
  RefSCMatrix m = dim()->create_matrix(dim());
  m->assign(0.0);
  m->accumulate(this);
  m->print(title, out, i);
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

void
SCVector::save(StateOut&s)
{
  int nr = n();
  s.put(nr);
  for (int i=0; i<nr; i++) {
      s.put(get_element(i));
    }
}

void
SCVector::restore(StateIn& s)
{
  int nrt, nr = n();
  s.get(nrt);
  if (nrt != nr) {
      cerr << "SCVector::restore(): bad dimension" << endl;
      abort();
    }
  for (int i=0; i<nr; i++) {
      double tmp;
      s.get(tmp);
      set_element(i, tmp);
    }
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
SCVector::randomize()
{
  RefSCElementOp op = new SCElementRandomize();
  this->element_op(op);
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

void
SCVector::accumulate_product(SymmSCMatrix *m, SCVector *v)
{
  RefSCMatrix mrect = m->dim()->create_matrix(m->dim());
  mrect->assign(0.0);
  mrect->accumulate(m);
  accumulate_product(mrect.pointer(), v);
}
