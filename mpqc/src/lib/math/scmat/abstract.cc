
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
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
SCMatrixKit::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCMatrixKit::SCMatrixKit()
{
}

SCMatrixKit::SCMatrixKit(const RefKeyVal&)
{
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
// SCMatrix members

#define CLASSNAME SCMatrix
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
SCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCMatrix::SCMatrix(const RefSCDimension&a1, const RefSCDimension&a2,
                   SCMatrixKit*kit):
  d1(a1),
  d2(a2),
  kit_(kit)
{
}

SCMatrix::~SCMatrix()
{
}

void
SCMatrix::save(StateOut&s)
{
  int nr = nrow();
  int nc = ncol();
  s.put(nr);
  s.put(nc);
  int has_subblocks = 0;
  s.put(has_subblocks);
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
  int has_subblocks;
  s.get(has_subblocks);
  if (!has_subblocks) {
      for (int i=0; i<nr; i++) {
          for (int j=0; j<nc; j++) {
              double tmp;
              s.get(tmp);
              set_element(i,j, tmp);
            }
        }
    }
  else {
      cerr << "SCMatrix::restore(): matrix has subblocks--cannot restore"
           << endl;
      abort();
    }
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
SCMatrix::convert(SCMatrix*a)
{
  assign(0.0);
  convert_accumulate(a);
}

void
SCMatrix::convert_accumulate(SCMatrix*a)
{
  RefSCElementOp op = new SCElementAccumulateSCMatrix(a);
  element_op(op);
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
  return kit()->matrix(rowdim(),coldim());
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
  SCMatrix *brect = kit()->matrix(b->dim(),b->dim());
  brect->assign(0.0);
  brect->accumulate(b);
  accumulate_product(a,brect);
  delete brect;
}

void
SCMatrix::accumulate_product(SCMatrix*a,DiagSCMatrix*b)
{
  SCMatrix *brect = kit()->matrix(b->dim(),b->dim());
  brect->assign(0.0);
  brect->accumulate(b);
  accumulate_product(a,brect);
  delete brect;
}

/////////////////////////////////////////////////////////////////////////
// SymmSCMatrix member functions

#define CLASSNAME SymmSCMatrix
#define PARENTS public DescribedClass
#include <util/class/classia.h>

SymmSCMatrix::SymmSCMatrix(const RefSCDimension&a, SCMatrixKit *kit):
  d(a),
  kit_(kit)
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
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
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
SymmSCMatrix::convert(SymmSCMatrix*a)
{
  assign(0.0);
  convert_accumulate(a);
}

void
SymmSCMatrix::convert_accumulate(SymmSCMatrix*a)
{
  RefSCElementOp op = new SCElementAccumulateSymmSCMatrix(a);
  element_op(op);
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
  RefSCMatrix m = kit()->matrix(dim(),dim());
  m->assign(0.0);
  m->accumulate(this);
  m->print(title, out, i);
}

SymmSCMatrix*
SymmSCMatrix::clone()
{
  return kit()->symmmatrix(dim());
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
  RefSCMatrix m = kit()->matrix(a->rowdim(),a->rowdim());
  m->assign(0.0);
  m->accumulate_product(a, at.pointer());
  scale(2.0);
  accumulate_symmetric_sum(m.pointer());
  scale(0.5);
}

void
SymmSCMatrix::accumulate_transform(SCMatrix *a, SymmSCMatrix *b)
{
  RefSCMatrix m = kit()->matrix(a->rowdim(),a->rowdim());
  RefSCMatrix brect = kit()->matrix(b->dim(),b->dim());
  brect->assign(0.0);
  brect->accumulate(b);

  RefSCMatrix tmp = a->clone();
  tmp->assign(0.0);
  tmp->accumulate_product(a,brect);
  brect = 0;

  RefSCMatrix at = a->copy();
  at->assign(0.0);
  at->transpose_this();

  RefSCMatrix res = kit()->matrix(dim(),dim());
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
  RefSCMatrix m = kit()->matrix(a->rowdim(),a->rowdim());
  RefSCMatrix brect = kit()->matrix(b->dim(),b->dim());
  brect->assign(0.0);
  brect->accumulate(b);

  RefSCMatrix tmp = a->clone();
  tmp->assign(0.0);
  tmp->accumulate_product(a,brect);
  brect = 0;

  RefSCMatrix at = a->copy();
  at->transpose_this();

  RefSCMatrix res = kit()->matrix(dim(),dim());
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
  RefSCMatrix m = kit()->matrix(a->dim(),a->dim());
  m->assign(0.0);
  m->accumulate(a);
  accumulate_transform(m.pointer(),b);
}

void
SymmSCMatrix::accumulate_symmetric_outer_product(SCVector *v)
{
  RefSCMatrix m = kit()->matrix(dim(),dim());
  m->assign(0.0);
  m->accumulate_outer_product(v,v);

  scale(2.0);
  accumulate_symmetric_sum(m.pointer());
  scale(0.5);
}

double
SymmSCMatrix::scalar_product(SCVector *v)
{
  RefSCVector v2 = kit()->vector(dim());
  v2->assign(0.0);
  v2->accumulate_product(this,v);
  return v2->scalar_product(v);
}

/////////////////////////////////////////////////////////////////////////
// DiagSCMatrix member functions

#define CLASSNAME DiagSCMatrix
#define PARENTS public DescribedClass
#include <util/class/classia.h>

void *
DiagSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

DiagSCMatrix::DiagSCMatrix(const RefSCDimension&a, SCMatrixKit *kit):
  d(a),
  kit_(kit)
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
DiagSCMatrix::convert(DiagSCMatrix*a)
{
  assign(0.0);
  convert_accumulate(a);
}

void
DiagSCMatrix::convert_accumulate(DiagSCMatrix*a)
{
  RefSCElementOp op = new SCElementAccumulateDiagSCMatrix(a);
  element_op(op);
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
  RefSCMatrix m = kit()->matrix(dim(),dim());
  m->assign(0.0);
  m->accumulate(this);
  m->print(title, out, i);
}

DiagSCMatrix*
DiagSCMatrix::clone()
{
  return kit()->diagmatrix(dim());
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
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
SCVector::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCVector::SCVector(const RefSCDimension&a, SCMatrixKit *kit):
  d(a),
  kit_(kit)
{
}

SCVector::~SCVector()
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
SCVector::convert(SCVector*a)
{
  assign(0.0);
  convert_accumulate(a);
}

void
SCVector::convert_accumulate(SCVector*a)
{
  RefSCElementOp op = new SCElementAccumulateSCVector(a);
  element_op(op);
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
  return kit()->vector(dim());
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
  RefSCMatrix mrect = kit()->matrix(m->dim(),m->dim());
  mrect->assign(0.0);
  mrect->accumulate(m);
  accumulate_product(mrect.pointer(), v);
}
