
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/cmatrix.h>

#define CLASSNAME LocalSCDimension
#define PARENTS public SCDimension
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

#define CLASSNAME LocalSCMatrix
#define PARENTS public SCMatrix
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

#define CLASSNAME LocalSCVector
#define PARENTS public SCVector
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

/////////////////////////////////////////////////////////////////////////////
// LocalSCDimension member functions

LocalSCDimension::LocalSCDimension(int n): n_(n)
{
}

LocalSCDimension::LocalSCDimension(StateIn&s):
  SavableState(s,LocalSCDimension::class_desc_)
{
  s.get(n_);
}

LocalSCDimension::LocalSCDimension(const RefKeyVal&keyval)
{
  n_ = keyval->intvalue("n");
}

void
LocalSCDimension::save_data_state(StateOut&s)
{
  s.put(n_);
}

void *
LocalSCDimension::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCDimension::_castdown(cd);
  return do_castdowns(casts,cd);
}

LocalSCDimension::~LocalSCDimension()
{
}

int
LocalSCDimension::n()
{
  return n_;
}
SCMatrix*
LocalSCDimension::create_matrix(SCDimension*a)
{
  LocalSCDimension*coldim
    = LocalSCDimension::require_castdown(a,"LocalSCDimension::create_matrix");
  return new LocalSCMatrix(this,coldim);
}
SymmSCMatrix*
LocalSCDimension::create_symmmatrix()
{
  return new LocalSymmSCMatrix(this);
}
DiagSCMatrix*
LocalSCDimension::create_diagmatrix()
{
  return new LocalDiagSCMatrix(this);
}
SCVector*
LocalSCDimension::create_vector()
{
  return new LocalSCVector(this);
}

SavableState_REF_def(LocalSCDimension);

/////////////////////////////////////////////////////////////////////////////
// LocalSCMatrix member functions

static double **
init_rect_rows(double *data, int ni,int nj)
{
  double** r = new double*[ni];
  for (int i=0; i<ni; i++) r[i] = &data[i*nj];
  return r;
}

LocalSCMatrix::LocalSCMatrix(): rows(0)
{
}

LocalSCMatrix::LocalSCMatrix(LocalSCDimension*a,LocalSCDimension*b):
  d1(a),
  d2(b),
  rows(0)
{
  resize(a->n(),b->n());
}

LocalSCMatrix::LocalSCMatrix(StateIn&s):
  SCMatrix(s),
  SavableState(s,LocalSCMatrix::class_desc_)
{
  d1.restore_state(s);
  d2.restore_state(s);
  block.restore_state(s);
  rows = init_rect_rows(block->data,d1->n(),d2->n());
}

LocalSCMatrix::LocalSCMatrix(const RefKeyVal&keyval)
{
  d1 = keyval->describedclassvalue("rowdim");
  d2 = keyval->describedclassvalue("coldim");
  d1.require_nonnull();
  d2.require_nonnull();
  block = new SCMatrixRectBlock(0,d1->n(),0,d2->n());
  rows = init_rect_rows(block->data,d1->n(),d2->n());
  for (int i=0; i<nrow(); i++) {
      for (int j=0; j<ncol(); j++) {
          set_element(i,j,keyval->doublevalue("data",i,j));
        }
    }
}

void
LocalSCMatrix::save_data_state(StateOut&s)
{
  SCMatrix::save_data_state(s);
  d1.save_state(s);
  d2.save_state(s);
  block.save_state(s);
}

void *
LocalSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

LocalSCMatrix::~LocalSCMatrix()
{
  if (rows) delete[] rows;
}

int
LocalSCMatrix::compute_offset(int i,int j)
{
  if (i<0 || j<0 || i>=d1->n() || j>=d2->n()) {
      fprintf(stderr,"LocalSCMatrix: index out of bounds\n");
      abort();
    }
  return i*(d2->n()) + j;
}

void
LocalSCMatrix::resize(int nr, int nc)
{
  block = new SCMatrixRectBlock(0,nr,0,nc);
  if (rows) delete[] rows;
  rows = init_rect_rows(block->data,nr,nc);
}

RefSCDimension
LocalSCMatrix::rowdim()
{
  return d1;
}

RefSCDimension
LocalSCMatrix::coldim()
{
  return d2;
}

double
LocalSCMatrix::get_element(int i,int j)
{
  int off = compute_offset(i,j);
  return block->data[off];
}

void
LocalSCMatrix::set_element(int i,int j,double a)
{
  int off = compute_offset(i,j);
  block->data[off] = a;
}

void
LocalSCMatrix::accumulate_product(SCMatrix*a,SCMatrix*b)
{
  const char* name = "LocalSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = LocalSCMatrix::require_castdown(a,name);
  LocalSCMatrix* lb = LocalSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!(this->rowdim() == a->rowdim())
      || !(this->coldim() == b->coldim())
      || !(a->coldim() == b->rowdim())) {
      fprintf(stderr,"LocalSCMatrix::"
              "accumulate_product(SCMatrix*a,SCMatrix*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  cmat_mxm(la->rows, 0,
           lb->rows, 0,
           rows, 0,
           nrow(), la->ncol(), this->ncol(),
           1);
}

// does the outer product a x b.  this must have rowdim() == a->dim() and
// coldim() == b->dim()
void
LocalSCMatrix::accumulate_outer_product(SCVector*a,SCVector*b)
{
  const char* name = "LocalSCMatrix::accumulate_outer_product";
  // make sure that the arguments are of the correct type
  LocalSCVector* la = LocalSCVector::require_castdown(a,name);
  LocalSCVector* lb = LocalSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!(this->rowdim() == a->dim())
      || !(this->coldim() == b->dim())) {
      fprintf(stderr,"LocalSCMatrix::"
              "accumulate_outer_product(SCVector*a,SCVector*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nr = a->n();
  int nc = b->n();
  double* adat = la->block->data;
  double* bdat = lb->block->data;
  double** thisdat = rows;
  for (int i=0; i<nr; i++) {
      for (int j=0; j<nc; j++) {
          thisdat[i][j] += adat[i] * bdat[j];
        }
    }
}

void
LocalSCMatrix::accumulate_product(SCMatrix*a,SymmSCMatrix*b)
{
  const char* name = "LocalSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = LocalSCMatrix::require_castdown(a,name);
  LocalSymmSCMatrix* lb = LocalSymmSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!(this->rowdim() == a->rowdim())
      || !(this->coldim() == b->dim())
      || !(a->coldim() == b->dim())) {
      fprintf(stderr,"LocalSCMatrix::"
              "accumulate_product(SCMatrix*a,SymmSCMatrix*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  double **cd = rows;
  double **ad = la->rows;
  double **bd = lb->rows;
  int ni = a->rowdim().n();
  int njk = b->dim().n();
  for (int i=0; i<ni; i++) {
      for (int j=0; j<njk; j++) {
          for (int k=0; k<=j; k++) {
              cd[i][k] += ad[i][j]*bd[j][k];
            }
          for (; k<njk; k++) {
              cd[i][k] += ad[i][j]*bd[k][j];
            }
        }
    }
}

void
LocalSCMatrix::accumulate_product(SCMatrix*a,DiagSCMatrix*b)
{
  const char* name = "LocalSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = LocalSCMatrix::require_castdown(a,name);
  LocalDiagSCMatrix* lb = LocalDiagSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!(this->rowdim() == a->rowdim())
      || !(this->coldim() == b->dim())
      || !(a->coldim() == b->dim())) {
      fprintf(stderr,"LocalSCMatrix::"
              "accumulate_product(SCMatrix*a,DiagSCMatrix*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  double **cd = rows;
  double **ad = la->rows;
  double *bd = lb->block->data;
  int ni = a->rowdim().n();
  int nj = b->dim().n();
  for (int i=0; i<ni; i++) {
      for (int j=0; j<nj; j++) {
          cd[i][j] += ad[i][j]*bd[j];
        }
    }
}

void
LocalSCMatrix::accumulate(SCMatrix*a)
{
  // make sure that the arguments is of the correct type
  LocalSCMatrix* la
    = LocalSCMatrix::require_castdown(a,"LocalSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!(this->rowdim() == a->rowdim())
      || !(this->coldim() == a->coldim())) {
      fprintf(stderr,"LocalSCMatrix::"
              "accumulate(SCMatrix*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = this->ncol() * this->nrow();
  for (int i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

void
LocalSCMatrix::transpose_this()
{
  cmat_transpose_matrix(rows,nrow(),ncol());
  delete[] rows;
  rows = new double*[ncol()];
  cmat_matrix_pointers(rows,block->data,ncol(),nrow());
  RefLocalSCDimension tmp = d1;
  d1 = d2;
  d2 = tmp;
}

double
LocalSCMatrix::invert_this()
{
  if (nrow() != ncol()) {
      fprintf(stderr,"LocalSCMatrix::invert_this: matrix is not square\n");
      abort();
    }
  return cmat_invert(rows,0,nrow());
}

void
LocalSCMatrix::gen_invert_this()
{
  fprintf(stderr,"LocalSCMatrix::gen_invert_this: SVD not implemented yet");
  abort();
}

double
LocalSCMatrix::determ_this()
{
  if (nrow() != ncol()) {
    fprintf(stderr,"LocalSCMatrix::determ_this: matrix is not square\n");
    abort();
  }
  return cmat_determ(rows,0,nrow());
}

double
LocalSCMatrix::trace()
{
  if (nrow() != ncol()) {
    fprintf(stderr,"LocalSCMatrix::trace: matrix is not square\n");
    abort();
  }
  double ret=0;
  for (int i=0; i < nrow(); i++)
    ret += rows[i][i];
  return ret;
}

double
LocalSCMatrix::solve_this(SCVector*v)
{
  LocalSCVector* lv =
    LocalSCVector::require_castdown(v,"LocalSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!(this->rowdim() == v->dim())) {
      fprintf(stderr,"LocalSCMatrix::solve_this(SCVector*v):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  return cmat_solve_lin(rows,0,lv->block->data,nrow());
}

void
LocalSCMatrix::element_op(const RefSCElementOp& op)
{
  op->process(block.pointer());
}

void
LocalSCMatrix::element_op(const RefSCElementOp2& op,
                          SCMatrix* m)
{
  LocalSCMatrix *lm
      = LocalSCMatrix::require_castdown(m,"LocalSCMatrix::element_op");
  if (!lm || d1 != lm->d1 || d2 != lm->d2) {
      fprintf(stderr,"LocalSCMatrix: bad element_op\n");
      abort();
    }
  op->process(block.pointer(), lm->block.pointer());
}

void
LocalSCMatrix::element_op(const RefSCElementOp3& op,
                          SCMatrix* m,SCMatrix* n)
{
  LocalSCMatrix *lm
      = LocalSCMatrix::require_castdown(m,"LocalSCMatrix::element_op");
  LocalSCMatrix *ln
      = LocalSCMatrix::require_castdown(n,"LocalSCMatrix::element_op");
  if (!lm || !ln
      || d1 != lm->d1 || d2 != lm->d2 || d1 != ln->d1 || d2 != ln->d2) {
      fprintf(stderr,"LocalSCMatrix: bad element_op\n");
      abort();
    }
  op->process(block.pointer(), lm->block.pointer(), ln->block.pointer());
}

// from Ed Seidl at the NIH
void
LocalSCMatrix::print(const char *title, ostream& os, int prec)
{
  int ii,jj,kk,nn;
  int i,j;
  int lwidth,width;
  double max=this->maxabs();

  max=(max==0.0)?1.0:log10(max);
  if(max < 0.0) max=1.0;

  lwidth = prec+5+(int) max; width = 75/lwidth;

  os.setf(ios::fixed,ios::floatfield); os.precision(prec);
  os.setf(ios::right,ios::adjustfield);

  if(title) os << "\n" << title << "\n";
  else os << "\n";

  if(nrow()==0 || ncol()==0) { os << " empty matrix\n"; return; }

  for(ii=jj=0;;) {
    ii++; jj++; kk=width*jj;
    nn=(ncol()>kk)?kk:ncol();

 // print column indices
    for(i=ii; i <= nn; i++) { os.width(lwidth); os << i; }
    os << "\n";

 // print the rows
    for(i=0; i < nrow() ; i++) {
      os.width(5); os << i+1;
      for(j=ii-1; j < nn; j++) { os.width(lwidth); os << rows[i][j]; }
      os << "\n";
      }
    os << "\n";

    if(ncol()<=kk) { os.flush(); return; }
    ii=kk;
    }
}

/////////////////////////////////////////////////////////////////////////////
// LocalSCVector member functions

LocalSCVector::LocalSCVector()
{
}

LocalSCVector::LocalSCVector(LocalSCDimension*a):
  d(a)
{
  resize(a->n());
}

LocalSCVector::LocalSCVector(StateIn&s):
  SCVector(s),
  SavableState(s,LocalSCVector::class_desc_)
{
  d.restore_state(s);
  block.restore_state(s);
}

LocalSCVector::LocalSCVector(const RefKeyVal&keyval)
{
  d = keyval->describedclassvalue("dim");
  d.require_nonnull();
  block = new SCVectorSimpleBlock(0,d->n());
  for (int i=0; i<n(); i++) {
      set_element(i,keyval->doublevalue("data",i,i));
    }
}

void
LocalSCVector::resize(int n)
{
  block = new SCVectorSimpleBlock(0,n);
}

void
LocalSCVector::save_data_state(StateOut&s)
{
  SCVector::save_data_state(s);
  d.save_state(s);
  block.save_state(s);
}

void *
LocalSCVector::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCVector::_castdown(cd);
  return do_castdowns(casts,cd);
}

LocalSCVector::~LocalSCVector()
{
}

RefSCDimension
LocalSCVector::dim()
{
  return d;
}

double
LocalSCVector::get_element(int i)
{
  int size = block->iend - block->istart;
  if (i < 0 || i >= size) {
      fprintf(stderr,"LocalSCVector::get_element: bad index\n");
      abort();
    }
  return block->data[i];
}

void
LocalSCVector::set_element(int i,double a)
{
  int size = block->iend - block->istart;
  if (i < 0 || i >= size) {
      fprintf(stderr,"LocalSCVector::set_element: bad index\n");
      abort();
    }
  block->data[i] = a;
}

void
LocalSCVector::accumulate_product(SCMatrix*a,SCVector*b)
{
  const char* name = "LocalSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = LocalSCMatrix::require_castdown(a,name);
  LocalSCVector* lb = LocalSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!(this->dim() == a->rowdim())
      || !(a->coldim() == b->dim())) {
      fprintf(stderr,"LocalSCVector::"
              "accumulate_product(SCMatrix*a,SCVector*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  cmat_mxm(la->rows, 0,
           &(lb->block->data), 1,
           &(block->data), 1,
           n(), la->ncol(), 1,
           1);
}

void
LocalSCVector::accumulate_product(SymmSCMatrix*a,SCVector*b)
{
  const char* name = "LocalSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSymmSCMatrix* la = LocalSymmSCMatrix::require_castdown(a,name);
  LocalSCVector* lb = LocalSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!(this->dim() == a->dim())
      || !(a->dim() == b->dim())) {
      fprintf(stderr,"LocalSCVector::"
              "accumulate_product(SymmSCMatrix*a,SCVector*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  double* thisdat = block->data;
  double** adat = la->rows;
  double* bdat = lb->block->data;
  int n = dim()->n();
  for (int i=0; i<n; i++) {
      int j;
      double tmp = 0.0;
      for (j=0; j<=i; j++) {
          tmp += adat[i][j] * bdat[j];
        }
      for (; j<n; j++) {
          tmp += adat[j][i] * bdat[j];
        }
      thisdat[i] += tmp;
    }
}

void
LocalSCVector::accumulate(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::accumulate");

  // make sure that the dimensions match
  if (!(this->dim() == la->dim())) {
      fprintf(stderr,"LocalSCVector::"
              "accumulate(SCVector*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = d->n();
  for (int i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

void
LocalSCVector::assign(double a)
{
  int nelem = d->n();
  for (int i=0; i<nelem; i++) block->data[i] = a;
}

void
LocalSCVector::assign(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::assign");

  // make sure that the dimensions match
  if (!(this->dim() == la->dim())) {
      fprintf(stderr,"LocalSCVector::"
              "assign(SCVector*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = d->n();
  for (int i=0; i<nelem; i++) block->data[i] = la->block->data[i];
}

void
LocalSCVector::assign(const double*a)
{
  int nelem = d->n();
  for (int i=0; i<nelem; i++) block->data[i] = a[i];
}

double
LocalSCVector::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::scalar_product");

  // make sure that the dimensions match
  if (!(this->dim() == la->dim())) {
      fprintf(stderr,"LocalSCVector::"
              "scale_product(SCVector*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = d->n();
  double result = 0.0;
  for (int i=0; i<nelem; i++) result += block->data[i] * la->block->data[i];
  return result;
}

void
LocalSCVector::element_op(const RefSCElementOp& op)
{
  op->process(block.pointer());
}

void
LocalSCVector::element_op(const RefSCElementOp2& op,
                          SCVector* m)
{
  LocalSCVector *lm
      = LocalSCVector::require_castdown(m, "LocalSCVector::element_op");
  if (!lm || d != lm->d) {
      fprintf(stderr,"LocalSCVector: bad element_op\n");
      abort();
    }
  op->process(block.pointer(), lm->block.pointer());
}

void
LocalSCVector::element_op(const RefSCElementOp3& op,
                          SCVector* m,SCVector* n)
{
  LocalSCVector *lm
      = LocalSCVector::require_castdown(m, "LocalSCVector::element_op");
  LocalSCVector *ln
      = LocalSCVector::require_castdown(n, "LocalSCVector::element_op");
  if (!lm || !ln || d != lm->d || d != ln->d) {
      fprintf(stderr,"LocalSCVector: bad element_op\n");
      abort();
    }
  op->process(block.pointer(), lm->block.pointer(), ln->block.pointer());
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
LocalSCVector::print(const char *title, ostream& os, int prec)
{
  int i;
  int lwidth;
  double max=this->maxabs();

  max=(max==0.0)?1.0:log10(max);
  if(max < 0.0) max=1.0;

  lwidth = prec+5+(int) max;

  os.setf(ios::fixed,ios::floatfield); os.precision(prec);
  os.setf(ios::right,ios::adjustfield);

  if(title) os << "\n" << title << "\n";
  else os << "\n";

  if(n()==0) { os << " empty vector\n"; return; }

  for (i=0; i<n(); i++) {
      os.width(5); os << i+1;
      os.width(lwidth); os << block->data[i];
      os << "\n";
    }
  os << "\n";

  os.flush();
}
