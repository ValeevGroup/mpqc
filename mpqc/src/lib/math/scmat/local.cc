
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/cmatrix.h>

/////////////////////////////////////////////////////////////////////////////
// LocalSCDimension member functions

#define CLASSNAME LocalSCDimension
#define PARENTS public SCDimension
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

LocalSCDimension::LocalSCDimension(int n): n_(n)
{
}

LocalSCDimension::LocalSCDimension(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(n_);
}

LocalSCDimension::LocalSCDimension(KeyVal&keyval)
{
  n_ = keyval.intvalue("n");
}

void
LocalSCDimension::save_data_state(StateOut&s)
{
  s.put(n_);
}

void *
LocalSCDimension::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCDimension::_castdown(cd) };
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

#define CLASSNAME LocalSCMatrix
#define PARENTS public SCMatrix
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

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
  rows(0),
  d1(a),
  d2(b)
{
  resize(a->n(),b->n());
}

LocalSCMatrix::LocalSCMatrix(StateIn&s):
  SCMatrix(s),
  SavableState(s,class_desc_)
{
  d1.restore_state(s);
  d2.restore_state(s);
  block.restore_state(s);
  rows = init_rect_rows(block->data,d1->n(),d2->n());
}

LocalSCMatrix::LocalSCMatrix(KeyVal&keyval)
{
  d1 = keyval.describedclassvalue("rowdim");
  d2 = keyval.describedclassvalue("coldim");
  d1.require_nonnull();
  d2.require_nonnull();
  block = new SCMatrixRectBlock(0,d1->n(),0,d2->n());
  rows = init_rect_rows(block->data,d1->n(),d2->n());
  for (int i=0; i<nrow(); i++) {
      for (int j=0; j<ncol(); j++) {
          set_element(i,j,keyval.doublevalue("data",i,j));
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
  void* casts[] =  { SCMatrix::_castdown(cd) };
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
  return block->data[compute_offset(i,j)];
}

void
LocalSCMatrix::set_element(int i,int j,double a)
{
  block->data[compute_offset(i,j)] = a;
}

void
LocalSCMatrix::accumulate_product(SCMatrix*a,SCMatrix*b)
{
  const char* name = "LocalSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = LocalSCMatrix::require_castdown(a,name);
  LocalSCMatrix* lb = LocalSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (this->rowdim() != a->rowdim()
      || this->coldim() != b->coldim()
      || a->coldim() != b->rowdim()) {
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

void
LocalSCMatrix::accumulate(SCMatrix*a)
{
  // make sure that the arguments is of the correct type
  LocalSCMatrix* la
    = LocalSCMatrix::require_castdown(a,"LocalSCMatrix::accumulate");

  // make sure that the dimensions match
  if (this->rowdim() != a->rowdim()
      || this->coldim() != a->coldim()) {
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
LocalSCMatrix::element_op(RefSCRectElementOp& op)
{
  op->process(block);
}

// from Ed Seidl at the NIH
void
LocalSCMatrix::print(const char *title, ostream& os, int prec)
{
  int ii,jj,kk,nn,ll;
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
    ll=2*(nn-ii+1)+1;

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
// LocalSymmSCMatrix member functions

#define CLASSNAME LocalSymmSCMatrix
#define PARENTS public SymmSCMatrix
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

static double **
init_symm_rows(double *data, int n)
{
  double** r = new double*[n];
  for (int i=0; i<n; i++) r[i] = &data[(i*(i+1))/2];
  return r;
}

LocalSymmSCMatrix::LocalSymmSCMatrix(): rows(0)
{
}

LocalSymmSCMatrix::LocalSymmSCMatrix(LocalSCDimension*a):
  rows(0),
  d(a)
{
  resize(a->n());
}

LocalSymmSCMatrix::LocalSymmSCMatrix(StateIn&s):
  SymmSCMatrix(s),
  SavableState(s,class_desc_)
{
  d.restore_state(s);
  block.restore_state(s);
  if (rows) delete[] rows;
  rows = init_symm_rows(block->data,d->n());
}

LocalSymmSCMatrix::LocalSymmSCMatrix(KeyVal&keyval)
{
  d = keyval.describedclassvalue("dim");
  d.require_nonnull();
  block = new SCMatrixLTriBlock(0,d->n());
  rows = init_symm_rows(block->data,d->n());
  for (int i=0; i<n(); i++) {
      for (int j=0; j<=i; j++) {
          set_element(i,j,keyval.doublevalue("data",i,j));
        }
    }
}

void
LocalSymmSCMatrix::save_data_state(StateOut&s)
{
  SymmSCMatrix::save_data_state(s);
  d.save_state(s);
  block.save_state(s);
}

void *
LocalSymmSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SymmSCMatrix::_castdown(cd) };
  return do_castdowns(casts,cd);
}

LocalSymmSCMatrix::~LocalSymmSCMatrix()
{
  if (rows) delete[] rows;
}

int
LocalSymmSCMatrix::compute_offset(int i,int j)
{
  if (i<0 || j<0 || i>=d->n() || j>=d->n()) {
      fprintf(stderr,"LocalSymmSCMatrix: index out of bounds\n");
      abort();
    }
  if (i<j) {
      int tmp = j; j=i; i=tmp;
    }
  return (i*(i+1))/2 + j;
}

void
LocalSymmSCMatrix::resize(int n)
{
  block = new SCMatrixLTriBlock(0,n);
  rows = init_symm_rows(block->data,n);
}

RefSCDimension
LocalSymmSCMatrix::dim()
{
  return d;
}

double
LocalSymmSCMatrix::get_element(int i,int j)
{
  return block->data[compute_offset(i,j)];
}

void
LocalSymmSCMatrix::set_element(int i,int j,double a)
{
  block->data[compute_offset(i,j)] = a;
}

void
LocalSymmSCMatrix::accumulate(SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  LocalSymmSCMatrix* la
    = LocalSymmSCMatrix::require_castdown(a,"LocalSymmSCMatrix::accumulate");

  // make sure that the dimensions match
  if (this->dim() != la->dim()) {
      fprintf(stderr,"LocalSymmSCMatrix::"
              "accumulate(SCMatrix*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = (this->n() * (this->n() + 1))/2;
  for (int i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

double
LocalSymmSCMatrix::invert_this()
{
  return cmat_invert(rows,1,n());
}

void
LocalSymmSCMatrix::diagonalize(DiagSCMatrix*a,SCMatrix*b)
{
  const char* name = "LocalSymmSCMatrix::diagonalize";
  // make sure that the arguments is of the correct type
  LocalDiagSCMatrix* la = LocalDiagSCMatrix::require_castdown(a,name);
  LocalSCMatrix* lb = LocalSCMatrix::require_castdown(b,name);

  if (   (la&&(la->dim() != dim()))
      || (lb&&(lb->coldim() != dim() || lb->rowdim() != dim()))) {
      fprintf(stderr,"LocalSymmSCMatrix::"
              "diagonalize(DiagSCMatrix*a,SCMatrix*b): bad dims");
      abort();
    }

  double *eigvals;
  double **eigvecs;
  if (!la) {
      eigvals = new double[n()];
    }
  else {
      eigvals = la->block->data;
    }

  if (!lb) {
      eigvecs = cmat_new_square_matrix(n());
    }
  else {
      eigvecs = lb->rows;
    }

  cmat_diag(rows,eigvals,eigvecs,n(),1,1.0e-15);

  if (!la) delete[] eigvals;
  if (!lb) cmat_delete_matrix(eigvecs);
}

// computes this += a * a.t
void
LocalSymmSCMatrix::accumulate_symmetric_product(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  LocalSCMatrix* la
    = LocalSCMatrix::require_castdown(a,"LocalSymmSCMatrix::"
                                          "accumulate_symmetric_product");

  if (la->rowdim() != dim()) {
      fprintf(stderr,"LocalSymmSCMatrix::"
              "accumulate_symmetric_product(SCMatrix*a): bad dim");
      abort();
    }

  cmat_symmetric_mxm(rows,n(),la->rows,la->ncol(),1);
}

// this += a * b * transpose(a)
void
LocalSymmSCMatrix::accumulate_transform(SCMatrix*a,SymmSCMatrix*b)
{
  // do the necessary castdowns
  LocalSCMatrix*la
    = LocalSCMatrix::require_castdown(a,"%s::accumulate_transform",
                                      class_name());
  LocalSymmSCMatrix*lb = require_castdown(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (n() != la->nrow() || la->ncol() != lb->n()) {
      fprintf(stderr,"LocalSymmSCMatrix::accumulate_transform: bad dim\n");
      abort();
    }

  cmat_transform_symmetric_matrix(rows,n(),lb->rows,lb->n(),la->rows,1);
}

void
LocalSymmSCMatrix::element_op(RefSCSymmElementOp& op)
{
  op->process(block);
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
LocalSymmSCMatrix::print(const char *title, ostream& os, int prec)
{
  int ii,jj,kk,nn,ll;
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

  if(n()==0) { os << " empty matrix\n"; return; }

  for(ii=jj=0;;) {
    ii++; jj++; kk=width*jj;
    nn=(n()>kk)?kk:n();
    ll=2*(nn-ii+1)+1;

 // print column indices
    for(i=ii; i <= nn; i++) { os.width(lwidth); os << i; }
    os << "\n";

 // print the rows
    for(i=0; i < n() ; i++) {
      os.width(5); os << i+1;
      for(j=ii-1; j<nn && j<=i; j++) { os.width(lwidth); os << rows[i][j]; }
      os << "\n";
      }
    os << "\n";

    if(n()<=kk) { os.flush(); return; }
    ii=kk;
    }
}

/////////////////////////////////////////////////////////////////////////////
// LocalDiagSCMatrix member functions

#define CLASSNAME LocalDiagSCMatrix
#define PARENTS public DiagSCMatrix
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

LocalDiagSCMatrix::LocalDiagSCMatrix()
{
}

LocalDiagSCMatrix::LocalDiagSCMatrix(LocalSCDimension*a):
  d(a)
{
  resize(a->n());
}

LocalDiagSCMatrix::LocalDiagSCMatrix(StateIn&s):
  DiagSCMatrix(s),
  SavableState(s,class_desc_)
{
  d.restore_state(s);
  block.restore_state(s);
}

LocalDiagSCMatrix::LocalDiagSCMatrix(KeyVal&keyval)
{
  d = keyval.describedclassvalue("dim");
  d.require_nonnull();
  block = new SCMatrixDiagBlock(0,d->n());
  for (int i=0; i<n(); i++) {
      set_element(i,keyval.doublevalue("data",i));
    }
}

void
LocalDiagSCMatrix::save_data_state(StateOut&s)
{
  DiagSCMatrix::save_data_state(s);
  d.save_state(s);
  block.save_state(s);
}

void *
LocalDiagSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { DiagSCMatrix::_castdown(cd) };
  return do_castdowns(casts,cd);
}

LocalDiagSCMatrix::~LocalDiagSCMatrix()
{
}

void
LocalDiagSCMatrix::resize(int n)
{
  block = new SCMatrixDiagBlock(0,n);
}

RefSCDimension
LocalDiagSCMatrix::dim()
{
  return d;
}

double
LocalDiagSCMatrix::get_element(int i)
{
  return block->data[i];
}

void
LocalDiagSCMatrix::set_element(int i,double a)
{
  block->data[i] = a;
}

void
LocalDiagSCMatrix::accumulate(DiagSCMatrix*a)
{
  // make sure that the argument is of the correct type
  LocalDiagSCMatrix* la
    = LocalDiagSCMatrix::require_castdown(a,"LocalDiagSCMatrix::accumulate");

  // make sure that the dimensions match
  if (this->dim() != la->dim()) {
      fprintf(stderr,"LocalDiagSCMatrix::"
              "accumulate(SCMatrix*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = n();
  for (int i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

double
LocalDiagSCMatrix::invert_this()
{
  double det = 1.0;
  int nelem = n();
  double*data = block->data;
  for (int i=0; i<nelem; i++) {
      data[i] = 1.0/data[i];
      det *= data[i];
    }
  return det;
}

void
LocalDiagSCMatrix::element_op(RefSCDiagElementOp& op)
{
  op->process(block);
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
LocalDiagSCMatrix::print(const char *title, ostream& os, int prec)
{
  int i;
  int lwidth,width;
  double max=this->maxabs();

  max=(max==0.0)?1.0:log10(max);
  if(max < 0.0) max=1.0;

  lwidth = prec+5+(int) max; width = 75/lwidth;

  os.setf(ios::fixed,ios::floatfield); os.precision(prec);
  os.setf(ios::right,ios::adjustfield);

  if(title) os << "\n" << title << "\n";
  else os << "\n";

  if(n()==0) { os << " empty matrix\n"; return; }

  for (i=0; i<n(); i++) {
      os.width(5); os << i+1;
      os.width(lwidth); os << block->data[i];
      os << "\n";
    }
  os << "\n";

  os.flush();
}

/////////////////////////////////////////////////////////////////////////////
// LocalSCVector member functions

#define CLASSNAME LocalSCVector
#define PARENTS public SCVector
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

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
  SavableState(s,class_desc_)
{
  d.restore_state(s);
  block.restore_state(s);
}

LocalSCVector::LocalSCVector(KeyVal&keyval)
{
  d = keyval.describedclassvalue("dim");
  d.require_nonnull();
  block = new SCVectorSimpleBlock(0,d->n());
  for (int i=0; i<n(); i++) {
      set_element(i,keyval.doublevalue("data",i,i));
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
  void* casts[] =  { SCVector::_castdown(cd) };
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
  return block->data[i];
}

void
LocalSCVector::set_element(int i,double a)
{
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
  if (this->dim() != a->rowdim()
      || a->coldim() != b->dim()) {
      fprintf(stderr,"LocalSCVector::"
              "accumulate_product(SCMatrix*a,SCVector*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  cmat_mxm(la->rows, 0,
           &(lb->block->data), 0,
           &(block->data), 0,
           n(), la->ncol(), 1,
           1);
}

void
LocalSCVector::accumulate(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::accumulate");

  // make sure that the dimensions match
  if (this->dim() != la->dim()) {
      fprintf(stderr,"LocalSCVector::"
              "accumulate(SCVector*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = n();
  for (int i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

double
LocalSCVector::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::scalar_product");

  // make sure that the dimensions match
  if (this->dim() != la->dim()) {
      fprintf(stderr,"LocalSCVector::"
              "scale_product(SCVector*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = n();
  double result = 0.0;
  for (int i=0; i<nelem; i++) result += block->data[i] * la->block->data[i];
  return result;
}

void
LocalSCVector::element_op(RefSCVectorElementOp& op)
{
  op->process(block);
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
LocalSCVector::print(const char *title, ostream& os, int prec)
{
  int i;
  int lwidth,width;
  double max=this->maxabs();

  max=(max==0.0)?1.0:log10(max);
  if(max < 0.0) max=1.0;

  lwidth = prec+5+(int) max; width = 75/lwidth;

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
