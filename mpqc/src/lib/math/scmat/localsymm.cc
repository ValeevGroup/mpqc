
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// LocalSymmSCMatrix member functions

#define CLASSNAME LocalSymmSCMatrix
#define PARENTS public SymmSCMatrix
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LocalSymmSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SymmSCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

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
  d(a),
  rows(0)
{
  resize(a->n());
}

LocalSymmSCMatrix::LocalSymmSCMatrix(StateIn&s):
  SymmSCMatrix(s)
{
  d.restore_state(s);
  block.restore_state(s);
  // if (rows) delete[] rows;
  rows = init_symm_rows(block->data,d->n());
}

LocalSymmSCMatrix::LocalSymmSCMatrix(const RefKeyVal&keyval)
{
  d = keyval->describedclassvalue("dim");
  d.require_nonnull();
  block = new SCMatrixLTriBlock(0,d->n());
  rows = init_symm_rows(block->data,d->n());
  for (int i=0; i<n(); i++) {
      for (int j=0; j<=i; j++) {
          set_element(i,j,keyval->doublevalue("data",i,j));
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
LocalSymmSCMatrix::accumulate_element(int i,int j,double a)
{
  block->data[compute_offset(i,j)] += a;
}

SCMatrix *
LocalSymmSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    fprintf(stderr,"LocalSymmSCMatrix::get_subblock: trying to get too big a "
            "subblock (%d,%d) from (%d,%d)\n",nsrow,nscol,n(),n());
    abort();
  }
  
  RefLocalSCDimension dnrow = new LocalSCDimension(nsrow);
  RefLocalSCDimension dncol = new LocalSCDimension(nscol);

  SCMatrix * sb = dnrow->create_matrix(dncol.pointer());
  sb->assign(0.0);

  LocalSCMatrix *lsb = LocalSCMatrix::require_castdown(sb,
                                      "LocalSymmSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      lsb->rows[i][j] = get_element(i+br,j+bc);
      
  return sb;
}

SymmSCMatrix *
LocalSymmSCMatrix::get_subblock(int br, int er)
{
  int nsrow = er-br+1;

  if (nsrow > n()) {
    fprintf(stderr,"LocalSymmSCMatrix::get_subblock: trying to get too big a "
            "subblock (%d,%d) from (%d,%d)\n",nsrow,nsrow,n(),n());
    abort();
  }
  
  RefLocalSCDimension dnrow = new LocalSCDimension(nsrow);

  SymmSCMatrix * sb = dnrow->create_symmmatrix();
  sb->assign(0.0);

  LocalSymmSCMatrix *lsb = LocalSymmSCMatrix::require_castdown(sb,
                                      "LocalSymmSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      lsb->rows[i][j] = get_element(i+br,j+br);
      
  return sb;
}

void
LocalSymmSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  LocalSCMatrix *lsb = LocalSCMatrix::require_castdown(sb,
                                      "LocalSCMatrix::assign_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    fprintf(stderr,"LocalSymmSCMatrix::assign_subblock: "
            "trying to assign too big a "
            "subblock (%d,%d) to (%d,%d)\n",nsrow,nscol,n(),n());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      set_element(i+br,j+bc,lsb->rows[i][j]);
}

void
LocalSymmSCMatrix::assign_subblock(SymmSCMatrix*sb, int br, int er)
{
  LocalSymmSCMatrix *lsb = LocalSymmSCMatrix::require_castdown(sb,
                                      "LocalSymmSCMatrix::assign_subblock");

  int nsrow = er-br+1;

  if (nsrow > n()) {
    fprintf(stderr,"LocalSymmSCMatrix::assign_subblock: "
            "trying to assign too big a "
            "subblock (%d,%d) to (%d,%d)\n",nsrow,nsrow,n(),n());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      set_element(i+br,j+br,lsb->rows[i][j]);
}

void
LocalSymmSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  LocalSCMatrix *lsb = LocalSCMatrix::require_castdown(sb,
                                  "LocalSymmSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    fprintf(stderr,"LocalSymmSCMatrix::accumulate_subblock: trying to "
            "accumulate too big a subblock (%d,%d) to (%d,%d)\n",
            nsrow,nscol,n(),n());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      set_element(i+br,j+br,get_element(i+br,j+br)+lsb->rows[i][j]);
}

void
LocalSymmSCMatrix::accumulate_subblock(SymmSCMatrix*sb, int br, int er)
{
  LocalSCMatrix *lsb = LocalSCMatrix::require_castdown(sb,
                                  "LocalSymmSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;

  if (nsrow > n()) {
    fprintf(stderr,"LocalSymmSCMatrix::accumulate_subblock: trying to "
            "accumulate too big a subblock (%d,%d) to (%d,%d)\n",
            nsrow,nsrow,n(),n());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      set_element(i+br,j+br,get_element(i+br,j+br)+lsb->rows[i][j]);
}

SCVector *
LocalSymmSCMatrix::get_row(int i)
{
  if (i >= n()) {
    fprintf(stderr,"LocalSymmSCMatrix::get_row: trying to get invalid row"
            "%d max %d\n",i,n());
    abort();
  }
  
  SCVector * v = dim()->create_vector();

  LocalSCVector *lv = LocalSCVector::require_castdown(v,
                                               "LocalSymmSCMatrix::get_row");

  for (int j=0; j < n(); j++)
    lv->set_element(j,get_element(i,j));
      
  return v;
}

void
LocalSymmSCMatrix::assign_row(SCVector *v, int i)
{
  if (i >= n()) {
    fprintf(stderr,"LocalSymmSCMatrix::assign_row: trying to assign invalid "
            "row %d max %d\n",i,n());
    abort();
  }
  
  if (v->n() != n()) {
    fprintf(stderr,"LocalSymmSCMatrix::assign_row: vector is wrong size"
            "is %d, should be %d\n",v->n(),n());
    abort();
  }
  
  LocalSCVector *lv = LocalSCVector::require_castdown(v,
                                          "LocalSymmSCMatrix::assign_row");

  for (int j=0; j < n(); j++)
    set_element(i,j,lv->get_element(j));
}

void
LocalSymmSCMatrix::accumulate_row(SCVector *v, int i)
{
  if (i >= n()) {
    fprintf(stderr,"LocalSymmSCMatrix::accumulate_row: trying to assign "
                   "invalid row %d max %d\n",i,n());
    abort();
  }
  
  if (v->n() != n()) {
    fprintf(stderr,"LocalSymmSCMatrix::accumulate_row: vector is wrong size"
            "is %d, should be %d\n",v->n(),n());
    abort();
  }
  
  LocalSCVector *lv = LocalSCVector::require_castdown(v,
                                        "LocalSymmSCMatrix::accumulate_row");

  for (int j=0; j < n(); j++)
    set_element(i,j,get_element(i,j)+lv->get_element(j));
}

void
LocalSymmSCMatrix::accumulate(SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  LocalSymmSCMatrix* la
    = LocalSymmSCMatrix::require_castdown(a,"LocalSymmSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!(this->dim() == la->dim())) {
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

double
LocalSymmSCMatrix::determ_this()
{
  return cmat_determ(rows,1,n());
}

double
LocalSymmSCMatrix::trace()
{
  double ret=0;
  for (int i=0; i < n(); i++) ret += rows[i][i];
  return ret;
}

double
LocalSymmSCMatrix::solve_this(SCVector*v)
{
  LocalSCVector* lv =
    LocalSCVector::require_castdown(v,"LocalSymmSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!(this->dim() == v->dim())) {
      fprintf(stderr,"LocalSymmSCMatrix::solve_this(SCVector*v):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  return cmat_solve_lin(rows,1,lv->block->data,n());
}

void
LocalSymmSCMatrix::gen_invert_this()
{
  double *evals = new double[n()];
  double **evecs = cmat_new_square_matrix(n());
  
  cmat_diag(rows,evals,evecs,n(),1,1.0e-15);

  for (int i=0; i < n(); i++) {
    if (fabs(evals[i]) > 1.0e-8)
      evals[i] = 1.0/evals[i];
    else
      evals[i] = 0;
  }

  cmat_transform_diagonal_matrix(rows, n(), evals, n(), evecs, 0);
  
  delete[] evals;
  cmat_delete_matrix(evecs);  
}

void
LocalSymmSCMatrix::diagonalize(DiagSCMatrix*a,SCMatrix*b)
{
  const char* name = "LocalSymmSCMatrix::diagonalize";
  // make sure that the arguments is of the correct type
  LocalDiagSCMatrix* la = LocalDiagSCMatrix::require_castdown(a,name);
  LocalSCMatrix* lb = LocalSCMatrix::require_castdown(b,name);

  if (   (la&&(!(la->dim() == dim())))
      || (lb&&(!(lb->coldim() == dim()) || !(lb->rowdim() == dim())))) {
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

  if (!(la->rowdim() == dim())) {
      fprintf(stderr,"LocalSymmSCMatrix::"
              "accumulate_symmetric_product(SCMatrix*a): bad dim");
      abort();
    }

  cmat_symmetric_mxm(rows,n(),la->rows,la->ncol(),1);
}

// computes this += a + a.t
void
LocalSymmSCMatrix::accumulate_symmetric_sum(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  LocalSCMatrix* la
    = LocalSCMatrix::require_castdown(a,"LocalSymmSCMatrix::"
                                          "accumulate_symmetric_sum");

  if (!(la->rowdim() == dim()) || !(la->coldim() == dim())) {
      fprintf(stderr,"LocalSymmSCMatrix::"
              "accumulate_symmetric_sum(SCMatrix*a): bad dim");
      abort();
    }

  int n = dim().n();
  double** tdat = this->rows;
  double** adat = la->rows;
  for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
          tdat[i][j] += adat[i][j] + adat[j][i];
        }
    }
}

void
LocalSymmSCMatrix::accumulate_symmetric_outer_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSymmSCMatrix::"
                                      "accumulate_symmetric_outer_product");

  if (!(la->dim() == dim())) {
      fprintf(stderr,"LocalSymmSCMatrix::"
              "accumulate_symmetric_outer_product(SCMatrix*a): bad dim");
      abort();
    }

  int n = dim().n();
  double** tdat = this->rows;
  double* adat = la->block->data;
  for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
          tdat[i][j] += adat[i]*adat[j];
        }
    }
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

// this += a * b * transpose(a)
void
LocalSymmSCMatrix::accumulate_transform(SCMatrix*a,DiagSCMatrix*b)
{
  // do the necessary castdowns
  LocalSCMatrix*la
    = LocalSCMatrix::require_castdown(a,"%s::accumulate_transform",
                                      class_name());
  LocalDiagSCMatrix*lb
    = LocalDiagSCMatrix::require_castdown(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (n() != la->nrow() || la->ncol() != lb->n()) {
      fprintf(stderr,"LocalSymmSCMatrix::accumulate_transform: bad dim\n");
      abort();
    }

  cmat_transform_diagonal_matrix(rows,n(),lb->block->data,lb->n(),la->rows,1);
}

double
LocalSymmSCMatrix::scalar_product(SCVector*a)
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

  int nelem = n();
  double* adat = la->block->data;
  double result = 0.0;
  for (int i=0; i<nelem; i++) {
      for (int j=0; j<i; j++) {
          result += 2.0 * rows[i][j] * adat[i] * adat[j];
        }
      result += rows[i][i] * adat[i] * adat[i];
    }
  return result;
}

void
LocalSymmSCMatrix::element_op(const RefSCElementOp& op)
{
  op->process(block.pointer());
}

void
LocalSymmSCMatrix::element_op(const RefSCElementOp2& op,
                              SymmSCMatrix* m)
{
  LocalSymmSCMatrix *lm
      = LocalSymmSCMatrix::require_castdown(m,"LocalSymSCMatrix::element_op");
  if (!lm || d != lm->d) {
      fprintf(stderr,"LocalSymmSCMatrix: bad element_op\n");
      abort();
    }
  op->process(block.pointer(), lm->block.pointer());
}

void
LocalSymmSCMatrix::element_op(const RefSCElementOp3& op,
                              SymmSCMatrix* m,SymmSCMatrix* n)
{
  LocalSymmSCMatrix *lm
      = LocalSymmSCMatrix::require_castdown(m,"LocalSymSCMatrix::element_op");
  LocalSymmSCMatrix *ln
      = LocalSymmSCMatrix::require_castdown(n,"LocalSymSCMatrix::element_op");
  if (!lm || !ln || d != lm->d || d != ln->d) {
      fprintf(stderr,"LocalSymmSCMatrix: bad element_op\n");
      abort();
    }
  op->process(block.pointer(), lm->block.pointer(), ln->block.pointer());
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
LocalSymmSCMatrix::print(const char *title, ostream& os, int prec)
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

  if(n()==0) { os << " empty matrix\n"; return; }

  for(ii=jj=0;;) {
    ii++; jj++; kk=width*jj;
    nn=(n()>kk)?kk:n();

 // print column indices
    for(i=ii; i <= nn; i++) { os.width(lwidth); os << i; }
    os << "\n";

 // print the rows
    for(i=ii-1; i < n() ; i++) {
      os.width(5); os << i+1;
      for(j=ii-1; j<nn && j<=i; j++) { os.width(lwidth); os << rows[i][j]; }
      os << "\n";
      }
    os << "\n";

    if(n()<=kk) { os.flush(); return; }
    ii=kk;
    }
}
