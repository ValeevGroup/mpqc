
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

extern "C" {
    int sing_(double *q, int *lq, int *iq, double *s, double *p,
              int *lp, int *ip, double *a, int *la, int *m, int *n,
              double *w);
};

/////////////////////////////////////////////////////////////////////////////
// LocalSCMatrix member functions

#define CLASSNAME LocalSCMatrix
#define PARENTS public SCMatrix
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LocalSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

static double **
init_rect_rows(double *data, int ni,int nj)
{
  double** r = new double*[ni];
  int i;
  for (i=0; i<ni; i++) r[i] = &data[i*nj];
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
  SCMatrix(s)
{
  d1.restore_state(s);
  d2.restore_state(s);
  block.restore_state(s);
  rows = init_rect_rows(block->data,d1->n(),d2->n());
}

LocalSCMatrix::LocalSCMatrix(const RefKeyVal&keyval)
{
  int i, j;
  d1 = keyval->describedclassvalue("rowdim");
  d2 = keyval->describedclassvalue("coldim");
  d1.require_nonnull();
  d2.require_nonnull();
  block = new SCMatrixRectBlock(0,d1->n(),0,d2->n());
  rows = init_rect_rows(block->data,d1->n(),d2->n());
  for (i=0; i<nrow(); i++) {
      for (j=0; j<ncol(); j++) {
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
LocalSCMatrix::accumulate_element(int i,int j,double a)
{
  int off = compute_offset(i,j);
  block->data[off] += a;
}

SCMatrix *
LocalSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
    fprintf(stderr,"LocalSCMatrix::get_subblock: trying to get too big a "
            "subblock (%d,%d) from (%d,%d)\n",nsrow,nscol,nrow(),ncol());
    abort();
  }
  
  RefLocalSCDimension dnrow = new LocalSCDimension(nsrow);
  RefLocalSCDimension dncol = new LocalSCDimension(nscol);

  SCMatrix * sb = dnrow->create_matrix(dncol.pointer());
  sb->assign(0.0);

  LocalSCMatrix *lsb = LocalSCMatrix::require_castdown(sb,
                                      "LocalSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      lsb->rows[i][j] = rows[i+br][j+bc];
      
  return sb;
}

void
LocalSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                               int source_br, int source_bc)
{
  LocalSCMatrix *lsb = LocalSCMatrix::require_castdown(sb,
                                      "LocalSCMatrix::assign_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
    fprintf(stderr,"LocalSCMatrix::assign_subblock: trying to assign too big a"
            " subblock (%d,%d) to (%d,%d)\n",nsrow,nscol,nrow(),ncol());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      rows[i+br][j+bc] = lsb->rows[source_br + i][source_bc + j];
}

void
LocalSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                                   int source_br, int source_bc)
{
  LocalSCMatrix *lsb = LocalSCMatrix::require_castdown(sb,
                                      "LocalSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
    fprintf(stderr,"LocalSCMatrix::accumulate_subblock: trying to accumulate"
            " too big a subblock (%d,%d) to (%d,%d)\n",
            nsrow,nscol,nrow(),ncol());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      rows[i+br][j+bc] += lsb->rows[source_br + i][source_bc + j]; 
}

SCVector *
LocalSCMatrix::get_row(int i)
{
  if (i >= nrow()) {
    fprintf(stderr,"LocalSCMatrix::get_row: trying to get invalid row"
            "%d max %d\n",i,nrow());
    abort();
  }
  
  SCVector * v = coldim()->create_vector();

  LocalSCVector *lv = LocalSCVector::require_castdown(v,
                                               "LocalSCMatrix::get_row");

  for (int j=0; j < ncol(); j++)
    lv->set_element(j,rows[i][j]);
      
  return v;
}

void
LocalSCMatrix::assign_row(SCVector *v, int i)
{
  if (i >= nrow()) {
    fprintf(stderr,"LocalSCMatrix::assign_row: trying to assign invalid row"
            "%d max %d\n",i,nrow());
    abort();
  }
  
  if (v->n() != ncol()) {
    fprintf(stderr,"LocalSCMatrix::assign_row: vector is wrong size"
            "is %d, should be %d\n",v->n(),ncol());
    abort();
  }
  
  LocalSCVector *lv = LocalSCVector::require_castdown(v,
                                               "LocalSCMatrix::assign_row");

  for (int j=0; j < ncol(); j++)
    rows[i][j] = lv->get_element(j);
}

void
LocalSCMatrix::accumulate_row(SCVector *v, int i)
{
  if (i >= nrow()) {
    fprintf(stderr,"LocalSCMatrix::accumulate_row: trying to assign invalid "
                   "row %d max %d\n",i,nrow());
    abort();
  }
  
  if (v->n() != ncol()) {
    fprintf(stderr,"LocalSCMatrix::accumulate_row: vector is wrong size"
            "is %d, should be %d\n",v->n(),ncol());
    abort();
  }
  
  LocalSCVector *lv = LocalSCVector::require_castdown(v,
                                            "LocalSCMatrix::accumulate_row");

  for (int j=0; j < ncol(); j++)
    rows[i][j] += lv->get_element(j);
}

SCVector *
LocalSCMatrix::get_column(int i)
{
  if (i >= ncol()) {
    fprintf(stderr,"LocalSCMatrix::get_column: trying to get invalid column"
            "%d max %d\n",i,ncol());
    abort();
  }
  
  SCVector * v = rowdim()->create_vector();

  LocalSCVector *lv = LocalSCVector::require_castdown(v,
                                               "LocalSCMatrix::get_column");

  for (int j=0; j < nrow(); j++)
    lv->set_element(j,rows[j][i]);
      
  return v;
}

void
LocalSCMatrix::assign_column(SCVector *v, int i)
{
  if (i >= ncol()) {
    fprintf(stderr,"LocalSCMatrix::assign_column: trying to assign invalid "
            "column %d max %d\n",i,ncol());
    abort();
  }
  
  if (v->n() != nrow()) {
    fprintf(stderr,"LocalSCMatrix::assign_column: vector is wrong size"
            "is %d, should be %d\n",v->n(),nrow());
    abort();
  }
  
  LocalSCVector *lv = LocalSCVector::require_castdown(v,
                                               "LocalSCMatrix::assign_column");

  for (int j=0; j < nrow(); j++)
    rows[j][i] = lv->get_element(j);
}

void
LocalSCMatrix::accumulate_column(SCVector *v, int i)
{
  if (i >= ncol()) {
    fprintf(stderr,"LocalSCMatrix::accumulate_column: trying to assign invalid"
            " column %d max %d\n",i,ncol());
    abort();
  }
  
  if (v->n() != nrow()) {
    fprintf(stderr,"LocalSCMatrix::accumulate_column: vector is wrong size"
            "is %d, should be %d\n",v->n(),nrow());
    abort();
  }
  
  LocalSCVector *lv = LocalSCVector::require_castdown(v,
                                     "LocalSCMatrix::accumulate_column");

  for (int j=0; j < nrow(); j++)
    rows[j][i] += lv->get_element(j);
}

void
LocalSCMatrix::assign(double a)
{
  int n = d1->n() * d2->n();
  double *data = block->data;
  for (int i=0; i<n; i++) data[i] = a;
}

void
LocalSCMatrix::accumulate_product(SCMatrix*a,SCMatrix*b)
{
  const char* name = "LocalSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = LocalSCMatrix::require_castdown(a,name);
  LocalSCMatrix* lb = LocalSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->coldim()) ||
      !la->coldim()->equiv(lb->rowdim())) {
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
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(lb->dim())) {
      fprintf(stderr,"LocalSCMatrix::"
              "accumulate_outer_product(SCVector*a,SCVector*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nr = a->n();
  int nc = b->n();
  int i, j;
  double* adat = la->block->data;
  double* bdat = lb->block->data;
  double** thisdat = rows;
  for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
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
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
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
  int i, j, k;
  for (i=0; i<ni; i++) {
      for (j=0; j<njk; j++) {
          for (k=0; k<=j; k++) {
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
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
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
  int i, j;
  for (i=0; i<ni; i++) {
      for (j=0; j<nj; j++) {
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
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(la->coldim())) {
      fprintf(stderr,"LocalSCMatrix::accumulate(SCMatrix*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = this->ncol() * this->nrow();
  int i;
  for (i=0; i<nelem; i++) block->data[i] += la->block->data[i];
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
  int i;
  for (i=0; i < nrow(); i++)
    ret += rows[i][i];
  return ret;
}

void
LocalSCMatrix::svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V)
{
  LocalSCMatrix* lU =
    LocalSCMatrix::require_castdown(U,"LocalSCMatrix::svd_this");
  LocalSCMatrix* lV =
    LocalSCMatrix::require_castdown(V,"LocalSCMatrix::svd_this");
  LocalDiagSCMatrix* lsigma =
    LocalDiagSCMatrix::require_castdown(sigma,"LocalSCMatrix::svd_this");

  RefSCDimension mdim = rowdim();
  RefSCDimension ndim = coldim();
  int m = mdim.n();
  int n = ndim.n();

  RefSCDimension pdim;
  if (m == n && m == sigma->dim().n())
    pdim = sigma->dim();
  else if (m<n)
    pdim = mdim;
  else
    pdim = ndim;

  int p = pdim.n();

  if (!mdim->equiv(lU->rowdim()) ||
      !mdim->equiv(lU->coldim()) ||
      !ndim->equiv(lV->rowdim()) ||
      !ndim->equiv(lV->coldim()) ||
      !pdim->equiv(sigma->dim())) {
      fprintf(stderr,"LocalSCMatrix: svd_this: dimension mismatch\n");
      abort();
    }

  // form a fortran style matrix for the SVD routines
  double *dA = new double[m*n];
  double *dU = new double[m*m];
  double *dV = new double[n*n];
  double *dsigma = new double[n];
  double *w = new double[(3*p-1>m)?(3*p-1):m];

  int i,j;
  for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {
          dA[i + j*m] = this->block->data[i*n + j];
        }
    }

  int three = 3;

  sing_(dU, &m, &three, dsigma, dV, &n, &three, dA, &m, &m, &n, w);

  for (i=0; i<m; i++) {
      for (j=0; j<m; j++) {
          lU->block->data[i*m + j] = dU[i + j*m];
        }
    }

  for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
          lV->block->data[i*n + j] = dV[i + j*n];
        }
    }

  for (i=0; i<p; i++) {
      lsigma->block->data[i] = dsigma[i];
    }

  delete[] dA;
  delete[] dU;
  delete[] dV;
  delete[] dsigma;
  delete[] w;
}

double
LocalSCMatrix::solve_this(SCVector*v)
{
  LocalSCVector* lv =
    LocalSCVector::require_castdown(v,"LocalSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lv->dim())) {
      fprintf(stderr,"LocalSCMatrix::solve_this(SCVector*v):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  return cmat_solve_lin(rows,0,lv->block->data,nrow());
}

void
LocalSCMatrix::schmidt_orthog(SymmSCMatrix *S, int nc)
{
  LocalSymmSCMatrix* lS =
    LocalSymmSCMatrix::require_castdown(S,"LocalSCMatrix::schmidt_orthog");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lS->dim())) {
      fprintf(stderr,"LocalSCMatrix::schmidt_orthog():\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  cmat_schmidt(rows,lS->block->data,nrow(),nc);
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

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim())) {
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

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim()) ||
      !rowdim()->equiv(ln->rowdim()) || !coldim()->equiv(ln->coldim())) {
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

RefSCMatrixSubblockIter
LocalSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  RefSCMatrixSubblockIter iter
      = new SCMatrixSimpleSubblockIter(access, block.pointer());
  return iter;
}

RefSCMatrixSubblockIter
LocalSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      cerr << "LocalSCMatrix::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  return local_blocks(access);
}
