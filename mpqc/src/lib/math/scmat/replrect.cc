
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/repl.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// ReplSCMatrix member functions

#define CLASSNAME ReplSCMatrix
#define PARENTS public SCMatrix
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
ReplSCMatrix::_castdown(const ClassDesc*cd)
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

ReplSCMatrix::ReplSCMatrix(const RefSCDimension&a,const RefSCDimension&b,
                           ReplSCMatrixKit*k):
  SCMatrix(a,b,k)
{
  int nr = a->n();
  int nc = b->n();

  matrix = new double[nr*nc];

  rows = init_rect_rows(matrix,nr,nc);

  init_blocklist();
}

void
ReplSCMatrix::before_elemop()
{
  // zero out the blocks not in my block list
  int i, j, index;
  int nc = d2->n();
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  for (i=0, index=0; i<d1->blocks()->nblock(); i++) {
      for (j=0; j<d2->blocks()->nblock(); j++, index++) {
          if (index%nproc == me) continue;
          for (int ii=d1->blocks()->start(i);
               ii<d1->blocks()->fence(i); ii++) {
              for (int jj=d2->blocks()->start(j);
                   jj<d2->blocks()->fence(j); jj++) {
                  matrix[ii*nc + jj] = 0.0;
                }
            }
        }
    }
}

void
ReplSCMatrix::after_elemop()
{
  messagegrp()->sum(matrix, d1->n()*d2->n());
}

void
ReplSCMatrix::init_blocklist()
{
  int i, j, index;
  int nc = d2->n();
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  for (i=0, index=0; i<d1->blocks()->nblock(); i++) {
      for (j=0; j<d2->blocks()->nblock(); j++, index++) {
          if (index%nproc == me) {
              blocklist->insert(
                  new SCMatrixRectSubBlock(d1->blocks()->start(i),
                                           d1->blocks()->fence(i),
                                           nc,
                                           d2->blocks()->start(j),
                                           d2->blocks()->fence(j),
                                           matrix));
            }
        }
    }
}

ReplSCMatrix::~ReplSCMatrix()
{
  if (matrix) delete[] matrix;
  if (rows) delete[] rows;
}

int
ReplSCMatrix::compute_offset(int i,int j)
{
  if (i<0 || j<0 || i>=d1->n() || j>=d2->n()) {
      fprintf(stderr,"ReplSCMatrix: index out of bounds\n");
      abort();
    }
  return i*(d2->n()) + j;
}

double
ReplSCMatrix::get_element(int i,int j)
{
  int off = compute_offset(i,j);
  return matrix[off];
}

void
ReplSCMatrix::set_element(int i,int j,double a)
{
  int off = compute_offset(i,j);
  matrix[off] = a;
}

void
ReplSCMatrix::accumulate_element(int i,int j,double a)
{
  int off = compute_offset(i,j);
  matrix[off] += a;
}

SCMatrix *
ReplSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
    fprintf(stderr,"ReplSCMatrix::get_subblock: trying to get too big a"
            "subblock (%d,%d) from (%d,%d)\n",nsrow,nscol,nrow(),ncol());
    abort();
  }
  
  RefSCDimension dnrow = new SCDimension(nsrow);
  RefSCDimension dncol = new SCDimension(nscol);

  SCMatrix * sb = kit()->matrix(dnrow,dncol);
  sb->assign(0.0);

  ReplSCMatrix *lsb = ReplSCMatrix::require_castdown(sb,
                                      "ReplSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      lsb->rows[i][j] = rows[i+br][j+bc];
      
  return sb;
}

void
ReplSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                              int source_br, int source_bc)
{
  ReplSCMatrix *lsb = ReplSCMatrix::require_castdown(sb,
                                      "ReplSCMatrix::assign_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
    fprintf(stderr,"ReplSCMatrix::assign_subblock: trying to assign too big a"
            " subblock (%d,%d) to (%d,%d)\n",nsrow,nscol,nrow(),ncol());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      rows[i+br][j+bc] = lsb->rows[source_br + i][source_bc + j];
}

void
ReplSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                                  int source_br, int source_bc)
{
  ReplSCMatrix *lsb = ReplSCMatrix::require_castdown(sb,
                                      "ReplSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
    fprintf(stderr,"ReplSCMatrix::accumulate_subblock: trying to accumulate"
            " too big a subblock (%d,%d) to (%d,%d)\n",
            nsrow,nscol,nrow(),ncol());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      rows[i+br][j+bc] += lsb->rows[source_br + i][source_bc + j];
}

SCVector *
ReplSCMatrix::get_row(int i)
{
  if (i >= nrow()) {
    fprintf(stderr,"ReplSCMatrix::get_row: trying to get invalid row"
            "%d max %d\n",i,nrow());
    abort();
  }
  
  SCVector * v = kit()->vector(coldim());

  ReplSCVector *lv = ReplSCVector::require_castdown(v,
                                               "ReplSCMatrix::get_row");

  for (int j=0; j < ncol(); j++)
    lv->set_element(j,rows[i][j]);
      
  return v;
}

void
ReplSCMatrix::assign_row(SCVector *v, int i)
{
  if (i >= nrow()) {
    fprintf(stderr,"ReplSCMatrix::assign_row: trying to assign invalid row"
            "%d max %d\n",i,nrow());
    abort();
  }
  
  if (v->n() != ncol()) {
    fprintf(stderr,"ReplSCMatrix::assign_row: vector is wrong size"
            "is %d, should be %d\n",v->n(),ncol());
    abort();
  }
  
  ReplSCVector *lv = ReplSCVector::require_castdown(v,
                                               "ReplSCMatrix::assign_row");

  for (int j=0; j < ncol(); j++)
    rows[i][j] = lv->get_element(j);
}

void
ReplSCMatrix::accumulate_row(SCVector *v, int i)
{
  if (i >= nrow()) {
    fprintf(stderr,"ReplSCMatrix::accumulate_row: trying to assign invalid "
                   "row %d max %d\n",i,nrow());
    abort();
  }
  
  if (v->n() != ncol()) {
    fprintf(stderr,"ReplSCMatrix::accumulate_row: vector is wrong size"
            "is %d, should be %d\n",v->n(),ncol());
    abort();
  }
  
  ReplSCVector *lv = ReplSCVector::require_castdown(v,
                                            "ReplSCMatrix::accumulate_row");

  for (int j=0; j < ncol(); j++)
    rows[i][j] += lv->get_element(j);
}

SCVector *
ReplSCMatrix::get_column(int i)
{
  if (i >= ncol()) {
    fprintf(stderr,"ReplSCMatrix::get_column: trying to get invalid column"
            "%d max %d\n",i,ncol());
    abort();
  }
  
  SCVector * v = kit()->vector(rowdim());

  ReplSCVector *lv = ReplSCVector::require_castdown(v,
                                               "ReplSCMatrix::get_column");

  for (int j=0; j < nrow(); j++)
    lv->set_element(j,rows[j][i]);
      
  return v;
}

void
ReplSCMatrix::assign_column(SCVector *v, int i)
{
  if (i >= ncol()) {
    fprintf(stderr,"ReplSCMatrix::assign_column: trying to assign invalid "
            "column %d max %d\n",i,ncol());
    abort();
  }
  
  if (v->n() != nrow()) {
    fprintf(stderr,"ReplSCMatrix::assign_column: vector is wrong size"
            "is %d, should be %d\n",v->n(),nrow());
    abort();
  }
  
  ReplSCVector *lv = ReplSCVector::require_castdown(v,
                                               "ReplSCMatrix::assign_column");

  for (int j=0; j < nrow(); j++)
    rows[j][i] = lv->get_element(j);
}

void
ReplSCMatrix::accumulate_column(SCVector *v, int i)
{
  if (i >= ncol()) {
    fprintf(stderr,"ReplSCMatrix::accumulate_column: trying to assign invalid"
            " column %d max %d\n",i,ncol());
    abort();
  }
  
  if (v->n() != nrow()) {
    fprintf(stderr,"ReplSCMatrix::accumulate_column: vector is wrong size"
            "is %d, should be %d\n",v->n(),nrow());
    abort();
  }
  
  ReplSCVector *lv = ReplSCVector::require_castdown(v,
                                     "ReplSCMatrix::accumulate_column");

  for (int j=0; j < nrow(); j++)
    rows[j][i] += lv->get_element(j);
}

void
ReplSCMatrix::assign(double a)
{
  int n = d1->n() * d2->n();
  for (int i=0; i<n; i++) matrix[i] = a;
}

void
ReplSCMatrix::accumulate_product(SCMatrix*a,SCMatrix*b)
{
  const char* name = "ReplSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  ReplSCMatrix* la = ReplSCMatrix::require_castdown(a,name);
  ReplSCMatrix* lb = ReplSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->coldim()) ||
      !la->coldim()->equiv(lb->rowdim())) {
      fprintf(stderr,"ReplSCMatrix::"
              "accumulate_product(SCMatrix*a,SCMatrix*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  cmat_transpose_matrix(lb->rows, la->ncol(), this->ncol());
  double** btrans;
  btrans = new double*[this->ncol()];
  btrans[0] = lb->rows[0];
  cmat_matrix_pointers(btrans,btrans[0],this->ncol(),la->ncol());

  RefSCElementOp op = new SCElementDot(la->rows, btrans, la->ncol());
  element_op(op);

  cmat_transpose_matrix(btrans,this->ncol(),la->ncol());
  delete[] btrans;
}

// does the outer product a x b.  this must have rowdim() == a->dim() and
// coldim() == b->dim()
void
ReplSCMatrix::accumulate_outer_product(SCVector*a,SCVector*b)
{
  const char* name = "ReplSCMatrix::accumulate_outer_product";
  // make sure that the arguments are of the correct type
  ReplSCVector* la = ReplSCVector::require_castdown(a,name);
  ReplSCVector* lb = ReplSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(lb->dim())) {
      fprintf(stderr,"ReplSCMatrix::"
              "accumulate_outer_product(SCVector*a,SCVector*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nr = a->n();
  int nc = b->n();
  int i, j;
  double* adat = la->vector;
  double* bdat = lb->vector;
  double** thisdat = rows;
  for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
          thisdat[i][j] += adat[i] * bdat[j];
        }
    }
}

void
ReplSCMatrix::accumulate_product(SCMatrix*a,SymmSCMatrix*b)
{
  const char* name = "ReplSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  ReplSCMatrix* la = ReplSCMatrix::require_castdown(a,name);
  ReplSymmSCMatrix* lb = ReplSymmSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
      fprintf(stderr,"ReplSCMatrix::"
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
ReplSCMatrix::accumulate_product(SCMatrix*a,DiagSCMatrix*b)
{
  const char* name = "ReplSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  ReplSCMatrix* la = ReplSCMatrix::require_castdown(a,name);
  ReplDiagSCMatrix* lb = ReplDiagSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
      fprintf(stderr,"ReplSCMatrix::"
              "accumulate_product(SCMatrix*a,DiagSCMatrix*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  double **cd = rows;
  double **ad = la->rows;
  double *bd = lb->matrix;
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
ReplSCMatrix::accumulate(SCMatrix*a)
{
  // make sure that the arguments is of the correct type
  ReplSCMatrix* la
    = ReplSCMatrix::require_castdown(a,"ReplSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(la->coldim())) {
      fprintf(stderr,"ReplSCMatrix::accumulate(SCMatrix*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = this->ncol() * this->nrow();
  int i;
  for (i=0; i<nelem; i++) matrix[i] += la->matrix[i];
}

void
ReplSCMatrix::accumulate(SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  ReplSymmSCMatrix* la
    = ReplSymmSCMatrix::require_castdown(a,"ReplSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
      fprintf(stderr,"ReplSCMatrix::accumulate(SymmSCMatrix*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int n = this->ncol();
  double *dat = la->matrix;
  int i, j;
  for (i=0; i<n; i++) {
      for (j=0; j<i; j++) {
          double tmp = *dat;
          matrix[i*n+j] += tmp;
          matrix[j*n+i] += tmp;
          dat++;
        }
      matrix[i*n+i] += *dat++;
    }
}

void
ReplSCMatrix::accumulate(DiagSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  ReplDiagSCMatrix* la
    = ReplDiagSCMatrix::require_castdown(a,"ReplSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
      fprintf(stderr,"ReplSCMatrix::accumulate(DiagSCMatrix*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int n = this->ncol();
  double *dat = la->matrix;
  int i;
  for (i=0; i<n; i++) {
      matrix[i*n+i] += *dat++;
    }
}

void
ReplSCMatrix::accumulate(SCVector*a)
{
  // make sure that the arguments is of the correct type
  ReplSCVector* la
    = ReplSCVector::require_castdown(a,"ReplSCVector::accumulate");

  // make sure that the dimensions match
  if (!((rowdim()->equiv(la->dim()) && coldim()->n() == 1)
        || (coldim()->equiv(la->dim()) && rowdim()->n() == 1))) {
      fprintf(stderr,"ReplSCMatrix::accumulate(SCVector*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int n = this->ncol();
  int i;
  double *dat = la->vector;
  for (i=0; i<n; i++) {
      matrix[i*n+i] += dat[i];
    }
}

void
ReplSCMatrix::transpose_this()
{
  cmat_transpose_matrix(rows,nrow(),ncol());
  delete[] rows;
  rows = new double*[ncol()];
  cmat_matrix_pointers(rows,matrix,ncol(),nrow());
  RefSCDimension tmp = d1;
  d1 = d2;
  d2 = tmp;
  init_blocklist();
}

double
ReplSCMatrix::invert_this()
{
  if (nrow() != ncol()) {
      fprintf(stderr,"ReplSCMatrix::invert_this: matrix is not square\n");
      abort();
    }
  return cmat_invert(rows,0,nrow());
}

void
ReplSCMatrix::gen_invert_this()
{
  fprintf(stderr,"ReplSCMatrix::gen_invert_this: SVD not implemented yet");
  abort();
}

double
ReplSCMatrix::determ_this()
{
  if (nrow() != ncol()) {
    fprintf(stderr,"ReplSCMatrix::determ_this: matrix is not square\n");
    abort();
  }
  return cmat_determ(rows,0,nrow());
}

double
ReplSCMatrix::trace()
{
  if (nrow() != ncol()) {
    fprintf(stderr,"ReplSCMatrix::trace: matrix is not square\n");
    abort();
  }
  double ret=0;
  int i;
  for (i=0; i < nrow(); i++)
    ret += rows[i][i];
  return ret;
}

double
ReplSCMatrix::solve_this(SCVector*v)
{
  ReplSCVector* lv =
    ReplSCVector::require_castdown(v,"ReplSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lv->dim())) {
      fprintf(stderr,"ReplSCMatrix::solve_this(SCVector*v):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  return cmat_solve_lin(rows,0,lv->vector,nrow());
}

void
ReplSCMatrix::schmidt_orthog(SymmSCMatrix *S, int nc)
{
  ReplSymmSCMatrix* lS =
    ReplSymmSCMatrix::require_castdown(S,"ReplSCMatrix::schmidt_orthog");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lS->dim())) {
      fprintf(stderr,"ReplSCMatrix::schmidt_orthog():\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  cmat_schmidt(rows,lS->matrix,nrow(),nc);
}

void
ReplSCMatrix::element_op(const RefSCElementOp& op)
{
  if (op->has_side_effects()) before_elemop();
  SCMatrixBlockListIter i;
  for (i = blocklist->begin(); i != blocklist->end(); i++) {
      op->process(i.block());
    }
  if (op->has_side_effects()) after_elemop();
}

void
ReplSCMatrix::element_op(const RefSCElementOp2& op,
                          SCMatrix* m)
{
  ReplSCMatrix *lm
      = ReplSCMatrix::require_castdown(m,"ReplSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim())) {
      fprintf(stderr,"ReplSCMatrix: bad element_op\n");
      abort();
    }
  if (op->has_side_effects()) before_elemop();
  if (op->has_side_effects_in_arg()) lm->before_elemop();
  SCMatrixBlockListIter i, j;
  for (i = blocklist->begin(), j = lm->blocklist->begin();
       i != blocklist->end();
       i++, j++) {
      op->process(i.block(), j.block());
    }
  if (op->has_side_effects()) after_elemop();
  if (op->has_side_effects_in_arg()) lm->after_elemop();
}

void
ReplSCMatrix::element_op(const RefSCElementOp3& op,
                          SCMatrix* m,SCMatrix* n)
{
  ReplSCMatrix *lm
      = ReplSCMatrix::require_castdown(m,"ReplSCMatrix::element_op");
  ReplSCMatrix *ln
      = ReplSCMatrix::require_castdown(n,"ReplSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim()) ||
      !rowdim()->equiv(ln->rowdim()) || !coldim()->equiv(ln->coldim())) {
      fprintf(stderr,"ReplSCMatrix: bad element_op\n");
      abort();
    }
  if (op->has_side_effects()) before_elemop();
  if (op->has_side_effects_in_arg1()) lm->before_elemop();
  if (op->has_side_effects_in_arg2()) ln->before_elemop();
  SCMatrixBlockListIter i, j, k;
  for (i = blocklist->begin(),
           j = lm->blocklist->begin(),
           k = ln->blocklist->begin();
       i != blocklist->end();
       i++, j++, k++) {
      op->process(i.block(), j.block(), k.block());
    }
  if (op->has_side_effects()) after_elemop();
  if (op->has_side_effects_in_arg1()) lm->after_elemop();
  if (op->has_side_effects_in_arg2()) ln->after_elemop();
}

// from Ed Seidl at the NIH
void
ReplSCMatrix::print(const char *title, ostream& os, int prec)
{
  if (messagegrp()->me() != 0) return;

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
ReplSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new ReplSCMatrixListSubblockIter(access, blocklist,
                                          messagegrp(),
                                          matrix, d1->n()*d2->n());
}

RefSCMatrixSubblockIter
ReplSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      cerr << "ReplSCMatrix::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  RefSCMatrixBlockList allblocklist = new SCMatrixBlockList();
  allblocklist->insert(new SCMatrixRectSubBlock(0, d1->n(), d1->n(),
                                                0, d2->n(), matrix));
  return new ReplSCMatrixListSubblockIter(access, allblocklist,
                                          messagegrp(),
                                          matrix, d1->n()*d2->n());
}

RefReplSCMatrixKit
ReplSCMatrix::skit()
{
  return ReplSCMatrixKit::castdown(kit().pointer());
}

RefMessageGrp
ReplSCMatrix::messagegrp()
{
  return skit()->messagegrp();
}
