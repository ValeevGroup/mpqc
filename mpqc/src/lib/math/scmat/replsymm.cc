
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/repl.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// ReplSymmSCMatrix member functions

#define CLASSNAME ReplSymmSCMatrix
#define PARENTS public SymmSCMatrix
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
ReplSymmSCMatrix::_castdown(const ClassDesc*cd)
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

ReplSymmSCMatrix::ReplSymmSCMatrix(const RefSCDimension&a,ReplSCMatrixKit*k):
  SymmSCMatrix(a,k),
  rows(0)
{
  int n = d->n();

  matrix = new double[n*(n+1)>>1];
  rows = init_symm_rows(matrix,n);

  init_blocklist();
}

void
ReplSymmSCMatrix::before_elemop()
{
  // zero out the blocks not in my block list
  int i, j, index;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  for (i=0, index=0; i<d->blocks()->nblock(); i++) {
      for (j=0; j<=i; j++, index++) {
          if (index%nproc == me) continue;
          for (int ii=d->blocks()->start(i); ii<d->blocks()->fence(i); ii++) {
              for (int jj=d->blocks()->start(j);
                   jj < d->blocks()->fence(j) && jj <= ii;
                   jj++) {
                  matrix[(ii*(ii+1)>>1) + jj] = 0.0;
                }
            }
        }
    }
}

void
ReplSymmSCMatrix::after_elemop()
{
  messagegrp()->sum(matrix, d->n()*(d->n()+1)>>1);
}

void
ReplSymmSCMatrix::init_blocklist()
{
  int i, j, index;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  for (i=0, index=0; i<d->blocks()->nblock(); i++) {
      for (j=0; j<=i; j++, index++) {
          if (index%nproc != me) continue;
          blocklist->insert(
              new SCMatrixLTriSubBlock(d->blocks()->start(i),
                                       d->blocks()->fence(i),
                                       d->blocks()->start(j),
                                       d->blocks()->fence(j),
                                       matrix));
        }
    }
}

ReplSymmSCMatrix::~ReplSymmSCMatrix()
{
  if (matrix) delete[] matrix;
  if (rows) delete[] rows;
}

int
ReplSymmSCMatrix::compute_offset(int i,int j)
{
  if (i<0 || j<0 || i>=d->n() || j>=d->n()) {
      fprintf(stderr,"ReplSymmSCMatrix: index out of bounds\n");
      abort();
    }
  if (i<j) {
      int tmp = j; j=i; i=tmp;
    }
  return (i*(i+1))/2 + j;
}

double
ReplSymmSCMatrix::get_element(int i,int j)
{
  return matrix[compute_offset(i,j)];
}

void
ReplSymmSCMatrix::set_element(int i,int j,double a)
{
  matrix[compute_offset(i,j)] = a;
}

void
ReplSymmSCMatrix::accumulate_element(int i,int j,double a)
{
  matrix[compute_offset(i,j)] += a;
}

SCMatrix *
ReplSymmSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    fprintf(stderr,"ReplSymmSCMatrix::get_subblock: trying to get too big a "
            "subblock (%d,%d) from (%d,%d)\n",nsrow,nscol,n(),n());
    abort();
  }
  
  RefSCDimension dnrow = new SCDimension(nsrow);
  RefSCDimension dncol = new SCDimension(nscol);

  SCMatrix * sb = kit()->matrix(dnrow,dncol);
  sb->assign(0.0);

  ReplSCMatrix *lsb = ReplSCMatrix::require_castdown(sb,
                                      "ReplSymmSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      lsb->rows[i][j] = get_element(i+br,j+bc);
      
  return sb;
}

SymmSCMatrix *
ReplSymmSCMatrix::get_subblock(int br, int er)
{
  int nsrow = er-br+1;

  if (nsrow > n()) {
    fprintf(stderr,"ReplSymmSCMatrix::get_subblock: trying to get too big a "
            "subblock (%d,%d) from (%d,%d)\n",nsrow,nsrow,n(),n());
    abort();
  }
  
  RefSCDimension dnrow = new SCDimension(nsrow);

  SymmSCMatrix * sb = kit()->symmmatrix(dnrow);
  sb->assign(0.0);

  ReplSymmSCMatrix *lsb = ReplSymmSCMatrix::require_castdown(sb,
                                      "ReplSymmSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      lsb->rows[i][j] = get_element(i+br,j+br);
      
  return sb;
}

void
ReplSymmSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  ReplSCMatrix *lsb = ReplSCMatrix::require_castdown(sb,
                                      "ReplSCMatrix::assign_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    fprintf(stderr,"ReplSymmSCMatrix::assign_subblock: "
            "trying to assign too big a "
            "subblock (%d,%d) to (%d,%d)\n",nsrow,nscol,n(),n());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      set_element(i+br,j+bc,lsb->rows[i][j]);
}

void
ReplSymmSCMatrix::assign_subblock(SymmSCMatrix*sb, int br, int er)
{
  ReplSymmSCMatrix *lsb = ReplSymmSCMatrix::require_castdown(sb,
                                      "ReplSymmSCMatrix::assign_subblock");

  int nsrow = er-br+1;

  if (nsrow > n()) {
    fprintf(stderr,"ReplSymmSCMatrix::assign_subblock: "
            "trying to assign too big a "
            "subblock (%d,%d) to (%d,%d)\n",nsrow,nsrow,n(),n());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      set_element(i+br,j+br,lsb->rows[i][j]);
}

void
ReplSymmSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  ReplSCMatrix *lsb = ReplSCMatrix::require_castdown(sb,
                                  "ReplSymmSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    fprintf(stderr,"ReplSymmSCMatrix::accumulate_subblock: trying to "
            "accumulate too big a subblock (%d,%d) to (%d,%d)\n",
            nsrow,nscol,n(),n());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      set_element(i+br,j+br,get_element(i+br,j+br)+lsb->rows[i][j]);
}

void
ReplSymmSCMatrix::accumulate_subblock(SymmSCMatrix*sb, int br, int er)
{
  ReplSCMatrix *lsb = ReplSCMatrix::require_castdown(sb,
                                  "ReplSymmSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;

  if (nsrow > n()) {
    fprintf(stderr,"ReplSymmSCMatrix::accumulate_subblock: trying to "
            "accumulate too big a subblock (%d,%d) to (%d,%d)\n",
            nsrow,nsrow,n(),n());
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      set_element(i+br,j+br,get_element(i+br,j+br)+lsb->rows[i][j]);
}

SCVector *
ReplSymmSCMatrix::get_row(int i)
{
  if (i >= n()) {
    fprintf(stderr,"ReplSymmSCMatrix::get_row: trying to get invalid row"
            "%d max %d\n",i,n());
    abort();
  }
  
  SCVector * v = kit()->vector(dim());

  ReplSCVector *lv = ReplSCVector::require_castdown(v,
                                               "ReplSymmSCMatrix::get_row");

  for (int j=0; j < n(); j++)
    lv->set_element(j,get_element(i,j));
      
  return v;
}

void
ReplSymmSCMatrix::assign_row(SCVector *v, int i)
{
  if (i >= n()) {
    fprintf(stderr,"ReplSymmSCMatrix::assign_row: trying to assign invalid "
            "row %d max %d\n",i,n());
    abort();
  }
  
  if (v->n() != n()) {
    fprintf(stderr,"ReplSymmSCMatrix::assign_row: vector is wrong size"
            "is %d, should be %d\n",v->n(),n());
    abort();
  }
  
  ReplSCVector *lv = ReplSCVector::require_castdown(v,
                                          "ReplSymmSCMatrix::assign_row");

  for (int j=0; j < n(); j++)
    set_element(i,j,lv->get_element(j));
}

void
ReplSymmSCMatrix::accumulate_row(SCVector *v, int i)
{
  if (i >= n()) {
    fprintf(stderr,"ReplSymmSCMatrix::accumulate_row: trying to assign "
                   "invalid row %d max %d\n",i,n());
    abort();
  }
  
  if (v->n() != n()) {
    fprintf(stderr,"ReplSymmSCMatrix::accumulate_row: vector is wrong size"
            "is %d, should be %d\n",v->n(),n());
    abort();
  }
  
  ReplSCVector *lv = ReplSCVector::require_castdown(v,
                                        "ReplSymmSCMatrix::accumulate_row");

  for (int j=0; j < n(); j++)
    set_element(i,j,get_element(i,j)+lv->get_element(j));
}

void
ReplSymmSCMatrix::assign(double val)
{
  int n = (d->n()*(d->n()+1))/2;
  for (int i=0; i<n; i++) matrix[i] = val;
}

void
ReplSymmSCMatrix::accumulate(SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  ReplSymmSCMatrix* la
    = ReplSymmSCMatrix::require_castdown(a,"ReplSymmSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      fprintf(stderr,"ReplSymmSCMatrix::accumulate(SCMatrix*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = (this->n() * (this->n() + 1))/2;
  for (int i=0; i<nelem; i++) matrix[i] += la->matrix[i];
}

double
ReplSymmSCMatrix::invert_this()
{
  return cmat_invert(rows,1,n());
}

double
ReplSymmSCMatrix::determ_this()
{
  return cmat_determ(rows,1,n());
}

double
ReplSymmSCMatrix::trace()
{
  double ret=0;
  for (int i=0; i < n(); i++) ret += rows[i][i];
  return ret;
}

double
ReplSymmSCMatrix::solve_this(SCVector*v)
{
  ReplSCVector* lv =
    ReplSCVector::require_castdown(v,"ReplSymmSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!dim()->equiv(lv->dim())) {
      fprintf(stderr,"ReplSymmSCMatrix::solve_this(SCVector*v):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  return cmat_solve_lin(rows,1,lv->vector,n());
}

void
ReplSymmSCMatrix::gen_invert_this()
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
ReplSymmSCMatrix::diagonalize(DiagSCMatrix*a,SCMatrix*b)
{
  const char* name = "ReplSymmSCMatrix::diagonalize";
  // make sure that the arguments is of the correct type
  ReplDiagSCMatrix* la = ReplDiagSCMatrix::require_castdown(a,name);
  ReplSCMatrix* lb = ReplSCMatrix::require_castdown(b,name);

  if (!dim()->equiv(la->dim()) ||
      !dim()->equiv(lb->coldim()) || !dim()->equiv(lb->rowdim())) {
      fprintf(stderr,"ReplSymmSCMatrix::"
              "diagonalize(DiagSCMatrix*a,SCMatrix*b): bad dims");
      abort();
    }

  double *eigvals;
  double **eigvecs;
  if (!la) {
      eigvals = new double[n()];
    }
  else {
      eigvals = la->matrix;
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
ReplSymmSCMatrix::accumulate_symmetric_product(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  ReplSCMatrix* la
    = ReplSCMatrix::require_castdown(a,"ReplSymmSCMatrix::"
                                          "accumulate_symmetric_product");

  if (!dim()->equiv(la->rowdim())) {
      fprintf(stderr,"ReplSymmSCMatrix::"
              "accumulate_symmetric_product(SCMatrix*a): bad dim");
      abort();
    }

  cmat_symmetric_mxm(rows,n(),la->rows,la->ncol(),1);
}

// computes this += a + a.t
void
ReplSymmSCMatrix::accumulate_symmetric_sum(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  ReplSCMatrix* la
    = ReplSCMatrix::require_castdown(a,"ReplSymmSCMatrix::"
                                          "accumulate_symmetric_sum");

  if (!dim()->equiv(la->rowdim()) || !dim()->equiv(la->coldim())) {
      fprintf(stderr,"ReplSymmSCMatrix::"
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
ReplSymmSCMatrix::accumulate_symmetric_outer_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  ReplSCVector* la
    = ReplSCVector::require_castdown(a,"ReplSymmSCMatrix::"
                                      "accumulate_symmetric_outer_product");

  if (!dim()->equiv(la->dim())) {
      fprintf(stderr,"ReplSymmSCMatrix::"
              "accumulate_symmetric_outer_product(SCMatrix*a): bad dim");
      abort();
    }

  int n = dim().n();
  double** tdat = this->rows;
  double* adat = la->vector;
  for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
          tdat[i][j] += adat[i]*adat[j];
        }
    }
}

// this += a * b * transpose(a)
void
ReplSymmSCMatrix::accumulate_transform(SCMatrix*a,SymmSCMatrix*b)
{
  // do the necessary castdowns
  ReplSCMatrix*la
    = ReplSCMatrix::require_castdown(a,"%s::accumulate_transform",
                                      class_name());
  ReplSymmSCMatrix*lb = require_castdown(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
      fprintf(stderr,"ReplSymmSCMatrix::accumulate_transform: bad dim\n");
      abort();
    }

  cmat_transform_symmetric_matrix(rows,n(),lb->rows,lb->n(),la->rows,1);
}

// this += a * b * transpose(a)
void
ReplSymmSCMatrix::accumulate_transform(SCMatrix*a,DiagSCMatrix*b)
{
  // do the necessary castdowns
  ReplSCMatrix*la
    = ReplSCMatrix::require_castdown(a,"%s::accumulate_transform",
                                      class_name());
  ReplDiagSCMatrix*lb
    = ReplDiagSCMatrix::require_castdown(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
      fprintf(stderr,"ReplSymmSCMatrix::accumulate_transform: bad dim\n");
      abort();
    }

  cmat_transform_diagonal_matrix(rows,n(),lb->matrix,lb->n(),la->rows,1);
}

double
ReplSymmSCMatrix::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  ReplSCVector* la
    = ReplSCVector::require_castdown(a,"ReplSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      fprintf(stderr,"ReplSCVector::"
              "scale_product(SCVector*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = n();
  double* adat = la->vector;
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
ReplSymmSCMatrix::element_op(const RefSCElementOp& op)
{
  if (op->has_side_effects()) before_elemop();
  SCMatrixBlockListIter i;
  for (i = blocklist->begin(); i != blocklist->end(); i++) {
      op->process_base(i.block());
    }
  if (op->has_side_effects()) after_elemop();
}

void
ReplSymmSCMatrix::element_op(const RefSCElementOp2& op,
                              SymmSCMatrix* m)
{
  ReplSymmSCMatrix *lm
      = ReplSymmSCMatrix::require_castdown(m,"ReplSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim())) {
      fprintf(stderr,"ReplSymmSCMatrix: bad element_op\n");
      abort();
    }
  if (op->has_side_effects()) before_elemop();
  if (op->has_side_effects_in_arg()) lm->before_elemop();
  SCMatrixBlockListIter i, j;
  for (i = blocklist->begin(), j = lm->blocklist->begin();
       i != blocklist->end();
       i++, j++) {
      op->process_base(i.block(), j.block());
    }
  if (op->has_side_effects()) after_elemop();
  if (op->has_side_effects_in_arg()) lm->after_elemop();
}

void
ReplSymmSCMatrix::element_op(const RefSCElementOp3& op,
                              SymmSCMatrix* m,SymmSCMatrix* n)
{
  ReplSymmSCMatrix *lm
      = ReplSymmSCMatrix::require_castdown(m,"ReplSymSCMatrix::element_op");
  ReplSymmSCMatrix *ln
      = ReplSymmSCMatrix::require_castdown(n,"ReplSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      fprintf(stderr,"ReplSymmSCMatrix: bad element_op\n");
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
      op->process_base(i.block(), j.block(), k.block());
    }
  if (op->has_side_effects()) after_elemop();
  if (op->has_side_effects_in_arg1()) lm->after_elemop();
  if (op->has_side_effects_in_arg2()) ln->after_elemop();
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
ReplSymmSCMatrix::print(const char *title, ostream& os, int prec)
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

RefSCMatrixSubblockIter
ReplSymmSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new ReplSCMatrixListSubblockIter(access, blocklist,
                                          messagegrp(),
                                          matrix, (d->n()*(d->n()+1))/2);
}

RefSCMatrixSubblockIter
ReplSymmSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      cerr << "ReplSymmSCMatrix::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  RefSCMatrixBlockList allblocklist = new SCMatrixBlockList();
  allblocklist->insert(new SCMatrixLTriSubBlock(0, d->n(),
                                                0, d->n(), matrix));
  return new ReplSCMatrixListSubblockIter(access, allblocklist,
                                          messagegrp(),
                                          matrix, (d->n()*(d->n()+1))/2);
}

RefReplSCMatrixKit
ReplSymmSCMatrix::skit()
{
  return ReplSCMatrixKit::castdown(kit().pointer());
}

RefMessageGrp
ReplSymmSCMatrix::messagegrp()
{
  return skit()->messagegrp();
}
