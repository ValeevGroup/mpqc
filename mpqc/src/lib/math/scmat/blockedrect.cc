
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// BlockedSCMatrix member functions

#define CLASSNAME BlockedSCMatrix
#define PARENTS public SCMatrix
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BlockedSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
BlockedSCMatrix::resize(BlockedSCDimension *a, BlockedSCDimension *b)
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
  
  d1 = a;
  d2 = b;
  
  if (!a || !b || !a->nblocks() || !b->nblocks())
    return;

  if (a->nblocks() != b->nblocks()) {
    fprintf(stderr,"BlockedSCMatrix::resize: wrong number of blocks\n");
    abort();
  }

  mats_ = new RefSCMatrix[d1->nblocks()];
  for (int i=0; i < d1->nblocks(); i++)
    mats_[i] = d1->dim(i)->create_matrix(d2->dim(i).pointer());
}


BlockedSCMatrix::BlockedSCMatrix() :
  mats_(0)
{
}

BlockedSCMatrix::BlockedSCMatrix(BlockedSCDimension*a,BlockedSCDimension*b):
  mats_(0)
{
  resize(a,b);
}

BlockedSCMatrix::BlockedSCMatrix(StateIn&s):
  SCMatrix(s)
{
  d1.restore_state(s);
  d2.restore_state(s);
  mats_ = new RefSCMatrix[d1->nblocks()];
  for (int i=0; i < d1->nblocks(); i++)
    mats_[i].restore_state(s);
}

BlockedSCMatrix::BlockedSCMatrix(const RefKeyVal&keyval)
{
  abort();
}

void
BlockedSCMatrix::save_data_state(StateOut&s)
{
  SCMatrix::save_data_state(s);
  d1.save_state(s);
  d2.save_state(s);
  for (int i=0; i < d1->nblocks(); i++)
    mats_[i].save_state(s);
}

BlockedSCMatrix::~BlockedSCMatrix()
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
}

RefSCDimension
BlockedSCMatrix::rowdim()
{
  return d1;
}

RefSCDimension
BlockedSCMatrix::coldim()
{
  return d2;
}

RefSCDimension
BlockedSCMatrix::rowdim(int i)
{
  return d1->dim(i);
}

RefSCDimension
BlockedSCMatrix::coldim(int i)
{
  return d2->dim(i);
}

void
BlockedSCMatrix::assign(double v)
{
  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->assign(v);
}

double
BlockedSCMatrix::get_element(int i,int j)
{
  int block_i = d1->block(i);
  int block_j = d2->block(j);

  if (block_i != block_j)
    return 0;
  
  int elem_i = i - d1->first(block_i);
  int elem_j = j - d2->first(block_j);
  
  return mats_[block_i]->get_element(elem_i,elem_j);
}

void
BlockedSCMatrix::set_element(int i,int j,double a)
{
  int block_i = d1->block(i);
  int block_j = d2->block(j);

  if (block_i != block_j)
    return;
  
  int elem_i = i - d1->first(block_i);
  int elem_j = j - d2->first(block_j);
  
  mats_[block_i]->set_element(elem_i,elem_j,a);
}

void
BlockedSCMatrix::accumulate_element(int i,int j,double a)
{
  int block_i = d1->block(i);
  int block_j = d2->block(j);

  if (block_i != block_j)
    return;
  
  int elem_i = i - d1->first(block_i);
  int elem_j = j - d2->first(block_j);
  
  mats_[block_i]->accumulate_element(elem_i,elem_j,a);
}

SCMatrix *
BlockedSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  fprintf(stderr,"BlockedSCMatrix::get_subblock: cannot get subblock\n");
  abort();
  return 0;
}

void
BlockedSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                               int source_br, int source_bc)
{
  fprintf(stderr,"BlockedSCMatrix::assign_subblock: cannot assign subblock\n");
  abort();
}

void
BlockedSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                                   int source_br, int source_bc)
{
  fprintf(stderr,"BlockedSCMatrix::accumulate_subblock:"
          " cannot accumulate subblock\n");
  abort();
}

SCVector *
BlockedSCMatrix::get_row(int i)
{
  fprintf(stderr,"BlockedSCMatrix::get_row: cannot get row\n");
  abort();

  return 0;
}

void
BlockedSCMatrix::assign_row(SCVector *v, int i)
{
  fprintf(stderr,"BlockedSCMatrix::assign_row: cannot assign row\n");
  abort();
}

void
BlockedSCMatrix::accumulate_row(SCVector *v, int i)
{
  fprintf(stderr,"BlockedSCMatrix::accumulate_row: cannot accumulate row\n");
  abort();
}

SCVector *
BlockedSCMatrix::get_column(int i)
{
  fprintf(stderr,"BlockedSCMatrix::get_column: cannot get column\n");
  abort();

  return 0;
}

void
BlockedSCMatrix::assign_column(SCVector *v, int i)
{
  fprintf(stderr,"BlockedSCMatrix::assign_column: cannot assign column\n");
  abort();
}

void
BlockedSCMatrix::accumulate_column(SCVector *v, int i)
{
  fprintf(stderr,"BlockedSCMatrix::accumulate_column: cannot accumulate column\n");
  abort();
}

// does the outer product a x b.  this must have rowdim() == a->dim() and
// coldim() == b->dim()
void
BlockedSCMatrix::accumulate_outer_product(SCVector*a,SCVector*b)
{
  const char* name = "BlockedSCMatrix::accumulate_outer_product";
  // make sure that the arguments are of the correct type
  BlockedSCVector* la = BlockedSCVector::require_castdown(a,name);
  BlockedSCVector* lb = BlockedSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!(this->rowdim() == a->dim())
      || !(this->coldim() == b->dim())) {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate_outer_product(SCVector*a,SCVector*b):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->accumulate_outer_product(la->vecs_[i], lb->vecs_[i]);
}

void
BlockedSCMatrix::accumulate_product(SCMatrix*a,SCMatrix*b)
{
  const char* name = "BlockedSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = BlockedSCMatrix::require_castdown(a,name);
  BlockedSCMatrix* lb = BlockedSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!(this->rowdim() == a->rowdim())
      || !(this->coldim() == b->coldim())
      || !(a->coldim() == b->rowdim())) {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate_product(SCMatrix*a,SCMatrix*b):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->accumulate_product(la->mats_[i], lb->mats_[i]);
}

void
BlockedSCMatrix::accumulate_product(SCMatrix*a,SymmSCMatrix*b)
{
  const char* name = "BlockedSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = BlockedSCMatrix::require_castdown(a,name);
  BlockedSymmSCMatrix* lb = BlockedSymmSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!(this->rowdim() == a->rowdim())
      || !(this->coldim() == b->dim())
      || !(a->coldim() == b->dim())) {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate_product(SCMatrix*a,SymmSCMatrix*b):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->accumulate_product(la->mats_[i], lb->mats_[i]);
}

void
BlockedSCMatrix::accumulate_product(SCMatrix*a,DiagSCMatrix*b)
{
  const char* name = "BlockedSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = BlockedSCMatrix::require_castdown(a,name);
  BlockedDiagSCMatrix* lb = BlockedDiagSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!(this->rowdim() == a->rowdim())
      || !(this->coldim() == b->dim())
      || !(a->coldim() == b->dim())) {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate_product(SCMatrix*a,DiagSCMatrix*b):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->accumulate_product(la->mats_[i], lb->mats_[i]);
}

void
BlockedSCMatrix::accumulate(SCMatrix*a)
{
  // make sure that the arguments is of the correct type
  BlockedSCMatrix* la
    = BlockedSCMatrix::require_castdown(a,"BlockedSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!(this->rowdim() == a->rowdim())
      || !(this->coldim() == a->coldim())) {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate(SCMatrix*a):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->accumulate(la->mats_[i]);
}

void
BlockedSCMatrix::transpose_this()
{
  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->transpose_this();
  
  RefBlockedSCDimension tmp = d1;
  d1 = d2;
  d2 = tmp;
}

double
BlockedSCMatrix::invert_this()
{
  double res=1;
  
  for (int i=0; i < d1->nblocks(); i++)
    res *= mats_[i]->invert_this();

  return res;
}

void
BlockedSCMatrix::gen_invert_this()
{
  fprintf(stderr,"BlockedSCMatrix::gen_invert_this: SVD not implemented yet");
  abort();
}

double
BlockedSCMatrix::determ_this()
{
  double res=1;
  
  for (int i=0; i < d1->nblocks(); i++)
    res *= mats_[i]->determ_this();

  return res;
}

double
BlockedSCMatrix::trace()
{
  double ret=0;
  for (int i=0; i < d1->nblocks(); i++)
    ret += mats_[i]->trace();
  
  return ret;
}

void
BlockedSCMatrix::svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V)
{
  BlockedSCMatrix* lU =
    BlockedSCMatrix::require_castdown(U,"BlockedSCMatrix::svd_this");
  BlockedSCMatrix* lV =
    BlockedSCMatrix::require_castdown(V,"BlockedSCMatrix::svd_this");
  BlockedDiagSCMatrix* lsigma =
    BlockedDiagSCMatrix::require_castdown(sigma,"BlockedSCMatrix::svd_this");

  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->svd_this(lU->mats_[i].pointer(),
                       lsigma->mats_[i].pointer(),
                       lV->mats_[i].pointer());
}

double
BlockedSCMatrix::solve_this(SCVector*v)
{
  BlockedSCVector* lv =
    BlockedSCVector::require_castdown(v,"BlockedSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!(this->rowdim() == v->dim())) {
    fprintf(stderr,"BlockedSCMatrix::solve_this(SCVector*v):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->solve_this(lv->vecs_[i].pointer());
}

void
BlockedSCMatrix::schmidt_orthog(SymmSCMatrix *S, int nc)
{
  BlockedSymmSCMatrix* lS =
    BlockedSymmSCMatrix::require_castdown(S,"BlockedSCMatrix::schmidt_orthog");
  
  // make sure that the dimensions match
  if (!(this->rowdim() == S->dim())) {
    fprintf(stderr,"BlockedSCMatrix::schmidt_orthog():\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->schmidt_orthog(lS->mats_[i].pointer(), lS->dim(i).n());
}

void
BlockedSCMatrix::element_op(const RefSCElementOp& op)
{
  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->element_op(op);
}

void
BlockedSCMatrix::element_op(const RefSCElementOp2& op,
                          SCMatrix* m)
{
  BlockedSCMatrix *lm
      = BlockedSCMatrix::require_castdown(m,"BlockedSCMatrix::element_op");
  if (!lm || d1 != lm->d1 || d2 != lm->d2) {
    fprintf(stderr,"BlockedSCMatrix: bad element_op\n");
    abort();
  }

  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->element_op(op,lm->mats_[i].pointer());
}

void
BlockedSCMatrix::element_op(const RefSCElementOp3& op,
                          SCMatrix* m,SCMatrix* n)
{
  BlockedSCMatrix *lm
      = BlockedSCMatrix::require_castdown(m,"BlockedSCMatrix::element_op");
  BlockedSCMatrix *ln
      = BlockedSCMatrix::require_castdown(n,"BlockedSCMatrix::element_op");
  if (!lm || !ln
      || d1 != lm->d1 || d2 != lm->d2 || d1 != ln->d1 || d2 != ln->d2) {
    fprintf(stderr,"BlockedSCMatrix: bad element_op\n");
    abort();
  }

  for (int i=0; i < d1->nblocks(); i++)
    mats_[i]->element_op(op,lm->mats_[i].pointer(),
                            ln->mats_[i].pointer());
}

// from Ed Seidl at the NIH
void
BlockedSCMatrix::print(const char *title, ostream& os, int prec)
{
  char *newtitle = new char[strlen(title) + 80];

  for (int i=0; i < d1->nblocks(); i++) {
    sprintf(newtitle,"%s:  block %d",title,i+1);
    mats_[i]->print(newtitle, os, prec);
  }

  delete[] newtitle;
}
