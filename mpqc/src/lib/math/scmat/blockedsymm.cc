
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// BlockedSymmSCMatrix member functions

#define CLASSNAME BlockedSymmSCMatrix
#define PARENTS public SymmSCMatrix
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BlockedSymmSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SymmSCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
BlockedSymmSCMatrix::resize(BlockedSCDimension *a)
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
  
  d = a;
  
  mats_ = new RefSymmSCMatrix[d->nblocks()];
  for (int i=0; i < d->nblocks(); i++)
    mats_[i] = d->dim(i)->create_symmmatrix();
}

BlockedSymmSCMatrix::BlockedSymmSCMatrix() :
  mats_(0)
{
}

BlockedSymmSCMatrix::BlockedSymmSCMatrix(BlockedSCDimension*a) :
  mats_(0)
{
  resize(a);
}

BlockedSymmSCMatrix::BlockedSymmSCMatrix(StateIn&s):
  SymmSCMatrix(s)
{
  d.restore_state(s);
  mats_ = new RefSymmSCMatrix[d->nblocks()];
  for (int i=0; i < d->nblocks(); i++)
    mats_[i].restore_state(s);
}

BlockedSymmSCMatrix::BlockedSymmSCMatrix(const RefKeyVal&keyval)
{
  abort();
}

void
BlockedSymmSCMatrix::save_data_state(StateOut&s)
{
  SymmSCMatrix::save_data_state(s);
  d.save_state(s);
  for (int i=0; i < d->nblocks(); i++)
    mats_[i].save_state(s);
}

BlockedSymmSCMatrix::~BlockedSymmSCMatrix()
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
}

RefSCDimension
BlockedSymmSCMatrix::dim()
{
  return d;
}

RefSCDimension
BlockedSymmSCMatrix::dim(int i)
{
  return d->dim(i);
}

double
BlockedSymmSCMatrix::get_element(int i,int j)
{
  int block_i = d->block(i);
  int block_j = d->block(j);

  if (block_i != block_j)
    return 0;
  
  int elem_i = i - d->first(block_i);
  int elem_j = j - d->first(block_j);
  
  return mats_[block_i]->get_element(elem_i,elem_j);
}

void
BlockedSymmSCMatrix::set_element(int i,int j,double a)
{
  int block_i = d->block(i);
  int block_j = d->block(j);

  if (block_i != block_j)
    return;
  
  int elem_i = i - d->first(block_i);
  int elem_j = j - d->first(block_j);
  
  mats_[block_i]->set_element(elem_i,elem_j,a);
}

void
BlockedSymmSCMatrix::accumulate_element(int i,int j,double a)
{
  int block_i = d->block(i);
  int block_j = d->block(j);

  if (block_i != block_j)
    return;
  
  int elem_i = i - d->first(block_i);
  int elem_j = j - d->first(block_j);
  
  mats_[block_i]->accumulate_element(elem_i,elem_j,a);
}

SCMatrix *
BlockedSymmSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  fprintf(stderr,"BlockedSymmSCMatrix::get_subblock: cannot get subblock\n");
  abort();
  return 0;
}

SymmSCMatrix *
BlockedSymmSCMatrix::get_subblock(int br, int er)
{
  fprintf(stderr,"BlockedSymmSCMatrix::get_subblock: cannot get subblock\n");
  abort();
  return 0;
}

void
BlockedSymmSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  fprintf(stderr,"BlockedSymmSCMatrix::assign_subblock:"
          " cannot assign subblock\n");
  abort();
}

void
BlockedSymmSCMatrix::assign_subblock(SymmSCMatrix*sb, int br, int er)
{
  fprintf(stderr,"BlockedSymmSCMatrix::assign_subblock:"
          " cannot assign subblock\n");
  abort();
}

void
BlockedSymmSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  fprintf(stderr,"BlockedSymmSCMatrix::accumulate_subblock:"
          " cannot accumulate subblock\n");
  abort();
}

void
BlockedSymmSCMatrix::accumulate_subblock(SymmSCMatrix*sb, int br, int er)
{
  fprintf(stderr,"BlockedSymmSCMatrix::accumulate_subblock:"
          " cannot accumulate subblock\n");
  abort();
}

SCVector *
BlockedSymmSCMatrix::get_row(int i)
{
  fprintf(stderr,"BlockedSymmSCMatrix::get_row: cannot get row\n");
  abort();

  return 0;
}

void
BlockedSymmSCMatrix::assign_row(SCVector *v, int i)
{
  fprintf(stderr,"BlockedSymmSCMatrix::assign_row: cannot assign row\n");
  abort();
}

void
BlockedSymmSCMatrix::accumulate_row(SCVector *v, int i)
{
  fprintf(stderr,"BlockedSymmSCMatrix::accumulate_row:"
          " cannot accumulate row\n");
  abort();
}

double
BlockedSymmSCMatrix::invert_this()
{
  double res=1;
  
  for (int i=0; i < d->nblocks(); i++)
    res *= mats_[i]->invert_this();

  return res;
}

double
BlockedSymmSCMatrix::determ_this()
{
  double res=1;
  
  for (int i=0; i < d->nblocks(); i++)
    res *= mats_[i]->determ_this();

  return res;
}

double
BlockedSymmSCMatrix::trace()
{
  double res=0;
  
  for (int i=0; i < d->nblocks(); i++)
    res += mats_[i]->determ_this();

  return res;
}

double
BlockedSymmSCMatrix::solve_this(SCVector*v)
{
  double res=1;
  
  BlockedSCVector* lv =
    BlockedSCVector::require_castdown(v,"BlockedSymmSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!dim()->equiv(lv->dim())) {
    fprintf(stderr,"BlockedSymmSCMatrix::solve_this(SCVector*v):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    res *= mats_[i]->solve_this(lv->vecs_[i].pointer());

  return res;
}

void
BlockedSymmSCMatrix::gen_invert_this()
{
  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->gen_invert_this();
}

double
BlockedSymmSCMatrix::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  BlockedSCVector* la
    = BlockedSCVector::require_castdown(a,"BlockedSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    fprintf(stderr,"BlockedSCVector::scale_product(SCVector*a):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  double result = 0.0;
  for (int i=0; i < d->nblocks(); i++)
    result += mats_[i]->scalar_product(la->vecs_[i].pointer());

  return result;
}

void
BlockedSymmSCMatrix::diagonalize(DiagSCMatrix*a,SCMatrix*b)
{
  const char* name = "BlockedSymmSCMatrix::diagonalize";
  // make sure that the arguments is of the correct type
  BlockedDiagSCMatrix* la = BlockedDiagSCMatrix::require_castdown(a,name);
  BlockedSCMatrix* lb = BlockedSCMatrix::require_castdown(b,name);

  if (!dim()->equiv(la->dim()) ||
      !dim()->equiv(lb->coldim()) || !dim()->equiv(lb->rowdim())) {
    fprintf(stderr,"BlockedSymmSCMatrix::"
            "diagonalize(DiagSCMatrix*a,SCMatrix*b): bad dims");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->diagonalize(la->mats_[i].pointer(),lb->mats_[i].pointer());
}

void
BlockedSymmSCMatrix::accumulate(SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  BlockedSymmSCMatrix* la = BlockedSymmSCMatrix::require_castdown(a,
                                   "BlockedSymmSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    fprintf(stderr,"BlockedSymmSCMatrix::accumulate(SCMatrix*a):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->accumulate(la->mats_[i].pointer());
}

// computes this += a * a.t
void
BlockedSymmSCMatrix::accumulate_symmetric_product(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  BlockedSCMatrix* la
    = BlockedSCMatrix::require_castdown(a,"BlockedSymmSCMatrix::"
                                          "accumulate_symmetric_product");

  if (!dim()->equiv(la->rowdim())) {
    fprintf(stderr,"BlockedSymmSCMatrix::"
            "accumulate_symmetric_product(SCMatrix*a): bad dim");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->accumulate_symmetric_product(la->mats_[i].pointer());
}

// computes this += a + a.t
void
BlockedSymmSCMatrix::accumulate_symmetric_sum(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  BlockedSCMatrix* la
    = BlockedSCMatrix::require_castdown(a,"BlockedSymmSCMatrix::"
                                          "accumulate_symmetric_sum");

  if (!dim()->equiv(la->rowdim()) || !dim()->equiv(la->coldim())) {
    fprintf(stderr,"BlockedSymmSCMatrix::"
            "accumulate_symmetric_sum(SCMatrix*a): bad dim");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->accumulate_symmetric_sum(la->mats_[i].pointer());
}

void
BlockedSymmSCMatrix::accumulate_symmetric_outer_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  BlockedSCVector* la
    = BlockedSCVector::require_castdown(a,"BlockedSymmSCMatrix::"
                                      "accumulate_symmetric_outer_product");

  if (!dim()->equiv(la->dim())) {
    fprintf(stderr,"BlockedSymmSCMatrix::"
            "accumulate_symmetric_outer_product(SCMatrix*a): bad dim");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->accumulate_symmetric_outer_product(la->vecs_[i].pointer());
}

// this += a * b * transpose(a)
void
BlockedSymmSCMatrix::accumulate_transform(SCMatrix*a,SymmSCMatrix*b)
{
  // do the necessary castdowns
  BlockedSCMatrix*la
    = BlockedSCMatrix::require_castdown(a,"%s::accumulate_transform",
                                      class_name());
  BlockedSymmSCMatrix*lb = require_castdown(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
    fprintf(stderr,"BlockedSymmSCMatrix::accumulate_transform: bad dim\n");
    abort();
  }

  if (lb->d->nblocks() == 1 && d->nblocks() > 1) {
    for (int i=0; i < d->nblocks(); i++)
      mats_[i]->accumulate_transform(la->mats_[i].pointer(),
                                     lb->mats_[0].pointer());

  } else if (lb->d->nblocks() > 1 && d->nblocks() == 1) {
    for (int i=0; i < lb->d->nblocks(); i++)
      mats_[0]->accumulate_transform(la->mats_[i].pointer(),
                                     lb->mats_[i].pointer());

  } else if (lb->d->nblocks() == d->nblocks()) {
    for (int i=0; i < d->nblocks(); i++)
      mats_[i]->accumulate_transform(la->mats_[i].pointer(),
                                     lb->mats_[i].pointer());
  }

}

// this += a * b * transpose(a)
void
BlockedSymmSCMatrix::accumulate_transform(SCMatrix*a,DiagSCMatrix*b)
{
  // do the necessary castdowns
  BlockedSCMatrix*la
    = BlockedSCMatrix::require_castdown(a,"%s::accumulate_transform",
                                      class_name());
  BlockedDiagSCMatrix*lb
    = BlockedDiagSCMatrix::require_castdown(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
    fprintf(stderr,"BlockedSymmSCMatrix::accumulate_transform: bad dim\n");
    abort();
  }

  if (lb->d->nblocks() == 1 && d->nblocks() > 1) {
    for (int i=0; i < d->nblocks(); i++)
      mats_[i]->accumulate_transform(la->mats_[i].pointer(),
                                     lb->mats_[0].pointer());

  } else if (lb->d->nblocks() > 1 && d->nblocks() == 1) {
    for (int i=0; i < lb->d->nblocks(); i++)
      mats_[0]->accumulate_transform(la->mats_[i].pointer(),
                                     lb->mats_[i].pointer());

  } else if (lb->d->nblocks() == d->nblocks()) {
    for (int i=0; i < d->nblocks(); i++)
      mats_[i]->accumulate_transform(la->mats_[i].pointer(),
                                     lb->mats_[i].pointer());
  }
}

void
BlockedSymmSCMatrix::element_op(const RefSCElementOp& op)
{
  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op);

  for (int i=0; i < d->nblocks(); i++) {
    if (bop)
      bop->working_on(i);
    mats_[i]->element_op(op);
  }
}

void
BlockedSymmSCMatrix::element_op(const RefSCElementOp2& op,
                                SymmSCMatrix* m)
{
  BlockedSymmSCMatrix *lm = BlockedSymmSCMatrix::require_castdown(m,
                                          "BlockedSymSCMatrix::element_op");
  if (!dim()->equiv(lm->dim())) {
    fprintf(stderr,"BlockedSymmSCMatrix: bad element_op\n");
    abort();
  }

  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op);

  for (int i=0; i < d->nblocks(); i++) {
    if (bop)
      bop->working_on(i);
    mats_[i]->element_op(op,lm->mats_[i].pointer());
  }
}

void
BlockedSymmSCMatrix::element_op(const RefSCElementOp3& op,
                              SymmSCMatrix* m,SymmSCMatrix* n)
{
  BlockedSymmSCMatrix *lm = BlockedSymmSCMatrix::require_castdown(m,
                                        "BlockedSymSCMatrix::element_op");
  BlockedSymmSCMatrix *ln = BlockedSymmSCMatrix::require_castdown(n,
                                        "BlockedSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
    fprintf(stderr,"BlockedSymmSCMatrix: bad element_op\n");
    abort();
  }

  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op);

  for (int i=0; i < d->nblocks(); i++) {
    if (bop)
      bop->working_on(i);
    mats_[i]->element_op(op,lm->mats_[i].pointer(),
                            ln->mats_[i].pointer());
  }
}

void
BlockedSymmSCMatrix::print(const char *title, ostream& os, int prec)
{
  int len = (title) ? strlen(title) : 0;
  char *newtitle = new char[len + 80];

  for (int i=0; i < d->nblocks(); i++) {
    sprintf(newtitle,"%s:  block %d",title,i+1);
    mats_[i]->print(newtitle, os, prec);
  }

  delete[] newtitle;
}
