
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// BlockedSCVector member functions

#define CLASSNAME BlockedSCVector
#define PARENTS public SCVector
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BlockedSCVector::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCVector::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
BlockedSCVector::resize(BlockedSCDimension *bsd)
{
  if (vecs_) {
    delete[] vecs_;
    vecs_=0;
  }

  d = bsd;
  
  if (!bsd || !bsd->nblocks())
    return;
  
  vecs_ = new RefSCVector[d->nblocks()];

  for (int i=0; i < d->nblocks(); i++)
    if (d->n(i))
      vecs_[i] = d->dim(i)->create_vector();
}

BlockedSCVector::BlockedSCVector() :
  vecs_(0)
{
}

BlockedSCVector::BlockedSCVector(BlockedSCDimension*a) :
  vecs_(0)
{
  resize(a);
}

BlockedSCVector::BlockedSCVector(const RefKeyVal&keyval)
{
  abort();
}

BlockedSCVector::BlockedSCVector(StateIn&s):
  SCVector(s)
{
  d.restore_state(s);

  int nb = d->nblocks();
  vecs_ = new RefSCVector[nb];
  
  for (int i=0; i < nb; i++)
    vecs_[i].restore_state(s);
}

void
BlockedSCVector::save_data_state(StateOut&s)
{
  SCVector::save_data_state(s);
  d.save_state(s);
  for (int i=0; i < d->nblocks(); i++)
    vecs_[i].save_state(s);
}

BlockedSCVector::~BlockedSCVector()
{
  if (vecs_) {
    delete[] vecs_;
    vecs_=0;
  }
}

RefSCDimension
BlockedSCVector::dim()
{
  return d;
}

RefSCDimension
BlockedSCVector::dim(int i)
{
  return d->dim(i);
}

RefSCVector
BlockedSCVector::vector(int i)
{
  return vecs_[i];
}

void
BlockedSCVector::assign(double a)
{
  for (int i=0; i < d->nblocks(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->assign(a);
}

void
BlockedSCVector::assign(SCVector*a)
{
  // make sure that the argument is of the correct type
  BlockedSCVector* la
    = BlockedSCVector::require_castdown(a,"BlockedSCVector::assign");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    fprintf(stderr,"BlockedSCVector::assign(SCVector*a):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->assign(la->vecs_[i]);
}

void
BlockedSCVector::assign(const double*a)
{
  for (int i=0; i < d->nblocks(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->assign(a+d->first(i));
}

double
BlockedSCVector::get_element(int i)
{
  int size = d->n();
  if (i < 0 || i >= size) {
    fprintf(stderr,"BlockedSCVector::get_element: bad index\n");
    abort();
  }

  int bi = d->block(i);
  return vecs_[bi].get_element(i-d->first(bi));
}

void
BlockedSCVector::set_element(int i,double a)
{
  int size = d->n();
  if (i < 0 || i >= size) {
    fprintf(stderr,"BlockedSCVector::set_element: bad index\n");
    abort();
  }

  int bi = d->block(i);
  vecs_[bi].set_element(i-d->first(bi),a);
}

void
BlockedSCVector::accumulate_element(int i,double a)
{
  int size = d->n();
  if (i < 0 || i >= size) {
    fprintf(stderr,"BlockedSCVector::accumulate_element: bad index\n");
    abort();
  }

  int bi = d->block(i);
  vecs_[bi].accumulate_element(i-d->first(bi),a);
}

void
BlockedSCVector::accumulate_product(SCMatrix*a,SCVector*b)
{
  const char* name = "BlockedSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = BlockedSCMatrix::require_castdown(a,name);
  BlockedSCVector* lb = BlockedSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->rowdim()) || !la->coldim()->equiv(lb->dim())) {
    fprintf(stderr,"BlockedSCVector::"
            "accumulate_product(SCMatrix*a,SCVector*b):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->accumulate_product(la->mats_[i], lb->vecs_[i]);
}

void
BlockedSCVector::accumulate_product(SymmSCMatrix*a,SCVector*b)
{
  const char* name = "BlockedSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSymmSCMatrix* la = BlockedSymmSCMatrix::require_castdown(a,name);
  BlockedSCVector* lb = BlockedSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim()) || !la->dim()->equiv(lb->dim())) {
    fprintf(stderr,"BlockedSCVector::"
            "accumulate_product(SymmSCMatrix*a,SCVector*b):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->accumulate_product(la->mats_[i], lb->vecs_[i]);
}

void
BlockedSCVector::accumulate(SCVector*a)
{
  // make sure that the argument is of the correct type
  BlockedSCVector* la
    = BlockedSCVector::require_castdown(a,"BlockedSCVector::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    fprintf(stderr,"BlockedSCVector::accumulate(SCVector*a):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->accumulate(la->vecs_[i]);
}

double
BlockedSCVector::scalar_product(SCVector*a)
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

  double result=0;

  for (int i=0; i < d->nblocks(); i++)
    if (vecs_[i].nonnull())
      result += vecs_[i]->scalar_product(la->vecs_[i]);
  
  return result;
}

void
BlockedSCVector::element_op(const RefSCElementOp& op)
{
  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op);

  for (int i=0; i < d->nblocks(); i++) {
    if (bop)
      bop->working_on(i);
    if (vecs_[i].nonnull())
      vecs_[i]->element_op(op);
  }
}

void
BlockedSCVector::element_op(const RefSCElementOp2& op,
                          SCVector* m)
{
  BlockedSCVector *lm
      = BlockedSCVector::require_castdown(m, "BlockedSCVector::element_op");

  if (!dim()->equiv(lm->dim())) {
    fprintf(stderr,"BlockedSCVector: bad element_op\n");
    abort();
  }

  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op);

  for (int i=0; i < d->nblocks(); i++) {
    if (bop)
      bop->working_on(i);
    if (vecs_[i].nonnull())
      vecs_[i]->element_op(op, lm->vecs_[i]);
  }
}

void
BlockedSCVector::element_op(const RefSCElementOp3& op,
                          SCVector* m,SCVector* n)
{
  BlockedSCVector *lm
      = BlockedSCVector::require_castdown(m, "BlockedSCVector::element_op");
  BlockedSCVector *ln
      = BlockedSCVector::require_castdown(n, "BlockedSCVector::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
    fprintf(stderr,"BlockedSCVector: bad element_op\n");
    abort();
  }

  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op);

  for (int i=0; i < d->nblocks(); i++) {
    if (bop)
      bop->working_on(i);
    if (vecs_[i].nonnull())
      vecs_[i]->element_op(op, lm->vecs_[i], ln->vecs_[i]);
  }
}

void
BlockedSCVector::print(const char *title, ostream& os, int prec)
{
  int len = (title) ? strlen(title) : 0;
  char *newtitle = new char[len + 80];

  for (int i=0; i < d->nblocks(); i++) {
    if (vecs_[i].null())
      continue;
    
    sprintf(newtitle,"%s:  block %d",title,i+1);
    vecs_[i]->print(newtitle, os, prec);
  }

  delete[] newtitle;
}
