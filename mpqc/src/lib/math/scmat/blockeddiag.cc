
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// BlockedDiagSCMatrix member functions

#define CLASSNAME BlockedDiagSCMatrix
#define PARENTS public DiagSCMatrix
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BlockedDiagSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DiagSCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
BlockedDiagSCMatrix::resize(BlockedSCDimension *a)
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }

  d = a;

  mats_ = new RefDiagSCMatrix[d->nblocks()];
  for (int i=0; i < d->nblocks(); i++)
    mats_[i] = d->dim(i)->create_diagmatrix();
}

BlockedDiagSCMatrix::BlockedDiagSCMatrix() :
  mats_(0)
{
}

BlockedDiagSCMatrix::BlockedDiagSCMatrix(BlockedSCDimension*a) :
  mats_(0)
{
  resize(a);
}

BlockedDiagSCMatrix::BlockedDiagSCMatrix(StateIn&s):
  DiagSCMatrix(s)
{
  d.restore_state(s);
  mats_ = new RefDiagSCMatrix[d->nblocks()];
  for (int i=0; i < d->nblocks(); i++)
    mats_[i].restore_state(s);
}

BlockedDiagSCMatrix::BlockedDiagSCMatrix(const RefKeyVal&keyval)
{
  abort();
}

void
BlockedDiagSCMatrix::save_data_state(StateOut&s)
{
  DiagSCMatrix::save_data_state(s);
  d.save_state(s);
  for (int i=0; i < d->nblocks(); i++)
    mats_[i].save_state(s);
}

BlockedDiagSCMatrix::~BlockedDiagSCMatrix()
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
}

RefSCDimension
BlockedDiagSCMatrix::dim()
{
  return d;
}

RefSCDimension
BlockedDiagSCMatrix::dim(int i)
{
  return d->dim(i);
}

double
BlockedDiagSCMatrix::get_element(int i)
{
  int block_i = d->block(i);
  int elem_i = i - d->first(block_i);
  
  return mats_[block_i]->get_element(elem_i);
}

void
BlockedDiagSCMatrix::set_element(int i,double a)
{
  int block_i = d->block(i);
  int elem_i = i - d->first(block_i);
  mats_[block_i]->set_element(elem_i,a);
}

void
BlockedDiagSCMatrix::accumulate_element(int i,double a)
{
  int block_i = d->block(i);
  int elem_i = i - d->first(block_i);
  mats_[block_i]->accumulate_element(elem_i,a);
}

void
BlockedDiagSCMatrix::accumulate(DiagSCMatrix*a)
{
  // make sure that the argument is of the correct type
  BlockedDiagSCMatrix* la = BlockedDiagSCMatrix::require_castdown(a,
                               "BlockedDiagSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    fprintf(stderr,"BlockedDiagSCMatrix:: accumulate(SCMatrix*a):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->accumulate(la->mats_[i].pointer());
}

double
BlockedDiagSCMatrix::invert_this()
{
  double det = 1.0;

  for (int i=0; i < d->nblocks(); i++)
    det *= mats_[i]->invert_this();
  
  return det;
}

double
BlockedDiagSCMatrix::determ_this()
{
  double det = 1.0;

  for (int i=0; i < d->nblocks(); i++)
    det *= mats_[i]->determ_this();
  
  return det;
}

double
BlockedDiagSCMatrix::trace()
{
  double det = 0;

  for (int i=0; i < d->nblocks(); i++)
    det += mats_[i]->trace();
  
  return det;
}

void
BlockedDiagSCMatrix::gen_invert_this()
{
  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->gen_invert_this();
}

void
BlockedDiagSCMatrix::element_op(const RefSCElementOp& op)
{
  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->element_op(op);
}

void
BlockedDiagSCMatrix::element_op(const RefSCElementOp2& op,
                              DiagSCMatrix* m)
{
  BlockedDiagSCMatrix *lm = BlockedDiagSCMatrix::require_castdown(m,
                                    "BlockedDiagSCMatrix::element_op");
  if (!dim()->equiv(lm->dim())) {
    fprintf(stderr,"BlockedDiagSCMatrix: bad element_op\n");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->element_op(op,lm->mats_[i].pointer());
}

void
BlockedDiagSCMatrix::element_op(const RefSCElementOp3& op,
                              DiagSCMatrix* m,DiagSCMatrix* n)
{
  BlockedDiagSCMatrix *lm = BlockedDiagSCMatrix::require_castdown(m,
                                      "BlockedDiagSCMatrix::element_op");
  BlockedDiagSCMatrix *ln = BlockedDiagSCMatrix::require_castdown(n,
                                      "BlockedDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
    fprintf(stderr,"BlockedDiagSCMatrix: bad element_op\n");
    abort();
  }

  for (int i=0; i < d->nblocks(); i++)
    mats_[i]->element_op(op,lm->mats_[i].pointer(),ln->mats_[i].pointer());
}

void
BlockedDiagSCMatrix::print(const char *title, ostream& os, int prec)
{
  int len = (title) ? strlen(title) : 0;
  char *newtitle = new char[len + 80];

  for (int i=0; i < d->nblocks(); i++) {
    sprintf(newtitle,"%s:  block %d",title,i+1);
    mats_[i]->print(newtitle, os, prec);
  }

  delete[] newtitle;
}
