
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

#include <math/scmat/local.h>
#include <math/scmat/repl.h>

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

  if (a->nblocks() > 1 && b->nblocks() == 1) {
    nblocks_ = d1->nblocks();
    mats_ = new RefSCMatrix[d1->nblocks()];
    for (int i=0; i < d1->nblocks(); i++)
      mats_[i] = d1->dim(i)->create_matrix(d2->dim(0));

  } else if (a->nblocks() == 1 && b->nblocks() > 1) {
    nblocks_ = d2->nblocks();
    mats_ = new RefSCMatrix[d2->nblocks()];
    for (int i=0; i < d2->nblocks(); i++)
      mats_[i] = d1->dim(0)->create_matrix(d2->dim(i));

  } else if (a->nblocks() == b->nblocks()) {
    nblocks_ = d2->nblocks();
    mats_ = new RefSCMatrix[d1->nblocks()];
    for (int i=0; i < d1->nblocks(); i++)
      mats_[i] = d1->dim(i)->create_matrix(d2->dim(i).pointer());

  } else {
    fprintf(stderr,"BlockedSCMatrix::resize: wrong number of blocks\n");
    abort();
  }

}


BlockedSCMatrix::BlockedSCMatrix() :
  mats_(0), nblocks_(0)
{
}

BlockedSCMatrix::BlockedSCMatrix(BlockedSCDimension*a,BlockedSCDimension*b):
  mats_(0), nblocks_(0)
{
  resize(a,b);
}

BlockedSCMatrix::BlockedSCMatrix(StateIn&s):
  SCMatrix(s)
{
  d1.restore_state(s);
  d2.restore_state(s);
  s.get(nblocks_);
  mats_ = new RefSCMatrix[nblocks_];
  for (int i=0; i < nblocks_; i++)
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
  s.put(nblocks_);
  for (int i=0; i < nblocks_; i++)
    mats_[i].save_state(s);
}

BlockedSCMatrix::~BlockedSCMatrix()
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
  nblocks_=0;
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
  for (int i=0; i < nblocks_; i++)
    mats_[i]->assign(v);
}

double
BlockedSCMatrix::get_element(int i,int j)
{
  int block_i = d1->block(i);
  int block_j = d2->block(j);

  int elem_i = i - d1->first(block_i);
  int elem_j = j - d2->first(block_j);

  if (d1->nblocks() == 1 && d2->nblocks() > 1) {
    return mats_[block_j]->get_element(elem_i,elem_j);

  } else if (d1->nblocks() > 1 && d2->nblocks() == 1 ||
             d1->nblocks() == d2->nblocks()) {
    return mats_[block_i]->get_element(elem_i,elem_j);

  } else {
    return 0;
  }
}

void
BlockedSCMatrix::set_element(int i,int j,double a)
{
  int block_i = d1->block(i);
  int block_j = d2->block(j);

  int elem_i = i - d1->first(block_i);
  int elem_j = j - d2->first(block_j);
  
  if (d1->nblocks() == 1 && d2->nblocks() > 1) {
    mats_[block_j]->set_element(elem_i,elem_j,a);

  } else if (d1->nblocks() > 1 && d2->nblocks() == 1 ||
             d1->nblocks() == d2->nblocks()) {
    mats_[block_i]->set_element(elem_i,elem_j,a);
  }
}

void
BlockedSCMatrix::accumulate_element(int i,int j,double a)
{
  int block_i = d1->block(i);
  int block_j = d2->block(j);

  int elem_i = i - d1->first(block_i);
  int elem_j = j - d2->first(block_j);
  
  if (d1->nblocks() == 1 && d2->nblocks() > 1) {
    mats_[block_j]->accumulate_element(elem_i,elem_j,a);

  } else if (d1->nblocks() > 1 && d2->nblocks() == 1 ||
             d1->nblocks() == d2->nblocks()) {
    mats_[block_i]->accumulate_element(elem_i,elem_j,a);
  }
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
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(lb->dim())) {
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
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->coldim()) ||
      !la->coldim()->equiv(lb->rowdim())) {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate_product(SCMatrix*a,SCMatrix*b):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  int nba = la->d1->nblocks();
  int nbb = lb->d2->nblocks();
  
  if (nba == 1 && nbb > 1) {
    for (int i=0; i < d2->nblocks(); i++)
      mats_[i]->accumulate_product(la->mats_[0], lb->mats_[i]);

  } else if (nba > 1 && nbb == 1) {
    for (int i=0; i < d1->nblocks(); i++)
      mats_[i]->accumulate_product(la->mats_[i], lb->mats_[0]);

  } else if (nba == nbb) {
    for (int i=0; i < d1->nblocks(); i++)
      mats_[i]->accumulate_product(la->mats_[i], lb->mats_[i]);

  } else {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate_product(SCMatrix*a,SCMatrix*b):\n");
    fprintf(stderr,"blocks don't match\n");
    abort();
  }
}

void
BlockedSCMatrix::accumulate_product(SCMatrix*a,SymmSCMatrix*b)
{
  const char* name = "BlockedSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = BlockedSCMatrix::require_castdown(a,name);
  BlockedSymmSCMatrix* lb = BlockedSymmSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate_product(SCMatrix*a,SymmSCMatrix*b):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  int nba = la->d1->nblocks();
  int nbb = lb->d->nblocks();
  
  if (nba == 1 && nbb > 1) {
    for (int i=0; i < d2->nblocks(); i++)
      mats_[i]->accumulate_product(la->mats_[0], lb->mats_[i]);

  } else if (nba > 1 && nbb == 1) {
    for (int i=0; i < d1->nblocks(); i++)
      mats_[i]->accumulate_product(la->mats_[i], lb->mats_[0]);

  } else if (nba == nbb) {
    for (int i=0; i < d1->nblocks(); i++)
      mats_[i]->accumulate_product(la->mats_[i], lb->mats_[i]);

  } else {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate_product(SCMatrix*a,SymmSCMatrix*b):\n");
    fprintf(stderr,"blocks don't match\n");
    abort();
  }
}

void
BlockedSCMatrix::accumulate_product(SCMatrix*a,DiagSCMatrix*b)
{
  const char* name = "BlockedSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = BlockedSCMatrix::require_castdown(a,name);
  BlockedDiagSCMatrix* lb = BlockedDiagSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate_product(SCMatrix*a,DiagSCMatrix*b):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  int nba = la->d1->nblocks();
  int nbb = lb->d->nblocks();
  
  if (nba == 1 && nbb > 1) {
    for (int i=0; i < d2->nblocks(); i++)
      mats_[i]->accumulate_product(la->mats_[0], lb->mats_[i]);

  } else if (nba > 1 && nbb == 1) {
    for (int i=0; i < d1->nblocks(); i++)
      mats_[i]->accumulate_product(la->mats_[i], lb->mats_[0]);

  } else if (nba == nbb) {
    for (int i=0; i < d1->nblocks(); i++)
      mats_[i]->accumulate_product(la->mats_[i], lb->mats_[i]);

  } else {
    fprintf(stderr,"BlockedSCMatrix::"
            "accumulate_product(SCMatrix*a,DiagSCMatrix*b):\n");
    fprintf(stderr,"blocks don't match\n");
    abort();
  }
}

void
BlockedSCMatrix::accumulate(SCMatrix*a)
{
  // make sure that the arguments is of the correct type
  BlockedSCMatrix* la
    = BlockedSCMatrix::require_castdown(a,"BlockedSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(la->coldim())) {
    fprintf(stderr,"BlockedSCMatrix::accumulate(SCMatrix*a):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    mats_[i]->accumulate(la->mats_[i]);
}

void
BlockedSCMatrix::transpose_this()
{
  for (int i=0; i < nblocks_; i++)
    mats_[i]->transpose_this();
  
  RefBlockedSCDimension tmp = d1;
  d1 = d2;
  d2 = tmp;
}

// hack, hack, hack!  One day we'll get svd working everywhere.
double
BlockedSCMatrix::invert_this()
{
  int i;
  double res=1;

  // if this matrix is block diagonal, then give a normal inversion a shot
  if (d1->nblocks() == d2->nblocks()) {
    for (i=0; i < nblocks_; i++)
      res *= mats_[i]->invert_this();
    return res;
  }

  // ok, let's make sure that the matrix is at least square
  if (d1->n() != d2->n()) {
    fprintf(stderr,"BlockedSCMatrix::invert_this: SVD not implemented yet");
    abort();
  }

  if (d1->nblocks() == 1) {
    RefSCMatrix tdim = d1->dim(0)->create_matrix(d1->dim(0));

    for (i=0; i < d2->nblocks(); i++)
      tdim.assign_subblock(mats_[i], 0, d1->n()-1,
                                     d2->first(i), d2->last(i)-1);

    res = tdim->invert_this();
    transpose_this();

    for (i=0; i < d1->nblocks(); i++)
      mats_[i].assign(tdim.get_subblock(d1->first(i), d1->last(i)-1,
                                        0, d2->n()-1));
    
    return res;

  } else if (d2->nblocks() == 1) {
    RefSCMatrix tdim = d2->dim(0)->create_matrix(d2->dim(0));

    for (i=0; i < d1->nblocks(); i++)
      tdim.assign_subblock(mats_[i], d1->first(i), d1->last(i)-1,
                                     0, d2->n()-1);

    res = tdim->invert_this();
    transpose_this();

    for (i=0; i < d2->nblocks(); i++)
      mats_[i].assign(tdim.get_subblock(0, d1->n()-1,
                                        d2->first(i), d2->last(i)-1));
    
    return res;

  } else {
    fprintf(stderr,"BlockedSCMatrix::invert_this: SVD not implemented yet");
    abort();
  }
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
  
  for (int i=0; i < nblocks_; i++)
    res *= mats_[i]->determ_this();

  return res;
}

double
BlockedSCMatrix::trace()
{
  double ret=0;
  for (int i=0; i < nblocks_; i++)
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

  for (int i=0; i < nblocks_; i++)
    mats_[i]->svd_this(lU->mats_[i], lsigma->mats_[i], lV->mats_[i]);
}

double
BlockedSCMatrix::solve_this(SCVector*v)
{
  double res=1;
  
  BlockedSCVector* lv =
    BlockedSCVector::require_castdown(v,"BlockedSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lv->dim())) {
    fprintf(stderr,"BlockedSCMatrix::solve_this(SCVector*v):\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    res *= mats_[i]->solve_this(lv->vecs_[i]);

  return res;
}

void
BlockedSCMatrix::schmidt_orthog(SymmSCMatrix *S, int nc)
{
  BlockedSymmSCMatrix* lS =
    BlockedSymmSCMatrix::require_castdown(S,"BlockedSCMatrix::schmidt_orthog");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lS->dim())) {
    fprintf(stderr,"BlockedSCMatrix::schmidt_orthog():\n");
    fprintf(stderr,"dimensions don't match\n");
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    mats_[i]->schmidt_orthog(lS->mats_[i].pointer(), lS->dim(i).n());
}

void
BlockedSCMatrix::element_op(const RefSCElementOp& op)
{
  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op);

  for (int i=0; i < nblocks_; i++) {
    if (bop)
      bop->working_on(i);
    mats_[i]->element_op(op);
  }
}

void
BlockedSCMatrix::element_op(const RefSCElementOp2& op,
                          SCMatrix* m)
{
  BlockedSCMatrix *lm
    = BlockedSCMatrix::require_castdown(m,"BlockedSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim())) {
    fprintf(stderr,"BlockedSCMatrix: bad element_op\n");
    abort();
  }

  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op);

  for (int i=0; i < nblocks_; i++) {
    if (bop)
      bop->working_on(i);
    mats_[i]->element_op(op,lm->mats_[i].pointer());
  }
}

void
BlockedSCMatrix::element_op(const RefSCElementOp3& op,
                          SCMatrix* m,SCMatrix* n)
{
  BlockedSCMatrix *lm
    = BlockedSCMatrix::require_castdown(m,"BlockedSCMatrix::element_op");
  BlockedSCMatrix *ln
    = BlockedSCMatrix::require_castdown(n,"BlockedSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim()) ||
      !rowdim()->equiv(ln->rowdim()) || !coldim()->equiv(ln->coldim())) {
    fprintf(stderr,"BlockedSCMatrix: bad element_op\n");
    abort();
  }

  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op);

  for (int i=0; i < nblocks_; i++) {
    if (bop)
      bop->working_on(i);
    mats_[i]->element_op(op,lm->mats_[i].pointer(),
                            ln->mats_[i].pointer());
  }
}

void
BlockedSCMatrix::print(const char *title, ostream& os, int prec)
{
  int len = (title) ? strlen(title) : 0;
  char *newtitle = new char[len + 80];

  for (int i=0; i < nblocks_; i++) {
    sprintf(newtitle,"%s:  block %d",title,i+1);
    mats_[i]->print(newtitle, os, prec);
  }

  delete[] newtitle;
}

RefSCMatrix
BlockedSCMatrix::unblock() const
{
  int i;
  
  int nrow = d1->n();
  int ncol = d2->n();
  
  RefSCDimension allrow, allcol;
  
  // ugly hack
  if (LocalSCDimension::castdown(d1) && LocalSCDimension::castdown(d2)) {
    allrow = new LocalSCDimension(nrow);
    allcol = new LocalSCDimension(ncol);
    
  } else if (ReplSCDimension::castdown(d1) && ReplSCDimension::castdown(d2)) {
    ReplSCDimension *foo = ReplSCDimension::castdown(d1);
    allrow = new ReplSCDimension(nrow, foo->messagegrp());
    allcol = new ReplSCDimension(ncol, foo->messagegrp());

  } else {
    fprintf(stderr,"BlockedSCMatrix::unblock: cannot unblock this type\n");
    abort();
  }
  
  RefSCMatrix ret = allrow->create_matrix(allcol);
  ret.assign(0.0);

  if (d1->nblocks()==d2->nblocks()) {
    for (i=0; i < d1->nblocks(); i++)
      ret.assign_subblock(mats_[i], d1->first(i), d1->last(i)-1,
                                    d2->first(i), d2->last(i)-1);
  } else if (d1->nblocks()==1) {
    for (i=0; i < d2->nblocks(); i++)
      ret.assign_subblock(mats_[i], 0, d1->n()-1,
                                    d2->first(i), d2->last(i)-1);

  } else if (d2->nblocks()==1) {
    for (i=0; i < d1->nblocks(); i++)
      ret.assign_subblock(mats_[i], d1->first(i), d1->last(i)-1,
                                    0, d2->n()-1);
  } else {
    fprintf(stderr,"BlockedSCMatrix::unblock: cannot unblock this type\n");
    abort();
  }
  
  return ret;
}

void
BlockedSCMatrix::block(const RefSCMatrix& m)
{
  int i;
  
  if (m->rowdim()->n() != rowdim()->n() ||
      m->coldim()->n() != coldim()->n()) {
    fprintf(stderr,"BlockedSCMatrix::block: wrong dimensions\n");
    abort();
  }
  
  if (d1->nblocks()==d2->nblocks()) {
    for (i=0; i < d1->nblocks(); i++)
      mats_[i].assign(m->get_subblock(d1->first(i), d1->last(i)-1,
                                      d2->first(i), d2->last(i)-1));
  } else if (d1->nblocks()==1) {
    for (i=0; i < d2->nblocks(); i++)
      mats_[i].assign(m->get_subblock(0, d1->n()-1,
                                      d2->first(i), d2->last(i)-1));

  } else if (d2->nblocks()==1) {
    for (i=0; i < d1->nblocks(); i++)
      mats_[i].assign(m->get_subblock(d1->first(i), d1->last(i)-1,
                                      0, d2->n()-1));
  } else {
    fprintf(stderr,"BlockedSCMatrix::block: cannot block this matrix\n");
    abort();
  }
}
