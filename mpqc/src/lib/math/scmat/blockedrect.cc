
#include <math.h>

#include <util/misc/formio.h>
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
#include <util/class/classi.h>
void *
BlockedSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
BlockedSCMatrix::resize(SCDimension *a, SCDimension *b)
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
  
  d1 = a;
  d2 = b;
  
  if (!a || !b || !a->blocks()->nblock() || !b->blocks()->nblock())
    return;

  if (a->blocks()->nblock() > 1 && b->blocks()->nblock() == 1) {
    nblocks_ = d1->blocks()->nblock();
    mats_ = new RefSCMatrix[d1->blocks()->nblock()];
    for (int i=0; i < d1->blocks()->nblock(); i++)
      if (d1->blocks()->size(i) && d2->blocks()->size(0))
        mats_[i] = subkit->matrix(d1->blocks()->subdim(i),
                                  d2->blocks()->subdim(0));

  } else if (a->blocks()->nblock() == 1 && b->blocks()->nblock() > 1) {
    nblocks_ = d2->blocks()->nblock();
    mats_ = new RefSCMatrix[d2->blocks()->nblock()];
    for (int i=0; i < d2->blocks()->nblock(); i++)
      if (d2->blocks()->size(i) && d1->blocks()->size(0))
        mats_[i] = subkit->matrix(d1->blocks()->subdim(0),
                                  d2->blocks()->subdim(i));

  } else if (a->blocks()->nblock() == b->blocks()->nblock()) {
    nblocks_ = d2->blocks()->nblock();
    mats_ = new RefSCMatrix[d1->blocks()->nblock()];
    for (int i=0; i < d1->blocks()->nblock(); i++)
      if (d2->blocks()->size(i) && d1->blocks()->size(i))
        mats_[i] = subkit->matrix(d1->blocks()->subdim(i),
                                  d2->blocks()->subdim(i));

  } else {
    cerr << indent << "BlockedSCMatrix::resize: wrong number of blocks\n";
    abort();
  }

}

BlockedSCMatrix::BlockedSCMatrix(const RefSCDimension&a,
                                 const RefSCDimension&b,
                                 BlockedSCMatrixKit*k):
  SCMatrix(a,b,k),
  subkit(k->subkit()),
  mats_(0), nblocks_(0)
{
  resize(a,b);
}

BlockedSCMatrix::~BlockedSCMatrix()
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
  nblocks_=0;
}

void
BlockedSCMatrix::assign(double v)
{
  for (int i=0; i < nblocks_; i++)
    if (mats_[i].nonnull())
      mats_[i]->assign(v);
}

double
BlockedSCMatrix::get_element(int i,int j)
{
  int block_i, block_j;
  int elem_i, elem_j;

  d1->blocks()->elem_to_block(i,block_i,elem_i);
  d2->blocks()->elem_to_block(j,block_j,elem_j);

  if (d1->blocks()->nblock() == 1 && d2->blocks()->nblock() > 1) {
    return mats_[block_j]->get_element(elem_i,elem_j);

  } else if (d1->blocks()->nblock() > 1 && d2->blocks()->nblock() == 1 ||
             d1->blocks()->nblock() == d2->blocks()->nblock()) {
    return mats_[block_i]->get_element(elem_i,elem_j);

  } else {
    return 0;
  }
}

void
BlockedSCMatrix::set_element(int i,int j,double a)
{
  int block_i, block_j;
  int elem_i, elem_j;

  d1->blocks()->elem_to_block(i,block_i,elem_i);
  d2->blocks()->elem_to_block(j,block_j,elem_j);
  
  if (d1->blocks()->nblock() == 1 && d2->blocks()->nblock() > 1) {
    mats_[block_j]->set_element(elem_i,elem_j,a);

  } else if (d1->blocks()->nblock() > 1 && d2->blocks()->nblock() == 1 ||
             d1->blocks()->nblock() == d2->blocks()->nblock()) {
    mats_[block_i]->set_element(elem_i,elem_j,a);
  }
}

void
BlockedSCMatrix::accumulate_element(int i,int j,double a)
{
  int block_i, block_j;
  int elem_i, elem_j;

  d1->blocks()->elem_to_block(i,block_i,elem_i);
  d2->blocks()->elem_to_block(j,block_j,elem_j);
  
  if (d1->blocks()->nblock() == 1 && d2->blocks()->nblock() > 1) {
    mats_[block_j]->accumulate_element(elem_i,elem_j,a);

  } else if (d1->blocks()->nblock() > 1 && d2->blocks()->nblock() == 1 ||
             d1->blocks()->nblock() == d2->blocks()->nblock()) {
    mats_[block_i]->accumulate_element(elem_i,elem_j,a);
  }
}

SCMatrix *
BlockedSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  cerr << indent << "BlockedSCMatrix::get_subblock: cannot get subblock\n";
  abort();
  return 0;
}

void
BlockedSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                               int source_br, int source_bc)
{
  cerr << indent
       << "BlockedSCMatrix::assign_subblock: cannot assign subblock\n";
  abort();
}

void
BlockedSCMatrix::accumulate_subblock(SCMatrix*sb,
                                     int br, int er, int bc, int ec,
                                     int source_br, int source_bc)
{
  cerr << indent << "BlockedSCMatrix::accumulate_subblock:"
       << " cannot accumulate subblock\n";
  abort();
}

SCVector *
BlockedSCMatrix::get_row(int i)
{
  cerr << indent << "BlockedSCMatrix::get_row: cannot get row\n";
  abort();

  return 0;
}

void
BlockedSCMatrix::assign_row(SCVector *v, int i)
{
  cerr << indent << "BlockedSCMatrix::assign_row: cannot assign row\n";
  abort();
}

void
BlockedSCMatrix::accumulate_row(SCVector *v, int i)
{
  cerr << indent << "BlockedSCMatrix::accumulate_row: cannot accumulate row\n";
  abort();
}

SCVector *
BlockedSCMatrix::get_column(int i)
{
  cerr << indent << "BlockedSCMatrix::get_column: cannot get column\n";
  abort();

  return 0;
}

void
BlockedSCMatrix::assign_column(SCVector *v, int i)
{
  cerr << indent << "BlockedSCMatrix::assign_column: cannot assign column\n";
  abort();
}

void
BlockedSCMatrix::accumulate_column(SCVector *v, int i)
{
  cerr << indent
       << "BlockedSCMatrix::accumulate_column: cannot accumulate column\n";
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
    cerr << indent
         << "BlockedSCMatrix::accumulate_outer_product(SCVector*,SCVector*): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d1->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->accumulate_outer_product(la->vecs_[i], lb->vecs_[i]);
}

void
BlockedSCMatrix::accumulate_product(SCMatrix*a,SCMatrix*b)
{
  int i, zero = 0;

  const char* name = "BlockedSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = BlockedSCMatrix::require_castdown(a,name);
  BlockedSCMatrix* lb = BlockedSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->coldim()) ||
      !la->coldim()->equiv(lb->rowdim())) {
    cerr << indent
         << "BlockedSCMatrix::accumulate_product(SCMatrix*a,SCMatrix*b): "
         << "dimensions don't match\n";
    abort();
  }

  // find out the number of blocks we need to process.
  int mxnb = (nblocks_ > la->nblocks_) ? nblocks_ : la->nblocks_;
  
  int nrba = la->d1->blocks()->nblock();
  int ncba = la->d2->blocks()->nblock();
  int nrbb = lb->d1->blocks()->nblock();
  int ncbb = lb->d2->blocks()->nblock();
  
  int &mi = (nrba==1 && ncba > 1 && nrbb > 1 && ncbb==1) ? zero : i;
  int &ai = (nrba==1 && ncba==1) ? zero : i;
  int &bi = (nrbb==1 && ncbb==1) ? zero : i;

  for (i=0; i < mxnb; i++) {
    if (mats_[mi].null() || la->mats_[ai].null() || lb->mats_[bi].null())
      continue;
    mats_[mi]->accumulate_product(la->mats_[ai], lb->mats_[bi]);
  }
}

void
BlockedSCMatrix::accumulate_product(SCMatrix*a,SymmSCMatrix*b)
{
  int i, zero=0;
  
  const char* name = "BlockedSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = BlockedSCMatrix::require_castdown(a,name);
  BlockedSymmSCMatrix* lb = BlockedSymmSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
    cerr << indent
         << "BlockedSCMatrix::accumulate_product(SCMatrix*a,SymmSCMatrix*b): "
         << "dimensions don't match\n";
    abort();
  }

  int &bi = (lb->d->blocks()->nblock()==1) ? zero : i;
  
  for (i=0; i < nblocks_; i++) {
    if (mats_[i].null() || la->mats_[i].null() || lb->mats_[bi].null())
      continue;
    mats_[i]->accumulate_product(la->mats_[i], lb->mats_[bi]);
  }
}


void
BlockedSCMatrix::accumulate_product(SCMatrix*a,DiagSCMatrix*b)
{
  int i, zero=0;
  
  const char* name = "BlockedSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = BlockedSCMatrix::require_castdown(a,name);
  BlockedDiagSCMatrix* lb = BlockedDiagSCMatrix::require_castdown(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
    cerr << indent
         << "BlockedSCMatrix::accumulate_product(SCMatrix*a,DiagSCMatrix*b): "
         << "dimensions don't match\n";
    abort();
  }

  int &bi = (lb->d->blocks()->nblock()==1) ? zero : i;
  
  for (i=0; i < nblocks_; i++) {
    if (mats_[i].null() || la->mats_[i].null() || lb->mats_[bi].null())
      continue;
    mats_[i]->accumulate_product(la->mats_[i], lb->mats_[bi]);
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
    cerr << indent << "BlockedSCMatrix::accumulate(SCMatrix*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i].nonnull())
      mats_[i]->accumulate(la->mats_[i]);
}

void
BlockedSCMatrix::accumulate(SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  BlockedSymmSCMatrix* la
    = BlockedSymmSCMatrix::require_castdown(a,"BlockedSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
    cerr << indent << "BlockedSCMatrix::accumulate(SymmSCMatrix*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i].nonnull())
      mats_[i]->accumulate(la->mats_[i]);
}

void
BlockedSCMatrix::accumulate(DiagSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  BlockedDiagSCMatrix* la
    = BlockedDiagSCMatrix::require_castdown(a,"BlockedSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
    cerr << indent << "BlockedSCMatrix::accumulate(DiagSCMatrix*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i].nonnull())
      mats_[i]->accumulate(la->mats_[i]);
}

void
BlockedSCMatrix::accumulate(SCVector*a)
{
  // make sure that the arguments is of the correct type
  BlockedSCVector* la
    = BlockedSCVector::require_castdown(a,"BlockedSCVector::accumulate");

  // make sure that the dimensions match
  if (!((rowdim()->equiv(la->dim()) && coldim()->n() == 1)
        || (coldim()->equiv(la->dim()) && rowdim()->n() == 1))) {
    cerr << indent << "BlockedSCMatrix::accumulate(SCVector*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i].nonnull())
      mats_[i]->accumulate(la->vecs_[i]);
}

void
BlockedSCMatrix::transpose_this()
{
  for (int i=0; i < nblocks_; i++)
    if (mats_[i].nonnull())
      mats_[i]->transpose_this();
  
  RefSCDimension tmp = d1;
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
  if (d1->blocks()->nblock() == d2->blocks()->nblock()) {
    for (i=0; i < nblocks_; i++)
      res *= mats_[i]->invert_this();
    return res;
  }

  // ok, let's make sure that the matrix is at least square
  if (d1->n() != d2->n()) {
    cerr << indent
         << "BlockedSCMatrix::invert_this: SVD not implemented yet\n";
    abort();
  }

  if (d1->blocks()->nblock() == 1) {
    RefSCMatrix tdim = subkit->matrix(d1->blocks()->subdim(0),
                                      d1->blocks()->subdim(0));
    tdim->convert(this);
    res = tdim->invert_this();
    transpose_this();

    // d1 and d2 were swapped by now
    for (i=0; i < d1->blocks()->nblock(); i++)
      if (mats_[i].nonnull())
        mats_[i]->convert(tdim.get_subblock(d1->blocks()->start(i),
                                            d1->blocks()->fence(i)-1,
                                            0, d2->n()-1));
    
    return res;

  } else if (d2->blocks()->nblock() == 1) {
    RefSCMatrix tdim = subkit->matrix(d2->blocks()->subdim(0),
                                      d2->blocks()->subdim(0));

    tdim->convert(this);
    res = tdim->invert_this();
    transpose_this();

    // d1 and d2 were swapped by now
    for (i=0; i < d2->blocks()->nblock(); i++)
      if (mats_[i].nonnull())
        mats_[i]->convert(tdim.get_subblock(0, d1->n()-1,
                                            d2->blocks()->start(i),
                                            d2->blocks()->fence(i)-1));
    
    return res;

  } else {
    cerr << indent
         << "BlockedSCMatrix::invert_this: SVD not implemented yet\n";
    abort();
  }

  return 0.0;
}

void
BlockedSCMatrix::gen_invert_this()
{
  cerr << indent
       << "BlockedSCMatrix::gen_invert_this: SVD not implemented yet\n";
  abort();
}

double
BlockedSCMatrix::determ_this()
{
  double res=1;
  
  for (int i=0; i < nblocks_; i++)
    if (mats_[i].nonnull())
      res *= mats_[i]->determ_this();

  return res;
}

double
BlockedSCMatrix::trace()
{
  double ret=0;
  for (int i=0; i < nblocks_; i++)
    if (mats_[i].nonnull())
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
    if (mats_[i].nonnull())
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
    cerr << indent << "BlockedSCMatrix::solve_this(SCVector*v): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i].nonnull())
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
    cerr << indent << "BlockedSCMatrix::schmidt_orthog(): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i].nonnull())
      mats_[i]->schmidt_orthog(lS->mats_[i].pointer(),
                               lS->dim()->blocks()->subdim(i).n());
}

void
BlockedSCMatrix::element_op(const RefSCElementOp& op)
{
  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op.pointer());

  for (int i=0; i < nblocks_; i++) {
    if (i < nblocks_-1)
      op->defer_collect(1);
    else
      op->defer_collect(0);
        
    if (bop)
      bop->working_on(i);
    if (mats_[i].nonnull())
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
    cerr << indent << "BlockedSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp2 *bop = BlockedSCElementOp2::castdown(op.pointer());

  for (int i=0; i < nblocks_; i++) {
    if (i < nblocks_-1)
      op->defer_collect(1);
    else
      op->defer_collect(0);
        
    if (bop)
      bop->working_on(i);
    if (mats_[i].nonnull())
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
    cerr << indent << "BlockedSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp3 *bop = BlockedSCElementOp3::castdown(op.pointer());

  for (int i=0; i < nblocks_; i++) {
    if (i < nblocks_-1)
      op->defer_collect(1);
    else
      op->defer_collect(0);
        
    if (bop)
      bop->working_on(i);
    if (mats_[i].nonnull())
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
    if (mats_[i].null())
      continue;

    sprintf(newtitle,"%s:  block %d",title,i+1);
    mats_[i]->print(newtitle, os, prec);
  }

  delete[] newtitle;
}

RefSCDimension
BlockedSCMatrix::rowdim(int i)
{
  return d1->blocks()->subdim(i);
}

RefSCDimension
BlockedSCMatrix::coldim(int i)
{
  return d2->blocks()->subdim(i);
}

int
BlockedSCMatrix::nblocks() const
{
  return nblocks_;
}

RefSCMatrix
BlockedSCMatrix::block(int i)
{
  return mats_[i];
}

RefSCMatrixSubblockIter
BlockedSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  RefSCMatrixCompositeSubblockIter iter
      = new SCMatrixCompositeSubblockIter(access, nblocks());
  for (int i=0; i<nblocks(); i++) {
      if (block(i).null())
          iter->set_iter(i, new SCMatrixNullSubblockIter(access));
      else
          iter->set_iter(i, block(i)->local_blocks(access));
    }
  RefSCMatrixSubblockIter ret = iter.pointer();
  return ret;
}

RefSCMatrixSubblockIter
BlockedSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  RefSCMatrixCompositeSubblockIter iter
      = new SCMatrixCompositeSubblockIter(access, nblocks());
  for (int i=0; i<nblocks(); i++) {
      if (block(i).null())
          iter->set_iter(i, new SCMatrixNullSubblockIter(access));
      else
          iter->set_iter(i, block(i)->all_blocks(access));
    }
  RefSCMatrixSubblockIter ret = iter.pointer();
  return ret;
}

void
BlockedSCMatrix::save(StateOut&s)
{
  int nr = nrow();
  int nc = ncol();
  s.put(nr);
  s.put(nc);
  int has_subblocks = 1;
  s.put(has_subblocks);
  s.put(nblocks());
  for (int i=0; i<nblocks(); i++) {
      block(i).save(s);
    }
}

void
BlockedSCMatrix::restore(StateIn& s)
{
  int nrt, nr = nrow();
  int nct, nc = ncol();
  s.get(nrt);
  s.get(nct);
  if (nrt != nr || nct != nc) {
      cerr << indent << "BlockedSCMatrix::restore(): bad dimensions" << endl;
      abort();
    }
  int has_subblocks;
  s.get(has_subblocks);
  if (has_subblocks) {
      int nblock;
      s.get(nblock);
      if (nblock != nblocks()) {
          cerr << indent
               << "BlockedSCMatrix::restore(): nblock differs\n" << endl;
          abort();
        }
      for (int i=0; i<nblocks(); i++) {
          block(i).restore(s);
        }
    }
  else {
      cerr << indent
           << "BlockedSCMatrix::restore(): no subblocks--cannot restore"
           << endl;
      abort();
    }
}
