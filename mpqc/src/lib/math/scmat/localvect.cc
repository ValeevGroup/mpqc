
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// LocalSCVector member functions

#define CLASSNAME LocalSCVector
#define PARENTS public SCVector
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LocalSCVector::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCVector::_castdown(cd);
  return do_castdowns(casts,cd);
}

LocalSCVector::LocalSCVector()
{
}

LocalSCVector::LocalSCVector(LocalSCDimension*a):
  d(a)
{
  resize(a->n());
}

LocalSCVector::LocalSCVector(StateIn&s):
  SCVector(s)
{
  d.restore_state(s);
  block.restore_state(s);
}

LocalSCVector::LocalSCVector(const RefKeyVal&keyval)
{
  d = keyval->describedclassvalue("dim");
  d.require_nonnull();
  block = new SCVectorSimpleBlock(0,d->n());
  int i;
  for (i=0; i<n(); i++) {
      set_element(i,keyval->doublevalue("data",i,i));
    }
}

void
LocalSCVector::resize(int n)
{
  block = new SCVectorSimpleBlock(0,n);
}

void
LocalSCVector::save_data_state(StateOut&s)
{
  SCVector::save_data_state(s);
  d.save_state(s);
  block.save_state(s);
}

LocalSCVector::~LocalSCVector()
{
}

RefSCDimension
LocalSCVector::dim()
{
  return d;
}

double
LocalSCVector::get_element(int i)
{
  int size = block->iend - block->istart;
  if (i < 0 || i >= size) {
      fprintf(stderr,"LocalSCVector::get_element: bad index\n");
      abort();
    }
  return block->data[i];
}

void
LocalSCVector::set_element(int i,double a)
{
  int size = block->iend - block->istart;
  if (i < 0 || i >= size) {
      fprintf(stderr,"LocalSCVector::set_element: bad index\n");
      abort();
    }
  block->data[i] = a;
}

void
LocalSCVector::accumulate_element(int i,double a)
{
  int size = block->iend - block->istart;
  if (i < 0 || i >= size) {
      fprintf(stderr,"LocalSCVector::accumulate_element: bad index\n");
      abort();
    }
  block->data[i] += a;
}

void
LocalSCVector::accumulate_product(SCMatrix*a,SCVector*b)
{
  const char* name = "LocalSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = LocalSCMatrix::require_castdown(a,name);
  LocalSCVector* lb = LocalSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->rowdim()) || !la->coldim()->equiv(lb->dim())) {
      fprintf(stderr,"LocalSCVector::"
              "accumulate_product(SCMatrix*a,SCVector*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  cmat_mxm(la->rows, 0,
           &(lb->block->data), 1,
           &(block->data), 1,
           n(), la->ncol(), 1,
           1);
}

void
LocalSCVector::accumulate_product(SymmSCMatrix*a,SCVector*b)
{
  const char* name = "LocalSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSymmSCMatrix* la = LocalSymmSCMatrix::require_castdown(a,name);
  LocalSCVector* lb = LocalSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim()) || !la->dim()->equiv(lb->dim())) {
      fprintf(stderr,"LocalSCVector::"
              "accumulate_product(SymmSCMatrix*a,SCVector*b):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  double* thisdat = block->data;
  double** adat = la->rows;
  double* bdat = lb->block->data;
  double tmp;
  int n = dim()->n();
  int i, j;
  for (i=0; i<n; i++) {
      tmp = 0.0;
      for (j=0; j<=i; j++) {
          tmp += adat[i][j] * bdat[j];
        }
      for (; j<n; j++) {
          tmp += adat[j][i] * bdat[j];
        }
      thisdat[i] += tmp;
    }
}

void
LocalSCVector::accumulate(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      fprintf(stderr,"LocalSCVector::accumulate(SCVector*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

void
LocalSCVector::assign(double a)
{
  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) block->data[i] = a;
}

void
LocalSCVector::assign(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::assign");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      fprintf(stderr,"LocalSCVector::assign(SCVector*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) block->data[i] = la->block->data[i];
}

void
LocalSCVector::assign(const double*a)
{
  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) block->data[i] = a[i];
}

double
LocalSCVector::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      fprintf(stderr,"LocalSCVector::scalar_product(SCVector*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = d->n();
  int i;
  double result = 0.0;
  for (i=0; i<nelem; i++) result += block->data[i] * la->block->data[i];
  return result;
}

void
LocalSCVector::element_op(const RefSCElementOp& op)
{
  op->process(block.pointer());
}

void
LocalSCVector::element_op(const RefSCElementOp2& op,
                          SCVector* m)
{
  LocalSCVector *lm
      = LocalSCVector::require_castdown(m, "LocalSCVector::element_op");

  if (!dim()->equiv(lm->dim())) {
      fprintf(stderr,"LocalSCVector: bad element_op\n");
      abort();
    }
  op->process(block.pointer(), lm->block.pointer());
}

void
LocalSCVector::element_op(const RefSCElementOp3& op,
                          SCVector* m,SCVector* n)
{
  LocalSCVector *lm
      = LocalSCVector::require_castdown(m, "LocalSCVector::element_op");
  LocalSCVector *ln
      = LocalSCVector::require_castdown(n, "LocalSCVector::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      fprintf(stderr,"LocalSCVector: bad element_op\n");
      abort();
    }
  op->process(block.pointer(), lm->block.pointer(), ln->block.pointer());
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
LocalSCVector::print(const char *title, ostream& os, int prec)
{
  int i;
  int lwidth;
  double max=this->maxabs();

  max=(max==0.0)?1.0:log10(max);
  if(max < 0.0) max=1.0;

  lwidth = prec+5+(int) max;

  os.setf(ios::fixed,ios::floatfield); os.precision(prec);
  os.setf(ios::right,ios::adjustfield);

  if(title) os << "\n" << title << "\n";
  else os << "\n";

  if(n()==0) { os << " empty vector\n"; return; }

  for (i=0; i<n(); i++) {
      os.width(5); os << i+1;
      os.width(lwidth); os << block->data[i];
      os << "\n";
    }
  os << "\n";

  os.flush();
}

RefSCMatrixSubblockIter
LocalSCVector::local_blocks(SCMatrixSubblockIter::Access access)
{
  RefSCMatrixSubblockIter iter
      = new SCMatrixSimpleSubblockIter(access, block.pointer());
  return iter;
}

RefSCMatrixSubblockIter
LocalSCVector::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      cerr << "LocalVectorSCMatrix::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  return local_blocks(access);
}
