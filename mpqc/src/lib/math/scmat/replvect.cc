
#include <iostream.h>
#include <iomanip.h>

#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/repl.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// ReplSCVector member functions

#define CLASSNAME ReplSCVector
#define PARENTS public SCVector
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
ReplSCVector::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCVector::_castdown(cd);
  return do_castdowns(casts,cd);
}

ReplSCVector::ReplSCVector(const RefSCDimension&a,ReplSCMatrixKit*k):
  SCVector(a,k)
{
  vector = new double[a->n()];
  init_blocklist();
}

void
ReplSCVector::before_elemop()
{
  // zero out the blocks not in my block list
  int i;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  for (i=0; i<d->blocks()->nblock(); i++) {
      if (i%nproc == me) continue;
      memset(&vector[d->blocks()->start(i)], 0,
             sizeof(double)*(d->blocks()->fence(i)
                             - d->blocks()->start(i)));
    }
}

void
ReplSCVector::after_elemop()
{
  messagegrp()->sum(vector, d->n());
}

void
ReplSCVector::init_blocklist()
{
  int i;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  for (i=0; i<d->blocks()->nblock(); i++) {
      if (i%nproc != me) continue;
      blocklist->insert(
          new SCVectorSimpleSubBlock(d->blocks()->start(i),
                                     d->blocks()->fence(i),
                                     d->blocks()->start(i),
                                     vector));
    }
}

ReplSCVector::~ReplSCVector()
{
}

double
ReplSCVector::get_element(int i)
{
  return vector[i];
}

void
ReplSCVector::set_element(int i,double a)
{
  vector[i] = a;
}

void
ReplSCVector::accumulate_element(int i,double a)
{
  vector[i] += a;
}

void
ReplSCVector::accumulate_product(SCMatrix*a,SCVector*b)
{
  const char* name = "ReplSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  ReplSCMatrix* la = ReplSCMatrix::require_castdown(a,name);
  ReplSCVector* lb = ReplSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->rowdim()) || !la->coldim()->equiv(lb->dim())) {
      cerr << indent
           << "ReplSCVector::accumulate_product(SCMatrix*a,SCVector*b): "
           << "dimensions don't match\n";
      abort();
    }

  cmat_mxm(la->rows, 0,
           &lb->vector, 1,
           &vector, 1,
           n(), la->ncol(), 1,
           1);
}

void
ReplSCVector::accumulate_product(SymmSCMatrix*a,SCVector*b)
{
  const char* name = "ReplSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  ReplSymmSCMatrix* la = ReplSymmSCMatrix::require_castdown(a,name);
  ReplSCVector* lb = ReplSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim()) || !la->dim()->equiv(lb->dim())) {
      cerr << indent
           << "ReplSCVector::accumulate_product(SymmSCMatrix*a,SCVector*b): "
           << "dimensions don't match\n";
      abort();
    }

  double** adat = la->rows;
  double* bdat = lb->vector;
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
      vector[i] += tmp;
    }
}

void
ReplSCVector::accumulate(SCVector*a)
{
  // make sure that the argument is of the correct type
  ReplSCVector* la
    = ReplSCVector::require_castdown(a,"ReplSCVector::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      cerr << indent << "ReplSCVector::accumulate(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) vector[i] += la->vector[i];
}

void
ReplSCVector::accumulate(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  ReplSCMatrix *la
    = ReplSCMatrix::require_castdown(a,"ReplSCVector::accumulate");

  // make sure that the dimensions match
  if (!((la->rowdim()->equiv(dim()) && la->coldim()->n() == 1)
        || (la->coldim()->equiv(dim()) && la->rowdim()->n() == 1))) {
      cerr << indent << "ReplSCVector::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) vector[i] += la->matrix[i];
}

void
ReplSCVector::assign(double a)
{
  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) vector[i] = a;
}

void
ReplSCVector::assign(SCVector*a)
{
  // make sure that the argument is of the correct type
  ReplSCVector* la
    = ReplSCVector::require_castdown(a,"ReplSCVector::assign");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      cerr << indent << "ReplSCVector::assign(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) vector[i] = la->vector[i];
}

void
ReplSCVector::assign(const double*a)
{
  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) vector[i] = a[i];
}

double
ReplSCVector::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  ReplSCVector* la
    = ReplSCVector::require_castdown(a,"ReplSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      cerr << indent << "ReplSCVector::scalar_product(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  double result = 0.0;
  for (i=0; i<nelem; i++) result += vector[i] * la->vector[i];
  return result;
}

void
ReplSCVector::element_op(const RefSCElementOp& op)
{
  if (op->has_side_effects()) before_elemop();
  SCMatrixBlockListIter i;
  for (i = blocklist->begin(); i != blocklist->end(); i++) {
      op->process_base(i.block());
    }
  if (op->has_side_effects()) after_elemop();
  if (op->has_collect()) op->collect(messagegrp());
}

void
ReplSCVector::element_op(const RefSCElementOp2& op,
                          SCVector* m)
{
  ReplSCVector *lm
      = ReplSCVector::require_castdown(m, "ReplSCVector::element_op");

  if (!dim()->equiv(lm->dim())) {
      cerr << indent << "ReplSCVector: bad element_op\n";
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
  if (op->has_collect()) op->collect(messagegrp());
}

void
ReplSCVector::element_op(const RefSCElementOp3& op,
                          SCVector* m,SCVector* n)
{
  ReplSCVector *lm
      = ReplSCVector::require_castdown(m, "ReplSCVector::element_op");
  ReplSCVector *ln
      = ReplSCVector::require_castdown(n, "ReplSCVector::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      cerr << indent << "ReplSCVector: bad element_op\n";
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
  if (op->has_collect()) op->collect(messagegrp());
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
ReplSCVector::vprint(const char *title, ostream& os, int prec)
{
  int i;
  int lwidth;
  double max=this->maxabs();

  if (messagegrp()->me() != 0) return;

  max = (max==0.0) ? 1.0 : log10(max);
  if (max < 0.0) max=1.0;

  lwidth = prec + 5 + (int) max;

  os.setf(ios::fixed,ios::floatfield); os.precision(prec);
  os.setf(ios::right,ios::adjustfield);

  if (title)
    os << endl << indent << title << endl;
  else
    os << endl;

  if (n()==0) {
    os << indent << "empty vector\n";
    return;
  }

  for (i=0; i < n(); i++)
    os << indent << setw(5) << i+1 << setw(lwidth) << vector[i] << endl;
  os << endl;

  os.flush();
}

RefSCMatrixSubblockIter
ReplSCVector::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new ReplSCMatrixListSubblockIter(access, blocklist,
                                          messagegrp(),
                                          vector, d->n());
}

RefSCMatrixSubblockIter
ReplSCVector::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      cerr << indent << "ReplSCVector::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  RefSCMatrixBlockList allblocklist = new SCMatrixBlockList();
  allblocklist->insert(new SCVectorSimpleSubBlock(0, d->n(), 0, vector));
  return new ReplSCMatrixListSubblockIter(access, allblocklist,
                                          messagegrp(),
                                          vector, d->n());
}

RefReplSCMatrixKit
ReplSCVector::skit()
{
  return ReplSCMatrixKit::castdown(kit().pointer());
}
