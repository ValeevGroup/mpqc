
#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/repl.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// ReplDiagSCMatrix member functions

#define CLASSNAME ReplDiagSCMatrix
#define PARENTS public DiagSCMatrix
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
ReplDiagSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DiagSCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

ReplDiagSCMatrix::ReplDiagSCMatrix(ReplSCDimension*a):
  d(a)
{
  matrix = new double[a->n()];
  init_blocklist();
}

void
ReplDiagSCMatrix::before_elemop()
{
  // zero out the blocks not in my block list
  int i;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  for (i=0; i<d->nblock(); i++) {
      if (i%nproc == me) continue;
      memset(&matrix[d->blockstart(i)], 0,
             sizeof(double)*(d->blockfence(i) - d->blockstart(i)));
    }
}

void
ReplDiagSCMatrix::after_elemop()
{
  messagegrp()->sum(matrix, d->n());
}

void
ReplDiagSCMatrix::init_blocklist()
{
  int i;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  for (i=0; i<d->nblock(); i++) {
      if (i%nproc != me) continue;
      blocklist->insert(new SCMatrixDiagSubBlock(d->blockstart(i),
                                                 d->blockfence(i),
                                                 d->blockstart(i),
                                                 matrix));
    }
}

ReplDiagSCMatrix::~ReplDiagSCMatrix()
{
}

RefSCDimension
ReplDiagSCMatrix::dim()
{
  return d;
}

double
ReplDiagSCMatrix::get_element(int i)
{
  return matrix[i];
}

void
ReplDiagSCMatrix::set_element(int i,double a)
{
  matrix[i] = a;
}

void
ReplDiagSCMatrix::accumulate_element(int i,double a)
{
  matrix[i] += a;
}

void
ReplDiagSCMatrix::assign(double val)
{
  int n = d->n();
  for (int i=0; i<n; i++) matrix[i] = val;
}

void
ReplDiagSCMatrix::accumulate(DiagSCMatrix*a)
{
  // make sure that the argument is of the correct type
  ReplDiagSCMatrix* la
    = ReplDiagSCMatrix::require_castdown(a,"ReplDiagSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      fprintf(stderr,"ReplDiagSCMatrix::"
              "accumulate(SCMatrix*a):\n");
      fprintf(stderr,"dimensions don't match\n");
      abort();
    }

  int nelem = n();
  for (int i=0; i<nelem; i++) matrix[i] += matrix[i];
}

double
ReplDiagSCMatrix::invert_this()
{
  double det = 1.0;
  int nelem = n();
  for (int i=0; i<nelem; i++) {
      det *= matrix[i];
      matrix[i] = 1.0/matrix[i];
    }
  return det;
}

double
ReplDiagSCMatrix::determ_this()
{
  double det = 1.0;
  int nelem = n();
  for (int i=0; i < nelem; i++) {
    det *= matrix[i];
  }
  return det;
}

double
ReplDiagSCMatrix::trace()
{
  double tr = 0;
  int nelem = n();
  for (int i=0; i < nelem; i++) {
    tr += matrix[i];
  }
  return tr;
}

void
ReplDiagSCMatrix::gen_invert_this()
{
  int nelem = n();
  for (int i=0; i < nelem; i++) {
    if (fabs(matrix[i]) > 1.0e-8)
      matrix[i] = 1.0/matrix[i];
    else
      matrix[i] = 0;
  }
}

void
ReplDiagSCMatrix::element_op(const RefSCElementOp& op)
{
  if (op->has_side_effects()) before_elemop();
  SCMatrixBlockListIter i;
  for (i = blocklist->begin(); i != blocklist->end(); i++) {
      op->process(i.block());
    }
  if (op->has_side_effects()) after_elemop();
}

void
ReplDiagSCMatrix::element_op(const RefSCElementOp2& op,
                              DiagSCMatrix* m)
{
  ReplDiagSCMatrix *lm
      = ReplDiagSCMatrix::require_castdown(m,"ReplDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim())) {
      fprintf(stderr,"ReplDiagSCMatrix: bad element_op\n");
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
ReplDiagSCMatrix::element_op(const RefSCElementOp3& op,
                              DiagSCMatrix* m,DiagSCMatrix* n)
{
  ReplDiagSCMatrix *lm
      = ReplDiagSCMatrix::require_castdown(m,"ReplDiagSCMatrix::element_op");
  ReplDiagSCMatrix *ln
      = ReplDiagSCMatrix::require_castdown(n,"ReplDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      fprintf(stderr,"ReplDiagSCMatrix: bad element_op\n");
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

// from Ed Seidl at the NIH (with a bit of hacking)
void
ReplDiagSCMatrix::print(const char *title, ostream& os, int prec)
{
  if (messagegrp()->me() != 0) return;

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

  if(n()==0) { os << " empty matrix\n"; return; }

  for (i=0; i<n(); i++) {
      os.width(5); os << i+1;
      os.width(lwidth); os << matrix[i];
      os << "\n";
    }
  os << "\n";

  os.flush();
}
