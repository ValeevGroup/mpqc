
#include <math.h>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>

static void
ortho(const RefIntegral& ints, const RefGaussianBasisSet&t,
      const RefSCMatrix&or, const RefSCMatrix&orinv)
{
  int n = t->nbasis();
  RefSCDimension dim = or.rowdim();
  if (dim.null()) {
      dim = orinv.coldim();
    }
  if (dim.n() != n) {
      fprintf(stderr,"chemistry/qc/integral/ortho:ortho: dim.n() != n\n");
      abort();
    }

  RefSCElementOp overlap = ints->overlap_op(t);
  
  RefSymmSCMatrix ov(dim);
  ov.element_op(overlap);

  RefSCMatrix trans(dim,dim);
  RefDiagSCMatrix eigval(dim);

  ov.diagonalize(eigval,trans);

  RefSCElementOp squareroot = new SCElementSquareRoot;
  eigval.element_op(squareroot);

  if (orinv.nonnull()) {
      orinv.assign(trans * eigval * trans.t());
    }

  if (or.nonnull()) {
      RefSCElementOp invert = new SCElementInvert;
      eigval.element_op(invert);
      or.assign(trans * eigval * trans.t());
    }
}

void
GaussianBasisSet::ortho(const RefIntegral& ints, const RefSCMatrix& or)
{
  RefSCMatrix orinv;
  ::ortho(ints,this,or,orinv);
}

void
GaussianBasisSet::ortho(const RefIntegral& ints, const RefSCMatrix& or,
                        const RefSCMatrix& orinv)
{
  ::ortho(ints,this,or,orinv);
}
