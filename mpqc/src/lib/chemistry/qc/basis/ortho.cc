
#include <math.h>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/obint.h>
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

  ints->set_basis(t);
  RefSCElementOp overlap = new OneBodyIntOp(ints->overlap());
  
  RefSymmSCMatrix ov(dim, t->matrixkit());
  ov.element_op(overlap);
  overlap=0;

  RefSCMatrix trans(dim,dim,t->matrixkit());
  RefDiagSCMatrix eigval(dim,t->matrixkit());

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
