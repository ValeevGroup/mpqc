
#include <math.h>
#include <math/scmat/matrix.h>
#include "integralv2.h"

void ortho(const GaussianBasisSet*t,RefSCMatrix&or,RefSCMatrix&orinv)
{
  int n = t->nbasis();
  RefSCDimension dim = or.rowdim();
  if (dim.n() != n) {
      fprintf(stderr,"chemistry/qc/integral/ortho:ortho: dim.n() != n\n");
      abort();
    }

  RefSCSymmElementOp overlap = new GaussianOverlapIntv2(t);
  
  RefSymmSCMatrix ov(dim);
  ov.element_op(overlap);

  RefSCMatrix trans(dim,dim);
  RefDiagSCMatrix eigval(dim);

  ov.diagonalize(eigval,trans);

  RefSCVectorElementOp squareroot = new SCElementSquareRoot;
  eigval.element_op(squareroot);

  if (orinv) {
      orinv.assign(trans * eigval * trans.t());
    }

  RefSCVectorElementOp invert = new SCElementInvert;
  eigval.element_op(invert);

  or.assign(trans * eigval * trans.t());
}

void GaussianBasisSet::ortho(RefSCMatrix&or) const
{
  RefSCMatrix orinv;
  ::ortho(this,or,orinv);
}

void GaussianBasisSet::ortho(RefSCMatrix&or,RefSCMatrix&orinv) const
{
  ::ortho(this,or,orinv);
}
