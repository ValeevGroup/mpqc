
#include <math.h>
#include <math/newmat7/newmatap.h>
#include "integralv2.h"

void ortho(const GaussianBasisSet*t,Matrix&or,Matrix*orinv)
{
  GaussianOverlapIntv2 overlap(t);
  int n = t->nbasis();
  SymmetricMatrix ov(n);
  overlap.compute(ov);

  Matrix trans(n,n);
  DiagonalMatrix eigval(n);
  EigenValues(ov,eigval,trans);

  or.ReDimension(n,n);
  or = trans * eigval * trans.t();

  for (int i=0; i<n; i++) eigval.element(i) = sqrt(eigval.element(i));
  if (orinv) {
      orinv->ReDimension(n,n);
      *orinv = trans * eigval * trans.t();
    }
  for (i=0; i<n; i++) eigval.element(i) = 1.0/eigval.element(i);

  or = trans * eigval * trans.t();
}

void GaussianBasisSet::ortho(Matrix&or) const
{
  ::ortho(this,or,0);
}

void GaussianBasisSet::ortho(Matrix&or,Matrix&orinv) const
{
  ::ortho(this,or,&orinv);
}
