
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>

void matrixtest(RefSCMatrixKit, RefKeyVal,
                RefSCDimension d1,RefSCDimension d2,RefSCDimension d3);

main()
{
  int i;

  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/matrixtest.in");
  RefSCMatrixKit kit = new LocalSCMatrixKit;

  RefSCDimension d1(kit->dimension(keyval->intvalue("n1")));
  RefSCDimension d2(kit->dimension(keyval->intvalue("n2")));
  RefSCDimension d3(kit->dimension(keyval->intvalue("n3")));

  matrixtest(kit,keyval,d1,d2,d3);

  // SVD is tested here since its not implemented for other specializations

  RefSCDimension m(kit->dimension(keyval->intvalue("n1")));
  RefSCDimension n(kit->dimension(keyval->intvalue("n2")));
  RefSCDimension p = ((m.n() < n.n()) ? m:n);
  RefSCMatrix A(m,n);
  RefSCMatrix U(m,m);
  RefSCMatrix V(n,n);
  RefDiagSCMatrix sigma(p);

  A.randomize();
  A.svd(U,sigma,V);

  A.print("A");
  U.print("U");
  (U*U.t()).print("U*U.t()");
  (U.t()*U).print("U.t()*U");
  sigma.print("sigma");
  V.print("V");
  (V*V.t()).print("V*V.t()");
  (V.t()*V).print("V.t()*V");
  RefSCMatrix sigmamat(m,n);
  sigmamat.assign(0.0);
  for (i=0; i<p.n(); i++) sigmamat(i,i) = sigma(i);
  (U*sigmamat*V.t()).print("U*sigmamat*V.t()");

  return 0;
}
