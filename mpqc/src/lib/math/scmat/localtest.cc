
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>

void matrixtest(RefSCMatrixKit kit, RefKeyVal keyval,
                RefSCDimension d1,RefSCDimension d2,RefSCDimension d3);

main()
{
  int i;

  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/matrixtest.in");
  RefSCMatrixKit kit = new LocalSCMatrixKit;

  matrixtest(kit,keyval,0,0,0);

  // SVD is tested here since its not implemented for other specializations

  RefSCDimension m(keyval->describedclassvalue("d1"));
  RefSCDimension n(keyval->describedclassvalue("d2"));
  RefSCDimension p = ((m.n() < n.n()) ? m:n);
  RefSCMatrix A(m,n,kit);
  RefSCMatrix U(m,m,kit);
  RefSCMatrix V(n,n,kit);
  RefDiagSCMatrix sigma(p,kit);

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
  RefSCMatrix sigmamat(m,n,kit);
  sigmamat.assign(0.0);
  for (i=0; i<p.n(); i++) sigmamat(i,i) = sigma(i);
  (U*sigmamat*V.t()).print("U*sigmamat*V.t()");

  return 0;
}
