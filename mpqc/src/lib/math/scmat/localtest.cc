
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>

void matrixtest(RefSCMatrixKit, RefKeyVal,
                RefSCDimension d1,RefSCDimension d2,RefSCDimension d3);

main()
{
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
  RefSCMatrix U(m,p);
  RefSCMatrix V(n,p);
  RefDiagSCMatrix sigma(p);

  A.randomize();
  A.svd(U,sigma,V);

  A.print("A");
  U.print("U");
  sigma.print("sigma");
  V.print("V");
  (U*sigma*V.t()).print("U*sigma*V.t()");

  return 0;
}
