
#include <math/scmat/local.h>

void matrixtest(RefSCDimension d1,RefSCDimension d2,RefSCDimension d3);

main()
{
  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/matrixtest.in");
  RefSCMatrixKit kit = new LocalSCMatrixKit;

  RefSCDimension d1(kit->dimension(keyval->intvalue("n1")));
  RefSCDimension d2(kit->dimension(keyval->intvalue("n2")));
  RefSCDimension d3(kit->dimension(keyval->intvalue("n3")));

  matrixtest(kit,keyval,d1,d2,d3);

  return 0;
}
