
#include <math/scmat/local.h>

void matrixtest(RefSCDimension d1,RefSCDimension d2,RefSCDimension d3);

main()
{
  RefSCDimension d1(new LocalSCDimension(2));
  RefSCDimension d2(new LocalSCDimension(3));
  RefSCDimension d3(new LocalSCDimension(4));
  matrixtest(d1,d2,d3);
  return 0;
}
