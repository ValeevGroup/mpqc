
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>

void matrixtest(RefSCMatrixKit, RefKeyVal,
                RefSCDimension d1,RefSCDimension d2,RefSCDimension d3);

main()
{
  int i;
  int nblks;
  int *blks1, *blks2, *blks3;

  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/matrixtest.in");
  RefBlockedSCMatrixKit kit = new BlockedSCMatrixKit;

  nblks = keyval->intvalue("nblocks");
  if (!nblks)
    nblks=3;
  
  blks1 = new int[nblks];
  blks2 = new int[nblks];
  blks3 = new int[nblks];
  
  for (i=0; i < nblks; i++) {
    blks1[i] = keyval->intvalue("n1");
    blks2[i] = keyval->intvalue("n2");
    blks3[i] = keyval->intvalue("n3");
  }
  
  RefSCDimension d1(kit->dimension(nblks,blks1));
  RefSCDimension d2(kit->dimension(nblks,blks2));
  RefSCDimension d3(kit->dimension(nblks,blks3));

  matrixtest(kit,keyval,d1,d2,d3);

  d1 = kit->dimension(blks1[0]);
  d2 = kit->dimension(blks2[0]);
  d3 = kit->dimension(blks3[0]);

  matrixtest(kit,keyval,d1,d2,d3);
  d1=0;

  return 0;
}
