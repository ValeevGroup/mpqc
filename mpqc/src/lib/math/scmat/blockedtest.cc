
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/blocked.h>

void matrixtest(RefSCMatrixKit, RefKeyVal,
                RefSCDimension d1,RefSCDimension d2,RefSCDimension d3);

main()
{
  int i;
  int nblks;
  int *blks1, *blks2, *blks3;

  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/matrixtest.in");
  RefSCMatrixKit subkit = new LocalSCMatrixKit;
  RefBlockedSCMatrixKit kit = new BlockedSCMatrixKit(subkit);

  nblks = keyval->intvalue("nblocks");
  if (!nblks)
    nblks=3;
  
  blks1 = new int[nblks];
  blks2 = new int[nblks];
  blks3 = new int[nblks];
  
  RefSCDimension sd1(keyval->describedclassvalue("d1"));
  RefSCDimension sd2(keyval->describedclassvalue("d2"));
  RefSCDimension sd3(keyval->describedclassvalue("d2"));
  
  RefSCDimension d1(new SCDimension(sd1.n()*nblks,nblks));
  RefSCDimension d2(new SCDimension(sd2.n()*nblks,nblks));
  RefSCDimension d3(new SCDimension(sd3.n()*nblks,nblks));

  for (i=0; i < nblks; i++) {
    d1->blocks()->set_subdim(i,sd1);
    d2->blocks()->set_subdim(i,sd2);
    d3->blocks()->set_subdim(i,sd3);
  }

  matrixtest(kit,keyval,d1,d2,d3);
  
  d1 = new SCDimension(sd1.n(),1);
  d2 = new SCDimension(sd2.n(),1);
  d3 = new SCDimension(sd3.n(),1);

  d1->blocks()->set_subdim(0,sd1);
  d2->blocks()->set_subdim(0,sd2);
  d3->blocks()->set_subdim(0,sd3);

  matrixtest(kit,keyval,d1,d2,d3);
  d1=0;

  return 0;
}
