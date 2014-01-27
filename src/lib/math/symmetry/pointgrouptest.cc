
#include <math/symmetry/pointgrp.h>
#include <math/symmetry/corrtab.h>

using namespace sc;

int
main(int argc, char *argv[])
{
  Ref<PointGroup> pg = new PointGroup(argv[1]);
  Ref<PointGroup> pg2;
  if (argc > 2) pg2 = new PointGroup(argv[2]);

  //pg.char_table().print();

  CharacterTable ct = pg->char_table();
  CharacterTable ct2;
  ct2=ct;

  ct.print();

  if (pg2) {
      pg2->char_table().print();
      if (argc <= 3) {
          CorrelationTable corrtab(pg,pg2);
          corrtab.print();
        }
    }

  // test given axis rearrangements
  for (int i=3; i<argc; i++) {
      pg2->symm_frame().zero();
      for (int j=0; j<3; j++) {
          if (argv[i][j] == 'x') pg2->symm_frame()(j,0) = 1.0;
          else if (argv[i][j] == 'y') pg2->symm_frame()(j,1) = 1.0;
          else if (argv[i][j] == 'z') pg2->symm_frame()(j,2) = 1.0;
        }
      CorrelationTable corrtab(pg,pg2);
      corrtab.print();
    }
}
