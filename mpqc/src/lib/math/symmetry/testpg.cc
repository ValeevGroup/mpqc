
#include <math/symmetry/pointgrp.h>

main(int argc, char *argv[])
{
  PointGroup pg(argv[1]);

  //pg.char_table().print();

  CharacterTable ct = pg.char_table();
  CharacterTable ct2;
  ct2=ct;

  ct.print();
}
