
#include "pointgrp.h"

main(int argc, char *argv[])
{
  PointGroup pg(argv[1]);

  pg.char_table().print();
}
