
#include <util/misc/regtime.h>

int
main()
{
  int i;
  RefRegionTimer tim = new RegionTimer("top", 1, 1);
  tim->enter("main");

  tim->enter("x");
  double x = 0.0;
  for (i=0; i<10000000; i++) {
      x += 0.0001;
    }
  tim->enter("subx");
  sleep(2);
  tim->exit("subx");
  cout << " x = " << x << endl;
  tim->exit("x");
  tim->enter("a");
  double a = 0.0;
  for (i=0; i<10000000; i++) {
      a += 0.0001;
    }
  tim->enter("subx");
  sleep(1);
  tim->exit("subx");
  cout << " a = " << a << endl;
  tim->exit("a");
  tim->enter("y");
  double y = 0.0;
  for (i=0; i<10000000; i++) {
      y += 0.0001;
    }
  cout << " y = " << y << endl;
  tim->change("z", "y");
  double z = 0.0;
  for (i=0; i<10000000; i++) {
      z += 0.0001;
    }
  cout << " z = " << z << endl;
  tim->exit();
  tim->exit("main");

  tim->print();

  return 0;
}
