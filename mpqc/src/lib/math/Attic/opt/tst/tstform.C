/*--------------------*/
/* Standard libraries */
/*--------------------*/

#define WANT_STREAM
#include "nlp.h"

main ()
{
  int iter;
  int NMAX = 10;
  int n = 5;
  double x,y,z;

  iter = 1;
  x = 1.0;
  y = 2.0;
  z = 3.0;
  printf("using printf\n");
  printf("x = %e, y = %e\n",x,y);
  cout << "using cout\n";
  cout << form("x = %e, y = %e\n",x,y);
  cout.flush();
}
