/*--------------------*/
/* Standard libraries */
/*--------------------*/

#define WANT_STREAM
#define WANT_MATH
#include "include.h"
#include "newmatap.h"
#include "nlp.h"
#include "cblas.h"

main ()
{
  int i;
  int NMAX = 10;
  int n = 3;
  double xnorm, ynorm, xdoty;
  ColumnVector x(n), y(n);
  SymmetricMatrix Hk(n);
  LowerTriangularMatrix L(n);


  ProgramException::SetAction(1);
  DataException::SetAction(1);

  x = 1.0;
  y = 2.0;

  printf("Enter dimension size:");
  scanf("%d",&n);

#ifdef tstvec /* Sat Jan  2 22:46:01 1993 */
/**/  for (i=1; i<= n; ++i)
/**/    printf("%d  %e  %e\n", i,x(i),y(i));
/**/  xnorm = Norm2(x);
/**/  ynorm = Norm2(y);
/**/  printf("xnorm = %e\n",xnorm);
/**/  printf("ynorm = %e\n",ynorm);
/**/  xnorm = x.Norm2();
/**/  ynorm = y.Norm2();
/**/  printf("xnorm = %e\n",xnorm);
/**/  printf("ynorm = %e\n",ynorm);
/**/  xdoty = x.Dot(y);
/**/  printf("xdoty = %e\n",xdoty);
/**/  xdoty = Dot(x,y);
/**/  printf("xdoty = %e\n",xdoty);
/**/
#endif /* tstvec Sat Jan  2 22:46:01 1993 */


  Hk(1,1) = 1.0;
  Hk(2,1) = 1.0;
  Hk(3,1) = 2.0;
  Hk(2,2) = 1.0+1.e-10;
  Hk(3,2) = 3.0;
  Hk(3,3) = 1.0;
  Print(Hk);

  L = MCholesky(Hk);
  Print(Hk);
  Print(L);

}
