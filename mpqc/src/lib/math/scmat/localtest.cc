
#include <iostream.h>
#include <math.h>
#include <math/scmat/local.h>

main()
{
  int i;
  int j;
  
  RefSCDimension d1(new LocalSCDimension(2));
  RefSCDimension d2(new LocalSCDimension(3));
  RefSCDimension d3(new LocalSCDimension(4));

  RefSCVector v(d3);

  v.print("Vector v");

  RefSCMatrix a(d1,d2);
  RefSCMatrix b(d2,d3);
  RefSCMatrix c(d1,d3);

  cout << "a(" << a.nrow() << "," << a.ncol() << ")\n";

  a.assign(0.0);
  b.assign(1.0);
  c.assign(2.0);

  a.print("a");
  b.print("b");
  c.print("c");

  RefSCMatrix d = c * b.t();

  d.print("d");

  RefSymmSCMatrix e(d3);

  e.assign(1.0);
  e.eigvals().print("e.eigvals()");
  e.eigvecs().print("e.eigvecs()");

  RefSymmSCMatrix f(d3);
  for (i=0; i<d3->n(); i++) {
      for (j=0; j<=i; j++) {
          f(i,j) = i + sqrt((double)j);
        }
    }
  f.print("f");
  f.eigvals().print("f.eigvals()");
  f.i().print("f.i()");
  //(f * f.i()).print("f * f.i()");

  RefSCMatrix g(d3,d3);
  for (i=0; i<d3->n(); i++) {
      for (j=0; j<d3->n(); j++) {
          if (i>j) g(i,j) = i + sqrt((double)j);
          else g(i,j) = j + sqrt((double)i);
        }
    }
  g.print("g");
  g.i().print("g.i()");
  (g * g.i()).print("g * g.i()");

}
