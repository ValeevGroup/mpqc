
#include <iostream.h>
#include <math.h>
#include <stdio.h>
#include <util/misc/libmisc.h>
#include <util/group/picl.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/elemop.h>

class Prod3: public SCElementOp3 {
  private:
  public:
    void process(SCMatrixBlockIter& i1,
                 SCMatrixBlockIter& i2,
                 SCMatrixBlockIter& i3) {
        for (i1.reset(),i2.reset(),i3.reset();
             i1;
             ++i1,++i2,++i3) {
            i1.set(i2.get()*i3.get());
          }
      }
    int has_side_effects() { return 1; }
};

void
randomize(RefSCMatrix&m)
{
  for (int i=0; i<m.rowdim().n(); i++) {
      for (int j=0; j<m.coldim().n(); j++) {
          m.set_element(i,j,drand48());
        }
    }
}

void
randomize(RefSymmSCMatrix&m)
{
  for (int i=0; i<m.dim().n(); i++) {
      for (int j=0; j<=i; j++) {
          m.set_element(i,j,drand48());
        }
    }
}

void
randomize(RefSCVector&m)
{
  for (int i=0; i<m.dim().n(); i++) {
      m.set_element(i,drand48());
    }
}

// test abstract matrices
void
matrixtest(
    RefSCMatrixKit kit, RefKeyVal keyval,
    RefSCDimension d1,RefSCDimension d2,RefSCDimension d3)
{
  // The tim_enter routines require PICL
  int numproc, me, host;
  open0(&numproc, &me, &host);

  tim_enter("matrixtest");
  int i;
  int j;

  // seed the random number generator
  srand48(0);

  RefSCMatrix a(d1,d2);
  RefSCMatrix a2(d1,d2);
  RefSCMatrix a3(d1,d2);
  RefSCMatrix b(d2,d3);
  RefSCMatrix c(d1,d3);

  cout << "a(" << a.nrow() << "," << a.ncol() << ")\n";

  a.assign(7.0);
  a2.assign(5.0);
  a3.assign(3.0);
  RefSCElementOp3 op3 = new Prod3;
  a.element_op(op3,a2,a3);
  a.print("a");
  a2.print("a2");
  a3.print("a3");

  /////////////////////////////////
  
  RefSymmSCMatrix sa(d3);
  RefSymmSCMatrix sa2(d3);
  RefSymmSCMatrix sa3(d3);

  sa.assign(7.0);
  sa2.assign(5.0);
  sa3.assign(3.0);
  sa.element_op(op3,sa2,sa3);
  sa.print("sa");
  sa2.print("sa2");
  sa3.print("sa3");

  /////////////////////////////////
  
  RefDiagSCMatrix da(d3);
  RefDiagSCMatrix da2(d3);
  RefDiagSCMatrix da3(d3);

  da.assign(7.0);
  da2.assign(5.0);
  da3.assign(3.0);
  da.element_op(op3,da2,da3);
  da.print("da");
  da2.print("da2");
  da3.print("da3");

  /////////////////////////////////
  
  RefSCVector vva(d3);
  RefSCVector vva2(d3);
  RefSCVector vva3(d3);

  vva.assign(7.0);
  vva2.assign(5.0);
  vva3.assign(3.0);
  vva.element_op(op3,vva2,vva3);
  vva.print("vva");
  vva2.print("vva2");
  vva3.print("vva3");

  ////////////////////////////////

  a.assign(0.0);
  b.assign(1.0);
  c.assign(2.0);

  a.print("a");
  b.print("b");
  c.print("c");

  tim_enter("mxm");
  RefSCMatrix d = c * b.t();
  tim_exit("mxm");

  d.print("d");

  int nd4 = keyval->intvalue("n4");
  if (!nd4) nd4 = 1;
  cout << "n4 = " << nd4 << endl;
  RefSCDimension d4 = kit->dimension(nd4);
  RefSCMatrix aaa(d4,d4);
  RefSCMatrix bbb(d4,d4);
  aaa.assign(1.0);
  bbb.assign(2.0);
  tim_enter("mxm2");
  RefSCMatrix ccc = aaa*bbb;
  tim_exit("mxm2");

  d.print("d later");

  RefSymmSCMatrix e(d3);

  e.assign(1.0);
  e.print("e");
  e.eigvals().print("e.eigvals()");
  e.eigvecs().print("e.eigvecs()");

  RefSymmSCMatrix f(d3);
  for (i=0; i<d3.n(); i++) {
      for (j=0; j<=i; j++) {
          f(i,j) = i + sqrt((double)j);
        }
    }
  f.print("f");
  f.eigvals().print("f.eigvals()");
  f.i().print("f.i()");

  RefSymmSCMatrix h(d3);
  for (i=0; i<d3.n(); i++) {
    for (j=0; j<=i; j++) {
      h(i,j) = f(i,j);
    }
  }
  h.print("h should be equal to f");

  RefSCMatrix g(d3,d3);
  for (i=0; i<d3.n(); i++) {
      for (j=0; j<d3.n(); j++) {
          if (i>j) g(i,j) = i + sqrt((double)j);
          else g(i,j) = j + sqrt((double)i);
        }
    }
  g.print("g");
  g.i().print("g.i()");
  (g * g.i()).print("g * g.i()");

  RefSCVector v(d3);
  for (i=0; i<d3.n(); i++) {
      v(i) = 1.0/(i+1);
    }
  v.print("Vector v");

  RefSCVector wa(d3);
  RefSCMatrix ma(d1,d3);
  randomize(ma);
  randomize(wa);
  RefSCVector va = ma * wa;
  ma.print("Matrix ma");
  va.print("Vector va");
  wa.print("Vector wa");

  RefSCVector wb(d1);
  RefSCMatrix mb(d3,d1);
  randomize(mb);
  randomize(wb);
  RefSCVector vb = mb * wb;
  ma.print("Matrix mb");
  va.print("Vector vb");
  wa.print("Vector wb");

  RefSymmSCMatrix bmbt(d3);
  RefSCMatrix redundant_ortho(d2,d3);
  RefSymmSCMatrix bmbt_fixed(d2);
  RefSCMatrix bmbt_fix_red(d2,d3);
  bmbt.assign(0.0);
  randomize(redundant_ortho);
  randomize(bmbt_fixed);
  randomize(bmbt_fix_red);
  bmbt.accumulate_transform(redundant_ortho.t(), bmbt_fixed);
  bmbt.accumulate_symmetric_sum(redundant_ortho.t() * bmbt_fix_red);
  cout << "bmbt:\n";
  bmbt.print();

  RefSCMatrix bmbt_test;
  RefSCMatrix bmbt_fixed_test(d2,d2);
  for (i=0; i<d2.n(); i++) {
      for (j=0; j<=i; j++) {
          bmbt_fixed_test(i,j) = bmbt_fixed(i,j);
          bmbt_fixed_test(j,i) = bmbt_fixed(i,j);
        }
    }
  RefSCMatrix tmp = redundant_ortho.t() * bmbt_fix_red;
  bmbt_test =  redundant_ortho.t() * bmbt_fixed_test * redundant_ortho
             + tmp + tmp.t();
  cout << "bmbt_test\n";
  bmbt_test.print();

  tim_exit("matrixtest");
  tim_print(0);

  close0(0);
}
