//
// matrixtest.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <iostream>
#include <math.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/elemop.h>
#include <util/misc/regtime.h>

using namespace std;
using namespace sc;

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
matrixtest(Ref<SCMatrixKit> kit, Ref<KeyVal> keyval,
           RefSCDimension d1,RefSCDimension d2,RefSCDimension d3,
           bool have_svd)
{
  if (d1 == 0) d1 << keyval->describedclassvalue("d1");
  if (d2 == 0) d2 << keyval->describedclassvalue("d2");
  if (d3 == 0) d3 << keyval->describedclassvalue("d3");

  d1.print();
  d2.print();
  d3.print();

  Timer tim("matrixtest");
  int i;
  int j;

  // seed the random number generator
  srand48(0);

  RefSCMatrix a(d1,d2,kit);
  RefSCMatrix a2(d1,d2,kit);
  RefSCMatrix a3(d1,d2,kit);
  RefSCMatrix b(d2,d3,kit);
  RefSCMatrix c(d1,d3,kit);

  cout << "a(" << a.nrow() << "," << a.ncol() << ")\n";

  a.assign(7.0);
  a2.assign(5.0);
  a3.assign(3.0);
  Ref<SCElementOp3> op3 = new Prod3;
  a.element_op(op3,a2,a3);
  a.print("a");
  a2.print("a2");
  a3.print("a3");

  /////////////////////////////////
  
  RefSymmSCMatrix sa(d3,kit);
  RefSymmSCMatrix sa2(d3,kit);
  RefSymmSCMatrix sa3(d3,kit);

  sa.assign(7.0);
  sa2.assign(5.0);
  sa3.assign(3.0);
  sa.element_op(op3,sa2,sa3);
  sa.print("sa");
  sa2.print("sa2");
  sa3.print("sa3");

  /////////////////////////////////
  
  RefDiagSCMatrix da(d3,kit);
  RefDiagSCMatrix da2(d3,kit);
  RefDiagSCMatrix da3(d3,kit);

  da.assign(7.0);
  da2.assign(5.0);
  da3.assign(3.0);
  da.element_op(op3,da2,da3);
  da.print("da");
  da2.print("da2");
  da3.print("da3");

  /////////////////////////////////
  
  RefSCVector vva(d3,kit);
  RefSCVector vva2(d3,kit);
  RefSCVector vva3(d3,kit);

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

  tim.enter("mxm");
  RefSCMatrix d = c * b.t();
  tim.exit("mxm");

  d.print("d");

  RefSCDimension d4; d4 << keyval->describedclassvalue("d4");
  int nd4 = d4->n();
  cout << "n4 = " << nd4 << endl;
  d4.print();
  RefSCMatrix aaa(d4,d4,kit);
  RefSCMatrix bbb(d4,d4,kit);
  aaa.assign(1.0);
  bbb.assign(2.0);
  tim.enter("mxm2");
  RefSCMatrix ccc = aaa*bbb;
  tim.exit("mxm2");

  d.print("d later");

  RefSymmSCMatrix e(d3,kit);

  e.assign(1.0);
  e.print("e");
  e.eigvals().print("e.eigvals()");
  e.eigvecs().print("e.eigvecs()");

  RefSymmSCMatrix f(d3,kit);
  for (i=0; i<d3.n(); i++) {
      for (j=0; j<=i; j++) {
          f(i,j) = i + sqrt((double)j);
        }
    }
  f.print("f");
  f.eigvals().print("f.eigvals()");
  f.i().print("f.i()");

  RefSymmSCMatrix h(d3,kit);
  for (i=0; i<d3.n(); i++) {
    for (j=0; j<=i; j++) {
      h(i,j) = f(i,j);
    }
  }
  h.print("h should be equal to f");

  RefSCMatrix g(d3,d3,kit);
  for (i=0; i<d3.n(); i++) {
      for (j=0; j<d3.n(); j++) {
          if (i>j) g(i,j) = i + sqrt((double)j);
          else g(i,j) = j + sqrt((double)i);
        }
    }
  g.print("g");
  g.i().print("g.i()");
  (g * g.i()).print("g * g.i()");

  if (have_svd) {
      g.gi().print("g.gi()");
      (g * g.gi()).print("g * g.gi()");
      (g.gi() * g).print("g.gi() * g");
    }

  RefSCVector v(d3,kit);
  for (i=0; i<d3.n(); i++) {
      v(i) = 1.0/(i+1);
    }
  v.print("Vector v");

  RefSCVector wa(d3,kit);
  RefSCMatrix ma(d1,d3,kit);
  randomize(ma);
  randomize(wa);
  RefSCVector va = ma * wa;
  ma.print("Matrix ma");
  va.print("Vector va");
  wa.print("Vector wa");

  if (have_svd) {
      ma.gi().print("ma.gi()");
      (ma * ma.gi()).print("ma * ma.gi()");
      (ma.gi() * ma).print("ma.gi() * ma");
    }

  RefSCVector wb(d1,kit);
  RefSCMatrix mb(d3,d1,kit);
  randomize(mb);
  randomize(wb);
  RefSCVector vb = mb * wb;
  ma.print("Matrix mb");
  va.print("Vector vb");
  wa.print("Vector wb");

  if (have_svd) {
      mb.gi().print("mb.gi()");
      (mb * mb.gi()).print("mb * mb.gi()");
      (mb.gi() * mb).print("mb.gi() * mb");
    }

  RefSymmSCMatrix bmbt(d3,kit);
  RefSCMatrix redundant_ortho(d2,d3,kit);
  RefSymmSCMatrix bmbt_fixed(d2,kit);
  RefSCMatrix bmbt_fix_red(d2,d3,kit);
  bmbt.assign(0.0);
  randomize(redundant_ortho);
  randomize(bmbt_fixed);
  randomize(bmbt_fix_red);
  bmbt.accumulate_transform(redundant_ortho.t(), bmbt_fixed);
  bmbt.print("bmbt");

  bmbt.assign(0.0);
  bmbt.accumulate_transform(redundant_ortho, bmbt_fixed,
                            SCMatrix::TransposeTransform);
  bmbt.print("bmbt (2)");

  bmbt.accumulate_symmetric_sum(redundant_ortho.t() * bmbt_fix_red);
  bmbt.print("bmbt (symmetric_sum)");

  RefSCMatrix bmbt_test;
  RefSCMatrix bmbt_fixed_test(d2,d2,kit);
  for (i=0; i<d2.n(); i++) {
      for (j=0; j<=i; j++) {
          bmbt_fixed_test(i,j) = bmbt_fixed(i,j);
          bmbt_fixed_test(j,i) = bmbt_fixed(i,j);
        }
    }
  RefSCMatrix tmp = redundant_ortho.t() * bmbt_fix_red;
  bmbt_test =  redundant_ortho.t() * bmbt_fixed_test * redundant_ortho
             + tmp + tmp.t();
  bmbt_test.print("bmbt_test");

  tim.exit("matrixtest");
  tim.print();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
