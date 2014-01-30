//
// opttest.cc
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
#include <util/state/stateio.h>
#include <math/optimize/function.h>
#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
#include <math/optimize/linkage.h>

using namespace std;
using namespace sc;

class Quadratic: public Function
{
  private:
    RefSCVector x0;
    RefSCVector g0;
    RefSymmSCMatrix h0;
    RefSymmSCMatrix hguess;
  public:
    Quadratic(StateIn&);
    Quadratic(const Ref<KeyVal>&);
    void save_data_state(StateOut&);
    void compute();
    void guess_hessian(RefSymmSCMatrix&);
};
static ClassDesc Quadratic_cd(
  typeid(Quadratic),"Quadratic",1,"public Function",
  0, create<Quadratic>, create<Quadratic>);
Quadratic::Quadratic(StateIn&s):
  SavableState(s),
  Function(s)
{
  x0 = matrixkit_->vector(dim_);
  x0.restore(s);
  g0 = matrixkit_->vector(dim_);
  g0.restore(s);
  h0 = matrixkit_->symmmatrix(dim_);
  h0.restore(s);
}
void
Quadratic::save_data_state(StateOut&s)
{
  Function::save_data_state(s);
  x0.save(s);
  g0.save(s);
  h0.save(s);
}
Quadratic::Quadratic(const Ref<KeyVal>&keyval):
  Function(keyval)
{
  set_dimension(new SCDimension(keyval->count("x0")));
  x0 = matrixkit_->vector(dim_);
  g0 = matrixkit_->vector(dim_);
  h0 = matrixkit_->symmmatrix(dim_);
  hguess = matrixkit_->symmmatrix(dim_);
  hguess.assign(0.0);
  Ref<SCElementOp> op(new SCElementShiftDiagonal(1.0));
  hguess.element_op(op);
  
  int dim = dimension()->n();
  for (int i=0; i<dim; i++) {
      x0(i) = keyval->doublevalue("x0",i);
      g0(i) = keyval->doublevalue("g0",i);
      for (int j=0; j<=i; j++) {
          h0(i,j) = keyval->doublevalue("h0",i,j);
          hguess(i,j) = keyval->doublevalue("hguess",i,j);
        }
    }
}
// this computes everything, whether or not it was requested
void
Quadratic::compute()
{
  cout << "Quadratic::compute(): entered\n";
  
  // compute the displacement from x0
  RefSCVector d = x_ - x0;

  // compute h * d
  RefSCVector h0d = h0 * d;
//   RefSCVector h0d(h0.dim());
//   int n=h0.dim().n();
//   for (int i=0; i<n; i++) {
//       double tmp = 0;
//       for (int j=0; j<n; j++) {
//           tmp += h0(i,j) * d(j);
//         }
//       h0d(i) = tmp;
//     }

  // compute the value
  value_.result_noupdate() =   d.scalar_product(g0)
                             + 0.5 * d.scalar_product(h0d);
  value_.computed() = 1;

  // compute the gradient
  gradient_.result_noupdate() = g0 + h0d;
  gradient_.computed() = 1;

  // compute the hessian
  hessian_.result_noupdate() = h0;
  hessian_.computed() = 1;
}
void
Quadratic::guess_hessian(RefSymmSCMatrix&gh)
{
  gh.assign(hguess);
}

int
main(int argc, char* argv[])
{
  Ref<KeyVal> kv = new ParsedKeyVal( SRCDIR "/opttest.in");
  Ref<KeyVal> pkv = new PrefixKeyVal(kv,"opt");

  for (int i=0; i<pkv->count(); i++) {
      Ref<Optimize> opt; opt << pkv->describedclassvalue(i);
      if (opt) {
          RefSCVector oldx = opt->function()->get_x();
          opt->optimize();
          // restore the orginal x, in case the function is used again
          opt->function()->set_x(oldx);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
