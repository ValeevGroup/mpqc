//
// opttest.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#include <iostream.h>
#include <math/optimize/function.h>
#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>

class Quadratic: public Function
{
#   define CLASSNAME Quadratic
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>    
  private:
    RefSCVector x0;
    RefSCVector g0;
    RefSymmSCMatrix h0;
    RefSymmSCMatrix hguess;
  public:
    Quadratic(StateIn&);
    Quadratic(const RefKeyVal&);
    void save_data_state(StateOut&);
    void compute();
    void guess_hessian(RefSymmSCMatrix&);
};
#define CLASSNAME Quadratic
#define PARENTS public Function
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
Quadratic::Quadratic(StateIn&s):
  SavableState(s,class_desc_),
  Function(s)
{
  x0.restore_state(s);
  g0.restore_state(s);
  h0.restore_state(s);
}
void
Quadratic::save_data_state(StateOut&s)
{
  Function::save_data_state(s);
  x0.save_state(s);
  g0.save_state(s);
  h0.save_state(s);
}
void *
Quadratic::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { Function::_castdown(cd) };
  return do_castdowns(casts,cd);
}
Quadratic::Quadratic(const RefKeyVal&keyval):
  Function(new LocalSCDimension(keyval->count("x0")))
{
  x0 = dimension()->create_vector();
  g0 = dimension()->create_vector();
  h0 = dimension()->create_symmmatrix();
  hguess = dimension()->create_symmmatrix();
  hguess.assign(0.0);
  RefSCElementOp op(new SCElementShiftDiagonal(1.0));
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
  RefSCVector d = _x - x0;

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
  _value.result_noupdate() =   d.scalar_product(g0)
                             + 0.5 * d.scalar_product(h0d);
  _value.computed() = 1;

  // compute the gradient
  _gradient.result_noupdate() = g0 + h0d;
  _gradient.computed() = 1;

  // compute the hessian
  _hessian.result_noupdate() = h0;
  _hessian.computed() = 1;
}
void
Quadratic::guess_hessian(RefSymmSCMatrix&gh)
{
  gh.assign(hguess);
}

main()
{
  RefKeyVal kv = new ParsedKeyVal( SRCDIR "/opttest.in");
  RefKeyVal pkv = new PrefixKeyVal("opt",*kv);

  for (int i=0; i<pkv->count(); i++) {
      RefOptimize opt(pkv->describedclassvalue(i));
      if (opt.nonnull()) {
          RefSCVector oldx = opt->function()->get_x().copy();
          opt->optimize();
          // restore the orginal x, in case the function is used again
          opt->function()->set_x(oldx);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
