//
// ortho.cc
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

#include <math.h>
#include <util/misc/formio.h>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/integral.h>

static void
ortho(const RefIntegral& ints, const RefGaussianBasisSet&t,
      const RefSCMatrix&or, const RefSCMatrix&orinv)
{
  int n = t->nbasis();
  RefSCDimension dim = or.rowdim();
  if (dim.null()) {
      dim = orinv.coldim();
    }
  if (dim.n() != n) {
      cerr << node0 << indent
           << "chemistry/qc/integral/ortho:ortho: dim.n() != n\n";
      abort();
    }

  ints->set_basis(t);
  RefSCElementOp overlap = new OneBodyIntOp(ints->overlap());
  
  RefSymmSCMatrix ov(dim, t->matrixkit());
  ov.element_op(overlap);
  overlap=0;

  RefSCMatrix trans(dim,dim,t->matrixkit());
  RefDiagSCMatrix eigval(dim,t->matrixkit());

  ov.diagonalize(eigval,trans);

  RefSCElementOp squareroot = new SCElementSquareRoot;
  eigval.element_op(squareroot);

  if (orinv.nonnull()) {
      orinv.assign(trans * eigval * trans.t());
    }

  if (or.nonnull()) {
      RefSCElementOp invert = new SCElementInvert;
      eigval.element_op(invert);
      or.assign(trans * eigval * trans.t());
    }
}

void
GaussianBasisSet::ortho(const RefIntegral& ints, const RefSCMatrix& or)
{
  RefSCMatrix orinv;
  ::ortho(ints,this,or,orinv);
}

void
GaussianBasisSet::ortho(const RefIntegral& ints, const RefSCMatrix& or,
                        const RefSCMatrix& orinv)
{
  ::ortho(ints,this,or,orinv);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
