//
// localtest.cc
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

#include <util/keyval/keyval.h>
#include <math/scmat/local.h>

void matrixtest(RefSCMatrixKit kit, RefKeyVal keyval,
                RefSCDimension d1,RefSCDimension d2,RefSCDimension d3);

main()
{
  int i;

  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/matrixtest.in");
  RefSCMatrixKit kit = new LocalSCMatrixKit;

  matrixtest(kit,keyval,0,0,0);

  // SVD is tested here since its not implemented for other specializations

  RefSCDimension m(keyval->describedclassvalue("d1"));
  RefSCDimension n(keyval->describedclassvalue("d2"));
  RefSCDimension p = ((m.n() < n.n()) ? m:n);
  RefSCMatrix A(m,n,kit);
  RefSCMatrix U(m,m,kit);
  RefSCMatrix V(n,n,kit);
  RefDiagSCMatrix sigma(p,kit);

  A.randomize();
  A.svd(U,sigma,V);

  A.print("A");
  U.print("U");
  (U*U.t()).print("U*U.t()");
  (U.t()*U).print("U.t()*U");
  sigma.print("sigma");
  V.print("V");
  (V*V.t()).print("V*V.t()");
  (V.t()*V).print("V.t()*V");
  RefSCMatrix sigmamat(m,n,kit);
  sigmamat.assign(0.0);
  for (i=0; i<p.n(); i++) sigmamat(i,i) = sigma(i);
  (U*sigmamat*V.t()).print("U*sigmamat*V.t()");

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
