//
// local.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// LocalSCMatrixKit member functions

#define CLASSNAME LocalSCMatrixKit
#define PARENTS public SCMatrixKit
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
LocalSCMatrixKit::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixKit::_castdown(cd);
  return do_castdowns(casts,cd);
}

LocalSCMatrixKit::LocalSCMatrixKit()
{
}

LocalSCMatrixKit::LocalSCMatrixKit(const RefKeyVal& keyval):
  SCMatrixKit(keyval)
{
}

LocalSCMatrixKit::~LocalSCMatrixKit()
{
}

SCMatrix*
LocalSCMatrixKit::matrix(const RefSCDimension&d1, const RefSCDimension&d2)
{
  return new LocalSCMatrix(d1,d2,this);
}

SymmSCMatrix*
LocalSCMatrixKit::symmmatrix(const RefSCDimension&d)
{
  return new LocalSymmSCMatrix(d,this);
}

DiagSCMatrix*
LocalSCMatrixKit::diagmatrix(const RefSCDimension&d)
{
  return new LocalDiagSCMatrix(d,this);
}

SCVector*
LocalSCMatrixKit::vector(const RefSCDimension&d)
{
  return new LocalSCVector(d,this);

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
