//
// repl.cc
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
#include <math/scmat/repl.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// ReplSCMatrixKit member functions

#define CLASSNAME ReplSCMatrixKit
#define PARENTS public SCMatrixKit
#define HAVE_KEYVAL_CTOR
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
ReplSCMatrixKit::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixKit::_castdown(cd);
  return do_castdowns(casts,cd);
}

ReplSCMatrixKit::ReplSCMatrixKit()
{
}

ReplSCMatrixKit::ReplSCMatrixKit(const RefKeyVal& keyval):
  SCMatrixKit(keyval)
{
}

ReplSCMatrixKit::~ReplSCMatrixKit()
{
}

SCMatrix*
ReplSCMatrixKit::matrix(const RefSCDimension&d1, const RefSCDimension&d2)
{
  return new ReplSCMatrix(d1,d2,this);
}

SymmSCMatrix*
ReplSCMatrixKit::symmmatrix(const RefSCDimension&d)
{
  return new ReplSymmSCMatrix(d,this);
}

DiagSCMatrix*
ReplSCMatrixKit::diagmatrix(const RefSCDimension&d)
{
  return new ReplDiagSCMatrix(d,this);
}

SCVector*
ReplSCMatrixKit::vector(const RefSCDimension&d)
{
  return new ReplSCVector(d,this);

}

/////////////////////////////////////////////////////////////////////////////
// ReplSCMatrixKit member functions

ReplSCMatrixListSubblockIter::ReplSCMatrixListSubblockIter(
    Access access,
    const RefSCMatrixBlockList &list,
    const RefMessageGrp &grp,
    double *data,
    int ndata
    ):
  SCMatrixListSubblockIter(access, list),
  grp_(grp),
  data_(data),
  ndata_(ndata)
{
  if (access == Write) {
      for (int i=0; i<ndata; i++) data[i] = 0.0;
    }
}

ReplSCMatrixListSubblockIter::~ReplSCMatrixListSubblockIter()
{
  if (access() == Write || access() == Accum) {
      grp_->sum(data_,ndata_);
    }
}

///////////////////////////////////////////////////////////////////////
// The static SCMatrixKit members.

static RefSCMatrixKit defaultmatrixkit;

SCMatrixKit*
SCMatrixKit::default_matrixkit()
{
  if (defaultmatrixkit.null()) defaultmatrixkit = new ReplSCMatrixKit;
  return defaultmatrixkit.pointer();
}

void
SCMatrixKit::set_default_matrixkit(const RefSCMatrixKit &k)
{
  defaultmatrixkit = k;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
