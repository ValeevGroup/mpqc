//
// blocked.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////
// BlockedSCMatrixKit member functions

DescribedClass_REF_def(BlockedSCMatrixKit);

#define CLASSNAME BlockedSCMatrixKit
#define PARENTS public SCMatrixKit
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
BlockedSCMatrixKit::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixKit::_castdown(cd);
  return do_castdowns(casts,cd);
}

BlockedSCMatrixKit::BlockedSCMatrixKit(const RefSCMatrixKit&subkit):
  subkit_(subkit)
{
}

BlockedSCMatrixKit::BlockedSCMatrixKit(const RefKeyVal& keyval):
  SCMatrixKit(keyval)
{
  subkit_ = keyval->describedclassvalue("subkit");
}

BlockedSCMatrixKit::~BlockedSCMatrixKit()
{
}

SCMatrix*
BlockedSCMatrixKit::matrix(const RefSCDimension&d1, const RefSCDimension&d2)
{
  int i;
  for (i=0; i<d1->blocks()->nblock(); i++) {
      if (d1->blocks()->subdim(i).null()) {
          ExEnv::err() << indent
               << "BlockedSCMatrixKit: given a dim without subdim info"
               << endl;
          abort();
        }
    }
  for (i=0; i<d2->blocks()->nblock(); i++) {
      if (d2->blocks()->subdim(i).null()) {
          ExEnv::err() << indent
               << "BlockedSCMatrixKit: given a dim without subdim info"
               << endl;
          abort();
        }
    }
  return new BlockedSCMatrix(d1,d2,this);
}

SymmSCMatrix*
BlockedSCMatrixKit::symmmatrix(const RefSCDimension&d)
{
  for (int i=0; i<d->blocks()->nblock(); i++) {
      if (d->blocks()->subdim(i).null()) {
          ExEnv::err() << indent
               << "BlockedSCMatrixKit: given a dim without subdim info"
               << endl;
          abort();
        }
    }
  return new BlockedSymmSCMatrix(d,this);
}

DiagSCMatrix*
BlockedSCMatrixKit::diagmatrix(const RefSCDimension&d)
{
  for (int i=0; i<d->blocks()->nblock(); i++) {
      if (d->blocks()->subdim(i).null()) {
          ExEnv::err() << indent
               << "BlockedSCMatrixKit: given a dim without subdim info"
               << endl;
          abort();
        }
    }
  return new BlockedDiagSCMatrix(d,this);
}

SCVector*
BlockedSCMatrixKit::vector(const RefSCDimension&d)
{
  for (int i=0; i<d->blocks()->nblock(); i++) {
      if (d->blocks()->subdim(i).null()) {
          ExEnv::err() << indent
               << "BlockedSCMatrixKit: given a dim without subdim info"
               << endl;
          abort();
        }
    }
  return new BlockedSCVector(d,this);
}

/////////////////////////////////////////////////////////////////////////////

#define CLASSNAME BlockedSCElementOp
#define PARENTS public SCElementOp
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
BlockedSCElementOp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}

BlockedSCElementOp::BlockedSCElementOp()
{
  current_block_=0;
}

void
BlockedSCElementOp::working_on(int b)
{
  current_block_ = b;
}

int
BlockedSCElementOp::current_block() const
{
  return current_block_;
}

#define CLASSNAME BlockedSCElementOp2
#define PARENTS public SCElementOp2
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
BlockedSCElementOp2::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp2::_castdown(cd);
  return do_castdowns(casts,cd);
}

BlockedSCElementOp2::BlockedSCElementOp2()
{
  current_block_=0;
}

void
BlockedSCElementOp2::working_on(int b)
{
  current_block_ = b;
}

int
BlockedSCElementOp2::current_block() const
{
  return current_block_;
}

#define CLASSNAME BlockedSCElementOp3
#define PARENTS public SCElementOp3
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
BlockedSCElementOp3::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp3::_castdown(cd);
  return do_castdowns(casts,cd);
}

BlockedSCElementOp3::BlockedSCElementOp3()
{
  current_block_=0;
}

void
BlockedSCElementOp3::working_on(int b)
{
  current_block_ = b;
}

int
BlockedSCElementOp3::current_block() const
{
  return current_block_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
