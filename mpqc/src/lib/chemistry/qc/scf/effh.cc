//
// effh.cc --- implementation of effective fock matrix builders
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

#include <chemistry/qc/scf/effh.h>

using namespace sc;

///////////////////////////////////////////////////////////////////////////
// AccumEffectiveH

AccumEffectiveH::AccumEffectiveH(SCF*s) :
  scf_(s)
{
}

AccumEffectiveH::~AccumEffectiveH()
{
}

int
AccumEffectiveH::index(int hindex, int shelli, int shellj)
{
  if (shellj > shelli) {
    int tmp = shelli;
    shelli = shellj;
    shellj = tmp;
  }

  return hindex * 9 + ((shelli+1)*shelli)/2 + shellj;
}

int
AccumEffectiveH::shell(double occ)
{
  if (occ==2.0) return 0;
  if (occ < 2.0 && occ > 0.0) return 1;
  return 2;
}

void
AccumEffectiveH::process(SCMatrixBlockIter&i,SCMatrixBlockIter&j)
{
  int ir=current_block();

  for (i.reset(),j.reset(); i; ++i,++j) {
    double occi = scf_->occupation(ir, i.i());
    double occj = scf_->occupation(ir, i.j());
    
    int ri = shell(occi);
    int rj = shell(occj);

    i.set(i.get() * coef_[index(0, ri, rj)]
          + j.get() * coef_[index(1, ri, rj)]);
  }
}

///////////////////////////////////////////////////////////////////////////
// GSGeneralEffH

void
GSGeneralEffH::init()
{
  coef(0,0,0) =  1.0;
  coef(0,1,0) =  2.0; coef(0,1,1) = 1.0;
  coef(0,2,0) =  1.0; coef(0,2,1) = 0.0; coef(0,2,2) = 1.0;

  coef(1,0,0) =  0.0;
  coef(1,1,0) = -1.0; coef(1,1,1) = 0.0;
  coef(1,2,0) =  0.0; coef(1,2,1) = 1.0; coef(1,2,2) = 0.0;
}

GSGeneralEffH::GSGeneralEffH(SCF *s) :
  AccumEffectiveH(s)
{
  init();
}

GSGeneralEffH::~GSGeneralEffH()
{
}

///////////////////////////////////////////////////////////////////////////
// GSHighSpinEffH

void
GSHighSpinEffH::init()
{
  coef(0,0,0) =  2.0;
  coef(0,1,0) =  2.0; coef(0,1,1) = 2.0;
  coef(0,2,0) =  1.0; coef(0,2,1) = 0.0; coef(0,2,2) = 2.0;

  coef(1,0,0) = -1.0;
  coef(1,1,0) = -1.0; coef(1,1,1) = -1.0;
  coef(1,2,0) =  0.0; coef(1,2,1) =  1.0; coef(1,2,2) = -1.0;
}

GSHighSpinEffH::GSHighSpinEffH(SCF* s) :
  AccumEffectiveH(s)
{
  init();
}

GSHighSpinEffH::~GSHighSpinEffH()
{
}

///////////////////////////////////////////////////////////////////////////
// TestEffH

void
TestEffH::init()
{
  coef(0,0,0) =  0.0;
  coef(0,1,0) =  2.0; coef(0,1,1) = 0.0;
  coef(0,2,0) =  1.0; coef(0,2,1) = 0.0; coef(0,2,2) = 0.0;

  coef(1,0,0) =  1.0;
  coef(1,1,0) = -1.0; coef(1,1,1) = 1.0;
  coef(1,2,0) =  0.0; coef(1,2,1) = 1.0; coef(1,2,2) = 1.0;
}

TestEffH::TestEffH(SCF* s) :
  AccumEffectiveH(s)
{
  init();
}

TestEffH::~TestEffH()
{
}

///////////////////////////////////////////////////////////////////////////
// PsiEffH

void
PsiEffH::init()
{
  coef(0,0,0) =  1.0;
  coef(0,1,0) =  2.0; coef(0,1,1) = 0.0;
  coef(0,2,0) =  1.0; coef(0,2,1) = 0.0; coef(0,2,2) = 0.0;

  coef(1,0,0) =  0.0;
  coef(1,1,0) = -1.0; coef(1,1,1) = 1.0;
  coef(1,2,0) =  0.0; coef(1,2,1) = 1.0; coef(1,2,2) = 1.0;
}

PsiEffH::PsiEffH(SCF*s) :
  AccumEffectiveH(s)
{
  init();
}

PsiEffH::~PsiEffH()
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
