//
// obint.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#include <chemistry/qc/cints/obintcints.h>

using namespace sc;

OneBodyIntCints::OneBodyIntCints(Integral* integral,
                           const Ref<GaussianBasisSet>&bs1,
                           const Ref<GaussianBasisSet>&bs2,
                           IntegralFunction ifunc):
  OneBodyInt(integral,bs1,bs2)
{
  bool need_overlap = true;
  bool need_coulomb = true;
  if (ifunc == &Int1eCints::nuclear) {
    need_overlap = false;
  }
  else if (ifunc == &Int1eCints::overlap || ifunc == &Int1eCints::kinetic ||
	   ifunc == &Int1eCints::edipole || ifunc == &Int1eCints::equadrupole) {
    need_coulomb = false;
  }

  int ntypes = 1;
  if (ifunc == &Int1eCints::edipole)
    ntypes = 3;
  if (ifunc == &Int1eCints::equadrupole)
    ntypes = 6;
    
  int1ecints_ = new Int1eCints(integral,bs1,bs2,0,need_overlap, need_coulomb, ntypes);
  intfunc_ = ifunc;
  buffer_ = int1ecints_->buffer();
}

OneBodyIntCints::~OneBodyIntCints()
{
}

void OneBodyIntCints::set_multipole_origin(const Ref<DipoleData>& origin)
{
  int1ecints_->set_multipole_origin(origin);
}

void OneBodyIntCints::set_EdotV_origin(const Ref<EfieldDotVectorData>& origin)
{
  int1ecints_->set_EdotV_origin(origin);
}

void OneBodyIntCints::set_Q_origin(const Ref<PointChargeData>& origin)
{
  int1ecints_->set_Q_origin(origin);
}

void
OneBodyIntCints::compute_shell(int i, int j)
{
  (int1ecints_.pointer()->*intfunc_)(i, j);
}

bool
OneBodyIntCints::cloneable() const
{
  return true;
}

Ref<OneBodyInt>
OneBodyIntCints::clone()
{
  Ref<OneBodyIntCints> obcints
      = new OneBodyIntCints(integral_, bs1_, bs2_, intfunc_);

  // make sure the full state of this object gets set up
  // in the clone
  obcints->set_multipole_origin(int1ecints_->multipole_origin());
  obcints->set_EdotV_origin(int1ecints_->EdotV_origin());
  obcints->set_Q_origin(int1ecints_->Q_origin());

  return obcints.pointer();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
