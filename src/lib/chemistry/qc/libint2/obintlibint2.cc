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

#include <chemistry/qc/libint2/obintlibint2.h>

using namespace sc;

OneBodyIntLibint2::OneBodyIntLibint2(Integral* integral,
                           const Ref<GaussianBasisSet>&bs1,
                           const Ref<GaussianBasisSet>&bs2,
                           IntegralFunction ifunc):
  OneBodyInt(integral,bs1,bs2)
{
  bool need_overlap = true;
  bool need_coulomb = true;
  if (ifunc == &Int1eLibint2::nuclear ||
      ifunc == &Int1eLibint2::efield ||
      ifunc == &Int1eLibint2::efield_grad) {
    need_overlap = false;
  }
  else if (ifunc == &Int1eLibint2::overlap || ifunc == &Int1eLibint2::kinetic ||
	   ifunc == &Int1eLibint2::edipole || ifunc == &Int1eLibint2::equadrupole ||
	   ifunc == &Int1eLibint2::p4) {
    need_coulomb = false;
  }

  int ntypes = 1;
  if (ifunc == &Int1eLibint2::edipole || ifunc == &Int1eLibint2::efield)
    ntypes = 3;
  if (ifunc == &Int1eLibint2::equadrupole || ifunc == &Int1eLibint2::efield_grad)
    ntypes = 6;
    
  int1elibint2_ = new Int1eLibint2(integral,bs1,bs2,0,need_overlap, need_coulomb, ntypes);
  intfunc_ = ifunc;
  buffer_ = int1elibint2_->buffer();
}

OneBodyIntLibint2::~OneBodyIntLibint2()
{
}

void OneBodyIntLibint2::set_params(const Ref<IntParams>& p)
{
  int1elibint2_->set_params(p);
}

void OneBodyIntLibint2::set_EdotV_origin(const Ref<EfieldDotVectorData>& origin)
{
  int1elibint2_->set_EdotV_origin(origin);
}

void OneBodyIntLibint2::set_Q_origin(const Ref<PointChargeData>& origin)
{
  int1elibint2_->set_Q_origin(origin);
}

void
OneBodyIntLibint2::compute_shell(int i, int j)
{
  (int1elibint2_.pointer()->*intfunc_)(i, j);
}

bool
OneBodyIntLibint2::cloneable() const
{
  return true;
}

Ref<OneBodyInt>
OneBodyIntLibint2::clone()
{
  Ref<OneBodyIntLibint2> oblibint2
      = new OneBodyIntLibint2(integral_, bs1_, bs2_, intfunc_);

  // make sure the full state of this object gets set up
  // in the clone
  oblibint2->set_params(int1elibint2_->params());
  oblibint2->set_EdotV_origin(int1elibint2_->EdotV_origin());
  oblibint2->set_Q_origin(int1elibint2_->Q_origin());

  return oblibint2.pointer();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
