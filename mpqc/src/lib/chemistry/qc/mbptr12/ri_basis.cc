//
// ri_basis.cc
//
// Copyright (C) 2004 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <sstream>

#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/mbptr12.h>

using namespace sc;

void
MBPT2_R12::construct_ri_basis_(bool safe)
{
  if (aux_basis_->equiv(basis())) {
    ri_basis_ = basis();
  }
  else {
    switch(abs_method_) {
      case LinearR12::ABS_KS:
	construct_ri_basis_ks_(safe);
	break;
      case LinearR12::ABS_KSPlus:
	construct_ri_basis_ksplus_(safe);
	break;
      case LinearR12::ABS_EV:
	construct_ri_basis_ev_(safe);
	break;
      case LinearR12::ABS_EVPlus:
	construct_ri_basis_evplus_(safe);
	break;
      default:
	throw std::runtime_error("MBPT2_R12::construct_ri_basis_ -- invalid ABS method");
    }
  }
}

void
MBPT2_R12::construct_ri_basis_ks_(bool safe)
{
  if (!abs_spans_obs_() && safe)
    throw std::runtime_error("MBPT2_R12::construct_ri_basis_ks_ -- auxiliary basis is not safe to use with the given orbital basis");
  ri_basis_ = aux_basis_;
}

void
MBPT2_R12::construct_ri_basis_ksplus_(bool safe)
{
  GaussianBasisSet& bs = *(basis().pointer());
  ri_basis_ = bs + aux_basis_;
}

void
MBPT2_R12::construct_ri_basis_ev_(bool safe)
{
  ri_basis_ = basis();
  if (!abs_spans_obs_() && safe)
    throw std::runtime_error("MBPT2_R12::construct_ri_basis_ev_ -- auxiliary basis is not safe to use with the given orbital basis");
  construct_ortho_comp_();
}

void
MBPT2_R12::construct_ri_basis_evplus_(bool safe)
{
  GaussianBasisSet& bs = *(basis().pointer());
  ri_basis_ = bs + aux_basis_;
  construct_ortho_comp_();
}


bool
MBPT2_R12::abs_spans_obs_()
{
  return true;
}

void
MBPT2_R12::construct_ortho_comp_()
{
  throw std::runtime_error("MBPT2_R12::construct_ortho_comp_ -- EV ABS method has not been implemented yet");
}
