//
// tawfn.cpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <chemistry/qc/wfn/tawfn.hpp>
#include <util/misc/scexception.h>

using namespace std;
using namespace mpqc;
using namespace mpqc::v3;
namespace TA = TiledArray;

sc::ClassDesc Wavefunction::class_desc_(
                typeid(mpqc::v3::Wavefunction),
                "v3.Wavefunction", 1, "public MolecularEnergy",
                0, 0, 0);

// mpqc::v3::Wavefunction::Wavefunction(sc::StateIn& s):
//     sc::SaveableState(s),
//     sc::MolecularEnergy(s),
//     overlap_(this),
//     hcore_(this)
// {
//     overlap_.compute() = 0;
//     overlap_.computed() = 0;
//     hcore_.compute() = 0;
//     hcore_.computed() = 0;
// }

mpqc::v3::Wavefunction::Wavefunction(const sc::Ref<sc::KeyVal>& kval) :
    sc::MolecularEnergy(kval)
{
}

mpqc::v3::Wavefunction::~Wavefunction()
{
}

//void mpqc::v3::Wavefunction::save_data_state(sc::StateOut& s) {
//}

double mpqc::v3::Wavefunction::total_charge() const {
  return molecule()->total_charge() - nelectron();
}

const TA::Array<double,2>&
mpqc::v3::Wavefunction::ao_density() {
  throw sc::FeatureNotImplemented("mpqc::v3::Wavefunction::ao_density() not yet ready", __FILE__, __LINE__);
}

const TA::Array<double,2>&
mpqc::v3::Wavefunction::ao_overlap() {
  throw sc::FeatureNotImplemented("mpqc::v3::Wavefunction::ao_overlap() not yet ready", __FILE__, __LINE__);
}

double
mpqc::v3::Wavefunction::magnetic_moment() const {

//  if (magnetic_moment_ > extent(osorange_)) // magnetic moment greater than the number of states means it has not been computed yet.
//    magnetic_moment_ = trace(alpha_density(), overlap()) -
//                       trace(beta_density(), overlap());
//  return magnetic_moment_;
  throw sc::FeatureNotImplemented("mpqc::v3::Wavefunction::magnetic_moment() not yet implemented", __FILE__, __LINE__);
}

bool mpqc::v3::Wavefunction::nonzero_efield_supported() const {
  // support efields in C1 symmetry only
  if (molecule()->point_group()->char_table().order() == 1)
    return true;
  return false;
}
