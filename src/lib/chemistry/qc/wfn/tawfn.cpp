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
#include <mpqc/interfaces/tiledarray/array_ints.hpp>
#include <mpqc/integrals/integralenginepool.hpp>

using namespace std;
using namespace mpqc;
using namespace mpqc::TA;

sc::ClassDesc Wavefunction::class_desc_(
                typeid(mpqc::TA::Wavefunction),
                "TA.Wavefunction", 1, "public MolecularEnergy",
                0, 0, 0);

// mpqc::TA::Wavefunction::Wavefunction(sc::StateIn& s):
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

mpqc::TA::Wavefunction::Wavefunction(const sc::Ref<sc::KeyVal>& kval) :
    sc::MolecularEnergy(kval), overlap_(this), rdm1_(this), rdm1_alpha_(this), rdm1_beta_(this)
{
  overlap_.compute() = 0;
  rdm1_.compute() = 0;
  rdm1_alpha_.compute() = 0;
  rdm1_beta_.compute() = 0;

  overlap_.computed() = 0;
  rdm1_.computed() = 0;
  rdm1_alpha_.computed() = 0;
  rdm1_beta_.computed() = 0;

  world_ << kval->describedclassvalue("world");
  if (world_.null())
    world_ = new mpqc::World;

  tbs_ << kval->describedclassvalue("basis");
  if (tbs_.null()) { // did the user give regular GaussianBasisSet? If so, try conversion
    sc::Ref<sc::GaussianBasisSet> bs;
    bs << kval->describedclassvalue("basis");
    if (bs.null())
      throw sc::InputError("missing \"basis\" keyword", __FILE__, __LINE__, "basis", "", this->class_desc());
    else {
      tbs_ = new TiledBasisSet(bs);
    }
  }

  integral_ << kval->describedclassvalue("integrals");
  if (integral_.null())
    integral_ = sc::Integral::get_default_integral()->clone();
  integral_->set_basis(tbs_);

  magnetic_moment_ = tbs_->nbasis() + 1;
}

mpqc::TA::Wavefunction::~Wavefunction()
{
}

//void mpqc::TA::Wavefunction::save_data_state(sc::StateOut& s) {
//}

double mpqc::TA::Wavefunction::total_charge() const {
  return molecule()->total_charge() - nelectron();
}

const mpqc::TA::Wavefunction::Matrix&
mpqc::TA::Wavefunction::rdm1() {
  if (not rdm1_.computed()) {
    if (rdm1_.result_noupdate().is_initialized() == false)
      rdm1_("a,b") = rdm1(sc::Alpha)("a,b") + rdm1(sc::Beta)("a,b");
    rdm1_.computed() = 1;
  }
  return rdm1_.result_noupdate();
}

const mpqc::TA::Wavefunction::Matrix&
mpqc::TA::Wavefunction::overlap() {

  if (not overlap_.computed()) {

    mpqc::IntegralEnginePool<sc::Ref<sc::OneBodyInt> > overlap_pool(integral_->overlap());
    overlap_ = mpqc::Integrals(*world_->madworld(), overlap_pool, tbs_);

    world_->madworld()->gop.fence();
    overlap_.computed() = 1;
  }

  return overlap_.result_noupdate();
}

double
mpqc::TA::Wavefunction::magnetic_moment() const {
  //if (magnetic_moment_ > extent(osorange_)) // magnetic moment greater than the number of states means it has not been computed yet.
  //  magnetic_moment_ = trace(rdm1(Alpha), overlap()) - trace(rdm1(Beta), overlap());
  //return magnetic_moment_;
  throw sc::FeatureNotImplemented("mpqc::v3::Wavefunction::magnetic_moment() not yet implemented", __FILE__, __LINE__);
}

bool mpqc::TA::Wavefunction::nonzero_efield_supported() const {
  // support efields in C1 symmetry only
  if (molecule()->point_group()->char_table().order() == 1)
    return true;
  return false;
}

void mpqc::TA::Wavefunction::obsolete() {
  magnetic_moment_ = tbs_->nbasis() + 1;
  MolecularEnergy::obsolete();
}

void mpqc::TA::Wavefunction::print(std::ostream& os) const {
  sc::MolecularEnergy::print(os);
  os << sc::indent << "Electronic basis:" << std::endl;
  os << sc::incindent;
  basis()->print_brief(os);
  os << sc::decindent;
  os << sc::indent << "Integral factory = " << integral()->class_name() << std::endl;
  //os << sc::indent << "magnetic moment = " << magnetic_moment() << std::endl;
}
