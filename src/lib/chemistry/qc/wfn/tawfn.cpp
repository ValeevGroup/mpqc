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
#include <tiledarray.h>
#include <mpqc/utility/mutex.hpp>
#include <util/misc/regtime.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/taskintegrals.hpp>
#include <chemistry/qc/basis/integralenginepool.hpp>

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

Wavefunction::Wavefunction(const sc::Ref<sc::KeyVal>& kval) :
    sc::MolecularEnergy(kval),
    overlap_(this),
    hcore_(this),
    rdm1_(this),
    rdm1_alpha_(this),
    rdm1_beta_(this)
{
  overlap_.compute() = 0;
  hcore_.compute() = 0;
  rdm1_.compute() = 0;
  rdm1_alpha_.compute() = 0;
  rdm1_beta_.compute() = 0;

  overlap_.computed() = 0;
  hcore_.computed() = 0;
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

Wavefunction::~Wavefunction(){}

//void mpqc::TA::Wavefunction::save_data_state(sc::StateOut& s) {
//}

const Wavefunction::TAMatrix&
Wavefunction::ao_overlap() {

  if (not overlap_.computed()) {

    sc::Timer tim("ao_overlap:");
    // Get integral pool
    integral_->set_basis(basis());
    auto overlap_pool =
        std::make_shared<IntegralEnginePool<sc::Ref<sc::OneBodyInt>>>
                                            (integral_->overlap());

    // Fill TAMatrix with integrals
    overlap_ = Integrals(*world_->madworld(), overlap_pool, tbs_);
    world_->madworld()->gop.fence();
    tim.exit("ao_overlap:");


    overlap_.computed() = 1; // Overlap_ is now computed
  }

  return overlap_.result_noupdate();
}

  return overlap_.result_noupdate();
}

const Wavefunction::TAMatrix&
Wavefunction::ao_hcore() {

  if (not hcore_.computed()) {

    sc::Timer tim("ao_hcore:");

    // Get integral pool
    integral_->set_basis(basis());
    auto hcore_pool =
      std::make_shared<IntegralEnginePool<sc::Ref<sc::OneBodyInt>>>
                                             (integral_->hcore());

    // Fill TAMatrix with integrals
    hcore_ = Integrals(*world_->madworld(), hcore_pool, tbs_);
    world_->madworld()->gop.fence();

    tim.exit("ao_hcore:");

    hcore_.computed() = 1; // hcore_ is now computed
  }

  return hcore_.result_noupdate();
}

Wavefunction::TAMatrixExpr
Wavefunction::rdm1_expr(std::string input) {
  return rdm1()(input);
}

Wavefunction::TAMatrixExpr
Wavefunction::ao_overlap_expr(std::string input){
  return ao_overlap()(input);
}

Wavefunction::TAMatrixExpr
Wavefunction::ao_hcore_expr(std::string input){
  return ao_hcore()(input);
}

double Wavefunction::magnetic_moment() const {
  //if (magnetic_moment_ > extent(osorange_)) // magnetic moment greater than the number of states means it has not been computed yet.
  //  magnetic_moment_ = trace(rdm1(Alpha), overlap()) - trace(rdm1(Beta), overlap());
  //return magnetic_moment_;
  throw sc::FeatureNotImplemented("mpqc::v3::Wavefunction::magnetic_moment() not yet implemented", __FILE__, __LINE__);
}

bool Wavefunction::nonzero_efield_supported() const {
  // support efields in C1 symmetry only
  if (molecule()->point_group()->char_table().order() == 1)
    return true;
  return false;
}

void Wavefunction::obsolete() {
  magnetic_moment_ = tbs_->nbasis() + 1;
  MolecularEnergy::obsolete();
}

void Wavefunction::print(std::ostream& os) const {
  sc::MolecularEnergy::print(os);
  os << sc::indent << "Electronic basis:" << std::endl;
  os << sc::incindent;
  basis()->print_brief(os);
  os << sc::decindent;
  os << sc::indent << "Integral factory = " << integral()->class_name() << std::endl;
  //os << sc::indent << "magnetic moment = " << magnetic_moment() << std::endl;
}

