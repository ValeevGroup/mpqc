//
// taclscf.cpp
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

#include <chemistry/qc/scf/taclscf.hpp>
#include <tiledarray.h>
#include <chemistry/qc/basis/integralenginepool.hpp>
#include <chemistry/qc/basis/taskintegrals.hpp>
#include <chemistry/qc/lcao/soad.h>
#include <mpqc/interfaces/tiledarray/symmscmat.hpp>
#include <math/elemental/eigensolver.hpp>

using namespace mpqc;
using namespace mpqc::TA;
using TAMatrix = mpqc::TA::CLSCF::TAMatrix;

sc::ClassDesc mpqc::TA::CLSCF::class_desc_(typeid(mpqc::TA::CLSCF), "TA.CLSCF",
                                           1, "public TA.SCF", 0, 0, 0);


CLSCF::CLSCF(const sc::Ref<sc::KeyVal>& kval) :
        SCF(kval) {
  if (nelectron() % 2 != 0) {
    throw sc::InputError("Number of electrons is not divisable by two",
      __FILE__, __LINE__, "", "", this->class_desc());
  }
  set_occupation(nelectron() / 2);
}

CLSCF::~CLSCF() { }

#warning "compute is not yet defined"
void CLSCF::compute() {
  MPQC_ASSERT(false);
}

const TAMatrix& CLSCF::rdm1() {
  // Check if computed, and if not then check that Evecs are at desired accuracy
  if(!rdm1_.computed()){

    if(!Coeff_.is_initialized()){
      Coeff_ = eigensolver_occ_Coeff(ao_fock(),
                                     ao_overlap(),
                                     occupation());
    }

    rdm1_.result_noupdate()("mu,nu") = Coeff_("mu,i") * Coeff_("nu,i");
    rdm1_.computed() = 1;

  }

  return rdm1_.result_noupdate();
}

const TAMatrix& CLSCF::rdm1(sc::SpinCase1){
  // To return a reference we have to store the data, just use alpha for
  // Closed shell systems
  if(!rdm1_alpha_.computed()){
    rdm1_alpha_.result_noupdate()("i,j") = 
                                0.5 * rdm1_.result_noupdate()("i,j");
    rdm1_alpha_.computed() = 1;
  }
  return rdm1_alpha_.result_noupdate();
}

// If we have not generated an intial guess for the density yet then do so now.
TAMatrix& CLSCF::ao_density() {

  // Check to see if data has been initialized if
  if (!Wavefunction::ao_density().is_initialized()) {
    // For constructing a SOAD object
    sc::Ref<sc::AssignedKeyVal> akv = new sc::AssignedKeyVal;
    akv->assign("molecule", molecule().pointer());
    akv->assign("basis", basis().pointer());
    akv->assign("integrals", integral().pointer());

    // Construct a density guess based on a SOAD object
    using Soad = sc::SuperpositionOfAtomicDensities;
    sc::Ref<Soad> guess = new Soad(sc::Ref<sc::KeyVal>(akv));

    // Copy the mpqc sc matrix into our tiledarray Matrix.
    Wavefunction::ao_density() = mpqc::SymmScMat_To_TiledArray(*world()->madworld(),
                                      guess->guess_density(basis(),integral()),
                                      basis()->trange1());

    world()->madworld()->gop.fence();
  }

  return Wavefunction::ao_density();
}

double CLSCF::scf_energy(){
  // E = \sum_{ij} \left( D_{ij} * (F_{ij} + H_{ij}) \right)
  return rdm1()("i,j").dot(ao_hcore()("i,j") + ao_fock()("i,j"))
                        + molecule()->nuclear_repulsion_energy();
}

double CLSCF::iter_energy() {
  // E = \sum_{ij} \left( D_{ij} * (F_{ij} + H_{ij}) \right)
  return ao_density()("i,j").dot(ao_hcore()("i,j") +
                                        scf_ao_fock_()("i,j"));
}

