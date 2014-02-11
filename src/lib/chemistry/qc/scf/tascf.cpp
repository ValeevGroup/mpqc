//
// tascf.cpp
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

#include <chemistry/qc/scf/tascf.hpp>

using namespace mpqc;
using namespace mpqc::v3;

sc::ClassDesc mpqc::v3::SCF::class_desc_(typeid(mpqc::v3::SCF), "v3.SCF",
                      1, "public v3.Wavefunction",
                      0,
                      sc::create<mpqc::v3::SCF>,
                      0);

mpqc::v3::SCF::SCF(const sc::Ref<sc::KeyVal>& kval) :
    Wavefunction(kval), tbints_()
{
    if(kval->exists("maxiter"))
        maxiter_= kval->intvalue("maxiter", sc::KeyValValueint(1000));
    if(kval->exists("miniter"))
        miniter_= kval->intvalue("miniter", sc::KeyValValueint(0));
}

mpqc::v3::SCF::~SCF(){
}

void mpqc::v3::SCF::compute() {
  MPQC_ASSERT(false);
}

const mpqc::v3::SCF::Matrix&
mpqc::v3::SCF::ao_fock() {
  MPQC_ASSERT(false);
}

const mpqc::v3::SCF::Matrix&
mpqc::v3::SCF::ao_density() {
  MPQC_ASSERT(false);
}

const mpqc::v3::SCF::Matrix&
mpqc::v3::SCF::ao_overlap() {
  MPQC_ASSERT(false);
}

double
mpqc::v3::SCF::scf_energy() {
  MPQC_ASSERT(false);
}

size_t
mpqc::v3::SCF::nelectron() const {

}
