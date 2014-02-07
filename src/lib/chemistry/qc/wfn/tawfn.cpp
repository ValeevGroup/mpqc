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
#include "tawfn.hpp"

using namespace std;
using namespace mpqc;
namespace TA = TiledArray;
using TAMat = TiledArrayWavefunction::TAMat;

// WORK ON SERIALIZATION LATER
//static sc::ClassDesc TiledArrayWavefunction_cd(
//                typeid(TiledArrayWavefunction)),
//                "TiledArrayWavefunction", 9, "public MolecularEnergy",
//                0, 0, 0);

// mpqc::TiledArrayWavefunction::TiledArrayWavefunction(sc::StateIn& s):
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

mpqc::TiledArrayWavefunction::TiledArrayWavefunction(
                const sc::Ref<sc::KeyVal>& kval) :
    sc::MolecularEnergy(kval), overlap_(this), hcore_(this)
{
    overlap_.compute() = 0;
    hcore_.compute() = 0;
    overlap_.computed() = 0;
    hcore_.computed() = 0;
}

//void mpqc::TiledArrayWavefunction::save_data_state(sc::StateOut& s) {
//}

double mpqc::TiledArrayWavefunction::total_charge() const {
    return 2.0;
}

TAMat mpqc::TiledArrayWavefunction::ao_density() {
    return 2.0;
}

TAMat mpqc::TiledArrayWavefunction::ao_overlap() {
    return 2.0;
}
