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

#include <mpqc/tiledarrayscf/tascf.hpp>

using namespace mpqc;
using TAMat = TiledArrayScf::TAMat;

static sc::ClassDesc TiledArrayScf_cd(typeid(TiledArrayScf), "TiledArrayScf",
                      7, "public TiledArrayWavefunction", 0,
                      sc::create<TiledArrayScf>,0);

mpqc::TiledArrayScf::TiledArrayScf(const sc::Ref<sc::KeyVal>& kval) :
    TiledArrayWavefunction(kval), tbints_()
{
    if(kval->exists("maxiter"))
        maxiter_= kval->intvalue("maxiter");
    if(kval->exists("miniter"))
        miniter_= kval->intvalue("miniter");
}

mpqc::TiledArrayScf::~TiledArrayScf(){};

void mpqc::TiledArrayScf::compute() {}

int mpqc::TiledArrayScf::nelectron() {
    return 2.0;
}

TAMat mpqc::TiledArrayScf::ao_fock() {
    return 2.0;
}

TAMat mpqc::TiledArrayScf::ao_density() {
    return 2.0;
}

TAMat mpqc::TiledArrayScf::ao_overlap() {
    return 2.0;
}

double mpqc::TiledArrayScf::scf_energy() {
    return 2.0;
}
