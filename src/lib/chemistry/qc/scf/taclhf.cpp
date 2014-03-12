//
// taclhf.cpp
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

#include <chemistry/qc/scf/taclhf.hpp>
#include <chemistry/qc/basis/integralenginepool.hpp>
#include <chemistry/qc/basis/taskintegrals.hpp>

using namespace mpqc;
using namespace mpqc::TA;
using Matrix = mpqc::TA::CLHF::Matrix;


sc::ClassDesc mpqc::TA::CLHF::class_desc_(typeid(mpqc::TA::CLHF), "TA.CLHF",
                      1, "public TA.CLSCF",
                      0,
                      sc::create<mpqc::TA::CLHF>,
                      0);

mpqc::TA::CLHF::CLHF(const sc::Ref<sc::KeyVal>& kval) :
    CLSCF(kval)
{
}

#warning "Gmat uses all four centered ints"
Matrix mpqc::TA::CLHF::Gmat(){
    std::shared_ptr<IntegralEnginePool<sc::Ref<sc::TwoBodyInt> > >
        eri_pool(new IntegralEnginePool<sc::Ref<sc::TwoBodyInt> >(
                        integral()->electron_repulsion()));
    ::TiledArray::Array<double, 4> eri = Integrals(*(world())->madworld(),
                                                eri_pool, basis());
    Matrix Gmat = 2*density()("r,s") * eri("m,r,n,s") -
                    density()("r,s") * eri("m,r,s,n");
    return Gmat;
}

void mpqc::TA::CLHF::minimize_energy() {
    MPQC_ASSERT(false);
}
