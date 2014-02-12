//
// tiledbasisset_test.cpp
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

#include "tiledbasisset.hpp"
#include <chemistry/qc/basis/integral.h>
#include <mpqc/integrals/integralenginepool.hpp>
#include <mpqc/interfaces/tiledarray/array_ints.hpp>

using namespace mpqc;
int main(int argc, char** argv){
    madness::World &world = madness::initialize(argc, argv);
    sc::Ref<sc::Molecule>  mol = new sc::Molecule;
    /*
    mol->add_atom( 6,     0,     0,    0);
    mol->add_atom( 9,    -1,    -1,    0);
    mol->add_atom( 1,   0.6,  -0.1,  0.9);
    mol->add_atom(17, -0.75,   1.5,    0);
    mol->add_atom(35,   1.1, -0.18, -1.5);
    */
    mol->add_atom( 1, 0, 1, -1);
    mol->add_atom( 1, 0, 1, 1);

    sc::Ref<sc::AssignedKeyVal> akv = new  sc::AssignedKeyVal;
    akv->assign("name", "3-21G");
    akv->assign("molecule", mol.pointer());
    sc::Ref<sc::GaussianBasisSet> bs = new sc::GaussianBasisSet(akv);
    sc::Ref<TiledBasisSet> tbasis =
                    new TiledBasisSet(bs, 2);

    sc::Ref<sc::Integral> Int_fac = sc::Integral::initial_integral(argc, argv);
    if(Int_fac)
        sc::Integral::set_default_integral(Int_fac);
    Int_fac = sc::Integral::get_default_integral()->clone();
    Int_fac->set_basis(tbasis);

    IntegralEnginePool<sc::Ref<sc::OneBodyInt> > overlap_pool(Int_fac->overlap());
    TiledArray::Array<double, 2> S = Integrals(world, overlap_pool,
                                               tbasis);

    world.gop.fence();
    std::cout << "S = \n" << S << std::endl;


    return 0;
}


