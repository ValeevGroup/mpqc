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

#include <util/madness/init.h>
#include <util/misc/regtime.h>
#include <util/madness/world.h>
#include <iostream>
#include <chemistry/qc/basis/tiledbasisset.hpp>
#include <chemistry/qc/basis/integral.h>
#include <mpqc/integrals/integralenginepool.hpp>
#include <mpqc/interfaces/tiledarray/array_ints.hpp>
#define BOOST_TEST_MODULE test_tabasis
#include <boost/test/included/unit_test.hpp>



using namespace boost::unit_test;
using namespace sc;
using namespace mpqc;

struct MADConfig {
    MADConfig() {
      mpqc::MADNESSRuntime::initialize();
    }
    ~MADConfig() {
      mpqc::MADNESSRuntime::finalize();
    }
};

BOOST_GLOBAL_FIXTURE( MADConfig );

BOOST_AUTO_TEST_CASE( construct_tiledbasisset_from_gaussianbasisset){
    Ref<World> world = new mpqc::World;

    Ref<Molecule>  mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 0.0);
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(10, 2.8, 0.0, 1.8);
    mol->add_atom(18, 0.0, 2.8, 1.8);
    mol->add_atom(36, 0.0, -2.8, 1.8);
    mol->add_atom(18, -2.8, 0.0, 1.8);

    Ref<AssignedKeyVal> akv = new  AssignedKeyVal;
    akv->assign("name", "STO-3G");
    akv->assign("molecule", mol.pointer());
    Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

    Ref<Integral> Int_fac = Integral::get_default_integral()->clone();

    for(size_t i = 1; i <= 6; ++i){
        Ref<TA::TiledBasisSet> tbasis = new TA::TiledBasisSet(bs, i);

        Int_fac->set_basis(tbasis);

        IntegralEnginePool<Ref<OneBodyInt> > overlap_pool(Int_fac->overlap());
        TiledArray::Array<double, 2> S = Integrals(*world->madworld(),
                                                   overlap_pool, tbasis);

        world->madworld()->gop.fence();
        std::cout << "S(" << i << ") = \n" << S << std::endl;
    }
}


