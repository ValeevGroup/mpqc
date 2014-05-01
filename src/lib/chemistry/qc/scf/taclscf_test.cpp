//
// taclscf_test.cpp
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
#include <iostream>
#include <chemistry/qc/scf/taclscf.hpp>
#define BOOST_TEST_MODULE test_taclscf
#include <boost/test/included/unit_test.hpp>

#include <chemistry/qc/libint2/linkage.h>

using namespace boost::unit_test;
using namespace sc;
using namespace mpqc;
using namespace mpqc::TA;

struct Mock_CLSCF : public CLSCF {
  Mock_CLSCF(const sc::Ref<sc::KeyVal> &kval) : CLSCF(kval) {}
  ~Mock_CLSCF(){}

};

struct MADConfig {
    MADConfig() {
      mpqc::MADNESSRuntime::initialize();
    }
    ~MADConfig() {
      mpqc::MADNESSRuntime::finalize();
    }
};

BOOST_GLOBAL_FIXTURE( MADConfig );

BOOST_AUTO_TEST_CASE( construct_clscf_programmatically ){

    // Make a molecule H2
    Ref<Molecule> mol = new Molecule;
    mol->add_atom(1,0,1,-1);
    mol->add_atom(1,0,1,1);

    // Make keyval
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("name", "3-21G");
    akv->assign("molecule", mol.pointer());
    Ref<GaussianBasisSet> basis =
                    new GaussianBasisSet(Ref<KeyVal>(akv));
    akv->assign("basis", basis.pointer());

    //Construct object
    Ref<Mock_CLSCF> tscf = new Mock_CLSCF(akv);
    tscf->print();
    std::cout << "Overlap \n" << tscf->overlap() << std::endl;
}











