//
// david_hack_file.cpp
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

#include <chemistry/qc/libint2/libint2.h>
#include <util/madness/init.h>
#include <util/madness/world.h>
#include <chemistry/qc/basis/tiledbasisset.hpp>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/integralenginepool.hpp>
#include <chemistry/qc/basis/taskintegrals.hpp>
#include <iostream>

using namespace sc;
using namespace mpqc;
using namespace mpqc::TA;

using Matrix = TiledArray::Array<double, 3, TiledArray::Tensor<double> >;

void try_main(int argc, char** argv){
  sc::ExEnv::init(argc, argv);
  mpqc::MADNESSRuntime::initialize();

  Ref<World> world = new World();

  const char *input = "./benzene_trimer.kv";
  Ref<KeyVal> kv = new ParsedKeyVal(input);

  Ref<Molecule> mol; mol << kv->describedclassvalue("benzene_trimer");
  int occ = mol->total_Z()/2;

  Ref<TiledBasisSet> tbs; tbs << kv->describedclassvalue("basis");
  Ref<TiledBasisSet> dftbs; dftbs << kv->describedclassvalue("dfbasis");
  Ref<IntegralLibint2> ints_fac; ints_fac << kv->describedclassvalue("integrals");
  Integral::set_default_integral(ints_fac);
  ints_fac->set_basis(tbs,tbs,dftbs);

  // Get pools
  using twobpool = IntegralEnginePool<Ref<TwoBodyThreeCenterInt> >;
  auto E_pool = std::make_shared<twobpool>(
          ints_fac->electron_repulsion3()->clone());
  world->madworld()->gop.fence();

  // Compute ints
  const double start_time = madness::wall_time();
  Matrix S = Integrals(*world->madworld(), E_pool, tbs, dftbs);
  world->madworld()->gop.fence();
  const double stop_time = madness::wall_time();

  if(world->madworld()->rank() == 0){
    std::cout << "Time was = " << stop_time - start_time << std::endl;
  }

  mpqc::MADNESSRuntime::finalize();
}

int main(int argc, char** argv){
  try{
    try_main(argc, argv);
  }
  catch(std::exception &e){
    std::cout << " Caught an exception " << e.what() << std::endl;
  }
  return 0;
}

