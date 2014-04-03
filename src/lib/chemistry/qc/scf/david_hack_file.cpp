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

#include <chemistry/qc/scf/cldfgengine.hpp>
#include <chemistry/qc/libint2/libint2.h>
#include <util/madness/init.h>
#include <iostream>

using namespace sc;
using namespace mpqc;
using namespace mpqc::TA;

using Matrix = ClDFGEngine::TAMatrix;

int main(int argc, char** argv){
  sc::ExEnv::init(argc, argv);
  mpqc::MADNESSRuntime::initialize();

  Ref<World> world = new World();

  const char *input = "./benzene_trimer.kv";
  Ref<KeyVal> kv = new ParsedKeyVal(input);

  Ref<Molecule> mol; mol << kv->describedclassvalue("benzene_trimer");


  Ref<TiledBasisSet> tbs; tbs << kv->describedclassvalue("basis");
  Ref<TiledBasisSet> dftbs; dftbs << kv->describedclassvalue("dfbasis");

  Ref<IntegralLibint2> ints_fac; ints_fac << kv->describedclassvalue("integrals");

  std::array<TiledArray::TiledRange1, 2>
        blocking{{tbs->trange1(), tbs->trange1()}};

  TiledArray::TiledRange trange(blocking.begin(), blocking.end());

  Matrix dens(*world->madworld(), trange);
  dens.set_all_local(1.0);

  ClDFGEngine G(ints_fac, tbs, dftbs, &dens, world);

  Matrix Gmat = G("i,j");
  world->madworld()->gop.fence();

  double t1 = madness::wall_time();
  Matrix Gmat2 = G("i,j");
  world->madworld()->gop.fence();
  double t2 = madness::wall_time();
  std::cout << "The formation of G took " << t2 - t1 << " seconds" << std::endl;


  mpqc::MADNESSRuntime::finalize();
  return 0;
}




