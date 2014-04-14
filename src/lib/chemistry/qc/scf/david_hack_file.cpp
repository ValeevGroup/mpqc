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
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/lcao/soad.h>
#include <chemistry/qc/libint2/libint2.h>
#include <util/group/pregtime.h>
#include <util/madness/init.h>
#include <math/elemental/eigensolver.hpp>
#include <TiledArray/algebra/diis.h>
#include <iostream>

using namespace sc;
using namespace mpqc;
using namespace mpqc::TA;

using Matrix = TiledArray::Array<double, 2, TiledArray::Tensor<double> >;
using DistMatrix = elem::DistMatrix<double>;

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
  ints_fac->set_basis(tbs);


  if(world->madworld()->rank()==0){
    mol->print();
    tbs->print();
    dftbs->print();
  }


  // Get pools
  using onebpool = IntegralEnginePool<Ref<OneBodyInt> >;
  auto S_pool = std::make_shared<onebpool>(ints_fac->overlap()->clone());
  auto H_pool = std::make_shared<onebpool>(ints_fac->hcore()->clone());

  // Compute ints
  Matrix S = Integrals(*world->madworld(), S_pool, tbs);
  world->madworld()->gop.fence();
  Matrix H = Integrals(*world->madworld(), H_pool, tbs);
  world->madworld()->gop.fence();

  std::array<TiledArray::TiledRange1, 2>
        blocking{{tbs->trange1(), tbs->trange1()}};

  TiledArray::TiledRange trange(blocking.begin(), blocking.end());

  // Intialize Density matrix
  Matrix dens(*world->madworld(), trange);
  dens.set_all_local(0.0);
  world->madworld()->gop.fence();

  // Guess for density
  double etimer0 = madness::wall_time();
  eigensolver_D(H, S, dens, occ);
  world->madworld()->gop.fence();
  double etimer1 = madness::wall_time();
  if(world->madworld()->rank()==0){
    std::cout << "Total time to form density = " << etimer1 - etimer0 << " s \n" << std::endl;
  }


  // Initialize G engine
  ClDFGEngine G(ints_fac, tbs, dftbs, &dens, world);

  // Compute integrals
  G("i,j");
  world->madworld()->gop.fence();

  // Make Fock matrix
  double t1 = madness::wall_time();
  Matrix F = H("i,j") + G("i,j");
  world->madworld()->gop.fence();
  double t2 = madness::wall_time();
  if(world->madworld()->rank()==0){
    std::cout << "The formation of F took " << t2 - t1 << " seconds" << std::endl;
  }

  // Perform SCF iterations
  double energyinit = 1;
  double energy = 0;
  int iter = 0;
  TiledArray::DIIS<Matrix> diis;
  world->madworld()->gop.fence();
  while(abs(energyinit - energy) >= 1e-9 && ++iter < 100){
    double scf0 = madness::wall_time();
    energyinit = energy;

    eigensolver_D(F, S, dens, occ);
    world->madworld()->gop.fence();
    auto f0 = madness::wall_time();
    F = H("i,j") + G("i,j");
    world->madworld()->gop.fence();
    auto f1 = madness::wall_time();

    Matrix gradient = 8 * (S("i,q") * dens("q,x") * F("x,j") -
                           F("i,q") * dens("q,x") * S("x,j"));

    diis.extrapolate(F, gradient);

    energy = TiledArray::expressions::dot((H("i,j") + F("i,j")), dens("i,j")) +
            mol->nuclear_repulsion_energy();
    world->madworld()->gop.fence();
    double scf1 = madness::wall_time();

    if(world->madworld()->rank()==0){
      std::cout << "SCF iteration " << iter << "\n\tTime = " << scf1 - scf0 <<
              " s\n\tFock build time = " << f1-f0 <<
              " s\n\tEnergy = " << energy << std::endl;
    }
  }

  std::cout << "Finished calculation energy = " <<  energy << std::endl;
  world->madworld()->gop.fence();

  Ref<AssignedKeyVal> kv_old = new AssignedKeyVal();
  kv_old->assign("basis", tbs.pointer());
  kv_old->assign("molecule", mol.pointer());
  kv_old->assign("integrals", ints_fac.pointer());

  Ref<SuperpositionOfAtomicDensities> soad_guess = new
                      SuperpositionOfAtomicDensities(kv_old);

  kv_old->assign("guess_wavefunction", soad_guess.pointer());

  // Try and allocate clhf object
  Ref<MessageGrp> msg;
  MessageGrp::initial_messagegrp(argc, argv);
  msg = MessageGrp::get_default_messagegrp();
  Ref<RegionTimer> regtim;
  regtim = new ParallelRegionTimer(msg, "SCF compare", 1,1);
  RegionTimer::set_default_regiontimer(regtim);
  Ref<sc::CLHF> old_hf = new sc::CLHF(kv_old);
  old_hf->print();

  Timer tim;
  tim.enter("old scf");
  double old_hf_energy = old_hf->energy();
  world->madworld()->gop.fence();
  tim.exit("old scf");
  tim.print(ExEnv::out0());

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

