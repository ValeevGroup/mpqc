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
#include <chemistry/qc/scf/taclhf.hpp>
#include <chemistry/qc/lcao/soad.h>
#include <chemistry/qc/libint2/libint2.h>
#include <util/group/pregtime.h>
#include <util/madness/init.h>
#include <util/madness/world.h>
#include <chemistry/qc/basis/tiledbasisset.hpp>
#include <chemistry/qc/basis/integral.h>
#include <math/elemental/eigensolver.hpp>
#include <chemistry/qc/basis/integralenginepool.hpp>
#include <chemistry/qc/basis/taskintegrals.hpp>
#include <TiledArray/algebra/diis.h>
#include <mpqc/interfaces/tiledarray/symmscmat.hpp>
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

  Ref<mpqc::TA::CLHF> tclhf = new TA::CLHF(kv);
  Ref<MessageGrp> msg;
  MessageGrp::initial_messagegrp(argc, argv);
  msg = MessageGrp::get_default_messagegrp();
  Ref<RegionTimer> regtim;
  regtim = new ParallelRegionTimer(msg, "SCF compare", 1,1);
  RegionTimer::set_default_regiontimer(regtim);

  world->madworld()->gop.fence();

  std::cout << "About to compute scf_energy\n";
  sc::Timer tim("Scf Energy time");
  double energy = tclhf->scf_energy();
  tim.exit("Scf Energy time");
  if(world->madworld()->rank() == 0){
    std::cout << "The Final Energy was " << energy << std::endl;
  }
  tim.print();

 #if 0
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

  Ref<AssignedKeyVal> akv = new AssignedKeyVal();
  akv->assign("molecule", mol.pointer());
  akv->assign("basis", tbs.pointer());
  akv->assign("dfbasis", dftbs.pointer());
  akv->assign("world", world.pointer());
  akv->assign("integrals", ints_fac.pointer());

  // Construct a density guess based on a SOAD object
  using Soad = sc::SuperpositionOfAtomicDensities;
  sc::Ref<Soad> guess = new Soad(sc::Ref<sc::KeyVal>(akv));

  // Copy the mpqc sc matrix into our tiledarray Matrix.
  dens = mpqc::SymmScMat_To_TiledArray(*world->madworld(),
                                    guess->guess_density(tbs,ints_fac),
                                    tbs->trange1());

  ClDFGEngine GC(static_cast<Ref<KeyVal> >(akv));
  GC.set_densities({&dens});
  world->madworld()->gop.fence();

  // Make Fock matrix
  double t1 = madness::wall_time();
  Matrix F = H("i,j") + GC("i,j");
  world->madworld()->gop.fence();
  double t2 = madness::wall_time();
  if(world->madworld()->rank()==0){
    std::cout << std::setprecision(12) << std::endl;
    std::cout << "The formation of F took " << t2 - t1 << " seconds" << std::endl;
  }

  // give coefficients to GC
  Matrix C;
  GC.set_coefficients({&C});
  world->madworld()->gop.fence();

  // Perform SCF iterations
  double energy = 0;
  double error = 10;
  int iter = 0;
  TiledArray::DIIS<Matrix> diis(1,6);
  double ftime = 0;
  double scftime = 0;
  world->madworld()->gop.fence();
  while(error >= 1e-8 && ++iter < 100){
    double scf0 = madness::wall_time();

    C = eigensolver_occ_Coeff(F, S, occ);
    dens = C("mu,i")*C("nu,i");
    world->madworld()->gop.fence();


    auto f0 = madness::wall_time();
    F = H("i,j") + GC("i,j");
    world->madworld()->gop.fence();
    auto f1 = madness::wall_time();
    ftime += (f1 - f0);

    Matrix gradient = 8 * (S("i,q") * dens("q,x") * F("x,j") -
                           F("i,q") * dens("q,x") * S("x,j"));

    error = TiledArray::expressions::norminf(gradient("i,j"));

    diis.extrapolate(F, gradient);

    energy = TiledArray::expressions::dot((H("i,j") + F("i,j")), dens("i,j")) +
            mol->nuclear_repulsion_energy();
    world->madworld()->gop.fence();
    double scf1 = madness::wall_time();
    scftime += scf1 - scf0;

    if(world->madworld()->rank()==0){
      std::cout << "SCF iteration " << iter << "\n\tTime = " << scf1 - scf0 <<
              " s\n\tFock build time = " << f1 - f0  <<
              " s\n\tEnergy = " << energy <<
              " s\n\tGradient Norm = " << error << std::endl;
    }
  }


  std::cout << "Average scf iter time = " << scftime/(double(iter)) << "s in " <<
          iter << " iterations " << std::endl;

  std::cout << "Average Fock build time = " << ftime/(double(iter)) << "s in " <<
          iter << " iterations " << std::endl;

  std::cout << "Finished calculation energy = " <<  energy << std::endl;
  world->madworld()->gop.fence();
  */

#endif
#if 0
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
#endif

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

