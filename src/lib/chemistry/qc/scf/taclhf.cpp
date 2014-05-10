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
#include <tiled_array.h>
#include <TiledArray/algebra/diis.h>
#include <chemistry/qc/basis/integralenginepool.hpp>
#include <chemistry/qc/basis/taskintegrals.hpp>
#include <chemistry/qc/scf/cldfgengine.hpp>
#include <math/elemental/eigensolver.hpp>
#include <util/misc/regtime.h>

using namespace mpqc;
using namespace mpqc::TA;
using TAMatrix = mpqc::TA::CLHF::TAMatrix;

sc::ClassDesc mpqc::TA::CLHF::class_desc_(typeid(mpqc::TA::CLHF), "TA.CLHF",
                      1, "public TA.CLSCF",
                      0,
                      0,
                      0);

mpqc::TA::CLHF::CLHF(const sc::Ref<sc::KeyVal>& kval) :
    CLSCF(kval), G_eng()
{
  G_eng << kval->describedclassvalue("GEngine");
  if(G_eng.null()){ // If no GEngine was given then use Density Fitting
    if(kval->exists("dfbasis")){ // check if density fitting basis was given
      sc::ExEnv::out0() <<
              "Constructing a default GEngine: using density fitting\n";
      G_eng = new ClDFGEngine(kval);
    }
    else{ // else if no dfbasis use default
      throw sc::InputError("No dfbasis was given to GEngine.");
    }
  }
}

GEngineBase::return_type
mpqc::TA::CLHF::G(const std::string s){
  // Prefer coefficent build, but only possible if we have a fock matrix
  if(!G_eng->coefficients_set() && SCF::scf_ao_fock_().is_initialized()){
    if(!Coeff_.is_initialized()){
      Coeff_ = eigensolver_occ_Coeff(scf_ao_fock_(), ao_overlap(), occupation());
    }
    G_eng->set_coefficients({&Coeff_});
  }
  else if(!G_eng->densities_set() ){
    G_eng->set_densities({&ao_density()});
  }

  return G_eng->operator()(s);
}

void CLHF::compute_ao_fock(double desired_accuracy) {
  sc::Timer tim("SCF iteration:");
  TAMatrix &F = scf_ao_fock_();
  const TAMatrix &S = ao_overlap();
  const TAMatrix &H = ao_hcore();
  TAMatrix &D = ao_density();

  double error_norminf = desired_accuracy + 1.0; // Just to ensure we do a loop
  TiledArray::DIIS<TAMatrix> diis;
  size_t iter = 0;

  std::cout << "Beginning scf iterations:" << std::endl;
  while(error_norminf > desired_accuracy && iter >= miniter() && iter < maxiter()){
    sc::Timer scf_tim;

    scf_tim.enter("Coeff_ eigensolver");
    Coeff_ = eigensolver_occ_Coeff(F, S, occupation());
    D = Coeff_("mu,i") * Coeff_("nu,i");
    scf_tim.exit("Coeff_ eigensolver");

    scf_tim.enter("Fock matrix Contraction");
    // Update the Fock matrix
    F("i,j") = H("i,j") + G("i,j");
    scf_tim.exit("Fock matrix Contraction");

    scf_tim.enter("Gradient and diis");
    // Compute the gradient and extrapolate with diis
    TAMatrix Grad = 8 * ( S("i,q") * D("q,x") * F("x,j") -
                          F("i,q") * D("q,x") * S("x,j") );
    diis.extrapolate(F,Grad);
    scf_tim.exit("Gradient and diis");

    // Compute the error as the largest element of the gradient.
    error_norminf = ::TiledArray::expressions::norminf(Grad("i,j"));
    world()->madworld()->gop.fence();
    if(world()->madworld()->rank() == 0){
      std::cout << "\tIteration " << iter++ << "\n";
      std::cout << "\t\tenergy = " << iter_energy() << "\n";
      std::cout << "\t\terror = " << error_norminf << "\n" << std::endl;
    }
  }
  tim.exit("SCF iteration:");
}

// if the base class SCF::scf_fock() matrix has not been initialized then
// we should construct it.  We call the base class to get to the tiledarray
TAMatrix& mpqc::TA::CLHF::scf_ao_fock_(){

  if(!SCF::scf_ao_fock_().is_initialized()){
    // Initialize fock with Gengine
    SCF::scf_ao_fock_() = ao_hcore()("i,j") + G("i,j");
    world()->madworld()->gop.fence();
  }

  return SCF::scf_ao_fock_();
}
