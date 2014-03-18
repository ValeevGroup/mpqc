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
    world()->madworld()->gop.fence();
    return Gmat;
}

void mpqc::TA::CLHF::minimize_energy() {
  Matrix &F = scf_fock();
  const Matrix &S = overlap();
  const Matrix &H = hcore();
  Matrix &D = density();

  size_t iter = 0;
  double error_norminf = 1.0;

  std::cout << "S = \n" << S << std::endl;
  std::cout << "H = \n" << H << std::endl;
  while(error_norminf > 1e-6 && iter < 100){
    Dguess(F);  // modifies D internally.
    std::cout << "D = \n" << D << std::endl;
    F("i,j") = H("i,j") + Gmat()("i,j");
    std::cout << "F = \n" << F << std::endl;
    Matrix Grad = 8 * ( S("i,q") * D("q,x") * F("x,j") -
                        F("i,q") * D("q,x") * S("x,j") );
    error_norminf = ::TiledArray::expressions::norminf(Grad("i,j"));
    //diis.extrapolate(F, Grad);
    world()->madworld()->gop.fence();
    std::cout << "Energy for iter " << ++iter << " = " << iter_energy() << std::endl;
    std::cout << "Error = " << error_norminf << std::endl;
  }
}

Matrix& mpqc::TA::CLHF::scf_fock(){
  if(!SCF::scf_fock().is_initialized()){
    Matrix& F = SCF::scf_fock();
    F = hcore()("i,j") + Gmat()("i,j");
    world()->madworld()->gop.fence();
  }
  return SCF::scf_fock();
}
