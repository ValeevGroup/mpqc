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
#include <mpqc/interfaces/tiledarray/array_ints.hpp>
#include <mpqc/integrals/integralenginepool.hpp>

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



const Matrix& mpqc::TA::CLHF::fock(){
    if(not fock_.computed()){
        fock_ = ao_fock();
        fock_.computed() = 1;
    }
    return fock_.result_noupdate();
}

Matrix mpqc::TA::CLHF::ao_fock() {
    Matrix F = hcore()("m,n") + Gmat()("m,n");
    return F;
}

#warning "Gmat uses all four centered ints"
Matrix mpqc::TA::CLHF::Gmat(){
    mpqc::IntegralEnginePool<sc::Ref<sc::TwoBodyInt> > eri_pool(
                                            integral_->electron_repulsion());
    ::TiledArray::Array<double, 4> eri = mpqc::Integrals(*world_->madworld(),
                                                eri_pool, tbs_);
    Matrix Gmat = 2*rdm1()("r,s") * eri("m,r,n,s") -
                    rdm1()("r,s") * eri("m,r,s,n");
    return Gmat;
}

void mpqc::TA::CLHF::minimize_energy() {
    // Get matrices we need.
    const Matrix& H = hcore();
    const Matrix& S = overlap();
    Matrix F = fock_.result_noupdate()("i,j");


    double error_norminf = 1.0;

    size_t iter = 0;
    while(error_norminf > 1e-8){
        // Generate new Matrices
        Matrix D = Dguess(F);
        tr_corr_purify(D);
        rdm1_.result_noupdate() = D("i,j"); // Set to current density

        F("m,n") = H("m,n") + Gmat()("m,n");
        Matrix gradient = 8 * (S("i,q") * D("q,x") * F("x,j") -
                               F("i,q") * D("q,x") * S("x,j"));
        error_norminf = ::TiledArray::expressions::norminf(gradient("i,j"));
        error_norminf = abs(error_norminf);
        diis.extrapolate(F,gradient);

        world_->madworld()->gop.fence();
    }
    fock_.result_noupdate() = F("i,j");
}
