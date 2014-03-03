//
// taclscf.cpp
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

#include <chemistry/qc/scf/taclscf.hpp>
#include <mpqc/interfaces/tiledarray/array_ints.hpp>
#include <mpqc/integrals/integralenginepool.hpp>
#include <Eigen/Dense> // To be removed once tiledarray has operations on diagonal.

using namespace mpqc;
using namespace mpqc::TA;
using Matrix = mpqc::TA::CLSCF::Matrix;

sc::ClassDesc mpqc::TA::CLSCF::class_desc_(typeid(mpqc::TA::CLSCF), "TA.CLSCF",
                      1, "public TA.SCF",
                      0,
                      sc::create<mpqc::TA::CLSCF>,
                      0);

mpqc::TA::CLSCF::CLSCF(const sc::Ref<sc::KeyVal>& kval) :
    SCF(kval), diis(), hcore_(this)
{
    hcore_.compute() = 0;
    hcore_.computed() = 0;

    if(nelectron() % 2 != 0){
        throw sc::InputError("Number of electrons is not divisable by two",
                             __FILE__, __LINE__, "", "", this->class_desc());
    }
    occupation = nelectron()/2;
}

mpqc::TA::CLSCF::~CLSCF() {}

#warning "compute is not yet defined"
void mpqc::TA::CLSCF::compute() {
    MPQC_ASSERT(false);
}

const Matrix& mpqc::TA::CLSCF::rdm1() {
    if(not rdm1_.computed()){
        Matrix temp = Dguess(hcore()); // Use hcore as guess till SOAD works
        rdm1_= temp("i,j");
        rdm1_.computed() = 1;
    }
    return rdm1_.result_noupdate();
}

Matrix& mpqc::TA::CLSCF::hcore(){
    if(not hcore_.computed()){
        mpqc::IntegralEnginePool<sc::Ref<sc::OneBodyInt> > hcore_pool(
                                                            integral_->hcore());
        hcore_ = mpqc::Integrals(*world_->madworld(), hcore_pool, tbs_);
        world_->madworld()->gop.fence();
        hcore_.computed() = 1;
    }
    return hcore_.result_noupdate();
}

double mpqc::TA::CLSCF::scf_energy() {
    // E = \sum_{ij} \left( D_{ij} * (F_{ij} + H_{ij}) \right)
    return ::TiledArray::expressions::dot(
                    hcore()("i,j") + fock()("i,j"), rdm1()("i,j"));
}

#warning "tr_corr_purify uses Eigen and is not production ready"
void mpqc::TA::CLSCF::tr_corr_purify(Matrix& P) {
    // Avoid Eigen in the future
    Eigen::MatrixXd Ep = ::TiledArray::array_to_eigen(P);

    // Purify if the matrix is not equal to it's square then purify
    while(Eigen::MatrixXd(Ep - Ep*Ep).lpNorm<Eigen::Infinity>() >= 1e-10){
        // If the trace of the matrix is too large shrink it else raise it
        Ep = (Ep.trace() >= occ()) ? Eigen::MatrixXd(Ep*Ep) :
                                            Eigen::MatrixXd(2*Ep - Ep * Ep);
    }
    P = ::TiledArray::eigen_to_array<Matrix>(*world_->madworld(), P.trange(),
                                             Ep);
}

#warning "Dguess uses Eigen and is not production ready"
Matrix mpqc::TA::CLSCF::Dguess(const Matrix& F){
    // Should be using Eval guesses instead of Frobenius norm, but this will work
    double Fnorm = ::TiledArray::expressions::norm2(F("i,j"));

    /* Needs to be changed to TiledArray only opps, but for now this is a
    *  stand in.
    */
    Eigen::MatrixXd Ef = ::TiledArray::array_to_eigen(F);

    // Shift spectrum of F
    for(size_t i = 0; i < Ef.rows(); ++i){
        Ef(i,i) = Fnorm - Ef(i,i);
    }
    Matrix D = ::TiledArray::eigen_to_array<Matrix>(*world_->madworld(),
                                                F.trange(),
                                                Ef);
    /* End part that needs to be replaced */

    // Scale Evals to the range (0,1)
    D("i,j") = D("i,j") * (1.0 / (2.0 * Fnorm));

    // Purifiy to idempotency
    tr_corr_purify(D);

    return D;
}
