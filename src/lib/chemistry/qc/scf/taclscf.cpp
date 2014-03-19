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
#include <chemistry/qc/basis/integralenginepool.hpp>
#include <chemistry/qc/basis/taskintegrals.hpp>
#include <chemistry/qc/lcao/soad.h>
#include <mpqc/interfaces/tiledarray/symmscmat.hpp>
#include <util/elemental/ta_interface.hpp>
#include <util/elemental/grid.hpp>
#include <Eigen/Dense> // To be removed once tiledarray has operations on diagonal.
using namespace mpqc;
using namespace mpqc::TA;
using Matrix = mpqc::TA::CLSCF::Matrix;

sc::ClassDesc mpqc::TA::CLSCF::class_desc_(typeid(mpqc::TA::CLSCF), "TA.CLSCF",
                                           1, "public TA.SCF", 0,
                                           sc::create<mpqc::TA::CLSCF>, 0);

mpqc::TA::CLSCF::CLSCF(const sc::Ref<sc::KeyVal>& kval) :
        SCF(kval) {
  if (nelectron() % 2 != 0) {
    throw sc::InputError("Number of electrons is not divisable by two",
    __FILE__,
                         __LINE__, "", "", this->class_desc());
  }
  occ() = nelectron() / 2;
}

mpqc::TA::CLSCF::~CLSCF() {
}

#warning "compute is not yet defined"
void mpqc::TA::CLSCF::compute() {
  MPQC_ASSERT(false);
}

const Matrix& mpqc::TA::CLSCF::rdm1() {
  MPQC_ASSERT(false);
}

double mpqc::TA::CLSCF::scf_energy() {
  // E = \sum_{ij} \left( D_{ij} * (F_{ij} + H_{ij}) \right)
  return ::TiledArray::expressions::dot(hcore()("i,j") + fock()("i,j"),
                                        rdm1()("i,j"));
}

#warning "tr_corr_purify uses Eigen and is not production ready"
void mpqc::TA::CLSCF::tr_corr_purify(Matrix &P) {
  // Avoid Eigen in the future
  Eigen::MatrixXd Ep = ::TiledArray::array_to_eigen(P);
  Eigen::MatrixXd Es = ::TiledArray::array_to_eigen(overlap());

  // Purify if the matrix is not equal to it's square then purify
  while (Eigen::MatrixXd(Ep - Ep*Es*Ep).lpNorm<Eigen::Infinity>() >= 1e-10) {
    // If the trace of the matrix is too large shrink it else raise it
    Ep = (Ep.trace() >= occupation()) ? Eigen::MatrixXd(Ep*Es*Ep) :
                                        Eigen::MatrixXd(2 * Ep - Ep*Es*Ep);
  }
  P = ::TiledArray::eigen_to_array < Matrix
          > (*(world())->madworld(), P.trange(), Ep);
}

// If we have not generated an intial guess for the density yet then do so now.
Matrix& mpqc::TA::CLSCF::density() {
  // Check to see if data has been initialized if
  if (!Wavefunction::density().is_initialized()) {
    // Get a reference to the array for ease of compuation
    Matrix &D = Wavefunction::density();

    // For constructing a SOAD object
    sc::Ref<sc::AssignedKeyVal> akv = new sc::AssignedKeyVal;
    akv->assign("molecule", molecule().pointer());
    akv->assign("basis", basis().pointer());

    // Construct a density guess based on a SOAD object
    using Soad = sc::SuperpositionOfAtomicDensities;
    sc::Ref<Soad> guess = new Soad(sc::Ref<sc::KeyVal>(akv));

    // Copy the mpqc sc matrix into our tiledarray Matrix.
    D = mpqc::SymmScMat_To_TiledArray(*world()->madworld(),
                                      guess->guess_density(basis(),integral()),
                                      basis()->trange1());
    world()->madworld()->gop.fence();
    std::cout << "D = \n" << D << std::endl;
    world()->madworld()->gop.fence();
  }
  return Wavefunction::density();
}

#warning "Dguess uses Eigen and is not production ready"
void mpqc::TA::CLSCF::Dguess(const Matrix& F) {

  // Grab density so we can work with it.
  Matrix& D = density();

  /* Needs to be changed to TiledArray only opps, but for now this is a
   *  stand in.
   */
  Eigen::MatrixXd Ef = ::TiledArray::array_to_eigen(F);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ef);
  double emin = es.eigenvalues().minCoeff();
  double emax = es.eigenvalues().maxCoeff();

  // Move to elem
  mpqc::Grid* grid = new mpqc::Grid();
  elem::DistMatrix<double> ElemF = array_to_distmat(F, grid->elemGrid());
  elem::Print(ElemF);
  std::cout << "F \n" << F << std::endl;

  // Shift spectrum of F
  for (size_t i = 0; i < Ef.rows(); ++i) {
    Ef(i, i) = emax - Ef(i, i);
  }
  D=::TiledArray::eigen_to_array<Matrix>(*world()->madworld(), F.trange(), Ef);
  /* End part that needs to be replaced */

  // Scale Evals to the range (0,1)
  D("i,j") = D("i,j") * (1.0 / (emax - emin));

  // Purifiy to idempotency
  tr_corr_purify(D);
}

double mpqc::TA::CLSCF::iter_energy() {
  // E = \sum_{ij} \left( D_{ij} * (F_{ij} + H_{ij}) \right)
  return ::TiledArray::expressions::dot(hcore()("i,j") + scf_fock()("i,j"),
                                        density()("i,j"));
}
