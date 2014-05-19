//
// tascf.cpp
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
#include <chemistry/qc/scf/tascf.hpp>
#include <tiledarray.h>
#include <chemistry/qc/basis/taskintegrals.hpp>
#include <chemistry/qc/basis/integralenginepool.hpp>
#include <math/elemental/eigensolver.hpp>

using namespace mpqc;
using namespace mpqc::TA;

sc::ClassDesc mpqc::TA::SCF::class_desc_(typeid(mpqc::TA::SCF), "TA.SCF",
                      1, "public TA.Wavefunction", 0, 0, 0);

mpqc::TA::SCF::SCF(const sc::Ref<sc::KeyVal>& kval) :
    Wavefunction(kval),
    MO_eigensystem_(this),
    ao_fock_(this)
{
    if(kval->exists("maxiter")){
        maxiter_= kval->intvalue("maxiter", sc::KeyValValueint(1000));
    }
    if(kval->exists("miniter")){
        miniter_= kval->intvalue("miniter", sc::KeyValValueint(0));
    }
    double desired_accuracy_ = 1e-6;
    if(kval->exists("accuracy")){
      desired_accuracy_ = kval->intvalue("accuracy",
                                           sc::KeyValValuedouble(1e-6));
    }

    ao_fock_.compute() = 0;
    ao_fock_.computed() = 0;
    ao_fock_.set_desired_accuracy(desired_accuracy_);

    MO_eigensystem_.compute() = 0;
    MO_eigensystem_.computed() = 0;
    MO_eigensystem_.set_desired_accuracy(desired_accuracy_);
}

SCF::~SCF(){}

SCF::TAMatrix&
SCF::ao_fock(double desired_accuracy){

  // Set the desired accuracy
  ao_fock_.set_desired_accuracy(desired_accuracy);

  if(!ao_fock_.computed_to_desired_accuracy()){
    compute_ao_fock(desired_accuracy);
    world()->madworld()->gop.fence();
    ao_fock_.computed() = 1;
    ao_fock_.set_actual_accuracy(desired_accuracy);
  }

  return ao_fock_.result_noupdate();
}

SCF::TAMatrixExpr
SCF::ao_fock_expr(std::string input){
  return ao_fock()(input);
}

SCF::ElemVector
SCF::MO_eigenvalues(double desired_accuracy){
  return MO_eigensystem(desired_accuracy).first;
}

SCF::TAMatrix
SCF::MO_eigenvectors(double desired_accuracy){
  return MO_eigensystem(desired_accuracy).second;
}

SCF::ElemTAEigenSystem
SCF::MO_eigensystem(double desired_accuracy){

  MO_eigensystem_.set_desired_accuracy(desired_accuracy);

  if(!MO_eigensystem_.computed_to_desired_accuracy()){
    MO_eigensystem_.result_noupdate() =
        eigensolver_full_Coeff(ao_fock(desired_accuracy),
                               ao_overlap());
    MO_eigensystem_.computed() = 1;
    MO_eigensystem_.set_actual_accuracy(desired_accuracy);
  }

  return MO_eigensystem_.result_noupdate();
}

#warning "nelectron has a temporary solution that assumes the charge is neutral
size_t mpqc::TA::SCF::nelectron() const {
    return molecule()->total_Z(); // Temporary solution.
}

// Add things to print later
void SCF::print(std::ostream &o) const {
  Wavefunction::print(o);
}
