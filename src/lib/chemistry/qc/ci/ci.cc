//
// ci.cc
//
// Copyright (C) 2012 Andrey Asadchev
//
// Author: Andrey Asadchev <asadchev@vt.edu>
// Maintainer: EFV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/ci/ci.h>

namespace sc {
    static ClassDesc CI_cd(typeid(CI), "CI", 1, "public ManyBodyWavefunction", 0,
                           create<CI>, create<CI>);
}

using namespace sc;

CI::CI(StateIn& s) :
    ManyBodyWavefunction(s)
{
}

CI::CI(const Ref<KeyVal> &kv)
    : ManyBodyWavefunction(kv), kv_(kv) {}

CI::~CI() {}

void
CI::save_data_state(StateOut& so) {
  ManyBodyWavefunction::save_data_state(so);
}

RefSymmSCMatrix CI::density() {
  return 0;
}

void CI::compute() {

  mpqc::ci::Config config;
  {
    size_t ne = refwfn()->nelectron();
    size_t no = refwfn()->basis()->nbasis();
    typedef KeyValValueint Int;
    config.core = kv_->intvalue("core", Int(0));
    config.orbitals = kv_->intvalue("orbitals", Int(no));
    config.alpha = kv_->intvalue("alpha", Int((ne + 1) / 2));
    config.beta = kv_->intvalue("beta", Int(ne / 2));
    config.level = kv_->intvalue("level", Int(0));
    config.max = kv_->intvalue("max", Int(30));
    config.collapse = kv_->intvalue("collapse", Int(config.collapse));
    config.cutoff = kv_->intvalue("cutoff", Int(config.cutoff));
    config.block = kv_->intvalue("block", Int(config.block));
    config.convergence = this->desired_value_accuracy();
    config.e_ref = molecule()->nuclear_repulsion_energy();
  }

  E_ = CI::compute(ManyBodyWavefunction::refwfn(), config);
  this->set_energy(E_.back() + config.e_ref);
}

int CI::value_implemented() const {
  return 1;
}

// /////////////////////////////////////////////////////////////////////////////

// // Local Variables:
// // mode: c++
// // c-file-style: "CLJ"
// // End:
