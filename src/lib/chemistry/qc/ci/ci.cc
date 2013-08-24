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
  int count;
  FromStateIn(config_, s, count);
}

CI::CI(const Ref<KeyVal> &kv)
    : ManyBodyWavefunction(kv) {
  {
    typedef KeyValValueint Int;

    /// only SD_RefWavefunction is tested for now
    Ref<SD_RefWavefunction> sd_refwfn; sd_refwfn << refwfn();
    if (sd_refwfn.null())
      throw InputError("only SD_RefWavefunction has been tested with CI",
                       __FILE__, __LINE__, "reference");

    const int charge =
        kv->intvalue("total_charge", Int(molecule()->total_Z() - refwfn()->nelectron()));
    const int magmom =
        kv->intvalue("magnetic_moment", Int(refwfn()->magnetic_moment()));
    {
        const int nelectron = molecule()->total_Z() - charge;
        if (nelectron%2 != magmom%2)
            throw InputError("charge and magnetic_moment inconsistent",
                             __FILE__, __LINE__);
    }

    config_.core = refwfn()->occ()->rank() - refwfn()->occ_act()->rank();
    config_.orbitals = refwfn()->occ_act()->rank() + refwfn()->uocc_act()->rank();

    const int nelectron = molecule()->total_Z() - charge - 2 * config_.core;
    config_.electrons.alpha = (nelectron + magmom ) / 2;
    config_.electrons.beta =  (nelectron - magmom ) / 2;

    config_.rank = kv->intvalue("rank", Int(0));
    // if (config_.rank != 0)
    //   throw FeatureNotImplemented("only max_ex_rank=0 currently supported (full CI)",
    //                               __FILE__, __LINE__);

    config_.max = kv->intvalue("max", Int(30));
    config_.collapse = kv->intvalue("collapse", Int(config_.collapse));
    config_.cutoff = kv->intvalue("cutoff", Int(config_.cutoff));
    config_.block = kv->intvalue("block", Int(config_.block));

    config_.convergence = this->desired_value_accuracy();
    config_.e_ref = molecule()->nuclear_repulsion_energy();
  }
}

CI::~CI() {}

void
CI::save_data_state(StateOut& so) {
  ManyBodyWavefunction::save_data_state(so);
  int count;
  ToStateOut(config_, so, count);
}

void
CI::print(std::ostream&o) const {
  using std::endl;
  o << indent << "CI:" << endl;
  o << incindent;
    config_.print(o);
    ManyBodyWavefunction::print(o);
  o << decindent;
}

RefSymmSCMatrix CI::density() {
  return 0;
}

void CI::compute() {
  E_ = CI::compute(ManyBodyWavefunction::refwfn(), config_);
  this->set_energy(E_.back() + config_.e_ref + config_.e_core);
}

int CI::value_implemented() const {
  return 1;
}

int CI::nelectron() {
  return config_.electrons.alpha + config_.electrons.beta;
}

double CI::magnetic_moment() const {
  return config_.electrons.alpha - config_.electrons.beta;
}

// /////////////////////////////////////////////////////////////////////////////

// // Local Variables:
// // mode: c++
// // c-file-style: "CLJ"
// // End:
