//
// ri_basis.cc
//
// Copyright (C) 2004 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#include <stdexcept>
#include <sstream>
#include <cassert>

#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/r12technology.h>
#include <math/scmat/svd.h>
#include <chemistry/qc/lcao/transform_factory.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <util/misc/print.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>

using namespace sc;
using namespace std;

void
R12WavefunctionWorld::construct_ri_basis_(bool safe)
{
  if (ribs_space_)
    return;

  Ref<GaussianBasisSet> obs = basis();
  const bool obs_eq_abs = bs_aux_->equiv(obs);
  const bool vbs_eq_abs = obs_eq_vbs_ ? obs_eq_abs : bs_aux_->equiv(basis_vir());
  if (obs_eq_abs) {
    bs_ri_ = obs;
  }
  else {

    switch(r12tech()->abs_method()) {

      case R12Technology::ABS_CABS:
        bs_ri_ = bs_aux_;
        break;

      case R12Technology::ABS_CABSPlus:
      {
        bs_ri_ = bs_aux_ + obs;
        if (!vbs_eq_abs && !obs_eq_vbs_)
          bs_ri_ = bs_ri_ + basis_vir();
      }
      break;

      default:
        throw std::logic_error("R12WavefunctionWorld::construct_ri_basis_ -- invalid abs_method");

    }
    construct_orthog_ri_();
  }

  {
    // also create AO space for RI basis
    Ref<OrbitalSpaceRegistry> idxreg = world()->tfactory()->orbital_registry();
    Ref<AOSpaceRegistry> aoidxreg = world()->tfactory()->ao_registry();
    Ref<Integral> localints = refwfn()->integral()->clone();
    if (!obs_eq_ribs()) { // RI-BS
      Ref<OrbitalSpace> mu = new AtomicOrbitalSpace("mu'", "RIBS(AO)", basis_ri(), localints);
      idxreg->add(make_keyspace_pair(mu));
      aoidxreg->add(mu->basis(),mu);
    }
  }
}

void
R12WavefunctionWorld::construct_cabs_()
{
  if (cabs_space_[Alpha])
    return;

  construct_ri_basis_(r12tech()->safety_check());

  Ref<GaussianBasisSet> obs = refwfn()->basis();
  if (bs_ri_->equiv(obs))
    throw std::logic_error("R12WavefunctionWorld::construct_cabs_ -- CABS and CABS+ methods can only be used when ABS != OBS");

  ref_acc_for_cabs_space_ = refwfn()->desired_value_accuracy();
  construct_ortho_comp_svd_();
}

void
R12WavefunctionWorld::construct_orthog_aux_()
{
  if (abs_space_)
    return;

  if (! this->basis()->equiv(bs_aux_) &&
      (orthog_method() == OverlapOrthog::Symmetric ||
       orthog_method() == OverlapOrthog::Canonical) &&
      this->r12tech()->abs_nlindep() >= 0) {  // if R12Technology has abs_nlindep, use it only with symmetric/canonical orthogonalization
    nlindep_aux_ = this->r12tech()->abs_nlindep();
  }
  else // will compute the number of linear dependencies automatically
    nlindep_aux_ = -1;
  abs_space_ = orthogonalize("p'","RIBS", bs_aux_, integral(), orthog_method(), r12tech()->abs_lindep_tol(), nlindep_aux_);
  if (bs_aux_ == bs_ri_) {
    ribs_space_ = abs_space_;
  }
}

/* WARNING R12WavefunctionWorld::construct_orthog_vir_() moved to r12wfnworld.o -- gcc 3.4.3 generates internal symbols
   for SpinSpaces. Bug? If yes and fixed -- will move the function back here.
*/

void
R12WavefunctionWorld::construct_orthog_ri_()
{
  if (bs_ri_.null())
    throw std::runtime_error("R12WavefunctionWorld::construct_orthog_ri_ -- RI basis has not been set yet");
  if (bs_aux_ == bs_ri_)
    construct_orthog_aux_();
  if (ribs_space_.null()) {

    if (! this->basis()->equiv(bs_ri_) &&
        (orthog_method() == OverlapOrthog::Symmetric ||
         orthog_method() == OverlapOrthog::Canonical) &&
        this->r12tech()->abs_nlindep() >= 0) {  // if R12Technology has abs_nlindep, use it only with symmetric/canonical orthogonalization
      nlindep_ri_ = this->r12tech()->abs_nlindep();
    }
    else  // will compute the number of linear dependencies automatically
      nlindep_ri_ = -1;

    ribs_space_ = orthogonalize("p'","RIBS", bs_ri_, integral(), orthog_method(), r12tech()->abs_lindep_tol(), nlindep_ri_);
  }
  const Ref<OrbitalSpaceRegistry> idxreg = this->world()->tfactory()->orbital_registry();
  idxreg->add(make_keyspace_pair(ribs_space_));
}

bool
R12WavefunctionWorld::abs_spans_obs_()
{
  construct_orthog_aux_();

  // Compute the bumber of linear dependencies in OBS+ABS
  GaussianBasisSet& abs = *(bs_aux_.pointer());
  Ref<GaussianBasisSet> ri_basis = abs + refwfn()->basis();
  int nlindep_ri;
  if (bs_ri_ && ri_basis->equiv(bs_ri_)) {
    construct_orthog_ri_();
    nlindep_ri = nlindep_ri_;
  }
  else {
    nlindep_ri = -1;  // will compute the number of linear dependencies automatically
    Ref<OrbitalSpace> ribs_space = orthogonalize("p+p'","OBS+ABS", ri_basis, integral(), orthog_method(), r12tech()->abs_lindep_tol(), nlindep_ri);
  }

  const int obs_rank = refwfn()->orbs(Alpha)->rank();
  if (nlindep_ri - nlindep_aux_ - obs_rank == 0)
    return true;
  else
    return false;
}

/////////////////////////////////////////////////////////////////////////////

void
R12WavefunctionWorld::construct_ortho_comp_svd_()
{
  const Ref<OrbitalSpaceRegistry> idxreg = this->world()->tfactory()->orbital_registry();

  const double tol = r12tech()->abs_lindep_tol();
  if (!refwfn()->spin_polarized()) {
    Ref<OrbitalSpace> tmp = orthog_comp(refwfn()->occ_sb(Alpha), ribs_space_, "p'-m", "CABS", tol);
    tmp = orthog_comp(refwfn()->uocc_sb(Alpha), tmp, "a'", "CABS", tol);
    cabs_space_[Alpha] = tmp;
    cabs_space_[Beta] = cabs_space_[Alpha];
    idxreg->add(make_keyspace_pair(cabs_space_[Alpha]));
  }
  else {
    Ref<OrbitalSpace> tmp_a = orthog_comp(refwfn()->occ_sb(Alpha), ribs_space_, "P'-M", "CABS (Alpha)", tol);
    Ref<OrbitalSpace> tmp_b = orthog_comp(refwfn()->occ_sb(Beta), ribs_space_, "p'-m", "CABS (Beta)", tol);
    if (USE_NEW_ORBITALSPACE_KEYS) {
      const std::string key_a = ParsedOrbitalSpaceKey::key(std::string("a'"),Alpha);
      const std::string key_b = ParsedOrbitalSpaceKey::key(std::string("a'"),Beta);
      cabs_space_[Alpha] = orthog_comp(refwfn()->uocc_sb(Alpha), tmp_a, key_a, "CABS (Alpha)", tol);
      cabs_space_[Beta] = orthog_comp(refwfn()->uocc_sb(Beta), tmp_b, key_b, "CABS (Beta)", tol);
    }
    else // old orbitalspace key no longer supported
      MPQC_ASSERT(false);

    idxreg->add(make_keyspace_pair(cabs_space_[Alpha]));
    idxreg->add(make_keyspace_pair(cabs_space_[Beta]));
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
