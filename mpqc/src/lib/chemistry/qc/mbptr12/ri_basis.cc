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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <sstream>

#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/r12technology.h>
#include <chemistry/qc/mbptr12/svd.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/mbptr12/orbitalspace_utils.h>

using namespace sc;
using namespace std;

void
R12WavefunctionWorld::construct_ri_basis_(bool safe)
{
  // RI basis is only needed if corrfactor != none
  Ref<R12Technology::NullCorrelationFactor> null_cf; null_cf << r12tech()->corrfactor();
  const bool ri_basis_not_needed = null_cf.nonnull();

  Ref<GaussianBasisSet> obs = basis();
  const bool obs_eq_abs = bs_aux_->equiv(obs);
  const bool vbs_eq_abs = obs_eq_vbs_ ? obs_eq_abs : bs_aux_->equiv(basis_vir());
  if (obs_eq_abs) {
    bs_ri_ = obs;
  }
  else {
    if (ri_basis_not_needed) {
      bs_ri_ = bs_aux_;
    }
    else {
      switch(r12tech()->abs_method()) {

        case R12Technology::ABS_ABS:
          bs_ri_ = bs_aux_;
          if (!abs_spans_obs_()) {
            ExEnv::out0() << endl << indent << "WARNING: the auxiliary basis is not safe to use with the given orbital basis" << endl << endl;
            if (safe)
              throw std::runtime_error("R12WavefunctionWorld::construct_ri_basis_abs_ -- auxiliary basis is not safe to use with the given orbital basis");
          }
          break;

        case R12Technology::ABS_CABS:
          bs_ri_ = bs_aux_;
          break;

        case R12Technology::ABS_ABSPlus:
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
  }
}

void
R12WavefunctionWorld::construct_cabs_()
{
  // CABS space is only needed if corrfactor != none ...
  Ref<R12Technology::NullCorrelationFactor> null_cf; null_cf << r12tech()->corrfactor();
  const bool ri_basis_not_needed = null_cf.nonnull();

  if (!ri_basis_not_needed) {

    Ref<GaussianBasisSet> obs = ref()->basis();
    if (bs_ri_->equiv(obs) &&
        (r12tech()->abs_method() == R12Technology::ABS_CABS ||
         r12tech()->abs_method() == R12Technology::ABS_CABSPlus
        )
       )
      throw std::logic_error("R12WavefunctionWorld::construct_cabs_ -- CABS and CABS+ methods can only be used when ABS != OBS");

    ref_acc_for_cabs_space_ = ref()->desired_value_accuracy();
    construct_ortho_comp_svd_();
  }
}

void
R12WavefunctionWorld::construct_orthog_aux_()
{
  if (abs_space_.nonnull())
    return;

  abs_space_ = orthogonalize("p'","RIBS", bs_aux_, integral(), orthog_method(), lindep_tol(), nlindep_aux_);
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
    ribs_space_ = orthogonalize("p'","RIBS", bs_ri_, integral(), orthog_method(), lindep_tol(), nlindep_ri_);
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
  Ref<GaussianBasisSet> ri_basis = abs + ref()->basis();
  int nlindep_ri = 0;
  if (bs_ri_.nonnull() && ri_basis->equiv(bs_ri_)) {
    construct_orthog_ri_();
    nlindep_ri = nlindep_ri_;
  }
  else {
    Ref<OrbitalSpace> ribs_space = orthogonalize("p+p'","OBS+ABS", ri_basis, integral(), orthog_method(), lindep_tol(), nlindep_ri);
  }

  const int obs_rank = ref()->orbs(Alpha)->rank();
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

  const double tol = lindep_tol();
  if (!ref()->spin_polarized()) {
    Ref<OrbitalSpace> tmp = orthog_comp(ref()->occ_sb(Alpha), ribs_space_, "p'-m", "CABS", tol);
    tmp = orthog_comp(ref()->uocc_sb(Alpha), tmp, "a'", "CABS", tol);
    cabs_space_[Alpha] = tmp;
    cabs_space_[Beta] = cabs_space_[Alpha];
    idxreg->add(make_keyspace_pair(cabs_space_[Alpha]));
  }
  else {
    Ref<OrbitalSpace> tmp_a = orthog_comp(ref()->occ_sb(Alpha), ribs_space_, "P'-M", "CABS (Alpha)", tol);
    Ref<OrbitalSpace> tmp_b = orthog_comp(ref()->occ_sb(Beta), ribs_space_, "p'-m", "CABS (Beta)", tol);
    if (USE_NEW_ORBITALSPACE_KEYS) {
      const std::string key_a = ParsedOrbitalSpaceKey::key(std::string("a'"),Alpha);
      const std::string key_b = ParsedOrbitalSpaceKey::key(std::string("a'"),Beta);
      cabs_space_[Alpha] = orthog_comp(ref()->uocc_sb(Alpha), tmp_a, key_a, "CABS (Alpha)", tol);
      cabs_space_[Beta] = orthog_comp(ref()->uocc_sb(Beta), tmp_b, key_b, "CABS (Beta)", tol);
    }
    else // old orbitalspace key no longer supported
      assert(false);

    idxreg->add(make_keyspace_pair(cabs_space_[Alpha]));
    idxreg->add(make_keyspace_pair(cabs_space_[Beta]));
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
