//
// r12wfnworld.cc
//
// Copyright (C) 2003 Edward Valeev
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <stdexcept>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <util/misc/string.h>
#include <util/class/scexception.h>
#include <math/scmat/local.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/orbitalspace.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/registry.timpl.h>
#include <chemistry/qc/mbptr12/ref.h>
#if HAVE_PSIMPQCIFACE
# include <chemistry/qc/psi/psiref.h>
#endif

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*---------------
  R12WavefunctionWorld
 ---------------*/
static ClassDesc R12WavefunctionWorld_cd(
  typeid(R12WavefunctionWorld),"R12WavefunctionWorld",11,"virtual public SavableState",
  0, 0, create<R12WavefunctionWorld>);

R12WavefunctionWorld::R12WavefunctionWorld(
    const Ref<KeyVal>& keyval,
    const Ref<RefWavefunction>& ref) :
    ref_(ref)
{
  // by default use spin-adapted algorithms for closed-shell
  spinadapted_ = keyval->booleanvalue("spinadapted",KeyValValueboolean(true));
  if (ref->spin_polarized()) spinadapted_ = false;

  bs_aux_ = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("aux_basis").pointer(),
      "R12Technology::R12Technology\n"
      );
  if (bs_aux_.pointer() == NULL)
      bs_aux_ = ref->basis();

  r12tech_ = new R12Technology(keyval,ref->basis(),ref->uocc_basis(),bs_aux_);
  // Make sure can use the integral factory for R12 calcs
  r12tech_->check_integral_factory(integral());
}

R12WavefunctionWorld::R12WavefunctionWorld(StateIn& si) : SavableState(si)
{
  ref_ << SavableState::restore_state(si);
  bs_aux_ << SavableState::restore_state(si);
  bs_ri_ << SavableState::restore_state(si);

  si.get(spinadapted_);

  abs_space_ << SavableState::restore_state(si);
  ribs_space_ << SavableState::restore_state(si);
}

R12WavefunctionWorld::~R12WavefunctionWorld()
{
}

void R12WavefunctionWorld::save_data_state(StateOut& so)
{
  SavableState::save_state(ref_.pointer(),so);
  SavableState::save_state(bs_aux_.pointer(),so);
  SavableState::save_state(bs_ri_.pointer(),so);
  so.put(spinadapted_);
  SavableState::save_state(abs_space_.pointer(),so);
  SavableState::save_state(ribs_space_.pointer(),so);
}

void
R12WavefunctionWorld::initialize()
{
  ref_->world()->initialize_ao_spaces();
  // provide hints to the factory about the likely use of transforms
  {
    // 1) if stdapprox is A' or A'' most transforms will never be reused
    // 2) otherwise use persistent transforms
    if (r12tech()->stdapprox() == R12Technology::StdApprox_Ap ||
        r12tech()->stdapprox() == R12Technology::StdApprox_App)
      world()->tfactory()->hints().data_persistent(false);
    else
      world()->tfactory()->hints().data_persistent(true);
  }

  nlindep_ri_ = nlindep_aux_ = -1;
  obs_eq_vbs_ = basis()->equiv(basis_vir());
  bs_ri_ = 0;
  ribs_space_ = 0;
  cabs_space_[Alpha] = cabs_space_[Beta] = 0;
}

void
R12WavefunctionWorld::obsolete() {
  ref_->obsolete();
  abs_space_ = 0;
  ribs_space_ = 0;
  cabs_space_[Alpha] = 0;
  cabs_space_[Beta] = 0;
  nlindep_ri_ = nlindep_aux_ = -1;
  bs_ri_ = 0;

//  this->initialize();   // can't initialize because there is no guarantee we are ready to compute again
}

const Ref<OrbitalSpace>&
R12WavefunctionWorld::ribs_space() const
{
  if (ribs_space_.null()) {
    R12WavefunctionWorld* this_nonconst_ptr = const_cast<R12WavefunctionWorld*>(this);
    this_nonconst_ptr->construct_ri_basis_(r12tech()->safety_check());
  }
  return ribs_space_;
}

const Ref<OrbitalSpace>&
R12WavefunctionWorld::cabs_space(const SpinCase1& S) const
{
  if (r12tech()->abs_method() == R12Technology::ABS_CABS ||
      r12tech()->abs_method() == R12Technology::ABS_CABSPlus) {
    if (ref_acc_for_cabs_space_ > ref()->desired_value_accuracy()) { // recompute if accuracy of reference has increased
      cabs_space_[Alpha] = 0;
      cabs_space_[Beta] = 0;
    }
    if (cabs_space_[S].null()) { // compute if needed
      R12WavefunctionWorld* this_nonconst_ptr = const_cast<R12WavefunctionWorld*>(this);
      this_nonconst_ptr->construct_cabs_();
    }
    return cabs_space_[S];
  }
  else
    throw ProgrammingError("CABS space requested by abs_method set to ABS/ABS+",__FILE__,__LINE__);
}

bool
R12WavefunctionWorld::sdref() const {

#define ALWAYS_USE_GENREF_ALGORITHM 0
#if !ALWAYS_USE_GENREF_ALGORITHM
  // only references based on OneBodyWavefunction are detected as single-determinant references!
  {
    Ref<SD_RefWavefunction> sd; sd << ref();
    if (sd.nonnull()) return true;
  }
#if HAVE_PSIMPQCIFACE
  {
    Ref<PsiSCF_RefWavefunction> sd; sd << ref();
    if (sd.nonnull()) return true;
  }
#endif
#endif
  return false;
}

void
R12WavefunctionWorld::print(std::ostream& o) const {

  o << indent << "R12WavefunctionWorld:" << endl;
  o << incindent;

  this->world()->print(o);

  if (!bs_aux_->equiv(basis())) {
      o << indent << "Auxiliary Basis Set (ABS):" << endl;
      o << incindent; bs_aux_->print(o); o << decindent << endl;
  }
  r12tech()->print(o);
  o << indent << "Spin-adapted algorithm: " << (spinadapted_ ? "true" : "false") << endl;
  o << decindent << endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
