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

#include <stdexcept>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <util/misc/string.h>
#include <util/misc/scexception.h>
#include <math/scmat/local.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <chemistry/qc/lcao/transform_factory.h>
#include <util/misc/registry.timpl.h>
#include <chemistry/qc/nbody/ref.h>
#ifdef HAVE_MADNESS
# include <util/madness/init.h>
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
    const Ref<RefWavefunction>& ref,
    Ref<OrbitalSpace> ri_space) :
    refwfn_(ref), ref_acc_for_cabs_space_(DBL_MAX),
    ribs_space_(ri_space), ribs_space_given_(ribs_space_.nonnull())
{
  // by default use spin-orbital algorithm
  spinadapted_ = keyval->booleanvalue("spinadapted",KeyValValueboolean(false));
  //refwfn_->set_spinfree(spinadapted_);

  if (ribs_space_.null()) {
    bs_aux_ = require_dynamic_cast<GaussianBasisSet*>(
        keyval->describedclassvalue("aux_basis").pointer(),
        "R12Technology::R12Technology\n"
        );
    if (bs_aux_.pointer() == NULL)
      bs_aux_ = ref->basis();
  }
  else {
    bs_ri_ = bs_aux_ = ribs_space_->basis();
  }

  r12tech_ = new R12Technology(keyval,ref->basis(),ref->uocc_basis(),bs_aux_);
  // Make sure can use the integral factory for R12 calcs
  r12tech_->check_integral_factory(integral());

  // boot up madness
#ifdef HAVE_MADNESS
  MADNESSRuntime::initialize();
#endif
}

R12WavefunctionWorld::R12WavefunctionWorld(StateIn& si) : SavableState(si)
{
  refwfn_ << SavableState::restore_state(si);
  bs_aux_ << SavableState::restore_state(si);
  bs_ri_ << SavableState::restore_state(si);


  si.get(spinadapted_);

  abs_space_ << SavableState::restore_state(si);
  ribs_space_ << SavableState::restore_state(si);

  // boot up madness
#ifdef HAVE_MADNESS
  MADNESSRuntime::initialize();
#endif
}

R12WavefunctionWorld::~R12WavefunctionWorld()
{
  // boot up madness
#ifdef HAVE_MADNESS
  MADNESSRuntime::finalize();
#endif
}

void R12WavefunctionWorld::save_data_state(StateOut& so)
{
  SavableState::save_state(refwfn_.pointer(),so);
  SavableState::save_state(bs_aux_.pointer(),so);
  SavableState::save_state(bs_ri_.pointer(),so);
  so.put(spinadapted_);
  SavableState::save_state(abs_space_.pointer(),so);
  SavableState::save_state(ribs_space_.pointer(),so);
  so.put(ribs_space_given_);
}

void
R12WavefunctionWorld::initialize()
{
  // annotate parameters
  Ref<ParamsRegistry> pregistry = ParamsRegistry::instance();
  Ref<R12Technology::CorrelationFactor> corrfactor = r12tech_->corrfactor();
  Ref<R12Technology::NullCorrelationFactor> nullcorrfactor;  nullcorrfactor << corrfactor;
  if (nullcorrfactor.null()) { // real correlation factor
    const unsigned int ncorrfunctions = corrfactor->nfunctions();
    for (unsigned int f = 0; f < ncorrfunctions; ++f) {
      std::ostringstream oss;
      oss << "[" << f << "]";
      const std::string pkey = oss.str();
      if (pregistry->key_exists(pkey) == false) {
        pregistry->add(pkey,
                       corrfactor->tbintdescr(this->integral(), f)->params());
      }
    }
    for (unsigned int f = 0; f < ncorrfunctions; ++f) {
      for (unsigned int g = 0; g < ncorrfunctions; ++g) {
        std::ostringstream oss;
        oss << "[" << f << "," << g << "]";
        const std::string pkey = oss.str();
        if (pregistry->key_exists(pkey) == false) {
          pregistry->add(
                         pkey,
                         corrfactor->tbintdescr(this->integral(), f, g)->params());
        }
      }
    }
  }

  refwfn_->world()->initialize_ao_spaces();
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
  obs_eq_vbs_ = basis_vir().null() || basis()->equiv(basis_vir());
  if (ribs_space_given_ == false) {
    bs_ri_ = 0;
    ribs_space_ = 0;
  }
  cabs_space_[Alpha] = cabs_space_[Beta] = 0;
}

void
R12WavefunctionWorld::obsolete() {
  refwfn_->obsolete();
  abs_space_ = 0;
  cabs_space_[Alpha] = 0;
  cabs_space_[Beta] = 0;

  if (ribs_space_given_ == false) {
    nlindep_ri_ = nlindep_aux_ = -1;
    ribs_space_ = 0;
    bs_ri_ = 0;
  }

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

bool
R12WavefunctionWorld::obs_eq_ribs() const {
  if (bs_ri_.null()) ribs_space();
  return basis()->equiv(bs_ri_);
}

const Ref<OrbitalSpace>&
R12WavefunctionWorld::cabs_space(const SpinCase1& S) const
{
  if (ref_acc_for_cabs_space_ > refwfn()->desired_value_accuracy()) { // recompute if accuracy of reference has increased
    cabs_space_[Alpha] = 0;
    cabs_space_[Beta] = 0;
  }
  if (cabs_space_[S].null()) { // compute if needed
    R12WavefunctionWorld* this_nonconst_ptr = const_cast<R12WavefunctionWorld*>(this);
    this_nonconst_ptr->construct_cabs_();
  }
  return cabs_space_[S];
}

bool
R12WavefunctionWorld::sdref() const {

#define ALWAYS_USE_GENREF_ALGORITHM 0
#if !ALWAYS_USE_GENREF_ALGORITHM
  return refwfn()->sdref();
#endif
  return false;
}

void
R12WavefunctionWorld::print(std::ostream& o) const {

  o << indent << "R12WavefunctionWorld:" << endl;
  o << incindent;

  this->world()->print(o);
  this->refwfn()->print(o);

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
