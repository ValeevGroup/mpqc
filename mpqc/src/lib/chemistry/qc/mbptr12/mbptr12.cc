//
// mbptr12.cc
//
// Copyright (C) 2001 Edward Valeev
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
#include <util/misc/string.h>
#include <util/class/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/hsosscf.h>

using namespace std;
using namespace sc;
using namespace sc::LinearR12;

/*--------------------------------
  MBPT2_R12
 --------------------------------*/

static ClassDesc MBPT2_R12_cd(
  typeid(MBPT2_R12),"MBPT2_R12",9,"public MBPT2",
  0, create<MBPT2_R12>, create<MBPT2_R12>);

MBPT2_R12::MBPT2_R12(StateIn& s):
  MBPT2(s)
{
  if (s.version(::class_desc<MBPT2_R12>()) < 9)
      throw InputError("Cannot use MBPT2_R12 class prior to version 9",__FILE__,__LINE__);

  r12eval_ << SavableState::restore_state(s);
  r12evalinfo_ << SavableState::restore_state(s);
  r12a_energy_ << SavableState::restore_state(s);
  r12ap_energy_ << SavableState::restore_state(s);
  r12app_energy_ << SavableState::restore_state(s);
  r12b_energy_ << SavableState::restore_state(s);
  r12c_energy_ << SavableState::restore_state(s);

  int spinadapted; s.get(spinadapted); spinadapted_ = (bool)spinadapted;
  int incmp1; s.get(incmp1); include_mp1_ = (bool)incmp1;
  int hylleraas; s.get(hylleraas); hylleraas_ = (bool)hylleraas;
  int unv; s.get(unv); new_energy_ = (bool)unv;
  s.get(mp2_corr_energy_);

  twopdm_grid_ << SavableState::restore_state(s);
  s.get(plot_pair_function_[0]);
  s.get(plot_pair_function_[1]);
}

MBPT2_R12::MBPT2_R12(const Ref<KeyVal>& keyval):
  MBPT2(keyval)
{
  // Verify that this is a closed-shell or high-spin open-shell system
  CLSCF* clscfref = dynamic_cast<CLSCF*>(ref().pointer());
  HSOSSCF* roscfref = dynamic_cast<HSOSSCF*>(ref().pointer());
  if (roscfref != 0) {
    ExEnv::out0() << indent << "WARNING: ROHF-based MBPT2-R12 method not completely tested yet" << endl;
  }
  const bool closedshell = (clscfref != 0);

  spinadapted_ = false;
  if (closedshell)
    // Default is to use spin-adapted algorithm
    spinadapted_ = keyval->booleanvalue("spinadapted",KeyValValueboolean((int)true));

  r12evalinfo_ = new R12IntEvalInfo(keyval,this,ref(),nfzcore(), nfzvirt(), spinadapted_, true);

  const bool diag = r12evalinfo()->ansatz()->diag();
  const bool fixedcoeff = r12evalinfo()->ansatz()->fixedcoeff();

  // Default is to not compute MP1 energy
  include_mp1_ = false;
  // only check if VBS != OBS
  if (!r12evalinfo()->basis_vir()->equiv(r12evalinfo()->basis()))
      include_mp1_ = keyval->booleanvalue("include_mp1",KeyValValueboolean((int)false));
  
  new_energy_ = diag ? true : false;
  new_energy_ = keyval->booleanvalue("new_energy",KeyValValueboolean((int)false));
  if((new_energy_==true) && (diag==false)) {
    ExEnv::out0() << indent << "Warning: The non diagonal ansatz is safer to be computed with the old version" << endl
                  << indent << "because the old version performs many security checks." << endl;
  }
  if((new_energy_==false) and (fixedcoeff==true)) {
    throw InputError("MBPT2_R12::MBPT2_R12 -- fixed coefficients are not implemented in the old version. Set new_energy to true in your input.",__FILE__,__LINE__);
  }
  
  const bool hyll_default = (fixedcoeff==true) ? true : false;
  hylleraas_ = keyval->booleanvalue("hylleraas",KeyValValueboolean((int)hyll_default));
  
  if((fixedcoeff==false) && (hylleraas_==true)){
    throw InputError("MBPT2_R12::MBPT2_R12 -- Hylleraas functional needn't be calculated if coefficients are optimized (gives the same result as non Hylleraas).",__FILE__,__LINE__);
  }

  twopdm_grid_ = require_dynamic_cast<TwoBodyGrid*>(
                   keyval->describedclassvalue("twopdm_grid").pointer(),
                   "MBPT2_R12::MBPT2_R12\n"
                 );
  if (twopdm_grid_.nonnull()) {
    plot_pair_function_[0] = static_cast<unsigned int>(keyval->intvalue("plot_pair_function", 0,
									KeyValValueint(0)
									));
    plot_pair_function_[1] = static_cast<unsigned int>(keyval->intvalue("plot_pair_function", 1,
									KeyValValueint(0)
									));
  }

  r12eval_ = 0;
  r12a_energy_ = 0;
  r12ap_energy_ = 0;
  r12app_energy_ = 0;
  r12b_energy_ = 0;
  r12c_energy_ = 0;
  mp2_corr_energy_ = 0.0;
  
  set_desired_value_accuracy(desired_value_accuracy());
}

MBPT2_R12::~MBPT2_R12()
{
}

void
MBPT2_R12::save_data_state(StateOut& s)
{
  MBPT2::save_data_state(s);
  SavableState::save_state(r12eval_.pointer(),s);
  SavableState::save_state(r12evalinfo_.pointer(),s);
  SavableState::save_state(r12a_energy_.pointer(),s);
  SavableState::save_state(r12ap_energy_.pointer(),s);
  SavableState::save_state(r12app_energy_.pointer(),s);
  SavableState::save_state(r12b_energy_.pointer(),s);
  SavableState::save_state(r12c_energy_.pointer(),s);

  s.put((int)spinadapted_);
  s.put((int)include_mp1_);
  s.put((int)hylleraas_);
  s.put((int)new_energy_);
  s.put(mp2_corr_energy_);

  SavableState::save_state(twopdm_grid_.pointer(),s);
  s.put((int)plot_pair_function_[0]);
  s.put((int)plot_pair_function_[1]);
}

void
MBPT2_R12::print(ostream&o) const
{
  o << indent << "MBPT2_R12:" << endl;
  o << incindent;

  o << indent << "Spin-adapted algorithm: " << (spinadapted_ ? "true" : "false") << endl;
  if (!r12evalinfo()->basis_vir()->equiv(r12evalinfo()->basis()))
    o << indent << "Compute MP1 energy: " << (include_mp1_ ? "true" : "false") << endl;
  o << indent << "Use new MP2R12Energy: " << (new_energy_ ? "true" : "false") << endl;
  o << indent << "Compute Hylleraas functional: " << (hylleraas_ ? "true" : "false") << endl;
  

  r12evalinfo()->r12tech()->print(o);

  MBPT2::print(o);
  o << decindent;
}

RefSymmSCMatrix
MBPT2_R12::density()
{
  ExEnv::out0() << "MBPT2_R12::density() is not implemented" << endl;
  abort();
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::compute()
{
  ExEnv::out0() << "switched off integral check in mbptr12.cc/compute" << endl;

/*  if (std::string(reference_->integral()->class_name()) != integral()->class_name()) {
      FeatureNotImplemented ex(
          "cannot use a reference with a different Integral specialization",
          __FILE__, __LINE__, class_desc());
      try {
          ex.elaborate()
              << "reference uses " << reference_->integral()->class_name()
              << " but this object uses " << integral()->class_name()
              << std::endl;
        }
      catch (...) {}
      throw ex;
    }
*/
  compute_energy_();
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::obsolete()
{
  r12eval_ = 0;
  r12a_energy_ = 0;
  r12ap_energy_ = 0;
  r12app_energy_ = 0;
  r12b_energy_ = 0;
  r12c_energy_ = 0;
  mp2_corr_energy_ = 0.0;
  MBPT2::obsolete();
}

//////////////////////////////////////////////////////////////////////////////

int
MBPT2_R12::gradient_implemented() const
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

int
MBPT2_R12::value_implemented() const
{
  return 1;
}

////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::corrfactor(const Ref<LinearR12::CorrelationFactor>& cf)
{
    if (!r12evalinfo()->r12tech()->corrfactor()->equiv(cf)) {
      r12evalinfo()->r12tech()->corrfactor(cf);
      obsolete();
  }
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2_R12::corr_energy()
{
  energy();
  return energy() - ref_energy();
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2_R12::r12_corr_energy()
{
  energy();
  return energy() - mp2_corr_energy_ - ref_energy();
}

////////////////////////////////////////////////////////////////////////////

bool
MBPT2_R12::hylleraas() const
{
  return hylleraas_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
