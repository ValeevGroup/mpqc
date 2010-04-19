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
#include <chemistry/qc/mbptr12/orbitalspace_utils.h>

using namespace std;
using namespace sc;


/*--------------------------------
  MBPT2_R12
 --------------------------------*/

static ClassDesc MBPT2_R12_cd(
  typeid(MBPT2_R12),"MBPT2_R12",9,"public MBPT2",
  0, create<MBPT2_R12>, create<MBPT2_R12>);

MBPT2_R12::MBPT2_R12(StateIn& s):
  MBPT2(s)
{
#if 0
  if (s.version(::class_desc<MBPT2_R12>()) < 9)
    throw InputError("Cannot use MBPT2_R12 class prior to version 9",__FILE__,__LINE__);
#endif

  r12eval_ << SavableState::restore_state(s);
  r12world_ << SavableState::restore_state(s);
  r12a_energy_ << SavableState::restore_state(s);
  r12ap_energy_ << SavableState::restore_state(s);
  r12app_energy_ << SavableState::restore_state(s);
  r12b_energy_ << SavableState::restore_state(s);
  r12c_energy_ << SavableState::restore_state(s);

  if((r12world()->r12tech()->ansatz()->orbital_product_GG()==R12Technology::OrbProdGG_pq) ||
     (r12world()->r12tech()->ansatz()->orbital_product_gg()==R12Technology::OrbProdgg_pq)) {
    throw InputError("MBPT2_R12::MBPT2_R12 -- pq Ansatz not allowed",__FILE__,__LINE__);
  }

  int spinadapted; s.get(spinadapted); spinadapted_ = (bool)spinadapted;
  int unv; s.get(unv); new_energy_ = (bool)unv;
  s.get(mp2_corr_energy_);
  s.get(cabs_singles_);
  s.get(cabs_singles_energy_);

  twopdm_grid_ << SavableState::restore_state(s);
  s.get(plot_pair_function_[0]);
  s.get(plot_pair_function_[1]);
}

MBPT2_R12::MBPT2_R12(const Ref<KeyVal>& keyval):
  MBPT2(keyval)
{
  // MP2-R12 calculations can use virtual orbitals expanded in a separate basis set
  Ref<GaussianBasisSet> bs_vir = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("vir_basis").pointer(),
      "MBPT2_R12::MBPT2_R12\n"
      );
  Ref<OrbitalSpace> vbs;
  if (bs_vir.nonnull()) {
    int nlindep_vir;
    vbs = orthogonalize("e(sym)", "VBS", bs_vir,
                        integral(),
                        this->orthog_method(), this->lindep_tol(),
                        nlindep_vir);
  }

  Ref<WavefunctionWorld> world = new WavefunctionWorld(keyval, this);
  Ref<RefWavefunction> refinfo = new SD_RefWavefunction(world, ref(), false,
                                                              nfzcore(), nfzvirt(),
                                                              vbs);
  r12world_ = new R12WavefunctionWorld(keyval, refinfo);
  Ref<R12Technology> r12tech = r12world_->r12tech();

  cabs_singles_ = keyval->booleanvalue("cabs_singles",KeyValValueboolean((int)false));

  const bool diag = r12tech->ansatz()->diag();
  const bool optimized_amplitudes = r12tech->ansatz()->amplitudes() == R12Technology::GeminalAmplitudeAnsatz_fullopt;

  const bool default_new_energy = diag ? true : false;
  new_energy_ = keyval->booleanvalue("new_energy",KeyValValueboolean(default_new_energy));
  if((new_energy_==true) && (diag==false)) {
    ExEnv::out0() << indent << "Warning: The non diagonal ansatz is safer to be computed with the old version" << endl
                  << indent << "because the old version performs many security checks." << endl;
  }
  if (optimized_amplitudes==false) { // fixed amplitudes can only be used with the new energy object
    new_energy_ = true;
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
  cabs_singles_energy_ = 0.0;

  this->set_desired_value_accuracy(desired_value_accuracy());
}

MBPT2_R12::~MBPT2_R12()
{
}

void
MBPT2_R12::save_data_state(StateOut& s)
{
  MBPT2::save_data_state(s);
  SavableState::save_state(r12eval_.pointer(),s);
  SavableState::save_state(r12world_.pointer(),s);
  SavableState::save_state(r12a_energy_.pointer(),s);
  SavableState::save_state(r12ap_energy_.pointer(),s);
  SavableState::save_state(r12app_energy_.pointer(),s);
  SavableState::save_state(r12b_energy_.pointer(),s);
  SavableState::save_state(r12c_energy_.pointer(),s);

  s.put((int)spinadapted_);
  s.put((int)new_energy_);
  s.put(mp2_corr_energy_);
  s.put(cabs_singles_);
  s.put(cabs_singles_energy_);

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
  o << indent << "Use new MP2R12Energy: " << (new_energy_ ? "true" : "false") << endl;
  o << indent << "Include CABS singles? : " << (cabs_singles_ ? "true" : "false") << endl;
  if (cabs_singles_) {
    o << indent << "  E(CABS singles) = " << scprintf("%25.15lf", cabs_singles_energy_)
                                          << endl;
  }

  r12world()->print(o);

  MBPT2::print(o);
  o << decindent;
}

void
MBPT2_R12::set_desired_value_accuracy(double acc)
{
  Function::set_desired_value_accuracy(acc);
  if (ref()->desired_value_accuracy_set_to_default()) {
    // reference should be computed to higher accuracy
    const double ref_acc = acc * ref_to_mp2r12_acc();
    ref()->set_desired_value_accuracy(ref_acc);
  }
}

//////////////////////////////////////////////////////////////////////////////

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
  cabs_singles_energy_ = 0.0;
  r12world_->world()->obsolete();
  r12world_->obsolete();
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
MBPT2_R12::corrfactor(const Ref<R12Technology::CorrelationFactor>& cf)
{
    if (!r12world()->r12tech()->corrfactor()->equiv(cf)) {
      r12world()->r12tech()->corrfactor(cf);
      obsolete();
  }
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2_R12::corr_energy()
{
  energy();
  return energy() - ref_energy() - cabs_singles_energy();
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2_R12::r12_corr_energy()
{
  energy();
  return energy() - mp2_corr_energy_ - ref_energy() - cabs_singles_energy();
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2_R12::cabs_singles_energy()
{
  energy();
  return cabs_singles_energy_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
