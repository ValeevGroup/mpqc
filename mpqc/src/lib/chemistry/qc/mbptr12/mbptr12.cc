//
// mbptr12.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/mbptr12/mbptr12.h>

using namespace std;
using namespace sc;

/*--------------------------------
  MBPT2_R12
 --------------------------------*/

static ClassDesc MBPT2_R12_cd(
  typeid(MBPT2_R12),"MBPT2_R12",1,"public MBPT2",
  0, create<MBPT2_R12>, create<MBPT2_R12>);

MBPT2_R12::MBPT2_R12(StateIn& s):
  MBPT2(s)
{
  r12eval_ << SavableState::restore_state(s);
  r12a_energy_ << SavableState::restore_state(s);
  r12ap_energy_ << SavableState::restore_state(s);
  r12b_energy_ << SavableState::restore_state(s);
  aux_basis_ << SavableState::restore_state(s);
  int stdapprox; s.get(stdapprox); stdapprox_ = (LinearR12::StandardApproximation) stdapprox;
  int spinadapted; s.get(spinadapted); spinadapted_ = (bool)spinadapted;
  s.getstring(r12ints_file_);
  s.get(mp2_corr_energy_);
  s.get(r12_corr_energy_);
}

MBPT2_R12::MBPT2_R12(const Ref<KeyVal>& keyval):
  MBPT2(keyval)
{
  // Make sure can use the integral factory for linear R12
  check_integral_factory_();
  
  // Check is this is a closed-shell reference
  CLHF* clhfref = dynamic_cast<CLHF*>(reference_.pointer());
  if (!clhfref) {
    ExEnv::err0() << "MBPT2_R12::MBPT2_R12: reference wavefunction is of non-CLHF type" << endl;
    abort();
  }

  aux_basis_ = require_dynamic_cast<GaussianBasisSet*>(
    keyval->describedclassvalue("aux_basis").pointer(),
    "MBPT2_R12::MBPT2_R12\n"
    );
  if (aux_basis_.pointer() == NULL)
    aux_basis_ = basis();

  // Default method is MBPT2-R12/A
  stdapprox_ = LinearR12::StdApprox_A;
  char* sa_string = 0;
  sa_string = keyval->pcharvalue("stdapprox");
  if ( !strcmp(sa_string,"A") ||
       !strcmp(sa_string,"a") ) {
    stdapprox_ = LinearR12::StdApprox_A;
  }
  else if ( !strcmp(sa_string,"Ap") ||
	    !strcmp(sa_string,"ap") ||
	    !strcmp(sa_string,"A'") ||
	    !strcmp(sa_string,"a'") ) {
    stdapprox_ = LinearR12::StdApprox_Ap;
  }
  else if ( !strcmp(sa_string,"B") ||
	    !strcmp(sa_string,"b") )
    throw std::runtime_error("MBPT2_R12: MP2-R12/B energy is not implemented yet");


  // Default is to use spin-adapted algorithm
  spinadapted_ = true;
  spinadapted_ = keyval->booleanvalue("spinadapted");

  // Get the filename to store the integrals
  r12ints_file_ = 0;
  r12ints_file_ = keyval->pcharvalue("r12ints_file");
  if (!r12ints_file_) {
    r12ints_file_ = strdup("/tmp/r12ints.dat");
  }

  r12eval_ = 0;
  r12a_energy_ = 0;
  r12ap_energy_ = 0;
  r12b_energy_ = 0;
  mp2_corr_energy_ = 0.0;
  r12_corr_energy_ = 0.0;
}

MBPT2_R12::~MBPT2_R12()
{
  r12eval_ = 0;
  r12a_energy_ = 0;
  r12ap_energy_ = 0;
  r12b_energy_ = 0;
  delete[] r12ints_file_;
}

void
MBPT2_R12::save_data_state(StateOut& s)
{
  MBPT2::save_data_state(s);
  SavableState::save_state(r12eval_.pointer(),s);
  SavableState::save_state(r12a_energy_.pointer(),s);
  SavableState::save_state(r12ap_energy_.pointer(),s);
  SavableState::save_state(r12b_energy_.pointer(),s);
  SavableState::save_state(aux_basis_.pointer(),s);
  s.put((int)stdapprox_);
  s.put((int)spinadapted_);
  s.putstring(r12ints_file_);

  s.put(mp2_corr_energy_);
  s.put(r12_corr_energy_);
}

void
MBPT2_R12::print(ostream&o) const
{
  o << indent << "MBPT2_R12:" << endl;
  o << incindent;
  switch (stdapprox_) {
    case LinearR12::StdApprox_A :
      o << indent << "Standard Approximation: A" << endl;
    break;
    case LinearR12::StdApprox_Ap :
      o << indent << "Standard Approximation: A'" << endl;
    break;
    case LinearR12::StdApprox_B :
      o << indent << "Standard Approximation: B" << endl;
    break;
  }
  o << indent << "Spin-adapted algorithm: " << (spinadapted_ ? "true" : "false") << endl;
  o << indent << "Transformed Integrals file: " << r12ints_file_ << endl << endl;
  o << indent << "Auxiliary Basis:" << endl;
  o << incindent; aux_basis_->print(o); o << decindent << endl;
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
  init_variables_();
  reference_->set_desired_value_accuracy(desired_value_accuracy()
                                         / ref_to_mp2r12_acc_);

  compute_energy_a_();
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::obsolete()
{
  r12eval_ = 0;
  r12a_energy_ = 0;
  r12ap_energy_ = 0;
  r12b_energy_ = 0;
  mp2_corr_energy_ = 0.0;
  r12_corr_energy_ = 0.0;
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

/////////////////////////////////////////////////////////////////////////////

Ref<GaussianBasisSet>
MBPT2_R12::aux_basis() const
{
  return aux_basis_;
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

/////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::init_variables_()
{
  MBPT2::init_variables();
}

/////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::check_integral_factory_()
{
  // Only IntegralCints can be used at the moment
  IntegralCints* r12intf = dynamic_cast<IntegralCints*>(integral().pointer());
  if (!r12intf) {
    ostringstream errmsg;
    errmsg << "Integral factory " << (integral().null() ? "null" : integral()->class_name())
           << " cannot be used in MBPT2_R12 class - try IntegralCints instead" << ends;
    throw runtime_error(errmsg.str());
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
