//
// hfacm.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#include <math.h>

#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <chemistry/qc/dft/functional.h>
#include <chemistry/qc/dft/integrator.h>
#include <chemistry/qc/dft/hfacm.h>

///////////////////////////////////////////////////////////////////////////
// HFACM

#define CLASSNAME HFACM
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public Wavefunction
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
HFACM::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Wavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

HFACM::HFACM(StateIn& s):
  Wavefunction(s)
  maybe_SavableState(s)
{
  abort();
}

HFACM::HFACM(const RefKeyVal& keyval):
  Wavefunction(keyval)
{
  integrator_ = keyval->describedclassvalue("integrator");
  if (integrator_.null()) integrator_ = new Murray93Integrator();

  functional_ = keyval->describedclassvalue("functional");
  if (functional_.null()) functional_ = new XalphaFunctional();

  a0_ = keyval->doublevalue("a0");
  if (keyval->error() != KeyVal::OK) a0_ = 0.0;

  scf_ = keyval->describedclassvalue("hf");
  if (scf_.null()) {
      cerr << "HFACM: requires \"hf\" keyword of type SCF" << endl;
      abort();
    }
}

HFACM::~HFACM()
{
}

void
HFACM::save_data_state(StateOut& s)
{
  cout << "HFACM: cannot save state" << endl;
  abort();
}

void
HFACM::compute()
{
  if (hessian_needed())
      set_desired_gradient_accuracy(desired_hessian_accuracy()/100.0);

  if (gradient_needed())
      set_desired_value_accuracy(desired_gradient_accuracy()/100.0);

  if (value_needed()) {
      cout << node0 << endl << indent
           << scprintf("HFACM::compute: energy accuracy = %10.7e\n",
                       desired_value_accuracy())
           << endl;
      double e = compute_energy();
      cout << node0 << endl << indent
           << scprintf("total hfacm energy = %20.15f", e)
           << endl;

      set_energy(e);
      set_actual_value_accuracy(desired_value_accuracy());
    }
  if (gradient_needed()) {
      cerr << "HFACM: gradient not implemented" << endl;
      abort();
    }
  if (hessian_needed()) {
      cerr << "HFACM: hessian not implemented" << endl;
      abort();
    }
}

double
HFACM::compute_energy()
{
  cout << incindent;
  scf_->set_desired_value_accuracy(desired_value_accuracy());
  double e1 = scf_->one_body_energy();
  double enn = molecule()->nuclear_repulsion_energy();
  double ec,ex;
  scf_->two_body_energy(ec,ex);
  double escf = scf_->energy();
  cout << decindent;

  integrator_->set_wavefunction(this);
  tim_enter("integrate");
  integrator_->integrate(functional_);
  tim_exit("integrate");

  cout << node0 << indent
       << scprintf("E1           = % 14.10f", e1) << endl
       << scprintf("Enn          = % 14.10f", enn) << endl
       << scprintf("Ec           = % 14.10f", ec) << endl
       << scprintf("Ex           = % 14.10f", ex) << endl
       << scprintf("E(rho)       = % 14.10f", integrator_->value()) << endl
       << scprintf("E1+Enn+Ec+Ex = % 14.10f", e1+enn+ec+ex) << endl
       << scprintf("Escf         = % 14.10f", escf) << endl;

  return e1 + enn + ec + a0_*ex + integrator_->value();
}

RefSymmSCMatrix
HFACM::density()
{
  return scf_->density();
}

int
HFACM::spin_polarized()
{
  return scf_->spin_polarized();
}

RefSymmSCMatrix
HFACM::alpha_density()
{
  return scf_->alpha_density();
}

RefSymmSCMatrix
HFACM::beta_density()
{
  return scf_->beta_density();
}

int
HFACM::value_implemented()
{
  return 1;
}

int
HFACM::gradient_implemented()
{
  return 0;
}

int
HFACM::nelectron()
{
  return scf_->nelectron();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
