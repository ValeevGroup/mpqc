//
// mp2extrap.cc
//
// Copyright (C) 1998 Limit Point Systems, Inc.
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

#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <chemistry/qc/mbpt/mbpt.h>
#include <chemistry/qc/mbpt/mp2extrap.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////
// MP2BasisExtrap

static ClassDesc MP2BasisExtrap_cd(
  typeid(MP2BasisExtrap),"MP2BasisExtrap",1,"public SumMolecularEnergy",
  0, create<MP2BasisExtrap>, create<MP2BasisExtrap>);

MP2BasisExtrap::MP2BasisExtrap(const Ref<KeyVal> &keyval):
  SumMolecularEnergy(keyval)
{
  if (n_ != 3) {
      ExEnv::out0() << "ERROR: MP2BasisExtrap: require exactly 3 energies"
           << endl;
      abort();
    }

  // the first row of the inverse of a gives the coefficients
  //a = [ 1, -1/81,  -1/243;
  //      1, -1/256, -1/1024;
  //      1, -1/625, -1/3125; ]
  if (!keyval->exists("coef",0)
      &&!keyval->exists("coef",1)
      &&!keyval->exists("coef",2)) {
    coef_[0] =  0.184090909090909;
    coef_[1] = -1.551515151515153;
    coef_[2] =  2.367424242424244;
    }

  MBPT2 *mbpt[3];
  if ((mbpt[0] = dynamic_cast<MBPT2*>(mole_[0].pointer())) == 0
      ||(mbpt[1] = dynamic_cast<MBPT2*>(mole_[1].pointer())) == 0
      ||(mbpt[2] = dynamic_cast<MBPT2*>(mole_[2].pointer())) == 0) {
      ExEnv::out0() << "ERROR: MP2BasisExtrap: need MBPT2 objects"
           << endl;
      abort();
    }
  if (strcmp(mbpt[0]->basis()->name(),"cc-pVDZ")
      ||strcmp(mbpt[1]->basis()->name(),"cc-pVTZ")
      ||strcmp(mbpt[2]->basis()->name(),"cc-pVQZ")) {
      ExEnv::out0() << "WARNING: MP2BasisExtrap:" << endl
           << "  given basis sets: "
           << mbpt[0]->basis()->name() << ", "
           << mbpt[1]->basis()->name() << ", "
           << mbpt[2]->basis()->name() << endl
           << "  but prefer cc-pVDZ, cc-pVTZ, cc-pVQZ" << endl;
    }
}

MP2BasisExtrap::MP2BasisExtrap(StateIn&s):
  SumMolecularEnergy(s)
{
}

void
MP2BasisExtrap::save_data_state(StateOut&s)
{
  SumMolecularEnergy::save_data_state(s);
}

MP2BasisExtrap::~MP2BasisExtrap()
{
}

void
MP2BasisExtrap::compute()
{
  int i;

  MBPT2 *mbpt2[3];
  mbpt2[0] = dynamic_cast<MBPT2*>(mole_[0].pointer());
  mbpt2[1] = dynamic_cast<MBPT2*>(mole_[1].pointer());
  mbpt2[2] = dynamic_cast<MBPT2*>(mole_[2].pointer());

  int *old_do_value = new int[n_];
  int *old_do_gradient = new int[n_];
  int *old_do_hessian = new int[n_];

  for (i=0; i<n_; i++)
      old_do_value[i] = mole_[i]->do_value(value_.compute());
  for (i=0; i<n_; i++)
      old_do_gradient[i]=mole_[i]->do_gradient(gradient_.compute());
  for (i=0; i<n_; i++)
      old_do_hessian[i] = mole_[i]->do_hessian(hessian_.compute());

  ExEnv::out0() << indent
       << "MP2BasisExtrap: compute" << endl;

  ExEnv::out0() << incindent;

  if (value_needed()) {
      double val = 0.0;
      double accuracy = 0.0;
      for (i=0; i<n_; i++) {
          val += coef_[i] * mbpt2[i]->corr_energy();
          if (mbpt2[i]->actual_value_accuracy() > accuracy)
              accuracy = mbpt2[i]->actual_value_accuracy();
        }
      val += mbpt2[2]->ref_energy();
      ExEnv::out0() << endl << indent
           << "MP2BasisExtrap =" << endl;
      for (i=0; i<n_; i++) {
          ExEnv::out0() << indent
               << scprintf("  %c % 16.12f * % 16.12f",
                           (i==0?' ':'+'),
                           coef_[i], mbpt2[i]->corr_energy())
               << endl;
        }
      ExEnv::out0() << indent
           << scprintf("  + % 16.12f",
                       mbpt2[2]->ref_energy())
           << endl;
      ExEnv::out0() << indent
           << scprintf("  = % 16.12f", val) << endl;
      set_energy(val);
      set_actual_value_accuracy(accuracy);
    }
  if (gradient_needed()) {
      RefSCVector gradientvec = matrixkit()->vector(moldim());
      gradientvec->assign(0.0);
      double accuracy = 0.0;
      for (i=0; i<n_; i++) {
          gradientvec.accumulate(coef_[i] * mbpt2[i]->corr_energy_gradient());
          if (mbpt2[i]->actual_gradient_accuracy() > accuracy)
              accuracy = mbpt2[i]->actual_gradient_accuracy();
        }
      gradientvec.accumulate(mbpt2[2]->ref_energy_gradient());
      print_natom_3(mbpt2[2]->gradient(),
                    "Total MP2 Gradient with Largest Basis Set");
      print_natom_3(gradientvec,"Total Extrapolated MP2 Gradient");
      set_gradient(gradientvec);
      set_actual_gradient_accuracy(accuracy);
    }
  if (hessian_needed()) {
    ExEnv::out0()
                 << "ERROR: MP2BasisExtrap: cannot do hessian" << endl;
    abort();
    }

  ExEnv::out0() << decindent;

  for (i=0; i<n_; i++) mole_[i]->do_value(old_do_value[i]);
  for (i=0; i<n_; i++) mole_[i]->do_gradient(old_do_gradient[i]);
  for (i=0; i<n_; i++) mole_[i]->do_hessian(old_do_hessian[i]);

  delete[] old_do_value;
  delete[] old_do_gradient;
  delete[] old_do_hessian;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
