//
// hfacm.h --- definition of the HFACM DFT class
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

#ifndef _chemistry_qc_dft_hfacm_h
#define _chemistry_qc_dft_hfacm_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/dft/integrator.h>
#include <chemistry/qc/dft/functional.h>

//. The \clsnm{HFACM} class performs adiabatic connection method
//calculations using Hartree-Fock orbitals.
class HFACM: public Wavefunction {
#   define CLASSNAME HFACM
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double a0_;
    RefDenIntegrator integrator_;
    RefDenFunctional functional_;
    RefSCF scf_;

    // implement the Compute::compute() function
    virtual void compute();
    double compute_energy();
  public:
    HFACM(StateIn& s);
    HFACM(const RefKeyVal& keyval);
    ~HFACM();
    void save_data_state(StateOut& s);

    int nelectron();
    int spin_polarized();
    RefSymmSCMatrix density();
    RefSymmSCMatrix alpha_density();
    RefSymmSCMatrix beta_density();

    int value_implemented();
    int gradient_implemented();
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
