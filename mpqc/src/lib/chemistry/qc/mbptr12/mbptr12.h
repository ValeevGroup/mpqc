//
// mbptr12.h
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

#ifndef _chemistry_qc_mbptr12_mbptr12_h
#define _chemistry_qc_mbptr12_mbptr12_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbpt/mbpt.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/vxb_eval.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/** The MBPT2_R12 class implements several linear R12 second-order perturbation theory
methods. */
class MBPT2_R12: public MBPT2 {

    Ref<R12IntEval> r12eval_;           // the R12 intermediates evaluator

    Ref<GaussianBasisSet> aux_basis_;
    Ref<SCVector> epair_0_, epair_1_;   // Singlet/triplet pair energies if spin-adapted
                                        // Alpha-beta/alpha-alpha pair energies if spin-orbital

#define ref_to_mp2r12_acc_ 100.0

    double mp2_corr_energy_;
    double r12_corr_energy_;
    LinearR12::StandardApproximation stdapprox_;
    char *r12ints_file_;
    bool spinadapted_;

    void init_variables_();

    // This checks if the integral factory is suitable for R12 calculations
    void check_integral_factory_();

    /* calculate the MP2-R12 energy in std approximations A and A' */
    void compute_energy_a_();

  protected:
    // implement the Compute::compute() function,
    // overrides MBPT2::compute()
    void compute();

  public:
    MBPT2_R12(StateIn&);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>stdapprox</tt><dd> This gives a string that must take on one
        of the values below.  The default is A.

        <dl>

          <dt><tt>A</tt><dd> Use second order M\o{}ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation A.
          Only energies can be computed with this method.

          <dt><tt>A'</tt><dd> Use second order M\o{}ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation A'.
          Only energies can be computed with this method.

          <dt><tt>B</tt><dd> Use second order M\o{}ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation B. 
	  This method is not implemented yet.

        </dl>

        <dt><tt>r12ints_file</tt><dd> This specifies which file to use to store transformed
	MO integrals if the multipass algorithm is used.  The default is /tmp/r12ints.dat.

        </dl> */
    MBPT2_R12(const Ref<KeyVal>&);
    ~MBPT2_R12();

    void save_data_state(StateOut&);

    Ref<GaussianBasisSet> aux_basis() const;

    double corr_energy();
    double r12_corr_energy();

    RefSymmSCMatrix density();

    void obsolete();
    int gradient_implemented() const;
    int value_implemented() const;

    void print(std::ostream&o=ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
