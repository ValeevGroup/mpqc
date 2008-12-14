//
// mbptr12.h
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

#ifndef _chemistry_qc_mbptr12_mbptr12_h
#define _chemistry_qc_mbptr12_mbptr12_h

#ifdef __GNUC__
#pragma interface
#endif

#include <string>
#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbpt/mbpt.h>
#include <chemistry/qc/mbptr12/linearr12.h>
//#include <chemistry/qc/mbptr12/vxb_eval.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/twobodygrid.h>
#include <chemistry/qc/mbptr12/ansatz.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

class R12Technology;
class R12IntEval;
class R12IntEvalInfo;
class MP2R12Energy;
  
/** The MBPT2_R12 class implements several linear R12 second-order perturbation theory
methods. */
class MBPT2_R12: public MBPT2 {

    bool spinadapted_;
    bool include_mp1_;
    
    bool hylleraas_;
    
    bool new_energy_;

    Ref<R12IntEval> r12eval_;           // the R12 intermediates evaluator
    Ref<R12IntEvalInfo> r12evalinfo_;   // parameters for r12eval_

    /** These are MP2-R12 energy objects for each MP2-R12 method, since several different energies
        can be evaluated with the same set of intermediates */
    Ref<MP2R12Energy> r12a_energy_;
    Ref<MP2R12Energy> r12ap_energy_;
    Ref<MP2R12Energy> r12app_energy_;
    Ref<MP2R12Energy> r12b_energy_;
    Ref<MP2R12Energy> r12c_energy_;

    Ref<TwoBodyGrid> twopdm_grid_;      // The set of 2 particle positions on which to compute values of pair function
    unsigned int plot_pair_function_[2];// Which pair function to plot

#define ref_to_mp2r12_acc_ 100.0

    double mp2_corr_energy_;

    // calculate the MP2-R12 energy (or energies, depending on which approximation is chosen)
    void compute_energy_();

  protected:
    // implement the Compute::compute() function,
    // overrides MBPT2::compute()
    void compute();

  public:
    MBPT2_R12(StateIn&);
    /** The KeyVal constructor uses keywords of MBPT2 and R12IntEvalInfo and the following keywords
        <dl>

	<dt><tt>spinadapted</tt><dd> This boolean specifies whether to compute spin-adapted
	or spin-orbital pair energies. Default is to compute spin-adapted energies for closed-shell
    systems and spin-orbital energies for open-shell systems. For some references, e.g. UHF, this keyword
    is not used.

    <dt><tt>twopdm_grid</tt><dd> This optional keyword specifies a TwoBodyGrid object on which to
    plot pair function given by <tt>plot_pair_function</tt>.

    <dt><tt>plot_pair_function</tt><dd> If <tt>twopdm_grid</tt> is given, this array of 2 MO indices
    specifies which pair function to plot.
    
    <dt><tt>new_energy</tt><dd> Use new version that has got improved diagonal Ansaetze preserving
    spin as well as a fixed coefficient version. For the non diagonal ansaetze the old version is to be 
    preferred since it has many more security checks, such as checks of positive definiteness of B etc.
    The default value of this variable is false.
    
    <dt><tt>hylleraas</tt><dd> If new_energy = true, and the ansatz sets diag and fixed_coeff to true, this variable
    tells the new version to compute the F12 energy according to the Hylleraas functional. The default
    is to use the Hylleraas functional if the ansatz defines fixedcoef = true.

    </dl> */
    MBPT2_R12(const Ref<KeyVal>&);
    ~MBPT2_R12();

    void save_data_state(StateOut&);

    const Ref<R12IntEvalInfo>& r12evalinfo() const { return r12evalinfo_; }
    const Ref<R12IntEval>& r12eval() const { return r12eval_; }
    /// this changes the correlation factor
    void corrfactor(const Ref<LinearR12::CorrelationFactor>&);

    // use Hylleraas functional?
    bool hylleraas() const;
    
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
