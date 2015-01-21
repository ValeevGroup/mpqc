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

#include <string>
#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbpt/mbpt.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/twobodygrid.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

class R12IntEval;
class R12WavefunctionWorld;
class MP2R12Energy;

/** The MBPT2_R12 class implements several R12 second-order Moeller-Plesset perturbation theory
methods. */
class MBPT2_R12: public MBPT2 {

    Ref<R12IntEval> r12eval_;           // the R12 intermediates evaluator
    Ref<R12WavefunctionWorld> r12world_;   // parameters for r12eval_

    /** These are MP2-R12 energy objects for each MP2-R12 method, since several different energies
        can be evaluated with the same set of intermediates */
    Ref<MP2R12Energy> r12a_energy_;
    Ref<MP2R12Energy> r12ap_energy_;
    Ref<MP2R12Energy> r12app_energy_;
    Ref<MP2R12Energy> r12b_energy_;
    Ref<MP2R12Energy> r12c_energy_;

    Ref<TwoBodyGrid> twopdm_grid_;      // The set of 2 particle positions on which to compute values of pair function
    unsigned int plot_pair_function_[2];// Which pair function to plot

    static double ref_to_mp2r12_acc() { return 0.01; }

    double mp2_corr_energy_;

    bool cabs_singles_;
    double cabs_singles_energy_;

    std::string occ_orbs_;
    std::string uocc_orbs_;
    double uocc_orbs_pno_truncate_threshold_;

    int gf2_orbital_; // orbital index to use for Green's function calcs, -1 = HOMO, +1 = LUMO, 0 = skip GF2

    // calculate the MP2-R12 energy (or energies, depending on which approximation is chosen)
    void compute_energy_();

    /// computes MP1 pair-natural orbitals of Meier, with or without consideration for the R12 factor
    /// @param spin the spincase
    /// @param deflate_geminal whether to subtract the geminal contribution from the MP1 density; this also causes to add
    ///        the occupied-virtual geminal contribution
    /// @param truncate_threshold omit all PNOs with eigenvalues below this number
    /// @return vector of pairs, each pair contains a coefficient matrix for the MP1 PNOs and the coefficient matrix for the CABS space;
    ///         the number of entries is the number of active orbital pairs
    ///         for this spincase; columns of each matrix are the coefficients of PNO/CABS in the corresponding AO basis
    std::vector< std::pair<RefSCMatrix,RefSCMatrix> >
    mp1_pno(SpinCase2 spin,
            bool deflate_geminal,
            double truncate_threshold = 1e-8);

  protected:
    // implement the Compute::compute() function,
    // overrides MBPT2::compute()
    void compute();

  public:
    MBPT2_R12(StateIn&);
    /** The KeyVal constructor uses keywords of MBPT2, WavefunctionWorld, and R12WavefunctionWorld, and the following keywords
        <dl>

    <dt><tt>twopdm_grid</tt><dd> This optional keyword specifies a TwoBodyGrid object on which to
    plot pair function given by <tt>plot_pair_function</tt>.

    <dt><tt>plot_pair_function</tt><dd> If <tt>twopdm_grid</tt> is given, this array of 2 MO indices
    specifies which pair function to plot.

    <dt><tt>cabs_singles</tt><dd> Evaluate the second-order energy contribution from
    CABS singles and include it into the MBPT(2)-R12 energy. The default is true, except if <tt>corr_factor=none</tt>.

    <dt><tt>occ_orbitals</tt><dd> The method used to provide the occupied orbitals. This can
        be <tt>pipek-mezey</tt> or <tt>canonical</tt> (the default). To be used by developers only.

    <dt><tt>uocc_orbitals</tt><dd> The method used to provide the unoccupied ("virtual") orbitals. This can
        be <tt>pno</tt>, <tt>pno-r12</tt>, or <tt>canonical</tt> (the default). To be used by developers only.

    <dt><tt>pno_truncate_threshold</tt><dd> If <tt>uocc_orbitals</tt> is set to <tt>pno</tt> or <tt>pno-r12</tt>,
    this keyword is used to determine the truncation threshold for the PNO space. The default is 1e-10.
    To be used by developers only.

    </dl> */
    MBPT2_R12(const Ref<KeyVal>&);
    ~MBPT2_R12();

    void save_data_state(StateOut&);

    const Ref<R12WavefunctionWorld>& r12world() const { return r12world_; }
    const Ref<R12IntEval>& r12eval() const { return r12eval_; }
    /// this changes the correlation factor
    void corrfactor(const Ref<R12Technology::CorrelationFactor>&);

    /// total MBPT(2)-R12 energy (does not include CABS singles)
    double corr_energy();
    /// R12 doubles contribution to the MBPT(2)-R12 energy
    double r12_corr_energy();
    /// CABS singles contribution to the total energy
    double cabs_singles_energy();

    RefSymmSCMatrix density();

    void obsolete();
    bool analytic_gradient_implemented() const;
    int value_implemented() const;
    void set_desired_value_accuracy(double acc);

    void print(std::ostream&o=ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
