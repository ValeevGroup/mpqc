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

class R12IntEval;
class R12IntEvalInfo;
class MP2R12Energy;
  
/** The MBPT2_R12 class implements several linear R12 second-order perturbation theory
methods. */
class MBPT2_R12: public MBPT2 {

    Ref<R12IntEval> r12eval_;           // the R12 intermediates evaluator

    /** These are MP2-R12 energy objects for each MP2-R12 method, since several different energies
        can be evaluated with the same set of intermediates */
    Ref<MP2R12Energy> r12a_energy_;
    Ref<MP2R12Energy> r12ap_energy_;
    Ref<MP2R12Energy> r12app_energy_;
    Ref<MP2R12Energy> r12b_energy_;
    Ref<MP2R12Energy> r12c_energy_;

    Ref<GaussianBasisSet> aux_basis_;   // This is the auxiliary basis set (ABS)
    Ref<GaussianBasisSet> vir_basis_;   // This is the virtuals basis set (VBS)
    Ref<SCVector> epair_0_, epair_1_;   // Singlet/triplet pair energies if spin-adapted
                                        // Alpha-beta/alpha-alpha pair energies if spin-orbital

    Ref<TwoBodyGrid> twopdm_grid_;      // The set of 2 particle positions on which to compute values of pair function
    unsigned int plot_pair_function_[2];// Which pair function to plot

#define ref_to_mp2r12_acc_ 100.0

    double mp2_corr_energy_;
    double r12_corr_energy_;
    Ref<LinearR12::CorrelationFactor> corrfactor_;
    LinearR12::StandardApproximation stdapprox_;
    Ref<LinearR12Ansatz> ansatz_;
    LinearR12::ABSMethod abs_method_;
    unsigned int maxnabs_;
    R12IntEvalInfo::StoreMethod::type r12ints_method_;
    std::string r12ints_file_;
    bool gbc_;
    bool ebc_;
    bool ks_ebcfree_;
    bool omit_P_;
    bool spinadapted_;
    bool include_mp1_;

    void init_variables_();

    // This checks if the integral factory is suitable for R12 calculations
    void check_integral_factory_();

    // calculate the MP2-R12 energy (or energies, depending on which approximation is chosen)
    void compute_energy_();

  protected:
    // implement the Compute::compute() function,
    // overrides MBPT2::compute()
    void compute();

  public:
    MBPT2_R12(StateIn&);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>corr_factor</tt><dd> This string specifies which correlation factor to use.
        Allowed values are "r12", "g12", and "none". The default is "r12".
      
        <dt><tt>corr_param</tt><dd> This keyword specifies optional parameters
        of the correlation factor. <tt>corr_param</tt> can be a single floating-point value
        an array of floating-point values, or an array of arrays of 2-element arrays of
        floating-point values. Single value specifies the parameter of the single
        correlation function. The 1-d array form specifies a set of primitive correlation functions
        characterized by the corresponding parameters. The 3-d array form specifies
        a set of contracted correlation functions. For example,
        <tt>corr_param = 3.0</tt> specifies a single correlation function
        with parameter 3.0. <tt>corr_param = [ 1.0 3.0 10.0 ]</tt> specifies
        3 correlation functions with parameters 1.0, 3.0 and 10.0.
        <tt>corr_param = [ [[1.0 0.35][3.0 0.65]]  [[10.0 1.0]] ]</tt>
        specifies 2 correlation functions, first composed of 2 primitive functions
        with parameters 1.0 and 3.0 combined linearly with coefficients
        0.35 and 0.65, and second primitive function with parameter 10.0 .

        This keyword has no meaning for some correlation factors, e.g., "r12" and "none",
        and is not used. There is no default.
      
        <dt><tt>stdapprox</tt><dd> This gives a string that must take on one
        of the values below.  The default is A'.

        <dl>

          <dt><tt>A</tt><dd> Use second order M&oslash;ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation A (MP2-R12/A).
          Only energies can be computed with the MP2-R12/A method.

          <dt><tt>A'</tt><dd> Use second order M&oslash;ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation A' (MP2-R12/A').
          This will cause MP2-R12/A energies to be computed also.
          Only energies can be computed with the MP2-R12/A' method.

          <dt><tt>A''</tt><dd> Use second order M&oslash;ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation A'' (MP2-R12/A'').
          Only energies can be computed with the MP2-R12/A'' method.

          <dt><tt>B</tt><dd> Use second order M&oslash;ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation B.
          This will cause A and A' energies to be computed also.
          Only energies can be computed with the MP2-R12/B method.

          <dt><tt>C</tt><dd> Use second order M&oslash;ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation C.
          Only energies can be computed with the MP2-R12/C method.

        </dl>
        
        <dt><tt>ansatz</tt><dd> This object specifies the ansatz (see LinearR12Ansatz).
        
	<dt><tt>spinadapted</tt><dd> This boolean specifies whether to compute spin-adapted
	or spin-orbital pair energies. Default is to compute spin-adapted energies for closed-shell
        systems and spin-orbital energies for open-shell systems. For some references, e.g. UHF, this keyword
        is not used.

        <dt><tt>gbc</tt><dd> This boolean specifies whether Generalized Brillouin
        Condition (GBC) is assumed to hold. The default is "true". This keyword is
        only valid if stdapprox=A'.
        The effect of setting this keyword to true is very small --
        hence it is not recommended to use this keyword.

        <dt><tt>ebc</tt><dd> This boolean specifies whether Extended Brillouin
        Condition (EBC) is assumed to hold. The default is "true". This keyword
        is only valid if stdapprox=A'.
        The effect of setting this keyword to true is small --
        hence it is not recommended to use this keyword.
       
        <dt><tt>maxnabs</tt><dd> This integer specifies the maximum number of ABS indices per integral.
        Valid values are between 1 and 2. The default is to include all terms necessary for a given method.
        For example, MP2-F12/B energy involves integrals with 2 ABS indices. Setting maxnabs to 1
        will leave out such terms.
       
	<dt><tt>aux_basis</tt><dd> This specifies the auxiliary AO basis to be used for the resolution
	of the identity. Default is to use the same basis as for the orbital expansion.

	<dt><tt>vir_basis</tt><dd> This specifies the AO basis to be used for the virtual orbitals.
	Default is to use the same basis as for the orbital expansion.

        <dt><tt>abs_method</tt><dd> This string specifies whether the old ABS method, introduced
        by Klopper and Samson, or the new ABS variant, CABS, introduced by Valeev, should be used.
	Valid values are "ABS" (Klopper and Samson), "ABS+", "CABS", and "CABS+", where the "+" labels
	a method where the union of OBS and ABS is used to construct the RI basis. The default is "ABS".
        The default in 2.3.0 and later will be "CABS+".

        <dt><tt>lindep_tol</tt><dd> The tolerance used to detect linearly
        dependent basis functions in the RI basis set.
        The precise meaning depends on the
        orthogonalization method.  The default value is 1e-8.

	<dt><tt>r12ints</tt><dd> This specifies how to store transformed MO integrals.
	Valid values are:

	<dl>

	  <dt><tt>mem-posix</tt><dd> Store integrals in memory for single-pass situations
	  and in a binary file on task 0's node using POSIX I/O for multipass situations.
	  <tt>posix</tt> is usually less efficient than <tt>mpi</tt> for distributed
	  parallel multipass runs since the I/O is performed by one task only. However, this method is guaranteed to
          work in all types of environments, hence <tt>mem-posix</tt> is the default.

	  <dt><tt>posix</tt><dd> Store integrals in a binary file on task 0's node using POSIX I/O.
	  This method is different from <tt>mem-posix</tt> in that it forces the integrals out to disk
          even if they could be stored in memory. <tt>posix</tt> should only be used for benchmarking
          and testing purposes.

	  <dt><tt>mem-mpi</tt><dd> Store integrals in memory for single-pass situations
	  and in a binary file using MPI-I/O for multipass situations. This method
	  assumes the availability of MPI-I/O. <tt>mem-mpi</tt> is the preferred choice
	  in distributed environments which have MPI-I/O available.

	  <dt><tt>mpi</tt><dd> Store integrals in a binary file using MPI-I/O. This method
	  is different from <tt>mem-mpi</tt> in that it forces the integrals out to disk
	  even if they could be stored in memory. <tt>mpi</tt> should only be used for benchmarking
	  and testing purposes.

	  <dt><tt>mem</tt><dd> Store integrals in memory. Can only be used with single-pass
	  transformations for MP2-R12/A and MP2-R12/A' methods.
          This method should only be used for testing purposes.

	</dl>

	If <tt>r12ints</tt> is not specified, then <tt>mem-posix</tt> method will be used.
	If user wishes to use MPI-I/O, pending its availability, for higher parallel efficiency,
	<tt>r12ints</tt> should be explicitly set to <tt>mem-mpi</tt>.

        <dt><tt>r12ints_file</tt><dd> This specifies the prefix for the transformed
	MO integrals file if <tt>r12ints</tt> is set to <tt>posix</tt>, <tt>mpi</tt>, <tt>mem-posix</tt>
        or <tt>mem-mpi</tt> is used.
	Default is "./<inputbasename>.r12ints", where <inputbasename> is the name of the input
	file without ".in". If MPI-I/O is used then it is user's responsibility to ensure
	that the file resides on a file system that supports MPI-I/O.

        <dt><tt>twopdm_grid</tt><dd> This optional keyword specifies a TwoBodyGrid object on which to
        plot pair function given by <tt>plot_pair_function</tt>.

        <dt><tt>plot_pair_function</tt><dd> If <tt>twopdm_grid</tt> is given, this array of 2 MO indices
        specifies which pair function to plot.

        </dl> */
    MBPT2_R12(const Ref<KeyVal>&);
    ~MBPT2_R12();

    void save_data_state(StateOut&);

    const Ref<LinearR12::CorrelationFactor>& corrfactor() const;
    Ref<GaussianBasisSet> aux_basis() const;
    Ref<GaussianBasisSet> vir_basis() const;
    unsigned int maxnabs() const;
    bool gbc() const;
    bool ebc() const;
    LinearR12::ABSMethod abs_method() const;
    LinearR12::StandardApproximation stdapprox() const;
    const Ref<LinearR12Ansatz>& ansatz() const;
    bool spinadapted() const;
    bool ks_ebcfree() const;
    bool omit_P() const;
    bool include_mp1() const;
    R12IntEvalInfo::StoreMethod::type r12ints_method() const;
    const std::string& r12ints_file() const;

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
