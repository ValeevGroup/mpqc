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
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>

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
    Ref<MP2R12Energy> r12b_energy_;

    Ref<GaussianBasisSet> aux_basis_;
    Ref<SCVector> epair_0_, epair_1_;   // Singlet/triplet pair energies if spin-adapted
                                        // Alpha-beta/alpha-alpha pair energies if spin-orbital

#define ref_to_mp2r12_acc_ 100.0

    double mp2_corr_energy_;
    double r12_corr_energy_;
    LinearR12::StandardApproximation stdapprox_;
    R12IntEvalInfo::StoreMethod r12ints_method_;
    char* r12ints_file_;
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
	  with linear R12 terms in standard approximation A (MP2-R12/A).
          Only energies can be computed with the MP2-R12/A method.

          <dt><tt>A'</tt><dd> Use second order M\o{}ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation A' (MP2-R12/A').
          This will cause MP2-R12/A energies to be computed also.
          Only energies can be computed with the MP2-R12/A' method.

          <dt><tt>B</tt><dd> Use second order M\o{}ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation B. 
	  This method is not implemented yet.

        </dl>

	<dt><tt>spinadapted</tt><dd> This specifies whether to compute spin-adapted
	or spin-orbital pair energies. Default is to compute spin-adapted energies.

	<dt><tt>aux_basis</tt><dd> This specifies the auxiliary basis to be used for the resolution
	of the identity. Default is to use the same basis as for the orbital expansion.

	<dt><tt>r12ints</tt><dd> This specifies how to store transformed MO integrals.
	Valid values are:

	<dl>

	  <dt><tt>mem-posix</tt><dd> Store integrals in memory for single-pass situations
	  and in a binary file on task 0's node using POSIX I/O for multipass situations.
	  <tt>posix</tt> is usually less efficient than <tt>mpi</tt> for distributed
	  parallel multipass runs since the I/O is performed by one task only. However, this method guaranteed to
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
	  transformations. This method should only be used for testing purposes

	</dl>

	If <tt>r12ints</tt> is not specified, then <tt>mem-posix</tt> method will be used.
	If user wishes to use MPI-I/O, pending its availability, for higher parallel efficiency,
	<tt>r12ints</tt> should be explicitly set to <tt>mem-mpi</tt>.

        <dt><tt>r12ints_file</tt><dd> This specifies which file to use to store transformed
	MO integrals if <tt>r12ints=posix-io</tt> or <tt>r12ints=mpi-io</tt> is used.
	Default is "./<inputbasename>.r12ints.dat", where <inputbasename> is the name of the input
	file without ".in". If MPI-I/O is used then it is user's responsibility to ensure
	that the file resides on a file system that supports MPI-I/O.

        </dl> */
    MBPT2_R12(const Ref<KeyVal>&);
    ~MBPT2_R12();

    void save_data_state(StateOut&);

    Ref<GaussianBasisSet> aux_basis() const;
    LinearR12::StandardApproximation stdapprox() const;
    bool spinadapted() const;
    R12IntEvalInfo::StoreMethod r12ints_method() const;
    char* r12ints_file() const;

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
