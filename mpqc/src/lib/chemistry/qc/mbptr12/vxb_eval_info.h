//
// vxb_eval_info.h
//
// Copyright (C) 2003 Edward Valeev
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _chemistry_qc_mbptr12_vxbevalinfo_h
#define _chemistry_qc_mbptr12_vxbevalinfo_h

#include <util/ref/ref.h>
#include <math/scmat/abstract.h>
#include <util/group/memory.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/moindexspace.h>

namespace sc {

class MBPT2_R12;

  /** Class R12IntEvalInfo contains information necessary for R12 intermediate
      evaluators */

class R12IntEvalInfo : virtual public SavableState {

public:

  /// Describes the method of storing transformed MO integrals. See MBPT2_R12.
  enum StoreMethod { mem_posix = 0, posix = 1, mem_mpi = 2, mpi = 3, mem_only = 4 };

private:

  Wavefunction* wfn_;     // Wavefunction that owns this
  Ref<SCF> ref_;
  Ref<Integral> integral_;
  Ref<GaussianBasisSet> bs_;
  Ref<GaussianBasisSet> bs_aux_;
  Ref<GaussianBasisSet> bs_ri_;
  Ref<SCMatrixKit> matrixkit_;
  Ref<MessageGrp> msg_;
  Ref<MemoryGrp> mem_;
  Ref<ThreadGrp> thr_;

  int nocc_;
  int nocc_act_;
  int nfzc_;
  int nfzv_;
  int noso_;

  size_t memory_;
  bool dynamic_;
  double print_percent_;
  int debug_;
  StoreMethod ints_method_;
  char* ints_file_;
  LinearR12::ABSMethod abs_method_;

  RefSCMatrix scf_vec_;
  RefDiagSCMatrix evals_;
  RefDiagSCMatrix occs_;
  int* orbsym_;
  RefSCMatrix orthog_aux_;
  int nlindep_aux_;
  RefSymmSCMatrix overlap_aux_;
  RefSCMatrix orthog_ri_;
  RefSCMatrix orthog_ri_so_;
  int nlindep_ri_;
  RefSymmSCMatrix overlap_ri_;
  
  Ref<MOIndexSpace> obs_space_;
  Ref<MOIndexSpace> abs_space_;
  Ref<MOIndexSpace> ribs_space_;
  Ref<MOIndexSpace> act_occ_space_;

  // construct the RI basis based on abs_method
  void construct_ri_basis_(bool safe);
  void construct_ri_basis_ks_(bool safe);
  void construct_ri_basis_ksplus_(bool safe);
  void construct_ri_basis_ev_(bool safe);
  void construct_ri_basis_evplus_(bool safe);
  // Uses ri_basis to construct a basis that spans the orthogonal complement to the OBS
  void construct_ortho_comp_();
  // Same as above, except for this function uses Curt's approach
  void construct_ortho_comp_jv_();
  // Same as above, except for it uses pure SVD
  void construct_ortho_comp_svd_();
  // Returns true if ABS spans OBS
  bool abs_spans_obs_();
  // OBSOLETE: Construct eigenvector and eigenvalues sorted by energy
  void eigen_(RefDiagSCMatrix& evals, RefSCMatrix& scf_vec, RefDiagSCMatrix& occs, int*& orbsym);
  // Construct eigenvector and eigenvalues sorted by energy
  void eigen2_(RefDiagSCMatrix& evals, RefSCMatrix& scf_vec, int*& orbsym);
  // Construct orthog_aux_
  void construct_orthog_aux_();
  // Construct orthog_ri_
  void construct_orthog_ri_();

public:
  R12IntEvalInfo(StateIn&);
  /// Constructs an R12IntEvalInfo object using data from the MBPT2_R12 object
  R12IntEvalInfo(MBPT2_R12*);
  ~R12IntEvalInfo();

  void save_data_state(StateOut&);

  /** Sets whether to use dynamic load balancing in parallel MO transformations.
      Default is no */
  void set_dynamic(bool dynamic) { dynamic_ = dynamic; };
  /// Sets how frequently updates of progress are printed out. Default is 10%
  void set_print_percent(double print_percent) { print_percent_ = print_percent; };
  /// Set debug level. Default is 0.
  void set_debug_level(int debug) { debug_ = debug; };
  /** Sets the method of storing transformed MO integrals. Default depends on
      how the object was constructed. */
  void set_ints_method(const StoreMethod method) { ints_method_ = method; };
  /** Sets name of the file used to store transformed integrals.
      Default depends on how the object was constructed. */
  void set_ints_file(const char* filename) { ints_file_ = strdup(filename); };
  /** Sets the amount of memory to use for the calculation. Default is
      determined by DEFAULT_SC_MEMORY. */
  void set_memory(size_t nbytes) { if (nbytes >= 0) memory_ = nbytes; };
  /** Sets the ABS approach to be used (ABS or CABS).
      Default depends on how the object was constructed. */
  void set_absmethod(LinearR12::ABSMethod abs_method);

  Wavefunction* wfn() const { return wfn_; };
  Ref<SCF> ref() const { return ref_; };
  Ref<Integral> integral() const { return integral_; };
  /// Returns the orbital basis set (OBS) object
  Ref<GaussianBasisSet> basis() const { return bs_; };
  /// Returns the resolution-of-the-identity basis set (RIBS) object
  Ref<GaussianBasisSet> basis_ri() const { return bs_ri_; };
  Ref<SCMatrixKit> matrixkit() const { return matrixkit_; };
  Ref<MemoryGrp> mem() const { return mem_;};
  Ref<MessageGrp> msg() const { return msg_;};
  Ref<ThreadGrp> thr() const { return thr_;};

  bool dynamic() const { return dynamic_; };
  double print_percent() const { return print_percent_; };
  int debug_level() const { return debug_; };
  const StoreMethod ints_method() const { return ints_method_; };
  char* ints_file() const;
  const size_t memory() const { return memory_; };

  const int nocc() const { return nocc_;};
  const int nocc_act() const { return nocc_act_;};
  const int noso() const { return noso_;};
  const int nfzc() const { return nfzc_;};
  const int nfzv() const { return nfzv_;};

  LinearR12::ABSMethod abs_method() const { return abs_method_; };

  RefSCMatrix scf_vec() const { return scf_vec_; };
  RefDiagSCMatrix evals() const { return evals_; };
  int *orbsym() const { return orbsym_; };
  RefSCMatrix orthog_ri();
  
  /// Returns the MOIndexSpace object for OBS
  Ref<MOIndexSpace> obs_space() const { return obs_space_; };
  /// Returns the MOIndexSpace object for the active occupied MOs
  Ref<MOIndexSpace> act_occ_space() const { return act_occ_space_; };
  /// Returns the MOIndexSpace object for ABS
  Ref<MOIndexSpace> abs_space() const { return abs_space_; };
  /// Returns the MOIndexSpace object for RI-BS
  Ref<MOIndexSpace> ribs_space() const { return ribs_space_; };

  /// Compute dipole and quadrupole moment matrices in active MO basis
  void compute_multipole_ints(RefSymmSCMatrix& MX,
			      RefSymmSCMatrix& MY,
			      RefSymmSCMatrix& MZ,
			      RefSymmSCMatrix& MXX,
			      RefSymmSCMatrix& MYY,
			      RefSymmSCMatrix& MZZ);
			      
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


