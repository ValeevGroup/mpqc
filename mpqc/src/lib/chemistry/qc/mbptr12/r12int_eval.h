//
// r12int_eval.h
//
// Copyright (C) 2004 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_r12inteval_h
#define _chemistry_qc_mbptr12_r12inteval_h

#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/r12_amps.h>

namespace sc {

  /** R12IntEval is the top-level class which computes intermediates occuring in linear R12 theories.
      This class is used by all Wavefunction classes that implement linear R12 methods.
  */

class R12IntEval : virtual public SavableState {

  bool evaluated_;

  // Calculation information (number of basis functions, R12 approximation, etc.)
  Ref<R12IntEvalInfo> r12info_;

  // Note that intermediate B is symmetric but is stored as a full matrix to simplify the code
  // that computes asymmetric form of B
  RefSCMatrix Vaa_, Vab_, Xaa_, Xab_, Baa_, Bab_, Aaa_, Aab_, T2aa_, T2ab_;
  RefSCMatrix Ac_aa_, Ac_ab_;
  RefSCMatrix Raa_, Rab_;    // Not sure if I'll compute and keep these explicitly later
  Ref<R12Amplitudes> Amps_;  // Amplitudes of various R12-contributed terms in pair functions
  RefSCVector emp2pair_aa_, emp2pair_ab_;
  RefSCDimension dim_ij_aa_, dim_ij_ab_, dim_ij_s_, dim_ij_t_;
  RefSCDimension dim_ab_aa_, dim_ab_ab_;

  bool gbc_;
  bool ebc_;
  LinearR12::ABSMethod abs_method_;
  LinearR12::StandardApproximation stdapprox_;
  bool spinadapted_;
  bool include_mp1_;
  /// should I follow Klopper-Samson approach in the intermediates formulation for the EBC-free method?
  bool follow_ks_ebcfree_;
  int debug_;

  // Map to TwoBodyMOIntsTransform objects that have been computed previously
  typedef std::map<std::string, Ref<TwoBodyMOIntsTransform> > TformMap;
  TformMap tform_map_;
  // Returns pointer to the appropriate transform.
  // If the transform is not found then throw runtime_error
  Ref<TwoBodyMOIntsTransform> get_tform_(const std::string&);

  /// Fock-weighted occupied space |i_f> = f_i^R |R>, where R is a function in RI-BS
  Ref<MOIndexSpace> focc_space_;
  /// Form Fock-weighted occupied space
  void form_focc_space_();
  /// Fock-weighted active occupied space |i_f> = f_i^R |R>, where R is a function in RI-BS
  Ref<MOIndexSpace> factocc_space_;
  /// Form Fock-weighted active occupied space
  void form_factocc_space_();
  /// Space of canonical virtual MOs
  Ref<MOIndexSpace> canonvir_space_;
  /// Form space of auxiliary virtuals
  void form_canonvir_space_();

  /// Initialize standard transforms
  void init_tforms_();
  /// Set intermediates to zero + add the "diagonal" contributions
  void init_intermeds_();
  /// Compute r^2 contribution to X
  void r2_contrib_to_X_orig_();
  /// Compute r^2 contribution to X using compute_r2_()
  void r2_contrib_to_X_new_();
  /// Compute <space1 space1|r_{12}^2|space1 space2> matrix
  RefSCMatrix compute_r2_(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2);
  /** Compute the Fock matrix between 2 spaces. scale_J and scale_K are used to scale Coulomb
      and exchange contributions */
  RefSCMatrix fock_(const Ref<MOIndexSpace>& occ_space, const Ref<MOIndexSpace>& bra_space,
                    const Ref<MOIndexSpace>& ket_space, double scale_J = 1.0, double scale_K = 1.0);
  /// Compute the coulomb matrix between 2 spaces
  RefSCMatrix coulomb_(const Ref<MOIndexSpace>& occ_space, const Ref<MOIndexSpace>& bra_space,
                       const Ref<MOIndexSpace>& ket_space);
  /// Compute the exchange matrix between 2 spaces
  RefSCMatrix exchange_(const Ref<MOIndexSpace>& occ_space, const Ref<MOIndexSpace>& bra_space,
                        const Ref<MOIndexSpace>& ket_space);
  /// Compute and print out neatly various matrix norms of A
  void compute_norms_(const RefSCMatrix& A, const std::string& label, std::ostream& os = ExEnv::out0());

  /// Checkpoint the top-level molecular energy
  void checkpoint_() const;

  /** Compute the number tasks which have access to the integrals.
      map_to_twi is a vector<int> each element of which will hold the
      number of tasks with access to the integrals of lower rank than that task
      (or -1 if the task doesn't have access to the integrals) */
  const int tasks_with_ints_(const Ref<R12IntsAcc> ints_acc, vector<int>& map_to_twi);
  
  /** Compute contribution to V, X, and B of the following form:
      0.5 * \bar{g}_{ij}^{pq} * \bar{r}_{pq}^{kl}, where p and q span mospace.
      tform_name is the name of the transform to be used to get the integrals.
      mospace is either space2() or space4() of that transform
  */
  void contrib_to_VXB_a_symm_(const std::string& tform_name);

  /** Compute contribution to V, X, and B of the following form:
      \bar{g}_{ij}^{am} * \bar{r}_{am}^{kl}, where m and a span space1 and space2, respectively.
      tform_name is the name of the transform to be used to get the integrals.
      mospace1 and mospace2 are space2() and space4() of that transform, respectively
  */
  void contrib_to_VXB_a_asymm_(const std::string& tform_name);

  /// Compute OBS contribution to V, X, and B (these contributions are independent of the method)
  void obs_contrib_to_VXB_gebc_vbseqobs_();

  /// Compute 1-ABS contribution to V, X, and B (these contributions are independent of the method)
  void abs1_contrib_to_VXB_gebc_();

  /// Equiv to the sum of above, except for this doesn't assume that VBS is the same as OBS
  void contrib_to_VXB_gebc_vbsneqobs_();

  /// Compute A using the "simple" formula obtained using direct substitution alpha'->a'
  void compute_A_simple_();

  /// Compute A using the standard commutator approach
  void compute_A_via_commutator_();

  /// Compute MP2 T2
  void compute_T2_();

  /// Compute R "intermediate" (r12 integrals in occ-pair/vir-pair basis)
  void compute_R_();

  /// Compute A*T2 contribution to V (needed if EBC is not assumed)
  void AT2_contrib_to_V_();

  /// Compute -2*A*R contribution to B (needed if EBC is not assumed)
  void AR_contrib_to_B_();

  /** Compute the first (r<sub>kl</sub>^<sup>AB</sup> f<sub>A</sub><sup>m</sup> r<sub>mB</sub>^<sup>ij</sup>)
      contribution to B that vanishes under GBC */
  void compute_B_gbc_1_();
  
  /** Compute the second (r<sub>kl</sub>^<sup>AB</sup> r<sub>AB</sub>^<sup>Kj</sup> f<sub>K</sub><sup>i</sup>)
      contribution to B that vanishes under GBC */
  void compute_B_gbc_2_();

  /// Compute dual-basis MP2 energy (contribution from doubles to MP2 energy)
  void compute_dualEmp2_();
  
  /// Compute dual-basis MP1 energy (contribution from singles to HF energy)
  void compute_dualEmp1_();
 
  /// This function computes T2 amplitudes
  void compute_T2_vbsneqobs_();

  /** New general function to compute <ij|r<sub>12</sub>|pq> integrals. ipjq_tform
      is the source of the integrals.*/
  void compute_R_vbsneqobs_(const Ref<TwoBodyMOIntsTransform>& ipjq_tform,
                            RefSCMatrix& Raa, RefSCMatrix& Rab);

  /** Initialize amplitude objects */
  void compute_amps_();

  /** Sum contributions to SCMatrix A from all nodes and broadcast so
      every node has the correct SCMatrix. If to_all_tasks is false, then
      collect all contributions to task 0 and zero out the matrix on other tasks,
      otherwise distribute the sum to every task. If to_average is true then
      each result is scaled by the inverse of the number of tasks. */
  void globally_sum_scmatrix_(RefSCMatrix& A, bool to_all_tasks = false, bool to_average = false);

  /** Sum contributions to SCVector A from all nodes and broadcast so
      every node has the correct SCVector. If to_all_tasks is false, then
      collect all contributions to task 0 and zero out the matrix on other tasks,
      otherwise distribute the sum to every task. If to_average is true then 
      each result is scaled by the inverse of the number of tasks. */
  void globally_sum_scvector_(RefSCVector& A, bool to_all_tasks = false, bool to_average = false);

  /** Sum contributions to the intermediates from all nodes and broadcast so
      every node has the correct matrices. If to_all_tasks is false, then
      collect all contributions to task 0 and zero out the matrix on other tasks,
      otherwise distribute the sum to every task. */
  void globally_sum_intermeds_(bool to_all_tasks = false);

public:
  R12IntEval(StateIn&);
  /** Constructs R12IntEval. If follow_ks_ebcfree is true then follow formalism of Klopper and Samson
      to compute EBC-free MP2-R12 energy. */
  R12IntEval(const Ref<R12IntEvalInfo>& info, bool gbc = true, bool ebc = true,
             LinearR12::ABSMethod abs_method = LinearR12::ABS_CABSPlus,
             LinearR12::StandardApproximation stdapprox = LinearR12::StdApprox_Ap,
             bool follow_ks_ebcfree = false);
  ~R12IntEval();

  void save_data_state(StateOut&);
  virtual void obsolete();

  void include_mp1(bool include_mp1);
  void set_debug(int debug);
  void set_dynamic(bool dynamic);
  void set_print_percent(double print_percent);
  void set_memory(size_t nbytes);

  const bool gbc() const { return gbc_; }
  const bool ebc() const { return ebc_; }
  const LinearR12::StandardApproximation stdapprox() const { return stdapprox_; }
  bool follow_ks_ebcfree() const { return follow_ks_ebcfree_; }

  Ref<R12IntEvalInfo> r12info() const;
  RefSCDimension dim_oo_aa() const;
  RefSCDimension dim_oo_ab() const;
  RefSCDimension dim_oo_s() const;
  RefSCDimension dim_oo_t() const;
  RefSCDimension dim_vv_aa() const;
  RefSCDimension dim_vv_ab() const;

  /// This function causes the intermediate matrices to be computed.
  virtual void compute();

  /// Returns alpha-alpha block of the V intermediate matrix.
  RefSCMatrix V_aa();
  /// Returns alpha-alpha block of the X intermediate matrix.
  RefSCMatrix X_aa();
  /// Returns alpha-alpha block of the B intermediate matrix.
  RefSymmSCMatrix B_aa();
  /// Returns alpha-alpha block of the A intermediate matrix. Returns 0 if EBC is assumed.
  RefSCMatrix A_aa();
  /// Returns alpha-alpha block of the A intermediate matrix. Returns 0 if EBC is assumed.
  RefSCMatrix Ac_aa();
  /// Returns alpha-alpha block of the MP2 T2 matrix. Returns 0 if EBC is assumed.
  RefSCMatrix T2_aa();
  /// Returns alpha-beta block of the V intermediate matrix.
  RefSCMatrix V_ab();
  /// Returns alpha-beta block of the X intermediate matrix.
  RefSCMatrix X_ab();
  /// Returns alpha-beta block of the B intermediate matrix.
  RefSymmSCMatrix B_ab();
  /// Returns alpha-beta block of the A intermediate matrix. Returns 0 if EBC is assumed
  RefSCMatrix A_ab();
  /// Returns alpha-beta block of the A intermediate matrix. Returns 0 if EBC is assumed
  RefSCMatrix Ac_ab();
  /// Returns alpha-beta block of the MP2 T2 matrix. Returns 0 if EBC is assumed
  RefSCMatrix T2_ab();
  /// Returns alpha-alpha MP2 pair energies.
  RefSCVector emp2_aa();
  /// Returns alpha-beta MP2 pair energies.
  RefSCVector emp2_ab();
  /// Returns amplitudes of pair correlation functions
  Ref<R12Amplitudes> amps();

  RefDiagSCMatrix evals() const;  
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


