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

#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/r12_amps.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/mbptr12/spin.h>

#ifndef _chemistry_qc_mbptr12_r12inteval_h
#define _chemistry_qc_mbptr12_r12inteval_h

namespace sc {
  
  class TwoBodyMOIntsTransform;
  class R12IntsAcc;

  /** R12IntEval is the top-level class which computes intermediates occuring in linear R12 theories.
      This class is used by all Wavefunction classes that implement linear R12 methods.
  */

class R12IntEval : virtual public SavableState {
  private:
  bool evaluated_;

  // Calculation information (number of basis functions, R12 approximation, etc.)
  Ref<R12IntEvalInfo> r12info_;

  RefSCMatrix V_[NSpinCases2];
  RefSCMatrix X_[NSpinCases2];
  // Note that intermediate B is symmetric but is stored as a full matrix to simplify the code
  // that computes asymmetric form of B
  RefSCMatrix B_[NSpinCases2];
  RefSCMatrix A_[NSpinCases2];
  RefSCMatrix Ac_[NSpinCases2];
  RefSCMatrix T2_[NSpinCases2];
  RefSCMatrix F12_[NSpinCases2];
  RefSCVector emp2pair_[NSpinCases2];
  RefSCDimension dim_oo_[NSpinCases2];
  RefSCDimension dim_vv_[NSpinCases2];
  RefSCDimension dim_f12_[NSpinCases2];
  
#if 0
  RefSCMatrix Raa_, Rab_;    // Not sure if I'll compute and keep these explicitly later
#endif
  //Ref<F12Amplitudes> Amps_;  // First-order amplitudes of various contributions to the pair functions
#if 0
  RefSCVector emp2pair_aa_, emp2pair_ab_, emp2pair_bb_;
  RefSCDimension dim_ij_aa_, dim_ij_ab_, dim_ij_bb_;
  RefSCDimension dim_ab_aa_, dim_ab_ab_, dim_ab_bb_;
#endif
  RefSCDimension dim_ij_s_, dim_ij_t_;

  bool gbc_;
  bool ebc_;
  Ref<LinearR12::CorrelationFactor> corrfactor_;
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

  /// Fock-weighted occupied space |i_f> = f_i^R |R>, where R is a function in RI-BS
  Ref<MOIndexSpace> focc_space_;
  /// Form Fock-weighted occupied space
  void form_focc_space_();
  /// Fock-weighted active occupied space |i_f> = f_i^R |R>, where R is a function in RI-BS
  Ref<MOIndexSpace> factocc_space_;
  /// Form Fock-weighted active occupied space
  void form_factocc_space_();
  /// Form space of auxiliary virtuals
  void form_canonvir_space_();

  /// Initialize standard transforms
  void init_tforms_();
  /// Set intermediates to zero + add the "diagonal" contributions
  void init_intermeds_();
  /// When F12=R12 number of simplifications occur so a specialized code is provided
  void init_intermeds_r12_();
  /// When F12 != R12 the following code is used
  void init_intermeds_g12_(SpinCase2 S);
  /// Compute r^2 contribution to X using compute_r2_()
  void r2_contrib_to_X_new_();
  /// Compute <space1 space2|r_{12}^2|space3 space4> matrix
  RefSCMatrix compute_r2_(const Ref<MOIndexSpace>& space1,
                          const Ref<MOIndexSpace>& space2,
                          const Ref<MOIndexSpace>& space3,
                          const Ref<MOIndexSpace>& space4);
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

  /// Checkpoint the top-level molecular energy
  void checkpoint_() const;

  /** Compute the number tasks which have access to the integrals.
      map_to_twi is a vector<int> each element of which will hold the
      number of tasks with access to the integrals of lower rank than that task
      (or -1 if the task doesn't have access to the integrals) */
  const int tasks_with_ints_(const Ref<R12IntsAcc> ints_acc, vector<int>& map_to_twi);

  /// New, all-encompassing version
  void contrib_to_VXB_a_new_(const Ref<MOIndexSpace>& ispace,
                             const Ref<MOIndexSpace>& xspace,
                             const Ref<MOIndexSpace>& jspace,
                             const Ref<MOIndexSpace>& yspace,
                             SpinCase2 spincase,
                             const Ref<LinearR12::TwoParticleContraction>& tpcontract);

  /// Compute MP2 pair energies of spin case S
  void compute_mp2_pair_energies_(SpinCase2 S);
  /// Compute MP2 pair energies of spin case S using < space1 space3|| space2 space4> integrals
  void compute_mp2_pair_energies_(RefSCVector& emp2pair,
                                  SpinCase2 S,
                                  const Ref<MOIndexSpace>& space1,
                                  const Ref<MOIndexSpace>& space2,
                                  const Ref<MOIndexSpace>& space3,
                                  const Ref<MOIndexSpace>& space4,
                                  const Ref<TwoBodyMOIntsTransform>& transform);
  /// Compute VXB intermeds when VBS is not the same as OBS
  void contrib_to_VXB_gebc_vbsneqobs_();

  /// Compute A using the "simple" formula obtained using direct substitution alpha'->a'
  void compute_A_simple_();

  /// Compute A using the standard commutator approach
  void compute_A_via_commutator_();

  /** Compute T2 amplitude in basis <space1, space3 | space2, space4>.
      AlphaBeta amplitudes are computed.
      If tform is not given (it should be!), this function will construct a generic
      transform. */
  void compute_T2_(RefSCMatrix& T2,
                   const Ref<MOIndexSpace>& space1,
                   const Ref<MOIndexSpace>& space2,
                   const Ref<MOIndexSpace>& space3,
                   const Ref<MOIndexSpace>& space4,
                   bool antisymmetrize,
                   const Ref<TwoBodyMOIntsTransform>& tform = 0);
  /** Compute F12 integrals in basis <space1, space3 | f12 | space2, space4>.
      Bra (rows) are blocked by correlation function index.
      AlphaBeta amplitudes are computed.
      If tform is not given (it should be!), this function will construct a generic
      transform. */
  void compute_F12_(RefSCMatrix& F12,
                   const Ref<MOIndexSpace>& space1,
                   const Ref<MOIndexSpace>& space2,
                   const Ref<MOIndexSpace>& space3,
                   const Ref<MOIndexSpace>& space4,
                   const std::vector< Ref<TwoBodyMOIntsTransform> >& transforms = std::vector< Ref<TwoBodyMOIntsTransform> >());

  /** Compute A intermediate using "direct" formula in basis <space1, space3 | f12 | space2, space4>.
      Bra (rows) are blocked by correlation function index.
      AlphaBeta amplitudes are computed.
      If tform is not given (it should be!), this function will construct a generic
      transform. */
  void compute_A_direct_(RefSCMatrix& A,
                         const Ref<MOIndexSpace>& space1,
                         const Ref<MOIndexSpace>& space2,
                         const Ref<MOIndexSpace>& space3,
                         const Ref<MOIndexSpace>& space4,
                         const Ref<MOIndexSpace>& rispace2,
                         const Ref<MOIndexSpace>& rispace4,
                         const std::vector< Ref<TwoBodyMOIntsTransform> >& transforms =
                               std::vector< Ref<TwoBodyMOIntsTransform> >() );

  /** compute_tbint_tensor computes a 2-body tensor T using integrals of type tbint_type.
      Computed tensor T is added to its previous contents.
      Class DataProcess defines a static function 'double I2T()' which processes the integrals.
      Set CorrFactorInBra to true if bra of target tensor depends on correlation function index.
   */
  template <typename DataProcess, bool CorrFactorInBra, bool CorrFactorInKet>
    void compute_tbint_tensor(RefSCMatrix& T,
                              int tbint_type,
                              const Ref<MOIndexSpace>& space1,
                              const Ref<MOIndexSpace>& space2,
                              const Ref<MOIndexSpace>& space3,
                              const Ref<MOIndexSpace>& space4,
                              bool antisymmetrize,
                              const std::vector< Ref<TwoBodyMOIntsTransform> >& tforms = 
                                std::vector< Ref<TwoBodyMOIntsTransform> >(),
                              const std::vector< Ref<TwoBodyIntDescr> >& tbintdescrs =
                                std::vector< Ref<TwoBodyIntDescr> >());
  
#if 0
  /// Compute R "intermediate" (r12 integrals in occ-pair/vir-pair basis)
  void compute_R_();
#endif

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

  /// Compute dual-basis MP1 energy (contribution from singles to HF energy)
  void compute_dualEmp1_();

#if 0
  /// This function computes T2 amplitudes
  void compute_T2_vbsneqobs_();
#endif

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
  R12IntEval(const Ref<R12IntEvalInfo>& info, const Ref<LinearR12::CorrelationFactor>& corrfactor,
             bool gbc = true, bool ebc = true,
             LinearR12::ABSMethod abs_method = LinearR12::ABS_CABSPlus,
             LinearR12::StandardApproximation stdapprox = LinearR12::StdApprox_Ap,
             bool follow_ks_ebcfree = false);
  /*R12IntEval(const Ref<R12IntEvalInfo>& info, bool gbc = true, bool ebc = true,
             LinearR12::ABSMethod abs_method = LinearR12::ABS_CABSPlus,
             LinearR12::StandardApproximation stdapprox = LinearR12::StdApprox_Ap,
             bool follow_ks_ebcfree = false);*/
  ~R12IntEval();

  void save_data_state(StateOut&);
  virtual void obsolete();

  void include_mp1(bool include_mp1);
  void set_debug(int debug);
  void set_dynamic(bool dynamic);
  void set_print_percent(double print_percent);
  void set_memory(size_t nbytes);

  const Ref<LinearR12::CorrelationFactor>& corrfactor() const { return corrfactor_; }
  bool spin_polarized() const { return r12info_->refinfo()->ref()->spin_polarized(); }
  bool gbc() const { return gbc_; }
  bool ebc() const { return ebc_; }
  LinearR12::StandardApproximation stdapprox() const { return stdapprox_; }
  bool follow_ks_ebcfree() const { return follow_ks_ebcfree_; }

  const Ref<R12IntEvalInfo>& r12info() const;
#if 0
  RefSCDimension dim_oo_aa() const;
  RefSCDimension dim_oo_ab() const;
  RefSCDimension dim_oo_bb() const;
  RefSCDimension dim_vv_aa() const;
  RefSCDimension dim_vv_ab() const;
  RefSCDimension dim_vv_bb() const;
#endif
  RefSCDimension dim_oo_s() const;
  RefSCDimension dim_oo_t() const;
  /// Dimension for active-occ/active-occ pairs of spin case S
  RefSCDimension dim_oo(SpinCase2 S) const;
  /// Dimension for active-vir/active-vir pairs of spin case S
  RefSCDimension dim_vv(SpinCase2 S) const;
  /// Dimension for geminal functions of spin case S
  RefSCDimension dim_f12(SpinCase2 S) const;

  /// This function causes the intermediate matrices to be computed.
  virtual void compute();

#if 0
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
  /// Returns beta-beta block of the V intermediate matrix.
  RefSCMatrix V_bb();
  /// Returns beta-beta block of the X intermediate matrix.
  RefSCMatrix X_bb();
  /// Returns beta-beta block of the B intermediate matrix.
  RefSymmSCMatrix B_bb();
  /// Returns beta-beta block of the A intermediate matrix. Returns 0 if EBC is assumed.
  RefSCMatrix A_bb();
  /// Returns beta-beta block of the MP2 T2 matrix. Returns 0 if EBC is assumed.
  RefSCMatrix T2_bb();
  /// Returns alpha-alpha MP2 pair energies.
  RefSCVector emp2_aa();
  /// Returns alpha-beta MP2 pair energies.
  RefSCVector emp2_ab();
  /// Returns beta-beta MP2 pair energies.
  RefSCVector emp2_bb();
#endif
  /// Returns amplitudes of pair correlation functions
  //Ref<F12Amplitudes> amps();
  
  /// Returns S block of intermediate V
  const RefSCMatrix& V(SpinCase2 S);
  /// Returns S block of intermediate X
  const RefSCMatrix& X(SpinCase2 S);
  /// Returns S block of intermediate B
  RefSymmSCMatrix B(SpinCase2 S);
  /// Returns S block of intermediate A
  const RefSCMatrix& A(SpinCase2 S);
  /// Returns S block of intermediate A computed using commutator method
  const RefSCMatrix& Ac(SpinCase2 S);
  /// Returns S block of intermediate T2
  const RefSCMatrix& T2(SpinCase2 S);
  /// Returns S block of intermediate F12
  const RefSCMatrix& F12(SpinCase2 S);
  /// Returns alpha-alpha MP2 pair energies
  const RefSCVector& emp2(SpinCase2 S);
  /// Returns the eigenvalues of spin case S
  const RefDiagSCMatrix& evals(SpinCase1 S) const;
  
  /// Returns the eigenvalues for the closed-shell case
  RefDiagSCMatrix evals() const;
  /// Returns the alpha eigenvalues
  RefDiagSCMatrix evals_a() const;
  /// Returns the beta eigenvalues
  RefDiagSCMatrix evals_b() const;
  
  // Returns the number of unique combinations of 2 spin cases
  int nspincases2() const { return ::sc::nspincases2(spin_polarized()); }
  /// Returns the act occ space for spin case S
  const Ref<MOIndexSpace>& occ_act(SpinCase1 S) const;
  /// Returns the act vir space for spin case S
  const Ref<MOIndexSpace>& vir_act(SpinCase1 S) const;
  
  /** Returns an already created transform.
      If the transform is not found then throw TransformNotFound */
  Ref<TwoBodyMOIntsTransform> get_tform_(const std::string&);
  /// Generates canonical id for transform. f12 is the index of the correlation function
  std::string transform_label(const Ref<MOIndexSpace>& space1,
                              const Ref<MOIndexSpace>& space2,
                              const Ref<MOIndexSpace>& space3,
                              const Ref<MOIndexSpace>& space4,
                              unsigned int f12) const;
  /// version of transform_label() applicable when left and right correlation factors differ
  std::string transform_label(const Ref<MOIndexSpace>& space1,
                              const Ref<MOIndexSpace>& space2,
                              const Ref<MOIndexSpace>& space3,
                              const Ref<MOIndexSpace>& space4,
                              unsigned int f12_left,
                              unsigned int f12_right) const;
  
};

class TransformNotFound: public ProgrammingError {
  public:
  TransformNotFound(const char *description=0, const char *file=0, int line=0) : ProgrammingError(description,file,line) {}
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


