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
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/mbptr12/spin.h>

#ifndef _chemistry_qc_mbptr12_r12inteval_h
#define _chemistry_qc_mbptr12_r12inteval_h

namespace sc {
  
  class TwoBodyMOIntsTransform;
  class R12IntsAcc;
  class F12Amplitudes;

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
  RefSCMatrix BB_[NSpinCases2];  // The difference between B intermediate of approximation B and A'
  RefSCMatrix A_[NSpinCases2];
  RefSCMatrix Ac_[NSpinCases2];
#if 0
  RefSCMatrix T2_[NSpinCases2];
  RefSCMatrix F12_[NSpinCases2];
#endif
  RefSCVector emp2pair_[NSpinCases2];
  RefSCDimension dim_oo_[NSpinCases2];
  RefSCDimension dim_vv_[NSpinCases2];
  RefSCDimension dim_f12_[NSpinCases2];
  
  Ref<F12Amplitudes> Amps_;  // First-order amplitudes of various contributions to the pair functions
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

  /// Returns the number of unique spin cases
  int nspincases1() const { return ::sc::nspincases1(spin_polarized()); }
  /// Returns the number of unique combinations of 2 spin cases
  int nspincases2() const { return ::sc::nspincases2(spin_polarized()); }
  /// "Spin-adapt" MO space id and name
  void spinadapt_mospace_labels(SpinCase1 spin, std::string& id, std::string& name) const;
  
  /// Fock-weighted occupied space |i_f> = f_i^R |R>, where R is a function in RI-BS
  Ref<MOIndexSpace> focc_space_[NSpinCases1];
  /// Fock-weighted active occupied space |i_f> = f_i^R |R>, where R is a function in RI-BS
  Ref<MOIndexSpace> factocc_space_[NSpinCases1];
  /// Exchange-weighted active occupied space |i_K> = K_i^R |R>, where R is a function in RI-BS
  Ref<MOIndexSpace> kactocc_space_[NSpinCases1];
  /// Fock-weighted active occupied space |a_f> = f_a^R |R>, where R is a function in RI-BS
  Ref<MOIndexSpace> factvir_space_[NSpinCases1];
  /// Exchange-weighted active virtual space |a_K> = K_a^R |R>, where R is a function in RI-BS
  Ref<MOIndexSpace> kactvir_space_[NSpinCases1];
  /// Exchange-weighted RI space |a'_K> = K_a'^A |A>, where A is a function in ABS
  Ref<MOIndexSpace> kribs_space_[NSpinCases1];
  /// Exchange-weighted (through OBS) active occupied space |i_K> = K_i^P |P>, where P is a function in OBS
  Ref<MOIndexSpace> kactocc_obs_space_[NSpinCases1];
  /// Exchange-weighted (through OBS) virtual space |a_K> = K_a^p |p>, where p is a function in OBS
  Ref<MOIndexSpace> kvir_obs_space_[NSpinCases1];
  /// Exchange-weighted (through RIBS) virtual space |a_K> = K_a^P |P>, where P is a function in RIBS
  Ref<MOIndexSpace> kvir_ribs_space_[NSpinCases1];
  /// Compute factocc and kactocc spaces, if needed
  void form_focc_act(SpinCase1 spin);
  /// Compute factvir and kactvir spaces, if needed
  void form_fvir_act(SpinCase1 spin);
  /// Compute kribs space, if needed
  void form_fribs(SpinCase1 spin);
  /// Compute kactocc_obs space, if needed
  void form_focc_act_obs(SpinCase1 spin);
  /// Compute kvir_obs space, if needed
  void form_fvir_obs(SpinCase1 spin);
  /// Compute kvir_ribs space, if needed
  void form_fvir_ribs(SpinCase1 spin);
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
  /** Compute the Fock matrix between bra_ and ket_ spaces of spin S.
      scale_J and scale_K are used to scale Coulomb
      and exchange contributions, T12IntEval::occ() is used for the occupied spaces.
      */
  RefSCMatrix fock_(const Ref<MOIndexSpace>& bra_space,
                    const Ref<MOIndexSpace>& ket_space,
                    SpinCase1 S = Alpha,
                    double scale_J = 1.0, double scale_K = 1.0);
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

#if 0
  /// Compute MP2 pair energies of spin case S
  void compute_mp2_pair_energies_(SpinCase2 S);
#endif
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
                         const Ref<MOIndexSpace>& rispace4);
  
  /** Compute A intermediate using "commutator" formula in basis <space1, space3 | f12 | space2, space4>.
      Bra (rows) are blocked by correlation function index.
      AlphaBeta amplitudes are computed.
      If tform is not given (it should be!), this function will construct a generic
      transform. */
  void compute_A_viacomm_(RefSCMatrix& A,
                          const Ref<MOIndexSpace>& space1,
                          const Ref<MOIndexSpace>& space2,
                          const Ref<MOIndexSpace>& space3,
                          const Ref<MOIndexSpace>& space4,
                          const std::vector< Ref<TwoBodyMOIntsTransform> >& tforms);
  
  /** compute_tbint_tensor computes a 2-body tensor T using integrals of type tbint_type.
      Computed tensor T is added to its previous contents.
      Class DataProcess defines a static function 'double I2T()' which processes the integrals.
      Set CorrFactorInBra to true if bra of target tensor depends on correlation function index.
      
      Of course, this ugliness should become constructor/member function of ManyBodyOperator
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
  
  /** contract_tbint_tensor computes a 2-body tensor T as a sum over mn : <ij|Tbra|mn> * <kl|Tket|mn>^t
      bra and ket integrals come from tforms_bra and tforms_ket.
      Computed tensor T is added to its previous contents.
      Class DataProcess_XXX defines a static function 'double I2T()' which processes the integrals.
      Set CorrFactorInBra to true if bra of target tensor depends on correlation function index.
      
      Semantics: 1) ranks of internal spaces must match, although the spaces don't have to be the same;
      2) if antisymmetrize is true, external (bra and ket) particle spaces
      should be strongly (identity) or weakly(rank) equivalent. If internal spaces do not match
      but antisymmetrization is requested
      
      Of course, this ugliness should become function/operator on 2 ManyBodyOperators
   */
  template <typename DataProcessBra,
            typename DataProcessKet,
            typename DataProcessBraKet,
            bool CorrFactorInBra,
            bool CorrFactorInKet,
            bool CorrFactorInInt>
    void contract_tbint_tensor(
           RefSCMatrix& T,
           unsigned int tbint_type_bra,
           unsigned int tbint_type_ket,
           const Ref<MOIndexSpace>& space1_bra,
           const Ref<MOIndexSpace>& space2_bra,
           const Ref<MOIndexSpace>& space1_intb,
           const Ref<MOIndexSpace>& space2_intb,
           const Ref<MOIndexSpace>& space1_ket,
           const Ref<MOIndexSpace>& space2_ket,
           const Ref<MOIndexSpace>& space1_intk,
           const Ref<MOIndexSpace>& space2_intk,
           const Ref<LinearR12::TwoParticleContraction>& tpcontract,
           bool antisymmetrize,
           const std::vector< Ref<TwoBodyMOIntsTransform> >& tforms_bra = 
             std::vector< Ref<TwoBodyMOIntsTransform> >(),
           const std::vector< Ref<TwoBodyMOIntsTransform> >& tforms_ket = 
             std::vector< Ref<TwoBodyMOIntsTransform> >(),
           const std::vector< Ref<TwoBodyIntDescr> >& tbintdescrs_bra =
             std::vector< Ref<TwoBodyIntDescr> >(),
           const std::vector< Ref<TwoBodyIntDescr> >& tbintdescrs_ket =
             std::vector< Ref<TwoBodyIntDescr> >());
  
  /** Compute X intermediate (F12 F12) in basis <bra1 bra2 | ket1 ket2>. sc2 specifies the spin case
      of particles 1 and 2. Resulting X is symmetric w.r.t bra-ket permutation, but will not be
      symmetric w.r.t. permutation of particles 1 and 2 if bra1 != bra2 || ket1 != ket2.
      If X is null, then allocate, otherwise check dimensions and accumulate result into X.
      
      Semantics: 1) sc2 == AlphaAlpha or BetaBeta means that X will be antisymmetric w.r.t
      permutations bra1<->bra2 or ket1<->ket2, maybe artificially so. So sc2 == AlphaBeta
      means "particles are not equivalent or different spin", whereas AlphaAlpha means "act like
      particles are equivalent".
  */
  void compute_X_(RefSCMatrix& X,
                  SpinCase2 sc2,
                  const Ref<MOIndexSpace>& bra1,
                  const Ref<MOIndexSpace>& bra2,
                  const Ref<MOIndexSpace>& ket1,
                  const Ref<MOIndexSpace>& ket2);
  /** Compute contraction <bra1 bra2|F12|intb1 int2> <intb1|x|intk1> <intk1 int2|F12|ket1 ket2> +
      <bra1 bra2|F12|int1 intb2> <intb2|x|intk2> <int1 intk2|F12|ket1 ket2>
      sc2 specifies the spin case of particles 1 and 2.
      Resulting X is symmetric w.r.t bra-ket permutation, but may not be
      symmetric w.r.t. permutation of particles 1 and 2.
      If FxF is null, then allocate, otherwise check dimensions and accumulate result into FxF.
      intkx1 is the "combined" space |intb1> <intb1|x|intk1>, etc.
  */
  void compute_FxF_(RefSCMatrix& FxF,
                    SpinCase2 sc2,
                    const Ref<MOIndexSpace>& bra1,
                    const Ref<MOIndexSpace>& bra2,
                    const Ref<MOIndexSpace>& ket1,
                    const Ref<MOIndexSpace>& ket2,
                    const Ref<MOIndexSpace>& int1,
                    const Ref<MOIndexSpace>& int2,
                    const Ref<MOIndexSpace>& intk1,
                    const Ref<MOIndexSpace>& intk2,
                    const Ref<MOIndexSpace>& intkx1,
                    const Ref<MOIndexSpace>& intkx2
                   );
  
  /// Compute A*T2 contribution to V (needed if EBC is not assumed)
  void AT2_contrib_to_V_();

  /// Compute -2*A*R contribution to B (needed if EBC is not assumed)
  void AF12_contrib_to_B_();

  /** Compute contributions to B that vanishe under GBC */
  void compute_B_gbc_();

  /** Compute contributions to B which vanish in standard approximation A' */
  void compute_BB_();

  /// Compute dual-basis MP1 energy (contribution from singles to HF energy)
  void compute_dualEmp1_();

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

  /// Returns amplitudes of pair correlation functions
  Ref<F12Amplitudes> amps();
  
  /// Returns S block of intermediate V
  const RefSCMatrix& V(SpinCase2 S);
  /// Returns S block of intermediate X
  const RefSCMatrix& X(SpinCase2 S);
  /// Returns S block of intermediate B
  RefSymmSCMatrix B(SpinCase2 S);
  /// Returns S block of intermediate BB
  RefSymmSCMatrix BB(SpinCase2 S);
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
  
  /// Returns the act occ space for spin case S
  const Ref<MOIndexSpace>& occ_act(SpinCase1 S) const;
  /// Returns the occ space for spin case S
  const Ref<MOIndexSpace>& occ(SpinCase1 S) const;
  /// Returns the act vir space for spin case S
  const Ref<MOIndexSpace>& vir_act(SpinCase1 S) const;
  /// Returns the vir space for spin case S
  const Ref<MOIndexSpace>& vir(SpinCase1 S) const;
  /// Form Fock-weighted occupied space for spin case S
  const Ref<MOIndexSpace>& focc(SpinCase1 S);
  /// Form Fock-weighted active occupied space for spin case S
  const Ref<MOIndexSpace>& focc_act(SpinCase1 S);
  /// Form exchange-weighted active occupied space for spin case S
  const Ref<MOIndexSpace>& kocc_act(SpinCase1 S);
  /// Form exchange-weighted (through OBS) active occupied space for spin case S
  const Ref<MOIndexSpace>& kocc_act_obs(SpinCase1 S);
  /// Form Fock-weighted active virtual space for spin case S
  const Ref<MOIndexSpace>& fvir_act(SpinCase1 S);
  /// Form exchange-weighted active virtual space for spin case S
  const Ref<MOIndexSpace>& kvir_act(SpinCase1 S);
  /// Form exchange-weighted RI space for spin case S
  const Ref<MOIndexSpace>& kribs(SpinCase1 S);
  /// Form exchange-weighted (through OBS) virtual space for spin case S
  const Ref<MOIndexSpace>& kvir_obs(SpinCase1 S);
  /// Form exchange-weighted (through RIBS) virtual space for spin case S
  const Ref<MOIndexSpace>& kvir_ribs(SpinCase1 S);
  
  /** Returns an already created transform.
      If the transform is not found then throw TransformNotFound */
  Ref<TwoBodyMOIntsTransform> get_tform_(const std::string&) const;
  /** Map transform T to label */
  void add_tform(const std::string& label,
                 const Ref<TwoBodyMOIntsTransform>& T);
  /// Generates canonical id for transform. no correlation function included
  std::string transform_label(const Ref<MOIndexSpace>& space1,
                              const Ref<MOIndexSpace>& space2,
                              const Ref<MOIndexSpace>& space3,
                              const Ref<MOIndexSpace>& space4) const;
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
                   bool antisymmetrize,
                   const std::vector< Ref<TwoBodyMOIntsTransform> >& transforms
                     = std::vector< Ref<TwoBodyMOIntsTransform> >(),
                   const std::vector< Ref<TwoBodyIntDescr> >& descrvec
                     = std::vector< Ref<TwoBodyIntDescr> >() );

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


