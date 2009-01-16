//
// r12int_eval.h
//
// Copyright (C) 2004 Edward Valeev
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
  // change to false to use the old fock builder
  static const bool USE_FOCKBUILD = true;

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

  RefSCVector emp2pair_[NSpinCases2];
  RefSCDimension dim_oo_[NSpinCases2];
  RefSCDimension dim_vv_[NSpinCases2];
  /// a is any index in a given basis
  RefSCDimension dim_aa_[NSpinCases2];
  RefSCDimension dim_f12_[NSpinCases2];
  RefSCDimension dim_xy_[NSpinCases2];

  Ref<F12Amplitudes> Amps_;  // First-order amplitudes of various contributions to the pair functions
  RefSCDimension dim_ij_s_, dim_ij_t_;
  double emp2_singles_;
  int debug_;

  // Map to TwoBodyMOIntsTransform objects that have been computed previously
  typedef std::map<std::string, Ref<TwoBodyMOIntsTransform> > TformMap;
  TformMap tform_map_;

  /// "Spin-adapt" MO space id and name
  void spinadapt_mospace_labels(SpinCase1 spin, std::string& id, std::string& name) const;

  /// Form space of auxiliary virtuals
  void form_canonvir_space_();

  /// This is the new way to generate needed spaces
  /// generates fock, h+J, or K weighted spaces
  void f_bra_ket(SpinCase1 spin,
		 bool make_F,
		 bool make_hJ,
		 bool make_K,
		 Ref<OrbitalSpace>& F,
		 Ref<OrbitalSpace>& hJ,
		 Ref<OrbitalSpace>& K,
		 const Ref<OrbitalSpace>& extspace,
		 const Ref<OrbitalSpace>& intspace
      );
  Ref<OrbitalSpace> hj_i_p_[NSpinCases1];
  Ref<OrbitalSpace> hj_i_A_[NSpinCases1];
  Ref<OrbitalSpace> hj_i_P_[NSpinCases1];
  Ref<OrbitalSpace> hj_i_m_[NSpinCases1];
  Ref<OrbitalSpace> hj_i_a_[NSpinCases1];
  Ref<OrbitalSpace> hj_p_p_[NSpinCases1];
  Ref<OrbitalSpace> hj_p_A_[NSpinCases1];
  Ref<OrbitalSpace> hj_p_P_[NSpinCases1];
  Ref<OrbitalSpace> hj_P_P_[NSpinCases1];
  Ref<OrbitalSpace> hj_p_m_[NSpinCases1];
  Ref<OrbitalSpace> hj_p_a_[NSpinCases1];
  Ref<OrbitalSpace> K_i_p_[NSpinCases1];
  Ref<OrbitalSpace> K_i_m_[NSpinCases1];
  Ref<OrbitalSpace> K_i_a_[NSpinCases1];
  Ref<OrbitalSpace> K_i_A_[NSpinCases1];
  Ref<OrbitalSpace> K_i_P_[NSpinCases1];
  Ref<OrbitalSpace> K_m_a_[NSpinCases1];
  Ref<OrbitalSpace> K_a_a_[NSpinCases1];
  Ref<OrbitalSpace> K_a_p_[NSpinCases1];
  Ref<OrbitalSpace> K_a_P_[NSpinCases1];
  Ref<OrbitalSpace> K_p_p_[NSpinCases1];
  Ref<OrbitalSpace> K_p_m_[NSpinCases1];
  Ref<OrbitalSpace> K_p_a_[NSpinCases1];
  Ref<OrbitalSpace> K_p_A_[NSpinCases1];
  Ref<OrbitalSpace> K_p_P_[NSpinCases1];
  Ref<OrbitalSpace> K_A_P_[NSpinCases1];
  Ref<OrbitalSpace> K_P_P_[NSpinCases1];
  Ref<OrbitalSpace> F_P_P_[NSpinCases1];
  Ref<OrbitalSpace> F_p_A_[NSpinCases1];
  Ref<OrbitalSpace> F_p_p_[NSpinCases1];
  Ref<OrbitalSpace> F_p_m_[NSpinCases1];
  Ref<OrbitalSpace> F_p_a_[NSpinCases1];
  Ref<OrbitalSpace> F_m_m_[NSpinCases1];
  Ref<OrbitalSpace> F_m_a_[NSpinCases1];
  Ref<OrbitalSpace> F_m_P_[NSpinCases1];
  Ref<OrbitalSpace> F_m_A_[NSpinCases1];
  Ref<OrbitalSpace> F_i_A_[NSpinCases1];
  Ref<OrbitalSpace> F_i_m_[NSpinCases1];
  Ref<OrbitalSpace> F_i_a_[NSpinCases1];
  Ref<OrbitalSpace> F_i_p_[NSpinCases1];
  Ref<OrbitalSpace> F_a_a_[NSpinCases1];
  Ref<OrbitalSpace> F_a_A_[NSpinCases1];

  /// Initialize standard transforms
  void init_tforms_();
  /// Set intermediates to zero + add the "diagonal" contributions
  void init_intermeds_();
  /// When F12=R12 number of simplifications occur so a specialized code is provided
  void init_intermeds_r12_();
  /// When F12 != R12 the following code is used
  void init_intermeds_g12_();
  /// Compute r^2 contribution to X using compute_r2_()
  void r2_contrib_to_X_new_();
  /// Compute <space1 space2|space3 space4> matrix
  RefSCMatrix compute_I_(const Ref<OrbitalSpace>& space1,
			 const Ref<OrbitalSpace>& space2,
			 const Ref<OrbitalSpace>& space3,
			 const Ref<OrbitalSpace>& space4);
  /// Compute <space1 space2|r_{12}^2|space3 space4> matrix
  RefSCMatrix compute_r2_(const Ref<OrbitalSpace>& space1,
                          const Ref<OrbitalSpace>& space2,
                          const Ref<OrbitalSpace>& space3,
                          const Ref<OrbitalSpace>& space4);
  /** Compute the relativistic hcore Hamiltonian using DKH2 and substract
      T, V and the mass-velocity term. Based on r12int_eval::fock. In file fock.cc
  */
  RefSCMatrix Delta_DKH_(const Ref<OrbitalSpace>& bra_space,
                         const Ref<OrbitalSpace>& ket_space,
                         SpinCase1 S = Alpha);
  // Computes T, V and the mass-velocity term in the momentum basis.
  // It's a modified version of Wavefunction::core_hamiltonian_dk
  RefSymmSCMatrix hcore_plus_massvelocity_(const Ref<GaussianBasisSet> &bas,
                                           const Ref<GaussianBasisSet> &p_bas,
                                           bool include_T = true,
                                           bool include_V = true,
                                           bool include_MV = true);
  // Computes the skalar Pauli-Hamiltonian (T + V + mass_velocity + Darwin),
  // with the mass-velocity term evaluated in the momentum basis.
  RefSymmSCMatrix pauli(const Ref<GaussianBasisSet> &bas,
                        const Ref<GaussianBasisSet> &pbas = 0,
                        const bool momentum=false);
  RefSymmSCMatrix pauli_realspace(const Ref<GaussianBasisSet> &bas);
  RefSymmSCMatrix pauli_momentumspace(const Ref<GaussianBasisSet> &bas,
                        const Ref<GaussianBasisSet> &p_bas);
  /// Compute the coulomb matrix between 2 spaces
  RefSCMatrix coulomb_(const Ref<OrbitalSpace>& occ_space, const Ref<OrbitalSpace>& bra_space,
                       const Ref<OrbitalSpace>& ket_space);
  /// Compute the exchange matrix between 2 spaces
  RefSCMatrix exchange_(const Ref<OrbitalSpace>& occ_space, const Ref<OrbitalSpace>& bra_space,
                        const Ref<OrbitalSpace>& ket_space);

  /// Checkpoint the top-level molecular energy
  void checkpoint_() const;

  /// New version which uses tensor contract functions
  void contrib_to_VXB_a_();
  /// New version which uses tensor contract functions
  void contrib_to_VXB_a_vbsneqobs_();

  /// Compute MP2 pair energies of spin case S using < space1 space3|| space2 space4> integrals
  void compute_mp2_pair_energies_(RefSCVector& emp2pair,
                                  SpinCase2 S,
                                  const Ref<OrbitalSpace>& space1,
                                  const Ref<OrbitalSpace>& space2,
                                  const Ref<OrbitalSpace>& space3,
                                  const Ref<OrbitalSpace>& space4,
                                  const std::string& tform_key);

  /** Compute A intermediate using "direct" formula in basis <space1, space3 | f12 | space2, space4>.
      Bra (rows) are blocked by correlation function index.
      AlphaBeta amplitudes are computed.
      If tform is not given (it should be!), this function will construct a generic
      transform. */
  void compute_A_direct_(RefSCMatrix& A,
                         const Ref<OrbitalSpace>& space1,
                         const Ref<OrbitalSpace>& space2,
                         const Ref<OrbitalSpace>& space3,
                         const Ref<OrbitalSpace>& space4,
                         const Ref<OrbitalSpace>& rispace2,
                         const Ref<OrbitalSpace>& rispace4,
                         bool antisymmetrize);

  /** compute_tbint_tensor computes a 2-body tensor T using integrals of type tbint_type.
      Computed tensor T is added to its previous contents.
      Class DataProcess defines a static function 'double I2T()' which processes the integrals.
      Set CorrFactorInBra to true if bra of target tensor depends on correlation function index.

      Of course, this ugliness should become constructor/member function of ManyBodyOperator
   */
  template <typename DataProcess, bool CorrFactorInBra, bool CorrFactorInKet>
    void compute_tbint_tensor(RefSCMatrix& T,
                              TwoBodyInt::tbint_type tbint_type,
                              const Ref<OrbitalSpace>& space1,
                              const Ref<OrbitalSpace>& space2,
                              const Ref<OrbitalSpace>& space3,
                              const Ref<OrbitalSpace>& space4,
                              bool antisymmetrize,
                              const std::vector<std::string>& tform_keys);

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
           TwoBodyInt::tbint_type tbint_type_bra,
           TwoBodyInt::tbint_type tbint_type_ket,
           const Ref<OrbitalSpace>& space1_bra,
           const Ref<OrbitalSpace>& space2_bra,
           const Ref<OrbitalSpace>& space1_intb,
           const Ref<OrbitalSpace>& space2_intb,
           const Ref<OrbitalSpace>& space1_ket,
           const Ref<OrbitalSpace>& space2_ket,
           const Ref<OrbitalSpace>& space1_intk,
           const Ref<OrbitalSpace>& space2_intk,
           const Ref<LinearR12::TwoParticleContraction>& tpcontract,
           bool antisymmetrize,
           const std::vector<std::string>& tformkeys_bra,
           const std::vector<std::string>& tformkeys_ket);

  /** Compute X intermediate (F12 F12) in basis <bra1 bra2 | ket1 ket2>. sc2 specifies the spin case
      of particles 1 and 2. Resulting X is symmetric w.r.t bra-ket permutation, but will not be
      symmetric w.r.t. permutation of particles 1 and 2 if bra1 != bra2 || ket1 != ket2.
      If X is null, then allocate, otherwise check dimensions and accumulate result into X.
      Setting F2_only to true will only leave F12^2 term.

      Semantics: 1) sc2 == AlphaAlpha or BetaBeta means that X will be antisymmetric w.r.t
      permutations bra1<->bra2 or ket1<->ket2, maybe artificially so. So sc2 == AlphaBeta
      means "particles are not equivalent or different spin", whereas AlphaAlpha means "act like
      particles are equivalent".
  */
  void compute_X_(RefSCMatrix& X,
                  SpinCase2 sc2,
                  const Ref<OrbitalSpace>& bra1,
                  const Ref<OrbitalSpace>& bra2,
                  const Ref<OrbitalSpace>& ket1,
                  const Ref<OrbitalSpace>& ket2,
                  bool F2_only = false);
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
                    const Ref<OrbitalSpace>& bra1,
                    const Ref<OrbitalSpace>& bra2,
                    const Ref<OrbitalSpace>& ket1,
                    const Ref<OrbitalSpace>& ket2,
                    const Ref<OrbitalSpace>& int1,
                    const Ref<OrbitalSpace>& int2,
                    const Ref<OrbitalSpace>& intk1,
                    const Ref<OrbitalSpace>& intk2,
                    const Ref<OrbitalSpace>& intkx1,
                    const Ref<OrbitalSpace>& intkx2
                   );

  /// Compute A*T2 contribution to V (needed if EBC is not assumed)
  void AT2_contrib_to_V_();

  /// Compute -2*A*R contribution to B (needed if EBC is not assumed)
  void AF12_contrib_to_B_();

  /** Compute contributions to B that vanish under GBC */
  void compute_B_gbc_();

  /** Compute contributions to B that vanish under BC */
  void compute_B_bc_();

  /** Compute contributions to B which appear in standard approximation B and not in A' */
  void compute_BB_();

  /** Compute B using standard approximation C */
  void compute_BC_();

  /** Compute B using standard approximation A'' -- exchange is dropped completely! */
  void compute_BApp_();

  /** Compute the mass-velocity contributions to B that appear when DKH-based R12 calculations are requested */
  void compute_B_DKH_();

  /// Compute singles contribution to the MP2 energy
  void compute_singles_emp2_();

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
  /** Constructs R12IntEval. */
  R12IntEval(const Ref<R12IntEvalInfo>& info);
  /*R12IntEval(const Ref<R12IntEvalInfo>& info, bool gbc = true, bool ebc = true,
             LinearR12::ABSMethod abs_method = LinearR12::ABS_CABSPlus,
             LinearR12::StandardApproximation stdapprox = LinearR12::StdApprox_Ap);*/
  ~R12IntEval();

  void save_data_state(StateOut&);
  virtual void obsolete();

  void set_debug(int debug);
  void set_dynamic(bool dynamic);
  void set_print_percent(double print_percent);
  void set_memory(size_t nbytes);

  /// Indicates whether Douglas-Kroll Hamiltonian is used. For output values see Wavefunction::dk().
  int dk() const;

  const Ref<LinearR12::CorrelationFactor>& corrfactor() const { return r12info()->corrfactor(); }
  LinearR12::ABSMethod abs_method() const { return r12info()->abs_method(); }
  const Ref<LinearR12Ansatz>& ansatz() const { return r12info()->ansatz(); }
  bool spin_polarized() const { return r12info()->refinfo()->spin_polarized(); }
  bool gbc() const { return r12info()->gbc(); }
  bool ebc() const { return r12info()->ebc(); }
  LinearR12::StandardApproximation stdapprox() const { return r12info()->stdapprox(); }
  bool omit_P() const { return r12info()->omit_P(); }

  const Ref<R12IntEvalInfo>& r12info() const;

  RefSCDimension dim_oo_s() const;
  RefSCDimension dim_oo_t() const;
  /// Dimension for active-occ/active-occ pairs of spin case S
  RefSCDimension dim_oo(SpinCase2 S) const;
  /// Dimension for active-vir/active-vir pairs of spin case S
  RefSCDimension dim_vv(SpinCase2 S) const;
  /// Dimension for any/any pairs of spin case S
  RefSCDimension dim_aa(SpinCase2 S) const;
  /// Dimension for geminal functions of spin case S = # of correlation factors x dim_xy
  RefSCDimension dim_f12(SpinCase2 S) const;
  /// Dimension of orbital product space used to generate geminal functions
  RefSCDimension dim_xy(SpinCase2 S) const;
  /// Returns the number of unique spin cases
  int nspincases1() const { return ::sc::nspincases1(spin_polarized()); }
  /// Returns the number of unique combinations of 2 spin cases
  int nspincases2() const { return ::sc::nspincases2(spin_polarized()); }

  /// This function causes the intermediate matrices to be computed.
  virtual void compute();

  /// Returns amplitudes of pair correlation functions
  Ref<F12Amplitudes> amps();

  /// Returns S block of intermediate V
  const RefSCMatrix& V(SpinCase2 S);
  /// Returns S block of intermediate X
  RefSymmSCMatrix X(SpinCase2 S);
  /// Returns S block of intermediate B
  RefSymmSCMatrix B(SpinCase2 S);
  /// Returns S block of the difference between intermediate B of approximations B and A'
  RefSymmSCMatrix BB(SpinCase2 S);
  /// Returns S block of intermediate A
  const RefSCMatrix& A(SpinCase2 S);
  /// Returns S block of intermediate T2
  const RefSCMatrix& T2(SpinCase2 S);
  /// Returns S block of intermediate F12
  const RefSCMatrix& F12(SpinCase2 S);

  /** Compute V = 1/2 g_{pq}^{\alpha\beta} R_{\alpha\beta}^{xy} */
  RefSCMatrix V(SpinCase2 spincase2,
                const Ref<OrbitalSpace>& p,
                const Ref<OrbitalSpace>& q);
  /// Compute P = RgR
  RefSymmSCMatrix P(SpinCase2 S);

  /// Returns the singles MP2 energy
  double emp2_singles();
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
  const Ref<OrbitalSpace>& occ_act(SpinCase1 S) const;
  /// Returns the occ space for spin case S
  const Ref<OrbitalSpace>& occ(SpinCase1 S) const;
  /// Returns the act vir space for spin case S
  const Ref<OrbitalSpace>& vir_act(SpinCase1 S) const;
  /// Returns the vir space for spin case S
  const Ref<OrbitalSpace>& vir(SpinCase1 S) const;
  /// Returns the geminal-generating orbital space for spin case S
  const Ref<OrbitalSpace>& xspace(SpinCase1 S) const;

  /// Form <P|h+J|x> space
  const Ref<OrbitalSpace>& hj_x_P(SpinCase1 S);
  /// Form <A|h+J|x> space
  const Ref<OrbitalSpace>& hj_x_A(SpinCase1 S);
  /// Form <p|h+J|x> space
  const Ref<OrbitalSpace>& hj_x_p(SpinCase1 S);
  /// Form <m|h+J|x> space
  const Ref<OrbitalSpace>& hj_x_m(SpinCase1 S);
  /// Form <a|h+J|x> space
  const Ref<OrbitalSpace>& hj_x_a(SpinCase1 S);
  /// Form <P|h+J|i> space
  const Ref<OrbitalSpace>& hj_i_P(SpinCase1 S);
  /// Form <A|h+J|i> space
  const Ref<OrbitalSpace>& hj_i_A(SpinCase1 S);
  /// Form <p|h+J|i> space
  const Ref<OrbitalSpace>& hj_i_p(SpinCase1 S);
  /// Form <m|h+J|i> space
  const Ref<OrbitalSpace>& hj_i_m(SpinCase1 S);
  /// Form <a|h+J|i> space
  const Ref<OrbitalSpace>& hj_i_a(SpinCase1 S);
  /// Form <P|h+J|p> space
  const Ref<OrbitalSpace>& hj_p_P(SpinCase1 S);
  /// Form <A|h+J|p> space
  const Ref<OrbitalSpace>& hj_p_A(SpinCase1 S);
  /// Form <p|h+J|p> space
  const Ref<OrbitalSpace>& hj_p_p(SpinCase1 S);
  /// Form <m|h+J|p> space
  const Ref<OrbitalSpace>& hj_p_m(SpinCase1 S);
  /// Form <a|h+J|p> space
  const Ref<OrbitalSpace>& hj_p_a(SpinCase1 S);
  /// Form <P|K|x> space
  const Ref<OrbitalSpace>& K_x_P(SpinCase1 S);
  /// Form <A|K|x> space
  const Ref<OrbitalSpace>& K_x_A(SpinCase1 S);
  /// Form <p|K|x> space
  const Ref<OrbitalSpace>& K_x_p(SpinCase1 S);
  /// Form <i|K|x> space
  const Ref<OrbitalSpace>& K_x_m(SpinCase1 S);
  /// Form <a|K|x> space
  const Ref<OrbitalSpace>& K_x_a(SpinCase1 S);
  /// Form <P|K|i> space
  const Ref<OrbitalSpace>& K_i_P(SpinCase1 S);
  /// Form <A|K|i> space
  const Ref<OrbitalSpace>& K_i_A(SpinCase1 S);
  /// Form <p|K|i> space
  const Ref<OrbitalSpace>& K_i_p(SpinCase1 S);
  /// Form <m|K|i> space
  const Ref<OrbitalSpace>& K_i_m(SpinCase1 S);
  /// Form <a|K|i> space
  const Ref<OrbitalSpace>& K_i_a(SpinCase1 S);
  /// Form <a|K|m> space
  const Ref<OrbitalSpace>& K_m_a(SpinCase1 S);
  /// Form <a|K|a> space
  const Ref<OrbitalSpace>& K_a_a(SpinCase1 S);
  /// Form <p|K|a> space
  const Ref<OrbitalSpace>& K_a_p(SpinCase1 S);
  /// Form <P|K|a> space
  const Ref<OrbitalSpace>& K_a_P(SpinCase1 S);
  /// Form <P|K|p> space
  const Ref<OrbitalSpace>& K_p_P(SpinCase1 S);
  /// Form <A|K|p> space
  const Ref<OrbitalSpace>& K_p_A(SpinCase1 S);
  /// Form <p|K|p> space
  const Ref<OrbitalSpace>& K_p_p(SpinCase1 S);
  /// Form <m|K|p> space
  const Ref<OrbitalSpace>& K_p_m(SpinCase1 S);
  /// Form <a|K|p> space
  const Ref<OrbitalSpace>& K_p_a(SpinCase1 S);
  /// Form <P|K|A> space
  const Ref<OrbitalSpace>& K_A_P(SpinCase1 S);
  /// Form <P|K|P> space
  const Ref<OrbitalSpace>& K_P_P(SpinCase1 S);
  /// Form <P|F|x> space
  const Ref<OrbitalSpace>& F_x_P(SpinCase1 S);
  /// Form <A|F|x> space
  const Ref<OrbitalSpace>& F_x_A(SpinCase1 S);
  /// Form <p|F|x> space
  const Ref<OrbitalSpace>& F_x_p(SpinCase1 S);
  /// Form <m|F|x> space
  const Ref<OrbitalSpace>& F_x_m(SpinCase1 S);
  /// Form <a|F|x> space
  const Ref<OrbitalSpace>& F_x_a(SpinCase1 S);
  /// Form <P|F|i> space
  const Ref<OrbitalSpace>& F_i_P(SpinCase1 S);
  /// Form <A|F|i> space
  const Ref<OrbitalSpace>& F_i_A(SpinCase1 S);
  /// Form <p|F|i> space
  const Ref<OrbitalSpace>& F_i_p(SpinCase1 S);
  /// Form <m|F|i> space
  const Ref<OrbitalSpace>& F_i_m(SpinCase1 S);
  /// Form <a|F|i> space
  const Ref<OrbitalSpace>& F_i_a(SpinCase1 S);
  /// Form <m|F|m> space
  const Ref<OrbitalSpace>& F_m_m(SpinCase1 S);
  /// Form <a|F|m> space
  const Ref<OrbitalSpace>& F_m_a(SpinCase1 S);
  /// Form <P|F|m> space
  const Ref<OrbitalSpace>& F_m_P(SpinCase1 S);
  /// Form <A|F|m> space
  const Ref<OrbitalSpace>& F_m_A(SpinCase1 S);
  /// Form <a|F|a> space
  const Ref<OrbitalSpace>& F_a_a(SpinCase1 S);
  /// Form <A|F|a> space
  const Ref<OrbitalSpace>& F_a_A(SpinCase1 S);
  /// Form <P|F|p> space
  const Ref<OrbitalSpace>& F_p_P(SpinCase1 S);
  /// Form <A|F|p> space
  const Ref<OrbitalSpace>& F_p_A(SpinCase1 S);
  /// Form <p|F|p> space
  const Ref<OrbitalSpace>& F_p_p(SpinCase1 S);
  /// Form <m|F|p> space
  const Ref<OrbitalSpace>& F_p_m(SpinCase1 S);
  /// Form <a|F|p> space
  const Ref<OrbitalSpace>& F_p_a(SpinCase1 S);
  /// Form <P|F|P> space
  const Ref<OrbitalSpace>& F_P_P(SpinCase1 S);

  /** Returns an already created transform.
      If the transform is not found then throw TransformNotFound */
  Ref<TwoBodyMOIntsTransform> get_tform_(const std::string&) const;
  /** Map transform T to label */
  void add_tform(const std::string& label,
                 const Ref<TwoBodyMOIntsTransform>& T);
  /// Generates canonical id for transform. no correlation function included
  std::string transform_label(const Ref<OrbitalSpace>& space1,
                              const Ref<OrbitalSpace>& space2,
                              const Ref<OrbitalSpace>& space3,
                              const Ref<OrbitalSpace>& space4,
                              const std::string& operator_label = std::string()) const;
  /// Generates canonical id for transform. f12 is the index of the correlation function
  std::string transform_label(const Ref<OrbitalSpace>& space1,
                              const Ref<OrbitalSpace>& space2,
                              const Ref<OrbitalSpace>& space3,
                              const Ref<OrbitalSpace>& space4,
                              unsigned int f12,
                              const std::string& operator_label = std::string()) const;
  /// version of transform_label() applicable when left and right correlation factors differ
  std::string transform_label(const Ref<OrbitalSpace>& space1,
                              const Ref<OrbitalSpace>& space2,
                              const Ref<OrbitalSpace>& space3,
                              const Ref<OrbitalSpace>& space4,
                              unsigned int f12_left,
                              unsigned int f12_right,
                              const std::string& operator_label = std::string()) const;

  /** Compute T2 amplitude in basis <space1, space3 | space2, space4>.
      AlphaBeta amplitudes are computed.
      If tform is not given (it should be!), this function will construct a generic
      transform. */
  void compute_T2_(RefSCMatrix& T2,
                   const Ref<OrbitalSpace>& space1,
                   const Ref<OrbitalSpace>& space2,
                   const Ref<OrbitalSpace>& space3,
                   const Ref<OrbitalSpace>& space4,
                   bool antisymmetrize,
                   const std::string& tform_key);
  /** Compute F12 integrals in basis <space1, space3 | f12 | space2, space4>.
      Bra (rows) are blocked by correlation function index.
      AlphaBeta amplitudes are computed.
      If tform is not given (it should be!), this function will construct a generic
      transform. */
  void compute_F12_(RefSCMatrix& F12,
                   const Ref<OrbitalSpace>& space1,
                   const Ref<OrbitalSpace>& space2,
                   const Ref<OrbitalSpace>& space3,
                   const Ref<OrbitalSpace>& space4,
                   bool antisymmetrize,
                   const std::vector<std::string>& transform_keys);

  /** Compute the Fock matrix between bra_ and ket_ spaces of spin S.
      scale_J and scale_K are used to scale Coulomb
      and exchange contributions, T12IntEval::occ() is used for the occupied spaces.
      */
  RefSCMatrix fock(const Ref<OrbitalSpace>& bra_space,
                   const Ref<OrbitalSpace>& ket_space, SpinCase1 S = Alpha,
                   double scale_J = 1.0, double scale_K = 1.0, double scale_H = 1.0);

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


