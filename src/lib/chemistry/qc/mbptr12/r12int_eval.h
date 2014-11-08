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

#ifndef _chemistry_qc_mbptr12_r12inteval_h
#define _chemistry_qc_mbptr12_r12inteval_h

#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/r12technology.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/wfn/spin.h>
#include <chemistry/qc/mbptr12/fixedcoefficient.h>
#include <chemistry/qc/mbptr12/twobodytensorinfo.h>

namespace sc {

  class TwoBodyMOIntsTransform;
  class DistArray4;
  class R12Amplitudes;
  class R12EnergyIntermediates;

  /** R12IntEval is the top-level class which computes intermediates occuring in R12 theories.
      This class is used by all Wavefunction classes that implement R12 methods.
  */

class R12IntEval : virtual public SavableState {
  private:
  // change to false to use the old fock builder
  static const bool USE_FOCKBUILD = true;


  bool evaluated_;

  // Calculation information (number of basis functions, R12 approximation, etc.)
  Ref<R12WavefunctionWorld> r12world_;

  RefSCMatrix V_[NSpinCases2];
  RefSCMatrix X_[NSpinCases2];
  // Note that intermediate B is symmetric but is stored as a full matrix to simplify the code
  // that computes asymmetric form of B
  RefSCMatrix B_[NSpinCases2];
  RefSCMatrix BB_[NSpinCases2];  // The difference between B intermediate of approximation B and A'
  RefSCMatrix A_[NSpinCases2];
  /// array of Ref<CuspConsistentGeminalCoefficient> for the generation of cusp consistent geminal coefficients
  Ref<CuspConsistentGeminalCoefficient> cuspconsistentgeminalcoefficient_[NSpinCases2];

  RefSCVector emp2pair_[NSpinCases2];
  RefSCDimension dim_oo_[NSpinCases2];
  RefSCDimension dim_vv_[NSpinCases2];
  /// a is any index in a given basis
  RefSCDimension dim_aa_[NSpinCases2];
  /// space of all geminal functions
  RefSCDimension dim_f12_[NSpinCases2];
  /// space of orbital products multiplied by f(r12) to produce the above
  RefSCDimension dim_GG_[NSpinCases2];
  /// space of orbital products from which geminal substitutions are allowed
  RefSCDimension dim_gg_[NSpinCases2];

  Ref<R12Amplitudes> Amps_;  // First-order amplitudes of various contributions to the pair functions
  double emp2_obs_singles_;
  double emp2_cabs_singles_;
  RefSCMatrix T1_cabs_[NSpinCases1];
  int debug_;

  mutable RefSymmSCMatrix ordm_[NSpinCases1];  //!< spin-orbital 1-RDM in MO basis

  /// "Spin-adapt" MO space id and name
  void spinadapt_mospace_labels(SpinCase1 spin, std::string& id, std::string& name) const;

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
  Ref<OrbitalSpace> hj_m_m_[NSpinCases1];
  Ref<OrbitalSpace> hj_m_p_[NSpinCases1];
  Ref<OrbitalSpace> hj_a_A_[NSpinCases1];
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
  Ref<OrbitalSpace> F_m_p_[NSpinCases1];
  Ref<OrbitalSpace> F_m_P_[NSpinCases1];
  Ref<OrbitalSpace> F_gg_P_[NSpinCases1];
  Ref<OrbitalSpace> F_m_A_[NSpinCases1];
  Ref<OrbitalSpace> F_i_A_[NSpinCases1];
  Ref<OrbitalSpace> F_i_m_[NSpinCases1];
  Ref<OrbitalSpace> F_i_a_[NSpinCases1];
  Ref<OrbitalSpace> F_i_p_[NSpinCases1];
  Ref<OrbitalSpace> F_i_P_[NSpinCases1];
  Ref<OrbitalSpace> F_a_a_[NSpinCases1];
  Ref<OrbitalSpace> F_a_A_[NSpinCases1];
  Ref<OrbitalSpace> h_P_P_[NSpinCases1];
  Ref<OrbitalSpace> J_i_p_[NSpinCases1];
  Ref<OrbitalSpace> J_i_P_[NSpinCases1];
  Ref<OrbitalSpace> J_P_P_[NSpinCases1];
  Ref<OrbitalSpace> F_A_A_[NSpinCases1];
  Ref<OrbitalSpace> F_p_P_[NSpinCases1];
  Ref<OrbitalSpace> gamma_p_p_[NSpinCases1];
  Ref<OrbitalSpace> gamma_m_m_[NSpinCases1];
  Ref<OrbitalSpace> gammaFgamma_p_p_[NSpinCases1];
  Ref<OrbitalSpace> Fgamma_p_P_[NSpinCases1];

  /// compute CABS space for spin s canonicalized by diagonalization of the Fock matrix with adjustable component weights
  Ref<OrbitalSpace> cabs_space_fockcanonical(SpinCase1 s,
                                             double scale_H,
                                             double scale_J,
                                             double scale_K);

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

#if 0
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
#endif

  // Computes the skalar Pauli-Hamiltonian (T + V + mass_velocity + Darwin),
  // with the mass-velocity term evaluated in the momentum basis.
  RefSymmSCMatrix pauli(const Ref<GaussianBasisSet> &bas,
                        const Ref<GaussianBasisSet> &pbas = 0,
                        const bool momentum=false);
  RefSymmSCMatrix pauli_realspace(const Ref<GaussianBasisSet> &bas);
  RefSymmSCMatrix pauli_momentumspace(const Ref<GaussianBasisSet> &bas,
                        const Ref<GaussianBasisSet> &p_bas);

  /// Compute the coulomb matrix between 2 spaces. Coulomb contribution computed w.r.t. density of spincase spin only
  RefSCMatrix coulomb_(const SpinCase1 &spin, const Ref<OrbitalSpace>& bra_space,
                       const Ref<OrbitalSpace>& ket_space);
  RefSCMatrix coulomb_(const Ref<OrbitalSpace>& occ_space, const Ref<OrbitalSpace>& bra_space,
                       const Ref<OrbitalSpace>& ket_space);
  /// Compute the exchange matrix between 2 spaces
  RefSCMatrix exchange_(const SpinCase1 &spin, const Ref<OrbitalSpace>& bra_space,
                        const Ref<OrbitalSpace>& ket_space);
  RefSCMatrix exchange_(const Ref<OrbitalSpace>& occ_space, const Ref<OrbitalSpace>& bra_space,
                        const Ref<OrbitalSpace>& ket_space);

  /// Checkpoint the top-level molecular energy
  void checkpoint_() const;

  /// New version which uses optimized contract functions
  void contrib_to_VXB_a_();
  void contrib_to_VXB_c_ansatz1_();
  void contrib_to_VX_GenRefansatz2_();
  void contrib_to_VX_GenRefansatz2_spinfree_();
  /// New version which uses tensor contract functions
  void contrib_to_VXB_a_vbsneqobs_();
  /// extra single-commutator contributions to B from relativistic terms
  void contrib_to_B_DKH_a_();

  /// Compute MP2 pair energies of spin case S using < space1 space3|| space2 space4> integrals
  void compute_mp2_pair_energies_(RefSCVector& emp2pair,
                                  SpinCase2 S,
                                  const Ref<OrbitalSpace>& space1,
                                  const Ref<OrbitalSpace>& space2,
                                  const Ref<OrbitalSpace>& space3,
                                  const Ref<OrbitalSpace>& space4,
                                  const std::string& tform_key);

  void compute_A(); //< computes A ("coupling") matrix for all spincases

  /** Compute A intermediate using "direct" formula in basis <space1, space3 | f12 | space2, space4>.
      Bra (rows) are blocked by correlation function index.
      AlphaBeta amplitudes are computed.
      If tform is not given (it should be!), this function will construct a generic
      transform.  */
  void compute_A_direct_(RefSCMatrix& A,
                         const Ref<OrbitalSpace>& space1,
                         const Ref<OrbitalSpace>& space2,
                         const Ref<OrbitalSpace>& space3,
                         const Ref<OrbitalSpace>& space4,
                         const Ref<OrbitalSpace>& rispace2,
                         const Ref<OrbitalSpace>& rispace4,
                         bool antisymmetrize);

  // make compute_tbint_tensor_ public for now
  public:
  /** compute_tbint_tensor computes a 2-body tensor T using integrals of type tbint_type.
      Computed tensor T is added to its previous contents.
      Class DataProcess defines a static function 'double I2T()' which processes the integrals.
      Set CorrFactorInBra to true if bra of target tensor depends on correlation function index.

      Of course, this ugliness should become constructor/member function of ManyBodyOperator

      template parameters:
      DataProcess                       -- classes in namespace ManyBodyTensors that describe
                                           the operator whose matrix
                                           elements are needed, e.g. Apply_H0minusE0<sign>
      CorrFactorInBra, CorrFactorInKet  -- 'true' if there is a correlation factor.

      function parameters:
      tbint_type                        -- type of tensor to be computed.
      space1, space2, space3, space4    -- index spaces of the tensor in chemist's (Mulliken) notation:
                                           ( space1 space2 | space3 space4 ).
      antisymmetrize                    -- if true, the computed tensor is antisymmetrized.
      tforms                            -- TwoBodyMOIntsTransform vector (dimension: number
                                           of correlation factors) that describes the transformation
                                           of the  tensor from the AO to the MO space.
      tbintdescrs                       -- integral descriptor of the tensor to be computed.

   */
  template <typename DataProcess, bool CorrFactorInBra, bool CorrFactorInKet>
  DEPRECATED void compute_tbint_tensor(RefSCMatrix& T,
                              TwoBodyOper::type tbint_type,
                              const Ref<OrbitalSpace>& space1,
                              const Ref<OrbitalSpace>& space2,
                              const Ref<OrbitalSpace>& space3,
                              const Ref<OrbitalSpace>& space4,
                              bool antisymmetrize,
                              const std::vector<std::string>& tform_keys);

  private:
  /**
     contract_tbint_tensor computes a 2-body tensor T as a sum over mn : <ij|Tbra|mn> * <kl|Tket|mn>^t
     bra and ket integrals come from tforms_bra and tforms_ket.
     Computed tensor T is added to its previous contents.
     Class DataProcess_XXX defines a static function 'double I2T()' which processes the integrals.
     The I2T() functions are implemented as members of classes in the namespace ManyBodyTensors.
     Set CorrFactorInBra to true if bra of target tensor depends on correlation function index.

     Semantics: 1) ranks of internal spaces must match, although the spaces don't have to be the same;
     2) if antisymmetrize is true, external (bra and ket) particle spaces
     should be strongly (identity) or weakly(rank) equivalent. If internal spaces do not match
     but antisymmetrization is requested, symmetrization w.r.t. particle swap will be performed,
      then antisymmetrization.

     Of course, this ugliness should become function/operator on 2 ManyBodyOperators

     template parameters:
     DataProcessBra, DataProcessKet,
     DataProcessBraKet                 -- classes in namespace ManyBodyTensors that describe the operator whose matrix
                                          elements are needed, e.g. Apply_H0minusE0<sign>
     CorrFactorInBra, CorrFactorInKet,
     CorrFactorInInt                   -- 'true' if there is a correlation factor.

     function parameters:
     tbint_type_bra and tbint_type_ket -- type of the first and second tensor.
     space1_bra and space2_bra         -- space of the first and second bra index of the first tensor - external indices.
     space1_intb and space2_intb       -- space of the first and second ket index of the first tensor - internal indices.
     space1_ket and space2_ket         -- space of the first and second ket index of the second tensor - external indices.
     space1_intk and space2_intk       -- space of the first and second bra index of the second tensor - internal indices.
     tpcontract                        -- provides information on how the two tensors are to be contracted - see
                                          document of class R12Technology::TwoParticleContraction and the classes derived
                                          from it.
     antisymmetrize                    -- antisymmetrize the final contracted product.
     tforms_bra and tforms_ket         -- TwoBodyMOIntsTransform vectors (dimension: number of correlation factors) that
                                          describe the transformation of the bra and ket tensor respectively from the AO
                                          MO space.
     intdescrs_bra and intdescrs_ket   -- integral descriptors belonging to the bra and ket tensor respectively.
     <space1_bra space2_bra|T1|space1_intb space2_intb> * <space1_intk space2_intk|T2|space1_ket space2_ket>
  */
  template <typename DataProcessBra,
            typename DataProcessKet,
            typename DataProcessBraKet,
            bool CorrFactorInBra,
            bool CorrFactorInKet,
            bool CorrFactorInInt>
  DEPRECATED void contract_tbint_tensor(
           RefSCMatrix& T,
           TwoBodyOper::type tbint_type_bra,
           TwoBodyOper::type tbint_type_ket,
           const Ref<OrbitalSpace>& space1_bra,
           const Ref<OrbitalSpace>& space2_bra,
           const Ref<OrbitalSpace>& space1_intb,
           const Ref<OrbitalSpace>& space2_intb,
           const Ref<OrbitalSpace>& space1_ket,
           const Ref<OrbitalSpace>& space2_ket,
           const Ref<OrbitalSpace>& space1_intk,
           const Ref<OrbitalSpace>& space2_intk,
           const Ref<mbptr12::TwoParticleContraction>& tpcontract,
           bool antisymmetrize,
           const std::vector<std::string>& tformkeys_bra,
           const std::vector<std::string>& tformkeys_ket);

  /** overload of the above when no pre- and post-processing is needed.
      this version is also much more efficient since it does contraction
      as a DGEMM (hence it tiles loads).

      \param antisymmetrize indicates whether the target tensor is antisymmetric w.r.t permutation of
      particles or now.
   */
  template <bool CorrFactorInBra,
            bool CorrFactorInKet>
    void DEPRECATED contract_tbint_tensor(
           RefSCMatrix& T,
           TwoBodyOper::type tbint_type_bra,
           TwoBodyOper::type tbint_type_ket,
           double scale,
           const Ref<OrbitalSpace>& space1_bra,
           const Ref<OrbitalSpace>& space2_bra,
           const Ref<OrbitalSpace>& space1_intb,
           const Ref<OrbitalSpace>& space2_intb,
           const Ref<OrbitalSpace>& space1_ket,
           const Ref<OrbitalSpace>& space2_ket,
           const Ref<OrbitalSpace>& space1_intk,
           const Ref<OrbitalSpace>& space2_intk,
           bool antisymmetrize,
           const std::vector<std::string>& tformkeys_bra,
           const std::vector<std::string>& tformkeys_ket);

  /** overload of the above when the result should be stored as a DistArray4

      \param antisymmetrize indicates whether the target tensor is antisymmetric w.r.t permutation of
      particles or now.
   */
  template <bool CorrFactorInBra,
            bool CorrFactorInKet>
  void
  contract_tbint_tensor(
           std::vector< Ref<DistArray4> >& results,
           TwoBodyOper::type tbint_type_bra,
           TwoBodyOper::type tbint_type_ket,
           double scale,
           const Ref<OrbitalSpace>& space1_bra,
           const Ref<OrbitalSpace>& space2_bra,
           const Ref<OrbitalSpace>& space1_intb,
           const Ref<OrbitalSpace>& space2_intb,
           const Ref<OrbitalSpace>& space1_ket,
           const Ref<OrbitalSpace>& space2_ket,
           const Ref<OrbitalSpace>& space1_intk,
           const Ref<OrbitalSpace>& space2_intk,
           bool antisymmetrize,
           const std::vector<std::string>& tformkeys_bra,
           const std::vector<std::string>& tformkeys_ket);

  // for now make it public
  public:
  /** <space1bra space1_intb |Tbra| space2_intb space3_intb> * <space2_intk space3_intk |Tket| space1ket space1_intk>
   *  contract_tbint_tensors_to_obtensor computes a 1-body tensor T as a sum over kmn : <ik|Tbra|mn><jk|Tket|mn>^t.
   *  The notation used here is analoguous to that of the routine contract_tbint_tensor.
   *  bra and ket integrals come from tforms_bra and tforms_ket.
   *  Computed tensor T is added to its previous contents.
   *  Class DataProcess_XXX defines a static function 'double I2T()' which processes the integrals.
   *  The I2T() functions are implemented as members of classes in the namespace ManyBodyTensors.
   *  Set CorrFactorInBra to true if bra of target tensor depends on correlation function index.
   *
   *  Semantics: ranks of internal spaces must match, although the spaces don't have to be the same.
   *
   *  Of course, this ugliness should become function/operator on 2 ManyBodyOperators
   *
   *  template parameters:
   *  DataProcessBra, DataProcessKet,   -- classes in namespace ManyBodyTensors that describe the operator whose matrix
   *                                       elements are needed, e.g. Apply_H0minusE0<sign>
   *  CorrFactorInBra, CorrFactorInKet,
   *  CorrFactorInInt                   -- 'true' if there is a correlation factor.
   *
   *  function parameters:
   *  tbtensor_type_bra,
   *  tbtensor_type_ket                 -- type of the first and second tensor.
   *  space1_bra                        -- space of the first bra index of the first tensor - external index.
   *  space1_intb, space2_intb,
   *  space3_intb                       -- space of the first bra and first and second ket index of the first tensor - internal indices.
   *  space1_ket                        -- space of the first ket index of the second tensor - external index.
   *  space1_intk, space2_intk,
   *  space3_intk                       -- space of the first bra and first and second ket index of the second tensor - internal indices.
   *  tpcontract                        -- provides information on how the two tensors are to be contracted - see
   *                                       document of class R12Technology::TwoParticleContraction and the classes derived
   *                                       from it.
   *  tforms_bra and tforms_ket         -- TwoBodyMOIntsTransform vectors (dimension: number of correlation factors) that
   *                                       describe the transformation of the bra and ket tensor respectively from the AO
   *                                       MO space.
   *  intdescrs_bra and intdescrs_ket   -- integral descriptors belonging to the bra and ket tensor respectively. */
  template <typename DataProcessBra,
            typename DataProcessKet,
            bool CorrFactorInBra,
            bool CorrFactorInKet,
            bool CorrFactorInInt>
    void contract_tbint_tensors_to_obtensor(RefSCMatrix& T,
                                            SpinCase2 pairspin,
                                            TwoBodyTensorInfo tbtensor_type_bra,
                                            TwoBodyTensorInfo tbtensor_type_ket,
                                            const Ref<OrbitalSpace>& space1_bra,
                                            const Ref<OrbitalSpace>& space1_intb,
                                            const Ref<OrbitalSpace>& space2_intb,
                                            const Ref<OrbitalSpace>& space3_intb,
                                            const Ref<OrbitalSpace>& space1_ket,
                                            const Ref<OrbitalSpace>& space1_intk,
                                            const Ref<OrbitalSpace>& space2_intk,
                                            const Ref<OrbitalSpace>& space3_intk,
                                            const Ref<mbptr12::TwoParticleContraction>& tpcontract,
                                            const std::vector<std::string>& tformkeys_bra,
                                            const std::vector<std::string>& tformkeys_ket);

  private:
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

      Example: term \f$\overline{r}^{p' q'}_{\bf p q} k^{s'}_{p'} \overline{r}^{\bf r s}_{s' q'} \f$.
      All tensors are generated from the integrals over AO's by means of the MO coefficients. In this
      case, \f$ \overline{r}^{p' q'}_{\bf p q} k^{s'}_{p'} = \overline{r}^{s'_k q'}_{\bf p q}\f$ as an
      ordinary tensor matrix element of \f$\overline{r}\f$, but with the new MO coefficient
      \f$C^{s'_k}_{\mu'} = C^{p'}_{\mu'} k^{s'}_{p'} \f$. compute_FxF_ then computes the product
      \f$\overline{r}^{s'_k q'}_{\bf p q} \overline{r}^{\bf r s}_{s' q'}\f$. The OrbitalSpace of
      \f$C^{s'_k}_{\mu'}\f$ is generated by the function K_p_p. For each internal index, the spin1
      and spin2 version is needed. Because of this there as six internal indices.
      <bra1 bra2|F12|intkx1 int2> * <intk1 int2|F12|ket1 ket2> + <bra1 bra2|F12|int1 intkx2> * <int1 intk2|F12|ket1 ket2>.
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

  /** Compute the fX contribution to B (only appears in stdapprox = A' or B).
      This includes contributions that depend on Brillouin condition. */
  void compute_B_fX_();

  /** Compute contributions to B which appear in standard approximation B and not in A' */
  void compute_BB_();

  /** Compute B using standard approximation C */
  void compute_BC_();
  void compute_BC_ansatz1_();
  void compute_BC_GenRefansatz2_();
  void compute_BC_GenRefansatz2_spinfree();


  /** Compute B using standard approximation C' */
  void compute_BCp_();

  /** Compute B using standard approximation A'' -- exchange is dropped completely! */
  void compute_BApp_();

  /** Compute the mass-velocity contributions to B that appear when DKH-based R12 calculations are requested */
  void compute_B_DKH_();

  /// computes MP2 energy for single-determinant reference
  const RefSCVector& compute_emp2(SpinCase2 s);
  /** Compute OBS singles contribution to the MP2 energy.
      If obs_singles is set to true, use OBS virtuals, else use correlating virtuals (these
      virtuals differ if VBS != OBS).
    */
  double compute_emp2_obs_singles(bool obs_singles);
#if 0
  /** Compute CABS singles contribution to the MP2 energy. This assumes occupied and virtual (OBS) orbitals to be
   *  canonical.
   *  */
  double compute_emp2_cabs_singles();
#endif
  /** Compute CABS singles contribution to the MP2 energy without
   * assuming canonical occupied or virtual orbitals.
   * \sa R12IntEval::compute_emp2_cabs_singles
   *
   * @param vir_cabs_coupling if true, will couple conventional and CABS T1
   */
  double compute_emp2_cabs_singles_noncanonical(bool vir_cabs_coupling);
  /** Compute CABS singles contribution to the CCSD energy without
   * assuming canonical occupied orbitals.
   * \sa R12IntEval::compute_emp2_cabs_singles_noncanonical */
  double compute_emp2_cabs_singles_noncanonical_ccsd(const RefSCMatrix& T1_ia_alpha,
                                                     const RefSCMatrix& T1_ia_beta);

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
  R12IntEval(const Ref<R12WavefunctionWorld>& r12w);
  ~R12IntEval();

  void save_data_state(StateOut&);
  virtual void obsolete();

  void debug(int d) { debug_ = d; }

  int dk() const { return r12world()->refwfn()->dk(); }

#if 1
  const Ref<R12Technology::CorrelationFactor>& corrfactor() const { return r12world()->r12tech()->corrfactor(); }
  R12Technology::ABSMethod abs_method() const { return r12world()->r12tech()->abs_method(); }
  const Ref<R12Technology::R12Ansatz>& ansatz() const { return r12world()->r12tech()->ansatz(); }
  bool spin_polarized() const { return r12world()->refwfn()->spin_polarized(); }
  bool gbc() const { return r12world()->r12tech()->gbc(); }
  bool ebc() const { return r12world()->r12tech()->ebc(); }
  bool coupling() const { return r12world()->r12tech()->coupling(); }
  bool compute_1rdm() const { return r12world()->r12tech()->compute_1rdm(); }
  R12Technology::StandardApproximation stdapprox() const { return r12world()->r12tech()->stdapprox(); }
  bool omit_P() const { return r12world()->r12tech()->omit_P(); }
  const Ref<MOIntsTransformFactory>& tfactory() const { return r12world()->world()->tfactory(); };
  const Ref<MOIntsRuntime>& moints_runtime() const { return r12world()->world()->moints_runtime(); };
  const Ref<TwoBodyFourCenterMOIntsRuntime>& moints_runtime4() const { return r12world()->world()->moints_runtime()->runtime_4c(); };
  const Ref<FockBuildRuntime>& fockbuild_runtime() const { return r12world()->world()->fockbuild_runtime(); };

#endif

  /// the R12World in which this object lives
  const Ref<R12WavefunctionWorld>& r12world() const { return r12world_; }

  /// Dimension for active-occ/active-occ pairs of spin case S
  RefSCDimension dim_oo(SpinCase2 S) const { return dim_oo_[S]; }
  /// Dimension for active-vir/active-vir pairs of spin case S
  RefSCDimension dim_vv(SpinCase2 S) const { return dim_vv_[S]; }
  /// Dimension for any/any pairs of spin case S
  RefSCDimension dim_aa(SpinCase2 S) const { return dim_aa_[S]; }
  /// Dimension for geminal functions of spin case S = # of correlation factors x dim_GG
  RefSCDimension dim_f12(SpinCase2 S) const { return dim_f12_[S]; }
  /// Dimension of orbital product space that multiply the correlation factor to produce geminal functions
  RefSCDimension dim_GG(SpinCase2 S) const { return dim_GG_[S]; }
  /// Dimension of orbital product space from which geminal substitutions are allowed
  RefSCDimension dim_gg(SpinCase2 S) const { return dim_gg_[S]; }

  /// Returns the number of unique spin cases
  int nspincases1() const {
    if (r12world()->spinadapted()) return 1;
    return ::sc::nspincases1(r12world()->refwfn()->spin_polarized());
  }
  /// Returns the number of unique combinations of 2 spin cases
  int nspincases2() const {
    if (r12world()->spinadapted()) return 1;
    return ::sc::nspincases2(r12world()->refwfn()->spin_polarized());
  }

  /// This function causes the intermediate matrices to be computed.
  virtual void compute();

  /// does Brillouin condition hold?
  bool bc() const;

  /// Returns amplitudes of pair correlation functions
  Ref<R12Amplitudes> amps();

  /// Returns S block of intermediate V
  const RefSCMatrix& V(SpinCase2 S);
  /// Returns spin-free intermediate V
  const RefSCMatrix& V();
  /// Returns S block of intermediate X
  RefSymmSCMatrix X(SpinCase2 S);
  /// Returns spin-free intermediate X
  RefSymmSCMatrix X();
  /// Returns S block of intermediate B
  RefSymmSCMatrix B(SpinCase2 S);
  /// Returns spin-free intermediate B
  RefSymmSCMatrix B();
  /// Returns S block of the difference between intermediate B of approximations B and A'
  RefSymmSCMatrix BB(SpinCase2 S);
  /// Returns S block of intermediate A
  const RefSCMatrix& A(SpinCase2 S);
  /// Returns S block of intermediate T2
  const RefSCMatrix& T2(SpinCase2 S);
  /// Returns S block of intermediate F12
  const RefSCMatrix& F12(SpinCase2 S);

  /** Compute \f$ V_{pq}^{xy} = \frac{1}{2} \bar{g}_{pq}^{\alpha\beta} \bar{R}_{\alpha\beta}^{xy}\f$ */
  RefSCMatrix V(SpinCase2 spincase2,
                const Ref<OrbitalSpace>& p,
                const Ref<OrbitalSpace>& q);
  /** Compute \f$ V_{pq}^{xy} = \frac{1}{2} \bar{g}_{pq}^{\alpha\beta} \bar{R}_{\alpha\beta}^{xy}\f$. \sa R12IntEval::V() */
  std::vector<Ref<DistArray4> > V_distarray4(SpinCase2 spincase2,
                               const Ref<OrbitalSpace>& p,
                               const Ref<OrbitalSpace>& q);
  /** Compute \f$ U_{pq}^{xy} = \frac{1}{2} \bar{R}_{p \alpha'}^{x l} \bar{g}_{q l}^{y \alpha'} \f$. */
  std::vector<Ref<DistArray4> > U_distarray4(SpinCase2 spincase2,
                                             const Ref<OrbitalSpace>& p,
                                             const Ref<OrbitalSpace>& q);
  /// Compute \f$ P_{uv}^{xy} = \frac{1}{4} \bar{R}^{\alpha\beta}_{uv} \bar{g}_{\alpha\beta}^{\gamma\delta} \bar{R}_{\gamma\delta}^{xy}\f$ P = RgR
  RefSymmSCMatrix P(SpinCase2 S);
  /** Returns the cusp consistent coefficient \f$C_{ij}^{kl}\f$. */
  double C_CuspConsistent(int i,int j,int k,int l,SpinCase2 pairspin);

  /// Returns the OBS singles correction to the MP2 energy
  double emp2_obs_singles();
  /**
   * Returns the CABS singles correction to the MP2 energy
   * @param vir_cabs_coupling if true, will couple conventional and CABS T1's,
   *        hence the OBS singles correction will not need to be computed separately
   */
  double emp2_cabs_singles(bool vir_cabs_coupling = true);
  /// Returns the CABS singles MP2 energy, with fixed conventional T1 amplitudes
  double emp2_cabs_singles(const RefSCMatrix& T1_ia_alpha,
                           const RefSCMatrix& T1_ia_beta);
  /// Returns alpha-alpha MP2 pair energies
  const RefSCVector& emp2(SpinCase2 S);

  /**
   * returns CABS singles amplitudes. This includes excitations into standard virtuals if
   * this is a correction to MP2, not CCSD. Must call emp2_cabs_singles() BEFORE calling this.
   * @param spin SpinCase
   * @return the amplitude matrix (MP2: occ by allvirt=vir+CABS; CCSD: occ by CABS)
   */
  const RefSCMatrix& T1_cabs(SpinCase1 spin) const;

  /// Returns the act occ space for spin case S
  const Ref<OrbitalSpace>& occ_act(SpinCase1 S) const;
  /// Returns the occ space for spin case S
  const Ref<OrbitalSpace>& occ(SpinCase1 S) const;
  /// Returns the act vir space for spin case S
  const Ref<OrbitalSpace>& vir_act(SpinCase1 S) const;
  /// Returns the vir space for spin case S
  const Ref<OrbitalSpace>& vir(SpinCase1 S) const;
  /// Returns the OBS space for spin case S
  const Ref<OrbitalSpace>& orbs(SpinCase1 S) const;
  /// Returns the geminal-generating orbital space for spin case S
  const Ref<OrbitalSpace>& GGspace(SpinCase1 S) const;
  /// Returns the space for spin case S from which geminal-generating substitutions are allowed
  const Ref<OrbitalSpace>& ggspace(SpinCase1 S) const;
  /// compute canonical CABS space for spin s
  const Ref<OrbitalSpace>& cabs_space_canonical(SpinCase1 s);
  /// compute CABS space for spin s canonicalized by diagonalization of core hamitonian
  const Ref<OrbitalSpace>& cabs_space_hcanonical(SpinCase1 s);
  /// Returns the 1-RDM for spin S in the ``MO'' basis (i.e. that provided by orbs(S) )
  RefSymmSCMatrix ordm(SpinCase1 S) const;
  /// Returns the total 1-RDM in the ``MO'' basis (i.e. that provided by orbs() )
  RefSymmSCMatrix ordm() const;
  /// Returns the average 1-RDM in the ``MO'' basis (i.e. that provided by orbs() ): (alpha-rdm1 + beta-rdm2)/2
  RefSymmSCMatrix ordm_av() const;
  RefSymmSCMatrix ordm_occ_av() const;
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
  /// Form <m|h+J|m> space
  const Ref<OrbitalSpace>& hj_m_m(SpinCase1 S);
  /// Form <p|h+J|m> space
  const Ref<OrbitalSpace>& hj_m_p(SpinCase1 S);
  /// Form <A|h+J|a> space
  const Ref<OrbitalSpace>& hj_a_A(SpinCase1 S);
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
  /// Form <P|h+J|P> space
  const Ref<OrbitalSpace>& hj_P_P(SpinCase1 S);
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
  /// Form <p|F|m> space
  const Ref<OrbitalSpace>& F_m_p(SpinCase1 S);
  /// Form <P|F|m> space
  const Ref<OrbitalSpace>& F_m_P(SpinCase1 S);
  /// Form <P|F|gg> space
  const Ref<OrbitalSpace>& F_gg_P(SpinCase1 S);
  /// Form <A|F|m> space
  const Ref<OrbitalSpace>& F_m_A(SpinCase1 S);
  /// Form <a|F|a> space
  const Ref<OrbitalSpace>& F_a_a(SpinCase1 S);
  /// Form <A|F|a> space
  const Ref<OrbitalSpace>& F_a_A(SpinCase1 S);
  /// Form <A|F|A> space
  const Ref<OrbitalSpace>& F_A_A(SpinCase1 S);
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
  /// Form <P|h|P> space
  const Ref<OrbitalSpace>& h_P_P(SpinCase1 S);
  /// Form <p|J|x> space
  const Ref<OrbitalSpace>& J_x_p(SpinCase1 S);
  /// Form <p|J|i> space
  const Ref<OrbitalSpace>& J_i_p(SpinCase1 S);
  /// Form <P|J|i> space
  const Ref<OrbitalSpace>& J_i_P(SpinCase1 S);
  /// Form <P|J|P> space
  const Ref<OrbitalSpace>& J_P_P(SpinCase1 S);

  /// Form <p|gamma|p> space
  const Ref<OrbitalSpace>& gamma_p_p(SpinCase1 S);
  /// Form spin-average <p|gamma|p> space using average (instead of total) rdm.
    const Ref<OrbitalSpace>& gamma_p_p_av();
    const Ref<OrbitalSpace>& gamma_m_m_av();
  /// Form <p|gammaFgamma|p> space
  const Ref<OrbitalSpace>& gammaFgamma_p_p(SpinCase1 S);
  const Ref<OrbitalSpace>& gammaFgamma_p_p();
  /// Form <P|Fgamma|p> space
  const Ref<OrbitalSpace>& Fgamma_p_P(SpinCase1 S);
  /// Form <P|Fgamma|p> space using spin-average rdm1
  const Ref<OrbitalSpace>& Fgamma_p_P();
//  const Ref<OrbitalSpace>& Fgamma_P_p(SpinCase1 S,const RefSymmSCMatrix &gamma,Ref<OrbitalSpace>& FGamma);
//  const Ref<OrbitalSpace>& gammaF_p_P(SpinCase1 S,const RefSymmSCMatrix &gamma,Ref<OrbitalSpace>& GammaF);
  /** Form <A|obtensor|p> space
   *  obtensor should have the dimension ncabs*norbs */
  Ref<OrbitalSpace> obtensor_p_A(const RefSCMatrix &obtensor,SpinCase1 S);

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
                   double scale_J = 1.0, double scale_K = 1.0, double scale_H = 1.0,
                   // leave at default to use R12Technology's pauli flag, else set to 0 or 1
                   int override_pauli = -1);


  /** Compute V intermediates using CABS/CABS+ approach */
  RefSCMatrix V_cabs(SpinCase2 spincase2,
                     const Ref<OrbitalSpace>& p,
                     const Ref<OrbitalSpace>& q);

  /** Compute V intermediates in SF-[2]R12 */
  RefSCMatrix V_genref_spinfree(const Ref<OrbitalSpace>& p,
                     const Ref<OrbitalSpace>& q);

#if defined(HAVE_MPQC3_RUNTIME)
  /**
   * TiledArray-based builder of closed-shell V intermediate
   */
  void V_diag_ta();

  void gf2_r12();

  // TiledArray-based MP2-F12 one-electron properties
  void compute_TA_mp2f12_1rdm();
#endif

  void compute_ccr12_1rdm(const RefSCMatrix& T1, const Ref<DistArray4> (&T2)[NSpinCases2]);

  /// returns the OrbitalSpaceRegistry object
  const Ref<OrbitalSpaceRegistry>& orbital_registry() const {
    return this->r12world()->world()->tfactory()->orbital_registry();
  }
  /// returns the AOSpaceRegistry object
  const Ref<AOSpaceRegistry>& ao_registry() const {
    return this->r12world()->world()->tfactory()->ao_registry();
  }

};

std::vector< Ref<DistArray4> >
A_distarray4(SpinCase2 spincase2, const Ref<R12IntEval>& r12eval);

// compute orbital Z-vector from F12 contribution
RefSCMatrix Onerdm_X_F12(SpinCase1 spin, const Ref<R12IntEval>& r12eval, int debug);

// compute orbital relaxation contributions from CABS Singles
RefSCMatrix Onerdm_X_CABS_Singles(SpinCase1 spin,
                                  const Ref<R12IntEval>& r12eval,
                                  const Ref<R12EnergyIntermediates>& r12intermediates,
                                  int debug);

} // end of namespace sc

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


