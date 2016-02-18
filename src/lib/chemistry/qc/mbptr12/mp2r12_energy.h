//
// mp2r12_energy.h
//
// Copyright (C) 2003 Edward Valeev
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

#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/r12technology.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/wfn/spin.h>
#include <chemistry/qc/mbptr12/twobodygrid.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/mbptr12/mp2r12_energy_util.h>

#ifndef _chemistry_qc_mbptr12_mp2r12energy_h
#define _chemistry_qc_mbptr12_mp2r12energy_h

#define MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION 1

namespace sc {
  /**
   * The class R12EnergyIntermediates is the front-end to
   * R12 intermediates.
   */
  class R12EnergyIntermediates : virtual public SavableState {
    private:
      R12Technology::StandardApproximation stdapprox_;
      Ref<R12IntEval> r12eval_;
      bool V_computed_;
      RefSCMatrix V_[NSpinCases2];
      bool X_computed_;
      RefSymmSCMatrix X_[NSpinCases2];
      bool B_computed_;
      RefSymmSCMatrix B_[NSpinCases2];
      bool A_computed_;
      RefSCMatrix A_[NSpinCases2];
      bool T1_cc_computed_;
      RefSCMatrix T1_cc_[NSpinCases1];
      bool T2_cc_computed_;
      Ref<DistArray4> T2_cc_[NSpinCases2];
      // lambda amplitudes
      bool L1_cc_computed_;
      RefSCMatrix L1_cc_[NSpinCases1];
      bool L2_cc_computed_;
      Ref<DistArray4> L2_cc_[NSpinCases2];
      // parameters for importing psi ccsd one-particle density
      bool Onerdm_cc_computed_;
      RefSCMatrix Onerdm_cc_[NSpinCases1];
      // parameters for orbial relaxation contribution to
      // ccsd_f12 one-particle density
      bool Onerdm_relax_computed_;
      RefSCMatrix Onerdm_relax_[NSpinCases1];
    public:
      typedef enum { V=0, X=1, B=2, A=3 } IntermediateType;
      R12EnergyIntermediates(const Ref<R12IntEval>& r12eval,
                             const R12Technology::StandardApproximation stdapp);
      R12EnergyIntermediates(StateIn &si);
      void save_data_state(StateOut &so);
      ~R12EnergyIntermediates(){}
      Ref<R12IntEval> r12eval() const;
      void set_r12eval(Ref<R12IntEval> &r12eval);
      R12Technology::StandardApproximation stdapprox() const;
      bool V_computed() const;
      bool X_computed() const;
      bool B_computed() const;
      bool A_computed() const;
      bool T1_cc_computed() const;
      bool T2_cc_computed() const;
      bool L1_cc_computed() const;
      bool L2_cc_computed() const;
      bool Onerdm_cc_computed() const;
      bool Onerdm_relax_computed() const;
      const RefSCMatrix& get_V(const SpinCase2 &spincase2) const;
      void assign_V(const SpinCase2 &spincase2, const RefSCMatrix& V);
      const RefSymmSCMatrix& get_X(const SpinCase2 &spincase2) const;
      void assign_X(const SpinCase2 &spincase2, const RefSymmSCMatrix& X);
      const RefSymmSCMatrix& get_B(const SpinCase2 &spincase2) const;
      void assign_B(const SpinCase2 &spincase2, const RefSymmSCMatrix& B);
      const RefSCMatrix& get_A(const SpinCase2 &spincase2) const;
      void assign_A(const SpinCase2 &spincase2, const RefSCMatrix& A);
      const RefSCMatrix& get_T1_cc(const SpinCase1 &spincase1) const;
      void assign_T1_cc(const SpinCase1 &spincase1, const RefSCMatrix& T1_cc);
      const Ref<DistArray4>& get_T2_cc(const SpinCase2 &spincase2) const;
      void assign_T2_cc(const SpinCase2 &spincase2, const Ref<DistArray4>& T2_cc);
      const RefSCMatrix& get_L1_cc(const SpinCase1 &spincase1) const;
      void assign_L1_cc(const SpinCase1 &spincase1, const RefSCMatrix& L1_cc);
      const Ref<DistArray4>& get_L2_cc(const SpinCase2 &spincase2) const;
      void assign_L2_cc(const SpinCase2 &spincase2, const Ref<DistArray4>& L2_cc);
      const RefSCMatrix& get_1rdm_cc(const SpinCase1 &spincase1) const;
      void assign_1rdm_cc(const SpinCase1 &spincase1, const RefSCMatrix& Onerdm_cc);
      const RefSCMatrix& get_1rdm_relax(const SpinCase1 &spincase1) const;
      void assign_1rdm_relax(const SpinCase1 &spincase1, const RefSCMatrix& Onerdm_relax);
  };

  /** Class MP2R12Energy is the object that computes and maintains MP2-R12 energies */
class MP2R12Energy : virtual public SavableState {
  private:
    /** Computes values of all 2-body products from
        space1 and space2 if electron 1 is at r1 and
        electron 2 is at r2. equiv specifies whether electrons
        are equivalent (same spin) or not */
    RefSCVector compute_2body_values_(bool equiv, const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                      const SCVector3& r1, const SCVector3& r2) const;

    // computes f12 contributions to the mp2 pair energies
    virtual void compute_ef12() =0;

  protected:
    Ref<R12EnergyIntermediates> r12intermediates_;
    Ref<R12IntEval> r12eval_;
    bool include_obs_singles_;
    int debug_;
    bool evaluated_;

    RefSCVector ef12_[NSpinCases2], emp2f12_[NSpinCases2];
    // The coefficients are stored xy by ij, where xy is the geminal-multiplied pair
    RefSCMatrix C_[NSpinCases2];

    // Initialize SCVectors and SCMatrices
    void init();

  public:

    MP2R12Energy(StateIn&);
    MP2R12Energy(const Ref<R12EnergyIntermediates>& r12intermediates,
                 bool include_obs_singles,
                 int debug);
    ~MP2R12Energy();

    void save_data_state(StateOut&);
    void obsolete();
    void print(std::ostream&o=ExEnv::out0()) const;

    Ref<R12IntEval> r12eval() const;
    const Ref<R12EnergyIntermediates>& r12intermediates() const;
    R12Technology::StandardApproximation stdapprox() const;
    bool include_obs_singles() const { return include_obs_singles_; }
    void set_debug(int debug);
    int get_debug() const;

    /// total correlation energy (including OBS singles if include_obs_singles() == true)
    double energy();
    const RefSCVector& ef12(SpinCase2 S);
    const RefSCVector& emp2f12(SpinCase2 S);
    double emp2f12tot(SpinCase2 S);
    double ef12tot(SpinCase2 S);
    void print_pair_energies(bool spinadapted,
                             double cabs_singles_energy,
                             std::ostream&so=ExEnv::out0());

#if MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION
    /** Computes values of pair function ij on tbgrid */
    void compute_pair_function(unsigned int i, unsigned int j, SpinCase2 spincase2,
                               const Ref<TwoBodyGrid>& tbgrid);
#endif

    /// Computes the first-order R12 wave function and MP2-R12 energy
    void compute();
    /** Returns the matrix of first-order amplitudes of r12-multiplied occupied orbital pairs.
      */
    RefSCMatrix C(SpinCase2 S);
    /** Returns the matrix of first-order amplitudes of conventional orbital pairs.
     */
    RefSCMatrix T2(SpinCase2 S);
};

/**
 * The class MP2R12Energy_SpinOrbital is the original implementation of MP2R12Energy
 * It supports only the standard orbital-invariant ansatz and the full set of features
 * of R12Technology.
 */
class MP2R12Energy_SpinOrbital : public MP2R12Energy
{
    void compute_ef12();

  public:
    MP2R12Energy_SpinOrbital(StateIn&);
    MP2R12Energy_SpinOrbital(Ref<R12EnergyIntermediates> &r12intermediates,
                             bool include_obs_singles,
                             int debug);
    ~MP2R12Energy_SpinOrbital();

    void save_data_state(StateOut&);
};

/**
 * The class MP2R12Energy_Diag is an implementation of MP2R12Energy that supports
 * Ten-no's diagonal orbital-invariant ansatz for closed and open-shells.
 */
class MP2R12Energy_Diag : public MP2R12Energy
{
    void compute_ef12();
    /// this computes correct alpha-beta energies for open-shell cases
    void compute_ef12_10132011();
    void activate_ints(const std::string&, const std::string&,
                       const std::string&, const std::string&,
                       const std::string&, Ref<TwoBodyFourCenterMOIntsRuntime>&,
                       Ref<DistArray4>&);
    // Label Y^ij_ij, Y^ij_ji, Y^ji_ij, Y^ji_ji
    enum idx_b1b2k1k2{ij_ij, ij_ji, ji_ij, ji_ji};
    void compute_Y(const int b1b2_k1k2, const double prefactor,
                   const unsigned int oper_idx,
                   Ref<DistArray4>& i1i2i1i2_ints, double* array_i1i2i1i2);
    void compute_YxF(const int b1b2_k1k2, const double prefactor,
                     const unsigned int oper1_idx, const unsigned int oper2_idx,
                     const Ref<DistArray4>& i1i2x1x2_ints, const Ref<DistArray4>& i2i1x1x2_ints,
                     double* array_ijij);
    void compute_FxT(const int b1b2_k1k2, const unsigned int f12_idx,
                     Ref<DistArray4>& F_ints, const double* Tiiaa,
                     double* V_coupling);
    void compute_VX(const int b1b2_k1k2, std::vector<std::string>& VX_output,
                    const unsigned int oper12_idx, Ref<DistArray4>& Y12_ints,
                    const unsigned int oper1_idx, const unsigned int oper2_idx,
                    const std::vector<Ref<DistArray4> >& Y_ints,
                    const std::vector<Ref<DistArray4> >& F_ints,
                    double* VX_array);
    void accumulate_P_YxF(std::vector<std::string>& P_output,
                          std::vector<int>& b1b2_k1k2, std::vector<double>& P_prefactor,
                          const unsigned int oper1_idx, const unsigned int oper2_idx,
                          std::vector<Ref<DistArray4> >& Y_ints,
                          std::vector<Ref<DistArray4> >& F_ints,
                          double* P);
    /// Compute U=FxG intermediate needed for CC V bar
    void compute_U(Ref<DistArray4> (&U)[NSpinCases2]);
    void contract_VT1(const Ref<DistArray4>& V,
                         const int b1b2_k1k2,  const bool swap_e12_V,
                         const double* const T1_array,
                         const int nv, const bool VT1_offset,
                         double* const VT1);
    //
    // compute the one electron density matrix for the diagonal ansatz
    //RefSCMatrix D_ccsdf12_[NSpinCases1];
    void compute_density_diag();

    // functions needed for compute_density_diag() function:
    void obtain_orbitals(const SpinCase2 spincase,
                         std::vector<Ref<OrbitalSpace> >& v_orbs1,
                         std::vector<Ref<OrbitalSpace> >& v_orbs2);

    // activate three f12_ints for computing X
    void activate_ints_X_f12(Ref<TwoBodyFourCenterMOIntsRuntime>& moints4_rtime, const std::string& index,
                             const std::vector< Ref<OrbitalSpace> >& v_orbs1,
                             const std::vector< Ref<OrbitalSpace> >& v_orbs2,
                             const std::string& descr_f12_key, std::vector<Ref<DistArray4> >& f12_ints);

    // test function for D^i_i
    void compute_Dii_test(const int nspincases1, const int nspincases2,
                          const int nocc_alpha, const int nocc_beta,
                          const double C_0, const double C_1);
    // functions for compute_Dii_test:
    // compute AlphaAlpha/BetaBeta: R^ij_ab R^ab_ij & R^ji_ab R^ab_ij
    void compute_RRii_ii(std::vector<std::string>& output,
                         const std::vector< Ref<OrbitalSpace> >& v_orbs1,
                         const std::vector< Ref<OrbitalSpace> >& v_orbs2,
                         double* RRij_ij, double* RRji_ij);
    // compute the openshell AlphaBeta:
    // R^ij_ab R^ab_ij, R^ji_ab R^ab_ij, R^ij_ab R^ab_ji, & R^ji_ab R^ab_ji
    void compute_Rii_ii(std::vector<std::string>& output,
                        const std::vector< Ref<OrbitalSpace> >& v_orbs1,
                        const std::vector< Ref<OrbitalSpace> >& v_orbs2,
                        double* RRi1i2_i1i2, double* RRi1i2_i2i1,
                        double* RRi2i1_i1i2, double* RRi2i1_i2i1);

    // compute D^m_i
    void compute_Dmi(const int nspincases1, const int nspincases2,
                     const double C_0, const double C_1,
                     const std::vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                     const std::vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                     double* const Dm_i_alpha, double* const Dm_i_beta);
    // functions needed for computing D^m_i :
    // compute R^ab_b1b2 R^k1k2 _ab (a, b represent complete virtual orbitals)
    // sum over a, b, and index 3
    // index 1,2: i, m; index 3: j  e.g.: R^ab_ij R^jm_ab
    void compute_RR_sum_abj2(const int RRb1b2_k1k2,
                             const int f12f12_idx, const int f12_idx,
                             const Ref<DistArray4>& f12f12_ints,
                             const std::vector< Ref<DistArray4> >& v_f12_ints1,
                             const std::vector< Ref<DistArray4> >& v_f12_ints2,
                             double* const RR_result);

    // compute D^m_i through R^ij_a'b' R_ij^a'b' + 2 R^ij_ab' R_ij^ab'
    void compute_Dmi_2(const int nspincases1, const int nspincases2,
                       const double C_0, const double C_1,
                       const std::vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                       const std::vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                       double* const Dm_i_alpha, double* const Dm_i_beta);

    // test function fo D^c_b
    void compute_Dcb(const int nspincases1, const int nspincases2,
                     const double C_0, const double C_1);

    // test function for D^c'_b'
    void compute_Dcpbp_a(const int nspincases1, const int nspincases2,
                         const double C_0, const double C_1);

    // test function for D^c'_b' sum over a'
    void compute_Dcpbp_ap(const int nspincases1, const int nspincases2,
                          const double C_0, const double C_1);

    // test function for D^a'_a RR part
    void compute_Dapa_RR(const int nspincases1, const int nspincases2,
                      const double C_0, const double C_1);

    // \bar{\tilde{R}}^31_ij \bar{\tilde{R}}^ij_32 with spin orbitals
    // which is for computing D^b'_c', D^a'_a
    void compute_RR31_32_spin(const int orbitals_label,
                              const int nspincases1, const int nspincases2,
                              const double C_0, const double C_1,
                              const std::vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                              const std::vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                              double* const D_alhpha, double* const D_beta);

    // function need for D^a'_a:
    // compute R * T2 (CC amplitudes)
    void compute_RT2_apa(const int nspincases1, const int nspincases2,
                         const double C_0, const double C_1,
                         const std::vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                         const std::vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                         double* RT2_alpha, double* RT2_beta);

    // compute MP2 T2 amplitude (non-antisymmetrized)
    void compute_T2_mp2(const std::vector< Ref<OrbitalSpace> >& v_orbs1,
                        const std::vector< Ref<OrbitalSpace> >& v_orbs2,
                        double* T2ab_ij);
    // compute R * T2 (mp2 amplitudes)
    void compute_RTmp2_apa(const int nspincases1, const int nspincases2,
                           const double C_0, const double C_1,
                           const std::vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                           const std::vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                           double* const RT2_alpha, double* const RT2_beta);

    RefSCMatrix compute_D_CABS(SpinCase1 spin);
    RefSCMatrix compute_D_CABS_test(SpinCase1 spin);
    // compute MP2 one-electron density matrix
    RefSCMatrix compute_1rdm_mp2(const SpinCase1 spin);
    RefSCMatrix compute_1rdm_mp2_test(const SpinCase1 spin);
    // MP2F12 one-electron density matrix
    // D_MP2F12 = T(MP2)T(MP2) + 2 T(MP2)T(F12) + T(F12)T(F12)
    void compute_1rdm_mp2f12(const int nspincases1, const int nspincases2,
                             const int C_0, const int C_1,
                             RefSCMatrix Dmp2f12[NSpinCases1]);
    // compute contribution for MP2F12 one-electron density matrix: D_MP2F12
    // compute T(MP2)T(MP2) or T(F12)T(F12)
    RefSCMatrix compute_1rdm_mp2part(const SpinCase1 spin,
                                     const int nocc1_act, const int nocc2_act,
                                     const int nvir1, const int nvir2,
                                     const double* const T2, const double* const T2_ab);
    // compute T(MP2)T(F12) or T(MP2)T(F12)
    RefSCMatrix compute_1rdm_mp2part(const SpinCase1 spin,
                                     const int nocc1_act, const int nocc2_act,
                                     const int nvir1, const int nvir2,
                                     const double* const T2_left,const double* const T2_right,
                                     const double* const T2_ab_left, const double* const T2_ab_right);
    // compute MP2 T2 amplitude (antisymmetrized), stored in array (a,b,i,j)
    void compute_T2_mp2(const SpinCase2 spincase,
                        double* const T2ab_ij);
    // compute F12 corrected T2 amplitude (antisymmetrized), stored in array (a,b,i,j)
    void compute_T2abij_f12corr(const SpinCase2 spincase,
                                const double C_0, const double C_1,
                                double* const T2ab_ij_f12corr);
    // compute MP2F12 T2 amplitude = T2 (MP2) + T2 (F12 corrected)
    void compute_T2abij_mp2f12(const int nocc1_act, const int nocc2_act,
                               const int nvir1, const int nvir2,
                               const double* const T2ab_ij_mp2,
                               const double* const T2ab_ij_f12corr,
                               double* const T2abij_mp2f12);
    // compute MP2F12 T2 amplitude, stored in array (a,b,i,j)
    void compute_T2abij_mp2f12(const SpinCase2 spincase,
                               const double C_0, const double C_1,
                               double* const T2ab_ij_mp2f12);
    // transform the one-particle density matrix into the MPQC ordering
    RefSCMatrix onepdm_transformed(const SpinCase1& spin, const bool frozen_core, const RefSCMatrix& D);
    // test function for computing CCSD dipole moment
    RefSCMatrix onepdm_transformed2(const SpinCase1& spin,const RefSCMatrix& D);
    // compute D^a'_i = R * T1 (CC)
    void compute_RT1_api(const int nspincases1, const int nspincases2,
                         const double C_0, const double C_1,
                         const std::vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                         const std::vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                         double* const D_alpha, double* const D_beta);

  public:
    MP2R12Energy_Diag(StateIn&);
    MP2R12Energy_Diag(Ref<R12EnergyIntermediates> &r12intermediates,
                      bool include_obs_singles,
                      int debug);
    ~MP2R12Energy_Diag();

    void save_data_state(StateOut&);
};

Ref<MP2R12Energy> construct_MP2R12Energy(Ref<R12EnergyIntermediates> &r12intermediates,
                                         bool include_obs_singles,
                                         int debug,
                                         bool diag);

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


