//
// mp2r12_energy_diag2.cc
//
// Copyright (C) 2010 Jinmei Zhang
//
// Author: Jinmei Zhang <jmzhang@vt.edu> and Edward Valeev <www.valeyev.net>
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

#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <util/misc/print.h>
#include <math/scmat/blas.h>
#include <cmath>

using namespace std;
using namespace sc;

/**********************************
 * class MP2R12Energy_Diag *
 **********************************/

void MP2R12Energy_Diag::compute_ef12_10132011() {
  if (evaluated_)
    return;

  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
  Ref<MessageGrp> msg = r12world->world()->msg();
  int me = msg->me();
  int ntasks = msg->n();
  const bool do_mp2 = r12intermediates_->T1_cc_computed() == false &&
                      r12intermediates_->T2_cc_computed() == false;

  // test for 1e density matrix
  double E_V_tot = 0.0;
  double E_Vcoupling_tot = 0.0;
  double E_Bcoupling_tot = 0.0;
  double E_X_tot = 0.0;
  double E_B_tot = 0.0;
  double E_X_noca_tot = 0.0;
  double E_VT_tot = 0.0;

  const Ref<R12Technology::CorrelationFactor> corrfactor =
      r12world->r12tech()->corrfactor();
  const bool obs_eq_ribs = r12world->obs_eq_ribs();
  const bool cabs_empty = obs_eq_ribs;
  // Only diagonal ansatz is supported
  const bool diag = r12world->r12tech()->ansatz()->diag();
  if (diag == false)
    throw ProgrammingError("only diagonal ansatz supported", __FILE__,
                           __LINE__);

  // obtain some preliminaries
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime =
      r12world->world()->moints_runtime4();
  Ref<TwoBodyIntDescr> descr_f12 =
      r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(), 0);
  Ref<TwoBodyIntDescr> descr_f12f12 =
      r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(), 0, 0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const std::string descr_f12f12_key = moints4_rtime->descr_key(descr_f12f12);

  //
  // Evaluate pair energies:
  // distribute workload among nodes by pair index
  //
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  if (debug_ >= DefaultPrintThresholds::N2)
    ExEnv::out0() << endl << indent << "num_unique_spincases2 = "
        << num_unique_spincases2 << endl;

  // Find the largest dimension (nocc1_act*nocc2_act) for all the ijij or ijji arrays
  // start from getting the number of alpha and beta electron in active space
  int nalpha = 0;
  int nbeta = 0;
  int nocc_max = 0;
  int nijij = 0;

  // Find the largest dimension for ij_ab array: nocc1_act*nocc2_act*nvir1_act*nvir2_act
  int nvir_act_alpha = 0;
  int nvir_act_beta = 0;
  int nvir_max = 0;
  int nijab = 0;

  // Alpha-alpha case
  const SpinCase2 spincase_AlphAlpha = static_cast<SpinCase2>(1);
  if (r12eval()->dim_oo(spincase_AlphAlpha).n() != 0) {

    const SpinCase1 spin1_alpha = case1(spincase_AlphAlpha);
    const Ref<OrbitalSpace>& occ_act_alpha = r12eval()->occ_act(spin1_alpha);
    const Ref<OrbitalSpace>& vir_act_alpha = r12eval()->vir_act(spin1_alpha);
    nalpha = occ_act_alpha->rank();
    nvir_act_alpha = vir_act_alpha->rank();
  }
  // Beta-beta case
  const SpinCase2 spincase_BetaBeta = static_cast<SpinCase2>(2);
  if (r12eval()->dim_oo(spincase_BetaBeta).n() != 0) {

    const SpinCase1 spin1_beta = case1(spincase_BetaBeta);
    const Ref<OrbitalSpace>& occ_act_beta = r12eval()->occ_act(spin1_beta);
    const Ref<OrbitalSpace>& vir_act_beta = r12eval()->vir_act(spin1_beta);
    nbeta = occ_act_beta->rank();
    nvir_act_beta = vir_act_beta->rank();
  }
  // Alpha-beta: if there is only no alpha-alpha or beta-beta case
  if (r12eval()->dim_oo(spincase_AlphAlpha).n() == 0
      && r12eval()->dim_oo(spincase_BetaBeta).n() == 0) {
    const SpinCase2 spincase_AlphaBeta = static_cast<SpinCase2>(0);
    const SpinCase1 spin1_alpha = case1(spincase_AlphaBeta);
    const SpinCase1 spin2_beta = case2(spincase_AlphaBeta);
    const Ref<OrbitalSpace>& occ_act_alpha = r12eval()->occ_act(spin1_alpha);
    const Ref<OrbitalSpace>& occ_act_beta = r12eval()->occ_act(spin2_beta);
    const Ref<OrbitalSpace>& vir_act_alpha = r12eval()->vir_act(spin1_alpha);
    const Ref<OrbitalSpace>& vir_act_beta = r12eval()->vir_act(spin2_beta);
    nalpha = occ_act_alpha->rank();
    nbeta = occ_act_beta->rank();
    nvir_act_alpha = vir_act_alpha->rank();
    nvir_act_beta = vir_act_beta->rank();
  }
  if (debug_ >= DefaultPrintThresholds::mostN2)
    ExEnv::out0() << indent << "#n of spin1(active) = " << nalpha
        << ";  #n of spin2(active) = " << nbeta << endl;

  // Obtain the larger number between alpha an beta electron: N
  // The dimension of all the intermediate matrixes is: NxN
  // This is needed for the computation of coupled cluster V contribution
  if (nalpha >= nbeta)
  nocc_max = nalpha;
  else
  nocc_max = nbeta;

  nijij = nocc_max * nocc_max;
  if (debug_ >= DefaultPrintThresholds::mostN2)
  ExEnv::out0() << indent << "The dimension of V X and B matrix: " << nijij << endl;

  // Compute the dimension of T^ij_ab
  if (nvir_act_alpha >= nvir_act_beta)
  nvir_max = nvir_act_alpha;
  else
  nvir_max = nvir_act_beta;

  nijab = nijij * nvir_max * nvir_max;
  if (debug_ >= DefaultPrintThresholds::mostN2)
  ExEnv::out0() << indent << "The dimension of T^ij_ab matrix: " << nijab << endl;

  if (nijij == 0)
  ExEnv::out0() << indent << "error: no electron" << endl;

  //
  // Obtain CC amplitudes from Psi
  //
  RefSCMatrix T1[NSpinCases1];
  Ref<DistArray4> T2[NSpinCases2]; // need to be outside for coupling term
  if (r12intermediates_->T1_cc_computed()
      && r12intermediates_->T2_cc_computed()) {
    // Obtain T1 amplitudes
    const int nspincases1 = r12eval()->nspincases1();
    for (int s = 0; s < nspincases1; ++s) {
      const SpinCase1 spin = static_cast<SpinCase1>(s);
      T1[spin] = r12intermediates_->get_T1_cc(spin);
    }
    if (nspincases1 == 1) {
      T1[Beta] = T1[Alpha];
    }

    // Obtain T2 amplitudes
    for (int s = 0; s < num_unique_spincases2; ++s) {
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      if (r12eval()->dim_oo(spincase2).n() == 0)
        continue;
      T2[spincase2] = r12intermediates_->get_T2_cc(spincase2);
    } // end of T2
  } // end of CC amplitudes

  //
  // geminal coefficients
  const double C_0 = 1.0 / 2.0;
  const double C_1 = 1.0 / 4.0;

  std::vector<double> ef12_VT[NSpinCases2];
  std::vector<double> ef12_TBT[NSpinCases2];
  for (int spin = 0; spin < num_unique_spincases2; spin++) {
    const SpinCase2 spincase = static_cast<SpinCase2>(spin);
    const int noo = r12eval()->dim_oo(spincase).n();
    ef12_VT[spin].resize(noo);
    ef12_TBT[spin].resize(noo);
  }

  //
  // Compute the MP2-R12 part
  //
  for (int spin = 0; spin < num_unique_spincases2; spin++) {

    const SpinCase2 spincase = static_cast<SpinCase2>(spin);
    if (r12eval()->dim_oo(spincase).n() == 0)
      continue;
    const SpinCase1 spin1 = case1(spincase);
    const SpinCase1 spin2 = case2(spincase);

    // Determine the spincase
    std::string spinletters = to_string(spincase);

    const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
    const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
    const Ref<OrbitalSpace>& vir1_act = r12eval()->vir_act(spin1);
    const Ref<OrbitalSpace>& vir2_act = r12eval()->vir_act(spin2);
    const int nocc1_act = occ1_act->rank();
    const int nocc2_act = occ2_act->rank();
    const int nvir1_act = vir1_act->rank();
    const int nvir2_act = vir2_act->rank();
    const int nocc12 = nocc1_act * nocc2_act;
    if (nocc12 == 0)
      continue; // skip spincase if no electron pairs of this kind
    const int nvir12_act = nvir1_act * nvir2_act;

    const Ref<OrbitalSpace>& occ1 = r12eval()->occ(spin1);
    const Ref<OrbitalSpace>& vir1 = r12eval()->vir(spin1);
    const Ref<OrbitalSpace>& orbs1 = r12eval()->orbs(spin1);
    const Ref<OrbitalSpace>& cabs1 = r12world->cabs_space(spin1);
    const Ref<OrbitalSpace>& ribs1 = r12world->ribs_space();
    const Ref<OrbitalSpace>& Kribs1 = r12eval()->K_P_P(spin1);
    const Ref<OrbitalSpace>& Fribs1 = r12eval()->F_P_P(spin1);
    const Ref<OrbitalSpace>& forbsp1 = r12eval()->F_p_p(spin1);
    const Ref<OrbitalSpace>& Focc1 = r12eval()->F_m_P(spin1);
    const Ref<OrbitalSpace>& focc1 = r12eval()->F_m_m(spin1);
    const Ref<OrbitalSpace>& forbsA1 = r12eval()->F_p_A(spin1);
    const Ref<OrbitalSpace>& fvir1_act = r12eval()->F_a_A(spin1);

    const Ref<OrbitalSpace>& occ2 = r12eval()->occ(spin2);
    const Ref<OrbitalSpace>& vir2 = r12eval()->vir(spin2);
    const Ref<OrbitalSpace>& orbs2 = r12eval()->orbs(spin2);
    const Ref<OrbitalSpace>& cabs2 = r12world->cabs_space(spin2);
    const Ref<OrbitalSpace>& ribs2 = r12world->ribs_space();
    const Ref<OrbitalSpace>& Kribs2 = r12eval()->K_P_P(spin2);
    const Ref<OrbitalSpace>& Fribs2 = r12eval()->F_P_P(spin2);
    const Ref<OrbitalSpace>& forbsp2 = r12eval()->F_p_p(spin2);
    const Ref<OrbitalSpace>& Focc2 = r12eval()->F_m_P(spin2);
    const Ref<OrbitalSpace>& focc2 = r12eval()->F_m_m(spin2);
    const Ref<OrbitalSpace>& forbsA2 = r12eval()->F_p_A(spin2);
    const Ref<OrbitalSpace>& fvir2_act = r12eval()->F_a_A(spin2);

    //
    // Allocate storage for V, X, and B intermediates
    //
    double* Vij_ij = new double[nocc12];
    fill_n(Vij_ij, nocc12, 0.0);
    double* Vij_ji = new double[nocc12];
    fill_n(Vij_ji, nocc12, 0.0);

    // Assign pointers to V coupling matrixes
    double* Vij_ij_coupling = NULL;
    double* Vij_ji_coupling = NULL;
    //double* Vji_ij_coupling = NULL;
    //double* Vji_ji_coupling = NULL;
    double* Tij_ab = NULL;
    if (this->r12eval()->coupling() == true
        || this->r12eval()->ebc() == false) {
      Vij_ij_coupling = new double[nocc12];
      fill_n(Vij_ij_coupling, nocc12, 0.0);
      Vij_ji_coupling = new double[nocc12];
      fill_n(Vij_ji_coupling, nocc12, 0.0);
      Tij_ab = new double[nocc12 * nvir12_act];
    }
    // in MP2 with coupling=true also need to compute a coupling contribution to B matrix
    double* Bij_ij_coupling = NULL;
    double* Bij_ji_coupling = NULL;
    double* Bji_ij_coupling = NULL;
    double* Bji_ji_coupling = NULL;
    if (this->r12eval()->coupling() == true && do_mp2) {
      Bij_ij_coupling = new double[nocc12];
      fill_n(Bij_ij_coupling, nocc12, 0.0);
      Bij_ji_coupling = new double[nocc12];
      fill_n(Bij_ji_coupling, nocc12, 0.0);
      Bji_ij_coupling = new double[nocc12];
      fill_n(Bji_ij_coupling, nocc12, 0.0);
      Bji_ji_coupling = new double[nocc12];
      fill_n(Bji_ji_coupling, nocc12, 0.0);
    }

    // stored as o1 x o2 matrix
    double* Xij_ij = new double[nocc12];
    fill_n(Xij_ij, nocc12, 0.0);
    double* Xij_ji = new double[nocc12];
    fill_n(Xij_ji, nocc12, 0.0);
    double* Bij_ij = new double[nocc12];
    fill_n(Bij_ij, nocc12, 0.0);
    double* Bij_ji = new double[nocc12];
    fill_n(Bij_ji, nocc12, 0.0);
    double* Pij_ij = new double[nocc12];
    fill_n(Pij_ij, nocc12, 0.0);
    double* Pij_ji = new double[nocc12];
    fill_n(Pij_ji, nocc12, 0.0);
    double* Qij_ij = new double[nocc12];
    fill_n(Qij_ij, nocc12, 0.0);
    double* Qij_ji = new double[nocc12];
    fill_n(Qij_ji, nocc12, 0.0);
    // these are only needed in alpha-beta open-shell
    // stored as o1 x o2 matrix
    double* Xji_ji = 0;
    double* Xji_ij = 0;
    // stored as o2 x o1 matrix
    double* Bji_ji = 0;
    double* Pji_ji = 0;
    double* Bji_ij = 0;
    double* Pji_ij = 0;
    if (spin1 != spin2 && num_unique_spincases2 == 3) {
      Xji_ji = new double[nocc12];
      fill_n(Xji_ji, nocc12, 0.0);
      Xji_ij = new double[nocc12];
      fill_n(Xji_ij, nocc12, 0.0);
      Bji_ji = new double[nocc12];
      fill_n(Bji_ji, nocc12, 0.0);
      Bji_ij = new double[nocc12];
      fill_n(Bji_ij, nocc12, 0.0);
      Pji_ji = new double[nocc12];
      fill_n(Pji_ji, nocc12, 0.0);
      Pji_ij = new double[nocc12];
      fill_n(Pji_ij, nocc12, 0.0);
    }

    //
    // compute intermediates V, X, B
    //
    const TwoBodyOper::type f12eri_type =
        r12world->r12tech()->corrfactor()->tbint_type_f12eri();
    const TwoBodyOper::type f12_type =
        r12world->r12tech()->corrfactor()->tbint_type_f12();
    const TwoBodyOper::type eri_type =
        r12world->r12tech()->corrfactor()->tbint_type_eri();
    const TwoBodyOper::type f12f12_type =
        r12world->r12tech()->corrfactor()->tbint_type_f12f12();
    const TwoBodyOper::type f12t1f12_type =
        r12world->r12tech()->corrfactor()->tbint_type_f12t1f12();
    const unsigned int f12eri_idx = descr_f12->intset(f12eri_type);
    const unsigned int f12_idx = descr_f12->intset(f12_type);
    const unsigned int eri_idx = descr_f12->intset(eri_type);
    const unsigned int f12f12_idx = descr_f12f12->intset(f12f12_type);
    const unsigned int f12t1f12_idx = descr_f12f12->intset(f12t1f12_type);

    // Get eigenvalues of Fock matrix
#define COMPUTE_ORBITALSPACE_EIGENVALUES 0
#if COMPUTE_ORBITALSPACE_EIGENVALUES
    RefDiagSCMatrix evals_i1;
    RefDiagSCMatrix evals_i2;
    {
      RefSCMatrix F_ii_1 = r12eval()->fock(occ1_act, occ1_act, spin1);
      evals_i1 = F_ii_1.kit()->diagmatrix(F_ii_1.rowdim());
      for(unsigned int o=0; o<nocc1_act; ++o) evals_i1.set_element(o, F_ii_1(o, o));
    }
    if (occ1_act != occ2_act) {
      RefSCMatrix F_ii_2 = r12eval()->fock(occ2_act, occ2_act, spin2);
      evals_i2 = F_ii_2.kit()->diagmatrix(F_ii_2.rowdim());
      for(unsigned int o=0; o<nocc2_act; ++o) evals_i2.set_element(o, F_ii_2(o, o));
    }
    else
    evals_i2 = evals_i1;
#else
    const RefDiagSCMatrix evals_i1 = occ1_act->evals();
    const RefDiagSCMatrix evals_i2 = occ2_act->evals();
#endif

    //
    // Compute the V intermediate matrix: V^ij_ij and V^ij_ji
    //
    // Alpha_beta V^ij_ij and V^ij_ji are antisymmetrized
    // Alpha_alpha or beta_beta antisymmetrized V = V^ij_ij - V^ij_ji
    // Note: '^' indicates bra, '_' indicates ket
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " V^ij_ij: " << endl;

    // this is necessary because how I wrote the compute_Y function
    fill_n(Vij_ij, nocc12, 0.0);

    // V^ij_ij = (f12/r12)^ij_ij -
    Ref<DistArray4> i1i2i1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), occ1_act->id(),
                  occ2_act->id(), descr_f12_key, moints4_rtime, i1i2i1i2_ints);

    // store all the ints
    std::vector<Ref<DistArray4> > f12_ij_ints;
    std::vector<std::string> VX_output;

    // V^ij_ij -= g^ij_pq f^pq_ij
    Ref<DistArray4> i1i2p1p2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), orbs1->id(), orbs2->id(),
                  descr_f12_key, moints4_rtime, i1i2p1p2_ints);
    f12_ij_ints.push_back(i1i2p1p2_ints);
    VX_output.push_back("diag-pq contribution");

    // V^ij_ij -= g^ij_ma' f^ma'_ij
    Ref<DistArray4> i1i2i1a2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), occ1->id(), cabs2->id(),
                  descr_f12_key, moints4_rtime, i1i2i1a2_ints);

    f12_ij_ints.push_back(i1i2i1a2_ints);
    VX_output.push_back("diag-pq-ma' contribution");

    // V^ij_ij -= g^ij_a'm f^a'm_ij
    // TODO: for RHF can simply scale previous contribution by 2 and "symmetrize" at the end V^ij_ij = 0.5 * (V^ij_ij + V^ji_ji)
    Ref<DistArray4> i1i2a1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), cabs1->id(), occ2->id(),
                  descr_f12_key, moints4_rtime, i1i2a1i2_ints);

    f12_ij_ints.push_back(i1i2a1i2_ints);
    VX_output.push_back("diag-pq-ma'-a'm contribution");

    compute_VX(ij_ij, VX_output, f12eri_idx, i1i2i1i2_ints, eri_idx, f12_idx,
               f12_ij_ints, f12_ij_ints, Vij_ij);

    //
    // Vij_ji = V^ij_ji = g^ij_ab f^ab_ji
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " V^ij_ji : " << endl;

    fill_n(Vij_ji, nocc12, 0.0);

    // V^ij_ji = (f12/r12)^ij_ji
    Ref<DistArray4> i1i2i2i1_ints;

    // V^ij_ji -= g^ij_pq f^pq_ji
    Ref<DistArray4> i2i1p1p2_ints;

    // V^ij_ji -= r^ij_ma' f^ma'_ji
    Ref<DistArray4> i2i1i1a2_ints;

    // V^ij_ji -= r^ij_a'm f^a'm_ji
    Ref<DistArray4> i2i1a1i2_ints;

    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ1_act->id(), occ2_act->id(), occ2_act->id(),
                    occ1_act->id(), descr_f12_key, moints4_rtime,
                    i1i2i2i1_ints);

      activate_ints(occ2_act->id(), occ1_act->id(), orbs1->id(), orbs2->id(),
                    descr_f12_key, moints4_rtime, i2i1p1p2_ints);

      activate_ints(occ2_act->id(), occ1_act->id(), occ1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime, i2i1i1a2_ints);

      activate_ints(occ2_act->id(), occ1_act->id(), cabs1->id(), occ2->id(),
                    descr_f12_key, moints4_rtime, i2i1a1i2_ints);
    } else {
      i1i2i2i1_ints = i1i2i1i2_ints;
      i2i1p1p2_ints = i1i2p1p2_ints;
      i2i1i1a2_ints = i1i2i1a2_ints;
      i2i1a1i2_ints = i1i2a1i2_ints;
    }

    std::vector<Ref<DistArray4> > f12_ji_ints;
    f12_ji_ints.push_back(i2i1p1p2_ints);
    f12_ji_ints.push_back(i2i1i1a2_ints);
    f12_ji_ints.push_back(i2i1a1i2_ints);

    compute_VX(ij_ji, VX_output, f12eri_idx, i1i2i2i1_ints, eri_idx, f12_idx,
               f12_ij_ints, f12_ji_ints, Vij_ji);

    if (debug_ >= DefaultPrintThresholds::N2 && spin1 == spin2)
      print_antisym_intermediate(spincase, "V^ij_ij", Vij_ij, Vij_ji, nocc1_act,
                                 nocc2_act);

    if (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2) {
      // Alpha-beta case
      print_intermediate(spincase, "V^ij_ij", Vij_ij, nocc1_act, nocc2_act);
      print_intermediate(spincase, "V^ij_ji", Vij_ji, nocc1_act, nocc2_act);
    }

    // Compute V coupling
    // Alpha-alpha or beta-beta:
    // sum(i<j,ab,a') C1* (R^ij_aa' f^a'_b T^ab_ij + R^ji_aa' f^a'_b T^ab_ji
    //                   - R^ji_aa' f^a'_b T^ab_ij - R^ij_aa' f^a'_b T^ab_ji)
    //
    // Alpha-beta: i:alpha, j:beta, a:alpha, b:beta
    // sum(i<j,a<b,a') (C0+C1)/2* (R^ij_a'(alpha)b f^a'_a T^ab_ij + R^ij_aa'(beta) f^a'_b T^ab_ij)
    //               + (C0-C1)/2* (R^ji_a'(alpha)b f^a'_a T^ab_ij + R^ji_aa'(beta) f^a'_b T^ab_ij)

    if (this->r12eval()->coupling() == true
        || this->r12eval()->ebc() == false) {
      fill_n(Vij_ij_coupling, nocc12, 0.0);
      fill_n(Vij_ji_coupling, nocc12, 0.0);

      // Get eigenvalues of Fock matrix
#if COMPUTE_ORBITALSPACE_EIGENVALUES
      RefDiagSCMatrix evals_a1;
      RefDiagSCMatrix evals_a2;
      {
        RefSCMatrix F_aa_1 = r12eval()->fock(vir1_act, vir1_act, spin1);
        evals_a1 = F_aa_1.kit()->diagmatrix(F_aa_1.rowdim());
        for(unsigned int o=0; o<nvir1_act; ++o) evals_a1.set_element(o, F_aa_1(o, o));
      }
      if (vir1_act != vir2_act) {
        RefSCMatrix F_aa_2 = r12eval()->fock(vir2_act, vir2_act, spin2);
        evals_a2 = F_aa_2.kit()->diagmatrix(F_aa_2.rowdim());
        for(unsigned int o=0; o<nvir2_act; ++o) evals_a2.set_element(o, F_aa_2(o, o));
      }
      else
      evals_a2 = evals_a1;
#else
      const RefDiagSCMatrix evals_a1 = vir1_act->evals();
      const RefDiagSCMatrix evals_a2 = vir2_act->evals();
#endif

      // Print out all the eigenvalues
#if 0
      if (debug_ >= DefaultPrintThresholds::mostN2) {
        ExEnv::out0() << endl << indent << "evals_i1: " << endl;
        for (int i1 = 0; i1 < nocc1_act; ++i1) {
          const double e_i1 = evals_i1(i1);
          ExEnv::out0() << indent << scprintf("%12.10f", e_i1) << endl;
        }

        ExEnv::out0() << endl << indent <<  "evals_i2: " << endl;
        for (int i2 = 0; i2 < nocc2_act; ++i2) {
          const double e_i2 = evals_i2(i2);
          ExEnv::out0() << indent << scprintf("%12.10f", e_i2) << endl;
        }

        ExEnv::out0() << endl << indent << "evals_a1: " << endl;
        for (int a1 = 0; a1 < nvir1_act; ++a1) {
          ExEnv::out0() << indent << evals_a1(a1) << endl;
        }

        ExEnv::out0() << endl << indent << "evals_a2: " << endl;
        for (int a2 = 0; a2 < nvir2_act; ++a2) {
          ExEnv::out0() << indent << evals_a2(a2) << endl;
        }
      }
#endif

      // Initialize all the integrals needed
      // Vij_ij_coupling: R^ij_aa' f^a'_b T^ab_ij
      Ref<DistArray4> i1i2a1AF2_ints = NULL;
      activate_ints(occ1_act->id(), occ2_act->id(), vir1_act->id(),
                    fvir2_act->id(), descr_f12_key, moints4_rtime,
                    i1i2a1AF2_ints);

#if 0
      /// print C, the bare coupling matrix
      {

        RefSCMatrix abF_mat = SCMatrixKit::default_matrixkit()->matrix(new SCDimension(nvir1_act),
            new SCDimension(nvir2_act));
        RefSCMatrix bFa_mat = SCMatrixKit::default_matrixkit()->matrix(new SCDimension(nvir2_act),
            new SCDimension(nvir1_act));
        for(int i1=0; i1<nocc1_act; ++i1) {
          for(int i2=0; i2<nocc2_act; ++i2) {
            const double* abF_blk = i1i2a1AF2_ints->retrieve_pair_block(i1, i2, f12_idx);
            const double* bFa_blk = i1i2a1AF2_ints->retrieve_pair_block(i2, i1, f12_idx);
            abF_mat.assign(abF_blk);
            bFa_mat.assign(bFa_blk);

            std::ostringstream oss;
            oss << "<i j | a b_F> + <i j | b_F a> integral (i = " << i1 << ", j = " << i2 << ")" << std::endl;
            (abF_mat + bFa_mat.t()).print(oss.str().c_str());

            i1i2a1AF2_ints->release_pair_block(i1, i2, f12_idx);
            i1i2a1AF2_ints->release_pair_block(i2, i1, f12_idx);
          }
        }

        i1i2a1a2_ints->deactivate();
      }
#endif

      // Vji_ji_coupling: R^ij_a'b f^a'_a T^ab_ij
      Ref<DistArray4> i1i2AF1a2_ints = NULL;
//      if (spin1 != spin2) {
        activate_ints(occ1_act->id(), occ2_act->id(), fvir1_act->id(),
                      vir2_act->id(), descr_f12_key, moints4_rtime,
                      i1i2AF1a2_ints);
//      }

      Ref<DistArray4> i2i1AF1a2_ints = NULL;
      Ref<DistArray4> i2i1a1AF2_ints = NULL;
      if (num_unique_spincases2 == 3 && spin1 != spin2) {

        // i:alpha  j:beta   a:alpha b:beta
        // Vij_ji_coupling: R^ji_a'b f^a'_a T^ab_ij
        // a':alpha
        activate_ints(occ2_act->id(), occ1_act->id(), fvir1_act->id(),
                      vir2_act->id(), descr_f12_key, moints4_rtime,
                      i2i1AF1a2_ints);

        // Vji_ij_coupling: R^ji_aa' f^a'_b T^ab_ij
        // a':beta
        activate_ints(occ2_act->id(), occ1_act->id(), vir1_act->id(),
                      fvir2_act->id(), descr_f12_key, moints4_rtime,
                      i2i1a1AF2_ints);
      } else {
        i2i1AF1a2_ints = i1i2AF1a2_ints;
        i2i1a1AF2_ints = i1i2a1AF2_ints;
      }

      if (r12intermediates_->T2_cc_computed()) {
        if (debug_ >= DefaultPrintThresholds::N2)
          ExEnv::out0() << endl << indent << "Coupled-cluster V coupling:"
              << endl;
        T2[spin]->activate();
        if (spin1 != spin2) {
          compute_YxF(ij_ij, 1.0, f12_idx, 0, i1i2a1AF2_ints, T2[spin],
                      Vij_ij_coupling);
          compute_YxF(ij_ij, 1.0, f12_idx, 0, i1i2AF1a2_ints, T2[spin],
                      Vij_ij_coupling);

          compute_YxF(ji_ij, 1.0, f12_idx, 0, i2i1AF1a2_ints, T2[spin],
                      Vij_ji_coupling);
          compute_YxF(ji_ij, 1.0, f12_idx, 0, i2i1a1AF2_ints, T2[spin],
                      Vij_ji_coupling);

        } else {
          compute_YxF(ij_ij, 1.0, f12_idx, 0, i1i2a1AF2_ints, T2[spin],
                      Vij_ij_coupling);
          compute_YxF(ji_ij, 1.0, f12_idx, 0, i1i2a1AF2_ints, T2[spin],
                      Vij_ji_coupling);
        }
        T2[spin]->deactivate();
      } // done with CC coupling contribution to V
      else if (do_mp2 && this->r12eval()->coupling() == true){ // MP2 coupling contibutions to V and B

        // Start computing MP2 V coupling
        if (debug_ >= DefaultPrintThresholds::N2)
          ExEnv::out0() << endl << indent << spinletters << "  MP2 V coupling:"
              << endl;
        // Compute T^ij_ab (not antisymmetrized)
        // g^ij_ab
        Ref<DistArray4> i1i2a1a2_ints;
        activate_ints(occ1_act->id(), occ2_act->id(), vir1_act->id(),
                      vir2_act->id(), descr_f12_key, moints4_rtime,
                      i1i2a1a2_ints);

        //   T^i(spin1)j(spin2)_a(spin1)b(spin2) 4-dimension matrix
        // = g^ij_ab / (e_i + e_j - e_a -e_b)
        double* iter_Tij_ab = Tij_ab;
        for (int i1 = 0; i1 < nocc1_act; ++i1) {
          for (int i2 = 0; i2 < nocc2_act; ++i2) {
            const double* gij_ab = i1i2a1a2_ints->retrieve_pair_block(i1, i2,
                                                                      eri_idx);
            double mp2_pair_energy = 0;

            for (int a1 = 0; a1 < nvir1_act; ++a1) {
              for (int a2 = 0; a2 < nvir2_act; ++a2, ++iter_Tij_ab, ++gij_ab) {
                *iter_Tij_ab =
                    (*gij_ab)
                        / (evals_i1(i1) + evals_i2(i2) - evals_a1(a1)
                            - evals_a2(a2));
              }
            }
            i1i2a1a2_ints->release_pair_block(i1, i2, eri_idx);
          }
        }
        i1i2a1a2_ints->deactivate();

        // Alpha-alpha, beta-beta and alpha-beta for closed shell:
        // Tji_ab matrix <=> Tij_ab matrix

        if (debug_ >= DefaultPrintThresholds::mostN2)
          ExEnv::out0() << endl << indent << spinletters
              << " Vij_ij_coupling: R^ij_aa' f^a'_b T^ab_ij" << endl;
        compute_FxT(ij_ij, f12_idx, i1i2a1AF2_ints, Tij_ab, Vij_ij_coupling);

        if (spin1 != spin2) {

          if (debug_ >= DefaultPrintThresholds::mostN2)
            ExEnv::out0() << indent << spinletters
                << " Vij_ij_coupling(a':alpha): + R^ij_a'b f^a'_a T^ab_ij"
                << endl;
          compute_FxT(ij_ij, f12_idx, i1i2AF1a2_ints, Tij_ab, Vij_ij_coupling);

          if (debug_ >= DefaultPrintThresholds::mostN2)
            ExEnv::out0() << indent << spinletters
                << " Vji_ij_coupling: R^ji_a'b f^a'_a T^ab_ij" << endl;
          compute_FxT(ji_ij, f12_idx, i2i1AF1a2_ints, Tij_ab, Vij_ji_coupling);

          if (debug_ >= DefaultPrintThresholds::mostN2)
            ExEnv::out0() << indent << spinletters
                << " Vji_ij_coupling(a':beta): + R^ji_aa' f^a'_b T^ab_ij"
                << endl;
          compute_FxT(ji_ij, f12_idx, i2i1a1AF2_ints, Tij_ab, Vij_ji_coupling);
        } else {
          // + R^ji_aa' f^a'_b T^ab_ji
          if (debug_ >= DefaultPrintThresholds::mostN2)
            ExEnv::out0() << indent << spinletters
                << " Vji_ji_coupling: + R^ji_aa' f^a'_b T^ab_ji" << endl;
          compute_FxT(ji_ji, f12_idx, i1i2a1AF2_ints, Tij_ab, Vij_ij_coupling);
          // - R^ij_aa' f^a'_b T^ab_ji
          if (debug_ >= DefaultPrintThresholds::mostN2)
            ExEnv::out0() << indent << spinletters
                << " Vij_ji_coupling: R^ij_aa' f^a'_b T^ab_ji" << endl;
          compute_FxT(ij_ji, f12_idx, i1i2a1AF2_ints, Tij_ab, Vij_ji_coupling);
          // - R^ji_aa' f^a'_b T^ab_ij
          if (debug_ >= DefaultPrintThresholds::mostN2)
            ExEnv::out0() << indent << spinletters
                << " Vji_ij_coupling: + R^ji_aa' f^a'_b T^ab_ij" << endl;
          compute_FxT(ji_ij, f12_idx, i1i2a1AF2_ints, Tij_ab, Vij_ji_coupling);
        }

        // Start computing MP2 B coupling
        //  1) compute A^i(spin1)j(spin2)_a(spin1)b(spin2)
        //     = C^ij_ab / (e_i + e_j - e_a -e_b)
        //  2) B_coupling = A^ij_ab C^ab_ij
        if (debug_ >= DefaultPrintThresholds::N2)
          ExEnv::out0() << endl << indent << spinletters << "  MP2 B coupling:"
              << endl;
        // compute A = - C / denom
        Ref<DistArray4> A_i1i2a1a2 = i1i2a1AF2_ints->clone();
        A_i1i2a1a2->activate();
        double* A_buf = new double[nvir12_act];
        for (int i1 = 0; i1 < nocc1_act; ++i1) {
          for (int i2 = 0; i2 < nocc2_act; ++i2) {
            const double* Cij_ab1 = i1i2a1AF2_ints->retrieve_pair_block(i1, i2,
                                                                        f12_idx);
            const double* Cij_ab2 = i1i2AF1a2_ints->retrieve_pair_block(i1, i2,
                                                                        f12_idx);

            for (int a1 = 0, a12 = 0; a1 < nvir1_act; ++a1) {
              for (int a2 = 0; a2 < nvir2_act; ++a2, ++a12) {
                A_buf[a12] =
                    (Cij_ab1[a12] + Cij_ab2[a12])
                        / (evals_i1(i1) + evals_i2(i2) - evals_a1(a1)
                            - evals_a2(a2));
              }
            }
            A_i1i2a1a2->store_pair_block(i1, i2, f12_idx, A_buf);
            i1i2a1AF2_ints->release_pair_block(i1, i2, f12_idx);
            i1i2AF1a2_ints->release_pair_block(i1, i2, f12_idx);
          }
        }

        compute_YxF(ij_ij, 1.0, f12_idx, f12_idx, A_i1i2a1a2, i1i2a1AF2_ints, Bij_ij_coupling);
        compute_YxF(ij_ij, 1.0, f12_idx, f12_idx, A_i1i2a1a2, i1i2AF1a2_ints, Bij_ij_coupling);

        compute_YxF(ij_ji, 1.0, f12_idx, f12_idx, A_i1i2a1a2, i2i1AF1a2_ints, Bij_ji_coupling);
        compute_YxF(ij_ji, 1.0, f12_idx, f12_idx, A_i1i2a1a2, i2i1a1AF2_ints, Bij_ji_coupling);

        A_i1i2a1a2->deactivate();
        A_i1i2a1a2 = 0;

        if (num_unique_spincases2 == 3 && spin1 != spin2) {

          // compute A = - C / denom
          Ref<DistArray4> A_i2i1a1a2 = i2i1a1AF2_ints->clone();
          A_i2i1a1a2->activate();
          double* A_buf = new double[nvir12_act];
          for (int i2 = 0; i2 < nocc2_act; ++i2) {
            for (int i1 = 0; i1 < nocc1_act; ++i1) {
              const double* Cji_ab1 = i2i1a1AF2_ints->retrieve_pair_block(i2, i1,
                                                                          f12_idx);
              const double* Cji_ab2 = i2i1AF1a2_ints->retrieve_pair_block(i2, i1,
                                                                          f12_idx);

              for (int a1 = 0, a12 = 0; a1 < nvir1_act; ++a1) {
                for (int a2 = 0; a2 < nvir2_act; ++a2, ++a12) {
                  A_buf[a12] =
                      (Cji_ab1[a12] + Cji_ab2[a12])
                          / (evals_i1(i1) + evals_i2(i2) - evals_a1(a1)
                              - evals_a2(a2));
                }
              }
              A_i2i1a1a2->store_pair_block(i2, i1, f12_idx, A_buf);
              i2i1a1AF2_ints->release_pair_block(i2, i1, f12_idx);
              i2i1AF1a2_ints->release_pair_block(i2, i1, f12_idx);
            }
          }

          compute_YxF(ji_ji, 1.0, f12_idx, f12_idx, A_i2i1a1a2, i2i1AF1a2_ints, Bji_ji_coupling);
          compute_YxF(ji_ji, 1.0, f12_idx, f12_idx, A_i2i1a1a2, i2i1a1AF2_ints, Bji_ji_coupling);

          compute_YxF(ji_ij, 1.0, f12_idx, f12_idx, A_i2i1a1a2, i1i2AF1a2_ints, Bji_ij_coupling);
          compute_YxF(ji_ij, 1.0, f12_idx, f12_idx, A_i2i1a1a2, i1i2a1AF2_ints, Bji_ij_coupling);

          A_i2i1a1a2->deactivate();
          A_i2i1a1a2 = 0;
        }

      } // end of MP2 V & B coupling computation
      i1i2a1AF2_ints->deactivate();
      i1i2AF1a2_ints->deactivate();
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        i2i1AF1a2_ints->deactivate();
        i2i1a1AF2_ints->deactivate();
      }

      if (debug_ >= DefaultPrintThresholds::N2) {
        if (spin1 == spin2) {
          print_antisym_intermediate(spincase, "V^ij_ij coupling",
                                     Vij_ij_coupling, Vij_ji_coupling,
                                     nocc1_act, nocc2_act);
          if (do_mp2 && this->r12eval()->coupling())
            print_antisym_intermediate(spincase, "B^ij_ij coupling",
                                       Bij_ij_coupling, Bij_ji_coupling,
                                       nocc1_act, nocc2_act);
        } else {
          print_intermediate(spincase, "V^ij_ij coupling", Vij_ij_coupling,
                             nocc1_act, nocc2_act);
          print_intermediate(spincase, "V^ij_ji coupling", Vij_ji_coupling,
                             nocc1_act, nocc2_act);

          if (do_mp2 && this->r12eval()->coupling()) {
            print_intermediate(spincase, "B^ij_ij coupling", Bij_ij_coupling, nocc1_act, nocc2_act);
            print_intermediate(spincase, "B^ij_ji coupling", Bij_ji_coupling, nocc1_act, nocc2_act);
            if (spin1 != spin2 && num_unique_spincases2 == 3) {
              print_intermediate(spincase, "B^ji_ij coupling", Bji_ij_coupling, nocc1_act, nocc2_act);
              print_intermediate(spincase, "B^ji_ji coupling", Bji_ji_coupling, nocc1_act, nocc2_act);
            }
          }
        }
      } // end of debug

    } // end of V coupling computation

    //
    // Compute X intermediate matrix: X^ij_ij, X^ij_ji, X^ji_ij, X^ji_ji (the latter two only needed for open-shell alpha-beta)
    //
    // Alpha_beta X^ij_ij and X^ij_ji are not antisymmetrized
    // Alpha_alpha or beta_beta antisymmetrized X = X^ij_ij - X^ij_ji

    // Xij_ij = X^ij_ij = f^ij_ab f^ab_ij
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " X^ij_ij : " << endl;
    fill_n(Xij_ij, nocc12, 0.0);

    // X^ij_ij += (f12f12)^ij_ij
    Ref<DistArray4> i1i2i1i2_f12f12_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), occ1_act->id(),
                  occ2_act->id(), descr_f12f12_key, moints4_rtime,
                  i1i2i1i2_f12f12_ints);

    compute_VX(ij_ij, VX_output, f12f12_idx, i1i2i1i2_f12f12_ints, f12_idx,
               f12_idx, f12_ij_ints, f12_ij_ints, Xij_ij);

    //
    // Xij_ji = X^ij_ji = f^ij_ab f^ab_ji
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " X^ij_ji : " << endl;
    fill_n(Xij_ji, nocc12, 0.0);

    // X^ij_ji = (f12f12)^ij_ji
    Ref<DistArray4> i1i2i2i1_f12f12_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ1_act->id(), occ2_act->id(), occ2_act->id(),
                    occ1_act->id(), descr_f12f12_key, moints4_rtime,
                    i1i2i2i1_f12f12_ints);
    } else {
      i1i2i2i1_f12f12_ints = i1i2i1i2_f12f12_ints;
    }

    compute_VX(ij_ji, VX_output, f12f12_idx, i1i2i2i1_f12f12_ints, f12_idx,
               f12_idx, f12_ij_ints, f12_ji_ints, Xij_ji);

    if (spin1 != spin2 && num_unique_spincases2 == 3) {

      // Xji_ij = X^ji_ij = f^ji_ab f^ab_ij
      if (debug_ >= DefaultPrintThresholds::mostN2)
        ExEnv::out0() << endl << indent << spinletters << " X^ji_ij : " << endl;

      // X^ji_ij = (f12f12)^ji_ij
      Ref<DistArray4> i2i1i1i2_f12f12_ints;
      activate_ints(occ2_act->id(), occ1_act->id(), occ1_act->id(),
                    occ2_act->id(), descr_f12f12_key, moints4_rtime,
                    i2i1i1i2_f12f12_ints);

      compute_VX(ji_ij, VX_output, f12f12_idx, i2i1i1i2_f12f12_ints, f12_idx,
                 f12_idx, f12_ji_ints, f12_ij_ints, Xji_ij);
      i2i1i1i2_f12f12_ints->deactivate();

      // Xji_ji = X^ji_ji = f^ji_ab f^ab_ji
      if (debug_ >= DefaultPrintThresholds::mostN2)
        ExEnv::out0() << endl << indent << spinletters << " X^ji_ji : " << endl;

      // X^ji_ji = (f12f12)^ji_ji
      Ref<DistArray4> i2i1i2i1_f12f12_ints;
      activate_ints(occ2_act->id(), occ1_act->id(), occ2_act->id(),
                    occ1_act->id(), descr_f12f12_key, moints4_rtime,
                    i2i1i2i1_f12f12_ints);

      compute_VX(ji_ji, VX_output, f12f12_idx, i2i1i2i1_f12f12_ints, f12_idx,
                 f12_idx, f12_ji_ints, f12_ji_ints, Xji_ji);
      i2i1i2i1_f12f12_ints->deactivate();
    }

    if (debug_ >= DefaultPrintThresholds::N2 && spin1 == spin2)
      print_antisym_intermediate(spincase, "X^ij_ij", Xij_ij, Xij_ji, nocc1_act,
                                 nocc2_act);

    if (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2) {
      // Alpha-beta case
      print_intermediate(spincase, "X^ij_ij", Xij_ij, nocc1_act, nocc2_act);
      print_intermediate(spincase, "X^ij_ji", Xij_ji, nocc1_act, nocc2_act);
      if (spin1 != spin2 && num_unique_spincases2 == 3) {
        print_intermediate(spincase, "X^ji_ij", Xji_ij, nocc1_act, nocc2_act);
        print_intermediate(spincase, "X^ji_ji", Xji_ji, nocc1_act, nocc2_act);
      }
    }

    //
    // Compute B intermediate matrix
    //
    // B^ij_ij
    // alpha beta case: not antisymmetrized
    // alpha alpha, beta beta case: B^ij_ij - B^ij_ji (antisymmetrized)
    //

    // This is for alpha apha, beta beta and alpha beta cases
    //   alpha beta case is indicated by the following:
    //   B^i(alpha)j(beta)_^i(alpha)j(beta)
    // = R^i(alpha)j(beta)_a(alpha)b(beta) * f^a(alpha)_a(alpha) * R^a(alpha)b(beta)_i(alpha)j(beta)

    // B^ij_ij += (f12t1f12)^ij_ij
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " B(diag) contribution"
          << endl;
    compute_Y(ij_ij, 1.0, f12t1f12_idx, i1i2i1i2_f12f12_ints, Bij_ij);
    if (spin1 != spin2 && num_unique_spincases2 == 3)
      compute_Y(ji_ji, 1.0, f12t1f12_idx, i1i2i1i2_f12f12_ints, Bji_ji);

    // P part of B intermediate

    // Store all the ints, prefactor, index for b1b2k1k2 and output
    // for each term
    std::vector<Ref<DistArray4> > Pijij_f12_ints;
    std::vector<Ref<DistArray4> > Pijij_fx12_ints;
    std::vector<double> P_prefactors;
    std::vector<int> Pijij_idx;
    std::vector<std::string> Pijij_output;

    // P +=  f^ij_PQ K^P_R f^RQ_ij    (antisymmetrized)

    //    =   f^ij_PQ K^P_R f^RQ_ij - f^ij_PQ K^P_R f^RQ_ji  =>B^ij_ij
    //      - f^ji_PQ K^P_R f^RQ_ij + f^ji_PQ K^P_R f^RQ_ji  =>B^ji_ji
    //
    // when spin1 != spin2: P+= f^ij_PQ K^P_R f^RQ_ij + f^ji_PQ K^P_R f^RQ_ji
    //                        = f^ij_PQ f^P_K Q _ij + f^ji_PQ f^P_K Q _ji ;

    Ref<DistArray4> b_i1i2P1P2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), ribs1->id(), ribs2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2P1P2_ints);

    Ref<DistArray4> b_i1i2PK1P2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), Kribs1->id(), ribs2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2PK1P2_ints);

    Ref<DistArray4> b_i2i1P2P1_ints;
    Ref<DistArray4> b_i2i1PK2P1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), ribs2->id(), ribs1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1P2P1_ints);
      activate_ints(occ2_act->id(), occ1_act->id(), Kribs2->id(), ribs1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1PK2P1_ints);
    } else {
      //spin1 == spin2 or closed shell
      b_i2i1P2P1_ints = b_i1i2P1P2_ints;
      b_i2i1PK2P1_ints = b_i1i2PK1P2_ints;
    }

    Pijij_output.push_back("P(f^ij_PQ K^P_R f^RQ_ij)");
    P_prefactors.push_back(1.0);
    Pijij_f12_ints.push_back(b_i1i2P1P2_ints);
    Pijij_fx12_ints.push_back(b_i1i2PK1P2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1P2P1_ints);
    Pijij_fx12_ints.push_back(b_i2i1PK2P1_ints);
    Pijij_idx.push_back(ji_ji);

    // P += f^ij_Pm f^P_F m _ij + f^ji_Pm f^P_F m _ji
    // P += f^ij_Pm f^P_F m _ij
    Ref<DistArray4> b_i1i2P1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), ribs1->id(), occ2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2P1i2_ints);

    Ref<DistArray4> b_i1i2PF1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), Fribs1->id(), occ2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2PF1i2_ints);

    // P += f^ji_Pm F^P_Q f^Qm_ji
    Ref<DistArray4> b_i2i1P2i1_ints;
    Ref<DistArray4> b_i2i1PF2i1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), ribs2->id(), occ1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1P2i1_ints);
      activate_ints(occ2_act->id(), occ1_act->id(), Fribs2->id(), occ1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1PF2i1_ints);
    } else {
      //spin1 == spin2 or closed shell
      b_i2i1P2i1_ints = b_i1i2P1i2_ints;
      b_i2i1PF2i1_ints = b_i1i2PF1i2_ints;
    }

    Pijij_output.push_back("P(including f^ij_Pm F^P_Q f^Qm_ij)");
    P_prefactors.push_back(1.0);
    Pijij_f12_ints.push_back(b_i1i2P1i2_ints);
    Pijij_fx12_ints.push_back(b_i1i2PF1i2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1P2i1_ints);
    Pijij_fx12_ints.push_back(b_i2i1PF2i1_ints);
    Pijij_idx.push_back(ji_ji);

    // P += f^ij_pe F^p_q f^qe_ij + f^ji_pe F^p_q f^qe_ji
    //
    //    = f^ij_pe f^p_F e _ij  + ...

    // P += f^ij_pe F^p_q f^qe_ij
    Ref<DistArray4> b_i1i2p1e2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), orbs1->id(), vir2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2p1e2_ints);

    Ref<DistArray4> b_i1i2pF1e2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), forbsp1->id(), vir2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2pF1e2_ints);

    // P += f^ji_pe F^p_q f^qe_ji
    Ref<DistArray4> b_i2i1p2e1_ints;
    Ref<DistArray4> b_i2i1pF2e1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), orbs2->id(), vir1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1p2e1_ints);
      activate_ints(occ2_act->id(), occ1_act->id(), forbsp2->id(), vir1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1pF2e1_ints);
    } else {
      //spin1 == spin2 or closed shell
      b_i2i1p2e1_ints = b_i1i2p1e2_ints;
      b_i2i1pF2e1_ints = b_i1i2pF1e2_ints;
    }

    Pijij_output.push_back("P(including f^ij_pe F^p_q f^qe_ij)");
    P_prefactors.push_back(1.0);
    Pijij_f12_ints.push_back(b_i1i2p1e2_ints);
    Pijij_fx12_ints.push_back(b_i1i2pF1e2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1p2e1_ints);
    Pijij_fx12_ints.push_back(b_i2i1pF2e1_ints);
    Pijij_idx.push_back(ji_ji);

    // P += f^ij_mA F^m_P f^PA_ij + f^ji_mA F^m_P f^PA_ji + h.c.
    //    = 2.0 * (f^ij_mA f^m_F A _ij + ...)

    // P += 2 f^ij_mA F^m_P f^PA_ij
    Ref<DistArray4> b_i1i2m1A2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), occ1->id(), cabs2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2m1A2_ints);

    Ref<DistArray4> b_i1i2mF1A2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), Focc1->id(), cabs2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2mF1A2_ints);

    // P += 2 f^ji_mA F^m_P f^PA_ji
    Ref<DistArray4> b_i2i1m2A1_ints;
    Ref<DistArray4> b_i2i1mF2A1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), occ2->id(), cabs1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1m2A1_ints);
      activate_ints(occ2_act->id(), occ1_act->id(), Focc2->id(), cabs1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1mF2A1_ints);
      //ExEnv::out0() << indent << "activate: b_i2i1m1A2_ints, b_i2i1mF1A2_key" << endl;
    } else {
      b_i2i1m2A1_ints = b_i1i2m1A2_ints;
      b_i2i1mF2A1_ints = b_i1i2mF1A2_ints;
    }

    Pijij_output.push_back("P(including 2 f^ij_mA F^m_P f^PA_ij)");
    P_prefactors.push_back(2.0);
    Pijij_f12_ints.push_back(b_i1i2m1A2_ints);
    Pijij_fx12_ints.push_back(b_i1i2mF1A2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1m2A1_ints);
    Pijij_fx12_ints.push_back(b_i2i1mF2A1_ints);
    Pijij_idx.push_back(ji_ji);

    // P += f^ij_Ae F^A_p f^pe_ij + f^ji_Ae F^A_p f^pe_ji + h.c.
    //    = 2 * (f^ij_Ae f^A_F e _ij ...)

    // P += 2 f^ij_Ae F^A_p f^pe_ij
    Ref<DistArray4> b_i1i2A1e2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), orbs1->id(), vir2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2A1e2_ints);

    Ref<DistArray4> b_i1i2AF1e2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), forbsA1->id(), vir2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2AF1e2_ints);

    // P += 2 f^ji_Ae F^A_p f^pe_ji
    Ref<DistArray4> b_i2i1A2e1_ints;
    Ref<DistArray4> b_i2i1AF2e1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), orbs2->id(), vir1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1A2e1_ints);
      activate_ints(occ2_act->id(), occ1_act->id(), forbsA2->id(), vir1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1AF2e1_ints);
    } else {
      //spin1 == spin2 or closed shell
      b_i2i1A2e1_ints = b_i1i2A1e2_ints;
      b_i2i1AF2e1_ints = b_i1i2AF1e2_ints;
    }

    Pijij_output.push_back("P(including 2 f^ij_Ae F^A_p f^pe_ij)");
    P_prefactors.push_back(2.0);
    Pijij_f12_ints.push_back(b_i1i2A1e2_ints);
    Pijij_fx12_ints.push_back(b_i1i2AF1e2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1A2e1_ints);
    Pijij_fx12_ints.push_back(b_i2i1AF2e1_ints);
    Pijij_idx.push_back(ji_ji);

    // P -= f^ij_mA F^m_n f^nA_ij + f^ji_mA F^m_n f^nA_ji
    //    = f^ij_mA f^m_F A _ij ...

    // P -= f^ij_mA F^m_n f^nA_ij
    // f^ij_mA is already computed in the previous step
    Ref<DistArray4> b_i1i2mf1A2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), focc1->id(), cabs2->id(),
                  descr_f12_key, moints4_rtime, b_i1i2mf1A2_ints);

    // P -= f^ji_mA F^m_n f^nA_ji
    // f^ji_mA is already computed in the previous step
    Ref<DistArray4> b_i2i1mf2A1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), focc2->id(), cabs1->id(),
                    descr_f12_key, moints4_rtime, b_i2i1mf2A1_ints);
      //ExEnv::out0() << indent << "activate: b_i2i1mf1A2_ints" << endl;
    } else {
      //spin1 == spin2 or closed shell
      b_i2i1mf2A1_ints = b_i1i2mf1A2_ints;
    }

    Pijij_output.push_back("P(including f^ij_mA F^m_n f^nA_ij)");
    P_prefactors.push_back(-1.0);
    Pijij_f12_ints.push_back(b_i1i2m1A2_ints);
    Pijij_fx12_ints.push_back(b_i1i2mf1A2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1m2A1_ints);
    Pijij_fx12_ints.push_back(b_i2i1mf2A1_ints);
    Pijij_idx.push_back(ji_ji);

    accumulate_P_YxF(Pijij_output, Pijij_idx, P_prefactors, f12_idx, f12_idx,
                     Pijij_f12_ints, Pijij_fx12_ints, Pij_ij);

    //
    // Q += 1/2[  (f12f12)^ij_Pj (f+k)^P_i - (f12f12)^ji_Pj (f+k)^P_i
    //          + (f12f12)^ij_iP (f+k)^P_j - (f12f12)^ji_iP (f+k)^P_j
    //          + h.c.]
    //    =  (f12f12)^ij_Pj (f+k)^P_i - (f12f12)^ji_Pj (f+k)^P_i
    //     + (f12f12)^ji_Pi (f+k)^P_j - (f12f12)^ij_Pi (f+k)^P_j
    //
    // spin1 != spin2: Q += (f12f12)^ij_Pj (f+k)^P_i + (f12f12)^ji_Pi (f+k)^P_j ;
    //
    // spin1 == spin2:  Q +=   (f12f12)^ij_Pj (f+k)^P_i - (f12f12)^ji_Pj (f+k)^P_i
    //                      - (f12f12)^ij_Pi (f+k)^P_j + (f12f12)^ji_Pi (f+k)^P_j
    //
    //                    =   (f12f12)^ij _i_hJ j - (f12f12)^ji _i_hJ j
    //                      - (f12f12)^ij _j_hJ i + (f12f12)^ji _j_hJ i
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " Q^ij_ij : " << endl;

    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << indent << "Q += (f12f12)^ij_Pj (f+k)^P_i " << endl;
    Ref<DistArray4> q_i1i2hi1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), r12eval_->hj_i_P(spin1)->id(),
                  occ2_act->id(), descr_f12f12_key, moints4_rtime,
                  q_i1i2hi1i2_ints);
    compute_Y(ij_ij, 1.0, f12f12_idx, q_i1i2hi1i2_ints, Qij_ij);

    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << indent << "Q += (f12f12)^ji_Pi (f+k)^P_j " << endl;
    Ref<DistArray4> q_i2i1hi2i1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(),
                    r12eval_->hj_i_P(spin2)->id(), occ1_act->id(),
                    descr_f12f12_key, moints4_rtime, q_i2i1hi2i1_ints);
    } else {
      // spin1 == spin2 or closed shell
      q_i2i1hi2i1_ints = q_i1i2hi1i2_ints;
    }
    compute_Y(ji_ji, 1.0, f12f12_idx, q_i2i1hi2i1_ints, Qij_ij);

    // B^ij_ij = B^ij_ij - P + Q
    for (int i1 = 0; i1 < nocc1_act; ++i1) {
      for (int i2 = 0; i2 < nocc2_act; ++i2) {
        const int i1i2 = i1 * nocc2_act + i2;
        Bij_ij[i1i2] += -Pij_ij[i1i2] + Qij_ij[i1i2];
      }
    }
    if (debug_ >= DefaultPrintThresholds::mostN2
        || (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2))
      print_intermediate(spincase, "B^ij_ij", Bij_ij, nocc1_act, nocc2_act);

    //
    // B^ij_ji
    //

    // B^ij_ji += (f12t1f12)^ij_ji
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " B(diag) contribution"
          << endl;
    compute_Y(ij_ji, 1.0, f12t1f12_idx, i1i2i2i1_f12f12_ints, Bij_ji);
    if (spin1 != spin2 && num_unique_spincases2 == 3)
      compute_Y(ji_ij, 1.0, f12t1f12_idx, i1i2i2i1_f12f12_ints, Bji_ij);

    // P part of B intermediate

    // Store all the ints, prefactor, index for b1b2k1k2 and output
    // for each term
    // Pijji_f12_ints =  Pijij_f12_ints;
    // Pijji_prefactors = P_prefactors;
    std::vector<Ref<DistArray4> > Pijji_fx12_ints;
    std::vector<int> Pijji_idx;
    std::vector<std::string> Pijji_output;

    Ref<DistArray4> b_i2i1PK1P2_ints;
    Ref<DistArray4> b_i1i2PK2P1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), Kribs1->id(), ribs2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1PK1P2_ints);
      activate_ints(occ1_act->id(), occ2_act->id(), Kribs2->id(), ribs1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2PK2P1_ints);
    }
    else {
      b_i2i1PK1P2_ints = b_i1i2PK1P2_ints;
      b_i1i2PK2P1_ints = b_i1i2PK1P2_ints;
    }

    Pijji_output.push_back("P(f^ij_PQ K^Q_R f^PR_ji)");
    // see above: Pijij_f12_ints.push_back(b_i1i2P1P2_ints)
    Pijji_fx12_ints.push_back(b_i2i1PK1P2_ints);
    Pijji_idx.push_back(ij_ji);
    // see above: Pijij_f12_ints.push_back(b_i2i1P2P1_ints);
    Pijji_fx12_ints.push_back(b_i1i2PK2P1_ints);
    Pijji_idx.push_back(ji_ij);

    // P += f^ji_Pm F^P_Q f^Qm_ji
    Ref<DistArray4> b_i2i1PF1i2_ints;
    Ref<DistArray4> b_i1i2PF2i1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), Fribs1->id(), occ2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1PF1i2_ints);
      activate_ints(occ1_act->id(), occ2_act->id(), Fribs2->id(), occ1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2PF2i1_ints);
    }
    else {
      b_i2i1PF1i2_ints = b_i1i2PF1i2_ints;
      b_i1i2PF2i1_ints = b_i1i2PF1i2_ints;
    }

    Pijji_output.push_back("P(including f^ij_Pm F^P_Q f^Qm_ji)");
    // see above: Pijij_f12_ints.push_back(b_i1i2P1i2_ints)
    Pijji_fx12_ints.push_back(b_i2i1PF1i2_ints);
    Pijji_idx.push_back(ij_ji);
    // see above: Pijij_f12_ints.push_back(b_i2i1P2i1_ints);
    Pijji_fx12_ints.push_back(b_i1i2PF2i1_ints);
    Pijji_idx.push_back(ji_ij);

    Ref<DistArray4> b_i2i1pF1e2_ints;
    Ref<DistArray4> b_i1i2pF2e1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), forbsp1->id(), vir2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1pF1e2_ints);
      activate_ints(occ1_act->id(), occ2_act->id(), forbsp2->id(), vir1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2pF2e1_ints);
    }
    else {
      b_i2i1pF1e2_ints = b_i1i2pF1e2_ints;
      b_i1i2pF2e1_ints = b_i1i2pF1e2_ints;
    }

    Pijji_output.push_back("P(including f^ij_pe F^p_q f^qe_ji)");
    // see above: Pijij_f12_ints.push_back(b_i1i2p1e2_ints)
    Pijji_fx12_ints.push_back(b_i2i1pF1e2_ints);
    Pijji_idx.push_back(ij_ji);
    // see above: Pijij_f12_ints.push_back(b_i2i1p2e1_ints);
    Pijji_fx12_ints.push_back(b_i1i2pF2e1_ints);
    Pijji_idx.push_back(ji_ij);

    Ref<DistArray4> b_i2i1mF1A2_ints;
    Ref<DistArray4> b_i1i2mF2A1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), Focc1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1mF1A2_ints);
      activate_ints(occ1_act->id(), occ2_act->id(), Focc2->id(), cabs1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2mF2A1_ints);
    }
    else {
      b_i2i1mF1A2_ints = b_i1i2mF1A2_ints;
      b_i1i2mF2A1_ints = b_i1i2mF1A2_ints;
    }

    Pijji_output.push_back("P(including 2 f^ij_mA F^m_P f^PA_ji)");
    // see above: Pijij_f12_ints.push_back(b_i1i2m1A2_ints)
    Pijji_fx12_ints.push_back(b_i2i1mF1A2_ints);
    Pijji_idx.push_back(ij_ji);
    // see above: Pijij_f12_ints.push_back(b_i2i1m2A1_ints);
    Pijji_fx12_ints.push_back(b_i1i2mF2A1_ints);
    Pijji_idx.push_back(ji_ij);

    Ref<DistArray4> b_i2i1AF1e2_ints;
    Ref<DistArray4> b_i1i2AF2e1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(), forbsA1->id(), vir2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1AF1e2_ints);
      activate_ints(occ1_act->id(), occ2_act->id(), forbsA2->id(), vir1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2AF2e1_ints);
    }
    else {
      b_i2i1AF1e2_ints = b_i1i2AF1e2_ints;
      b_i1i2AF2e1_ints = b_i1i2AF1e2_ints;
    }

    Pijji_output.push_back("P(including 2 f^ij_Ae F^A_p f^pe_ji)");
    // see above: Pijij_f12_ints.push_back(b_i1i2A1e2_ints)
    Pijji_fx12_ints.push_back(b_i2i1AF1e2_ints);
    Pijji_idx.push_back(ij_ji);
    // see above: Pijij_f12_ints.push_back(b_i2i1A2e1_ints);
    Pijji_fx12_ints.push_back(b_i1i2AF2e1_ints);
    Pijji_idx.push_back(ji_ij);

    Ref<DistArray4> b_i1i2mf2A1_ints;
    Ref<DistArray4> b_i2i1mf1A2_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ1_act->id(), occ2_act->id(), focc2->id(), cabs1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2mf2A1_ints);
      activate_ints(occ2_act->id(), occ1_act->id(), focc1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1mf1A2_ints);
    }
    else {
      b_i1i2mf2A1_ints = b_i1i2mf1A2_ints;
      b_i2i1mf1A2_ints = b_i1i2mf1A2_ints;
    }

    Pijji_output.push_back("P(including f^ij_mA F^m_n f^nA_ji)");
    // see above: Pijij_f12_ints.push_back(b_i1i2m1A2_ints)
    Pijji_fx12_ints.push_back(b_i2i1mf1A2_ints);
    Pijji_idx.push_back(ij_ji);
    // see above: Pijij_f12_ints.push_back(b_i2i1m2A1_ints)
    Pijji_fx12_ints.push_back(b_i1i2mf2A1_ints);
    Pijji_idx.push_back(ji_ij);

    accumulate_P_YxF(Pijji_output, Pijji_idx, P_prefactors, f12_idx, f12_idx,
                     Pijij_f12_ints, Pijji_fx12_ints, Pij_ji);

    // Qij_ji += (f12f12)^ij_Pi (f+k)^P_j + (f12f12)^ij_jP (f+k)^P_i
    //         = (f12f12)^ij _j_hJ i + (f12f12)^ji _i_hJ j
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " Q^ij_ji : " << endl;

    Ref<DistArray4> q_i1i2hi2i1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ1_act->id(), occ2_act->id(),
                    r12eval_->hj_i_P(spin2)->id(), occ1_act->id(),
                    descr_f12f12_key, moints4_rtime, q_i1i2hi2i1_ints);
    } else {
      q_i1i2hi2i1_ints = q_i2i1hi2i1_ints;
    }
    compute_Y(ij_ji, 1.0, f12f12_idx, q_i1i2hi2i1_ints, Qij_ji);

    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl;

    Ref<DistArray4> q_i2i1hi1i2_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(),
                    r12eval_->hj_i_P(spin1)->id(), occ2_act->id(),
                    descr_f12f12_key, moints4_rtime, q_i2i1hi1i2_ints);
    } else {
      q_i2i1hi1i2_ints = q_i1i2hi1i2_ints;
    }
    compute_Y(ji_ij, 1.0, f12f12_idx, q_i2i1hi1i2_ints, Qij_ji);

    // B^ij_ji = B^ij_ji - P^ij_ji + Q^ij_ji
    for (int i1 = 0; i1 < nocc1_act; ++i1) {
      for (int i2 = 0; i2 < nocc2_act; ++i2) {
        const int i1i2 = i1 * nocc2_act + i2;
        Bij_ji[i1i2] += -Pij_ji[i1i2] + Qij_ji[i1i2];
      }
    }
    if (debug_ >= DefaultPrintThresholds::mostN2
        || (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2))
      print_intermediate(spincase, "B^ij_ji", Bij_ji, nocc1_act, nocc2_act);

    //
    //  B^i(beta)j(alpha)_^i(beta)j(alpha)
    //
    //   R^i(alpha)j(beta)_a(alpha)b(beta) * f^b(beta)_b(beta) * R^a(alpha)b(beta)_i(alpha)j(beta)
    // =>  electron 1 <-> electron 2  i<->j  a<->b
    // = R^i(beta)j(alpha)_a(beta)r(alpha) * f^a(beta)_a(beta) * R^a(beta)b(alpha)_i(beta)j(alpha)
    //
    if (num_unique_spincases2 == 3 && spin1 != spin2) {

      // Store all the ints for b1b2k1k2 and output for each term
      // prefactor and  index are the same as Pij_ij
      std::vector<Ref<DistArray4> > Pjiji_f12_ints;
      std::vector<Ref<DistArray4> > Pjiji_fx12_ints;
      std::vector<std::string> Pjiji_output;

      // P +=  f^ij_PQ K^Q_R f^PR_ij    (antisymmetrized)

      //    =   f^ij_QP K^Q_R f^RP_ij - f^ij_QP K^Q_R f^RP_ji
      //      - f^ji_QP K^Q_R f^RP_ij + f^ji_QP K^Q_R f^RP_ji
      //
      // when spin1 != spin2: P+= f^ij_QP K^Q_R f^RP_ij + f^ji_QP K^Q_R f^RP_ji
      //                        = f^ij_QP f^Q_K P _ij + f^ji_QP f^Q_K P _ji ;

      // P += f^ji_QP f^Q_K P _ji
      Ref<DistArray4> b_i2i1P1P2_ints;
      activate_ints(occ2_act->id(), occ1_act->id(), ribs1->id(), ribs2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1P1P2_ints);
      Ref<DistArray4> b_i1i2P2P1_ints;
      activate_ints(occ1_act->id(), occ2_act->id(), ribs2->id(), ribs1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2P2P1_ints);

      Pjiji_output.push_back("P(f^ij_PQ K^Q_R f^PR_ij) (beta alpha case)");
      Pjiji_f12_ints.push_back(b_i2i1P1P2_ints);
      Pjiji_fx12_ints.push_back(b_i2i1PK1P2_ints);
      Pjiji_f12_ints.push_back(b_i1i2P2P1_ints);
      Pjiji_fx12_ints.push_back(b_i1i2PK2P1_ints);

      // P += f^ij_Pm f^P_F m _ij + f^ji_Pm f^P_F m _ji
      Ref<DistArray4> b_i2i1P1i2_ints;
      activate_ints(occ2_act->id(), occ1_act->id(), ribs1->id(), occ2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1P1i2_ints);
      Ref<DistArray4> b_i1i2P2i1_ints;
      activate_ints(occ1_act->id(), occ2_act->id(), ribs2->id(), occ1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2P2i1_ints);

      Pjiji_output.push_back(
          "P(including f^ij_Pm F^P_Q f^Qm_ij) (beta alpha case)");
      Pjiji_f12_ints.push_back(b_i2i1P1i2_ints);
      Pjiji_fx12_ints.push_back(b_i2i1PF1i2_ints);
      Pjiji_f12_ints.push_back(b_i1i2P2i1_ints);
      Pjiji_fx12_ints.push_back(b_i1i2PF2i1_ints);

      // P += f^ij_pe F^p_q f^qe_ij + f^ji_pe F^p_q f^qe_ji
      Ref<DistArray4> b_i2i1p1e2_ints;
      activate_ints(occ2_act->id(), occ1_act->id(), orbs1->id(), vir2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1p1e2_ints);
      Ref<DistArray4> b_i1i2p2e1_ints;
      activate_ints(occ1_act->id(), occ2_act->id(), orbs2->id(), vir1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2p2e1_ints);

      Pjiji_output.push_back(
          "P(including f^ij_pe F^p_q f^qe_ij) (beta alpha case)");
      Pjiji_f12_ints.push_back(b_i2i1p1e2_ints);
      Pjiji_fx12_ints.push_back(b_i2i1pF1e2_ints);
      Pjiji_f12_ints.push_back(b_i1i2p2e1_ints);
      Pjiji_fx12_ints.push_back(b_i1i2pF2e1_ints);

      // P += f^ij_mA F^m_P f^PA_ij + f^ji_mA F^m_P f^PA_ji + h.c.
      Ref<DistArray4> b_i2i1m1A2_ints;
      activate_ints(occ2_act->id(), occ1_act->id(), occ1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1m1A2_ints);
      Ref<DistArray4> b_i1i2m2A1_ints;
      activate_ints(occ1_act->id(), occ2_act->id(), occ2->id(), cabs1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2m2A1_ints);

      Pjiji_output.push_back(
          "P(including 2 f^ij_mA F^m_P f^PA_ij) (beta alpha case)");
      Pjiji_f12_ints.push_back(b_i2i1m1A2_ints);
      Pjiji_fx12_ints.push_back(b_i2i1mF1A2_ints);
      Pjiji_f12_ints.push_back(b_i1i2m2A1_ints);
      Pjiji_fx12_ints.push_back(b_i1i2mF2A1_ints);

      // P += f^ij_Ae F^A_p f^pe_ij + f^ji_Ae F^A_p f^pe_ji + h.c.
      Ref<DistArray4> b_i2i1A1e2_ints;
      activate_ints(occ2_act->id(), occ1_act->id(), orbs1->id(), vir2->id(),
                    descr_f12_key, moints4_rtime, b_i2i1A1e2_ints);
      Ref<DistArray4> b_i1i2A2e1_ints;
      activate_ints(occ1_act->id(), occ2_act->id(), orbs2->id(), vir1->id(),
                    descr_f12_key, moints4_rtime, b_i1i2A2e1_ints);

      Pjiji_output.push_back(
          "P(including 2 f^ij_Ae F^A_p f^pe_ij) (beta alpha case)");
      Pjiji_f12_ints.push_back(b_i2i1A1e2_ints);
      Pjiji_fx12_ints.push_back(b_i2i1AF1e2_ints);
      Pjiji_f12_ints.push_back(b_i1i2A2e1_ints);
      Pjiji_fx12_ints.push_back(b_i1i2AF2e1_ints);

      // P -= f^ij_mA F^m_n f^nA_ij + f^ji_mA F^m_n f^nA_ji
      Pjiji_output.push_back(
          "P(including f^ij_mA F^m_n f^nA_ij) (beta alpha case)");
      Pjiji_f12_ints.push_back(b_i2i1m1A2_ints);
      Pjiji_fx12_ints.push_back(b_i2i1mf1A2_ints);
      Pjiji_f12_ints.push_back(b_i1i2m2A1_ints);
      Pjiji_fx12_ints.push_back(b_i1i2mf2A1_ints);

      accumulate_P_YxF(Pjiji_output, Pijij_idx, P_prefactors, f12_idx,
                       f12_idx, Pjiji_f12_ints, Pjiji_fx12_ints,
                       Pji_ji);

      // B^ij_ij = B^ij_ij - P + Q
      for (int i2 = 0; i2 < nocc2_act; ++i2) {
        for (int i1 = 0; i1 < nocc1_act; ++i1) {
          const int i2i1 = i2 * nocc1_act + i1;
          const int i1i2 = i1 * nocc2_act + i2;
          Bji_ji[i2i1] += -Pji_ji[i2i1] + Qij_ij[i1i2];
        }
      }
      if (debug_ >= DefaultPrintThresholds::mostN2
          || (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2))
        print_intermediate(spincase, "Bji_ji:", Bji_ji,
                           nocc2_act, nocc1_act);

      Pji_ji = NULL;

      //
      // B^ij_ji
      // B^i(beta)j(alpha)_j(alpha)i(beta)

      // P part of B intermediate

      // Store all the ints, prefactor, index for b1b2k1k2 and output
      // for each term
      // Pjiij_f12_ints =  Pjiji_f12_ints;
      // Pjiij_prefactprs = P_prefactors;
      // Pjiij_idx = Pijji_idx
      std::vector<Ref<DistArray4> > Pjiij_fx12_ints;
      std::vector<std::string> Pjiij_output;

      Pjiij_output.push_back("P(f^ij_PQ K^Q_R f^PR_ji) (beta alpha case)");
      // see above: Pjiji_f12_ints.push_back(b_i2i1P1P2_ints);
      Pjiij_fx12_ints.push_back(b_i1i2PK1P2_ints);
      // Pijji_idx.push_back(ij_ji);
      // see above: Pjiji_f12_ints.push_back(b_i1i2P2P1_ints);
      Pjiij_fx12_ints.push_back(b_i2i1PK2P1_ints);
      // Pijji_idx.push_back(ji_ij);

      Pjiij_output.push_back(
          "P(including f^ij_Pm F^P_Q f^Qm_ji) (beta alpha case)");
      // see above Pjiji_f12_ints.push_back(b_i2i1P1i2_ints);
      Pjiij_fx12_ints.push_back(b_i1i2PF1i2_ints);
      // Pijji_idx.push_back(ij_ji);
      // see above Pjiji_f12_ints.push_back(b_i1i2P2i1_ints);
      Pjiij_fx12_ints.push_back(b_i2i1PF2i1_ints);
      // Pijji_idx.push_back(ji_ij);

      Pjiij_output.push_back(
          "P(including f^ij_pe F^p_q f^qe_ji) (beta alpha case)");
      // see above: Pjiji_f12_ints.push_back(b_i2i1p1e2_ints);
      Pjiij_fx12_ints.push_back(b_i1i2pF1e2_ints);
      // Pijji_idx.push_back(ij_ji);
      // see above: Pjiji_f12_ints.push_back(b_i1i2p2e1_ints);
      Pjiij_fx12_ints.push_back(b_i2i1pF2e1_ints);
      // Pijji_idx.push_back(ji_ij);

      Pjiij_output.push_back(
          "P(including 2 f^ij_mA F^m_P f^PA_ji) (beta alpha case)");
      // see above: Pjiji_f12_ints.push_back(b_i2i1m1A2_ints);
      Pjiij_fx12_ints.push_back(b_i1i2mF1A2_ints);
      // Pijji_idx.push_back(ij_ji);
      // Pjiji_f12_ints.push_back(b_i1i2m2A1_ints);
      Pjiij_fx12_ints.push_back(b_i2i1mF2A1_ints);
      // Pijji_idx.push_back(ji_ij);

      Pjiij_output.push_back(
          "P(including 2 f^ij_Ae F^A_p f^pe_ji) (beta alpha case)");
      // see above: Pjiji_f12_ints.push_back(b_i2i1A1e2_ints);
      Pjiij_fx12_ints.push_back(b_i1i2AF1e2_ints);
      // Pijji_idx.push_back(ij_ji);
      // see above: Pjiji_f12_ints.push_back(b_i1i2A2e1_ints);
      Pjiij_fx12_ints.push_back(b_i2i1AF2e1_ints);
      // Pijji_idx.push_back(ji_ij);

      Pjiij_output.push_back(
          "P(including f^ij_mA F^m_n f^nA_ji) (beta alpha case)");
      // see above: Pjiji_f12_ints.push_back(b_i2i1m1A2_ints);
      Pjiij_fx12_ints.push_back(b_i1i2mf1A2_ints);
      // Pijji_idx.push_back(ij_ji);
      // see above: Pjiji_f12_ints.push_back(b_i1i2m2A1_ints);
      Pjiij_fx12_ints.push_back(b_i2i1mf2A1_ints);
      // Pijji_idx.push_back(ji_ij);

      accumulate_P_YxF(Pjiij_output, Pijji_idx, P_prefactors, f12_idx,
                       f12_idx, Pjiji_f12_ints, Pjiij_fx12_ints,
                       Pji_ij);

      // B^ij_ji = B^ij_ji - P^ij_ji + Q^ij_ji
      for (int i2 = 0; i2 < nocc2_act; ++i2) {
        for (int i1 = 0; i1 < nocc1_act; ++i1) {
          const int i2i1 = i2 * nocc1_act + i1;
          const int i1i2 = i1 * nocc2_act + i2;
          Bji_ij[i2i1] += -Pji_ij[i2i1] + Qij_ji[i1i2];
        }
      }
      if (debug_ >= DefaultPrintThresholds::mostN2
          || (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2))
        print_intermediate(spincase, "Bji_ij:", Bji_ij,
                           nocc2_act, nocc1_act);

      Pji_ij = NULL;

      // Deactivate all the ints for Pji_ij
      for (std::vector<Ref<DistArray4> >::iterator it =
          Pjiji_f12_ints.begin(); it < Pjiji_f12_ints.end(); ++it) {
        (*it)->deactivate();
      }
      for (std::vector<Ref<DistArray4> >::iterator it =
          Pjiji_fx12_ints.begin(); it < Pjiji_fx12_ints.end(); ++it) {
        (*it)->deactivate();
      }


    } // end of alpha-beta open-shell case

    // code for testing 1e density matrix
    // compute X terms using non-canonical orbitals
#if 1
    double* Xioccj_ij = NULL;
    double* Xioccj_ji = NULL;
    double* Xijocc_ij = NULL;
    double* Xijocc_ji = NULL;
    double* Xjiocc_ij = NULL;
    double* Xjiocc_ji = NULL;
    double* Xjocci_ij = NULL;
    double* Xjocci_ji = NULL;
    if (this->r12eval()->compute_1rdm() ||
        r12intermediates_->Onerdm_cc_computed()) {

      const Ref<OrbitalSpace>& focc1_im = r12eval()->F_i_m(spin1);
      const Ref<OrbitalSpace>& focc2_im = r12eval()->F_i_m(spin2);
      Xioccj_ij = new double[nocc12];
      Xioccj_ji = new double[nocc12];
      fill_n(Xioccj_ij, nocc12, 0.0);
      fill_n(Xioccj_ji, nocc12, 0.0);

      // Xioccj_ij += (f12f12)^ioccj_ij
      Ref<DistArray4> ioccjij_f12f12_ints;
      activate_ints(focc1_im->id(), occ2_act->id(), occ1_act->id(),
                    occ2_act->id(), descr_f12f12_key, moints4_rtime,
                    ioccjij_f12f12_ints);
      // store all the ints
      std::vector<Ref<DistArray4> > f12_ioccj_ints;

      // Xioccj_ij -= f^ioccj_pq f^pq_ij
      Ref<DistArray4> ioccjpp_ints;
      activate_ints(focc1_im->id(), occ2_act->id(), orbs1->id(), orbs2->id(),
                    descr_f12_key, moints4_rtime, ioccjpp_ints);
      f12_ioccj_ints.push_back(ioccjpp_ints);

      // Xioccj_ij-= f^ioccj_ma' f^ma'_ij
      Ref<DistArray4> ioccjia_ints;
      activate_ints(focc1_im->id(), occ2_act->id(), occ1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime, ioccjia_ints);
      f12_ioccj_ints.push_back(ioccjia_ints);

      // Xioccj_ij-= f^ioccj_a'm f^a'm_ij
      Ref<DistArray4> ioccjai_ints;
      activate_ints(focc1_im->id(), occ2_act->id(), cabs1->id(), occ2->id(),
                    descr_f12_key, moints4_rtime, ioccjai_ints);
      f12_ioccj_ints.push_back(ioccjai_ints);

      compute_VX(ij_ij, VX_output, f12f12_idx, ioccjij_f12f12_ints, f12_idx,
                 f12_idx, f12_ioccj_ints, f12_ij_ints, Xioccj_ij);
//    print_intermediate(spincase, "X^ioccj_ij", Xioccj_ij, nocc1_act, nocc2_act);

      // Xioccj_ji
      Ref<DistArray4> ioccjji_f12f12_ints;
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        activate_ints(focc1_im->id(), occ2_act->id(),occ2_act->id(),
                      occ1_act->id(), descr_f12f12_key, moints4_rtime,
                      ioccjji_f12f12_ints);
      } else {
          ioccjji_f12f12_ints = ioccjij_f12f12_ints;
      }

      std::vector<Ref<DistArray4> > f12_ij_ints;
      Ref<DistArray4> i1i2p1p2_ints;
      activate_ints(occ1_act->id(), occ2_act->id(), orbs1->id(), orbs2->id(),
                    descr_f12_key, moints4_rtime, i1i2p1p2_ints);
      f12_ij_ints.push_back(i1i2p1p2_ints);

      Ref<DistArray4> i1i2i1a2_ints;
      activate_ints(occ1_act->id(), occ2_act->id(), occ1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime, i1i2i1a2_ints);
      f12_ij_ints.push_back(i1i2i1a2_ints);

      Ref<DistArray4> i1i2a1i2_ints;
      activate_ints(occ1_act->id(), occ2_act->id(), cabs1->id(), occ2->id(),
                    descr_f12_key, moints4_rtime, i1i2a1i2_ints);
      f12_ij_ints.push_back(i1i2a1i2_ints);

      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        activate_ints(occ2_act->id(), occ1_act->id(), orbs1->id(), orbs2->id(),
                      descr_f12_key, moints4_rtime, i2i1p1p2_ints);

        activate_ints(occ2_act->id(), occ1_act->id(), occ1->id(), cabs2->id(),
                      descr_f12_key, moints4_rtime, i2i1i1a2_ints);

        activate_ints(occ2_act->id(), occ1_act->id(), cabs1->id(), occ2->id(),
                      descr_f12_key, moints4_rtime, i2i1a1i2_ints);
      } else {
        i2i1p1p2_ints = i1i2p1p2_ints;
        i2i1i1a2_ints = i1i2i1a2_ints;
        i2i1a1i2_ints = i1i2a1i2_ints;
      }
      std::vector<Ref<DistArray4> > f12_ji_ints;
      f12_ji_ints.push_back(i2i1p1p2_ints);
      f12_ji_ints.push_back(i2i1i1a2_ints);
      f12_ji_ints.push_back(i2i1a1i2_ints);

      compute_VX(ij_ji, VX_output, f12f12_idx, ioccjji_f12f12_ints, f12_idx,
                 f12_idx, f12_ioccj_ints, f12_ji_ints, Xioccj_ji);
//    print_intermediate(spincase, "X^ioccj_ji", Xioccj_ji, nocc1_act, nocc2_act);

      ioccjij_f12f12_ints->deactivate();
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        ioccjji_f12f12_ints->deactivate();
      }
      for (std::vector<Ref<DistArray4> >::iterator it = f12_ioccj_ints.begin();
          it < f12_ioccj_ints.end(); ++it) {
        (*it)->deactivate();
      }

      Xijocc_ij = new double[nocc12];
      Xijocc_ji = new double[nocc12];
      fill_n(Xijocc_ij, nocc12, 0.0);
      fill_n(Xijocc_ji, nocc12, 0.0);

      // Xijocc_ij += (f12f12)^ijocc_ij
      Ref<DistArray4> ijoccij_f12f12_ints;
      activate_ints(occ1_act->id(), focc2_im->id(), occ1_act->id(),
                    occ2_act->id(), descr_f12f12_key, moints4_rtime,
                    ijoccij_f12f12_ints);
      // store all the ints
      std::vector<Ref<DistArray4> > f12_ijocc_ints;

      // Xijocc_ij -= f^ijocc_pq f^pq_ij
      Ref<DistArray4> ijoccpp_ints;
      activate_ints(occ1_act->id(), focc2_im->id(), orbs1->id(), orbs2->id(),
                    descr_f12_key, moints4_rtime, ijoccpp_ints);
      f12_ijocc_ints.push_back(ijoccpp_ints);

      // Xijocc_ij -= f^ijocc_ma' f^ma'_ij
      Ref<DistArray4> ijoccia_ints;
      activate_ints(occ1_act->id(), focc2_im->id(), occ1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime, ijoccia_ints);
      f12_ijocc_ints.push_back(ijoccia_ints);

      // Xijocc_ij -= f^ijocc_a'm f^a'm_ij
      Ref<DistArray4> ijoccai_ints;
      activate_ints(occ1_act->id(), focc2_im->id(), cabs1->id(), occ2->id(),
                    descr_f12_key, moints4_rtime, ijoccai_ints);
      f12_ijocc_ints.push_back(ijoccai_ints);

      compute_VX(ij_ij, VX_output, f12f12_idx, ijoccij_f12f12_ints, f12_idx,
                 f12_idx, f12_ijocc_ints, f12_ij_ints, Xijocc_ij);
//    print_intermediate(spincase, "X^ijocc_ij", Xijocc_ij, nocc1_act, nocc2_act);

      // Xijocc_ji
      Ref<DistArray4> ijoccji_f12f12_ints;
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        activate_ints(occ1_act->id(), focc2_im->id(), occ2_act->id(),
                      occ1_act->id(), descr_f12f12_key, moints4_rtime,
                      ijoccji_f12f12_ints);
      } else {
          ijoccji_f12f12_ints = ijoccij_f12f12_ints;
      }

      //ExEnv::out0() << endl << "f12_ji_ints size" <<
      compute_VX(ij_ji, VX_output, f12f12_idx, ijoccji_f12f12_ints, f12_idx,
                 f12_idx, f12_ijocc_ints, f12_ji_ints, Xijocc_ji);
//    print_intermediate(spincase, "X^ijocc_ji", Xijocc_ji, nocc1_act, nocc2_act);

      ijoccij_f12f12_ints->deactivate();
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        ijoccji_f12f12_ints->deactivate();
      }
      for (std::vector<Ref<DistArray4> >::iterator it = f12_ijocc_ints.begin();
          it < f12_ijocc_ints.end(); ++it) {
        (*it)->deactivate();
      }

      // these are only needed in alpha-beta open-shell
      if (spin1 != spin2 && num_unique_spincases2 == 3) {

        Xjiocc_ij = new double[nocc12];
        fill_n(Xjiocc_ij, nocc12, 0.0);
        Xjiocc_ji = new double[nocc12];
        fill_n(Xjiocc_ji, nocc12, 0.0);

        // X^jiocc_ij = (f12f12)^jiocc_ij
        Ref<DistArray4> j2iocc1ci1j2_f12f12_ints;
        activate_ints(occ2_act->id(), focc1_im->id(), occ1_act->id(),
                      occ2_act->id(), descr_f12f12_key, moints4_rtime,
                      j2iocc1ci1j2_f12f12_ints);

        // store all the ints
        std::vector<Ref<DistArray4> > f12_jiocc_ints;

        // Xjiocc_ij -= f^jiocc_pq f^pq_ij
        Ref<DistArray4> jioccpp_ints;
        activate_ints(occ2_act->id(), focc1_im->id(), orbs1->id(), orbs2->id(),
                      descr_f12_key, moints4_rtime, jioccpp_ints);
        f12_jiocc_ints.push_back(jioccpp_ints);

        // Xjiocc_ij-= f^jiocc_ma' f^ma'_ij
        Ref<DistArray4> jioccia_ints;
        activate_ints(occ2_act->id(), focc1_im->id(), occ1->id(), cabs2->id(),
                      descr_f12_key, moints4_rtime, jioccia_ints);
        f12_jiocc_ints.push_back(jioccia_ints);

        // Xjiocc_ij-= f^jiocc_a'm f^a'm_ij
        Ref<DistArray4> jioccai_ints;
        activate_ints(occ2_act->id(), focc1_im->id(), cabs1->id(), occ2->id(),
                      descr_f12_key, moints4_rtime, jioccai_ints);
        f12_jiocc_ints.push_back(jioccai_ints);

        compute_VX(ji_ij, VX_output, f12f12_idx, j2iocc1ci1j2_f12f12_ints, f12_idx,
                   f12_idx, f12_jiocc_ints, f12_ij_ints, Xjiocc_ij);
//      print_intermediate(spincase, "X^jiocc_ij", Xjiocc_ij, nocc1_act, nocc2_act);

        // X^jiocc_ji = (f12f12)^jiocc_ji
        Ref<DistArray4> j2iocc1cj2i1_f12f12_ints;
        activate_ints(occ2_act->id(), focc1_im->id(), occ2_act->id(),
                      occ1_act->id(), descr_f12f12_key, moints4_rtime,
                      j2iocc1cj2i1_f12f12_ints);

        compute_VX(ji_ji, VX_output, f12f12_idx, j2iocc1cj2i1_f12f12_ints, f12_idx,
                   f12_idx, f12_jiocc_ints, f12_ji_ints, Xjiocc_ji);
//      print_intermediate(spincase, "X^jiocc_ji", Xjiocc_ji, nocc1_act, nocc2_act);

        j2iocc1ci1j2_f12f12_ints->deactivate();
        j2iocc1cj2i1_f12f12_ints->deactivate();
        for (std::vector<Ref<DistArray4> >::iterator it = f12_jiocc_ints.begin();
            it < f12_jiocc_ints.end(); ++it) {
          (*it)->deactivate();
        }

        Xjocci_ij = new double[nocc12];
        fill_n(Xjocci_ij, nocc12, 0.0);
        Xjocci_ji = new double[nocc12];
        fill_n(Xjocci_ji, nocc12, 0.0);

        // Xjocci_ij += (f12f12)^ijocc_ij
        Ref<DistArray4> jocc2i1i1j2_f12f12_ints;
        activate_ints(focc2_im->id(), occ1_act->id(), occ1_act->id(),
                      occ2_act->id(), descr_f12f12_key, moints4_rtime,
                      jocc2i1i1j2_f12f12_ints);
        // store all the ints
        std::vector<Ref<DistArray4> > f12_jocci_ints;

        // Xjocci_ij -= f^jocci_pq f^pq_ij
        Ref<DistArray4> joccipp_ints;
        activate_ints(focc2_im->id(), occ1_act->id(), orbs1->id(), orbs2->id(),
                      descr_f12_key, moints4_rtime, joccipp_ints);
        f12_jocci_ints.push_back(joccipp_ints);

        // Xjocci_ij-= f^jocci_ma' f^ma'_ij
        Ref<DistArray4> jocciia_ints;
        activate_ints(focc2_im->id(), occ1_act->id(), occ1->id(), cabs2->id(),
                      descr_f12_key, moints4_rtime, jocciia_ints);
        f12_jocci_ints.push_back(jocciia_ints);

        // Xjocci_ij-= f^jocci_a'm f^a'm_ij
        Ref<DistArray4> jocciai_ints;
        activate_ints(focc2_im->id(), occ1_act->id(), cabs1->id(), occ2->id(),
                      descr_f12_key, moints4_rtime, jocciai_ints);
        f12_jocci_ints.push_back(jocciai_ints);

        compute_VX(ji_ij, VX_output, f12f12_idx, jocc2i1i1j2_f12f12_ints, f12_idx,
                   f12_idx, f12_jocci_ints, f12_ij_ints, Xjocci_ij);

        // Xjocci_ji += (f12f12)^jocci_ji
        Ref<DistArray4> jocc2i1i2j1_f12f12_ints;
        activate_ints(focc2_im->id(), occ1_act->id(), occ2_act->id(),
                      occ1_act->id(), descr_f12f12_key, moints4_rtime,
                      jocc2i1i2j1_f12f12_ints);

        compute_VX(ji_ji, VX_output, f12f12_idx, jocc2i1i2j1_f12f12_ints, f12_idx,
                   f12_idx, f12_jocci_ints, f12_ji_ints, Xjocci_ji);

        jocc2i1i1j2_f12f12_ints->deactivate();
        jocc2i1i2j1_f12f12_ints->deactivate();
        for (std::vector<Ref<DistArray4> >::iterator it = f12_jocci_ints.begin();
            it < f12_jocci_ints.end(); ++it) {
          (*it)->deactivate();
        }

//      print_intermediate(spincase, "X^jocci_ij", Xjocci_ij, nocc1_act, nocc2_act);
//      print_intermediate(spincase, "X^jocci_ji", Xjocci_ji, nocc1_act, nocc2_act);
      } // end of if (spin1 != spin2 && num_unique_spincases2 == 3)

      for (std::vector<Ref<DistArray4> >::iterator it = f12_ij_ints.begin();
        it < f12_ij_ints.end(); ++it) {
        (*it)->deactivate();
      }
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        for (std::vector<Ref<DistArray4> >::iterator it = f12_ji_ints.begin();
          it < f12_ji_ints.end(); ++it) {
          (*it)->deactivate();
        }
      }

      // print out the V coupling, X, and B contributions
      double E_V = 0.0;
      double E_Vcoupling = 0.0;
      double E_Bcoupling = 0.0;
      double E_X = 0.0;
      double E_B = 0.0;

      // X terms using non-canonical orbitals
      double E_X_noca = 0.0;

      if (spin1 == spin2) {
        // Alpha_alpha or beta_beta case
        double E_Vij = 0.0;
        double E_Vji = 0.0;
        double E_Bij = 0.0;
        double E_Bji = 0.0;
        for (int i1 = 0; i1 < nocc1_act; ++i1) {
          for (int i2 = i1 + 1; i2 < nocc2_act; ++i2) {
            const int ij = i1 * nocc2_act + i2;

            E_V += 2.0 * C_1 * (Vij_ij[ij] - Vij_ji[ij]);
            E_X += - C_1 * C_1* (evals_i1(i1) + evals_i2(i2))
                               * (Xij_ij[ij] - Xij_ji[ij]);
            E_B += C_1 * C_1 * (Bij_ij[ij] - Bij_ji[ij]);

            E_X_noca += - C_1 * C_1* (Xioccj_ij[ij] - Xioccj_ji[ij]
                                    + Xijocc_ij[ij] - Xijocc_ji[ij]);

            if (this->r12eval()->coupling() == true
                || this->r12eval()->ebc() == false) {
              E_Vcoupling += 2.0 * C_1
                  * (Vij_ij_coupling[ij] - Vij_ji_coupling[ij]);

              if (do_mp2 && this->r12eval()->coupling())
                E_Bcoupling += C_1 * C_1 * (Bij_ij_coupling[ij] - Bij_ji_coupling[ij]);
            }

          }
        }
        if (num_unique_spincases2 == 2) {
          E_V_tot += E_V * 2.0;
          E_X_tot += E_X * 2.0;
          E_B_tot += E_B * 2.0;
          E_X_noca_tot += E_X_noca * 2.0;
          E_Vcoupling_tot += E_Vcoupling * 2.0;
          E_Bcoupling_tot += E_Bcoupling * 2.0;
        } else {
            E_V_tot += E_V;
            E_X_tot += E_X;
            E_B_tot += E_B;
            E_X_noca_tot = +E_X_noca;
            E_Vcoupling_tot += E_Vcoupling;
            E_Bcoupling_tot += E_Bcoupling;
        }
      } else if (num_unique_spincases2 == 2) {
        // Alpha_beta case for closed shell
        double E_Vij = 0.0;
        double E_Vji = 0.0;
        double E_Bij = 0.0;
        double E_Bji = 0.0;
        for (int i1 = 0; i1 < nocc1_act; ++i1) {
          for (int i2 = 0; i2 < nocc2_act; ++i2) {
            int ij = i1 * nocc2_act + i2;

            E_V += 2.0  * (0.5 * (C_0 + C_1) * Vij_ij[ij]
                         + 0.5 * (C_0 - C_1) * Vij_ji[ij]);
            E_X += - pow(0.5 * (C_0 + C_1), 2)
                            * (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij]
                   - 0.25 * (C_0 * C_0 - C_1 * C_1)
                            * (evals_i1(i1) + evals_i2(i2)) * (2.0 * Xij_ji[ij])
                   - pow(0.5 * (C_0 - C_1), 2)
                            * (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij];
            E_B += pow(0.5 * (C_0 + C_1), 2) * Bij_ij[ij]
                   + 0.25 * (C_0 * C_0 - C_1 * C_1) * 2.0 * Bij_ji[ij]
                   + pow(0.5 * (C_0 - C_1), 2) * Bij_ij[ij] ;

            E_X_noca += - pow(0.5 * (C_0 + C_1), 2)
                                  * (Xioccj_ij[ij] + Xijocc_ij[ij])
                        - 0.25 * (C_0 * C_0 - C_1 * C_1)
                               * (Xioccj_ji[ij] + Xijocc_ji[ij]) * 2.0
                        - pow(0.5 * (C_0 - C_1), 2)
                                  * (Xioccj_ij[ij] + Xijocc_ij[ij]);

            E_Vij += Vij_ij[ij];
            E_Vji += Vij_ji[ij];
            E_Bij += Bij_ij[ij] - (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij];
            E_Bji += Bij_ji[ij] - (evals_i1(i1) + evals_i2(i2)) * Xij_ji[ij];

            if (this->r12eval()->coupling() == true
                || this->r12eval()->ebc() == false) {
              E_Vcoupling += 2.0
                  * (0.5 * (C_0 + C_1) * Vij_ij_coupling[ij]
                      + 0.5 * (C_0 - C_1) * Vij_ji_coupling[ij]);
              if (do_mp2 && this->r12eval()->coupling())
                E_Bcoupling += (pow(0.5 * (C_0 + C_1), 2)
                    * Bij_ij_coupling[ij]
                    + 0.25 * (C_0 * C_0 - C_1 * C_1)
                        * (2.0 * Bij_ji_coupling[ij])
                    + pow(0.5 * (C_0 - C_1), 2)
                        * (Bij_ij_coupling[ij]));
            }

          }
        }

        ExEnv::out0() << endl << "Individual contributions to the R12 energy (AlphaBeta):"
                      << endl << "E_Vij = " << scprintf("%12.10f", E_Vij)
                      << "     E_Vji = " << scprintf("%12.10f", E_Vji)
                      << endl << "E_Bij = " << scprintf("%12.10f", E_Bij)
                      << "     E_Bji = " << scprintf("%12.10f", E_Bji)
                      << endl << "E_VBX = " << scprintf("%12.10f", (E_Vij*5/4-E_Vji/4+E_Bij*7/32+E_Bji/32))
                      << endl << "E_VBX(new fixed T) = " << scprintf("%12.10f", (E_Vij*2-E_Vji+E_Bij/2-E_Bji/4))
                      << endl << endl;
        E_V_tot += E_V;
        E_X_tot += E_X;
        E_B_tot += E_B;
        E_X_noca_tot = +E_X_noca;
        E_Vcoupling_tot += E_Vcoupling;
        E_Bcoupling_tot += E_Bcoupling;
      } else {
        // Alpha_beta case for open shell
        for (int i1 = 0; i1 < nocc1_act; ++i1) {
          for (int i2 = 0; i2 < nocc2_act; ++i2) {
            double Hij_pair_energy;
            int ij = i1 * nocc2_act + i2;
            int ji = i2 * nocc1_act + i1;

            E_V += (2.0
                * (0.5 * (C_0 + C_1) * Vij_ij[ij]
                    + 0.5 * (C_0 - C_1) * Vij_ji[ij]));
            E_X += pow(0.5 * (C_0 + C_1), 2)
                * (- (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij])
                + 0.25 * (C_0 * C_0 - C_1 * C_1)
                    * (- (evals_i1(i1) + evals_i2(i2))
                            * (Xij_ji[ij] + Xji_ij[ij]))
                + pow(0.5 * (C_0 - C_1), 2)
                    * (- (evals_i1(i1) + evals_i2(i2)) * Xji_ji[ij]);
            E_B += pow(0.5 * (C_0 + C_1), 2) * Bij_ij[ij]
                   + 0.25 * (C_0 * C_0 - C_1 * C_1)
                         * (Bij_ji[ij] + Bji_ij[ji])
                   + pow(0.5 * (C_0 - C_1), 2) * Bji_ji[ji];

            E_X_noca += - pow(0.5 * (C_0 + C_1), 2)
                                  * (Xioccj_ij[ij] + Xijocc_ij[ij])
                        - 0.25 * (C_0 * C_0 - C_1 * C_1)
                               *  (Xioccj_ji[ij] + Xijocc_ji[ij] + Xjiocc_ij[ij] + Xjocci_ij[ij])
                        - pow(0.5 * (C_0 - C_1), 2)
                                  * (Xjiocc_ji[ij]+ Xjocci_ji[ij]);

            if (this->r12eval()->coupling() == true
                || this->r12eval()->ebc() == false) {
              E_Vcoupling += 2.0
                              * (0.5 * (C_0 + C_1) * Vij_ij_coupling[ij]
                               + 0.5 * (C_0 - C_1) * Vij_ji_coupling[ij]);
              if (do_mp2 && this->r12eval()->coupling())
                E_Bcoupling += (pow(0.5 * (C_0 + C_1), 2)
                                  * (Bij_ij_coupling[ij])
                                  + 0.25 * (C_0 * C_0 - C_1 * C_1)
                                      * ((Bij_ji_coupling[ij] + Bji_ij_coupling[ji])
                                        )
                                  + pow(0.5 * (C_0 - C_1), 2)
                                      * (Bji_ji_coupling[ji]));
            }
          }
        }
        E_V_tot += E_V;
        E_X_tot += E_X;
        E_B_tot += E_B;
        E_X_noca_tot = +E_X_noca;
        E_Vcoupling_tot += E_Vcoupling;
        E_Bcoupling_tot += E_Bcoupling;
      }
      delete[] Xioccj_ij;
      delete[] Xioccj_ji;
      delete[] Xijocc_ij;
      delete[] Xijocc_ji;
      if (spin1 != spin2 && num_unique_spincases2 == 3) {
        delete[] Xjiocc_ij;
        delete[] Xjiocc_ji;
        delete[] Xjocci_ij;
        delete[] Xjocci_ji;
      }

    }
#endif

    // Deactivate all the ints
    i1i2i1i2_ints->deactivate();
    for (std::vector<Ref<DistArray4> >::iterator it = f12_ij_ints.begin();
        it < f12_ij_ints.end(); ++it) {
      (*it)->deactivate();
    }
    i1i2i1i2_f12f12_ints->deactivate();

    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      i1i2i2i1_ints->deactivate();
      for (std::vector<Ref<DistArray4> >::iterator it = f12_ji_ints.begin();
          it < f12_ji_ints.end(); ++it) {
        (*it)->deactivate();
      }
      i1i2i2i1_f12f12_ints->deactivate();
    }
    // Deactivate ints for P
    for (std::vector<Ref<DistArray4> >::iterator it = Pijij_f12_ints.begin();
        it < Pijij_f12_ints.end(); ++it) {
      (*it)->deactivate();
      ++it;
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        (*it)->deactivate();
      }
    }
    for (std::vector<Ref<DistArray4> >::iterator it = Pijij_fx12_ints.begin();
        it < Pijij_fx12_ints.end(); ++it) {
      (*it)->deactivate();
      ++it;
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        (*it)->deactivate();
      }
    }
    // Deactivate ints for Q
    q_i1i2hi1i2_ints->deactivate();
    q_i2i1hi2i1_ints->deactivate();
    if (num_unique_spincases2 == 3 && spin1 != spin2) {
      q_i1i2hi2i1_ints->deactivate();
      q_i2i1hi1i2_ints->deactivate();
    }
    //
    // when spin1 != spin2: B^ji_ij = B^ij_ji, B^ji_ji = B^ij_ij
    //

    // Print the antisymmetrized Qij_ij and Bij_ij
    // for alpha_alpha and beta_beta case
    if (spin1 == spin2) {
      if (debug_ >= DefaultPrintThresholds::mostN2)
        print_antisym_intermediate(spincase, "Q^ij_ij:", Qij_ij, Qij_ji,
                                   nocc1_act, nocc2_act);
      if (debug_ >= DefaultPrintThresholds::N2)
        print_antisym_intermediate(spincase, "B^ij_ij:", Bij_ij, Bij_ji,
                                   nocc1_act, nocc2_act);
    }

    // Compute the f12 correction pair energy
    if (spin1 == spin2) {
      // Alpha_alpha or beta_beta case

      for (int i1 = 0; i1 < nocc1_act; ++i1) {
        for (int i2 = i1 + 1; i2 < nocc2_act; ++i2) {

          double Hij_pair_energy;
          const int ij = i1 * nocc2_act + i2;

          if (debug_ >= DefaultPrintThresholds::mostN2) {
            ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << ","
                << i2 << " e(VT) = " << (2.0 * C_1 * (Vij_ij[ij] - Vij_ji[ij]))
                << endl;
            ExEnv::out0()
                << indent
                << "Hij_pair_energy: ij = "
                << i1
                << ","
                << i2
                << " e(TBT) = "
                << (C_1 * C_1
                    * (Bij_ij[ij] - Bij_ji[ij]
                        - (evals_i1(i1) + evals_i2(i2))
                            * (Xij_ij[ij] - Xij_ji[ij]))) << endl;
          }

          Hij_pair_energy = 2.0 * C_1 * (Vij_ij[ij] - Vij_ji[ij])
              + C_1 * C_1
                  * (Bij_ij[ij] - Bij_ji[ij]
                      - (evals_i1(i1) + evals_i2(i2))
                          * (Xij_ij[ij] - Xij_ji[ij]));

          if (this->r12eval()->coupling() == true
              || this->r12eval()->ebc() == false) {

            if (debug_ >= DefaultPrintThresholds::mostN2) {
              ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << ","
                  << i2 << " e(couplingV) = "
                  << (2.0 * C_1 * (Vij_ij_coupling[ij] - Vij_ji_coupling[ij]))
                  << endl;
            }

            Hij_pair_energy += 2.0 * C_1
                * (Vij_ij_coupling[ij] - Vij_ji_coupling[ij]);

            if (do_mp2 && this->r12eval()->coupling()) {

              if (debug_ >= DefaultPrintThresholds::mostN2) {
                ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << ","
                    << i2 << " e(couplingB) = "
                    << (C_1 * C_1 * (Bij_ij_coupling[ij] - Bij_ji_coupling[ij]))
                    << endl;
              }

              Hij_pair_energy += C_1 * C_1 * (Bij_ij_coupling[ij] - Bij_ji_coupling[ij]);
            }
          }

          // Get indices for the lower triangle matrix
          // index = nrow*(nrow-1) + ncolumn
          // nrow > ncolumn (diagonal elements are excluded)
          const int i21 = i2 * (i2 - 1) / 2 + i1;
          ef12_[spin].set_element(i21, Hij_pair_energy);
        }
      }
    } else if (num_unique_spincases2 == 2) {
      // Alpha_beta case for closed shell
      for (int i1 = 0; i1 < nocc1_act; ++i1) {
        for (int i2 = 0; i2 < nocc2_act; ++i2) {
          double Hij_pair_energy;
          int ij = i1 * nocc2_act + i2;

          if (debug_ >= DefaultPrintThresholds::mostN2) {
            ExEnv::out0()
                << indent
                << "Hij_pair_energy: ij = "
                << i1
                << ","
                << i2
                << " e(VT) = "
                << (2.0
                    * (0.5 * (C_0 + C_1) * Vij_ij[ij]
                        + 0.5 * (C_0 - C_1) * Vij_ji[ij])) << endl;
            ExEnv::out0()
                << indent
                << "Hij_pair_energy: ij = "
                << i1
                << ","
                << i2
                << " e(TBT) = "
                << (pow(0.5 * (C_0 + C_1), 2)
                    * (Bij_ij[ij] - (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij])
                    + 0.25 * (C_0 * C_0 - C_1 * C_1)
                        * (2.0 * Bij_ji[ij]
                            - (evals_i1(i1) + evals_i2(i2)) * (2.0 * Xij_ji[ij]))
                    + pow(0.5 * (C_0 - C_1), 2)
                        * (Bij_ij[ij]
                            - (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij]))
                << endl;
          }

          Hij_pair_energy =
              2.0
                  * (0.5 * (C_0 + C_1) * Vij_ij[ij]
                      + 0.5 * (C_0 - C_1) * Vij_ji[ij])
                  + pow(0.5 * (C_0 + C_1), 2)
                      * (Bij_ij[ij] - (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij])
                  + 0.25 * (C_0 * C_0 - C_1 * C_1)
                      * (2.0 * Bij_ji[ij]
                          - (evals_i1(i1) + evals_i2(i2)) * (2.0 * Xij_ji[ij]))
                  + pow(0.5 * (C_0 - C_1), 2)
                      * (Bij_ij[ij] - (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij]);

          if (this->r12eval()->coupling() == true
              || this->r12eval()->ebc() == false) {
            if (debug_ >= DefaultPrintThresholds::mostN2) {
              ExEnv::out0()
                  << indent
                  << "Hij_pair_energy: ij = "
                  << i1
                  << ","
                  << i2
                  << " e(couplingV) = "
                  << (2.0
                      * (0.5 * (C_0 + C_1) * Vij_ij_coupling[ij]
                          + 0.5 * (C_0 - C_1) * Vij_ji_coupling[ij])) << endl;
            }
            Hij_pair_energy += 2.0
                * (0.5 * (C_0 + C_1) * Vij_ij_coupling[ij]
                    + 0.5 * (C_0 - C_1) * Vij_ji_coupling[ij]);

            if (do_mp2 && this->r12eval()->coupling()) {

              if (debug_ >= DefaultPrintThresholds::mostN2) {
                ExEnv::out0()
                    << indent
                    << "Hij_pair_energy: ij = "
                    << i1
                    << ","
                    << i2
                    << " e(couplingB) = "
                    << (pow(0.5 * (C_0 + C_1), 2)
                        * Bij_ij_coupling[ij]
                        + 0.25 * (C_0 * C_0 - C_1 * C_1)
                            * (2.0 * Bij_ji_coupling[ij])
                        + pow(0.5 * (C_0 - C_1), 2)
                            * (Bij_ij_coupling[ij])) << endl;
              }
              Hij_pair_energy += (pow(0.5 * (C_0 + C_1), 2)
                  * Bij_ij_coupling[ij]
                  + 0.25 * (C_0 * C_0 - C_1 * C_1)
                      * (2.0 * Bij_ji_coupling[ij])
                  + pow(0.5 * (C_0 - C_1), 2)
                      * (Bij_ij_coupling[ij]));
            }
          }

          const int i12 = i1 * nocc2_act + i2;
          ef12_[spin].set_element(i12, Hij_pair_energy);
        }
      }
    } else {
      // Alpha_beta case for open shell
      for (int i1 = 0; i1 < nocc1_act; ++i1) {
        for (int i2 = 0; i2 < nocc2_act; ++i2) {
          double Hij_pair_energy;
          int ij = i1 * nocc2_act + i2;
          int ji = i2 * nocc1_act + i1;

          double Hij_pair_energy_VT = (2.0
              * (0.5 * (C_0 + C_1) * Vij_ij[ij]
                  + 0.5 * (C_0 - C_1) * Vij_ji[ij]));
          double Hij_pair_energy_TBT = (pow(0.5 * (C_0 + C_1), 2)
              * (Bij_ij[ij]
                  - (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij])
              + 0.25 * (C_0 * C_0 - C_1 * C_1)
                  * ((Bij_ji[ij] + Bji_ij[ji])
                      - (evals_i1(i1) + evals_i2(i2))
                          * (Xij_ji[ij] + Xji_ij[ij]))
              + pow(0.5 * (C_0 - C_1), 2)
                  * (Bji_ji[ji]
                      - (evals_i1(i1) + evals_i2(i2)) * Xji_ji[ij]));
          if (debug_ >= DefaultPrintThresholds::mostN2) {
            ExEnv::out0()
                << indent
                << "Hij_pair_energy: ij = "
                << i1
                << ","
                << i2
                << " e(VT) = "
                << scprintf("%20.15lf",Hij_pair_energy_VT) << endl;
            ExEnv::out0()
                << indent
                << "Hij_pair_energy: ij = "
                << i1
                << ","
                << i2
                << " e(TBT) = "
                << scprintf("%20.15lf",Hij_pair_energy_TBT)
                << endl;
          }

          if (this->r12eval()->coupling() == true
              || this->r12eval()->ebc() == false) {
            if (debug_ >= DefaultPrintThresholds::mostN2) {
              ExEnv::out0()
                  << indent
                  << "Hij_pair_energy: ij = "
                  << i1
                  << ","
                  << i2
                  << " e(couplingV) = "
                  << (2.0
                      * (0.5 * (C_0 + C_1) * Vij_ij_coupling[ij]
                          + 0.5 * (C_0 - C_1) * Vij_ji_coupling[ij])) << endl;
            }
            Hij_pair_energy_VT += 2.0
                * (0.5 * (C_0 + C_1) * Vij_ij_coupling[ij]
                    + 0.5 * (C_0 - C_1) * Vij_ji_coupling[ij]);

            if (do_mp2 && this->r12eval()->coupling()) {
              if (debug_ >= DefaultPrintThresholds::mostN2) {
                ExEnv::out0()
                << indent
                << "Hij_pair_energy: ij = "
                << i1
                << ","
                << i2
                << " e(couplingB) = "
                << (pow(0.5 * (C_0 + C_1), 2)
                    * (Bij_ij_coupling[ij])
                    + 0.25 * (C_0 * C_0 - C_1 * C_1)
                        * ((Bij_ji_coupling[ij] + Bji_ij_coupling[ji])
                          )
                    + pow(0.5 * (C_0 - C_1), 2)
                        * (Bji_ji_coupling[ji])) << endl;
              }
              Hij_pair_energy_VT += (pow(0.5 * (C_0 + C_1), 2)
                  * (Bij_ij_coupling[ij])
                  + 0.25 * (C_0 * C_0 - C_1 * C_1)
                      * ((Bij_ji_coupling[ij] + Bji_ij_coupling[ji])
                        )
                  + pow(0.5 * (C_0 - C_1), 2)
                      * (Bji_ji_coupling[ji]));
            }
          }

          Hij_pair_energy = Hij_pair_energy_VT + Hij_pair_energy_TBT;

          const int i12 = i1 * nocc2_act + i2;
          ef12_[spin].set_element(i12, Hij_pair_energy);
          ef12_VT[spin][i12] = Hij_pair_energy_VT;
          ef12_TBT[spin][i12] = Hij_pair_energy_TBT;
        }
      }
      delete[] Bji_ij;
      delete[] Bji_ji;
      delete[] Xji_ij;
      delete[] Xji_ji;
      delete[] Pji_ij;
      delete[] Pji_ji;
    }

    // deallocate memory
    delete[] Vij_ij;
    delete[] Vij_ji;
    delete[] Vij_ij_coupling;
    delete[] Vij_ji_coupling;
    delete[] Tij_ab;
    delete[] Xij_ij;
    delete[] Xij_ji;
    delete[] Bij_ij;
    delete[] Bij_ji;
    delete[] Pij_ij;
    delete[] Pij_ji;
    delete[] Qij_ij;
    delete[] Qij_ji;
    
  } // end of spincase loop

  //
  // Coupled cluster V contribution
  //
  if (r12intermediates_->T1_cc_computed()
      && r12intermediates_->T2_cc_computed()) {

    if (debug_ >= DefaultPrintThresholds::N2)
      ExEnv::out0() << endl << indent << "Coupled-cluster V contribution"
          << endl;
    // Alpha-alpha, beta-beta, or close-shell case:

    // uncomment, if need to compute ALL terms linear in amplitudes + V . T2^2
//#define CCSD_2_R12_COMPLETE_TO_FIRST_ORDER
#ifdef CCSD_2_R12_COMPLETE_TO_FIRST_ORDER
      assert(num_unique_spincases2 == 2); // closed-shell only for now

      const double U2T2_contrib = U2T2_ta();
      ExEnv::out0() << "U2T2_contrib = " << scprintf("%20.15lf\n", U2T2_contrib);
      const double U1T1_plus_C1T1_contrib = U1T1_plus_C1T1_ta();
      ExEnv::out0() << "U1T1_plus_C1T1_contrib = " << scprintf("%20.15lf\n", U1T1_plus_C1T1_contrib);
      const double VTT_contrib = VTT_ta();
      ExEnv::out0() << "VTT_contrib = " << scprintf("%20.15lf\n", VTT_contrib);
      const double PT2R12_contrib = PT2R12_ta();
      ExEnv::out0() << "PT2R12_contrib = " << scprintf("%20.15lf\n", PT2R12_contrib);

      ef12_[0](0) = ef12_[0](0)
          + 2 * (U2T2_contrib + U1T1_plus_C1T1_contrib) // 2 = hermitian terms in CC lagrangian a la Hylleraas functional
          + VTT_contrib // only comes from Lambda [[H,T],R], not from R [[H,T],T]
      ;
#endif

    //  V^ij_ij(cc) = 1/2*C1 * V^ab_ij*T^ij_ab + C1 * V^ia_ij*T^j_a
    //                                         + C1 * V^aj_ij*T^i_a
    // Alpha-beta case:
    //  V^ij_ij(cc) = 1/2*[(C0+C1)/2 * V^ab_ij*T^ij_ab + (C0-C1)/2 * V^ab_ji*T^ij_ab]
    //              + (C0+C1)/2 * V^ia_ij*T^j_a+(C0-C1)/2 * V^ia_ji*T^j_a
    //              + (C0+C1)/2 * V^aj_ij*T^i_a+(C0-C1)/2 * V^aj_ji*T^i_a (=>V^ja_ij*T^i_a)

    double* const VT2ij_ij = new double[nijij];
    double* const VT2ij_ji = new double[nijij];
    double* const VT1ij_ij = new double[nijij];
    double* const VT1ij_ji = new double[nijij];
    for (int s = 0; s < num_unique_spincases2; s++) {
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      const SpinCase1 spin1 = case1(spincase2);
      const SpinCase1 spin2 = case2(spincase2);

      if (r12eval()->dim_oo(spincase2).n() == 0)
        continue;
      const Ref<OrbitalSpace>& v1 = r12eval()->vir_act(spin1);
      const Ref<OrbitalSpace>& v2 = r12eval()->vir_act(spin2);
      const int nv1 = v1->rank();
      const int nv2 = v2->rank();
      const Ref<OrbitalSpace>& o1 = r12eval()->occ_act(spin1);
      const Ref<OrbitalSpace>& o2 = r12eval()->occ_act(spin2);
      const int no1 = o1->rank();
      const int no2 = o2->rank();
      const Ref<OrbitalSpace>& p1 = r12eval()->orbs(spin1);
      const Ref<OrbitalSpace>& p2 = r12eval()->orbs(spin2);
      const unsigned int np1 = p1->rank();
      const unsigned int np2 = p2->rank();
      const blasint nv12 = nv1 * nv2;
      const blasint one = 1;

      // Vpq_ij
      std::vector<Ref<DistArray4> > Vpq_vec = r12eval()->V_distarray4(spincase2,
                                                                      p1, p2);
      MPQC_ASSERT(Vpq_vec.size() == 1);
      Ref<DistArray4> Vpq_ij = Vpq_vec[0];
      // V^qp_ij (V^p2p1_i1i2) for open-shell alpha-beta case
      Ref<DistArray4> Vqp_ij = NULL;

      // VT2ij_ij = V^ab_ij * T^ij_ab
      fill_n(VT2ij_ij, nijij, 0.0);
      // extract Vab_ij from Vpq_ij
      Ref<DistArray4> Vab_ij;
      map(Vpq_ij, o1, o2, p1, p2, Vab_ij, o1, o2, v1, v2);
      //_print(spincase2, Vab_ij, prepend_spincase(spincase2,"Vab_ij matrix").c_str());

      Vab_ij->activate();
      T2[s]->activate();
      compute_YxF(ij_ij, 1.0, 0, 0, Vab_ij, T2[s], VT2ij_ij);
      if (debug_ >= DefaultPrintThresholds::N2)
        print_intermediate(spincase2, "V^ab_ij * T^ij_ab", VT2ij_ij, no1, no2);

      // VT2ij_ji = V^ab_ji * T^ij_ab
      if (spincase2 == AlphaBeta) {
        fill_n(VT2ij_ji, nijij, 0.0);

        if (num_unique_spincases2 == 3) {
          // V^qp_ij (V^p2p1_i1i2)
          std::vector<Ref<DistArray4> > Vqp_vec = r12eval()->V_distarray4(
              spincase2, p2, p1);
          MPQC_ASSERT(Vqp_vec.size() == 1);
          Vqp_ij = Vqp_vec[0];

          // extract Vba_ij from Vqp_ij
          Ref<DistArray4> Vba_ij;
          map(Vqp_ij, o1, o2, p2, p1, Vba_ij, o1, o2, v2, v1);

          Vba_ij->activate();
          double* const Vab_blk = new double[nv12];
          fill(Vab_blk, Vab_blk + nv12, 0.0);
          for (int i = 0; i < no1; ++i) {
            for (int j = 0; j < no2; ++j) {
              const double* Vba_blk = Vba_ij->retrieve_pair_block(i, j, 0);
              swap_e12(Vba_blk, nv2, nv1, Vab_blk);

              const double* T2ab_blk = T2[s]->retrieve_pair_block(i, j, 0);
              const double vt2ijji = F77_DDOT(&nv12, Vab_blk, &one, T2ab_blk,
                                              &one);
              const int ij = i * no2 + j;
              VT2ij_ji[ij] = vt2ijji;

              Vba_ij->release_pair_block(i, j, 0);
              T2[s]->release_pair_block(i, j, 0);
            }
          }
          delete[] Vab_blk;
          Vba_ij->deactivate();
          //print_intermediate(spincase2,VT2ij_ji,no1,no2);
        } else {
          compute_YxF(ji_ij, 1.0, 0, 0, Vab_ij, T2[s], VT2ij_ji);
        }
        if (debug_ >= DefaultPrintThresholds::N2)
          print_intermediate(spincase2, "V^ab_ji * T^ij_ab", VT2ij_ji, no1,
                             no2);
      }
      Vab_ij->deactivate();
      T2[s]->deactivate();

      //
      // VT1 contraction
      const bool swap_e12_V = true;
      const bool nonswap_e12_V = false;
      const bool T1_ja = 0;
      const bool T1_ia = 1;

      // VT1ij_ij = Via_ij * T1^j_a
      fill_n(VT1ij_ij, nijij, 0.0);

      // extract Via_ij from Vpq_ij
      Ref<DistArray4> Via_ij;
      map(Vpq_ij, o1, o2, p1, p2, Via_ij, o1, o2, o1, v2);
      //_print(spincase2, Via_ij, prepend_spincase(spincase2,"Via_ij matrix").c_str());

      // Convert T1^j_a (spin2) to an array
      const int no2v2 = no2 * nv2;
      double* const raw_T1_ja = new double[no2v2];
      T1[spin2].convert(raw_T1_ja);

      Via_ij->activate();
      contract_VT1(Via_ij, ij_ij, nonswap_e12_V, raw_T1_ja, nv2, T1_ja,
                   VT1ij_ij);
      if (debug_ >= DefaultPrintThresholds::N2)
        print_intermediate(spincase2, "VT1ij_ij = Via_ij * T1^j_a", VT1ij_ij,
                           no1, no2);

      // VT1ij_ji = V^ia_ji * T^ij_ab
      if (spincase2 == AlphaBeta) {
        fill_n(VT1ij_ji, nijij, 0.0);

        if (num_unique_spincases2 == 3) {
          // extract Vai_ij from Vqp_ij
          Ref<DistArray4> Vai_ij;
          map(Vqp_ij, o1, o2, p2, p1, Vai_ij, o1, o2, v2, o1);

          Vai_ij->activate();
          contract_VT1(Vai_ij, ij_ij, swap_e12_V, raw_T1_ja, nv2, T1_ja,
                       VT1ij_ji);
          Vai_ij->deactivate();
        } else {
          contract_VT1(Via_ij, ji_ij, nonswap_e12_V, raw_T1_ja, nv2, T1_ja,
                       VT1ij_ji);
        }
        if (debug_ >= DefaultPrintThresholds::N2)
          print_intermediate(spincase2, "VT1ij_ji = V^ia_ji * T1^j_a", VT1ij_ji,
                             no1, no2);
      }
      Via_ij->deactivate();
      delete[] raw_T1_ja;

      // VT1ij_ij += V^aj_ij * T1^i_a
      //
      // extract Vaj_ij from Vpq_ij
      Ref<DistArray4> Vaj_ij;
      map(Vpq_ij, o1, o2, p1, p2, Vaj_ij, o1, o2, v1, o2);
      //_print(spincase2, Vaj_ij, prepend_spincase(spincase2,"Vaj_ij matrix").c_str());

      // Convert T1^i_a (spin1) to an array
      const int no1v1 = no1 * nv1;
      double* const raw_T1_ia = new double[no1v1];
      T1[spin1].convert(raw_T1_ia);

      Vaj_ij->activate();
      contract_VT1(Vaj_ij, ij_ij, swap_e12_V, raw_T1_ia, nv1, T1_ia, VT1ij_ij);
      Vaj_ij->deactivate();
      if (debug_ >= DefaultPrintThresholds::N2)
        print_intermediate(spincase2, "VT1ij_ij += V^aj_ij * T1^i_a", VT1ij_ij,
                           no1, no2);

      // VT1ij_ji += V^aj_ji * T1^i_a = V^ja_ij * T1^i_a
      if (spincase2 == AlphaBeta) {

        Ref<DistArray4> Vja_ij;
        if (num_unique_spincases2 == 3) {
          // extract Vja_ij from Vqp_ij
          map(Vqp_ij, o1, o2, p2, p1, Vja_ij, o1, o2, o2, v1);
        } else {
          Vja_ij = Via_ij;
        }

        Vja_ij->activate();
        contract_VT1(Vja_ij, ij_ij, nonswap_e12_V, raw_T1_ia, nv1, T1_ia,
                     VT1ij_ji);
        Vja_ij->deactivate();
        if (debug_ >= DefaultPrintThresholds::N2)
          print_intermediate(spincase2, "VT1ij_ji += V^ja_ij * T1^i_a",
                             VT1ij_ji, no1, no2);
      }
      delete[] raw_T1_ia;

      // Add VT2 & VT1 contribution to ef12
      if (spin1 == spin2) {
        for (int i1 = 0; i1 < no1; ++i1) {
          for (int i2 = i1 + 1; i2 < no2; ++i2) {
            const int ij = i1 * no2 + i2;
            const int i21 = i2 * (i2 - 1) / 2 + i1;
            // 2.0 from Hylleraas functional
            const double Hij_pair_energy = 2.0 * C_1
                * (0.5 * VT2ij_ij[ij] + VT1ij_ij[ij]);
//          const double Hij_pair_energy = VT2ij_ij[ij];
            if (debug_ >= DefaultPrintThresholds::mostN2) {
              ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(CC) = "
                  << Hij_pair_energy << endl;
            }
            ef12_[s].accumulate_element(i21, Hij_pair_energy);
          }
        }
        ef12_[s].print("f12 pair energies");
      } else {
        // Alpha_beta case
        for (int i1 = 0; i1 < no1; ++i1) {
          for (int i2 = 0; i2 < no2; ++i2) {
            const int ij = i1 * no2 + i2;
            const double Hij_pair_energy = 2.0
                * (0.5 * (C_0 + C_1) * VT2ij_ij[ij]
                    + 0.5 * (C_0 - C_1) * VT2ij_ji[ij]
                    + 0.5 * (C_0 + C_1) * VT1ij_ij[ij]
                    + 0.5 * (C_0 - C_1) * VT1ij_ji[ij]);
            if (debug_ >= DefaultPrintThresholds::mostN2) {
              ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(CC) = "
                  << Hij_pair_energy << endl;
            }
            ef12_[s].accumulate_element(ij, Hij_pair_energy);
            ef12_VT[s][ij] += Hij_pair_energy;
          }
        }
        ef12_[s].print("f12 pair energies");
        for(int i12=0; i12<ef12_VT[s].size(); ++i12) {
          ExEnv::out0() << indent << scprintf("%20.15lf",ef12_VT[s][i12]) << endl;
          ExEnv::out0() << indent << scprintf("%20.15lf",ef12_TBT[s][i12]) << endl;
        }
      }

      if (this->r12eval()->compute_1rdm() == true ||
          r12intermediates_->Onerdm_cc_computed()) {
        double E_VT = 0.0;

        if (spin1 == spin2) {
          for (int i1 = 0; i1 < no1; ++i1) {
            for (int i2 = i1 + 1; i2 < no2; ++i2) {
              const int ij = i1 * no2 + i2;

              // 2.0 from Hylleraas functional
              const double Hij_pair_energy = 2.0 * C_1
                    * (0.5 * VT2ij_ij[ij] + VT1ij_ij[ij]);
              E_VT += Hij_pair_energy;
            }
          }
          if (num_unique_spincases2 == 2) {
            E_VT_tot += E_VT * 2.0 ;
          } else {
              E_VT_tot += E_VT;
          }
        } else {
          // Alpha_beta case
          for (int i1 = 0; i1 < no1; ++i1) {
            for (int i2 = 0; i2 < no2; ++i2) {
              const int ij = i1 * no2 + i2;
              const double Hij_pair_energy = 2.0
                       * (0.5 * (C_0 + C_1) * VT2ij_ij[ij]
                        + 0.5 * (C_0 - C_1) * VT2ij_ji[ij]
                        + 0.5 * (C_0 + C_1) * VT1ij_ij[ij]
                        + 0.5 * (C_0 - C_1) * VT1ij_ji[ij]);
              E_VT += Hij_pair_energy;
            }
          }
          E_VT_tot += E_VT;
        }
    }

    } // end of spin iteration
    delete[] VT2ij_ij;
    delete[] VT2ij_ji;
    delete[] VT1ij_ij;
    delete[] VT1ij_ji;
  } // end of CC V contribution

  if (this->r12eval()->compute_1rdm() == true ||
      r12intermediates_->Onerdm_cc_computed()) {
    ExEnv::out0() << endl << "Individual contributions to the R12 energy:"
                  << endl << "E_V = " << scprintf("%12.10f", E_V_tot)
                  << endl << "E_Vcoupling = " << scprintf("%12.10f", E_Vcoupling_tot)
                  << endl << "E_Bcoupling = " << scprintf("%12.10f", E_Bcoupling_tot)
                  << endl << "E_X = " << scprintf("%12.10f", E_X_tot)
                  << endl << "E_X (non canonical orbitals) = " << scprintf("%12.10f", E_X_noca_tot)
                  << endl << "E_B = " << scprintf("%12.10f", E_B_tot)
                  << endl << "E_VT = " << scprintf("%12.10f", E_VT_tot)
                  << endl << endl;
  }

  // Set beta-beta energies to alpha-alpha for closed-shell
  if (!r12world->refwfn()->spin_polarized()) {
    C_[BetaBeta] = C_[AlphaAlpha];
    emp2f12_[BetaBeta] = emp2f12_[AlphaAlpha];
    ef12_[BetaBeta] = ef12_[AlphaAlpha];
  }

  evaluated_ = true;

  return;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
