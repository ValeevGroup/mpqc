//
// mp2r12_energy_diag.cc
//
// Copyright (C) 2010 Jinmei Zhang
//
// Author: Jinmei Zhang <jmzhang@vt.edu>
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
#pragma implementation
#endif

#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <math.h>

using namespace std;
using namespace sc;

/**********************************
 * class MP2R12Energy_Diag *
 **********************************/

// Activate the computation of integrals
void MP2R12Energy_Diag::activate_ints(const std::string& occ1_act_id, const std::string& occ2_act_id,
                                      const std::string& orbs1_id, const std::string& orbs2_id,
                                      const std::string& descr_key, Ref<TwoBodyFourCenterMOIntsRuntime>& moints4_rtime,
                                      Ref<DistArray4>& i1i2i1i2_ints)
{
  const std::string i1i2i1i2_key = ParsedTwoBodyFourCenterIntKey::key(occ1_act_id, occ2_act_id,
                                                                      orbs1_id, orbs2_id,
                                                                      descr_key,
                                                                      TwoBodyIntLayout::b1b2_k1k2);
  Ref<TwoBodyMOIntsTransform> i1i2i1i2_tform = moints4_rtime->get(i1i2i1i2_key);
  i1i2i1i2_tform->compute();
  i1i2i1i2_ints = i1i2i1i2_tform->ints_acc();
  i1i2i1i2_ints->activate();
}

// Compute Y^ij_ij, Y^ij_ji, Y^ji_ij, Y^ji_ji (indicated by b1b2_k1k2)
// i.e. Y^ij_ij <=> ij_ij, Y^ij_ji <=> ij_ji
void MP2R12Energy_Diag::compute_Y(const int b1b2_k1k2, const double prefactor,
                                  const unsigned int oper_idx,
                                  Ref<DistArray4>& i1i2i1i2_ints, double* array_i1i2i1i2)
{

  // ExEnv::out0() << indent << "b1b2_k1k2:" << b1b2_k1k2 << endl;

  int nocc1_act = i1i2i1i2_ints->ni();
  int nocc2_act = i1i2i1i2_ints->nj();;
  if (b1b2_k1k2 == ji_ij || b1b2_k1k2 == ji_ji) {
    nocc1_act = i1i2i1i2_ints->nj();
    nocc2_act = i1i2i1i2_ints->ni();
  }

  // get the pair block ij or ji based on k1k2
  int i=0;
  int j=0;
  int* blk_i1 = &i;
  int* blk_i2 = &j;
  // for Y^ji_ij and Y^ji_ji, the ji pair block is needed
  if (b1b2_k1k2 == ji_ij || b1b2_k1k2 == ji_ji) {
    blk_i1 = &j;
    blk_i2 = &i;
  }

  // Now only V^ij_ij V^ij_ji X^ij_ij X^ij_ji B^ij_ji are computed,
  // so the index for array which stores V^ij_ij ... is always ij
  // ij = i * nocc2_act + j which increase one for each loop

  // increase one each time
  int idx_ij = 0;
  // increase nocc1_act each time
  int idx_ji = 0;

  int* array_idx = &idx_ij;
  int* blk_idx = &idx_ij;
  if (b1b2_k1k2 == ij_ji || b1b2_k1k2 == ji_ji)
    blk_idx = &idx_ji;

  for( ; i<nocc1_act; ++i) {
    idx_ji = i;

    for(j=0; j<nocc2_act; ++j) {

      // based on the index of bra, get the ij or ji pair block
      const double* blk = i1i2i1i2_ints->retrieve_pair_block(*blk_i1, *blk_i2, oper_idx);

      array_i1i2i1i2[*array_idx] += prefactor * blk[*blk_idx];
      ExEnv::out0() << indent << scprintf("%12.10f", array_i1i2i1i2[*array_idx]) << " ";

      ++idx_ij;
      idx_ji += nocc1_act;

      i1i2i1i2_ints->release_pair_block(*blk_i1, *blk_i2, oper_idx);
    }
    ExEnv::out0() << endl;
  }
}

// Compute (YxF)^ij_ij or (YxF)^ij_ji
// Y can be g or f
void MP2R12Energy_Diag::compute_YxF(const int b1b2_k1k2, const double prefactor,
                                    const unsigned int oper1_idx, const unsigned int oper2_idx,
                                    Ref<DistArray4>& i1i2x1x2_ints, Ref<DistArray4>& j1j2x1x2_ints,
                                    double* array_i1i2i1i2)
{
  const int norbs12 = i1i2x1x2_ints->nx() * i1i2x1x2_ints->ny();

  int nocc1_act = i1i2x1x2_ints->ni();
  int nocc2_act = i1i2x1x2_ints->nj();;
  if (b1b2_k1k2 == ji_ij || b1b2_k1k2 == ji_ji) {
    nocc1_act = i1i2x1x2_ints->nj();
    nocc2_act = i1i2x1x2_ints->ni();
  }

  // get the right pair block
  int i=0;
  int j=0;

  // Select block indices
  int* blk1_i1 = &i;
  int* blk1_i2 = &j;
  int* blk2_i1 = &i;
  int* blk2_i2 = &j;

  switch (b1b2_k1k2) {
  case ij_ij:
    // default value
  break;

  case ij_ji:
    blk2_i1 = &j;
    blk2_i2 = &i;
  break;

  case ji_ij:
    blk1_i1 = &j;
    blk1_i2 = &i;
  break;

  case ji_ji:
    blk1_i1 = &j;
    blk1_i2 = &i;
    blk2_i1 = &j;
    blk2_i2 = &i;
  break;

  default:
    ExEnv::out0() << "There is no such index";
    break;
  }

  // Now only V^ij_ij V^ij_ji X^ij_ij X^ij_ji B^ij_ji are computed,
  // so the index for array which stores V^ij_ij ... is always ij
  // ij = i * nocc2_act + j which increase one for each loop
  int array_idx = 0;

  for( ; i<nocc1_act; ++i) {
    for(j=0; j<nocc2_act; ++j) {
      const int one = 1;

      const double* blk1 = i1i2x1x2_ints->retrieve_pair_block(*blk1_i1, *blk1_i2, oper1_idx);
      const double* blk2 = j1j2x1x2_ints->retrieve_pair_block(*blk2_i1, *blk2_i2, oper2_idx);
      const double i1i2i1i2 = prefactor * F77_DDOT(&norbs12, blk1, &one, blk2, &one);

      array_i1i2i1i2[array_idx] += i1i2i1i2 ;
      ExEnv::out0() << indent << scprintf("%12.10f",array_i1i2i1i2[array_idx]) << "  ";

      i1i2x1x2_ints->release_pair_block(*blk1_i1, *blk1_i2, oper1_idx);
      j1j2x1x2_ints->release_pair_block(*blk2_i1, *blk2_i2, oper2_idx);

      ++array_idx;
    }
    ExEnv::out0() << endl;
  }
}

void MP2R12Energy_Diag::compute_VX(const int b1b2_k1k2, std::vector<std::string>& VX_output,
                                   const unsigned int oper12_idx, Ref<DistArray4>& Y12_ints,
                                   const unsigned int oper1_idx, const unsigned int oper2_idx,
                                   std::vector<Ref<DistArray4> >& Y_ints,
                                   std::vector<Ref<DistArray4> >& F_ints,
                                   double* VX_array)
{
  ExEnv::out0() << indent << "(diag) contribution" << endl;
  compute_Y(b1b2_k1k2, 1.0,
            oper12_idx,
            Y12_ints, VX_array);

  const int nterms = Y_ints.size();
  for (int i=0; i != nterms; ++i) {
    ExEnv::out0() << indent << VX_output[i] << endl;
    compute_YxF(b1b2_k1k2, -1.0,
                oper1_idx, oper2_idx,
                Y_ints[i], F_ints[i],
                VX_array);
  }
  ExEnv::out0() << endl;
}

void MP2R12Energy_Diag::accumulate_P_YxF(std::vector<std::string>& P_output,
                                         std::vector<int>& b1b2_k1k2, std::vector<double>& P_prefactor,
                                         const unsigned int oper1_idx, const unsigned int oper2_idx,
                                         std::vector<Ref<DistArray4> >& Y_ints,
                                         std::vector<Ref<DistArray4> >& F_ints,
                                         double* P)
{
  const int nterms = Y_ints.size();
  if (Y_ints.size() != F_ints.size())
    ExEnv::out0() <<  indent << "input error" << endl;

  std::vector<std::string>::iterator output = P_output.begin();
  std::vector<double>::iterator prefactor = P_prefactor.begin();

  for (int i=0; i!=nterms; ++i, ++output, ++prefactor) {
    ExEnv::out0() << indent << *output << endl;

    // P += f^ij_xy F^x_z f^zy_ij
    //    = ()^ij_ij
    compute_YxF(b1b2_k1k2[i], *prefactor,
                oper1_idx, oper2_idx,
                Y_ints[i], F_ints[i],
                P);
    ExEnv::out0() << endl;

    // P += f^ji_xy F^x_z f^zy_ji
    //    = ()^ji_ji
    ++i;
    compute_YxF(b1b2_k1k2[i], *prefactor,
                oper1_idx, oper2_idx,
                Y_ints[i], F_ints[i],
                P);
    ExEnv::out0() << endl;
  }
}


void MP2R12Energy_Diag::compute_ef12() {
  if (evaluated_)
    return;

  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
  Ref<MessageGrp> msg = r12world->world()->msg();
  int me = msg->me();
  int ntasks = msg->n();

  const Ref<R12Technology::CorrelationFactor> corrfactor = r12world->r12tech()->corrfactor();
  const bool obs_eq_ribs = r12world->obs_eq_ribs();
  const bool cabs_empty = obs_eq_ribs;
  // Only diagonal ansatz is supported
  const bool diag = r12world->r12tech()->ansatz()->diag();
  if (diag == false)
    throw ProgrammingError("only diagonal ansatz supported",__FILE__,__LINE__);

  // obtain some preliminaries
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  Ref<TwoBodyIntDescr> descr_f12f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0,0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const std::string descr_f12f12_key = moints4_rtime->descr_key(descr_f12f12);

  //
  // Evaluate pair energies:
  // distribute workload among nodes by pair index
  //
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  ExEnv::out0() << indent << "num_unique_spincases2 = " << num_unique_spincases2 << endl;

  double f12_energy_aa = 0;   // f12 correlation energy for alpha alpha part
  double f12_energy_bb = 0;
  double f12_energy_ab = 0;   // alpha beta part

  // find the largest dimension (nxn) of the matrixes
  // n is the larger number between the #n of alpha and beta electron
  // start from getting the number of alpha and beta electron in active space
  int num_alpha = 0;
  int num_beta = 0;
  int dim_matrix = 0;
  const SpinCase2 spincase_AlphAlpha = static_cast<SpinCase2> (1);      // alpha alpha case
  if (r12eval()->dim_oo(spincase_AlphAlpha).n() != 0) {

      const SpinCase1 spin1_alpha = case1(spincase_AlphAlpha);
      const Ref<OrbitalSpace>& occ_act_alpha = r12eval()->occ_act(spin1_alpha);
      num_alpha = occ_act_alpha->rank();
  }

  const SpinCase2 spincase_BetaBeta = static_cast<SpinCase2> (2);
  if (r12eval()->dim_oo(spincase_BetaBeta).n() != 0) {

      const SpinCase1 spin1_beta = case1(spincase_BetaBeta);
      const Ref<OrbitalSpace>& occ_act_beta = r12eval()->occ_act(spin1_beta);
      num_beta = occ_act_beta->rank();
  }
  // ExEnv::out0() << indent << "#n of spin1(active) = " << num_alpha << "  #n of spin2(active) = " << num_beta << endl;

  // obtain the larger number between alpha an beta electron: N
  // the dimension of all the intermediate matrixes is: NxN
  if (num_alpha >= num_beta)
    dim_matrix = num_alpha * num_alpha;
  else
    dim_matrix = num_beta * num_beta;
  ExEnv::out0() << indent << "The dimension of the V X and B matrix is " << dim_matrix << endl;

  if (dim_matrix == 0)
    ExEnv::out0() << indent << "error: no electron" << endl;

  // assign pointers to each matrix
  double* Vij_ij = new double[dim_matrix];
  double* Vij_ji = new double[dim_matrix];
  double* Xij_ij = new double[dim_matrix];
  double* Xij_ji = new double[dim_matrix];
  double* Bij_ij = new double[dim_matrix];
  double* Pij_ij = new double[dim_matrix];
  double* Qij_ij = new double[dim_matrix];
  double* Bij_ji = new double[dim_matrix];
  double* Pij_ji = new double[dim_matrix];
  double* Qij_ji = new double[dim_matrix];

  for (int spin = 0; spin < num_unique_spincases2; spin++) {

    const SpinCase2 spincase = static_cast<SpinCase2> (spin);
    if (r12eval()->dim_oo(spincase).n() == 0)
      continue;
    const SpinCase1 spin1 = case1(spincase);
    const SpinCase1 spin2 = case2(spincase);
    // Determine the spincase
    std::string spinletters = to_string(spincase);

    ExEnv::out0() << indent << "spin1 = " << spin1 << "  spin2 = " << spin2
                            << "  "<< spinletters<< endl;

    const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
    const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
    const int nocc1_act = occ1_act->rank();
    const int nocc2_act = occ2_act->rank();
    const int nocc12 = nocc1_act * nocc2_act;

    const Ref<OrbitalSpace>& vir1_act = r12eval()->vir_act(spin1);
    const Ref<OrbitalSpace>& xspace1 = r12eval()->xspace(spin1);
    const Ref<OrbitalSpace>& occ1 = r12eval()->occ(spin1);
    const Ref<OrbitalSpace>& vir1 = r12eval()->vir(spin1);
    const Ref<OrbitalSpace>& orbs1 = r12eval()->orbs(spin1);
    const Ref<OrbitalSpace>& cabs1 = r12world->cabs_space(spin1);
    const Ref<OrbitalSpace>& ribs1 = r12world->ribs_space();
    const Ref<OrbitalSpace>& Kribs1 = r12eval()->K_P_P(spin1);
    const Ref<OrbitalSpace>& Fribs1 = r12eval()->F_P_P(spin1);
    const Ref<OrbitalSpace>& fribs1 = r12eval()->F_p_p(spin1);
    const Ref<OrbitalSpace>& Focc1 = r12eval()->F_m_P(spin1);
    const Ref<OrbitalSpace>& focc1 = r12eval()->F_m_m(spin1);
    const Ref<OrbitalSpace>& forbs1 = r12eval()->F_p_A(spin1);

    const Ref<OrbitalSpace>& vir2_act = r12eval()->vir_act(spin2);
    const Ref<OrbitalSpace>& xspace2 = r12eval()->xspace(spin2);
    const Ref<OrbitalSpace>& occ2 = r12eval()->occ(spin2);
    const Ref<OrbitalSpace>& vir2 = r12eval()->vir(spin2);
    const Ref<OrbitalSpace>& orbs2 = r12eval()->orbs(spin2);
    const Ref<OrbitalSpace>& cabs2 = r12world->cabs_space(spin2);
    const Ref<OrbitalSpace>& ribs2 = r12world->ribs_space();
    const Ref<OrbitalSpace>& Kribs2 = r12eval()->K_P_P(spin2);
    const Ref<OrbitalSpace>& Fribs2 = r12eval()->F_P_P(spin2);
    const Ref<OrbitalSpace>& fribs2 = r12eval()->F_p_p(spin2);
    const Ref<OrbitalSpace>& Focc2 = r12eval()->F_m_P(spin2);
    const Ref<OrbitalSpace>& focc2 = r12eval()->F_m_m(spin2);
    const Ref<OrbitalSpace>& forbs2 = r12eval()->F_p_A(spin2);

    //
    // compute intermediates V, X, B
    //
    const TwoBodyOper::type f12eri_type = r12world->r12tech()->corrfactor()->tbint_type_f12eri();
    const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
    const TwoBodyOper::type eri_type = r12world->r12tech()->corrfactor()->tbint_type_eri();
    const TwoBodyOper::type f12f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12f12();
    const TwoBodyOper::type f12t1f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12t1f12();
    const unsigned int f12eri_idx = descr_f12->intset(f12eri_type);
    const unsigned int f12_idx = descr_f12->intset(f12_type);
    const unsigned int eri_idx = descr_f12->intset(eri_type);
    const unsigned int f12f12_idx = descr_f12f12->intset(f12f12_type);
    const unsigned int f12t1f12_idx = descr_f12f12->intset(f12t1f12_type);

    //
    // Compute the V intermediate matrix: V^ij_ij and V^ij_ji
    //
    // Alpha_beta V^ij_ij and V^ij_ji are antisymmetrized
    // Alpha_alpha or beta_beta antisymmetrized V = V^ij_ij - V^ij_ji
    // Note: '^' indicates bra, '_' indicates ket

    // this is necessary because how I wrote the compute_Y function
    std::fill(Vij_ij, Vij_ij + nocc12, 0.0);
    ExEnv::out0() << indent << spinletters << " V^ij_ij: " << endl;

    // V^ij_ij = (f12/r12)^ij_ij -
    Ref<DistArray4> i1i2i1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  occ1_act->id(), occ2_act->id(),
                  descr_f12_key, moints4_rtime,
                  i1i2i1i2_ints);

#if 0
    ExEnv::out0() << indent << "V(diag) contribution" << endl;

    for(int i1=0; i1<nocc1_act; ++i1) {
      for(int i2=0; i2<nocc2_act; ++i2) {

        const double* V_f12eri_blk = i1i2i1i2_ints->retrieve_pair_block(i1, i2, f12eri_idx);
        const int V_i1i2 = i1 * nocc2_act + i2;
        ExEnv::out0() << indent << V_i1i2 << " "<<  scprintf("%12.10f",V_f12eri_blk[V_i1i2]) << "  ";

        i1i2i1i2_ints->release_pair_block(i1, i2, f12eri_idx);
      }
      ExEnv::out0() << endl;
    }
#endif

    // store all the ints
    std::vector<Ref<DistArray4> > f12_ij_ints;
    std::vector<std::string> VX_output;

    // V^ij_ij -= g^ij_pq f^pq_ij
    Ref<DistArray4> i1i2p1p2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  orbs1->id(), orbs2->id(),
                  descr_f12_key, moints4_rtime,
                  i1i2p1p2_ints);

    f12_ij_ints.push_back(i1i2p1p2_ints);
    VX_output.push_back("diag-pq contribution");

    // V^ij_ij -= g^ij_ma' f^ma'_ij
    Ref<DistArray4> i1i2i1a2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  occ1->id(), cabs2->id(),
                  descr_f12_key, moints4_rtime,
                  i1i2i1a2_ints);

    f12_ij_ints.push_back(i1i2i1a2_ints);
    VX_output.push_back("diag-pq-ma' contribution");

    // V^ij_ij -= g^ij_a'm f^a'm_ij
    Ref<DistArray4> i1i2a1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  cabs1->id(), occ2->id(),
                  descr_f12_key, moints4_rtime,
                  i1i2a1i2_ints);

    f12_ij_ints.push_back(i1i2a1i2_ints);
    VX_output.push_back("diag-pq-ma'-a'm contribution");

    compute_VX(ij_ij, VX_output,
               f12eri_idx, i1i2i1i2_ints,
               eri_idx, f12_idx,
               f12_ij_ints, f12_ij_ints,
               Vij_ij);

    //
    // Vij_ji = V^ij_ji = g^ij_ab f^ab_ji
    ExEnv::out0() << indent << spinletters << " V^ij_ji : " << endl;
    std::fill(Vij_ji, Vij_ji + nocc12, 0.0);

    // V^ij_ji = (f12/r12)^ij_ji
    Ref<DistArray4> i1i2i2i1_ints;

    // V^ij_ji -= g^ij_pq f^pq_ji
    Ref<DistArray4> i2i1p1p2_ints;

    // V^ij_ji -= r^ij_ma' f^ma'_ji
    Ref<DistArray4> i2i1i1a2_ints;

    // V^ij_ji -= r^ij_a'm f^a'm_ji
    Ref<DistArray4> i2i1a1i2_ints;

    if (num_unique_spincases2 == 3 && spin1 != spin2){
      activate_ints(occ1_act->id(), occ2_act->id(),
                    occ2_act->id(), occ1_act->id(),
                    descr_f12_key, moints4_rtime,
                    i1i2i2i1_ints);

      activate_ints(occ2_act->id(), occ1_act->id(),
                    orbs1->id(), orbs2->id(),
                    descr_f12_key, moints4_rtime,
                    i2i1p1p2_ints);

      activate_ints(occ2_act->id(), occ1_act->id(),
                    occ1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime,
                    i2i1i1a2_ints);

      activate_ints(occ2_act->id(), occ1_act->id(),
                    cabs1->id(), occ2->id(),
                    descr_f12_key, moints4_rtime,
                    i2i1a1i2_ints);
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

    compute_VX(ij_ji, VX_output,
               f12eri_idx, i1i2i2i1_ints,
               eri_idx, f12_idx,
               f12_ij_ints, f12_ji_ints,
               Vij_ji);

    //
    // Compute X intermediate matrix: X^ij_ij and X^ij_ji
    //
    // Alpha_beta X^ij_ij and X^ij_ji are antisymmetrized
    // Alpha_alpha or beta_beta antisymmetrized X = X^ij_ij - X^ij_ji

    // Xij_ij = X^ij_ij = f^ij_ab f^ab_ij
    std::fill(Xij_ij, Xij_ij + nocc12, 0.0);
    ExEnv::out0() << indent << spinletters << " X^ij_ij : " << endl;

    // X^ij_ij += (f12f12)^ij_ij
    Ref<DistArray4> x_i1i2i1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  occ1_act->id(), occ2_act->id(),
                  descr_f12f12_key, moints4_rtime,
                  x_i1i2i1i2_ints);

    compute_VX(ij_ij, VX_output,
               f12f12_idx, x_i1i2i1i2_ints,
               f12_idx, f12_idx,
               f12_ij_ints, f12_ij_ints,
               Xij_ij);

    //
    // Xij_ji = X^ij_ji = f^ij_ab f^ab_ji
    std::fill(Xij_ji, Xij_ji + nocc12, 0.0);
    ExEnv::out0() << indent << spinletters << " X^ij_ji : " << endl;

    // X^ij_ji = (f12f12)^ij_ji
    Ref<DistArray4> x_i1i2i2i1_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2){
      activate_ints(occ1_act->id(), occ2_act->id(),
                    occ2_act->id(), occ1_act->id(),
                    descr_f12f12_key, moints4_rtime,
                    x_i1i2i2i1_ints);
    } else {
      x_i1i2i2i1_ints = x_i1i2i1i2_ints;
    }

    compute_VX(ij_ji,  VX_output,
               f12f12_idx, x_i1i2i2i1_ints,
               f12_idx, f12_idx,
               f12_ij_ints, f12_ji_ints,
               Xij_ji);

    //
    // Compute B intermediate matrix
    //
    // B^ij_ij
    // alpha beta case: antisymmetrized
    // alpha alpha, beta beta case: B^ij_ij - B^ij_ji (antisymmetrized)
    //

    // This is for alpha apha, beta beta and alpha beta cases
    //   alpha beta case is indicated by the following:
    //   B^i(alpha)j(beta)_^i(alpha)j(beta)
    // = R^i(alpha)j(beta)_a(alpha)b(beta) * f^a(alpha)_a(alpha) * R^a(alpha)b(beta)_i(alpha)j(beta)
    std::fill(Bij_ij, Bij_ij + nocc12, 0.0);
    ExEnv::out0() << indent << spinletters << " B^ij_ij: " << endl;

    // B^ij_ij += (f12t1f12)^ij_ij
    ExEnv::out0() << indent << spinletters << " B(diag) contribution" << endl;

    Ref<DistArray4> b_i1i2i1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  occ1_act->id(), occ2_act->id(),
                  descr_f12f12_key, moints4_rtime,
                  b_i1i2i1i2_ints);
    compute_Y(ij_ij, 1.0,
              f12t1f12_idx,
              b_i1i2i1i2_ints, Bij_ij);
    ExEnv::out0() << endl;

    // P part of B intermediate
    std::fill(Pij_ij, Pij_ij + nocc12, 0.0);

    // Store all the ints, prefactor, index for b1b2k1k2 and output
    // for each term
    std::vector<Ref<DistArray4> > Pijij_f12_ints;
    std::vector<Ref<DistArray4> > Pijij_fx12_ints;
    std::vector<double> P_prefactors;
    std::vector<int> Pijij_idx;
    std::vector<std::string> Pijij_output;

    // P +=  f^ij_PQ K^Q_R f^PR_ij    (antisymmetrized)

    //    =   f^ij_QP K^Q_R f^RP_ij - f^ij_QP K^Q_R f^RP_ji
    //      - f^ji_QP K^Q_R f^RP_ij + f^ji_QP K^Q_R f^RP_ji
    //
    // when spin1 != spin2: P+= f^ij_QP K^Q_R f^RP_ij + f^ji_QP K^Q_R f^RP_ji
    //                        = f^ij_QP f^Q_K P _ij + f^ji_QP f^Q_K P _ji ;

    Ref<DistArray4> b_i1i2P1P2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  ribs1->id(), ribs2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2P1P2_ints);

    Ref<DistArray4> b_i1i2PK1P2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  Kribs1->id(), ribs2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2PK1P2_ints);

    // P += f^ji_QP f^Q_K P _ji
    Ref<DistArray4> b_i2i1P1P2_ints;
    Ref<DistArray4> b_i2i1PK1P2_ints;
    if(num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(),
                    ribs1->id(), ribs2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1P1P2_ints);
      activate_ints(occ2_act->id(), occ1_act->id(),
                    Kribs1->id(), ribs2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1PK1P2_ints);
      //ExEnv::out0() << indent << "activate: b_i2i1P1P2_ints, b_i2i1PK1P2_ints" << endl;
    } else {
      //spin1 == spin2 or close shell
      b_i2i1P1P2_ints = b_i1i2P1P2_ints;
      b_i2i1PK1P2_ints = b_i1i2PK1P2_ints;

    }

    Pijij_output.push_back("P(f^ij_PQ K^Q_R f^PR_ij)");
    P_prefactors.push_back(1.0);
    Pijij_f12_ints.push_back(b_i1i2P1P2_ints);
    Pijij_fx12_ints.push_back(b_i1i2PK1P2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1P1P2_ints);
    Pijij_fx12_ints.push_back(b_i2i1PK1P2_ints);
    Pijij_idx.push_back(ji_ji);

    // P += f^ij_Pm f^P_F m _ij + f^ji_Pm f^P_F m _ji
    // P += f^ij_Pm f^P_F m _ij
    Ref<DistArray4> b_i1i2P1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  ribs1->id(), occ2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2P1i2_ints);

    Ref<DistArray4> b_i1i2PF1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  Fribs1->id(), occ2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2PF1i2_ints);

    // P += f^ji_Pm F^P_Q f^Qm_ji
    Ref<DistArray4> b_i2i1P1i2_ints;
    Ref<DistArray4> b_i2i1PF1i2_ints;
    if(num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(),
                    ribs1->id(), occ2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1P1i2_ints);
      activate_ints(occ2_act->id(), occ1_act->id(),
                    Fribs1->id(), occ2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1PF1i2_ints);
      //ExEnv::out0() << indent << "activate: b_i2i1P1i2_ints, b_i2i1PF1i2_ints" << endl;
    } else {
      //spin1 == spin2 or close shell
      b_i2i1P1i2_ints = b_i1i2P1i2_ints;
      b_i2i1PF1i2_ints = b_i1i2PF1i2_ints;
    }

    Pijij_output.push_back("P(including f^ij_Pm F^P_Q f^Qm_ij)");
    P_prefactors.push_back(1.0);
    Pijij_f12_ints.push_back(b_i1i2P1i2_ints);
    Pijij_fx12_ints.push_back(b_i1i2PF1i2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1P1i2_ints);
    Pijij_fx12_ints.push_back(b_i2i1PF1i2_ints);
    Pijij_idx.push_back(ji_ji);


    // P += f^ij_pe F^p_q f^qe_ij + f^ji_pe F^p_q f^qe_ji
    //
    //    = f^ij_pe f^p_F e _ij  + ...

    // P += f^ij_pe F^p_q f^qe_ij
    Ref<DistArray4> b_i1i2p1e2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  orbs1->id(), vir2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2p1e2_ints);

    Ref<DistArray4> b_i1i2pF1e2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  fribs1->id(), vir2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2pF1e2_ints);

    // P += f^ji_pe F^p_q f^qe_ji
    Ref<DistArray4> b_i2i1p1e2_ints;
    Ref<DistArray4> b_i2i1pF1e2_ints;
    if(num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(),
                    orbs1->id(), vir2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1p1e2_ints);
      activate_ints(occ2_act->id(), occ1_act->id(),
                    fribs1->id(), vir2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1pF1e2_ints);
      //ExEnv::out0() << indent << "activate: b_i2i1p1e2_ints, b_i2i1pF1e2_ints" << endl;
    } else {
      //spin1 == spin2 or close shell
      b_i2i1p1e2_ints = b_i1i2p1e2_ints;
      b_i2i1pF1e2_ints = b_i1i2pF1e2_ints;
    }

    Pijij_output.push_back("P(including f^ij_pe F^p_q f^qe_ij)");
    P_prefactors.push_back(1.0);
    Pijij_f12_ints.push_back(b_i1i2p1e2_ints);
    Pijij_fx12_ints.push_back(b_i1i2pF1e2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1p1e2_ints);
    Pijij_fx12_ints.push_back(b_i2i1pF1e2_ints);
    Pijij_idx.push_back(ji_ji);

    // P += f^ij_mA F^m_P f^PA_ij + f^ji_mA F^m_P f^PA_ji + h.c.
    //    = 2.0 * (f^ij_mA f^m_F A _ij + ...)

    // P += 2 f^ij_mA F^m_P f^PA_ij
    Ref<DistArray4> b_i1i2m1A2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  occ1->id(), cabs2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2m1A2_ints);

    Ref<DistArray4> b_i1i2mF1A2_key;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  Focc1->id(), cabs2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2mF1A2_key);

    // P += 2 f^ji_mA F^m_P f^PA_ji
    Ref<DistArray4> b_i2i1m1A2_ints;
    Ref<DistArray4> b_i2i1mF1A2_key;
    if(num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(),
                    occ1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1m1A2_ints);
      activate_ints(occ2_act->id(), occ1_act->id(),
                    Focc1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1mF1A2_key);
      //ExEnv::out0() << indent << "activate: b_i2i1m1A2_ints, b_i2i1mF1A2_key" << endl;
    } else {
      b_i2i1m1A2_ints = b_i1i2m1A2_ints;
      b_i2i1mF1A2_key = b_i1i2mF1A2_key;
    }

    Pijij_output.push_back("P(including 2 f^ij_mA F^m_P f^PA_ij)");
    P_prefactors.push_back(2.0);
    Pijij_f12_ints.push_back(b_i1i2m1A2_ints);
    Pijij_fx12_ints.push_back(b_i1i2mF1A2_key);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1m1A2_ints);
    Pijij_fx12_ints.push_back(b_i2i1mF1A2_key);
    Pijij_idx.push_back(ji_ji);

    // P += f^ij_Ae F^A_p f^pe_ij + f^ji_Ae F^A_p f^pe_ji + h.c.
    //    = 2 * (f^ij_Ae f^A_F e _ij ...)

    // P += 2 f^ij_Ae F^A_p f^pe_ij
    Ref<DistArray4> b_i1i2A1e2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  orbs1->id(), vir2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2A1e2_ints);

    Ref<DistArray4> b_i1i2AF1e2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  forbs1->id(), vir2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2AF1e2_ints);

    // P += 2 f^ji_Ae F^A_p f^pe_ji
    Ref<DistArray4> b_i2i1A1e2_ints;
    Ref<DistArray4> b_i2i1AF1e2_ints;
    if(num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(),
                    orbs1->id(), vir2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1A1e2_ints);
      activate_ints(occ2_act->id(), occ1_act->id(),
                    forbs1->id(), vir2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1AF1e2_ints);
      //ExEnv::out0() << indent << "activate: b_i2i1A1e2_ints, b_i2i1AF1e2_ints" << endl;
    } else {
      //spin1 == spin2 or close shell
      b_i2i1A1e2_ints = b_i1i2A1e2_ints;
      b_i2i1AF1e2_ints = b_i1i2AF1e2_ints;
    }

    Pijij_output.push_back("P(including 2 f^ij_Ae F^A_p f^pe_ij)");
    P_prefactors.push_back(2.0);
    Pijij_f12_ints.push_back(b_i1i2A1e2_ints);
    Pijij_fx12_ints.push_back(b_i1i2AF1e2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1A1e2_ints);
    Pijij_fx12_ints.push_back(b_i2i1AF1e2_ints);
    Pijij_idx.push_back(ji_ji);

    // P -= f^ij_mA F^m_n f^nA_ij + f^ji_mA F^m_n f^nA_ji
    //    = f^ij_mA f^m_F A _ij ...

    // P -= f^ij_mA F^m_n f^nA_ij
    // f^ij_mA is already computed in the previous step
    Ref<DistArray4> b_i1i2mf1A2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  focc1->id(), cabs2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2mf1A2_ints);

    // P -= f^ji_mA F^m_n f^nA_ji
    // f^ji_mA is already computed in the previous step
    Ref<DistArray4> b_i2i1mf1A2_ints;
    if(num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(),
                    focc1->id(), cabs2->id(),
                    descr_f12_key, moints4_rtime,
                    b_i2i1mf1A2_ints);
      //ExEnv::out0() << indent << "activate: b_i2i1mf1A2_ints" << endl;
    } else {
      //spin1 == spin2 or close shell
      b_i2i1mf1A2_ints = b_i1i2mf1A2_ints;
    }

    Pijij_output.push_back("P(including f^ij_mA F^m_n f^nA_ij)");
    P_prefactors.push_back(-1.0);
    Pijij_f12_ints.push_back(b_i1i2m1A2_ints);
    Pijij_fx12_ints.push_back(b_i1i2mf1A2_ints);
    Pijij_idx.push_back(ij_ij);

    Pijij_f12_ints.push_back(b_i2i1m1A2_ints);
    Pijij_fx12_ints.push_back(b_i2i1mf1A2_ints);
    Pijij_idx.push_back(ji_ji);

    accumulate_P_YxF(Pijij_output,
                     Pijij_idx, P_prefactors,
                     f12_idx, f12_idx,
                     Pijij_f12_ints, Pijij_fx12_ints,
                     Pij_ij);
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

    ExEnv::out0() << indent << spinletters << " Q^ij_ij : " << endl;
    std::fill(Qij_ij, Qij_ij + nocc12, 0.0);

    ExEnv::out0() << indent << "Q += (f12f12)^ij_Pj (f+k)^P_i " << endl;
    Ref<DistArray4> q_i1i2hi1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  r12eval_->hj_i_P(spin1)->id(), occ2_act->id(),
                  descr_f12f12_key, moints4_rtime,
                  q_i1i2hi1i2_ints);
    compute_Y(ij_ij, 1.0,
              f12f12_idx,
              q_i1i2hi1i2_ints, Qij_ij);

    ExEnv::out0() << indent << "Q += (f12f12)^ji_Pi (f+k)^P_j " << endl;
    Ref<DistArray4> q_i2i1hi2i1_ints;
    if(num_unique_spincases2 == 3 && spin1 != spin2) {
      activate_ints(occ2_act->id(), occ1_act->id(),
                    r12eval_->hj_i_P(spin2)->id(), occ1_act->id(),
                    descr_f12f12_key, moints4_rtime,
                    q_i2i1hi2i1_ints);
    } else {
      // spin1 == spin2 or close shell
      q_i2i1hi2i1_ints = q_i1i2hi1i2_ints;
    }
    compute_Y(ji_ji, 1.0,
              f12f12_idx,
              q_i2i1hi2i1_ints, Qij_ij);
    ExEnv::out0() << endl;

      // B^ij_ij = B^ij_ij - P + Q
      ExEnv::out0() << indent << spinletters << " Bij_ij : " << endl;
#if 1
      for(int i1=0; i1<nocc1_act; ++i1) {
        for(int i2=0; i2<nocc2_act; ++i2) {

          const int i1i2 = i1 * nocc2_act + i2;
          Bij_ij[i1i2] += - Pij_ij[i1i2] + Qij_ij[i1i2];
          ExEnv::out0() << indent << scprintf("%12.10f",Bij_ij[i1i2])<< "  ";
        }
        ExEnv::out0() << endl;
      }
      ExEnv::out0() << endl;
#endif

      //
      // B^ij_ji
      std::fill(Bij_ji, Bij_ji + nocc12, 0.0);
      ExEnv::out0() << indent << spinletters << " B^ij_ji : " << endl;

      // B^ij_ji += (f12t1f12)^ij_ji
      ExEnv::out0() << indent << spinletters << " B(diag) contribution" << endl;
      Ref<DistArray4> b_i1i2i2i1_ints;
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        activate_ints(occ1_act->id(), occ2_act->id(),
                      occ2_act->id(), occ1_act->id(),
                      descr_f12f12_key, moints4_rtime,
                      b_i1i2i2i1_ints);
      } else {
        b_i1i2i2i1_ints = b_i1i2i1i2_ints;;
      }
      compute_Y(ij_ji, 1.0,
                f12t1f12_idx,
                b_i1i2i2i1_ints, Bij_ji);
      ExEnv::out0() << endl;

      // P part of B intermediate
      std::fill(Pij_ji, Pij_ji + nocc12, 0.0);

      // Store all the ints, prefactor, index for b1b2k1k2 and output
      // for each term
      // Pijji_f12_ints =  Pijij_f12_ints;
      // Pijji_prefactprs = P_prefactors;
      std::vector<Ref<DistArray4> > Pijji_fx12_ints;
      std::vector<int> Pijji_idx;
      std::vector<std::string> Pijji_output;

      Pijji_output.push_back("P(f^ij_PQ K^Q_R f^PR_ji)");
      Pijji_fx12_ints.push_back(b_i2i1PK1P2_ints);
      Pijji_idx.push_back(ij_ji);
      Pijji_fx12_ints.push_back(b_i1i2PK1P2_ints);
      Pijji_idx.push_back(ji_ij);

      Pijji_output.push_back("P(including f^ij_Pm F^P_Q f^Qm_ji)");
      Pijji_fx12_ints.push_back(b_i2i1PF1i2_ints);
      Pijji_idx.push_back(ij_ji);
      Pijji_fx12_ints.push_back(b_i1i2PF1i2_ints);
      Pijji_idx.push_back(ji_ij);

      Pijji_output.push_back("P(including f^ij_pe F^p_q f^qe_ji)");
      Pijji_fx12_ints.push_back(b_i2i1pF1e2_ints);
      Pijji_idx.push_back(ij_ji);
      Pijji_fx12_ints.push_back(b_i1i2pF1e2_ints);
      Pijji_idx.push_back(ji_ij);

      Pijji_output.push_back("P(including 2 f^ij_mA F^m_P f^PA_ji)");
      Pijji_fx12_ints.push_back(b_i2i1mF1A2_key);
      Pijji_idx.push_back(ij_ji);
      Pijji_fx12_ints.push_back(b_i1i2mF1A2_key);
      Pijji_idx.push_back(ji_ij);

      Pijji_output.push_back("P(including 2 f^ij_Ae F^A_p f^pe_ji)");
      Pijji_fx12_ints.push_back(b_i2i1AF1e2_ints);
      Pijji_idx.push_back(ij_ji);
      Pijji_fx12_ints.push_back(b_i1i2AF1e2_ints);
      Pijji_idx.push_back(ji_ij);

      Pijji_output.push_back("P(including f^ij_mA F^m_n f^nA_ji)");
      Pijji_fx12_ints.push_back(b_i2i1mf1A2_ints);
      Pijji_idx.push_back(ij_ji);
      Pijji_fx12_ints.push_back(b_i1i2mf1A2_ints);
      Pijji_idx.push_back(ji_ij);

      accumulate_P_YxF(Pijji_output,
                       Pijji_idx, P_prefactors,
                       f12_idx, f12_idx,
                       Pijij_f12_ints, Pijji_fx12_ints,
                       Pij_ji);

      // Qij_ji += (f12f12)^ij_Pi (f+k)^P_j + (f12f12)^ij_jP (f+k)^P_i
      //         = (f12f12)^ij _j_hJ i + (f12f12)^ji _i_hJ j
      ExEnv::out0() << indent << spinletters << " Q^ij_ji : " << endl;
      std::fill(Qij_ji, Qij_ji + nocc12, 0.0);

      Ref<DistArray4> q_i1i2hi2i1_ints;
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        activate_ints(occ1_act->id(), occ2_act->id(),
                      r12eval_->hj_i_P(spin2)->id(), occ1_act->id(),
                      descr_f12f12_key, moints4_rtime,
                      q_i1i2hi2i1_ints);
      } else {
        q_i1i2hi2i1_ints = q_i2i1hi2i1_ints;
      }
      compute_Y(ij_ji, 1.0,
                f12f12_idx,
                q_i1i2hi2i1_ints, Qij_ji);

      Ref<DistArray4> q_i2i1hi1i2_ints;
      if(num_unique_spincases2 == 3 && spin1 != spin2) {
        activate_ints(occ2_act->id(), occ1_act->id(),
                      r12eval_->hj_i_P(spin1)->id(), occ2_act->id(),
                      descr_f12f12_key, moints4_rtime,
                      q_i2i1hi1i2_ints);
      } else {
        q_i2i1hi1i2_ints = q_i1i2hi1i2_ints;
      }
      compute_Y(ji_ij, 1.0,
                f12f12_idx,
                q_i2i1hi1i2_ints, Qij_ji);
      ExEnv::out0() << endl;

      // B^ij_ji = B^ij_ji - P^ij_ji + Q^ij_ji
      ExEnv::out0() << indent << spinletters << " Bij_ji : " << endl;
#if 1
      for(int i1=0; i1<nocc1_act; ++i1) {
        for(int i2=0; i2<nocc2_act; ++i2) {

          const int i1i2 = i1 * nocc2_act + i2;
          Bij_ji[i1i2] += - Pij_ji[i1i2] + Qij_ji[i1i2];
          ExEnv::out0() << indent << scprintf("%12.10f",Bij_ji[i1i2])<< "  ";
        }
        ExEnv::out0() << endl;
      }
      ExEnv::out0() << endl;
#endif

      //
      //  B^i(beta)j(alpha)_^i(beta)j(alpha)
      //
      //   R^i(alpha)j(beta)_a(alpha)b(beta) * f^b(beta)_b(beta) * R^a(alpha)b(beta)_i(alpha)j(beta)
      // =>  electron 1 <-> electron 2  i<->j  a<->b
      // = R^i(beta)j(alpha)_a(beta)r(alpha) * f^a(beta)_a(beta) * R^a(beta)b(alpha)_i(beta)j(alpha)
      //
      double* Bij_ij_beta = NULL;
      double* Bij_ji_beta = NULL;

      if (num_unique_spincases2 == 3 && spin1 != spin2) {

          ExEnv::out0() << indent << spinletters << " B^ij_ij (beta alpha case): " << endl;
          Bij_ij_beta = new double[nocc12];
          std::fill(Bij_ij_beta, Bij_ij_beta + nocc12, 0.0);

          // B^ij_ij += (f12t1f12)^ij_ij
          ExEnv::out0() << indent << spinletters << " B(diag) contribution (beta alpha case)" << endl;

          Ref<DistArray4> b_i2i1i2i1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        occ2_act->id(), occ1_act->id(),
                        descr_f12f12_key, moints4_rtime,
                        b_i2i1i2i1_ints);
          compute_Y(ij_ij, 1.0,
                    f12t1f12_idx,
                    b_i2i1i2i1_ints, Bij_ij_beta);
          ExEnv::out0() << endl;

          // P part of B intermediate
          double* Pij_ij_beta = new double[nocc12];
          std::fill(Pij_ij_beta, Pij_ij_beta + nocc12, 0.0);

          // Store all the ints for b1b2k1k2 and output for each term
          // prefactor and  index are the same as Pij_ij (alpha case)
          std::vector<Ref<DistArray4> > Pijij_beta_f12_ints;
          std::vector<Ref<DistArray4> > Pijij_beta_fx12_ints;
          std::vector<std::string> Pijij_beta_output;

          // P +=  f^ij_PQ K^Q_R f^PR_ij    (antisymmetrized)

          //    =   f^ij_QP K^Q_R f^RP_ij - f^ij_QP K^Q_R f^RP_ji
          //      - f^ji_QP K^Q_R f^RP_ij + f^ji_QP K^Q_R f^RP_ji
          //
          // when spin1 != spin2: P+= f^ij_QP K^Q_R f^RP_ij + f^ji_QP K^Q_R f^RP_ji
          //                        = f^ij_QP f^Q_K P _ij + f^ji_QP f^Q_K P _ji ;

          Ref<DistArray4> b_i2i1P2P1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        ribs2->id(), ribs1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1P2P1_ints);

          Ref<DistArray4> b_i2i1PK2P1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        Kribs2->id(), ribs1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1PK2P1_ints);

          // P += f^ji_QP f^Q_K P _ji
          Ref<DistArray4> b_i1i2P2P1_ints;
          Ref<DistArray4> b_i1i2PK2P1_ints;
          activate_ints(occ1_act->id(), occ2_act->id(),
                        ribs2->id(), ribs1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2P2P1_ints);
          activate_ints(occ1_act->id(), occ2_act->id(),
                        Kribs2->id(), ribs1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2PK2P1_ints);

          Pijij_beta_output.push_back("P(f^ij_PQ K^Q_R f^PR_ij) (beta alpha case)");
          Pijij_beta_f12_ints.push_back(b_i2i1P2P1_ints);
          Pijij_beta_fx12_ints.push_back(b_i2i1PK2P1_ints);

          Pijij_beta_f12_ints.push_back(b_i1i2P2P1_ints);
          Pijij_beta_fx12_ints.push_back(b_i1i2PK2P1_ints);

          // P += f^ij_Pm f^P_F m _ij + f^ji_Pm f^P_F m _ji
          // P += f^ij_Pm f^P_F m _ij
          Ref<DistArray4> b_i2i1P2i1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        ribs2->id(), occ1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1P2i1_ints);

          Ref<DistArray4> b_i2i1PF2i1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        Fribs2->id(), occ1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1PF2i1_ints);

          // P += f^ji_Pm F^P_Q f^Qm_ji
          Ref<DistArray4> b_i1i2P2i1_ints;
          Ref<DistArray4> b_i1i2PF2i1_ints;
          activate_ints(occ1_act->id(), occ2_act->id(),
                        ribs2->id(), occ1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2P2i1_ints);
          activate_ints(occ1_act->id(), occ2_act->id(),
                        Fribs2->id(), occ1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2PF2i1_ints);

          Pijij_beta_output.push_back("P(including f^ij_Pm F^P_Q f^Qm_ij) (beta alpha case)");
          Pijij_beta_f12_ints.push_back(b_i2i1P2i1_ints);
          Pijij_beta_fx12_ints.push_back(b_i2i1PF2i1_ints);

          Pijij_beta_f12_ints.push_back(b_i1i2P2i1_ints);
          Pijij_beta_fx12_ints.push_back(b_i1i2PF2i1_ints);

          // P += f^ij_pe F^p_q f^qe_ij + f^ji_pe F^p_q f^qe_ji
          //
          //    = f^ij_pe f^p_F e _ij  + ...

          // P += f^ij_pe F^p_q f^qe_ij
          Ref<DistArray4> b_i2i1p2e1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        orbs2->id(), vir1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1p2e1_ints);

          Ref<DistArray4> b_i2i1pF2e1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        fribs2->id(), vir1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1pF2e1_ints);

          // P += f^ji_pe F^p_q f^qe_ji
          Ref<DistArray4> b_i1i2p2e1_ints;
          Ref<DistArray4> b_i1i2pF2e1_ints;
          activate_ints(occ1_act->id(), occ2_act->id(),
                        orbs2->id(), vir1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2p2e1_ints);
          activate_ints(occ1_act->id(), occ2_act->id(),
                        fribs2->id(), vir1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2pF2e1_ints);

          Pijij_beta_output.push_back("P(including f^ij_pe F^p_q f^qe_ij) (beta alpha case)");
          Pijij_beta_f12_ints.push_back(b_i2i1p2e1_ints);
          Pijij_beta_fx12_ints.push_back(b_i2i1pF2e1_ints);

          Pijij_beta_f12_ints.push_back(b_i1i2p2e1_ints);
          Pijij_beta_fx12_ints.push_back(b_i1i2pF2e1_ints);

          // P += f^ij_mA F^m_P f^PA_ij + f^ji_mA F^m_P f^PA_ji + h.c.
          //    = 2.0 * (f^ij_mA f^m_F A _ij + ...)

          // P += 2 f^ij_mA F^m_P f^PA_ij
          Ref<DistArray4> b_i2i1m2A1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        occ2->id(), cabs1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1m2A1_ints);

          Ref<DistArray4> b_i2i1mF2A1_key;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        Focc2->id(), cabs1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1mF2A1_key);

          // P += 2 f^ji_mA F^m_P f^PA_ji
          Ref<DistArray4> b_i1i2m2A1_ints;
          Ref<DistArray4> b_i1i2mF2A1_key;
          activate_ints(occ1_act->id(), occ2_act->id(),
                        occ2->id(), cabs1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2m2A1_ints);
          activate_ints(occ1_act->id(), occ2_act->id(),
                        Focc2->id(), cabs1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2mF2A1_key);

          Pijij_beta_output.push_back("P(including 2 f^ij_mA F^m_P f^PA_ij) (beta alpha case)");
          Pijij_beta_f12_ints.push_back(b_i2i1m2A1_ints);
          Pijij_beta_fx12_ints.push_back(b_i2i1mF2A1_key);

          Pijij_beta_f12_ints.push_back(b_i1i2m2A1_ints);
          Pijij_beta_fx12_ints.push_back(b_i1i2mF2A1_key);

          // P += f^ij_Ae F^A_p f^pe_ij + f^ji_Ae F^A_p f^pe_ji + h.c.
          //    = 2 * (f^ij_Ae f^A_F e _ij ...)

          // P += 2 f^ij_Ae F^A_p f^pe_ij
          Ref<DistArray4> b_i2i1A2e1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        orbs2->id(), vir1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1A2e1_ints);

          Ref<DistArray4> b_i2i1AF2e1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        forbs2->id(), vir1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1AF2e1_ints);

          // P += 2 f^ji_Ae F^A_p f^pe_ji
          Ref<DistArray4> b_i1i2A2e1_ints;
          Ref<DistArray4> b_i1i2AF2e1_ints;
          activate_ints(occ1_act->id(), occ2_act->id(),
                        orbs2->id(), vir1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2A2e1_ints);
          activate_ints(occ1_act->id(), occ2_act->id(),
                        forbs2->id(), vir1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2AF2e1_ints);

          Pijij_beta_output.push_back("P(including 2 f^ij_Ae F^A_p f^pe_ij) (beta alpha case)");
          Pijij_beta_f12_ints.push_back(b_i2i1A2e1_ints);
          Pijij_beta_fx12_ints.push_back(b_i2i1AF2e1_ints);

          Pijij_beta_f12_ints.push_back(b_i1i2A2e1_ints);
          Pijij_beta_fx12_ints.push_back(b_i1i2AF2e1_ints);

          // P -= f^ij_mA F^m_n f^nA_ij + f^ji_mA F^m_n f^nA_ji
          //    = f^ij_mA f^m_F A _ij ...

          // P -= f^ij_mA F^m_n f^nA_ij
          // f^ij_mA is already computed in the previous step
          Ref<DistArray4> b_i2i1mf2A1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        focc2->id(), cabs1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i2i1mf2A1_ints);

          // P -= f^ji_mA F^m_n f^nA_ji
          // f^ji_mA is already computed in the previous step
          Ref<DistArray4> b_i1i2mf2A1_ints;
          activate_ints(occ1_act->id(), occ2_act->id(),
                        focc2->id(), cabs1->id(),
                        descr_f12_key, moints4_rtime,
                        b_i1i2mf2A1_ints);

          Pijij_beta_output.push_back("P(including f^ij_mA F^m_n f^nA_ij) (beta alpha case)");
          Pijij_beta_f12_ints.push_back(b_i2i1m2A1_ints);
          Pijij_beta_fx12_ints.push_back(b_i2i1mf2A1_ints);

          Pijij_beta_f12_ints.push_back(b_i1i2m2A1_ints);
          Pijij_beta_fx12_ints.push_back(b_i1i2mf2A1_ints);

          accumulate_P_YxF(Pijij_beta_output,
                           Pijij_idx, P_prefactors,
                           f12_idx, f12_idx,
                           Pijij_beta_f12_ints, Pijij_beta_fx12_ints,
                           Pij_ij_beta);
          //
          // Q += 1/2[  (f12f12)^ij_Pj (f+k)^P_i - (f12f12)^ji_Pj (f+k)^P_i
          //          + (f12f12)^ij_iP (f+k)^P_j - (f12f12)^ji_iP (f+k)^P_j
          //          + h.c.]
          //    =  (f12f12)^ij_Pj (f+k)^P_i - (f12f12)^ji_Pj (f+k)^P_i
          //     + (f12f12)^ji_Pi (f+k)^P_j - (f12f12)^ij_Pi (f+k)^P_j
          //
          // spin1 != spin2: Q += (f12f12)^ij_Pj (f+k)^P_i + (f12f12)^ji_Pi (f+k)^P_j ;
          //
          // Qij_ij (alpha beta case) = Qji_ji (beta alpha case)
          ExEnv::out0() << indent << spinletters
                                  << " Q^ij_ij (beta alpha case) = Qji_ji (alpha beta case) "
                                  << endl << endl;

#if 0
          // The following codes is for testing
          ExEnv::out0() << indent << "Q^ij_ij (beta alpha case): " << endl;
          double* Qij_ij_beta = new double[nocc12];
          std::fill(Qij_ij_beta, Qij_ij_beta + nocc12, 0.0);

          ExEnv::out0() << indent << "Q += (f12f12)^ij_Pj (f+k)^P_i (beta alpha case)" << endl;
          Ref<DistArray4> q_i2i1hi2i1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        r12eval_->hj_i_P(spin2)->id(), occ1_act->id(),
                        descr_f12f12_key, moints4_rtime,
                        q_i2i1hi2i1_ints);
          compute_Y(ij_ij, 1.0,
                    f12f12_idx,
                    q_i2i1hi2i1_ints, Qij_ij_beta);

          ExEnv::out0() << indent << "Q += (f12f12)^ji_Pi (f+k)^P_j (beta alpha case)" << endl;
          Ref<DistArray4> q_i1i2hi1i2_ints;
          activate_ints(occ1_act->id(), occ2_act->id(),
                        r12eval_->hj_i_P(spin1)->id(), occ2_act->id(),
                        descr_f12f12_key, moints4_rtime,
                        q_i1i2hi1i2_ints);
          compute_Y(ji_ji, 1.0,
                    f12f12_idx,
                    q_i1i2hi1i2_ints, Qij_ij_beta);
          ExEnv::out0() << endl;
#endif

            // B^ij_ij = B^ij_ij - P + Q
            ExEnv::out0() << indent << spinletters << " Bij_ij (beta alpha case): " << endl;
      #if 1
            for(int i1=0; i1<nocc2_act; ++i1) {
              for(int i2=0; i2<nocc1_act; ++i2) {

                const int i1i2 = i1 * nocc1_act + i2;
                const int i2i1 = i2 * nocc2_act + i1;
                Bij_ij_beta[i1i2] += - Pij_ij_beta[i1i2] + Qij_ij[i2i1];
                ExEnv::out0() << indent << scprintf("%12.10f",Bij_ij_beta[i1i2])<< "  ";
              }
              ExEnv::out0() << endl;
            }
            ExEnv::out0() << endl;
      #endif
            delete[] Pij_ij_beta;
            Pij_ij_beta = NULL;

            //
            // B^ij_ji
            // B^i(beta)j(alpha)_j(alpha)i(beta)
            ExEnv::out0() << indent << spinletters << " B^ij_ji (beta alpha case): " << endl;
            Bij_ji_beta = new double[nocc12];
            std::fill(Bij_ji_beta, Bij_ji_beta + nocc12, 0.0);

            // B^ij_ji += (f12t1f12)^ij_ji
            ExEnv::out0() << indent << spinletters << " B(diag) contribution (beta alpha case):" << endl;
            Ref<DistArray4> b_i2i1i1i2_ints;
            activate_ints(occ2_act->id(), occ1_act->id(),
                          occ1_act->id(), occ2_act->id(),
                          descr_f12f12_key, moints4_rtime,
                          b_i2i1i1i2_ints);
            compute_Y(ij_ji, 1.0,
                      f12t1f12_idx,
                      b_i2i1i1i2_ints, Bij_ji_beta);
            ExEnv::out0() << endl;

            // P part of B intermediate
            double* Pij_ji_beta = new double[nocc12];
            std::fill(Pij_ji_beta, Pij_ji_beta + nocc12, 0.0);

            // Store all the ints, prefactor, index for b1b2k1k2 and output
            // for each term
            // Pijji_beta_f12_ints =  Pijij_beta_f12_ints;
            // Pijji_beta_prefactprs = P_prefactors;
            // Pijji_beta_idx = Pijji_idx
            std::vector<Ref<DistArray4> > Pijji_beta_fx12_ints;
            std::vector<std::string> Pijji_beta_output;

            Pijji_beta_output.push_back("P(f^ij_PQ K^Q_R f^PR_ji) (beta alpha case)");
            Pijji_beta_fx12_ints.push_back(b_i1i2PK2P1_ints);
            // Pijji_idx.push_back(ij_ji);
            Pijji_beta_fx12_ints.push_back(b_i2i1PK2P1_ints);
            // Pijji_idx.push_back(ji_ij);

            Pijji_beta_output.push_back("P(including f^ij_Pm F^P_Q f^Qm_ji) (beta alpha case)");
            Pijji_beta_fx12_ints.push_back(b_i1i2PF2i1_ints);
            // Pijji_idx.push_back(ij_ji);
            Pijji_beta_fx12_ints.push_back(b_i2i1PF2i1_ints);
            // Pijji_idx.push_back(ji_ij);

            Pijji_beta_output.push_back("P(including f^ij_pe F^p_q f^qe_ji) (beta alpha case)");
            Pijji_beta_fx12_ints.push_back(b_i1i2pF2e1_ints);
            // Pijji_idx.push_back(ij_ji);
            Pijji_beta_fx12_ints.push_back(b_i2i1pF2e1_ints);
            // Pijji_idx.push_back(ji_ij);

            Pijji_beta_output.push_back("P(including 2 f^ij_mA F^m_P f^PA_ji) (beta alpha case)");
            Pijji_beta_fx12_ints.push_back(b_i1i2mF2A1_key);
            // Pijji_idx.push_back(ij_ji);
            Pijji_beta_fx12_ints.push_back(b_i2i1mF2A1_key);
            // Pijji_idx.push_back(ji_ij);

            Pijji_beta_output.push_back("P(including 2 f^ij_Ae F^A_p f^pe_ji) (beta alpha case)");
            Pijji_beta_fx12_ints.push_back(b_i1i2AF2e1_ints);
            // Pijji_idx.push_back(ij_ji);
            Pijji_beta_fx12_ints.push_back(b_i2i1AF2e1_ints);
            // Pijji_idx.push_back(ji_ij);

            Pijji_beta_output.push_back("P(including f^ij_mA F^m_n f^nA_ji) (beta alpha case)");
            Pijji_beta_fx12_ints.push_back(b_i1i2mf2A1_ints);
            // Pijji_idx.push_back(ij_ji);
            Pijji_beta_fx12_ints.push_back(b_i2i1mf2A1_ints);
            // Pijji_idx.push_back(ji_ij);

            accumulate_P_YxF(Pijji_beta_output,
                             Pijji_idx, P_prefactors,
                             f12_idx, f12_idx,
                             Pijij_beta_f12_ints, Pijji_beta_fx12_ints,
                             Pij_ji_beta);

            // Qij_ji += (f12f12)^ij_Pi (f+k)^P_j + (f12f12)^ij_jP (f+k)^P_i
            //         = (f12f12)^ij _j_hJ i + (f12f12)^ji _i_hJ j
            // Qij_ji (beta alpha case) = Qji_ij (alpha beta case)
            ExEnv::out0() << indent << spinletters
                                    << " Q^ij_ji (beta alpha case) = Qji_ij (alpha beta case) "
                                    << endl << endl;

#if 0
            // test
            ExEnv::out0() << indent << "Q^ij_ji (beta alpha case): " << endl;
            double* Qij_ji_beta = new double[nocc12];
            std::fill(Qij_ji_beta, Qij_ji_beta + nocc12, 0.0);

            Ref<DistArray4> q_i2i1hi1i2_ints;
            activate_ints(occ2_act->id(), occ1_act->id(),
                          r12eval_->hj_i_P(spin1)->id(), occ2_act->id(),
                          descr_f12f12_key, moints4_rtime,
                          q_i2i1hi1i2_ints);
            compute_Y(ij_ji, 1.0,
                      f12f12_idx,
                      q_i2i1hi1i2_ints, Qij_ji_beta);
            ExEnv::out0() << endl;

            Ref<DistArray4> q_i1i2hi2i1_ints;
            activate_ints(occ1_act->id(), occ2_act->id(),
                          r12eval_->hj_i_P(spin2)->id(), occ1_act->id(),
                          descr_f12f12_key, moints4_rtime,
                          q_i1i2hi2i1_ints);
            compute_Y(ji_ij, 1.0,
                      f12f12_idx,
                      q_i1i2hi2i1_ints, Qij_ji_beta);
            ExEnv::out0() << endl;
#endif

            // B^ij_ji = B^ij_ji - P^ij_ji + Q^ij_ji
            ExEnv::out0() << indent << spinletters << " Bij_ji (beta alpha case): " << endl;
      #if 1
            for(int i1=0; i1<nocc2_act; ++i1) {
              for(int i2=0; i2<nocc1_act; ++i2) {

                const int i1i2 = i1 * nocc1_act + i2;
                const int i2i1 = i2 * nocc2_act + i1;
                Bij_ji_beta[i1i2] += - Pij_ji_beta[i1i2] + Qij_ji[i2i1];
                ExEnv::out0() << indent << scprintf("%12.10f",Bij_ji_beta[i1i2])<< "  ";
              }
              ExEnv::out0() << endl;
            }
            ExEnv::out0() << endl;
      #endif
            delete[] Pij_ji_beta;
            Pij_ji_beta = NULL;
      }

      //
      // when spin1 != spin2: B^ji_ij = B^ij_ji, B^ji_ji = B^ij_ij
      //

      // Print the antisymmetrized Vij_ij, Xij_ij and Bij_ij
      // for alpha_alpha and beta_beta case
      if (spin1 == spin2){

          ExEnv::out0() << indent << spinletters << " antisymmetrized V"<< endl;
          for(int i1=0; i1<nocc1_act; ++i1) {
              for(int i2=i1+1; i2<nocc2_act; ++i2) {
                  int i1i2 = i1 * nocc2_act + i2;
                  ExEnv::out0() << indent << Vij_ij[i1i2] - Vij_ji[i1i2] << "  ";
              }
              ExEnv::out0() << endl;
          }
          ExEnv::out0() << endl;

          ExEnv::out0() << indent << spinletters << " antisymmetrized X"<< endl;
          for(int i1=0; i1<nocc1_act; ++i1) {
              for(int i2=i1+1; i2<nocc2_act; ++i2) {
                  int i1i2 = i1 * nocc2_act + i2;
                  ExEnv::out0() << indent << Xij_ij[i1i2] - Xij_ji[i1i2] << "  ";
              }
              ExEnv::out0() << endl;
          }
          ExEnv::out0() << endl;

          ExEnv::out0() << indent << spinletters << " antisymmetrized Q"<< endl;
          for(int i1=0; i1<nocc1_act; ++i1) {
              for(int i2=i1+1; i2<nocc2_act; ++i2) {
                  int i1i2 = i1 * nocc2_act + i2;
                  ExEnv::out0() << indent << Qij_ij[i1i2] - Qij_ji[i1i2] << "  ";
              }
              ExEnv::out0() << endl;
          }
          ExEnv::out0() << endl;

          ExEnv::out0() << indent << spinletters << " antisymmetrized B"<< endl;
          for(int i1=0; i1<nocc1_act; ++i1) {
              for(int i2=i1+1; i2<nocc2_act; ++i2) {
                  int i1i2 = i1 * nocc2_act + i2;
                  ExEnv::out0() << indent << Bij_ij[i1i2] - Bij_ji[i1i2] << "  ";
              }
              ExEnv::out0() << endl;
          }
          ExEnv::out0() << endl;
      }

      //
      // for each ij pair compute its contribution to the Hylleraas functional for second-order R12 energy
      const RefDiagSCMatrix evals1 = occ1_act->evals();
      const RefDiagSCMatrix evals2 = occ2_act->evals();
      ExEnv::out0() << indent << "evals1: " << endl;
      for(int i1=0; i1<nocc1_act; ++i1) {
        ExEnv::out0() << indent << evals1(i1) << endl;
      }
      ExEnv::out0() << indent << "evals2: " << endl;
      for(int i2=0; i2<nocc2_act; ++i2) {
        ExEnv::out0() << indent << evals2(i2) << endl;
      }

      // compute coefficients
      double C_0;  // singlet
      double C_1;  // triplet
      if (R12Technology::STG(corrfactor->geminaldescriptor())) {
        const double gamma = R12Technology::single_slater_exponent(corrfactor->geminaldescriptor());
        C_0=-1.0/(2.0*gamma);
        C_1=-1.0/(4.0*gamma);
      }
      else if (R12Technology::R12(corrfactor->geminaldescriptor())) {
        C_0=1.0/2.0;
        C_1=1.0/4.0;
      } else {
        throw ProgrammingError("firstorder_cusp_coefficients -- geminal coefficients can be only determined for STG or R12 geminals",__FILE__,__LINE__);
      }

      // Compute the f12 correction pair energy
      ExEnv::out0() << indent << "i " << " j  " << spinletters
                              << " Hij" << endl;
      if (spin1 == spin2) {
          // Alpha_alpha or beta_beta case
        double f12_energy = 0;

        for(int i1=0; i1<nocc1_act; ++i1) {
          for(int i2=i1+1; i2<nocc2_act; ++i2) {

            double Hij_pair_energy;
            int ij = i1 * nocc2_act + i2;

            Hij_pair_energy =  2.0 * C_1 * (Vij_ij[ij] - Vij_ji[ij])
                                   + C_1*C_1 * (Bij_ij[ij] - Bij_ji[ij] - (evals1(i1) + evals2(i2)) * (Xij_ij[ij] - Xij_ji[ij]));
            ExEnv::out0() << indent << i1 << "  " << i2 << "  "
                                    << scprintf("%12.10f", Hij_pair_energy) << endl;

            f12_energy += Hij_pair_energy;
          }
        }

        if (spin1 == 0 && spin2 == 0) {
          f12_energy_aa = f12_energy;
          ExEnv::out0() << indent << "alpha-alpha f12 correction energy is "
                                  << scprintf("%12.10f", f12_energy_aa) << endl;
        } else {
          f12_energy_bb = f12_energy;
          ExEnv::out0() << indent << "beta-beta f12 correction energy is "
                                  << scprintf("%12.10f", f12_energy_bb) << endl;
        }

      }  else if (num_unique_spincases2 == 2) {
          // Alpha_beta case for close shell
          for(int i1=0; i1<nocc1_act; ++i1) {
            for(int i2=0; i2<nocc2_act; ++i2) {
              double Hij_pair_energy;
              int ij = i1 * nocc2_act + i2;

              Hij_pair_energy =   2.0 * ( 0.5*(C_0+C_1) * Vij_ij[ij] + 0.5*(C_0-C_1) * Vij_ji[ij])
                                  + pow(0.5*(C_0+C_1), 2) * (Bij_ij[ij]- (evals1(i1) + evals2(i2)) * Xij_ij[ij])
                                  + 0.25*(C_0*C_0 - C_1*C_1)*2.0*(Bij_ji[ij] - (evals1(i1) + evals2(i2)) * Xij_ji[ij])
                                  + pow(0.5*(C_0-C_1), 2) * (Bij_ij[ij] - (evals1(i1) + evals2(i2)) * Xij_ij[ij]);
              ExEnv::out0() << indent << i1 << "  " << i2 << "  "
                  << scprintf("%12.10f", Hij_pair_energy) << endl;
              f12_energy_ab += Hij_pair_energy;
            }
          }
          ExEnv::out0() << indent << "unique spincases: "<< num_unique_spincases2 << "  alpha beta f12 correction energy is "
              << scprintf("%12.10f", f12_energy_ab) << endl;
      } else {
          // Alpha_beta case for open shell
          for(int i1=0; i1<nocc1_act; ++i1) {
            for(int i2=0; i2<nocc2_act; ++i2) {
              double Hij_pair_energy;
              int ij = i1 * nocc2_act + i2;
              int ji = i2 * nocc1_act + i1;

              Hij_pair_energy =   2.0 * ( 0.5*(C_0+C_1) * Vij_ij[ij] + 0.5*(C_0-C_1) * Vij_ji[ij])
                                + pow(0.5*(C_0+C_1), 2) * ((Bij_ij[ij] +Bij_ij_beta[ji])*0.5 - (evals1(i1) + evals2(i2)) * Xij_ij[ij])
                                + 0.25*(C_0*C_0 - C_1*C_1)*2.0*((Bij_ji[ij] +Bij_ji_beta[ji])*0.5 - (evals1(i1) + evals2(i2)) * Xij_ji[ij])
                                + pow(0.5*(C_0-C_1), 2) * ((Bij_ij[ij] +Bij_ij_beta[ji])*0.5 - (evals1(i1) + evals2(i2)) * Xij_ij[ij]);
              ExEnv::out0() << indent << i1 << "  " << i2 << "  "
                                      << scprintf("%12.10f", Hij_pair_energy) << endl;
              f12_energy_ab += Hij_pair_energy;
            }
          }
          ExEnv::out0() << indent << "unique spincases: "<< num_unique_spincases2 << "  alpha beta f12 correction energy is "
                                  << scprintf("%12.10f", f12_energy_ab) << endl;
          // deallocate memory for B_beta (it will be used once)
          delete[] Bij_ij_beta;
          delete[] Bij_ji_beta;
      }

  } // end of spincase loop

  // deallocate memory
  delete[] Vij_ij;
  delete[] Vij_ji;
  delete[] Xij_ij;
  delete[] Xij_ji;
  delete[] Bij_ij;
  delete[] Bij_ji;
  delete[] Pij_ij;
  delete[] Pij_ji;
  delete[] Qij_ij;
  delete[] Qij_ji;

  // Set beta-beta energies to alpha-alpha for closed-shell
  if (!r12world->ref()->spin_polarized()) {
    C_[BetaBeta] = C_[AlphaAlpha];
    emp2f12_[BetaBeta] = emp2f12_[AlphaAlpha];
    ef12_[BetaBeta] = ef12_[AlphaAlpha];
    f12_energy_bb = f12_energy_aa;
  }

  double f12_energy = f12_energy_aa + f12_energy_bb + f12_energy_ab;
  ExEnv::out0() << indent << "f12 correction energy is "
                          << scprintf("%12.10f", f12_energy) << endl;

  evaluated_ = true;

  return;
}

