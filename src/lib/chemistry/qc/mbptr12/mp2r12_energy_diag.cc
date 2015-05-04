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

#include <cassert>
#include <cmath>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#ifdef MPQC_NEW_FEATURES
#include <chemistry/qc/mbptr12/sr_r12intermediates.h>
#endif // MPQC_NEW_FEATURES
#include <util/misc/print.h>
#include <math/scmat/blas.h>

using namespace std;
using namespace sc;

/**********************************
 * class MP2R12Energy_Diag *
 **********************************/

namespace {
  void _print(SpinCase2 spin,
             const Ref<DistArray4>& mat,
             const char* label) {
    if (mat->msg()->me() == 0) {
      const size_t nij = (spin != AlphaBeta && mat->ni() == mat->nj()) ? mat->ni() * (mat->ni()-1) / 2 : mat->ni() * mat->nj();
      const size_t nxy = (spin != AlphaBeta && mat->nx() == mat->ny()) ? mat->nx() * (mat->nx()-1) / 2 : mat->nx() * mat->ny();
      RefSCMatrix scmat = SCMatrixKit::default_matrixkit()->matrix(new SCDimension(nij), new SCDimension(nxy));
      scmat << mat;
      scmat.print(label);
    }
  }

  // Swap electron 1 and 2 in the integral
  void swap_e12(const double* ints_ijba,
                  const int nb, const int na,
                  double* const ints_jiab) {

    const double* iter_ints_ijba = ints_ijba;
    const int offset = nb;

    for (int b=0; b<nb; ++b) {
      double* iter_ints_jiab = ints_jiab + b;

      for (int a=0; a<na; ++a) {
        const int ab = a*nb+b;
        *iter_ints_jiab = *iter_ints_ijba;
        ++iter_ints_ijba;
        iter_ints_jiab += offset;
      }
    }
  }

  void print_intermediate(SpinCase2 spin, const string& label,
                          const double* const array,
                          const int no1, const int no2) {

    string spinletters = to_string(spin);
    ExEnv::out0() << indent <<  spinletters << "  " << label << endl;
    if (spin == AlphaBeta) {
        const double* iter = array;
        for (int i=0 ; i < no1; ++i) {
          for (int j = 0; j < no2; ++j, ++iter) {
              ExEnv::out0() << indent << scprintf("%12.10f", *iter) << "  ";
          }
          ExEnv::out0() << endl;
        }
    } else if (no1 != 1) {
        const int unique_ij = no1*(no1-1)/2;
        double* const array_unique_ij = new double[unique_ij];
        for (int i=0 ; i < no1; ++i) {
          for (int j = i+1; j < no2; ++j) {
              const int ij = j*(j-1)/2 + i;
              const int i1i2 = i*no2 +j;
              array_unique_ij[ij] = array[i1i2];

          }
        }

        const double* iter = array_unique_ij;
        const int i_end = (no1-1)/2;
        for (int i=0; i<i_end; ++i) {
          for (int j=0 ; j<no2; ++j, ++iter) {
            ExEnv::out0() << indent << scprintf("%12.10f", *iter);
          }
          ExEnv::out0() << endl;
        }
        delete[] array_unique_ij;
    } else {
        // There is only one element
        ExEnv::out0() << indent << scprintf("%12.10f", array[0]) << "  ";
    }
    ExEnv::out0() << endl;
  }

  void print_antisym_intermediate(SpinCase2 spin, const string& label,
                                  const double* const array_ij, const double* const array_ji,
                                  const int no1, const int no2) {

    string spinletters = to_string(spin);
    ExEnv::out0() << indent <<  spinletters << "  " << label << endl;
    if (spin == AlphaBeta) {
        ExEnv::out0() << indent << "error: this is not only needed for alpha-beta cases" << endl;
    } else {
        const int unique_ij = no1*(no1-1)/2;
        double* const array1_unique_ij = new double[unique_ij];
        double* const array2_unique_ij = new double[unique_ij];
        for (int i=0 ; i < no1; ++i) {
          for (int j = i+1; j < no2; ++j) {
              const int ij = j*(j-1)/2 + i;
              const int i1i2 = i*no2 +j;
              array1_unique_ij[ij] = array_ij[i1i2];
              array2_unique_ij[ij] = array_ji[i1i2];

          }
        }

        const double* iter1 = array1_unique_ij;
        const double* iter2 = array2_unique_ij;
        const int i_end = (no1-1)/2;
        for (int i=0; i<i_end; ++i) {
          for (int j=0 ; j<no2; ++j, ++iter1, ++iter2) {
            ExEnv::out0() << indent << scprintf("%12.10f", *iter1-*iter2);
          }
          ExEnv::out0() << endl;
        }
        delete[] array1_unique_ij;
        delete[] array2_unique_ij;
    }
    ExEnv::out0() << endl;
  }
}

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
  i1i2i1i2_ints = i1i2i1i2_tform->ints_distarray4();
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
  int nocc2_act = i1i2i1i2_ints->nj();
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
      if (debug_ >= DefaultPrintThresholds::mostN2)
        ExEnv::out0() << indent << scprintf("%12.10f", array_i1i2i1i2[*array_idx]) << " ";

      ++idx_ij;
      idx_ji += nocc1_act;

      i1i2i1i2_ints->release_pair_block(*blk_i1, *blk_i2, oper_idx);
    }
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl;
  }
}

// Compute (YxF)^ij_ij or (YxF)^ij_ji
// Y can be g or f
void MP2R12Energy_Diag::compute_YxF(const int b1b2_k1k2, const double prefactor,
                                    const unsigned int oper1_idx, const unsigned int oper2_idx,
                                    const Ref<DistArray4>& i1i2x1x2_ints, const Ref<DistArray4>& j1j2x1x2_ints,
                                    double* array_i1i2i1i2)
{
  const blasint norbs12 = i1i2x1x2_ints->nx() * i1i2x1x2_ints->ny();

  int nocc1_act = i1i2x1x2_ints->ni();
  int nocc2_act = i1i2x1x2_ints->nj();
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

  // the index for array which stores V^ij_ij, V^ij_ji, V^ji_ij, and V^ji_ji is always ij
  // ij = i * nocc2_act + j which increase one for each loop
  int array_idx = 0;

  for( ; i<nocc1_act; ++i) {
    for(j=0; j<nocc2_act; ++j) {
      const blasint one = 1;

      const double* blk1 = i1i2x1x2_ints->retrieve_pair_block(*blk1_i1, *blk1_i2, oper1_idx);
      const double* blk2 = j1j2x1x2_ints->retrieve_pair_block(*blk2_i1, *blk2_i2, oper2_idx);
      const double i1i2i1i2 = prefactor * F77_DDOT(&norbs12, blk1, &one, blk2, &one);

      array_i1i2i1i2[array_idx] += i1i2i1i2 ;
      if (debug_ >= DefaultPrintThresholds::mostN2)
        ExEnv::out0() << indent << scprintf("%12.10f",array_i1i2i1i2[array_idx]) << " ";

      i1i2x1x2_ints->release_pair_block(*blk1_i1, *blk1_i2, oper1_idx);
      j1j2x1x2_ints->release_pair_block(*blk2_i1, *blk2_i2, oper2_idx);

      ++array_idx;
    }
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl;
  }
}

// Compute (FxT)^ij_ij, (FxT)^ij_ji, (FxT)^ji_ij, (FxT)^ji_ji for V coupling
// (FxT)^ij_ij = R^ij_ab T^ab_ij
void MP2R12Energy_Diag::compute_FxT(const int b1b2_k1k2, const unsigned int f12_idx,
                                    Ref<DistArray4>& iiPFa_ints, const double* Ti1i2ab,
                                    double* V_coupling_array)
{
  // iiPFa_int dimension: nocc1_act*nocc2_act*nvir1_act*nvir1_act
  const blasint nvir12_act = iiPFa_ints->nx() * iiPFa_ints->ny();

  // ij case
  int nocc1_act = iiPFa_ints->ni();
  int nocc2_act = iiPFa_ints->nj();

  // get the right pair block from iiPFa_ints
  int i=0;
  int j=0;
  // Select block indices
  int* F_blk_i1 = &i;
  int* F_blk_i2 = &j;

  if (b1b2_k1k2 == ji_ij || b1b2_k1k2 == ji_ji) {
    nocc1_act = iiPFa_ints->nj();
    nocc2_act = iiPFa_ints->ni();
    F_blk_i1 = &j;
    F_blk_i2 = &i;
  }

  // Index V_coupling_array by ij: increase one each loop
  int idx_ij = 0;

  // Get the right start iterator for Tij or Tji
  int iter_Tij_start = 0;
  int iter_Tji_start = 0;

  // Get the right offset for Tij or Tji in the inner loop
  int offset_Tij = nvir12_act;
  // Tji only exist for alpha-alpha, beta-beta, close-shell
  // so it does not matter * nocc1_act or nocc2_act
  int offset_Tji = nvir12_act * nocc1_act;

  int* iter_T_start = &iter_Tij_start;
  int* offset_T = &offset_Tij;
  if (b1b2_k1k2 == ij_ji || b1b2_k1k2 == ji_ji) {
    offset_T = &offset_Tji;
    iter_T_start = &iter_Tji_start;
  }


  for( ; i<nocc1_act; ++i) {
    iter_Tij_start = i * nvir12_act * nocc2_act;
    iter_Tji_start = i * nvir12_act;
    const double* iter_T = Ti1i2ab + *iter_T_start;

    for(j=0; j<nocc2_act; ++j, ++idx_ij, iter_T += *offset_T) {
      const blasint one = 1;
      const double* F_blk = iiPFa_ints->retrieve_pair_block(*F_blk_i1, *F_blk_i2, f12_idx);

      // Get the right block of Tij or Tji
      const double i1i2i1i2 = F77_DDOT(&nvir12_act, F_blk, &one, iter_T, &one);

      V_coupling_array[idx_ij] += i1i2i1i2 ;
      if (debug_ >= DefaultPrintThresholds::mostN2)
        ExEnv::out0() << indent << scprintf("%12.10f",V_coupling_array[idx_ij]) << "  ";

      iiPFa_ints->release_pair_block(*F_blk_i1, *F_blk_i2, f12_idx);
    }
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl;
  }
}

void MP2R12Energy_Diag::compute_VX(const int b1b2_k1k2, std::vector<std::string>& VX_output,
                                   const unsigned int oper12_idx, Ref<DistArray4>& Y12_ints,
                                   const unsigned int oper1_idx, const unsigned int oper2_idx,
                                   const std::vector<Ref<DistArray4> >& Y_ints,
                                   const std::vector<Ref<DistArray4> >& F_ints,
                                   double* VX_array)
{
  if (debug_ >= DefaultPrintThresholds::mostN2)
    ExEnv::out0() << indent << "(diag) contribution" << endl;

  compute_Y(b1b2_k1k2, 1.0,
            oper12_idx,
            Y12_ints, VX_array);

  const int nterms = Y_ints.size();
  for (int i=0; i != nterms; ++i) {
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << indent << VX_output[i] << endl;

    compute_YxF(b1b2_k1k2, -1.0,
                oper1_idx, oper2_idx,
                Y_ints[i], F_ints[i],
                VX_array);
  }
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
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << *output << endl;

    // P += f^ij_xy F^x_z f^zy_ij
    //    = ()^ij_ij
    compute_YxF(b1b2_k1k2[i], *prefactor,
                oper1_idx, oper2_idx,
                Y_ints[i], F_ints[i],
                P);

    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl;

    // P += f^ji_xy F^x_z f^zy_ji
    //    = ()^ji_ji
    ++i;
    compute_YxF(b1b2_k1k2[i], *prefactor,
                oper1_idx, oper2_idx,
                Y_ints[i], F_ints[i],
                P);
  }
}

// Compute Via * T1^j_a or Vaj * T1^i_a
void MP2R12Energy_Diag::contract_VT1(const Ref<DistArray4>& V,
                                     const int b1b2_k1k2,  const bool swap_e12_V,
                                     const double* const T1_array,
                                     const int nv, const bool T1_ia,
                                     double* const VT1) {
   int no1 = V->ni();
   int no2 = V->nj();
   if (b1b2_k1k2 == ji_ij) {
     no1 = V->nj();
     no2 = V->ni();
   }
   const int norbs1 = V->nx();
   const int norbs2 = V->ny();
   const int norbs12 = norbs1 * norbs2;

   // Get the pair block ij or ji based on k1k2
   int i=0;
   int j=0;
   int* blk_i1 = &i;
   int* blk_i2 = &j;
   if (b1b2_k1k2 == ji_ij) {
     blk_i1 = &j;
     blk_i2 = &i;
   }
   // Get the offset for V and T1 in order to obtain Va and T1a block
   int* offset_V = &i;
   int* offset_T1 = &j;      // T^j_a (beta)
   if (T1_ia) {
     offset_V = &j;
     offset_T1 = &i;   // T^i_a (alpha)
   }

   const blasint one = 1;
   // Allocate memory for V_blk which contains the swapped electron 1 and 2 integrals
   double* Vma_blk = NULL;
   if (swap_e12_V) {
     Vma_blk = new double[norbs12];
     fill_n(Vma_blk, norbs12, 0.0);
   }

   double* iter_VT1 = VT1;
   for (; i<no1; ++i) {
     for (j=0; j<no2; ++j) {
       const double* const V_blk = V->retrieve_pair_block(*blk_i1, *blk_i2, 0);

       const double* iter_V = NULL;
       if (swap_e12_V) {
           swap_e12(V_blk,norbs1,norbs2,Vma_blk);
           iter_V = Vma_blk + (*offset_V) * nv;
       } else
         iter_V = V_blk + (*offset_V) * nv;

       const double* iter_T1 = T1_array + (*offset_T1) * nv;

       const blasint nvf = nv;
       const double vt1_ij = F77_DDOT(&nvf, iter_V, &one, iter_T1, &one);
       const int ij = i*no2 + j;
       *iter_VT1 += vt1_ij;
       ++iter_VT1;

       V->release_pair_block(*blk_i1, *blk_i2, 0);
     }
   }
   delete[] Vma_blk;
 }


void MP2R12Energy_Diag::compute_ef12() {

//  if (this->r12eval()->compute_1rdm()) {
//    // tests for open-shell systems using MP2-F12 method
//    if (r12intermediates_->Onerdm_relax_computed()) {
//      const int nspincases1 = r12eval()->nspincases1();
//      for(int s = 0; s < nspincases1; ++s) {
//        const SpinCase1 spin = static_cast<SpinCase1>(s);
//        RefSCMatrix Xf12_alpha = Onerdm_X_F12(spin, r12eval_, debug_);
//
//        // test: print the Z-vector from CABS Singles contribution
//        RefSCMatrix X_cabs = Onerdm_X_CABS_Singles(spin, r12eval_, r12intermediates_, debug_);
//        X_cabs.print(prepend_spincase(spin,"CABS_Singles Z-vector X:").c_str());
//      }
//    }
//    ExEnv::out0() << std:: endl << "Onerdm_cc_computed() " << r12intermediates_->Onerdm_cc_computed() << std::endl;
//    compute_density_diag();
//  }

  // switch to new implementation that should work correctly for alpha-beta contributions in open-shell molecules
  return this->compute_ef12_10132011();

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
  if (debug_ >= DefaultPrintThresholds::N2)
    ExEnv::out0() <<endl << indent << "num_unique_spincases2 = " << num_unique_spincases2 << endl;

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
  const SpinCase2 spincase_AlphAlpha = static_cast<SpinCase2> (1);
  if (r12eval()->dim_oo(spincase_AlphAlpha).n() != 0) {

      const SpinCase1 spin1_alpha = case1(spincase_AlphAlpha);
      const Ref<OrbitalSpace>& occ_act_alpha = r12eval()->occ_act(spin1_alpha);
      const Ref<OrbitalSpace>& vir_act_alpha = r12eval()->vir_act(spin1_alpha);
      nalpha = occ_act_alpha->rank();
      nvir_act_alpha = vir_act_alpha->rank();
  }
  // Beta-beta case
  const SpinCase2 spincase_BetaBeta = static_cast<SpinCase2> (2);
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
      const SpinCase2 spincase_AlphaBeta = static_cast<SpinCase2> (0);
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
  Ref<DistArray4> T2[NSpinCases2];   // need to be outside for coupling term
  if (r12intermediates_->T1_cc_computed() &&
      r12intermediates_->T2_cc_computed() ) {
    // Obtain T1 amplitudes
    const int nspincases1 = r12eval()->nspincases1();
    for(int s=0; s<nspincases1; ++s) {
      const SpinCase1 spin = static_cast<SpinCase1>(s);
      T1[spin] = r12intermediates_->get_T1_cc(spin);
    }
    if (nspincases1 == 1) {
      T1[Beta] = T1[Alpha];
    }

    // Obtain T2 amplitudes
    for(int s=0; s<num_unique_spincases2; ++s) {
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      if (r12eval()->dim_oo(spincase2).n() == 0)
        continue;
      T2[spincase2] = r12intermediates_->get_T2_cc(spincase2);
    } // end of T2
  } // end of CC amplitudes

  //
  // geminal coefficients
  const double C_0=1.0/2.0;
  const double C_1=1.0/4.0;

  //
  // Compute the MP2-R12 part
  //
  // Assign pointers to each matrix
  double* Vij_ij = new double[nijij];
  double* Vij_ji = new double[nijij];

  // Assign pointers to V coupling matrixes
  double* Vij_ij_coupling = NULL;
  double* Vij_ji_coupling = NULL;
  //double* Vji_ij_coupling = NULL;
  //double* Vji_ji_coupling = NULL;
  double* Tij_ab = NULL;
  if (this->r12eval()->coupling() == true
      || this->r12eval()->ebc() == false) {
    Vij_ij_coupling = new double[nijij];
    Vij_ji_coupling = new double[nijij];
    Tij_ab = new double[nijab];
  }

  double* Xij_ij = new double[nijij];
  double* Xij_ji = new double[nijij];

  double* Bij_ij = new double[nijij];
  double* Pij_ij = new double[nijij];
  double* Qij_ij = new double[nijij];
  double* Bij_ji = new double[nijij];
  double* Pij_ji = new double[nijij];
  double* Qij_ji = new double[nijij];

  for (int spin = 0; spin < num_unique_spincases2; spin++) {

    const SpinCase2 spincase = static_cast<SpinCase2> (spin);
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
    const int nvir12_act = nvir1_act * nvir2_act;

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
    const Ref<OrbitalSpace>& fvir1_act = r12eval()->F_a_A(spin1);

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
    const Ref<OrbitalSpace>& fvir2_act = r12eval()->F_a_A(spin2);

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
    if (occ1_act != occ2_act){
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
    activate_ints(occ1_act->id(), occ2_act->id(),
                  occ1_act->id(), occ2_act->id(),
                  descr_f12_key, moints4_rtime,
                  i1i2i1i2_ints);

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
    // TODO: for RHF can simply scale previous contribution by 2 and "symmetrize" at the end V^ij_ij = 0.5 * (V^ij_ij + V^ji_ji)
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

    if (debug_ >= DefaultPrintThresholds::N2 && spin1 == spin2)
      print_antisym_intermediate(spincase,"V^ij_ij",Vij_ij, Vij_ji, nocc1_act,nocc2_act);

    if (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2) {
      // Alpha-beta case
      print_intermediate(spincase,"V^ij_ij",Vij_ij,nocc1_act,nocc2_act);
      print_intermediate(spincase,"V^ij_ji",Vij_ji,nocc1_act,nocc2_act);
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
    if (vir1_act != vir2_act){
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
#if 1
      if (debug_ >= DefaultPrintThresholds::mostN2) {
        ExEnv::out0() << endl << indent << "evals_i1: " << endl;
        for(int i1=0; i1<nocc1_act; ++i1) {
          ExEnv::out0() << indent << evals_i1(i1) << endl;
        }

        ExEnv::out0() << endl << indent << "evals_i2: " << endl;
        for(int i2=0; i2<nocc2_act; ++i2) {
          ExEnv::out0() << indent << evals_i2(i2) << endl;
        }

        ExEnv::out0() << endl << indent << "evals_a1: " << endl;
        for(int a1=0; a1<nvir1_act; ++a1) {
          ExEnv::out0() << indent << evals_a1(a1) << endl;
        }

        ExEnv::out0() << endl << indent << "evals_a2: " << endl;
        for(int a2=0; a2<nvir2_act; ++a2) {
          ExEnv::out0() << indent << evals_a2(a2) << endl;
        }
      }
#endif

      // Initialize all the integrals needed
      // Vij_ij_coupling: R^ij_aa' f^a'_b T^ab_ij
      Ref<DistArray4> i1i2a1AF2_ints = NULL;
      activate_ints(occ1_act->id(), occ2_act->id(),
                    vir1_act->id(), fvir2_act->id(),
                    descr_f12_key, moints4_rtime,
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
      if (spin1 != spin2) {
      activate_ints(occ1_act->id(), occ2_act->id(),
                    fvir1_act->id(), vir2_act->id(),
                    descr_f12_key, moints4_rtime,
                    i1i2AF1a2_ints);
      }

      Ref<DistArray4> i2i1AF1a2_ints = NULL;
      Ref<DistArray4> i2i1a1AF2_ints = NULL;
      if (num_unique_spincases2 == 3 && spin1 != spin2){

          // i:alpha  j:beta   a:alpha b:beta
          // Vij_ji_coupling: R^ji_a'b f^a'_a T^ab_ij
          // a':alpha
          activate_ints(occ2_act->id(), occ1_act->id(),
                        fvir1_act->id(), vir2_act->id(),
                        descr_f12_key, moints4_rtime,
                        i2i1AF1a2_ints);

          // Vji_ij_coupling: R^ji_aa' f^a'_b T^ab_ij
          // a':beta
          activate_ints(occ2_act->id(), occ1_act->id(),
                        vir1_act->id(), fvir2_act->id(),
                        descr_f12_key, moints4_rtime,
                        i2i1a1AF2_ints);
      } else {
          i2i1AF1a2_ints = i1i2AF1a2_ints;
          i2i1a1AF2_ints = i1i2a1AF2_ints;
      }

      if(r12intermediates_->T2_cc_computed()) {
        if (debug_ >= DefaultPrintThresholds::N2)
          ExEnv::out0() << endl << indent << "Coupled-cluster V coupling:" << endl;
        T2[spin]->activate();
        if (spin1 != spin2) {
          compute_YxF(ij_ij, 1.0,
                      f12_idx, 0,
                      i1i2a1AF2_ints, T2[spin],
                      Vij_ij_coupling);
          compute_YxF(ij_ij, 1.0,
                      f12_idx, 0,
                      i1i2AF1a2_ints, T2[spin],
                      Vij_ij_coupling);

          compute_YxF(ji_ij, 1.0,
                      f12_idx, 0,
                      i2i1AF1a2_ints, T2[spin],
                      Vij_ji_coupling);
          compute_YxF(ji_ij, 1.0,
                      f12_idx, 0,
                      i2i1a1AF2_ints, T2[spin],
                      Vij_ji_coupling);

        } else {
            compute_YxF(ij_ij, 1.0,
                        f12_idx, 0,
                        i1i2a1AF2_ints, T2[spin],
                        Vij_ij_coupling);
            compute_YxF(ji_ij, 1.0,
                        f12_idx, 0,
                        i1i2a1AF2_ints, T2[spin],
                        Vij_ji_coupling);
        }
        T2[spin]->deactivate();
      } else {
          // Start computing MP2 V coupling
          if (debug_ >= DefaultPrintThresholds::N2)
            ExEnv::out0() << endl << indent << spinletters << "  MP2 V coupling:" << endl;
          // Compute T^ij_ab (not antisymmetrized)
          // g^ij_ab
          Ref<DistArray4> i1i2a1a2_ints;
          activate_ints(occ1_act->id(), occ2_act->id(),
                        vir1_act->id(), vir2_act->id(),
                        descr_f12_key, moints4_rtime,
                        i1i2a1a2_ints);

          //   T^i(spin1)j(spin2)_a(spin1)b(spin2) 4-dimention matrix
          // = g^ij_ab / (e_i + e_j - e_a -e_b)
          double* iter_Tij_ab = Tij_ab;
          for (int i1=0; i1<nocc1_act; ++i1){
            for (int i2=0; i2<nocc2_act; ++i2){
              const double* gij_ab = i1i2a1a2_ints->retrieve_pair_block(i1, i2, eri_idx);
              double mp2_pair_energy = 0;

              for (int a1=0; a1<nvir1_act; ++a1){
                for (int a2=0; a2<nvir2_act; ++a2, ++iter_Tij_ab, ++gij_ab){
                  *iter_Tij_ab = (*gij_ab) /(evals_i1(i1)+evals_i2(i2)-evals_a1(a1)-evals_a2(a2));
                }
              }
              i1i2a1a2_ints->release_pair_block(i1, i2, eri_idx);
            }
          }
          i1i2a1a2_ints->deactivate();

          // Alpha-alpha, beta-beta and alpha-beta for close shell:
          // Tji_ab matrix <=> Tij_ab matrix

          if (debug_ >= DefaultPrintThresholds::mostN2)
            ExEnv::out0() << endl << indent << spinletters << " Vij_ij_coupling: R^ij_aa' f^a'_b T^ab_ij" << endl;
          compute_FxT(ij_ij, f12_idx,
                      i1i2a1AF2_ints, Tij_ab,
                      Vij_ij_coupling);

          if (spin1 != spin2){

              if (debug_ >= DefaultPrintThresholds::mostN2)
                ExEnv::out0() << indent << spinletters << " Vij_ij_coupling(a':alpha): + R^ij_a'b f^a'_a T^ab_ij" << endl;
              compute_FxT(ij_ij, f12_idx,
                          i1i2AF1a2_ints, Tij_ab,
                          Vij_ij_coupling);

              if (debug_ >= DefaultPrintThresholds::mostN2)
                ExEnv::out0() << indent << spinletters << " Vji_ij_coupling: R^ji_a'b f^a'_a T^ab_ij" << endl;
              compute_FxT(ji_ij, f12_idx,
                          i2i1AF1a2_ints, Tij_ab,
                          Vij_ji_coupling);

              if (debug_ >= DefaultPrintThresholds::mostN2)
                ExEnv::out0() << indent << spinletters << " Vji_ij_coupling(a':beta): + R^ji_aa' f^a'_b T^ab_ij" << endl;
              compute_FxT(ji_ij, f12_idx,
                          i2i1a1AF2_ints, Tij_ab,
                          Vij_ji_coupling);
          } else {
              // + R^ji_aa' f^a'_b T^ab_ji
              if (debug_ >= DefaultPrintThresholds::mostN2)
                ExEnv::out0() << indent << spinletters << " Vji_ji_coupling: + R^ji_aa' f^a'_b T^ab_ji" << endl;
              compute_FxT(ji_ji, f12_idx,
                          i1i2a1AF2_ints, Tij_ab,
                          Vij_ij_coupling);
              // - R^ij_aa' f^a'_b T^ab_ji
              if (debug_ >= DefaultPrintThresholds::mostN2)
                ExEnv::out0() << indent << spinletters << " Vij_ji_coupling: R^ij_aa' f^a'_b T^ab_ji" << endl;
              compute_FxT(ij_ji, f12_idx,
                          i1i2a1AF2_ints, Tij_ab,
                          Vij_ji_coupling);
              // - R^ji_aa' f^a'_b T^ab_ij
              if (debug_ >= DefaultPrintThresholds::mostN2)
                ExEnv::out0() << indent << spinletters << " Vji_ij_coupling: + R^ji_aa' f^a'_b T^ab_ij" << endl;
              compute_FxT(ji_ij, f12_idx,
                          i1i2a1AF2_ints, Tij_ab,
                          Vij_ji_coupling);
          }
      } // end of V coupling computation
      i1i2a1AF2_ints->deactivate();
      if (spin1 != spin2)
        i1i2AF1a2_ints->deactivate();
      if (num_unique_spincases2 == 3 && spin1 != spin2) {
        i2i1AF1a2_ints->deactivate();
        i2i1a1AF2_ints->deactivate();
      }

      if (debug_ >= DefaultPrintThresholds::N2) {
        if (spin1 == spin2) {
            print_antisym_intermediate(spincase,"V^ij_ij coupling",Vij_ij_coupling,Vij_ji_coupling,nocc1_act,nocc2_act);
        } else {
            print_intermediate(spincase,"V^ij_ij coupling",Vij_ij_coupling,nocc1_act,nocc2_act);
            print_intermediate(spincase,"V^ij_ji coupling",Vij_ji_coupling,nocc1_act,nocc2_act);
        }
      } // end of debug

    } // end of V coupling computation

    //
    // Compute X intermediate matrix: X^ij_ij, X^ij_ji, X^ji_ij, X^ji_ji,
    //
    // Alpha_beta X^ij_ij and X^ij_ji are antisymmetrized
    // Alpha_alpha or beta_beta antisymmetrized X = X^ij_ij - X^ij_ji

    // Xij_ij = X^ij_ij = f^ij_ab f^ab_ij
    if (debug_ >= DefaultPrintThresholds::mostN2)
    ExEnv::out0() << endl << indent << spinletters << " X^ij_ij : " << endl;
    fill_n(Xij_ij, nocc12, 0.0);

    // X^ij_ij += (f12f12)^ij_ij
    Ref<DistArray4> i1i2i1i2_f12f12_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  occ1_act->id(), occ2_act->id(),
                  descr_f12f12_key, moints4_rtime,
                  i1i2i1i2_f12f12_ints);

    compute_VX(ij_ij, VX_output,
               f12f12_idx, i1i2i1i2_f12f12_ints,
               f12_idx, f12_idx,
               f12_ij_ints, f12_ij_ints,
               Xij_ij);

    //
    // Xij_ji = X^ij_ji = f^ij_ab f^ab_ji
    if (debug_ >= DefaultPrintThresholds::mostN2)
    ExEnv::out0() << endl << indent << spinletters << " X^ij_ji : " << endl;
    fill_n(Xij_ji, nocc12, 0.0);

    // X^ij_ji = (f12f12)^ij_ji
    Ref<DistArray4> i1i2i2i1_f12f12_ints;
    if (num_unique_spincases2 == 3 && spin1 != spin2){
      activate_ints(occ1_act->id(), occ2_act->id(),
                    occ2_act->id(), occ1_act->id(),
                    descr_f12f12_key, moints4_rtime,
                    i1i2i2i1_f12f12_ints);
    } else {
        i1i2i2i1_f12f12_ints = i1i2i1i2_f12f12_ints;
    }

    compute_VX(ij_ji,  VX_output,
               f12f12_idx, i1i2i2i1_f12f12_ints,
               f12_idx, f12_idx,
               f12_ij_ints, f12_ji_ints,
               Xij_ji);

    double* Xji_ij = NULL;
    double* Xji_ji = NULL;
    if (spin1 != spin2) {
      Xji_ij = new double[nocc12];
      Xji_ji = new double[nocc12];

      // Xji_ij = X^ji_ij = f^ji_ab f^ab_ij
      if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " X^ji_ij : " << endl;
      fill_n(Xji_ij, nocc12, 0.0);

      // X^ji_ij = (f12f12)^ji_ij
      Ref<DistArray4> i2i1i1i2_f12f12_ints;
      if (num_unique_spincases2 == 3){
        activate_ints(occ2_act->id(), occ1_act->id(),
                      occ1_act->id(), occ2_act->id(),
                      descr_f12f12_key, moints4_rtime,
                      i2i1i1i2_f12f12_ints);
      } else {
          i2i1i1i2_f12f12_ints = i1i2i1i2_f12f12_ints;
      }

      compute_VX(ji_ij,  VX_output,
                 f12f12_idx, i2i1i1i2_f12f12_ints,
                 f12_idx, f12_idx,
                 f12_ji_ints, f12_ij_ints,
                 Xji_ij);
      if (num_unique_spincases2 == 3)
        i2i1i1i2_f12f12_ints->deactivate();

      // Xji_ji = X^ji_ji = f^ji_ab f^ab_ji
      if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " X^ji_ji : " << endl;
      fill_n(Xji_ji, nocc12, 0.0);

      // X^ji_ji = (f12f12)^ji_ji
      Ref<DistArray4> i2i1i2i1_f12f12_ints;
      if (num_unique_spincases2 == 3){
        activate_ints(occ2_act->id(), occ1_act->id(),
                      occ2_act->id(), occ1_act->id(),
                      descr_f12f12_key, moints4_rtime,
                      i2i1i2i1_f12f12_ints);
      } else {
          i2i1i2i1_f12f12_ints = i1i2i1i2_f12f12_ints;
      }

      compute_VX(ji_ji,  VX_output,
                 f12f12_idx, i2i1i2i1_f12f12_ints,
                 f12_idx, f12_idx,
                 f12_ji_ints, f12_ji_ints,
                 Xji_ji);
      if (num_unique_spincases2 == 3)
        i2i1i2i1_f12f12_ints->deactivate();
    }

    if (debug_ >= DefaultPrintThresholds::N2 && spin1 == spin2)
      print_antisym_intermediate(spincase,"X^ij_ij",Xij_ij,Xij_ji,nocc1_act,nocc2_act);

    if (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2) {
      // Alpha-beta case
      print_intermediate(spincase,"X^ij_ij",Xij_ij,nocc1_act,nocc2_act);
      print_intermediate(spincase,"X^ij_ji",Xij_ji,nocc1_act,nocc2_act);
      print_intermediate(spincase,"X^ji_ij",Xji_ij,nocc1_act,nocc2_act);
      print_intermediate(spincase,"X^ji_ji",Xji_ji,nocc1_act,nocc2_act);
    }

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
    fill_n(Bij_ij, nocc12, 0.0);

    // B^ij_ij += (f12t1f12)^ij_ij
    if (debug_ >= DefaultPrintThresholds::mostN2)
    ExEnv::out0() << endl << indent << spinletters << " B(diag) contribution" << endl;
    compute_Y(ij_ij, 1.0,
              f12t1f12_idx,
              i1i2i1i2_f12f12_ints, Bij_ij);

    // P part of B intermediate
    fill_n(Pij_ij, nocc12, 0.0);

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
    activate_ints(occ1_act->id(), occ2_act->id(),
                  ribs1->id(), ribs2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2P1P2_ints);

    Ref<DistArray4> b_i1i2PK1P2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  Kribs1->id(), ribs2->id(),
                  descr_f12_key, moints4_rtime,
                  b_i1i2PK1P2_ints);

    // P += f^ji_PQ K^P_K Q _ji
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
    } else {
      //spin1 == spin2 or close shell
      b_i2i1P1P2_ints = b_i1i2P1P2_ints;
      b_i2i1PK1P2_ints = b_i1i2PK1P2_ints;

    }

    Pijij_output.push_back("P(f^ij_PQ K^P_R f^RQ_ij)");
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
    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << endl << indent << spinletters << " Q^ij_ij : " << endl;
    fill_n(Qij_ij, nocc12, 0.0);

    if (debug_ >= DefaultPrintThresholds::mostN2)
      ExEnv::out0() << indent << "Q += (f12f12)^ij_Pj (f+k)^P_i " << endl;
    Ref<DistArray4> q_i1i2hi1i2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(),
                  r12eval_->hj_i_P(spin1)->id(), occ2_act->id(),
                  descr_f12f12_key, moints4_rtime,
                  q_i1i2hi1i2_ints);
    compute_Y(ij_ij, 1.0,
              f12f12_idx,
              q_i1i2hi1i2_ints, Qij_ij);

    if (debug_ >= DefaultPrintThresholds::mostN2)
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

    // B^ij_ij = B^ij_ij - P + Q
    for(int i1=0; i1<nocc1_act; ++i1) {
      for(int i2=0; i2<nocc2_act; ++i2) {
        const int i1i2 = i1 * nocc2_act + i2;
        Bij_ij[i1i2] += - Pij_ij[i1i2] + Qij_ij[i1i2];
      }
    }
    if (debug_ >= DefaultPrintThresholds::mostN2
        || (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2))
      print_intermediate(spincase,"B^ij_ij",Bij_ij,nocc1_act,nocc2_act);

      //
      // B^ij_ji
      fill_n(Bij_ji, nocc12, 0.0);

      // B^ij_ji += (f12t1f12)^ij_ji
      if (debug_ >= DefaultPrintThresholds::mostN2)
        ExEnv::out0() << endl << indent << spinletters << " B(diag) contribution" << endl;
      compute_Y(ij_ji, 1.0,
                f12t1f12_idx,
                i1i2i2i1_f12f12_ints, Bij_ji);

      // P part of B intermediate
      fill_n(Pij_ji, nocc12, 0.0);

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
      if (debug_ >= DefaultPrintThresholds::mostN2)
        ExEnv::out0() << endl << indent << spinletters << " Q^ij_ji : " << endl;
      fill_n(Qij_ji, nocc12, 0.0);

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

      if (debug_ >= DefaultPrintThresholds::mostN2)
        ExEnv::out0() << endl;

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

      // B^ij_ji = B^ij_ji - P^ij_ji + Q^ij_ji
      for(int i1=0; i1<nocc1_act; ++i1) {
        for(int i2=0; i2<nocc2_act; ++i2) {
            const int i1i2 = i1 * nocc2_act + i2;
            Bij_ji[i1i2] += - Pij_ji[i1i2] + Qij_ji[i1i2];
        }
      }
      if (debug_ >= DefaultPrintThresholds::mostN2
          || (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2))
        print_intermediate(spincase,"B^ij_ji",Bij_ji,nocc1_act,nocc2_act);

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

          Bij_ij_beta = new double[nocc12];
          fill_n(Bij_ij_beta, nocc12, 0.0);

          // B^ij_ij += (f12t1f12)^ij_ij
          if (debug_ >= DefaultPrintThresholds::mostN2)
            ExEnv::out0() << endl << indent << spinletters << " B(diag) contribution (beta alpha case)" << endl;

          Ref<DistArray4> b_i2i1i2i1_ints;
          activate_ints(occ2_act->id(), occ1_act->id(),
                        occ2_act->id(), occ1_act->id(),
                        descr_f12f12_key, moints4_rtime,
                        b_i2i1i2i1_ints);
          compute_Y(ij_ij, 1.0,
                    f12t1f12_idx,
                    b_i2i1i2i1_ints, Bij_ij_beta);

          // P part of B intermediate
          double* Pij_ij_beta = new double[nocc12];
          fill_n(Pij_ij_beta, nocc12, 0.0);

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
          if (debug_ >= DefaultPrintThresholds::mostN2) {
            ExEnv::out0() << endl << indent << spinletters
                                    << " Q^ij_ij (beta alpha case) = Qji_ji (alpha beta case) "
                                    << endl;
          }

#if 0
          // The following codes is for testing
          ExEnv::out0() << indent << "Q^ij_ij (beta alpha case): " << endl;
          double* Qij_ij_beta = new double[nocc12];
          fill_n(Qij_ij_beta, nocc12, 0.0);

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
          for(int i2=0; i2<nocc2_act; ++i2) {
            for(int i1=0; i1<nocc1_act; ++i1) {
              const int i2i1 = i2 * nocc1_act + i1;
              const int i1i2 = i1 * nocc2_act + i2;
              Bij_ij_beta[i2i1] += - Pij_ij_beta[i2i1] + Qij_ij[i1i2];
            }
          }
          if (debug_ >= DefaultPrintThresholds::mostN2
              || (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2))
            print_intermediate(spincase,"Bij_ij (beta alpha case):",Bij_ij_beta,nocc2_act,nocc1_act);

            delete[] Pij_ij_beta;
            Pij_ij_beta = NULL;

            //
            // B^ij_ji
            // B^i(beta)j(alpha)_j(alpha)i(beta)
            Bij_ji_beta = new double[nocc12];
            fill_n(Bij_ji_beta, nocc12, 0.0);

            // B^ij_ji += (f12t1f12)^ij_ji
            if (debug_ >= DefaultPrintThresholds::mostN2)
              ExEnv::out0() << endl << indent << spinletters << " B(diag) contribution (beta alpha case):" << endl;
            Ref<DistArray4> b_i2i1i1i2_ints;
            activate_ints(occ2_act->id(), occ1_act->id(),
                          occ1_act->id(), occ2_act->id(),
                          descr_f12f12_key, moints4_rtime,
                          b_i2i1i1i2_ints);
            compute_Y(ij_ji, 1.0,
                      f12t1f12_idx,
                      b_i2i1i1i2_ints, Bij_ji_beta);

            // P part of B intermediate
            double* Pij_ji_beta = new double[nocc12];
            fill_n(Pij_ji_beta, nocc12, 0.0);

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
            if (debug_ >= DefaultPrintThresholds::mostN2) {
              ExEnv::out0() << endl << indent << spinletters
                                      << " Q^ij_ji (beta alpha case) = Qji_ij (alpha beta case) "
                                      << endl;
            }

#if 0
            // test
            ExEnv::out0() << indent << "Q^ij_ji (beta alpha case): " << endl;
            double* Qij_ji_beta = new double[nocc12];
            fill_n(Qij_ji_beta, nocc12, 0.0);

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
            for(int i2=0; i2<nocc2_act; ++i2) {
              for(int i1=0; i1<nocc1_act; ++i1) {
                const int i2i1 = i2 * nocc1_act + i1;
                const int i1i2 = i1 * nocc2_act + i2;
                Bij_ji_beta[i2i1] += - Pij_ji_beta[i2i1] + Qij_ji[i1i2];
              }
            }
            if (debug_ >= DefaultPrintThresholds::mostN2
                || (debug_ == DefaultPrintThresholds::N2 && spin1 != spin2))
              print_intermediate(spincase,"Bij_ji (beta alpha case):",Bij_ji_beta,nocc2_act,nocc1_act);

            delete[] Pij_ji_beta;
            Pij_ji_beta = NULL;

            // Deactive ints for B_beta
            b_i2i1i2i1_ints->deactivate();
            b_i2i1i1i2_ints->deactivate();
            // Deactive all the ints for Pij_ij_beta
            for (std::vector<Ref<DistArray4> >::iterator it = Pijij_beta_f12_ints.begin();
                it < Pijij_beta_f12_ints.end();++it) {
                (*it)->deactivate();
            }
            for (std::vector<Ref<DistArray4> >::iterator it = Pijij_beta_fx12_ints.begin();
                it < Pijij_beta_fx12_ints.end();++it) {
                (*it)->deactivate();
            }
      }

      // Deactive all the ints
      i1i2i1i2_ints->deactivate();
      for (std::vector<Ref<DistArray4> >::iterator it = f12_ij_ints.begin();
          it < f12_ij_ints.end();++it) {
          (*it)->deactivate();
      }
      i1i2i1i2_f12f12_ints->deactivate();

      if (num_unique_spincases2 == 3 && spin1 != spin2) {
          i1i2i2i1_ints->deactivate();
          for (std::vector<Ref<DistArray4> >::iterator it = f12_ji_ints.begin();
              it < f12_ji_ints.end();++it) {
              (*it)->deactivate();
          }
          i1i2i2i1_f12f12_ints->deactivate();
      }
      // Deactive ints for P
      for (std::vector<Ref<DistArray4> >::iterator it = Pijij_f12_ints.begin();
          it < Pijij_f12_ints.end();++it) {
          (*it)->deactivate();
          ++it;
          if (num_unique_spincases2 == 3 && spin1 != spin2) {
              (*it)->deactivate();
              }
          }
      for (std::vector<Ref<DistArray4> >::iterator it = Pijij_fx12_ints.begin();
          it < Pijij_fx12_ints.end();++it) {
          (*it)->deactivate();
          ++it;
          if (num_unique_spincases2 == 3 && spin1 != spin2) {
              (*it)->deactivate();
              }
          }
      // Deactive ints for Q
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
          print_antisym_intermediate(spincase,"Q^ij_ij:",Qij_ij,Qij_ji,nocc1_act,nocc2_act);
        if (debug_ >= DefaultPrintThresholds::N2)
          print_antisym_intermediate(spincase,"B^ij_ij:",Bij_ij,Bij_ji,nocc1_act,nocc2_act);
      }

      // Compute the f12 correction pair energy
      if (spin1 == spin2) {
        // Alpha_alpha or beta_beta case

        for(int i1=0; i1<nocc1_act; ++i1) {
          for(int i2=i1+1; i2<nocc2_act; ++i2) {

            double Hij_pair_energy;
            const int ij = i1 * nocc2_act + i2;

            if (debug_ >= DefaultPrintThresholds::mostN2) {
              ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(VT) = "
                            << (2.0 * C_1 * (Vij_ij[ij] - Vij_ji[ij]))
                            << endl;
              ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(TBT) = "
                            << (C_1*C_1 * (Bij_ij[ij] - Bij_ji[ij] - (evals_i1(i1) + evals_i2(i2)) * (Xij_ij[ij] - Xij_ji[ij])))
                            << endl;
            }

            Hij_pair_energy =  2.0 * C_1 * (Vij_ij[ij] - Vij_ji[ij])
                                   + C_1*C_1 * (Bij_ij[ij] - Bij_ji[ij] - (evals_i1(i1) + evals_i2(i2)) * (Xij_ij[ij] - Xij_ji[ij]));

            if (this->r12eval()->coupling() == true
                || this->r12eval()->ebc() == false) {

              if (debug_ >= DefaultPrintThresholds::mostN2) {
                ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(coupling) = "
                              << (2.0 * C_1 * ( Vij_ij_coupling[ij] - Vij_ji_coupling[ij]))
                              << endl;
              }

              Hij_pair_energy += 2.0 * C_1 * ( Vij_ij_coupling[ij] - Vij_ji_coupling[ij]);

            }

            // Get indices for the lower triangle matrix
            // index = nrow*(nrow-1) + ncolumn
            // nrow > ncolumn (diagonal elements are excluded)
            const int i21 = i2 * (i2-1)/2 + i1;
            ef12_[spin].set_element(i21, Hij_pair_energy);
          }
        }
      }  else if (num_unique_spincases2 == 2) {
          // Alpha_beta case for closed shell
          for(int i1=0; i1<nocc1_act; ++i1) {
            for(int i2=0; i2<nocc2_act; ++i2) {
              double Hij_pair_energy;
              int ij = i1 * nocc2_act + i2;

              if (debug_ >= DefaultPrintThresholds::mostN2) {
                ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(VT) = "
                              << (2.0 * ( 0.5*(C_0+C_1) * Vij_ij[ij] + 0.5*(C_0-C_1) * Vij_ji[ij]))
                              << endl;
                ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(TBT) = "
                              << (pow(0.5*(C_0+C_1), 2) * (Bij_ij[ij]- (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij])
                                  + 0.25*(C_0*C_0 - C_1*C_1)*(2.0*Bij_ji[ij] - (evals_i1(i1) + evals_i2(i2)) * (Xij_ji[ij] +Xji_ij[ij]))
                                  + pow(0.5*(C_0-C_1), 2) * (Bij_ij[ij] - (evals_i1(i1) + evals_i2(i2)) * Xji_ji[ij]))
                              << endl;
              }

              Hij_pair_energy =   2.0 * ( 0.5*(C_0+C_1) * Vij_ij[ij] + 0.5*(C_0-C_1) * Vij_ji[ij])
                                  + pow(0.5*(C_0+C_1), 2) * (Bij_ij[ij]- (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij])
                                  + 0.25*(C_0*C_0 - C_1*C_1)*(2.0*Bij_ji[ij] - (evals_i1(i1) + evals_i2(i2)) * (Xij_ji[ij] +Xji_ij[ij]))
                                  + pow(0.5*(C_0-C_1), 2) * (Bij_ij[ij] - (evals_i1(i1) + evals_i2(i2)) * Xji_ji[ij]);

              if (this->r12eval()->coupling() == true
                  || this->r12eval()->ebc() == false) {
                if (debug_ >= DefaultPrintThresholds::mostN2) {
                  ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(coupling) = "
                                << (2.0 * ( 0.5*(C_0+C_1) * Vij_ij_coupling[ij]
                                          + 0.5*(C_0-C_1) * Vij_ji_coupling[ij]))
                                << endl;
                }
                  Hij_pair_energy += 2.0 * ( 0.5*(C_0+C_1) * Vij_ij_coupling[ij]
                                           + 0.5*(C_0-C_1) * Vij_ji_coupling[ij]);
              }

              const int i12 = i1 * nocc2_act + i2;
              ef12_[spin].set_element(i12, Hij_pair_energy);
            }
          }
      } else {
          // Alpha_beta case for open shell
          for(int i1=0; i1<nocc1_act; ++i1) {
            for(int i2=0; i2<nocc2_act; ++i2) {
              double Hij_pair_energy;
              int ij = i1 * nocc2_act + i2;
              int ji = i2 * nocc1_act + i1;

              if (debug_ >= DefaultPrintThresholds::mostN2) {
                ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(VT) = "
                    << (2.0 * ( 0.5*(C_0+C_1) * Vij_ij[ij] + 0.5*(C_0-C_1) * Vij_ji[ij]))
                    << endl;
                ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(TBT) = "
                    << (pow(0.5*(C_0+C_1), 2) * ((Bij_ij[ij] +Bij_ij_beta[ji])*0.5 - (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij])
                        + 0.25*(C_0*C_0 - C_1*C_1)*((Bij_ji[ij] +Bij_ji_beta[ji]) - (evals_i1(i1) + evals_i2(i2)) * (Xij_ji[ij] +Xji_ij[ij]))
                        + pow(0.5*(C_0-C_1), 2) * ((Bij_ij[ij] +Bij_ij_beta[ji])*0.5 - (evals_i1(i1) + evals_i2(i2)) * Xji_ji[ij]))
                        << endl;
              }
              Hij_pair_energy =   2.0 * ( 0.5*(C_0+C_1) * Vij_ij[ij] + 0.5*(C_0-C_1) * Vij_ji[ij])
                                + pow(0.5*(C_0+C_1), 2) * ((Bij_ij[ij] +Bij_ij_beta[ji])*0.5 - (evals_i1(i1) + evals_i2(i2)) * Xij_ij[ij])
                                + 0.25*(C_0*C_0 - C_1*C_1)*((Bij_ji[ij] +Bij_ji_beta[ji]) - (evals_i1(i1) + evals_i2(i2)) * (Xij_ji[ij] +Xji_ij[ij]))
                                + pow(0.5*(C_0-C_1), 2) * ((Bij_ij[ij] +Bij_ij_beta[ji])*0.5 - (evals_i1(i1) + evals_i2(i2)) * Xji_ji[ij]);

              if (this->r12eval()->coupling() == true
                  || this->r12eval()->ebc() == false) {
                if (debug_ >= DefaultPrintThresholds::mostN2) {
                  ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(coupling) = "
                      << (2.0 * ( 0.5*(C_0+C_1) * Vij_ij_coupling[ij]
                                + 0.5*(C_0-C_1) * Vij_ji_coupling[ij]))
                      << endl;
                }
                Hij_pair_energy += 2.0 * ( 0.5*(C_0+C_1) * Vij_ij_coupling[ij]
                                         + 0.5*(C_0-C_1) * Vij_ji_coupling[ij]);
              }

              const int i12 = i1 * nocc2_act + i2;
              ef12_[spin].set_element(i12, Hij_pair_energy);
            }
          }
          // deallocate memory for B_beta (it will be used once)
          delete[] Bij_ij_beta;
          delete[] Bij_ji_beta;
          delete[] Xji_ij;
          delete[] Xji_ji;
      }

  } // end of spincase loop

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

  //
  // Coupled cluster V contribution
  //
  if (r12intermediates_->T1_cc_computed() &&
      r12intermediates_->T2_cc_computed() ) {

      if (debug_ >= DefaultPrintThresholds::N2)
        ExEnv::out0() << endl << indent << "Coupled-cluster V contribution" << endl;
    // Alpha-alpha, beta-beta, or close-shell case:
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
    for (int s=0; s<num_unique_spincases2; s++) {
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
      const blasint nv12 = nv1 *nv2;
      const blasint one = 1;

      // Vpq_ij
      std::vector< Ref<DistArray4> > Vpq_vec = r12eval()->V_distarray4(spincase2, p1, p2);
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
      compute_YxF(ij_ij, 1.0,
                  0, 0,
                  Vab_ij, T2[s],
                  VT2ij_ij);
      if (debug_ >= DefaultPrintThresholds::N2)
        print_intermediate(spincase2,"V^ab_ij * T^ij_ab",VT2ij_ij,no1,no2);

      // VT2ij_ji = V^ab_ji * T^ij_ab
      if (spincase2 == AlphaBeta) {
        fill_n(VT2ij_ji, nijij, 0.0);

        if (num_unique_spincases2 == 3) {
          // V^qp_ij (V^p2p1_i1i2)
          std::vector< Ref<DistArray4> > Vqp_vec = r12eval()->V_distarray4(spincase2, p2, p1);
          MPQC_ASSERT(Vqp_vec.size() == 1);
          Vqp_ij = Vqp_vec[0];

          // extract Vba_ij from Vqp_ij
          Ref<DistArray4> Vba_ij;
          map(Vqp_ij, o1, o2, p2, p1, Vba_ij, o1, o2, v2, v1);

          Vba_ij->activate();
          double* const Vab_blk = new double[nv12];
          fill(Vab_blk, Vab_blk+nv12, 0.0);
          for (int i=0; i<no1; ++i) {
            for (int j=0; j<no2; ++j) {
              const double*  Vba_blk = Vba_ij->retrieve_pair_block(i, j, 0);
              swap_e12(Vba_blk,nv2,nv1,Vab_blk);

              const double* T2ab_blk = T2[s]->retrieve_pair_block(i, j, 0);
              const double vt2ijji = F77_DDOT(&nv12, Vab_blk, &one, T2ab_blk, &one);
              const int ij = i*no2+j;
              VT2ij_ji[ij] = vt2ijji;

              Vba_ij->release_pair_block(i, j, 0);
              T2[s]->release_pair_block(i, j, 0);
            }
          }
          delete[] Vab_blk;
          Vba_ij->deactivate();
          //print_intermediate(spincase2,VT2ij_ji,no1,no2);
        } else {
            compute_YxF(ji_ij, 1.0,
                        0, 0,
                        Vab_ij, T2[s],
                        VT2ij_ji);
        }
        if (debug_ >= DefaultPrintThresholds::N2)
          print_intermediate(spincase2,"V^ab_ji * T^ij_ab",VT2ij_ji,no1,no2);
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
      contract_VT1(Via_ij,
                   ij_ij, nonswap_e12_V,
                   raw_T1_ja,
                   nv2, T1_ja,
                   VT1ij_ij);
      if (debug_ >= DefaultPrintThresholds::N2)
        print_intermediate(spincase2,"VT1ij_ij = Via_ij * T1^j_a",VT1ij_ij,no1,no2);

      // VT1ij_ji = V^ia_ji * T^ij_ab
      if (spincase2 == AlphaBeta) {
        fill_n(VT1ij_ji, nijij, 0.0);

        if (num_unique_spincases2 == 3) {
          // extract Vai_ij from Vqp_ij
          Ref<DistArray4> Vai_ij;
          map(Vqp_ij, o1, o2, p2, p1, Vai_ij, o1, o2, v2, o1);

          Vai_ij->activate();
          contract_VT1(Vai_ij,
                       ij_ij, swap_e12_V,
                       raw_T1_ja,
                       nv2, T1_ja,
                       VT1ij_ji);
          Vai_ij->deactivate();
        } else {
            contract_VT1(Via_ij,
                         ji_ij, nonswap_e12_V,
                         raw_T1_ja,
                         nv2, T1_ja,
                         VT1ij_ji);
        }
        if (debug_ >= DefaultPrintThresholds::N2)
          print_intermediate(spincase2,"VT1ij_ji = V^ia_ji * T1^j_a",VT1ij_ji,no1,no2);
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
      contract_VT1(Vaj_ij,
                   ij_ij, swap_e12_V,
                   raw_T1_ia,
                   nv1, T1_ia,
                   VT1ij_ij);
      Vaj_ij->deactivate();
      if (debug_ >= DefaultPrintThresholds::N2)
        print_intermediate(spincase2,"VT1ij_ij += V^aj_ij * T1^i_a",VT1ij_ij,no1,no2);

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
        contract_VT1(Vja_ij,
                     ij_ij, nonswap_e12_V,
                     raw_T1_ia,
                     nv1, T1_ia,
                     VT1ij_ji);
        Vja_ij->deactivate();
        if (debug_ >= DefaultPrintThresholds::N2)
          print_intermediate(spincase2,"VT1ij_ji += V^ja_ij * T1^i_a",VT1ij_ji,no1,no2);
      }
      delete[] raw_T1_ia;

      // Add VT2 & VT1 contribution to ef12
      if (spin1 == spin2) {
        for(int i1=0; i1<no1; ++i1) {
          for(int i2=i1+1; i2<no2; ++i2) {
            const int ij = i1 * no2 + i2;
            const int i21 = i2 * (i2-1)/2 + i1;
            // 2.0 from Hylleraas functional
            const double Hij_pair_energy = 2.0*C_1*(0.5*VT2ij_ij[ij] + VT1ij_ij[ij]);
//          const double Hij_pair_energy = VT2ij_ij[ij];
            if (debug_ >= DefaultPrintThresholds::N2) {
              ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(CC) = "
                            << Hij_pair_energy << endl;
            }
            ef12_[s].accumulate_element(i21, Hij_pair_energy);
          }
        }
      } else {
          // Alpha_beta case
          for(int i1=0; i1<no1; ++i1) {
            for(int i2=0; i2<no2; ++i2) {
              const int ij = i1 * no2 + i2;
              const double Hij_pair_energy = 2.0*( 0.5*(C_0+C_1) * VT2ij_ij[ij] + 0.5*(C_0-C_1) * VT2ij_ji[ij]
                                                 + 0.5*(C_0+C_1) * VT1ij_ij[ij] + 0.5*(C_0-C_1) * VT1ij_ji[ij]
                                                   );
              if (debug_ >= DefaultPrintThresholds::N2) {
                ExEnv::out0() << indent << "Hij_pair_energy: ij = " << i1 << "," << i2 << " e(CC) = "
                              << Hij_pair_energy << endl;
              }
              ef12_[s].accumulate_element(ij, Hij_pair_energy);
            }
          }
      }
    } // end of spin iteration
    delete[] VT2ij_ij;
    delete[] VT2ij_ji;
    delete[] VT1ij_ij;
    delete[] VT1ij_ji;
  } // end of CC V contribution

  // Set beta-beta energies to alpha-alpha for closed-shell
  if (!r12world->refwfn()->spin_polarized()) {
    C_[BetaBeta] = C_[AlphaAlpha];
    emp2f12_[BetaBeta] = emp2f12_[AlphaAlpha];
    ef12_[BetaBeta] = ef12_[AlphaAlpha];
  }

  evaluated_ = true;

  return;
}
