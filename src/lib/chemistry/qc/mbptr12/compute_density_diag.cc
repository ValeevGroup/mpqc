/*
 * compute_density_diag.cc
 *
 *  Created on: Nov 8, 2011
 *      Author: jinmei
 */

#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>
#include <chemistry/qc/psi/psiwfn.h>
#include <util/misc/print.h>
#include <math/scmat/blas.h>

using namespace std;
using namespace sc;

// case1: R^idx1 idx3_ij * R^ij_idx2 idx3, eg: R^ba'_ij R^ij_ca'
// case2: R^idx3 idx1_ij * R^ij_idx2 idx3, eg: R^a'b_ij R^ij_ca'
// case3: R^idx1 idx3_ij * R^ij_idx3 idx2, eg: R^ba'_ij R^ij_a'c
// case4: R^idx3 idx1_ij * R^ij_idx3 idx2, eg: R^a'b_ij R^ij_a'c
enum RRidx_cases{RR13_23, RR31_23, RR13_32, RR31_32};
enum RT2idx_cases{RT13_23, RT31_23, RT13_32, RT31_32};

// orbs1_orbs2_orbs3: D^orbs1_orbs2 sum over orbs2
enum orbitals_cases{vir_vir_cabs, cabs_cabs_vir, cabs_cabs_cabs, cabs_vir_cabs};

namespace {
  // print tension with 1 index
  void print_intermediate(const string& spinlabel, const string& label,
                          const double* const array, const int size_idx)
  {
     ExEnv::out0() << indent << spinlabel << " " << label << endl;

     const int col_print = 6;

     // print the column label
     const int size_col = (size_idx > col_print? col_print : size_idx);
     for (int i = 1; i <= size_col; ++i) {
         ExEnv::out0() << indent << scprintf("%12d", i) << "  ";
     }
//     ExEnv::out0() << endl << indent << 1 << " " ;

     const double* iter = array;
     for (int idx = 1 ; idx <= size_idx; ++idx) {
       ExEnv::out0() << indent << scprintf("%12.10f", *iter) << "  ";

       // print the column label for the next 5 elements
       const int left_col = idx % col_print;
       if (left_col == 0) {
         ExEnv::out0() << endl ;
         const int size_col2 = ((idx + col_print) > size_idx? size_idx : idx + col_print);
         for (int i = idx + 1; i <= size_col2; ++i) {
           ExEnv::out0() << indent << scprintf("%12d", i) << "  ";
         }
         ExEnv::out0() << endl << indent << 1 << " ";
       }

       ++iter;
     } // end of loop
     ExEnv::out0() << endl << endl;
  }
  // end of print function

  // print tension with 2 indices
  void print_intermediate(const string& spinlabel, const string& label,
                          const double* const array,
                          const int size_idx1, const int size_idx2)
  {
     ExEnv::out0() << indent << spinlabel << " " << label << endl;

     const int col_print = 6;
     const int size_col1 = (size_idx2 > col_print? col_print : size_idx2);

     if (size_idx2 <= col_print) {
       // print the column label
       for (int i = 1; i <= size_col1; ++i) {
         ExEnv::out0() << indent << scprintf("%12d", i) << "  ";
       }
       // print element
       const double* iter = array;
       for (int idx1 = 1 ; idx1 <= size_idx1; ++idx1) {
         ExEnv::out0() << endl << indent << idx1 << " ";

         for (int idx2 = 1; idx2 <= size_idx2; ++idx2) {
           ExEnv::out0() << indent << scprintf("%12.10f", *iter) << "  ";
           ++iter;
         }
         ExEnv::out0() << endl;
       }

       ExEnv::out0() << endl;
     } else {
         int left_col = size_idx2;
         int iteration_idx = 0;

         while (left_col > 0) {
           const int size_col = (left_col > col_print? col_print : left_col);

           // print the column label
           const int start_idx2 = iteration_idx * col_print;
           for (int idx2 = 1; idx2 <= size_col; ++idx2) {
             ExEnv::out0() << indent << scprintf("%12d", idx2 + start_idx2) << "  ";
           }
           ExEnv::out0() << endl;

           // print element
           for (int idx1 = 1 ; idx1 <= size_idx1; ++idx1) {
             ExEnv::out0() << indent << idx1 << " ";

             for (int idx2 = 1; idx2 <= size_col; ++idx2) {
               const int index = (idx1 - 1) * size_idx2 + (idx2 + start_idx2 - 1);
               ExEnv::out0() << indent << scprintf("%12.10f", array[index]) << "  ";
             }
             ExEnv::out0() << endl;
           }

           ++iteration_idx;
           left_col = size_idx2 - iteration_idx * col_print;
         } // end of while

     }

     ExEnv::out0() << endl;
  }
  // end of print function

  // functions needed for compute_Dii :
  // compute (f12f12)^b1b2_k1k2 sum over index 3
  // index 1 and 2 are i, index 3 is j: eg, (f12f12)^ij_ji
  void compute_F12F12_sum_idx3(const int RRb1b2_k1k2,
                               const unsigned int f12f12_idx, Ref<DistArray4>& f12f12_ints,
                               double* result)
  {
    // index1, 2: i, index3: j
    int index1 = 0;
    int index3 = 0;
    int size_idx1 = f12f12_ints->ni();
    int size_idx3 = f12f12_ints->nj();
    const int size_k2 = f12f12_ints->ny();

    int* f12f12_ints_idx1 = &index1;
    int* f12f12_ints_idx2 = &index3;

    int* f12f12_blk_idx1 = &index1;
    int* f12f12_blk_idx2 = &index3;

    switch (RRb1b2_k1k2) {
    case RR13_23:
      // default value
    break;

    case RR31_23:
      size_idx1 = f12f12_ints->nj();
      size_idx3 = f12f12_ints->ni();
      f12f12_ints_idx1 = &index3;
      f12f12_ints_idx2 = &index1;
    break;

    case RR13_32:
      f12f12_blk_idx1 = &index3;
      f12f12_blk_idx2 = &index1;
    break;

    case RR31_32:
      size_idx1 = f12f12_ints->nj();
      size_idx3 = f12f12_ints->ni();
      f12f12_ints_idx1 = &index3;
      f12f12_ints_idx2 = &index1;
      f12f12_blk_idx1 = &index3;
      f12f12_blk_idx2 = &index1;
    break;

    default:
      ExEnv::out0() << "There is no such index";
      break;
    }

    // test
    //ExEnv::out0() << "f12f12 part in compute_F12F12_sum_idx3" << endl;

    double* iter_result = result;
    for(index1 = 0; index1 < size_idx1; ++index1) {

      double f12f12_sum_idx3 = 0;
      for(index3 = 0; index3 < size_idx3; ++index3) {

        const double* f12f12_blk = f12f12_ints->retrieve_pair_block(*f12f12_ints_idx1, *f12f12_ints_idx2, f12f12_idx);
        const int element = (*f12f12_blk_idx1) * size_k2 + *f12f12_blk_idx2;
        f12f12_sum_idx3 += f12f12_blk[element];

        // test
        //ExEnv::out0() << indent << scprintf("%12.10f", f12f12_blk[element]) << " ";

        f12f12_ints->release_pair_block(*f12f12_ints_idx1, *f12f12_ints_idx2, f12f12_idx);
      }
      *iter_result = f12f12_sum_idx3;
      ++iter_result;

      // test
      //ExEnv::out0() << endl;
    }
  }
  // end of function: compute_F12F12_sum_idx3

  // functions needed for computing D^m_i :
  // compute (f12f12)^b1b2_k1k2 sum over index 3
  // index 1 and 2 are i, m, index 3 is j: eg, (f12f12)^ij_mj
  void compute_F12F12_sum_idx3_2(const int RRb1b2_k1k2,
                                 const unsigned int f12f12_idx, const Ref<DistArray4>& f12f12_ints,
                                 double* const result)
  {
    // index1, 2: i m, index3: j
    int index1 = 0;
    int index2 = 0;
    int index3 = 0;
    int size_idx1 = f12f12_ints->ni();
    int size_idx2 = f12f12_ints->nx();
    int size_idx3 = f12f12_ints->nj();
    // size of ket2
    const int size_k2 = f12f12_ints->ny();

    int* f12f12_ints_idx1 = &index1;
    int* f12f12_ints_idx2 = &index3;

    int* f12f12_blk_idx1 = &index2;
    int* f12f12_blk_idx2 = &index3;

    switch (RRb1b2_k1k2) {
    case RR13_23:
      // default value
    break;

    case RR31_23:
      size_idx1 = f12f12_ints->nj();
      size_idx3 = f12f12_ints->ni();
      f12f12_ints_idx1 = &index3;
      f12f12_ints_idx2 = &index1;
    break;

    case RR13_32:
      size_idx2 = f12f12_ints->ny();
      f12f12_blk_idx1 = &index3;
      f12f12_blk_idx2 = &index2;
    break;

    case RR31_32:
      size_idx1 = f12f12_ints->nj();
      size_idx2 = f12f12_ints->ny();
      size_idx3 = f12f12_ints->ni();
      f12f12_ints_idx1 = &index3;
      f12f12_ints_idx2 = &index1;
      f12f12_blk_idx1 = &index3;
      f12f12_blk_idx2 = &index2;
    break;

    default:
      ExEnv::out0() << "There is no such index";
      break;
    }

    // test
    //ExEnv::out0() << "f12f12 part in compute_F12F12_sum_idx3" << endl;

    double* iter_result = result;
    for(index1 = 0; index1 < size_idx1; ++index1) {
      for(index2 = 0; index2 < size_idx2; ++index2) {

        double f12f12_sum_idx3 = 0;
        for(index3 = 0; index3 < size_idx3; ++index3) {
          const double* f12f12_blk = f12f12_ints->retrieve_pair_block(*f12f12_ints_idx1, *f12f12_ints_idx2, f12f12_idx);
          const int element = (*f12f12_blk_idx1) * size_k2 + *f12f12_blk_idx2;
          f12f12_sum_idx3 += f12f12_blk[element];

          // test
          //ExEnv::out0() << indent << scprintf("%12.10f", f12f12_blk[element]) << " ";
          f12f12_ints->release_pair_block(*f12f12_ints_idx1, *f12f12_ints_idx2, f12f12_idx);
        }
        *iter_result = f12f12_sum_idx3;
        ++iter_result;

        // test
        //ExEnv::out0() << endl;
      }
    }
  }
  // end of function: compute_F12F12_sum_idx3_2

  // compute R^b1b2_ab * R^ab_k1k2 contracted over a, b, and j
  // b1, b2, k1, k2 can be i or j (only two indices)
  void compute_RR1_sum_3idx(const int RRb1b2_k1k2, const unsigned int f12_idx,
                           Ref<DistArray4>& f12_ints1, Ref<DistArray4>& f12_ints2,
                           double* f12f12_array)
  {
    // idx1: external indices
    // idx3: the other contracted index
    // different from compute_RR2_sum_3idx, there is no index2
    // eg: R^ij_ab R^ab_ij: idx1 is i, and idx3 is j
    int idx1 = 0;
    int idx3 = 0;

    int size_idx1 = f12_ints1->ni();
    int size_idx3 = f12_ints1->nj();
    // get the number of the contracted indices
    // eg: R^ij_ab R^ab_ij: sum_idx is ab
    const blasint size_sum_idx = f12_ints1->nx() * f12_ints1->ny();

    // get block indices
    int* f12_blk1_idx1 = &idx1;
    int* f12_blk1_idx2 = &idx3;
    int* f12_blk2_idx1 = &idx1;
    int* f12_blk2_idx2 = &idx3;

    switch (RRb1b2_k1k2) {
    case RR13_23:
      // default value
    break;

    case RR31_23:
      size_idx1 = f12_ints1->nj();
      size_idx3 = f12_ints1->ni();
      f12_blk1_idx1 = &idx3;
      f12_blk1_idx2 = &idx1;
    break;

    case RR13_32:
      f12_blk2_idx1 = &idx3;
      f12_blk2_idx2 = &idx1;
    break;

    case RR31_32:
      size_idx1 = f12_ints1->nj();
      size_idx3 = f12_ints1->ni();
      f12_blk1_idx1 = &idx3;
      f12_blk1_idx2 = &idx1;
      f12_blk2_idx1 = &idx3;
      f12_blk2_idx2 = &idx1;
    break;

    default:
      ExEnv::out0() << "There is no such index case";
      break;
    }

    const blasint one = 1; // for F77_DDOT
    double* iter_array = f12f12_array;

    for (idx1 = 0; idx1 < size_idx1; ++idx1) {

      double f12f12_sum_idx3 = 0;
      for(idx3 = 0; idx3 < size_idx3; ++idx3) {
        const double* f12_blk1 = f12_ints1->retrieve_pair_block((*f12_blk1_idx1), (*f12_blk1_idx2), f12_idx);
        const double* f12_blk2 = f12_ints2->retrieve_pair_block((*f12_blk2_idx1), (*f12_blk2_idx2), f12_idx);

        const double f12f12_sum_idx = F77_DDOT(&size_sum_idx, f12_blk1, &one, f12_blk2, &one);

        f12f12_sum_idx3 += f12f12_sum_idx;

        f12_ints1->release_pair_block((*f12_blk1_idx1), (*f12_blk1_idx2), f12_idx);
        f12_ints2->release_pair_block((*f12_blk2_idx1), (*f12_blk2_idx2), f12_idx);
      }
      *iter_array = f12f12_sum_idx3;
      ++iter_array;
    }
  }
  // end of compute_RR1_sum_3dix

  // function needed for compute_Dii_test
  void compute_Dii_spin(const SpinCase1 spin, const double C_0, const double C_1,
                        int nocc_alpha, int nocc_beta,
                        const double* RRij_ij, const double* RRji_ij,
                        const double* RRi1i2_i1i2, const double* RRi1i2_i2i1,
                        const double* RRi2i1_i1i2, const double* RRi2i1_i2i1)
  {
    int* nocc = &nocc_alpha;

    // get the right element of AlphaBeta R*R for alpha/beta spin
    int* size_j_ab = &nocc_beta;

    int start_alpha = 0;
    int start_beta = 0;
    int* start_RR = &start_alpha;

    int offset_alpha = 1;
    int offset_beta = nocc_beta;
    int* offset_RR = &offset_alpha;

    if (spin == Beta) {
      nocc = &nocc_beta;
      size_j_ab = &nocc_alpha;
      start_RR = &start_beta;
      offset_RR = &offset_beta;
    }

    ExEnv::out0() << "C_0: "  << C_0 << " C_1: " << C_1 << endl;
    // spin=0 (Alpha) => AlphaAlpha case (1)
    // spin=1 (Beta) => BetaBeta case (2)
    const int s = spin;
    const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
    string spinletters = to_string(spincase);

    const double* iter_RRij_ij = RRij_ij;
    const double* iter_RRji_ij = RRji_ij;
    for (int i = 0;  i < (*nocc); ++i) {
      double dii_12 = 0.0;

      // AlphaAlpha/BetaBeta part
      double dii_1 = 0.0;
      for (int j = 0; j < (*nocc); ++j) {
        dii_1 += *iter_RRij_ij - *iter_RRji_ij;

//        ExEnv::out0() << "j = " << j << " dii_1 = "
//                      << scprintf("%12.10f", dii_1) << endl;

        ++iter_RRij_ij;
        ++iter_RRji_ij;
      }
      dii_12 += C_1 * C_1 * dii_1;
      ExEnv::out0() << spinletters << " part d^" << i << "_" << i << " = "
                    << scprintf("%12.10f", dii_12) << endl;

      // AlphaBeta part
      if (nocc_alpha != 0 && nocc_beta != 0) {

        start_alpha = i * nocc_beta;
        start_beta = i;
        const double* iter_RRi1i2_i1i2 = RRi1i2_i1i2 + (*start_RR);
        const double* iter_RRi2i1_i1i2 = RRi2i1_i1i2 + (*start_RR);
        const double* iter_RRi1i2_i2i1 = RRi1i2_i2i1 + (*start_RR);
        const double* iter_RRi2i1_i2i1 = RRi2i1_i2i1 + (*start_RR);

        double RRi1i2i1i2 = 0.0;
        double RRi2i1i1i2 = 0.0;
        double RRi1i2i2i1 = 0.0;
        double RRi2i1i2i1 = 0.0;
        for (int j = 0; j < (*size_j_ab); ++j) {
          RRi1i2i1i2 += *iter_RRi1i2_i1i2;
          RRi2i1i1i2 += *iter_RRi2i1_i1i2;
          RRi1i2i2i1 += *iter_RRi1i2_i2i1;
          RRi2i1i2i1 += *iter_RRi2i1_i2i1;

          iter_RRi1i2_i1i2 += (*offset_RR);
          iter_RRi2i1_i1i2 += (*offset_RR);
          iter_RRi1i2_i2i1 += (*offset_RR);
          iter_RRi2i1_i2i1 += (*offset_RR);
        }
//        ExEnv::out0() << "RRi1i2i1i2: "  << scprintf("%12.10f", RRi1i2i1i2)
//                      << "  RRi2i1i1i2: " << scprintf("%12.10f", RRi2i1i1i2)
//                      << "  RRi1i2i2i1: "  << scprintf("%12.10f", RRi1i2i2i1)
//                      << "  RRi2i1i2i1: " << scprintf("%12.10f", RRi2i1i2i1)
//                      << endl;
        dii_12 += pow(0.5 * (C_0 + C_1), 2) * RRi1i2i1i2
               + 0.25 * (C_0 * C_0 - C_1 * C_1) * (RRi2i1i1i2 + RRi1i2i2i1)
               + pow(0.5 * (C_0 - C_1), 2) * RRi2i1i2i1;
        ExEnv::out0() << "AlphaBeta part: d^"  << i << "_" << i << " = "
                      << scprintf("%12.10f", dii_12) << endl;
      } // end of summing AlphaBeta X

    }
  }
  // end of function: compute_Dii_spin

  // compute R * R contracted over 3 indices: i j index 3 or a b index 3
  void compute_RR2_sum_3idx(const int RRb1b2_k1k2, const unsigned int f12_idx,
                           const Ref<DistArray4>& f12_ints1, const Ref<DistArray4>& f12_ints2,
                           double* const f12f12_array)
  {
    // idx1 and idx2 are label indices, idx3 is the other contracted index
    // eg: R^a'b_ij R^ij_a'c: idx1 is b, idx2 is c, and idx3 is a'
    int idx1 = 0;
    int idx2 = 0;
    int idx3 = 0;

    int size_idx1 = f12_ints1->ni();
    int size_idx2 = f12_ints2->ni();
    int size_idx3 = f12_ints1->nj();
    // get the number of the contracted indices
    // eg: R^a'b_ij R^ij_a'c: sum_idx is ij
    const blasint size_sum_idx = f12_ints1->nx() * f12_ints1->ny();

    // get block indices
    int* f12_blk1_idx1 = &idx1;
    int* f12_blk1_idx2 = &idx3;
    int* f12_blk2_idx1 = &idx2;
    int* f12_blk2_idx2 = &idx3;

    switch (RRb1b2_k1k2) {
    case RR13_23:
      // default value
    break;

    case RR31_23:
      size_idx1 = f12_ints1->nj();
      size_idx3 = f12_ints1->ni();
      f12_blk1_idx1 = &idx3;
      f12_blk1_idx2 = &idx1;
    break;

    case RR13_32:
      size_idx2 = f12_ints2->nj();
      f12_blk2_idx1 = &idx3;
      f12_blk2_idx2 = &idx2;
    break;

    case RR31_32:
      size_idx1 = f12_ints1->nj();
      size_idx2 = f12_ints2->nj();
      size_idx3 = f12_ints1->ni();
      f12_blk1_idx1 = &idx3;
      f12_blk1_idx2 = &idx1;
      f12_blk2_idx1 = &idx3;
      f12_blk2_idx2 = &idx2;
    break;

    default:
      ExEnv::out0() << "There is no such index case";
      break;
    }
    // test
//    ExEnv::out0() << "size of indices: " << size_idx1 << size_idx2 << size_idx3 << endl;

    const blasint one = 1; // for F77_DDOT
    double* iter_array = f12f12_array;

      for (idx1 = 0; idx1 < size_idx1; ++idx1) {
          // can be change to idx2 = idx1 for d^b_c and d^b'_c'
        for(idx2 = 0; idx2 < size_idx2; ++idx2) {

          double f12f12_sum_idx3 = 0;
          for(idx3 = 0; idx3 < size_idx3; ++idx3) {
            const double* f12_ij_blk1 = f12_ints1->retrieve_pair_block((*f12_blk1_idx1), (*f12_blk1_idx2), f12_idx);
            const double* f12_ij_blk2 = f12_ints2->retrieve_pair_block((*f12_blk2_idx1), (*f12_blk2_idx2), f12_idx);

            const double f12f12_sum_idx = F77_DDOT(&size_sum_idx, f12_ij_blk1, &one, f12_ij_blk2, &one);

            f12f12_sum_idx3 += f12f12_sum_idx;

            f12_ints1->release_pair_block((*f12_blk1_idx1), (*f12_blk1_idx2), f12_idx);
            f12_ints2->release_pair_block((*f12_blk2_idx1), (*f12_blk2_idx2), f12_idx);
          }
          *iter_array = f12f12_sum_idx3;
//          ExEnv::out0() << *iter_array << " ";
          ++iter_array;
        }
//        ExEnv::out0() << endl;
      }
  }
  // end of compute_RR2_sum_3dix

  // compute R * T2 (CC amplitudes, stored as ij, ab) contracted over 3 indices
  // here they are: i, j, idx3
  void compute_RT2_sum_3idx(const int RTb1b2_k1k2, const unsigned int f12_idx,
                            Ref<DistArray4>& F12_ints, Ref<DistArray4>& T2_ints,
                            double* const RT2)
  {
    const int size_i = F12_ints->ni();
    const int size_j = F12_ints->nj();

    int size_ap = F12_ints->nx();
    int size_a = T2_ints->nx();
    int size_b = F12_ints->ny();

    int Rapb_blk_start = 0;
    int Rbap_blk_start = 0;
    int T2ab_blk_start = 0;
    int T2ba_blk_start = 0;
    int* F12_blk_start = &Rapb_blk_start;
    int* T2_blk_start = &T2ab_blk_start;

    // R^ij_a'b * T2^ab_ij
    int* F12_blk_offset = &size_b;
    int* T2_blk_offset = &size_b;
    int one = 1;

    switch (RTb1b2_k1k2) {
    case RT13_23:
      // default value
    break;

    case RT31_23:
      // R^ij_ba' * T2^ab_ij
      size_ap = F12_ints->ny();
      size_b = F12_ints->nx();
      F12_blk_start = &Rbap_blk_start;
      F12_blk_offset = &one;
    break;

    case RT13_32:
      // R^ij_a'b * T2^ba_ij
      size_a = T2_ints->ny();
      T2_blk_start = &T2ba_blk_start;
      T2_blk_offset = &one;
    break;

    case RT31_32:
      // R^ij_ba' * T2^ba_ij
      size_ap = F12_ints->ny();
      size_b = F12_ints->nx();
      size_a = T2_ints->ny();
      F12_blk_start = &Rbap_blk_start;
      F12_blk_offset = &one;
      T2_blk_start = &T2ba_blk_start;
      T2_blk_offset = &one;
    break;

    default:
      ExEnv::out0() << "There is no such index case";
      break;
    }

//    ExEnv::out0() << "compute_RT2_sum_3idx: ap a b " << size_ap << size_a << size_b << endl;

    const int ncabs_vir = size_ap * size_a;
    fill_n(RT2, ncabs_vir, 0.0);

    for (int i = 0; i < size_i; ++i) {
      for(int j = 0; j < size_j; ++j) {
        const double* const F12_blk = F12_ints->retrieve_pair_block(i, j, f12_idx);
        const double* const T2_blk = T2_ints->retrieve_pair_block(i, j, 0);

        for (int b = 0; b < size_b; ++b) {
          T2ab_blk_start = b;                // T^ab_ij
          T2ba_blk_start = b * size_a;       // T^ba_ij

          Rapb_blk_start = b;                // R^ij_a'b
          Rbap_blk_start = b * size_ap;      // R^ij_ba'
          const double* iter_F12_blk = F12_blk + *F12_blk_start;

          double* iter_RT2 = RT2;
          for(int ap = 0; ap < size_ap; ++ap) {

            const double* iter_T2_blk = T2_blk + *T2_blk_start;
            for (int a = 0; a < size_a; ++a) {
              // eg. R^ij_a'b * T^ab_ij (sum over b)
              const double RTapa = (*iter_F12_blk) * (*iter_T2_blk);
              *iter_RT2 += RTapa;
              ++iter_RT2;

              iter_T2_blk += *T2_blk_offset;    // T^ab_ij: size_b; T^ba_ij: one
            } // end of a loop

            iter_F12_blk += *F12_blk_offset;   // R^ij_a'b: size_b; R^ij_ba': one
          } // end of a' loop
        } // end of b loop

        F12_ints->release_pair_block(i, j, f12_idx);
        T2_ints->release_pair_block(i, j, 0);
      } // end of j loop
    } // end of i loop

  }
  // end of compute_RT2_sum_3idx

  // compute R * T2 (mp2 amplitudes) contracted over 3 indices (eg: R^ij_a'b * T^ab_ij)
  // here they are: i, j, idx3
  void compute_RTmp2_sum_3idx(const int RTb1b2_k1k2, const unsigned int f12_idx,
                              const int size_idx2,
                              const Ref<DistArray4>& f12_ints, const double* const T2_ints,
                              double* const RT2)
  {
    // idx1 and idx2 are label indices, idx3 is the other contracted index
    int idx1 = 0;    // a'
    int idx2 = 0;    // a
    int idx3 = 0;    // b

    int size_idx1 = f12_ints->ni();
    int size_idx3 = f12_ints->nj();
    const blasint nij = f12_ints->nx() * f12_ints->ny();

    // get block indices
    int* f12_blk_idx1 = &idx1;
    int* f12_blk_idx2 = &idx3;

    // Get the right start iterator for T^ab_ij or T^ba_ij
    int iter_T2ab_start = 0;
    int iter_T2ba_start = 0;
    int* iter_T2_start = &iter_T2ab_start;

    // Get the right offset for T^ab_ij or T^ba_ij in the inner loop
    int offset_T2ab = nij;
    int offset_T2ba = size_idx2 * nij;
    int* offset_T2 = &offset_T2ab;

    switch (RTb1b2_k1k2) {
    case RT13_23:
      // default value
    break;

    case RT31_23:
      size_idx1 = f12_ints->nj();
      size_idx3 = f12_ints->ni();
      f12_blk_idx1 = &idx3;
      f12_blk_idx2 = &idx1;
    break;

    case RT13_32:
      iter_T2_start = &iter_T2ba_start;
      offset_T2 = &offset_T2ba;
    break;

    case RT31_32:
      size_idx1 = f12_ints->nj();
      size_idx3 = f12_ints->ni();
      f12_blk_idx1 = &idx3;
      f12_blk_idx2 = &idx1;
      iter_T2_start = &iter_T2ba_start;
      offset_T2 = &offset_T2ba;
    break;

    default:
      ExEnv::out0() << "There is no such index case";
      break;
    }

    const blasint one = 1; // for F77_DDOT
    double* iter_RT2 = RT2;

    // index 1: a', index 2: a, index 3: b
      for (idx1 = 0; idx1 < size_idx1; ++idx1) {   // a'

        for(idx2 = 0; idx2 < size_idx2; ++idx2) {  // a
          iter_T2ab_start = idx2 * size_idx3 * nij;
          iter_T2ba_start = idx2 * nij;
          const double* iter_T2 = T2_ints + *iter_T2_start;

          double f12T2_sum_idx3 = 0;
          for(idx3 = 0; idx3 < size_idx3; ++idx3) {   // b
            const double* f12_ij_blk = f12_ints->retrieve_pair_block((*f12_blk_idx1), (*f12_blk_idx2), f12_idx);
            const double f12T2_sum_ij = F77_DDOT(&nij, f12_ij_blk, &one, iter_T2, &one);
            f12T2_sum_idx3 += f12T2_sum_ij;

            iter_T2 += *offset_T2;
            f12_ints->release_pair_block((*f12_blk_idx1), (*f12_blk_idx2), f12_idx);
          }
//          ExEnv::out0() << indent << scprintf("%12.10f", f12T2_sum_idj) << " ";
          *iter_RT2 = f12T2_sum_idx3;
          ++iter_RT2;
        }
//        ExEnv::out0() << endl;
      }
  }
  // end of compute_RTmp2_sum_3idx

  // print out integrals for test propose
  void print_f12_ints(const string& spinlabel, const string& label,
                      const int f12_idx, Ref<DistArray4>& f12_ints)
  {
    const int size_idx1 = f12_ints->ni();
    const int size_idx2 = f12_ints->nj();
    const int size_idx3 = f12_ints->nx();
    const int size_idx4 = f12_ints->ny();

    ExEnv::out0() << indent << spinlabel << " " << label << endl;
    ExEnv::out0() << indent << "size idx1: " << size_idx1
                           <<"  size idx2: "<< size_idx2 << endl;
    for (int idx1 = 0; idx1 < size_idx1; ++idx1) {
      for(int idx2 = 0; idx2 < size_idx2; ++idx2) {
//        ExEnv::out0() << indent << "index1: " << idx1 <<"  index2:"<< idx2 << endl;

        const double* f12_ij_blk = f12_ints->retrieve_pair_block(idx1, idx2, f12_idx);
        for (int idx3 = 0; idx3 < size_idx3; ++idx3) {
          for(int idx4 = 0; idx4 < size_idx4; ++idx4) {

            ExEnv::out0() << indent << scprintf("%12.10f", *f12_ij_blk) << " ";
            ++f12_ij_blk;
          }
          ExEnv::out0() << " * ";
        }
        f12_ints->release_pair_block(idx1, idx2, f12_idx);
        ExEnv::out0() << endl;
      }
      ExEnv::out0() << endl;
    }
  } // end of print_ints

  void print_T2ijab_mp2(const string& spinlabel, const string& label,
                        const int no1, const int no2,
                        const int nv1, const int nv2,
                        const double* const T2)
  {
    const int nv12 = nv1 * nv2;
    const int no2v12 = no2 * nv12;

    ExEnv::out0() << indent << spinlabel << " " << label << endl
                  << indent << "a b" << "    i j" << endl;

    for (int a = 0 ; a < nv1; ++a) {
      for (int b = 0; b < nv2; ++b) {
        ExEnv::out0() << indent << a << " " << b << ":  " ;

        for (int i = 0 ; i < no1; ++i) {

          const double* iter_T2 = T2 + i * no2v12 + a * nv2 + b;
          for (int j = 0; j < no2; ++j) {
            ExEnv::out0() << indent << scprintf("%12.10f", *iter_T2) << "  ";

            iter_T2 += nv12;
          }
        }
        ExEnv::out0() << endl;
      }

      ExEnv::out0() << endl;
    }
  }

  void print_T2abij_mp2(const string& spinlabel, const string& label,
                        const int no1, const int no2,
                        const int nv1, const int nv2,
                        const double* const T2)
  {
     ExEnv::out0() << indent << spinlabel << " " << label << endl
                   << indent << "a b" << "    i j" << endl;

     const double* iter_T2 = T2;
     for (int a = 0 ; a < nv1; ++a) {
       for (int b = 0; b < nv2; ++b) {
         ExEnv::out0() << indent << a << " " << b << ":  ";

         for (int i = 0 ; i < no1; ++i) {
           for (int j = 0; j < no2; ++j) {
             ExEnv::out0() << indent << scprintf("%12.10f", *iter_T2) << "  ";

             ++iter_T2;
           }
         }
         ExEnv::out0() << endl;
       }
       ExEnv::out0() << endl;
     }
  }

  // print ints^b1b2_k1k2
  void print_ints(const string& label, Ref<DistArray4>& ints,
                  const unsigned int ints_idx)
  {
    const int size_b1 = ints->ni();
    const int size_b2 = ints->nj();
    const int size_k1 = ints->nx();
    const int size_k2 = ints->ny();

    ExEnv::out0() << label << endl;

    for(int idx_b1 = 0; idx_b1 < size_b1; ++idx_b1) {
      for (int idx_b2 = 0; idx_b2 < size_b2; ++idx_b2) {
        ExEnv::out0() << indent << "idx1 = " << idx_b1 << " idx2 = " << idx_b2 << ": ";

        const double* ints_blk = ints->retrieve_pair_block(idx_b1, idx_b2, ints_idx);
        print_intermediate("", "", ints_blk, size_k1, size_k2);

        ints->release_pair_block(idx_b1, idx_b2, ints_idx);
      }
      ExEnv::out0() << endl;
    }
  }

  // print R^ij_a'b or R^ij_ba' or R^ji_a'b or  R^ji_ba'
  void print_f12_ints(const string& label, const SpinCase1 spin,
                      Ref<DistArray4>& f12_ints, const unsigned int f12_idx,
                      const int nap, const int nb, const int ni, const int nj)
  {
    int ap = 0;
    int b = 0;
    int* f12_blk_idx1 = &ap;
    int* f12_blk_idx2 = &b;
    if (spin == Beta) {
      f12_blk_idx1 = &b;
      f12_blk_idx2 = &ap;
    }

    ExEnv::out0() << label << endl;

    for(int ap = 0; ap < nap; ++ap) {
      for (int b = 0; b < nb; ++b) {
        ExEnv::out0() << indent << "a' = " << ap << " b = " << b << ": ";

        const double* f12_blk = f12_ints->retrieve_pair_block(*f12_blk_idx1, *f12_blk_idx2, f12_idx);
        const double* iter_f12_blk = f12_blk;

        for (int i = 0 ; i < ni; ++i) {
          for (int j = 0 ; j < nj; ++j) {
          ExEnv::out0() << indent << scprintf("%12.10f", *iter_f12_blk) << " ";

          ++iter_f12_blk;
          }
        }
        ExEnv::out0() << "**";
        f12_ints->release_pair_block(*f12_blk_idx1, *f12_blk_idx2, f12_idx);
      }
      ExEnv::out0() << endl;
    }
  }

  // compute D^a'_a matrix
  void compute_Dapa_spin(const SpinCase1 spin, const double C_0, const double C_1,
                         const int nocc_alpha, const int nocc_beta, const int ncabs,
                         const double* RT2ij, const double* RT2ji,
                         const double* RT2ij_ab, const double* RT2ji_ab,
                         RefSCMatrix Dapa)
  {
    const int nvir = (spin == Alpha? nocc_alpha : nocc_beta);
    // determine the size of index b for the calculation of AlphaBeta contribution
    const int* size_b = &nocc_beta;
    if (spin == Beta) size_b = &nocc_alpha;

    for (int ap = 0; ap < ncabs; ++ap) {
      for (int a = 0; a < nvir; ++a) {

        // AlphaAlpha/BetaBeta part
        const int nab = nvir * nvir;
        const int iter_start = ap * nab + a;
        const double* iter_RT2ij = RT2ij + iter_start;
        const double* iter_RT2ji = RT2ji + iter_start;
        double rt2ij = 0.0;
        double rt2ji = 0.0;

        double dap_a = 0.0;
        for (int b = 0; b < nvir; ++b) {
          rt2ij += *iter_RT2ij;
          rt2ji += *iter_RT2ji;

          iter_RT2ij += nvir;
          iter_RT2ji += nvir;
        }
        dap_a += 0.5 * C_1 * (rt2ij - rt2ji);
        ExEnv::out0() << "spincase1 part d^" << ap
                    << "_" << a << " = " << scprintf("%12.10f", dap_a) << endl;

        // AlphaBeta part
        if (nocc_alpha != 0 && nocc_beta != 0) {
          const int nab = nvir * (*size_b);
          const int iter_start_ab = ap * nab + a;
          const double* iter_RT2ij_ab = RT2ij_ab + iter_start;
          const double* iter_RT2ji_ab = RT2ji_ab + iter_start;

          double rt2ij_ab = 0;
          double rt2ji_ab = 0;
          for (int b = 0; b < (*size_b); ++b) {
            rt2ij_ab += *iter_RT2ij_ab;
            rt2ji_ab += *iter_RT2ji_ab;

            iter_RT2ij_ab += nvir;
            iter_RT2ji_ab += nvir;
          }
          dap_a += 0.5 * (C_0 + C_1) * rt2ij_ab
                 + 0.5 * (C_0 - C_1) * rt2ji_ab;
        }
        ExEnv::out0() << "d^"  << ap << "_" << a << " = " << scprintf("%12.10f", dap_a) << endl;
        Dapa.set_element(ap, a, dap_a);
      }
    }
  } // end of function: compute_Dapa

  // transform T2(i,j,a,b) to T2(a,b,i,j)
  void transform_T2ijab_T2abij (const Ref<DistArray4>& T2ij_ab, double* const T2ab_ij )
  {
    const int size_i = T2ij_ab->ni();
    const int size_j = T2ij_ab->nj();
    const int size_a = T2ij_ab->nx();
    const int size_b = T2ij_ab->ny();

    for (int i = 0; i != size_i; ++i) {
      for (int j = 0; j != size_j; ++j) {
        const double* const T2_ab_blk = T2ij_ab->retrieve_pair_block(i, j, 0);

//        const double* iter_T2_ab_blk = T2_ab_blk;
        for (int a = 0; a != size_a; ++a) {
          for (int b = 0; b != size_b; ++b) {
            const int T2abij_idx = ((a * size_b + b) * size_i + i) * size_j + j;
            const int T2_ab_blk_idx = a * size_b + b;
            T2ab_ij[T2abij_idx] = T2_ab_blk[T2_ab_blk_idx];

//              ++iter_T2_ab_blk;
          }
        }

      } // end of j loop
   } // end of i loop

  } // end of transform_T2ijab_T2abij function

} // end of namespace

// obtain vectors of spin1 and spin2 orbitals
// in order of occ_act, vir, orbs, cabs, and occ orbitals
void MP2R12Energy_Diag::obtain_orbitals(const SpinCase2 spincase,
                                        vector<Ref<OrbitalSpace> >& v_orbs1,
                                        vector<Ref<OrbitalSpace> >& v_orbs2)
{
  const SpinCase1 spin1 = case1(spincase);
  const SpinCase1 spin2 = case2(spincase);

  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();

  const Ref<OrbitalSpace> occ1_act = r12eval()->occ_act(spin1);
  const Ref<OrbitalSpace> vir1 = r12eval()->vir(spin1);
  const Ref<OrbitalSpace> orbs1 = r12eval()->orbs(spin1);
  const Ref<OrbitalSpace> cabs1 = r12world->cabs_space(spin1);
  // use canonical cabs for test propose
  //const Ref<OrbitalSpace> cabs1 = r12eval()->cabs_space_hcanonical(spin1);

  const Ref<OrbitalSpace> occ1 = r12eval()->occ(spin1);

  const Ref<OrbitalSpace> occ2_act = r12eval()->occ_act(spin2);
  const Ref<OrbitalSpace> vir2 = r12eval()->vir(spin2);
  const Ref<OrbitalSpace> orbs2 = r12eval()->orbs(spin2);
  const Ref<OrbitalSpace> cabs2 = r12world->cabs_space(spin2);
  // use canonical cabs for test propose
  //const Ref<OrbitalSpace> cabs2 = r12eval()->cabs_space_hcanonical(spin2);
  const Ref<OrbitalSpace> occ2 = r12eval()->occ(spin2);

  v_orbs1.push_back(occ1_act);
  v_orbs1.push_back(vir1);
  v_orbs1.push_back(orbs1);
  v_orbs1.push_back(cabs1);
  v_orbs1.push_back(occ1);

  v_orbs2.push_back(occ2_act);
  v_orbs2.push_back(vir2);
  v_orbs2.push_back(orbs2);
  v_orbs2.push_back(cabs2);
  v_orbs2.push_back(occ2);
}
// end of obtain_orbitals

// activate the three f12_ints for X
void MP2R12Energy_Diag::activate_ints_X_f12(Ref<TwoBodyFourCenterMOIntsRuntime>& moints4_rtime, const string& index,
                                            const vector< Ref<OrbitalSpace> >& v_orbs1,
                                            const vector< Ref<OrbitalSpace> >& v_orbs2,
                                            const string& descr_f12_key, vector<Ref<DistArray4> >& f12_ints)
{
  // obtain orbitals
  const Ref<OrbitalSpace> occ1_act = v_orbs1[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2[0];
  const Ref<OrbitalSpace> vir1 = v_orbs1[1];
  const Ref<OrbitalSpace> vir2 = v_orbs2[1];
  const Ref<OrbitalSpace> orbs1 = v_orbs1[2];
  const Ref<OrbitalSpace> orbs2 = v_orbs2[2];
  const Ref<OrbitalSpace> cabs1 = v_orbs1[3];
  const Ref<OrbitalSpace> cabs2 = v_orbs2[3];
  const Ref<OrbitalSpace> occ1 = v_orbs1[4];
  const Ref<OrbitalSpace> occ2 = v_orbs2[4];

  // indices for (f12)^IJ or (f12)_IJ
  const Ref<OrbitalSpace>* bk1 = &occ1_act;
  const Ref<OrbitalSpace>* bk2 = &occ2_act;

 if (index == "ji") {
    bk1 = &occ2_act;
    bk2 = &occ1_act;
 }
  // f12^b1b2_PQ or f12^PQ_k1k2
  Ref<DistArray4> p1p2_ints;
  activate_ints((*bk1)->id(), (*bk2)->id(), orbs1->id(), orbs2->id(),
                descr_f12_key, moints4_rtime, p1p2_ints);
  f12_ints.push_back(p1p2_ints);

  // f12^b1b2_MA' or f12^MA'_k1k2
  Ref<DistArray4> i1a2_ints;
  activate_ints((*bk1)->id(), (*bk2)->id(), occ1->id(), cabs2->id(),
                descr_f12_key, moints4_rtime, i1a2_ints);
  f12_ints.push_back(i1a2_ints);

  // f12^b1b2_A'M or f12^A'M_k1k2
  Ref<DistArray4> a1i2_ints;
  activate_ints((*bk1)->id(), (*bk2)->id(), cabs1->id(), occ2->id(),
                descr_f12_key, moints4_rtime, a1i2_ints);
  f12_ints.push_back(a1i2_ints);
}
// end of function: activate_ints_X_f12

// functions needed for compute_Dii_test:
// 1) compute AlphaAlpha/BetaBeta: R^ij_ab R^ab_ij & R^ji_ab R^ab_ij
void MP2R12Energy_Diag::compute_RRii_ii(vector<string>& output,
                                        const vector< Ref<OrbitalSpace> >& v_orbs1,
                                        const vector< Ref<OrbitalSpace> >& v_orbs2,
                                        double* RRij_ij, double* RRji_ij)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  Ref<TwoBodyIntDescr> descr_f12f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(), 0, 0);
  const string descr_f12f12_key = moints4_rtime->descr_key(descr_f12f12);

  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);
  const TwoBodyOper::type f12f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12f12();
  const unsigned int f12f12_idx = descr_f12f12->intset(f12f12_type);

  // activate_ints
  Ref<DistArray4> ijij_f12f12_ints;
  vector<Ref<DistArray4> > ij_f12_ints;

  const Ref<OrbitalSpace> occ1_act = v_orbs1[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2[0];
  activate_ints(occ1_act->id(), occ2_act->id(), occ1_act->id(), occ2_act->id(),
                descr_f12f12_key, moints4_rtime, ijij_f12f12_ints);

  activate_ints_X_f12(moints4_rtime, "ij",
                      v_orbs1, v_orbs2,
                      descr_f12_key, ij_f12_ints);
  // compute X^ij_ij
  compute_VX(ij_ij, output, f12f12_idx, ijij_f12f12_ints, f12_idx,
             f12_idx, ij_f12_ints, ij_f12_ints, RRij_ij);
  // compute X^ji_ij
  // AlphaAlpha or BetaBeta: ji_f12_ints = ij_f12_ints
  compute_VX(ji_ij, output, f12f12_idx, ijij_f12f12_ints, f12_idx,
             f12_idx, ij_f12_ints, ij_f12_ints, RRji_ij);

  // deactivate integrals
  ijij_f12f12_ints->deactivate();
  for (int i=0; i<ij_f12_ints.size(); ++i) {
    ij_f12_ints[i]->deactivate();
  }

}
// end of function: compute_RRii_ii

// 2) compute the openshell AlphaBeta:
// R^ij_ab R^ab_ij, R^ji_ab R^ab_ij, R^ij_ab R^ab_ji, and R^ji_ab R^ab_ji
void MP2R12Energy_Diag::compute_Rii_ii(vector<string>& output,
                                       const vector< Ref<OrbitalSpace> >& v_orbs1,
                                       const vector< Ref<OrbitalSpace> >& v_orbs2,
                                       double* RRi1i2_i1i2, double* RRi1i2_i2i1,
                                       double* RRi2i1_i1i2, double* RRi2i1_i2i1)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  Ref<TwoBodyIntDescr> descr_f12f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(), 0, 0);
  const string descr_f12f12_key = moints4_rtime->descr_key(descr_f12f12);

  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);
  const TwoBodyOper::type f12f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12f12();
  const unsigned int f12f12_idx = descr_f12f12->intset(f12f12_type);

  // activate_ints
  // f12f12_ints
  const Ref<OrbitalSpace> occ1_act = v_orbs1[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2[0];

  Ref<DistArray4> ijij_f12f12_ints;
  Ref<DistArray4> jiij_f12f12_ints;
  Ref<DistArray4> ijji_f12f12_ints;
  Ref<DistArray4> jiji_f12f12_ints;

  activate_ints(occ1_act->id(), occ2_act->id(), occ1_act->id(), occ2_act->id(),
                descr_f12f12_key, moints4_rtime, ijij_f12f12_ints);

  activate_ints(occ2_act->id(), occ1_act->id(), occ1_act->id(), occ2_act->id(),
                descr_f12f12_key, moints4_rtime, jiij_f12f12_ints);

  activate_ints(occ1_act->id(), occ2_act->id(), occ2_act->id(), occ1_act->id(),
                descr_f12f12_key, moints4_rtime, ijji_f12f12_ints);

  activate_ints(occ2_act->id(), occ1_act->id(), occ2_act->id(), occ1_act->id(),
                descr_f12f12_key, moints4_rtime, jiji_f12f12_ints);

  // f12_ints
  vector<Ref<DistArray4> > ij_f12_ints;
  vector<Ref<DistArray4> > ji_f12_ints;
  activate_ints_X_f12(moints4_rtime, "ij",
                      v_orbs1, v_orbs2,
                      descr_f12_key, ij_f12_ints);
  activate_ints_X_f12(moints4_rtime, "ji",
                      v_orbs1, v_orbs2,
                      descr_f12_key, ji_f12_ints);
  // X^IJ_IJ
  compute_VX(ij_ij, output, f12f12_idx, ijij_f12f12_ints, f12_idx,
             f12_idx, ij_f12_ints, ij_f12_ints, RRi1i2_i1i2);
  // X^JI_IJ
  compute_VX(ji_ij, output, f12f12_idx, jiij_f12f12_ints, f12_idx,
             f12_idx, ji_f12_ints, ij_f12_ints, RRi2i1_i1i2);
  // X^IJ_JI
  compute_VX(ij_ji, output, f12f12_idx, ijji_f12f12_ints, f12_idx,
             f12_idx, ij_f12_ints, ji_f12_ints, RRi1i2_i2i1);
  // X^JI_JI
  compute_VX(ji_ji, output, f12f12_idx, jiji_f12f12_ints, f12_idx,
             f12_idx, ji_f12_ints, ji_f12_ints, RRi2i1_i2i1);

  // deactivate integrals
  ijij_f12f12_ints->deactivate();
  jiij_f12f12_ints->deactivate();
  ijji_f12f12_ints->deactivate();
  jiji_f12f12_ints->deactivate();
  for (int i=0; i<ij_f12_ints.size(); ++i) {
    ij_f12_ints[i]->deactivate();
    ji_f12_ints[i]->deactivate();
  }
}
// end of function: compute_Rii_ii

// 3) function for testing D^i_i
void MP2R12Energy_Diag::compute_Dii_test(const int nspincases1, const int nspincases2,
                                         const int nocc_alpha, const int nocc_beta,
                                        const double C_0, const double C_1)
{
  // output used for printing the intermediates in R * R (X) computation
  vector<string> output;
  output.push_back("diag-pq contribution");
  output.push_back("diag-pq-ma' contribution");
  output.push_back("diag-pq-ma'-a'm contribution");

  // AlphaBeta R * R (X)
  double* RRi1i2_i1i2 = 0;
  double* RRi1i2_i2i1 = 0;
  double* RRi2i1_i1i2 = 0;
  double* RRi2i1_i2i1 = 0;
  bool ab_RRii_ii_evaluated = false;

  for(int s = 0; s < nspincases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);

    // spin=0 (Alpha) => AlphaAlpha case (1)
    // spin=1 (Beta) => BetaBeta case (2)
    const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
    string spinletters = to_string(spincase);

    // obtain occ, vir, orbs, and cabs orbitals in that order
    vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
    vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
    obtain_orbitals(spincase, v_orbs1, v_orbs2);

    const int nocc_act = v_orbs1[0]->rank();
    if (debug_ >= DefaultPrintThresholds::N2)
      ExEnv::out0() << endl << spinletters << " D^i_i (test):" << endl
                    << "number of occupied orbital: " << nocc_act << endl;

    // skip spincase if no electron of this kind
    if (nocc_act == 0)
      continue;
    const int nocc12 = nocc_act * nocc_act;

    // calculate AlphaAlpha/BetaBeta part for D^i_i:
    // R^ij_ab R^ab_ij & R^ji_ab R^ab_ij contracted over ab
    double* RRij_ij = new double[nocc12];
    double* RRji_ij = new double[nocc12];
    fill_n(RRij_ij, nocc12, 0.0);
    fill_n(RRji_ij, nocc12, 0.0);

    compute_RRii_ii(output,
                    v_orbs1, v_orbs2,
                    RRij_ij, RRji_ij);
    if (debug_ >= DefaultPrintThresholds::N2) {
      print_intermediate(spinletters, "RRij_ij", RRij_ij, nocc_act, nocc_act);
      print_intermediate(spinletters, "RRji_ij", RRji_ij, nocc_act, nocc_act);
    }

    // calculate AlphaBeta part for D^i_i:
    if (nocc_alpha != 0 && nocc_beta != 0) {

      if (nspincases2 == 3 && !ab_RRii_ii_evaluated) {
        // only need to compute once, same intermediates will be used for Beta D^i_i

        // obtain occ, vir, orbs, and cabs orbitals
        vector< Ref<OrbitalSpace> > v_orbs1;
        vector< Ref<OrbitalSpace> > v_orbs2;
        obtain_orbitals(AlphaBeta, v_orbs1, v_orbs2);

        // R^ij_ab * R^ab_ij, R^ji_ab * R^ab_ij, R^ij_ab * R^ab_ji, R^ji_ab * R^ab_ji
        const int nocc1_act = v_orbs1[0]->rank();
        const int nocc2_act = v_orbs2[0]->rank();
        const int nocc12 = nocc1_act* nocc2_act;
        if (debug_ >= DefaultPrintThresholds::N2)
          ExEnv::out0() << endl << "AlphaBeta D^i_i (test)"
                        << endl << "Number of alpha occupied orbital: " << nocc1_act
                        << endl << "Number of beta occupied orbital: " << nocc2_act << endl;

        RRi1i2_i1i2 = new double[nocc12];
        RRi1i2_i2i1 = new double[nocc12];
        RRi2i1_i1i2 = new double[nocc12];
        RRi2i1_i2i1 = new double[nocc12];
        fill_n(RRi1i2_i1i2, nocc12, 0.0);
        fill_n(RRi1i2_i2i1, nocc12, 0.0);
        fill_n(RRi2i1_i1i2, nocc12, 0.0);
        fill_n(RRi2i1_i2i1, nocc12, 0.0);

        compute_Rii_ii(output,
                       v_orbs1, v_orbs2,
                       RRi1i2_i1i2, RRi1i2_i2i1,
                       RRi2i1_i1i2, RRi2i1_i2i1);
        ab_RRii_ii_evaluated = true;

        if (debug_ >= DefaultPrintThresholds::N2) {
          print_intermediate("AlphaBeta", "RRi1i2_i1i2", RRi1i2_i1i2, nocc_alpha, nocc_beta);
          print_intermediate("AlphaBeta", "RRi1i2_i2i1", RRi1i2_i2i1, nocc_alpha, nocc_beta);
          print_intermediate("AlphaBeta", "RRi2i1_i1i2", RRi2i1_i1i2, nocc_alpha, nocc_beta);
          print_intermediate("AlphaBeta", "RRi2i1_i2i1", RRi2i1_i2i1, nocc_alpha, nocc_beta);
        }

      } else if (nspincases2 == 2) {
          RRi1i2_i1i2 = RRij_ij;
          RRi1i2_i2i1 = RRji_ij;
          RRi2i1_i1i2 = RRji_ij;
          RRi2i1_i2i1 = RRij_ij;
     }

    }// end of AlphaBeta part for Di_i

    // calculate D^i_i
    compute_Dii_spin(spin, C_0, C_1,
                     nocc_alpha, nocc_beta,
                     RRij_ij, RRji_ij,
                     RRi1i2_i1i2, RRi1i2_i2i1, RRi2i1_i1i2, RRi2i1_i2i1);

     delete[] RRij_ij;
     delete[] RRji_ij;
  } // end of spincase1 loop

  if (nspincases2 == 3) {
    delete[] RRi1i2_i1i2;
    delete[] RRi1i2_i2i1;
    delete[] RRi2i1_i1i2;
    delete[] RRi2i1_i2i1;
  }
}
// end of compute_Dii_test

// compute R^ab_b1b2 R^k1k2 _ab (a, b represent complete virtual orbitals)
// sum over a, b, and index 3
// index 1,2: i, m; index 3: j  e.g.: R^ab_ij R^jm_ab
void MP2R12Energy_Diag::compute_RR_sum_abj2(const int RRb1b2_k1k2,
                                            const int f12f12_idx, const int f12_idx,
                                            const Ref<DistArray4>& f12f12_ints,
                                            const vector< Ref<DistArray4> >& v_f12_ints1,
                                            const vector< Ref<DistArray4> >& v_f12_ints2,
                                            double* const RR_result)
{
  int nocc_act = f12f12_ints->ni();
  if (RRb1b2_k1k2 == RR31_23 || RRb1b2_k1k2 == RR31_32) {
    nocc_act = f12f12_ints->nj();
  }
  const int nocc12 = nocc_act * nocc_act;

  // R^ab_b1b2 R^k1k2 _ab = (f12f12)^b1b2_k1k2
  //                       - f12^b1b2_pq f12^pq_k1k2
  //                       - f12^b1b2_pq_ma' f12^ma'_k1k2
  //                       - f12^b1b2_pq_a'm f12^a'm_k1k2

  // add (f12f12)^b1b2_k1k2
  compute_F12F12_sum_idx3_2(RRb1b2_k1k2,
                            f12f12_idx, f12f12_ints,
                            RR_result);
//  print_intermediate("test", "f12f12 part", RR_result, nocc_act, nocc_act);

  if (v_f12_ints1.size() != 3)
    ExEnv::out0() << "Then number of integrals in computing R^ab_b1b2 R^k1k2 _ab are wrong"
                  << endl;

  // subtract f12^b1b2_pq f12^pq_k1k2, f12^b1b2_pq_ma' f12^ma'_k1k2,
  // and f12^b1b2_pq_a'm f12^a'm_k1k2
  for (int i = 0; i != 3; ++i) {
    const Ref<DistArray4> f12_ints1 = v_f12_ints1[i];
    const Ref<DistArray4> f12_ints2 = v_f12_ints2[i];

    double* RR_part = new double[nocc12];
    fill_n(RR_part, nocc12, 0.0);
    compute_RR2_sum_3idx(RRb1b2_k1k2, f12_idx,
                         f12_ints1, f12_ints2,
                         RR_part);
//    print_intermediate("test", "f12 part", RR_part, nocc_act, nocc_act);

    double* iter_RR_result = RR_result;
    const double* iter_RR_part = RR_part;
    for (int j = 0; j != nocc_act; ++j) {
      for (int m = 0; m != nocc_act; ++m) {
        *iter_RR_result -= *iter_RR_part;
//        ExEnv::out0() << scprintf("%12.10f", *iter_RR_result) << " ";
        ++iter_RR_result;
        ++iter_RR_part;
      }
//      ExEnv::out0() << endl;
    }
//    ExEnv::out0() << endl;

    delete[] RR_part;
    RR_part = NULL;
//    print_intermediate("test", "plus f12 part", RR_result, nocc_act, nocc_act);
  }

}
// end of compute_RR_sum_abj_2

// compute D^m_i
void MP2R12Energy_Diag::compute_Dmi(const int nspincases1, const int nspincases2,
                                    const double C_0, const double C_1,
                                    const vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                                    const vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                                    double* const Dm_i_alpha, double* const Dm_i_beta)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  Ref<TwoBodyIntDescr> descr_f12f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(), 0, 0);
  const string descr_f12f12_key = moints4_rtime->descr_key(descr_f12f12);

  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);
  const TwoBodyOper::type f12f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12f12();
  const unsigned int f12f12_idx = descr_f12f12->intset(f12f12_type);

  // preliminaries for AlphaBeta parts
  const Ref<OrbitalSpace> occ1_act = v_orbs1_ab[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2_ab[0];
  const int nocc_alpha = occ1_act->rank();
  const int nocc_beta = occ2_act->rank();
  // test
  if (debug_ >= DefaultPrintThresholds::N2) {
    ExEnv::out0() << endl << "Number of alpha occupied orbital: " << nocc_alpha
                  << endl << "Number of beta occupied orbital: " << nocc_beta << endl;
  }

  // activate integrals for R^i1i2_a1a2 R^a1a2_i1i2, R^i2i1_a1a2 R^a1a2_i1i2
  //                        R^i1i2_a1a2 R^a1a2_i2i1, R^i2i2_a1a2 R^a1a2_i2i1
  // i: occupied orbital, a: complete virtual orbitals
  // 1: alpha orbital, 2: beta orbital

  // R^b1b2_ab R^ab_k1k2 = (f12f12)^b1b2_k1k2
  //                     - f12^b1b2_pq f12^pq_k1k2
  //                     - f12^b1b2_pq_ma' f12^ma'_k1k2
  //                     - f12^b1b2_pq_a'm f12^a'm_k1k2
  // which means (f12f12) and (f12) ints are needed

  // (f12f12)^ii_ii integrals
  Ref<DistArray4> i1i2_f12f12_i1i2;
  Ref<DistArray4> i2i1_f12f12_i1i2;
  Ref<DistArray4> i1i2_f12f12_i2i1;
  Ref<DistArray4> i2i1_f12f12_i2i1;

  // vectors for (f12)^ii_k1k2 integrals
  // b1 and b2: occupied orbitals
  vector<Ref<DistArray4> > f12_ints_i1i2;
  vector<Ref<DistArray4> > f12_ints_i2i1;

  if (nspincases2 == 3) {
    // activate (f12f12) ints
     activate_ints(occ1_act->id(), occ2_act->id(), occ1_act->id(), occ2_act->id(),
                   descr_f12f12_key, moints4_rtime, i1i2_f12f12_i1i2);

     activate_ints(occ2_act->id(), occ1_act->id(), occ1_act->id(), occ2_act->id(),
                   descr_f12f12_key, moints4_rtime, i2i1_f12f12_i1i2);

     activate_ints(occ1_act->id(), occ2_act->id(), occ2_act->id(), occ1_act->id(),
                   descr_f12f12_key, moints4_rtime, i1i2_f12f12_i2i1);

     activate_ints(occ2_act->id(), occ1_act->id(), occ2_act->id(), occ1_act->id(),
                   descr_f12f12_key, moints4_rtime, i2i1_f12f12_i2i1);
     // activate f12 ints
     activate_ints_X_f12(moints4_rtime, "ij",
                         v_orbs1_ab, v_orbs2_ab,
                         descr_f12_key, f12_ints_i1i2);
     activate_ints_X_f12(moints4_rtime, "ji",
                         v_orbs1_ab, v_orbs2_ab,
                         descr_f12_key, f12_ints_i2i1);
     // test for f12 ints
#if 0
     {
       // obtain orbitals
       const Ref<OrbitalSpace> occ1_act = v_orbs1_ab[0];
       const Ref<OrbitalSpace> occ2_act = v_orbs2_ab[0];
       const Ref<OrbitalSpace> vir1 = v_orbs1_ab[1];
       const Ref<OrbitalSpace> vir2 = v_orbs2_ab[1];
       const Ref<OrbitalSpace> orbs1 = v_orbs1_ab[2];
       const Ref<OrbitalSpace> orbs2 = v_orbs2_ab[2];
       const Ref<OrbitalSpace> cabs1 = v_orbs1_ab[3];
       const Ref<OrbitalSpace> cabs2 = v_orbs2_ab[3];
       const Ref<OrbitalSpace> occ1 = v_orbs1_ab[4];
       const Ref<OrbitalSpace> occ2 = v_orbs2_ab[4];

       // f12^b1b2_PQ
       Ref<DistArray4> p1p2_ints1;
       activate_ints(occ1_act->id(), occ2_act->id(), orbs1->id(), orbs2->id(),
                     descr_f12_key, moints4_rtime, p1p2_ints1);
       f12_ints_i1i2.push_back(p1p2_ints1);

       // f12^b1b2_MA'
       Ref<DistArray4> i1a2_ints1;
       activate_ints(occ1_act->id(), occ2_act->id(), occ1->id(), cabs2->id(),
                     descr_f12_key, moints4_rtime, i1a2_ints1);
       f12_ints_i1i2.push_back(i1a2_ints1);

       // f12^b1b2_A'M
       Ref<DistArray4> a1i2_ints1;
       activate_ints(occ1_act->id(), occ2_act->id(), cabs1->id(), occ2->id(),
                     descr_f12_key, moints4_rtime, a1i2_ints1);
       f12_ints_i1i2.push_back(a1i2_ints1);


       // f12^PQ_k1k2
       Ref<DistArray4> p1p2_ints2;
       activate_ints(occ2_act->id(), occ1_act->id(), orbs1->id(), orbs2->id(),
                     descr_f12_key, moints4_rtime, p1p2_ints2);
       f12_ints_i2i1.push_back(p1p2_ints2);

       // f12^MA'_k1k2
       Ref<DistArray4> i1a2_ints2;
       activate_ints(occ2_act->id(), occ1_act->id(), occ1->id(), cabs2->id(),
                     descr_f12_key, moints4_rtime, i1a2_ints2);
       f12_ints_i2i1.push_back(i1a2_ints2);

       // f12^A'M_k1k2
       Ref<DistArray4> a1i2_ints2;
       activate_ints(occ2_act->id(), occ1_act->id(), cabs1->id(), occ2->id(),
                     descr_f12_key, moints4_rtime, a1i2_ints2);
       f12_ints_i2i1.push_back(a1i2_ints2);
     }
#endif
  } // end of activating ints for alphabeta part

  for(int s = 0; s < nspincases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);

    // spin=0 (Alpha) => AlphaAlpha case (1)
    // spin=1 (Beta) => BetaBeta case (2)
    const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
    string spinletters = to_string(spincase);

    // obtain occ, vir, orbs, and cabs orbitals in that order
    vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
    vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
    obtain_orbitals(spincase, v_orbs1, v_orbs2);

    const Ref<OrbitalSpace> occ_act = v_orbs1[0];
    const int nocc_act = occ_act->rank();
    if (debug_ >= DefaultPrintThresholds::N2) {
      ExEnv::out0() << endl << spinletters << " D^m_i " << endl
                    << "number of occupied orbital: " << nocc_act << endl;
    }

    // skip spincase if no electron of this kind
    if (nocc_act == 0)
      continue;

    // calculate AlphaAlpha/BetaBeta part for D^m_i:
    // R^ij_ab R^ab_mj & R^ji_ab R^ab_mj
    // which are contracted over a, b, and j
    const int nocc_occ = nocc_act * nocc_act;
    double* RRij_mj = new double[nocc_occ];
    double* RRji_mj = new double[nocc_occ];
    fill_n(RRij_mj, nocc_occ, 0.0);
    fill_n(RRji_mj, nocc_occ, 0.0);

    // R^ij_ab R^ab_mj = (f12f12)^ij_mj
    //                 - f12^ij_pq f12^pq_mj
    //                 - f12^ij_ma' f12^ma'_mj
    //                 - f12^ij_a'm f12^a'm_mj

    // activate_ints: (f12f12)^ij_mj & (f12)^ij_k1k2/(f12)^mj_k1k2 ints
    Ref<DistArray4> f12f12_ints;
    vector<Ref<DistArray4> > v_f12_ints;

    activate_ints(occ_act->id(), occ_act->id(), occ_act->id(), occ_act->id(),
                  descr_f12f12_key, moints4_rtime, f12f12_ints);
    activate_ints_X_f12(moints4_rtime, "ij",
                        v_orbs1, v_orbs2,
                        descr_f12_key, v_f12_ints);
    // test f12_ints
#if 0
    {
      const Ref<OrbitalSpace> vir = v_orbs1[1];
      const Ref<OrbitalSpace> orbs = v_orbs1[2];
      const Ref<OrbitalSpace> cabs = v_orbs1[3];
      const Ref<OrbitalSpace> occ = v_orbs1[4];

      // f12^b1b2_PQ
      Ref<DistArray4> p1p2_ints;
      activate_ints(occ_act->id(), occ_act->id(), orbs->id(), orbs->id(),
                    descr_f12_key, moints4_rtime, p1p2_ints);
      v_f12_ints.push_back(p1p2_ints);

      // f12^b1b2_MA'
      Ref<DistArray4> i1a2_ints;
      activate_ints(occ_act->id(), occ_act->id(), occ->id(), cabs->id(),
                    descr_f12_key, moints4_rtime, i1a2_ints);
      v_f12_ints.push_back(i1a2_ints);

      // f12^b1b2_A'M
      Ref<DistArray4> a1i2_ints;
      activate_ints(occ_act->id(), occ_act->id(), cabs->id(), occ->id(),
                    descr_f12_key, moints4_rtime, a1i2_ints);
      v_f12_ints.push_back(a1i2_ints);
    }
#endif

    // for RR^ij_mj:
    // 1st, 2nd, 3rd index: i, m, j
    compute_RR_sum_abj2(RR13_23, f12f12_idx, f12_idx,
                        f12f12_ints, v_f12_ints, v_f12_ints,
                        RRij_mj);

    // R^ji_ab R^ab_ij:
    compute_RR_sum_abj2(RR31_23, f12f12_idx, f12_idx,
                        f12f12_ints, v_f12_ints, v_f12_ints,
                        RRji_mj);

    f12f12_ints->deactivate();
    for (int i = 0; i != v_f12_ints.size(); ++i) {
      v_f12_ints[i]->deactivate();
    }

    if (debug_ >= DefaultPrintThresholds::N2) {
      print_intermediate(spinletters, "R^ij_ab R^ab_mj", RRij_mj, nocc_act, nocc_act);
      print_intermediate(spinletters, "R^ji_ab R^ab_mj", RRji_mj, nocc_act, nocc_act);
    }

    // calculate AlphaBeta part for D^i_i:
    double* RR1 = NULL;
    double* RR2 = NULL;
    double* RR3 = NULL;
    double* RR4 = NULL;
    if (nspincases2 == 3) {

      RR1 = new double[nocc_occ];
      RR2 = new double[nocc_occ];
      RR3 = new double[nocc_occ];
      RR4 = new double[nocc_occ];
      fill_n(RR1, nocc_occ, 0.0);
      fill_n(RR2, nocc_occ, 0.0);
      fill_n(RR3, nocc_occ, 0.0);
      fill_n(RR4, nocc_occ, 0.0);

      if (spin == Alpha) {
        // R^i1j2_a1b2 R^a1b2_m1j2
//        ExEnv::out0() << "R^i1j2_a1b2 R^a1b2_m1j2" <<endl;
        compute_RR_sum_abj2(RR13_23, f12f12_idx, f12_idx,
                            i1i2_f12f12_i1i2, f12_ints_i1i2, f12_ints_i1i2,
                            RR1);
        // R^j2i1_a1b2 R^a1b2_m1j2
//        ExEnv::out0() << "R^j2i1_a1b2 R^a1b2_m1j2" <<endl;
        compute_RR_sum_abj2(RR31_23, f12f12_idx, f12_idx,
                            i2i1_f12f12_i1i2, f12_ints_i2i1, f12_ints_i1i2,
                            RR2);
        // R^i1j2_a1b2 R^a1b2_j2m1
//        ExEnv::out0() << "R^i1j2_a1b2 R^a1b2_j2m1" <<endl;
        compute_RR_sum_abj2(RR13_32, f12f12_idx, f12_idx,
                            i1i2_f12f12_i2i1, f12_ints_i1i2, f12_ints_i2i1,
                            RR3);
        // R^j2i1_a1b2 R^a1b2_j2m1
//        ExEnv::out0() << "R^j2i1_a1b2 R^a1b2_j2m1" <<endl;
        compute_RR_sum_abj2(RR31_32, f12f12_idx, f12_idx,
                            i2i1_f12f12_i2i1, f12_ints_i2i1, f12_ints_i2i1,
                            RR4);
      } else {
          // R^j1i2_a1b2 R^a1b2_j1m2
          compute_RR_sum_abj2(RR31_32, f12f12_idx, f12_idx,
                              i1i2_f12f12_i1i2, f12_ints_i1i2, f12_ints_i1i2,
                              RR1);
          // R^i2j1_a1b2 R^a1b2_j1m2
          compute_RR_sum_abj2(RR13_32, f12f12_idx, f12_idx,
                              i2i1_f12f12_i1i2, f12_ints_i2i1, f12_ints_i1i2,
                              RR2);
          // R^j1i2_a1b2 R^a1b2_m2j1
          compute_RR_sum_abj2(RR31_23, f12f12_idx, f12_idx,
                              i1i2_f12f12_i2i1, f12_ints_i1i2, f12_ints_i2i1,
                              RR3);
          // R^i2j1_a1b2 R^a1b2_m2j1
          compute_RR_sum_abj2(RR13_23, f12f12_idx, f12_idx,
                              i2i1_f12f12_i2i1, f12_ints_i2i1, f12_ints_i2i1,
                              RR4);
      }

      if (debug_ >= DefaultPrintThresholds::N2) {
        const string spinlabel = (spin == Alpha? "alpha" : "beta");
        const string RR1_label = (spin == Alpha? "R^i1j2_a1b2 R^a1b2_m1j2"
                                               : "R^j1i2_a1b2 R^a1b2_j1m2");
        const string RR2_label = (spin == Alpha? "R^j2i1_a1b2 R^a1b2_m1j2"
                                               : "R^i2j1_a1b2 R^a1b2_j1m2");
        const string RR3_label = (spin == Alpha? "R^i1j2_a1b2 R^a1b2_j2m1"
                                               : "R^j1i2_a1b2 R^a1b2_m2j1");
        const string RR4_label = (spin == Alpha? "R^j2i1_a1b2 R^a1b2_j2m1"
                                               : "R^i2j1_a1b2 R^a1b2_m2j1");
        print_intermediate(spinlabel, RR1_label, RR1, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR2_label, RR2, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR3_label, RR3, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR4_label, RR4, nocc_act, nocc_act);
      }

    } else if (nspincases2 == 2) {
        RR1 = RRij_mj;
        RR2 = RRji_mj;
        RR3 = RRji_mj;
        RR4 = RRij_mj;
   }


    // calculate D^m_i
    double* iter_D = (spin == Alpha? Dm_i_alpha : Dm_i_beta);

    const double* iter_RRij_mj = RRij_mj;
    const double* iter_RRji_mj = RRji_mj;
    const double* iter_RR1 = RR1;
    const double* iter_RR2 = RR2;
    const double* iter_RR3 = RR3;
    const double* iter_RR4 = RR4;

    for (int m = 0;  m < nocc_act; ++m) {
      for (int i = 0;  i < nocc_act; ++i) {

        // AlphaAlpha/BetaBeta part: C1^2 * (R^IJ_AB R^AB_MJ - R^JI_AB R^AB_MJ)
        double dmi_12 = C_1 * C_1 * (*iter_RRij_mj - *iter_RRji_mj);
//        ExEnv::out0() << spinletters << " part d^" << m << "_" << i << " = "
//                    << scprintf("%12.10f", dmi_12) << endl;

        ++iter_RRij_mj;
        ++iter_RRji_mj;

        // AlphaBeta part:
        // Alpha D^m_i: [(C0+C1)/2]^2 * R^IJ_AB R^AB_MJ
        //               + (C0+C1)/2*(C0-C1)/2 * (R^JI_AB R^AB_MJ + R^IJ_AB R^AB_JM)
        //               + [(C0-C1)/2]^2 * R^JI_AB R^AB_JM
        // Beta D^m_i: [(C0+C1)/2]^2 * R^JI_AB R^AB_JM
        //             + (C0+C1)/2*(C0-C1)/2 * (R^IJ_AB R^AB_JM + R^IJ_AB R^AB_MJ)
        //             + [(C0-C1)/2]^2 * R^IJ_AB R^AB_MJ
        if (nocc_alpha != 0 && nocc_beta != 0) {
//        ExEnv::out0() << "RR1: "  << scprintf("%12.10f", RR1)
//                      << "  RR2: " << scprintf("%12.10f", RR2)
//                      << "  RR3: "  << scprintf("%12.10f", RR3)
//                      << "  RR4: " << scprintf("%12.10f", RR4)
//                      << endl;

           dmi_12 += pow(0.5 * (C_0 + C_1), 2) * (*iter_RR1)
                  + 0.25 * (C_0 * C_0 - C_1 * C_1) * (*iter_RR2 + *iter_RR3)
                  + pow(0.5 * (C_0 - C_1), 2) * (*iter_RR4);
//         ExEnv::out0() << "AlphaBeta part: d^"  << m << "_" << i << " = "
//                        << scprintf("%12.10f", dmi_12) << endl;

           ++iter_RR1;
           ++iter_RR2;
           ++iter_RR3;
           ++iter_RR4;
         } // end of AlphaBeta part
        *iter_D = - dmi_12;
        ++iter_D;
      }
    }

    const double* iter_Dmi = (spin == Alpha? Dm_i_alpha : Dm_i_beta);
    print_intermediate(spinletters, "D^m_i:", iter_Dmi, nocc_act, nocc_act);

     delete[] RRij_mj;
     delete[] RRji_mj;
     RRij_mj = NULL;
     RRji_mj = NULL;
     if (nspincases2 == 3) {
       delete[] RR1;
       delete[] RR2;
       delete[] RR3;
       delete[] RR4;
       RR1 = NULL;
       RR2 = NULL;
       RR3 = NULL;
       RR4 = NULL;
     }

  } // end of spincase1 loop

  if (nspincases2 == 3) {
    i1i2_f12f12_i1i2->deactivate();
    i2i1_f12f12_i1i2->deactivate();
    i1i2_f12f12_i2i1->deactivate();
    i2i1_f12f12_i2i1->deactivate();

    for (int i = 0; i != f12_ints_i1i2.size(); ++i) {
      f12_ints_i1i2[i]->deactivate();
      f12_ints_i2i1[i]->deactivate();
    }
  }
}
// end of computation of D^m_i

// compute D^m_i through R^ij_a'b' R_ij^a'b' + 2 R^ij_ab' R_ij^ab'
void MP2R12Energy_Diag::compute_Dmi_2(const int nspincases1, const int nspincases2,
                                    const double C_0, const double C_1,
                                    const vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                                    const vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                                    double* const Dm_i_alpha, double* const Dm_i_beta)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  Ref<TwoBodyIntDescr> descr_f12f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(), 0, 0);
  const string descr_f12f12_key = moints4_rtime->descr_key(descr_f12f12);

  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);
  const TwoBodyOper::type f12f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12f12();
  const unsigned int f12f12_idx = descr_f12f12->intset(f12f12_type);

  // preliminaries for AlphaBeta parts
  const Ref<OrbitalSpace> occ1_act = v_orbs1_ab[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2_ab[0];
  const Ref<OrbitalSpace> vir1 = v_orbs1_ab[1];
  const Ref<OrbitalSpace> vir2 = v_orbs2_ab[1];
  const Ref<OrbitalSpace> cabs1 = v_orbs1_ab[3];
  const Ref<OrbitalSpace> cabs2 = v_orbs2_ab[3];
  const int nocc_alpha = occ1_act->rank();
  const int nocc_beta = occ2_act->rank();
  // test
//  if (debug_ >= DefaultPrintThresholds::mostN2) {
//    ExEnv::out0() << endl << "Number of alpha occupied orbital: " << nocc_alpha
//                  << endl << "Number of beta occupied orbital: " << nocc_beta << endl;
//  }

  // activate integrals: R^i1i2_a'1a'2 or R^a'1a'2_i1i2, R^i2i1_a'1a'2 or R^a'1a'2_i1i2
  //                     R^i1i2_a1a'2 or R^a1a'2_i2i1, R^i2i1_a1a'2 or R^a1a'2_i2i1
  //                     R^i1i2_a'1a2 or R^a'1a'1_i2i1, R^i2i1_a'1a'1 R^a'1a2_i2i1
  // 1: alpha orbital, 2: beta orbital
  Ref<DistArray4> i1i2ap1ap2_ints;
  Ref<DistArray4> i2i1ap1ap2_ints;
  Ref<DistArray4> i1i2a1ap2_ints;
  Ref<DistArray4> i2i1a1ap2_ints;
  Ref<DistArray4> i1i2ap1a2_ints;
  Ref<DistArray4> i2i1ap1a2_ints;

  if (nspincases2 == 3) {
    // activate ints
     activate_ints(occ1_act->id(), occ2_act->id(), cabs1->id(), cabs2->id(),
                   descr_f12_key, moints4_rtime, i1i2ap1ap2_ints);

     activate_ints(occ2_act->id(), occ1_act->id(), cabs1->id(), cabs2->id(),
                   descr_f12_key, moints4_rtime, i2i1ap1ap2_ints);

     activate_ints(occ1_act->id(), occ2_act->id(), vir1->id(), cabs2->id(),
                   descr_f12_key, moints4_rtime, i1i2a1ap2_ints);

     activate_ints(occ2_act->id(), occ1_act->id(), vir1->id(), cabs2->id(),
                   descr_f12_key, moints4_rtime, i2i1a1ap2_ints);

     activate_ints(occ1_act->id(), occ2_act->id(), cabs1->id(), vir2->id(),
                   descr_f12_key, moints4_rtime, i1i2ap1a2_ints);

     activate_ints(occ2_act->id(), occ1_act->id(), cabs1->id(), vir2->id(),
                   descr_f12_key, moints4_rtime, i2i1ap1a2_ints);
  } // end of activating ints for alphabeta part

  for(int s = 0; s < nspincases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);

    // spin=0 (Alpha) => AlphaAlpha case (1)
    // spin=1 (Beta) => BetaBeta case (2)
    const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
    string spinletters = to_string(spincase);

    // obtain occ, vir, orbs, and cabs orbitals in that order
    vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
    vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
    obtain_orbitals(spincase, v_orbs1, v_orbs2);

    const Ref<OrbitalSpace> occ_act = v_orbs1[0];
    const Ref<OrbitalSpace> vir = v_orbs1[1];
    const Ref<OrbitalSpace> cabs = v_orbs1[3];
    const int nocc_act = occ_act->rank();
    if (debug_ >= DefaultPrintThresholds::mostN2) {
      ExEnv::out0() << endl << spinletters << " D^m_i " << endl
                    << "number of occupied orbital: " << nocc_act << endl;
    }

    // skip spincase if no electron of this kind
    if (nocc_act == 0)
      continue;

    // calculate AlphaAlpha/BetaBeta part for D^m_i:
    const int nocc_occ = nocc_act * nocc_act;

    // R^ij_a'b' R^a'b'_mj & R^ji_a'b' R^a'b'_mj;
    double* RRij_mj_ap = new double[nocc_occ];
    double* RRji_mj_ap = new double[nocc_occ];
    fill_n(RRij_mj_ap, nocc_occ, 0.0);
    fill_n(RRji_mj_ap, nocc_occ, 0.0);

    // activate_ints
    Ref<DistArray4> iiapap_ints;
    activate_ints(occ_act->id(), occ_act->id(), cabs->id(), cabs->id(),
                  descr_f12_key, moints4_rtime, iiapap_ints);

    // R^ij_a'b' R^a'b'_mj
    compute_RR2_sum_3idx(RR13_23, f12_idx,
                         iiapap_ints, iiapap_ints,
                         RRij_mj_ap);
    // R^ji_a'b' R^a'b'_mj
    compute_RR2_sum_3idx(RR31_23, f12_idx,
                         iiapap_ints, iiapap_ints,
                         RRji_mj_ap);

    iiapap_ints->deactivate();

    if (debug_ >= DefaultPrintThresholds::mostN2) {
      print_intermediate(spinletters, "R^ij_a'b' R^a'b'_mj", RRij_mj_ap, nocc_act, nocc_act);
      print_intermediate(spinletters, "R^ji_a'b' R^a'b'_mj", RRji_mj_ap, nocc_act, nocc_act);
    }

#if 0
    // test Tr D^m_i ?= Tr D^b_c + Tr D^b'_c'
    double trace = 0;
    for (int i = 0; i != nocc_act; ++i) {
      const int idx = i * nocc_act + i;
      trace += RRij_mj_ap[idx] - RRji_mj_ap[idx];
    }
    ExEnv::out0() << endl << spinletters << " trace of R^ij_a'b' R^a'b'_mj: " << scprintf("%12.10f", trace) << endl;
#endif

    // R^ij_ab' R^ab'_mj & R^ji_ab' R^ab'_mj
    double* RRij_mj_a = new double[nocc_occ];
    double* RRji_mj_a = new double[nocc_occ];
    fill_n(RRij_mj_a, nocc_occ, 0.0);
    fill_n(RRji_mj_a, nocc_occ, 0.0);

    // activate_ints
    Ref<DistArray4> iiaap_ints;
    activate_ints(occ_act->id(), occ_act->id(), vir->id(), cabs->id(),
                  descr_f12_key, moints4_rtime, iiaap_ints);

    // R^ij_ab' R^ab'_mj
    compute_RR2_sum_3idx(RR13_23, f12_idx,
                         iiaap_ints, iiaap_ints,
                         RRij_mj_a);
    // R^ji_ab' R^ab'_mj
    compute_RR2_sum_3idx(RR31_23, f12_idx,
                         iiaap_ints, iiaap_ints,
                         RRji_mj_a);

    // R^ij_ab' R^ab'_jm & R^ji_ab' R^ab'_jm
    double* RRij_jm_a = new double[nocc_occ];
    double* RRji_jm_a = new double[nocc_occ];
    fill_n(RRij_jm_a, nocc_occ, 0.0);
    fill_n(RRji_jm_a, nocc_occ, 0.0);

    // R^ij_ab' R^ab'_jm
    compute_RR2_sum_3idx(RR13_32, f12_idx,
                         iiaap_ints, iiaap_ints,
                         RRij_jm_a);
    // R^ji_ab' R^ab'_jm
    compute_RR2_sum_3idx(RR31_32, f12_idx,
                         iiaap_ints, iiaap_ints,
                         RRji_jm_a);

    iiaap_ints->deactivate();

    if (debug_ >= DefaultPrintThresholds::mostN2) {
      print_intermediate(spinletters, "R^ij_ab' R^ab'_mj", RRij_mj_a, nocc_act, nocc_act);
      print_intermediate(spinletters, "R^ji_ab' R^ab'_mj", RRji_mj_a, nocc_act, nocc_act);
      print_intermediate(spinletters, "R^ij_ab' R^ab'_jm", RRij_jm_a, nocc_act, nocc_act);
      print_intermediate(spinletters, "R^ji_ab' R^ab'_jm", RRji_jm_a, nocc_act, nocc_act);
    }

#if 0
    // test Tr D^m_i ?= Tr D^b_c + Tr D^b'_c'
    double trace = 0;
    for (int i = 0; i != nocc_act; ++i) {
      const int idx = i * nocc_act + i;
      trace += RRij_mj_a[idx] - RRji_mj_a[idx];
    }
    ExEnv::out0() << endl << spinletters << " trace of R^ij_ab' R^ab'_mj: " << scprintf("%12.10f", trace) << endl;
#endif

    // calculate AlphaBeta part for D^i_i:
    double* RR1_ap = NULL;
    double* RR2_ap = NULL;
    double* RR3_ap = NULL;
    double* RR4_ap = NULL;

    // sum over a in alpha
    double* RR1_a_alpha = NULL;
    double* RR2_a_alpha = NULL;
    double* RR3_a_alpha = NULL;
    double* RR4_a_alpha = NULL;

    // sum over a in beta
    double* RR1_a_beta = NULL;
    double* RR2_a_beta = NULL;
    double* RR3_a_beta = NULL;
    double* RR4_a_beta = NULL;
    if (nspincases2 == 3) {

      RR1_ap = new double[nocc_occ];
      RR2_ap = new double[nocc_occ];
      RR3_ap = new double[nocc_occ];
      RR4_ap = new double[nocc_occ];
      fill_n(RR1_ap, nocc_occ, 0.0);
      fill_n(RR2_ap, nocc_occ, 0.0);
      fill_n(RR3_ap, nocc_occ, 0.0);
      fill_n(RR4_ap, nocc_occ, 0.0);

      RR1_a_alpha = new double[nocc_occ];
      RR2_a_alpha = new double[nocc_occ];
      RR3_a_alpha = new double[nocc_occ];
      RR4_a_alpha = new double[nocc_occ];
      fill_n(RR1_a_alpha, nocc_occ, 0.0);
      fill_n(RR2_a_alpha, nocc_occ, 0.0);
      fill_n(RR3_a_alpha, nocc_occ, 0.0);
      fill_n(RR4_a_alpha, nocc_occ, 0.0);

      RR1_a_beta = new double[nocc_occ];
      RR2_a_beta = new double[nocc_occ];
      RR3_a_beta = new double[nocc_occ];
      RR4_a_beta = new double[nocc_occ];
      fill_n(RR1_a_beta, nocc_occ, 0.0);
      fill_n(RR2_a_beta, nocc_occ, 0.0);
      fill_n(RR3_a_beta, nocc_occ, 0.0);
      fill_n(RR4_a_beta, nocc_occ, 0.0);

      if (spin == Alpha) {
        // R^i1j2_a'1b'2 R^a'1b'2_m1j2
        compute_RR2_sum_3idx(RR13_23, f12_idx,
                             i1i2ap1ap2_ints, i1i2ap1ap2_ints,
                             RR1_ap);
        // R^j2i1_a'1b'2 R^a'1b'2_m1j2
        compute_RR2_sum_3idx(RR31_23, f12_idx,
                             i2i1ap1ap2_ints, i1i2ap1ap2_ints,
                             RR2_ap);
        // R^i1j2_a'1b'2 R^a'1b'2_j2m1
        compute_RR2_sum_3idx(RR13_32, f12_idx,
                             i1i2ap1ap2_ints, i2i1ap1ap2_ints,
                             RR3_ap);
        // R^j2i1_a'1b'2 R^a'1b'2_j2m1
        compute_RR2_sum_3idx(RR31_32, f12_idx,
                             i2i1ap1ap2_ints, i2i1ap1ap2_ints,
                             RR4_ap);

        // R^i1j2_a1b'2 R^a1b'2_m1j2
        compute_RR2_sum_3idx(RR13_23, f12_idx,
                             i1i2a1ap2_ints, i1i2a1ap2_ints,
                             RR1_a_alpha);
        // R^j2i1_a1b'2 R^a1b'2_m1j2
        compute_RR2_sum_3idx(RR31_23, f12_idx,
                             i2i1a1ap2_ints, i1i2a1ap2_ints,
                             RR2_a_alpha);
        // R^i1j2_a1b'2 R^a1b'2_j2m1
        compute_RR2_sum_3idx(RR13_32, f12_idx,
                             i1i2a1ap2_ints, i2i1a1ap2_ints,
                             RR3_a_alpha);
        // R^j2i1_a1b'2 R^a1b'2_j2m1
        compute_RR2_sum_3idx(RR31_32, f12_idx,
                             i2i1a1ap2_ints, i2i1a1ap2_ints,
                             RR4_a_alpha);

        // R^i1j2_b'1a2 R^b'1a2_m1j2
        compute_RR2_sum_3idx(RR13_23, f12_idx,
                             i1i2ap1a2_ints, i1i2ap1a2_ints,
                             RR1_a_beta);
        // R^j2i1_b'1a2 R^b'1a2_m1j2
        compute_RR2_sum_3idx(RR31_23, f12_idx,
                             i2i1ap1a2_ints, i1i2ap1a2_ints,
                             RR2_a_beta);
        // R^i1j2_b'1a2 R^b'1a2_j2m1
        compute_RR2_sum_3idx(RR13_32, f12_idx,
                             i1i2ap1a2_ints, i2i1ap1a2_ints,
                             RR3_a_beta);
        // R^j2i1_b'1a2 R^b'1a2_j2m1
        compute_RR2_sum_3idx(RR31_32, f12_idx,
                             i2i1ap1a2_ints, i2i1ap1a2_ints,
                             RR4_a_beta);

      } else {
          // R^j1i2_a'1b'2 R^a'1b'2_j1m2
          compute_RR2_sum_3idx(RR31_32, f12_idx,
                               i1i2ap1ap2_ints, i1i2ap1ap2_ints,
                               RR1_ap);
          // R^i2j1_a'1b'2 R^a'1b'2_j1m2
          compute_RR2_sum_3idx(RR13_32, f12_idx,
                               i2i1ap1ap2_ints, i1i2ap1ap2_ints,
                               RR2_ap);
          // R^j1i2_a'1b'2 R^a'1b'2_m2j1
          compute_RR2_sum_3idx(RR31_23, f12_idx,
                               i1i2ap1ap2_ints, i2i1ap1ap2_ints,
                               RR3_ap);
          // R^i2j1_a'1b'2 R^a'1b'2_m2j1
          compute_RR2_sum_3idx(RR13_23, f12_idx,
                               i2i1ap1ap2_ints, i2i1ap1ap2_ints,
                               RR4_ap);

          // R^j1i2_a1b'2 R^a1b'2_j1m2
          compute_RR2_sum_3idx(RR31_32, f12_idx,
                               i1i2a1ap2_ints, i1i2a1ap2_ints,
                               RR1_a_alpha);
          // R^i2j1_a1b'2 R^a1b'2_j1m2
          compute_RR2_sum_3idx(RR13_32, f12_idx,
                               i2i1a1ap2_ints, i1i2a1ap2_ints,
                               RR2_a_alpha);
          // R^j1i2_a1b'2 R^a1b'2_m2j1
          compute_RR2_sum_3idx(RR31_23, f12_idx,
                               i1i2a1ap2_ints, i2i1a1ap2_ints,
                               RR3_a_alpha);
          // R^i2j1_a1b'2 R^a1b'2_m2j1
          compute_RR2_sum_3idx(RR13_23, f12_idx,
                               i2i1a1ap2_ints, i2i1a1ap2_ints,
                               RR4_a_alpha);

          // R^j1i2_b'1a2 R^b'1a2_j1m2
          compute_RR2_sum_3idx(RR31_32, f12_idx,
                               i1i2ap1a2_ints, i1i2ap1a2_ints,
                               RR1_a_beta);
          // R^i2j1_b'1a2 R^b'1a2_j1m2
          compute_RR2_sum_3idx(RR13_32, f12_idx,
                               i2i1ap1a2_ints, i1i2ap1a2_ints,
                               RR2_a_beta);
          // R^j1i2_b'1a2 R^b'1a2_m2j1
          compute_RR2_sum_3idx(RR31_23, f12_idx,
                               i1i2ap1a2_ints, i2i1ap1a2_ints,
                               RR3_a_beta);
          // R^i2j1_b'1a2 R^b'1a2_m2j1
          compute_RR2_sum_3idx(RR13_23, f12_idx,
                               i2i1ap1a2_ints, i2i1ap1a2_ints,
                               RR4_a_beta);
      }

      if (debug_ >= DefaultPrintThresholds::mostN2) {
        const string spinlabel = (spin == Alpha? "alpha" : "beta");
        const string RR1_ap_label = (spin == Alpha? "R^i1j2_a'1b'2 R^a'1b'2_m1j2"
                                               : "R^j1i2_a'1b'2 R^a'1b'2_j1m2");
        const string RR2_ap_label = (spin == Alpha? "R^j2i1_a'1b'2 R^a'1b'2_m1j2"
                                               : "R^i2j1_a'1b'2 R^a'1b'2_j1m2");
        const string RR3_ap_label = (spin == Alpha? "R^i1j2_a'1b'2 R^a'1b'2_j2m1"
                                               : "R^j1i2_a'1b'2 R^a'1b'2_m2j1");
        const string RR4_ap_label = (spin == Alpha? "R^j2i1_a'1b'2 R^a'1b'2_j2m1"
                                               : "R^i2j1_a'1b'2 R^a'1b'2_m2j1");
        print_intermediate(spinlabel, RR1_ap_label, RR1_ap, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR2_ap_label, RR2_ap, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR3_ap_label, RR3_ap, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR4_ap_label, RR4_ap, nocc_act, nocc_act);

        const string RR1_a_alpha_label = (spin == Alpha? "R^i1j2_a1b'2 R^a1b'2_m1j2"
                                               : "R^j1i2_a1b'2 R^a1b'2_j1m2");
        const string RR2_a_alpha_label = (spin == Alpha? "R^j2i1_a1b'2 R^a1b'2_m1j2"
                                               : "R^i2j1_a1b'2 R^a1b'2_j1m2");
        const string RR3_a_alpha_label = (spin == Alpha? "R^i1j2_a1b'2 R^a1b'2_j2m1"
                                               : "R^j1i2_a1b'2 R^a1b'2_m2j1");
        const string RR4_a_alpha_label = (spin == Alpha? "R^j2i1_a1b'2 R^a1b'2_j2m1"
                                               : "R^i2j1_a1b'2 R^a1b'2_m2j1");
        print_intermediate(spinlabel, RR1_a_alpha_label, RR1_a_alpha, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR2_a_alpha_label, RR2_a_alpha, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR3_a_alpha_label, RR3_a_alpha, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR4_a_alpha_label, RR4_a_alpha, nocc_act, nocc_act);

        const string RR1_a_beta_label = (spin == Alpha? "R^i1j2_b'1a2 R^b'1a2_m1j2"
                                               : "R^j1i2_b'1a2 R^b'1a2_j1m2");
        const string RR2_a_beta_label = (spin == Alpha? "R^j2i1_b'1a2 R^b'1a2_m1j2"
                                               : "R^i2j1_b'1a2 R^b'1a2_j1m2");
        const string RR3_a_beta_label = (spin == Alpha? "R^i1j2_b'1a2 R^b'1a2_j2m1"
                                               : "R^j1i2_b'1a2 R^b'1a2_m2j1");
        const string RR4_a_beta_label = (spin == Alpha? "R^j2i1_b'1a2 R^b'1a2_j2m1"
                                               : "R^i2j1_b'1a2 R^b'1a2_m2j1");
        print_intermediate(spinlabel, RR1_a_beta_label, RR1_a_beta, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR2_a_beta_label, RR2_a_beta, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR3_a_beta_label, RR3_a_beta, nocc_act, nocc_act);
        print_intermediate(spinlabel, RR4_a_beta_label, RR4_a_beta, nocc_act, nocc_act);
      }

      // test Tr D^m_i ?= Tr D^b_c + Tr D^b'_c'
#if 0
    double trace_RR1_alpha = 0;
    double trace_RR2_alpha = 0;
    double trace_RR3_alpha = 0;
    double trace_RR4_alpha = 0;
    double trace_RR1_beta = 0;
    double trace_RR2_beta = 0;
    double trace_RR3_beta = 0;
    double trace_RR4_beta = 0;
    for (int i = 0; i != nocc_act; ++i) {
      const int idx = i * nocc_act + i;
      trace_RR1_alpha += RR1_a_alpha[idx];
      trace_RR2_alpha += RR2_a_alpha[idx];
      trace_RR3_alpha += RR3_a_alpha[idx];
      trace_RR4_alpha += RR4_a_alpha[idx];

      trace_RR1_beta += RR1_a_beta[idx];
      trace_RR2_beta += RR2_a_beta[idx];
      trace_RR3_beta += RR3_a_beta[idx];
      trace_RR4_beta += RR4_a_beta[idx];
    }
    ExEnv::out0() << endl << spinletters << " RR1_a_alpha: " << scprintf("%12.10f", trace_RR1_alpha) << endl;
    ExEnv::out0() << endl << spinletters << " RR2_a_alpha: " << scprintf("%12.10f", trace_RR2_alpha) << endl;
    ExEnv::out0() << endl << spinletters << " RR3_a_alpha: " << scprintf("%12.10f", trace_RR3_alpha) << endl;
    ExEnv::out0() << endl << spinletters << " RR4_a_alpha: " << scprintf("%12.10f", trace_RR4_alpha) << endl;

    ExEnv::out0() << endl << spinletters << " RR1_a_beta: " << scprintf("%12.10f", trace_RR1_beta) << endl;
    ExEnv::out0() << endl << spinletters << " RR2_a_beta: " << scprintf("%12.10f", trace_RR2_beta) << endl;
    ExEnv::out0() << endl << spinletters << " RR3_a_beta: " << scprintf("%12.10f", trace_RR3_beta) << endl;
    ExEnv::out0() << endl << spinletters << " RR4_a_beta: " << scprintf("%12.10f", trace_RR4_beta) << endl;
#endif

    } else if (nspincases2 == 2) {
        RR1_ap = RRij_mj_ap;
        RR2_ap = RRji_mj_ap;
        RR3_ap = RRji_mj_ap;
        RR4_ap = RRij_mj_ap;

        RR1_a_alpha = RRij_mj_a;
        RR2_a_alpha = RRji_mj_a;
        RR3_a_alpha = RRij_jm_a;
        RR4_a_alpha = RRji_jm_a;

        RR1_a_beta = RRji_jm_a;
        RR2_a_beta = RRij_jm_a;
        RR3_a_beta = RRji_mj_a;
        RR4_a_beta = RRij_mj_a;
   }


    // calculate D^m_i
    double* iter_D = (spin == Alpha? Dm_i_alpha : Dm_i_beta);

    const double* iter_RRij_mj_ap = RRij_mj_ap;
    const double* iter_RRji_mj_ap = RRji_mj_ap;
    const double* iter_RR1_ap = RR1_ap;
    const double* iter_RR2_ap = RR2_ap;
    const double* iter_RR3_ap = RR3_ap;
    const double* iter_RR4_ap = RR4_ap;

    const double* iter_RRij_mj_a = RRij_mj_a;
    const double* iter_RRji_mj_a = RRji_mj_a;
    const double* iter_RRij_jm_a = RRij_jm_a;
    const double* iter_RRji_jm_a = RRji_jm_a;
    const double* iter_RR1_a_alpha = RR1_a_alpha;
    const double* iter_RR2_a_alpha = RR2_a_alpha;
    const double* iter_RR3_a_alpha = RR3_a_alpha;
    const double* iter_RR4_a_alpha = RR4_a_alpha;
    const double* iter_RR1_a_beta = RR1_a_beta;
    const double* iter_RR2_a_beta = RR2_a_beta;
    const double* iter_RR3_a_beta = RR3_a_beta;
    const double* iter_RR4_a_beta = RR4_a_beta;

    for (int m = 0;  m < nocc_act; ++m) {
      for (int i = 0;  i < nocc_act; ++i) {

        // AlphaAlpha/BetaBeta part: C1^2 * (R^IJ_AB R^AB_MJ - R^JI_AB R^AB_MJ)
        double dmi_12 = C_1 * C_1 * ( (*iter_RRij_mj_ap - *iter_RRji_mj_ap)
                                     + (*iter_RRij_mj_a - *iter_RRji_mj_a
                                      - *iter_RRij_jm_a + *iter_RRji_jm_a)
                                    );

//        double dmi_12 = C_1 * C_1 * (*iter_RRij_mj_ap - *iter_RRji_mj_ap);

//          double dmi_12 = 2* C_1 * C_1 * (*iter_RRij_mj_a - *iter_RRji_mj_a);

        ++iter_RRij_mj_ap;
        ++iter_RRji_mj_ap;
        ++iter_RRij_mj_a;
        ++iter_RRji_mj_a;
        ++iter_RRij_jm_a;
        ++iter_RRji_jm_a;

        // AlphaBeta part:
        // Alpha D^m_i: [(C0+C1)/2]^2 * R^IJ_AB R^AB_MJ
        //               + (C0+C1)/2*(C0-C1)/2 * (R^JI_AB R^AB_MJ + R^IJ_AB R^AB_JM)
        //               + [(C0-C1)/2]^2 * R^JI_AB R^AB_JM
        // Beta D^m_i: [(C0+C1)/2]^2 * R^JI_AB R^AB_JM
        //             + (C0+C1)/2*(C0-C1)/2 * (R^IJ_AB R^AB_JM + R^IJ_AB R^AB_MJ)
        //             + [(C0-C1)/2]^2 * R^IJ_AB R^AB_MJ
        if (nocc_alpha != 0 && nocc_beta != 0) {

           dmi_12 += pow(0.5 * (C_0 + C_1), 2) * (*iter_RR1_ap
                                                  + (*iter_RR1_a_alpha + *iter_RR1_a_beta)
                                                  )
                  + 0.25 * (C_0 * C_0 - C_1 * C_1) * ((*iter_RR2_ap + *iter_RR3_ap)
                                                     + (*iter_RR2_a_alpha + *iter_RR3_a_alpha
                                                      + *iter_RR2_a_beta + *iter_RR3_a_beta)
                                                     )
                  + pow(0.5 * (C_0 - C_1), 2) * (*iter_RR4_ap
                                                 + (*iter_RR4_a_alpha + *iter_RR4_a_beta)
                                                 );
//             dmi_12 += pow(0.5 * (C_0 + C_1), 2) * (*iter_RR1_ap)
//                    + 0.25 * (C_0 * C_0 - C_1 * C_1) * (*iter_RR2_ap + *iter_RR3_ap)
//                    + pow(0.5 * (C_0 - C_1), 2) * (*iter_RR4_ap);

//             dmi_12 += pow(0.5 * (C_0 + C_1), 2) *  (*iter_RR1_a_alpha + *iter_RR1_a_beta)
//                    + 0.25 * (C_0 * C_0 - C_1 * C_1) *  (*iter_RR2_a_alpha + *iter_RR3_a_alpha
//                                                       + *iter_RR2_a_beta + *iter_RR3_a_beta)
//                    + pow(0.5 * (C_0 - C_1), 2) *  (*iter_RR4_a_alpha + *iter_RR4_a_beta);

           ++iter_RR1_ap;
           ++iter_RR2_ap;
           ++iter_RR3_ap;
           ++iter_RR4_ap;
           ++iter_RR1_a_alpha;
           ++iter_RR2_a_alpha;
           ++iter_RR3_a_alpha;
           ++iter_RR4_a_alpha;
           ++iter_RR1_a_beta;
           ++iter_RR2_a_beta;
           ++iter_RR3_a_beta;
           ++iter_RR4_a_beta;
         } // end of AlphaBeta part
        *iter_D = - dmi_12;
        ++iter_D;
      }
    }

     delete[] RRij_mj_ap;
     delete[] RRji_mj_ap;
     RRij_mj_ap = NULL;
     RRji_mj_ap = NULL;
     delete[] RRij_mj_a;
     delete[] RRji_mj_a;
     RRij_mj_a = NULL;
     RRji_mj_a = NULL;
     if (nspincases2 == 3) {
       delete[] RR1_ap;
       delete[] RR2_ap;
       delete[] RR3_ap;
       delete[] RR4_ap;
       RR1_ap = NULL;
       RR2_ap = NULL;
       RR3_ap = NULL;
       RR4_ap = NULL;
       delete[] RR1_a_alpha;
       delete[] RR2_a_alpha;
       delete[] RR3_a_alpha;
       delete[] RR4_a_alpha;
       RR1_a_alpha = NULL;
       RR2_a_alpha = NULL;
       RR3_a_alpha = NULL;
       RR4_a_alpha = NULL;
       delete[] RR1_a_beta;
       delete[] RR2_a_beta;
       delete[] RR3_a_beta;
       delete[] RR4_a_beta;
       RR1_a_beta = NULL;
       RR2_a_beta = NULL;
       RR3_a_beta = NULL;
       RR4_a_beta = NULL;
     }

  } // end of spincase1 loop

  if (nspincases2 == 3) {
    i1i2ap1ap2_ints->deactivate();
    i2i1ap1ap2_ints->deactivate();
    i1i2a1ap2_ints->deactivate();
    i2i1a1ap2_ints->deactivate();
  }
}
// end of computation of D^m_i_2


// compute \bar{\tilde{R}}^31_ij \bar{\tilde{R}}^ij_32 with spin orbitals:
// alpha:  D^1_2 =
// index 1,2,3,i,j: alpha orbitals
//              + C1^2 * (R^31_ij R^32_ij - R^13_ij R^32_ij)
// index 1,2,i: alpha orbitals; index 3,j: beta orbitals
//              + [(C0+C1)/2]^2 * R^13_ij R^23_ij
//              + (C0+C1)/2*(C0-C1)/2 * (R^31_ij R^23_ij + R^13_ij R^32_ij)
//              + [(C0-C1)/2]^2 * R^31_ij R^32_ij
//
// beta:  D^1_2 =
// index 1,2,3,i,j: beta orbitals
//              + C1^2 * (R^31_ij R^32_ij - R^13_ij R^32_ij)
// index 1,2,j: beta orbitals; index 3,i: alpha orbitals
//              + [(C0+C1)/2]^2 * R^31_ij R^32_ij
//              + (C0+C1)/2*(C0-C1)/2 * (R^13_ij R^23_ij + R^31_ij R^23_ij)
//              + [(C0-C1)/2]^2 * R^13_ij R^23_ij
void MP2R12Energy_Diag::compute_RR31_32_spin(const int orbitals_label,
                                             const int nspincases1, const int nspincases2,
                                             const double C_0, const double C_1,
                                             const vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                                             const vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                                             double* const D_alpha, double* const D_beta)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);

  const Ref<OrbitalSpace> occ1_act = v_orbs1_ab[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2_ab[0];
  const Ref<OrbitalSpace> vir1 = v_orbs1_ab[1];
  const Ref<OrbitalSpace> vir2 = v_orbs2_ab[1];
  const Ref<OrbitalSpace> cabs1 = v_orbs1_ab[3];
  const Ref<OrbitalSpace> cabs2 = v_orbs2_ab[3];

  const int nocc_alpha = occ1_act->rank();
  const int nocc_beta = occ2_act->rank();

  // test propose
//  if (debug_ >= DefaultPrintThresholds::mostN2) {
//    ExEnv::out0() << endl << "number of alpha active occupied orbital: " << occ1_act->rank() << endl
//                  <<"number of beta active occupied orbital: " << occ2_act->rank() << endl
//                  << "number of alpha virtula orbital: " << vir1->rank() << endl
//                  <<"number of beta virtual orbital: " << vir2->rank() << endl
//                  << "number of alpha cabs: " << cabs1->rank() << endl
//                  <<"number of beta cabs: " << cabs2->rank() << endl;
//  }

  // can not use Ref<OrbitalSpace>&
  Ref<OrbitalSpace> idx1_orbs1 = vir1;
  Ref<OrbitalSpace> idx1_orbs2 = vir2;
  Ref<OrbitalSpace> idx2_orbs1 = vir1;
  Ref<OrbitalSpace> idx2_orbs2 = vir2;
  Ref<OrbitalSpace> idx3_orbs1 = cabs1;
  Ref<OrbitalSpace> idx3_orbs2 = cabs2;

  switch (orbitals_label) {
  case vir_vir_cabs:
    // default value
  break;

  case cabs_cabs_vir:
    idx1_orbs1 = cabs1;
    idx1_orbs2 = cabs2;
    idx2_orbs1 = cabs1;
    idx2_orbs2 = cabs2;
    idx3_orbs1 = vir1;
    idx3_orbs2 = vir2;
  break;

  case cabs_cabs_cabs:
    idx1_orbs1 = cabs1;
    idx1_orbs2 = cabs2;
    idx2_orbs1 = cabs1;
    idx2_orbs2 = cabs2;
  break;

  case cabs_vir_cabs:
    idx1_orbs1 = cabs1;
    idx1_orbs2 = cabs2;
  break;

  default:
    ExEnv::out0() << "There is no such case";
    break;
  }

  // test
  if (debug_ >= DefaultPrintThresholds::mostN2) {
    const int nidx1_orbs1 = idx1_orbs1->rank();
    const int nidx1_orbs2 = idx1_orbs2->rank();
    const int nidx2_orbs1 = idx2_orbs1->rank();
    const int nidx2_orbs2 = idx2_orbs2->rank();
    const int nidx3_orbs1 = idx3_orbs1->rank();
    const int nidx3_orbs2 = idx3_orbs2->rank();
    ExEnv::out0() << endl << "AlphaBeta test: "
                  << endl << "Number of index1 alpha orbital: " << nidx1_orbs1
                  << endl << "Number of index1 beta orbital: " << nidx1_orbs2
                  << endl << "Number of index2 alpha orbital: " << nidx2_orbs1
                  << endl << "Number of index2 beta orbital: " << nidx2_orbs2
                  << endl << "number of index3 alpha orbital: " << nidx3_orbs1
                  << endl << "number of index3 beta orbital: " << nidx3_orbs2 << endl;
  }

  for(int s = 0; s < nspincases1; ++s) {
     const SpinCase1 spin = static_cast<SpinCase1>(s);

     // spin=0 (Alpha) => AlphaAlpha case (1)
     // spin=1 (Beta) => BetaBeta case (2)
     const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
     string spinletters = to_string(spincase);

     // obtain occ_act, vir, orbs, cabs, and occ orbitals in that order
     vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
     vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
     obtain_orbitals(spincase, v_orbs1, v_orbs2);

     const Ref<OrbitalSpace> occ_act = v_orbs1[0];
     const int nocc_act = occ_act->rank();

     // skip spincase if no electron of this kind
     if (nocc_act == 0)
       continue;

     const Ref<OrbitalSpace> vir = v_orbs1[1];
     const Ref<OrbitalSpace> cabs = v_orbs1[3];
     const int nvir = vir->rank();
     const int ncabs = cabs->rank();

     int size_idx1 = nvir;
     int size_idx2 = nvir;
     // can not use Ref<OrbitalSpace>&, change the values in 'switch statement'
     Ref<OrbitalSpace> idx1_orbs = vir;
     Ref<OrbitalSpace> idx2_orbs = vir;
     Ref<OrbitalSpace> idx3_orbs = cabs;

     switch (orbitals_label) {
     case vir_vir_cabs:
       // default value
     break;

     case cabs_cabs_vir:
       size_idx1 = ncabs;
       size_idx2 = ncabs;
       idx1_orbs = cabs;
       idx2_orbs = cabs;
       idx3_orbs = vir;
     break;

     case cabs_cabs_cabs:
       size_idx1 = ncabs;
       size_idx2 = ncabs;
       idx1_orbs = cabs;
       idx2_orbs = cabs;
     break;

     case cabs_vir_cabs:
       size_idx1 = ncabs;
       idx1_orbs = cabs;
     break;
     }

     if (debug_ >= DefaultPrintThresholds::mostN2) {
       int size_idx3 = ncabs;
       if (orbitals_label == cabs_cabs_vir) size_idx3 = nvir;
       ExEnv::out0() << endl << spinletters << endl
                     << "number of index1 orbital: " << idx1_orbs->rank() << endl
                     << "number of index2 orbital: " << idx2_orbs->rank() << endl
                     << "number of index3 orbital: " << idx3_orbs->rank() << endl;
     }

     // test AlphaAlpha/BetaBeta part
#if 0
    {
       Ref<DistArray4> f12_ints1;
       activate_ints(cabs->id(), cabs->id(), occ_act->id(), occ_act->id(),
                     descr_f12_key, moints4_rtime, f12_ints1);
//       print_f12_ints(spinletters, "R^A'B_IJ ints", f12_idx, f12_ints1);

       Ref<DistArray4> f12_ints2;
       activate_ints(cabs->id(), vir->id(), occ_act->id(), occ_act->id(),
                     descr_f12_key, moints4_rtime, f12_ints2);


       const string RR_label1 = "testing result for R^B'A'_IJ R^IJ_B'A";
       const string RR_label2 = "testing result for R^A'B'_IJ R^IJ_B'A";

       const blasint one = 1;
       const blasint nocc12 = nocc_act * nocc_act;
       const int size_idx3 = ncabs;

       double* RR = new double[size_idx1 * size_idx2];
       double* iter_RR = RR;
       for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
         for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

           double RR_sum_idx3 = 0;
           for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
             const double* blk1 = f12_ints1->retrieve_pair_block(idx3, idx1, f12_idx);
             const double* blk2 = f12_ints2->retrieve_pair_block(idx3, idx2, f12_idx);
             const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

             RR_sum_idx3 += blk12;
             f12_ints1->release_pair_block(idx3, idx1, f12_idx);
             f12_ints2->release_pair_block(idx3, idx2, f12_idx);

           }
           *iter_RR = RR_sum_idx3;
           ++iter_RR;
         }
       }
       print_intermediate(spinletters, RR_label1, RR, size_idx1, size_idx2);

       iter_RR = RR;
       for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
         for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

           double RR_sum_idx3 = 0;
           for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
             const double* blk1 = f12_ints1->retrieve_pair_block(idx1, idx3, f12_idx);
             const double* blk2 = f12_ints2->retrieve_pair_block(idx3, idx2, f12_idx);
             const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

             RR_sum_idx3 += blk12;
             f12_ints1->release_pair_block(idx1, idx3, f12_idx);
             f12_ints2->release_pair_block(idx3, idx2, f12_idx);

           }
           *iter_RR = RR_sum_idx3;
           ++iter_RR;
         }
       }
       print_intermediate(spinletters, RR_label2, RR, size_idx1, size_idx2);
       delete[] RR;

       f12_ints1->deactivate();
       f12_ints2->deactivate();
     } // end of testing ints
#endif

     // calculate AlphaAlpha/BetaBeta part for D^idx1_idx2:
     // R^31_ij R^ij_32 & R^13_ij R^ij_32 contracted over i, j, and index 3
     // 1 and 2 are external indices: in same orbitals
     const int size_idx12 = size_idx1 * size_idx2;
     double* RR31_32_array = new double[size_idx12];
     double* RR13_32_array = new double[size_idx12];
     fill_n(RR31_32_array, size_idx12, 0.0);
     fill_n(RR13_32_array, size_idx12, 0.0);

     // activate integrals
     Ref<DistArray4> R31_ii_ints;
     Ref<DistArray4> R13_ii_ints;
     Ref<DistArray4> R32_ii_ints;

     activate_ints(idx3_orbs->id(), idx1_orbs->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, R31_ii_ints);

     switch (orbitals_label) {
     case vir_vir_cabs:
       activate_ints(idx1_orbs->id(), idx3_orbs->id(), occ_act->id(), occ_act->id(),
                     descr_f12_key, moints4_rtime, R13_ii_ints);
       R32_ii_ints = R31_ii_ints;
     break;

     case cabs_cabs_vir:
       activate_ints(idx1_orbs->id(), idx3_orbs->id(), occ_act->id(), occ_act->id(),
                     descr_f12_key, moints4_rtime, R13_ii_ints);
       R32_ii_ints = R31_ii_ints;
     break;

     case cabs_cabs_cabs:
       R13_ii_ints = R31_ii_ints;
       R32_ii_ints = R31_ii_ints;
     break;

     case cabs_vir_cabs:
       R13_ii_ints = R31_ii_ints;
       activate_ints(idx3_orbs->id(), idx2_orbs->id(), occ_act->id(), occ_act->id(),
                     descr_f12_key, moints4_rtime, R32_ii_ints);
     break;
     }

     // R^31_ij R^ij_32:
     compute_RR2_sum_3idx(RR31_32, f12_idx,
                          R31_ii_ints, R32_ii_ints,
                          RR31_32_array);
     // R^13_ij R^ij_32:
     compute_RR2_sum_3idx(RR13_32, f12_idx,
                          R13_ii_ints, R32_ii_ints,
                          RR13_32_array);

     R31_ii_ints->deactivate();
     switch (orbitals_label) {
     case vir_vir_cabs:
       R13_ii_ints->deactivate();
     break;

     case cabs_cabs_vir:
       R13_ii_ints->deactivate();
     break;

     case cabs_cabs_cabs:
     break;

     case cabs_vir_cabs:
       R32_ii_ints->deactivate();
     break;
     }

     if (debug_ >= DefaultPrintThresholds::mostN2) {
       string RR31_32_label;
       string RR13_32_label;
       switch (orbitals_label) {
       case vir_vir_cabs:
         RR31_32_label = "R^a'b_ij R^ij_a'c";
         RR13_32_label = "R^ba'_ij R^ij_a'c";
       break;

       case cabs_cabs_vir:
         RR31_32_label = "R^ab'_ij R^ij_ac'";
         RR13_32_label = "R^b'a_ij R^ij_ac'";
       break;

       case cabs_cabs_cabs:
         RR31_32_label = "R^a'b'_ij R^ij_a'c'";
         RR13_32_label = "R^b'a'_ij R^ij_a'c'";
       break;

       case cabs_vir_cabs:
         RR31_32_label = "R^b'a'_ij R^ij_b'a";
         RR13_32_label = "R^a'b'_ij R^ij_b'a";
       break;
       }

       print_intermediate(spinletters, RR31_32_label, RR31_32_array, size_idx1, size_idx2);
       print_intermediate(spinletters, RR13_32_label, RR13_32_array, size_idx1, size_idx2);
     }

     // test Tr D^m_i ?= Tr D^b_c + Tr D^b'_c'
#if 0
     double Tr_RR = 0;
     for (int i = 0; i != size_idx1; ++i) {
       const int idx = i * size_idx2 + i;
       Tr_RR += RR31_32_array[idx] - RR13_32_array[idx];
     }
     ExEnv::out0() << endl << spinletters << " trace of RR31_32: " << scprintf("%12.10f", Tr_RR) << endl;
#endif

     // calculate AlphaBeta part for D^idx1_idx2:
     //            alpha                      beta
     // RR1:   R^13_ij R^ij_23   or   R^31_ij R^ij_32
     // RR2:   R^31_ij R^ij_23   or   R^13_ij R^ij_32
     // RR3:   R^13_ij R^ij_32   or   R^31_ij R^ij_23
     // RR4:   R^31_ij R^ij_32   or   R^13_ij R^ij_23
     double* RR1 = NULL;
     double* RR2 = NULL;
     double* RR3 = NULL;
     double* RR4 = NULL;

     if (nspincases2 == 3) {

       RR1 = new double[size_idx12];
       RR2 = new double[size_idx12];
       RR3 = new double[size_idx12];
       RR4 = new double[size_idx12];
       fill_n(RR1, size_idx12, 0.0);
       fill_n(RR2, size_idx12, 0.0);
       fill_n(RR3, size_idx12, 0.0);
       fill_n(RR4, size_idx12, 0.0);

       if (spin == Alpha) {
         Ref<DistArray4> R13ii_1212_ints;   // R^1(e1)3(e2)_i(e1)j(e2)
         Ref<DistArray4> R31ii_2112_ints;   // R^3(e2)1(e1)_i(e1)j(e2)
         Ref<DistArray4> R23ii_1212_ints;   // R^2(e1)3(e2)_i(e1)j(e2)
         Ref<DistArray4> R32ii_2112_ints;   // R^3(e2)2(e1)_i(e1)j(e2)

         activate_ints(idx1_orbs1->id(), idx3_orbs2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, R13ii_1212_ints);
         activate_ints(idx3_orbs2->id(), idx1_orbs1->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, R31ii_2112_ints);

         if (orbitals_label == cabs_vir_cabs) {
           activate_ints(idx2_orbs1->id(), idx3_orbs2->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, R23ii_1212_ints);
           activate_ints(idx3_orbs2->id(), idx2_orbs1->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, R32ii_2112_ints);
          } else {
              R23ii_1212_ints = R13ii_1212_ints;
              R32ii_2112_ints = R31ii_2112_ints;
          }

         // R^1(e1)3(e2)_i(e1)j(e2) * R^i(e1)j(e2)_2(e1)3(e2)
         compute_RR2_sum_3idx(RR13_23, f12_idx,
                              R13ii_1212_ints, R23ii_1212_ints,
                              RR1);
         // R^3(e2)1(e1)_i(e1)j(e2) * R^i(e1)j(e2)_2(e1)3(e2)
         compute_RR2_sum_3idx(RR31_23, f12_idx,
                              R31ii_2112_ints, R23ii_1212_ints,
                              RR2);
         // R^1(e1)3(e2)_i(e1)j(e2) * R^i(e1)j(e2)_3(e2)2(e1)
         compute_RR2_sum_3idx(RR13_32, f12_idx,
                              R13ii_1212_ints, R32ii_2112_ints,
                              RR3);
         // R^3(e2)1(e1)_i(e1)j(e2) * R^i(e1)j(e2)_3(e2)2(e1)
         compute_RR2_sum_3idx(RR31_32, f12_idx,
                              R31ii_2112_ints, R32ii_2112_ints,
                              RR4);

         R13ii_1212_ints->deactivate();
         R31ii_2112_ints->deactivate();
         if (orbitals_label == cabs_vir_cabs) {
           R23ii_1212_ints->deactivate();
           R32ii_2112_ints->deactivate();
         }

         // test for AlphaBeta of alpha spin
#if 0
        {
           Ref<DistArray4> f12_ints1;
           activate_ints(cabs1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, f12_ints1);

           Ref<DistArray4> f12_ints2;
           activate_ints(cabs2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, f12_ints2);

            Ref<DistArray4> f12_ints21;
            activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                          descr_f12_key, moints4_rtime, f12_ints21);

            Ref<DistArray4> f12_ints22;
            activate_ints(cabs2->id(), vir1->id(), occ1_act->id(), occ2_act->id(),
                          descr_f12_key, moints4_rtime, f12_ints22);

           const string RR_label1 = "testing result for R^A'B'_IJ R^IJ_AB'";
           const string RR_label2 = "testing result for R^B'A'_IJ R^IJ_AB'";
           const string RR_label3 = "testing result for R^A'B'_IJ R^IJ_B'A";
           const string RR_label4 = "testing result for R^B'A'_IJ R^IJ_B'A";

           const int size_idx1 = cabs1->rank();
           const int size_idx2 = vir1->rank();
           const int size_idx3 = cabs2->rank();


           const blasint one = 1;
           const blasint nocc12 = nocc_alpha * nocc_beta;
           double* RR = new double[size_idx1 * size_idx2];

           double* iter_RR = RR;
           for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
             for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

               double RR_sum_idx3 = 0;
               for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                 const double* blk1 = f12_ints1->retrieve_pair_block(idx1, idx3, f12_idx);
                 const double* blk2 = f12_ints21->retrieve_pair_block(idx2, idx3, f12_idx);
                 const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                 RR_sum_idx3 += blk12;
                 f12_ints1->release_pair_block(idx1, idx3, f12_idx);
                 f12_ints21->release_pair_block(idx2, idx3, f12_idx);

               }
               *iter_RR = RR_sum_idx3;
               ++iter_RR;
             }
           }
           print_intermediate(spinletters, RR_label1, RR, size_idx1, size_idx2);

           iter_RR = RR;
           for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
             for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

               double RR_sum_idx3 = 0;
               for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                 const double* blk1 = f12_ints2->retrieve_pair_block(idx3, idx1, f12_idx);
                 const double* blk2 = f12_ints21->retrieve_pair_block(idx2, idx3, f12_idx);
                 const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                 RR_sum_idx3 += blk12;
                 f12_ints2->release_pair_block(idx3, idx1, f12_idx);
                 f12_ints21->release_pair_block(idx2, idx3, f12_idx);

               }
               *iter_RR = RR_sum_idx3;
               ++iter_RR;
             }
           }
           print_intermediate(spinletters, RR_label2, RR, size_idx1, size_idx2);

           iter_RR = RR;
           for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
             for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

               double RR_sum_idx3 = 0;
               for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                 const double* blk1 = f12_ints1->retrieve_pair_block(idx1, idx3, f12_idx);
                 const double* blk2 = f12_ints22->retrieve_pair_block(idx3, idx2, f12_idx);
                 const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                 RR_sum_idx3 += blk12;
                 f12_ints1->release_pair_block(idx1, idx3, f12_idx);
                 f12_ints22->release_pair_block(idx3, idx2, f12_idx);
               }
               *iter_RR = RR_sum_idx3;
               ++iter_RR;
             }
           }
           print_intermediate(spinletters, RR_label3, RR, size_idx1, size_idx2);

           iter_RR = RR;
            for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
              for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

                double RR_sum_idx3 = 0;
                for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                  const double* blk1 = f12_ints2->retrieve_pair_block(idx3, idx1, f12_idx);
                  const double* blk2 = f12_ints22->retrieve_pair_block(idx3, idx2, f12_idx);
                  const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                  RR_sum_idx3 += blk12;
                  f12_ints2->release_pair_block(idx3, idx1, f12_idx);
                  f12_ints22->release_pair_block(idx3, idx2, f12_idx);

                }
                *iter_RR = RR_sum_idx3;
                ++iter_RR;
              }
            }
            print_intermediate(spinletters, RR_label4, RR, size_idx1, size_idx2);

           delete[] RR;
           f12_ints1->deactivate();
           f12_ints2->deactivate();
           f12_ints21->deactivate();
           f12_ints22->deactivate();
         } // end of testing ints
#endif

       } else {
           // test for AlphaBeta of beta spin
#if 0
          {
             Ref<DistArray4> f12_ints1;
             activate_ints(cabs1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                           descr_f12_key, moints4_rtime, f12_ints1);

             Ref<DistArray4> f12_ints2;
             activate_ints(cabs2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                           descr_f12_key, moints4_rtime, f12_ints2);

             Ref<DistArray4> f12_ints21;
             activate_ints(cabs1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                           descr_f12_key, moints4_rtime, f12_ints21);

             Ref<DistArray4> f12_ints22;
             activate_ints(vir2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                           descr_f12_key, moints4_rtime, f12_ints22);

             const string RR_label1 = "testing result for R^B'A'_IJ R^IJ_B'A";
             const string RR_label2 = "testing result for R^A'B'_IJ R^IJ_B'A";
             const string RR_label3 = "testing result for R^B'A'_IJ R^IJ_AB'";
             const string RR_label4 = "testing result for R^A'B'_IJ R^IJ_AB'";

             const int size_idx1 = cabs2->rank();
             const int size_idx2 = vir2->rank();
             const int size_idx3 = cabs1->rank();

             const blasint one = 1;
             const blasint nocc12 = nocc_alpha * nocc_beta;
             double* RR = new double[size_idx1 * size_idx2];

             double* iter_RR = RR;
             for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
               for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

                 double RR_sum_idx3 = 0;
                 for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                   const double* blk1 = f12_ints1->retrieve_pair_block(idx3, idx1, f12_idx);
                   const double* blk2 = f12_ints21->retrieve_pair_block(idx3, idx2, f12_idx);
                   const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                   RR_sum_idx3 += blk12;
                   f12_ints1->release_pair_block(idx3, idx1, f12_idx);
                   f12_ints21->release_pair_block(idx3, idx2, f12_idx);

                 }
                 *iter_RR = RR_sum_idx3;
                 ++iter_RR;
               }
             }
             print_intermediate(spinletters, RR_label1, RR, size_idx1, size_idx2);

             iter_RR = RR;
             for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
               for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

                 double RR_sum_idx3 = 0;
                 for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                   const double* blk1 = f12_ints2->retrieve_pair_block(idx1, idx3, f12_idx);
                   const double* blk2 = f12_ints21->retrieve_pair_block(idx3, idx2, f12_idx);
                   const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                   RR_sum_idx3 += blk12;
                   f12_ints2->release_pair_block(idx1, idx3, f12_idx);
                   f12_ints21->release_pair_block(idx3, idx2, f12_idx);

                 }
                 *iter_RR = RR_sum_idx3;
                 ++iter_RR;
               }
             }
             print_intermediate(spinletters, RR_label2, RR, size_idx1, size_idx2);

             iter_RR = RR;
             for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
               for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

                 double RR_sum_idx3 = 0;
                 for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                   const double* blk1 = f12_ints1->retrieve_pair_block(idx3, idx1, f12_idx);
                   const double* blk2 = f12_ints22->retrieve_pair_block(idx2, idx3, f12_idx);
                   const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                   RR_sum_idx3 += blk12;
                   f12_ints1->release_pair_block(idx3, idx1, f12_idx);
                   f12_ints22->release_pair_block(idx2, idx3, f12_idx);
                 }
                 *iter_RR = RR_sum_idx3;
                 ++iter_RR;
               }
             }
             print_intermediate(spinletters, RR_label3, RR, size_idx1, size_idx2);

             iter_RR = RR;
              for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
                for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

                  double RR_sum_idx3 = 0;
                  for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                    const double* blk1 = f12_ints2->retrieve_pair_block(idx1, idx3, f12_idx);
                    const double* blk2 = f12_ints22->retrieve_pair_block(idx2, idx3, f12_idx);
                    const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                    RR_sum_idx3 += blk12;
                    f12_ints2->release_pair_block(idx1, idx3, f12_idx);
                    f12_ints22->release_pair_block(idx2, idx3, f12_idx);

                  }
                  *iter_RR = RR_sum_idx3;
                  ++iter_RR;
                }
              }
              print_intermediate(spinletters, RR_label4, RR, size_idx1, size_idx2);

             delete[] RR;
             f12_ints1->deactivate();
             f12_ints2->deactivate();
             f12_ints21->deactivate();
             f12_ints22->deactivate();
           } // end of testing ints
#endif
           Ref<DistArray4> R31ii_1212_ints;   // R^3(e1)1(e2)_i(e1)j(e2)
           Ref<DistArray4> R13ii_2112_ints;   // R^1(e2)3(e1)_i(e1)j(e2)
           Ref<DistArray4> R32ii_1212_ints;   // R^3(e1)2(e2)_i(e1)j(e2)
           Ref<DistArray4> R23ii_2112_ints;   // R^2(e2)3(e1)_i(e1)j(e2)

           activate_ints(idx3_orbs1->id(), idx1_orbs2->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, R31ii_1212_ints);
           activate_ints(idx1_orbs2->id(), idx3_orbs1->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, R13ii_2112_ints);

           if (orbitals_label == cabs_vir_cabs) {
             activate_ints(idx3_orbs1->id(), idx2_orbs2->id(), occ1_act->id(), occ2_act->id(),
                           descr_f12_key, moints4_rtime, R32ii_1212_ints);
             activate_ints(idx2_orbs2->id(), idx3_orbs1->id(), occ1_act->id(), occ2_act->id(),
                           descr_f12_key, moints4_rtime, R23ii_2112_ints);
           } else {
               R32ii_1212_ints = R31ii_1212_ints;
               R23ii_2112_ints = R13ii_2112_ints;
           }

           // R^3(e1)1(e2)_i(e1)j(e2) * R^i(e1)j(e2)_3(e1)2(e2)
           compute_RR2_sum_3idx(RR31_32, f12_idx,
                                R31ii_1212_ints, R32ii_1212_ints,
                                RR1);
           // R^1(e2)3(e1)_i(e1)j(e2) * R^i(e1)j(e2)_3(e1)2(e2)
           compute_RR2_sum_3idx(RR13_32, f12_idx,
                                R13ii_2112_ints, R32ii_1212_ints,
                                RR2);
           // R^3(e1)1(e2)_i(e1)j(e2) * R^i(e1)j(e2)_2(e2)3(e1)
           compute_RR2_sum_3idx(RR31_23, f12_idx,
                                R31ii_1212_ints, R23ii_2112_ints,
                                RR3);
           // R^1(e2)3(e1)_i(e1)j(e2) * R^i(e1)j(e2)_3(e1)2(e2)
           compute_RR2_sum_3idx(RR13_23, f12_idx,
                                R13ii_2112_ints, R23ii_2112_ints,
                                RR4);

           R31ii_1212_ints->deactivate();
           R13ii_2112_ints->deactivate();
           if (orbitals_label == cabs_vir_cabs) {
             R32ii_1212_ints->deactivate();
             R23ii_2112_ints->deactivate();
           }
       }

       if (debug_ >= DefaultPrintThresholds::mostN2) {
         const string spin_label = (spin == Alpha? "alpha" : "beta");
         string RR1_label;
         string RR2_label;
         string RR3_label;
         string RR4_label;

         switch (orbitals_label) {
         case vir_vir_cabs:
           RR1_label = (spin == Alpha? "R^ba'_ij R^ij_ca'" : "R^a'b_ij R^ij_a'c");
           RR2_label = (spin == Alpha? "R^a'b_ij R^ij_ca'" : "R^ba'_ij R^ij_a'c");
           RR3_label = (spin == Alpha? "R^ba'_ij R^ij_a'c" : "R^a'b_ij R^ij_ca'");
           RR4_label = (spin == Alpha? "R^a'b_ij R^ij_a'c" : "R^ba'_ij R^ij_ca'");
         break;

         case cabs_cabs_vir:
           RR1_label = (spin == Alpha? "R^b'a_ij R^ij_c'a" : "R^ab'_ij R^ij_ac'");
           RR2_label = (spin == Alpha? "R^ab'_ij R^ij_c'a" : "R^b'a_ij R^ij_ac'");
           RR3_label = (spin == Alpha? "R^b'a_ij R^ij_ac'" : "R^ab'_ij R^ij_c'a");
           RR4_label = (spin == Alpha? "R^ab'_ij R^ij_ac'" : "R^b'a_ij R^ij_c'a");
         break;

         case cabs_cabs_cabs:
           RR1_label = (spin == Alpha? "R^b'a'_ij R^ij_c'a'" : "R^a'b'_ij R^ij_a'c'");
           RR2_label = (spin == Alpha? "R^a'b'_ij R^ij_c'a'" : "R^b'a'_ij R^ij_a'c'");
           RR3_label = (spin == Alpha? "R^b'a'_ij R^ij_a'c'" : "R^a'b'_ij R^ij_c'a'");
           RR4_label = (spin == Alpha? "R^a'b'_ij R^ij_a'c'" : "R^b'a'_ij R^ij_c'a'");
         break;

         case cabs_vir_cabs:
           RR1_label = (spin == Alpha? "R^a'b'_ij R^ij_ab'" : "R^b'a'_ij R^ij_b'a");
           RR2_label = (spin == Alpha? "R^b'a'_ij R^ij_ab'" : "R^a'b'_ij R^ij_b'a");
           RR3_label = (spin == Alpha? "R^a'b'_ij R^ij_b'a" : "R^b'a'_ij R^ij_ab'");
           RR4_label = (spin == Alpha? "R^b'a'_ij R^ij_b'a" : "R^a'b'_ij R^ij_ab'");
         break;
         }

         print_intermediate(spin_label, RR1_label, RR1, size_idx1, size_idx2);
         print_intermediate(spin_label, RR2_label, RR2, size_idx1, size_idx2);
         print_intermediate(spin_label, RR3_label, RR3, size_idx1, size_idx2);
         print_intermediate(spin_label, RR4_label, RR4, size_idx1, size_idx2);
       }

       // test Tr D^m_i ?= Tr D^b_c + Tr D^b'_c'
#if 0
    double trace_RR1 = 0;
    double trace_RR2 = 0;
    double trace_RR3 = 0;
    double trace_RR4 = 0;
    for (int i = 0; i != size_idx1; ++i) {
      const int idx = i * size_idx2 + i;
      trace_RR1 += RR1[idx];
      trace_RR2 += RR2[idx];
      trace_RR3 += RR3[idx];
      trace_RR4 += RR4[idx];
    }
    ExEnv::out0() << endl << spinletters << " RR1: " << scprintf("%12.10f", trace_RR1) << endl;
    ExEnv::out0() << endl << spinletters << " RR2: " << scprintf("%12.10f", trace_RR2) << endl;
    ExEnv::out0() << endl << spinletters << " RR3: " << scprintf("%12.10f", trace_RR3) << endl;
    ExEnv::out0() << endl << spinletters << " RR4: " << scprintf("%12.10f", trace_RR4) << endl;
#endif

     } else if (nspincases2 == 2) {
         RR1 = RR31_32_array;
         RR2 = RR13_32_array;
         RR3 = RR13_32_array;
         RR4 = RR31_32_array;
    }// end of AlphaBeta part for D

     // calculate D[s]
     const double* iter_RR31_32_array = RR31_32_array;
     const double* iter_RR13_32_array = RR13_32_array;
     const double* iter_RR1 = RR1;
     const double* iter_RR2 = RR2;
     const double* iter_RR3 = RR3;
     const double* iter_RR4 = RR4;

     double* iter_D = (spin == Alpha? D_alpha : D_beta);
     for (int idx1 = 0;  idx1 < size_idx1; ++idx1) {
       for (int idx2 = 0; idx2 < size_idx2; ++idx2) {

         // AlphaAlpha/BetaBeta part
         double d_12 = C_1 * C_1 * (*iter_RR31_32_array - *iter_RR13_32_array);
//         ExEnv::out0() << spinletters << " part d^" << idx1 << "_" << idx2 << " = "
//                       << scprintf("%12.10f", d_12) << endl;

         ++iter_RR31_32_array;
         ++iter_RR13_32_array;

         // AlphaBeta part
         if (nocc_alpha != 0 && nocc_beta != 0) {
   //        ExEnv::out0() << "RR1: "  << scprintf("%12.10f", RR1)
   //                      << "  RR2: " << scprintf("%12.10f", RR2)
   //                      << "  RR3: "  << scprintf("%12.10f", RR3)
   //                      << "  RR4: " << scprintf("%12.10f", RR4)
   //                      << endl;

           d_12 += pow(0.5 * (C_0 + C_1), 2) * (*iter_RR1)
                  + 0.25 * (C_0 * C_0 - C_1 * C_1) * (*iter_RR2 + *iter_RR3)
                  + pow(0.5 * (C_0 - C_1), 2) * (*iter_RR4);
//           ExEnv::out0() << "AlphaBeta part: d^"  << idx1 << "_" << idx2 << " = "
//                         << scprintf("%12.10f", d_12) << endl;

           ++iter_RR1;
           ++iter_RR2;
           ++iter_RR3;
           ++iter_RR4;
         } // end of AlphaBeta part

         *iter_D = d_12;
         ++iter_D;
       } // end of looping over c
     } // end of calculating D^b_c[s]

     delete[] RR31_32_array;
     delete[] RR13_32_array;
     if (nspincases2 == 3) {
       delete[] RR1;
       delete[] RR2;
       delete[] RR3;
       delete[] RR4;
     }
  } // end of spincase1 loop

}
// end of computation of compute_RR31_32_spin

// test function for D^c_b
void MP2R12Energy_Diag::compute_Dcb(const int nspincases1, const int nspincases2,
                                    const double C_0, const double C_1)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);

  // obtain occ, vir, orbs, and cabs orbitals for AlphaBeta case
//  vector< Ref<OrbitalSpace> > v_orbs1_ab;
//  vector< Ref<OrbitalSpace> > v_orbs2_ab;
//  obtain_orbitals(AlphaBeta, v_orbs1_ab, v_orbs2_ab);

  const SpinCase1 spin1 = case1(AlphaBeta);
  const SpinCase1 spin2 = case2(AlphaBeta);

  const Ref<OrbitalSpace> occ1_act = r12eval()->occ_act(spin1);
  const Ref<OrbitalSpace> vir1 = r12eval()->vir(spin1);
  const Ref<OrbitalSpace> orbs1 = r12eval()->orbs(spin1);
  const Ref<OrbitalSpace> cabs1 = r12eval()->cabs_space_hcanonical(spin1);

  const Ref<OrbitalSpace> occ2_act = r12eval()->occ_act(spin2);
  const Ref<OrbitalSpace> vir2 = r12eval()->vir(spin2);
  const Ref<OrbitalSpace> orbs2 = r12eval()->orbs(spin2);
  const Ref<OrbitalSpace> cabs2 = r12eval()->cabs_space_hcanonical(spin2);

  const int nvir1 = vir1->rank();
  const int nvir2 = vir2->rank();
  const int ncabs1 = cabs1->rank();
  const int ncabs2 = cabs2->rank();

  const int nocc_alpha = occ1_act->rank();
  const int nocc_beta = occ2_act->rank();

  if (debug_ >= DefaultPrintThresholds::N2)
    ExEnv::out0() << endl << "AlphaBeta D^b_c"
                  << endl << "Number of alpha virtual orbital: " << nvir1
                  << endl << "Number of beta virtual orbital: " << nvir2
                  << endl << "number of alpha cabs orbital: " << ncabs1
                  << endl << "number of beta cabs orbital: " << ncabs2 << endl;

  for(int s = 0; s < nspincases1; ++s) {
     const SpinCase1 spin = static_cast<SpinCase1>(s);

     // spin=0 (Alpha) => AlphaAlpha case (1)
     // spin=1 (Beta) => BetaBeta case (2)
     const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
     string spinletters = to_string(spincase);

     // obtain occ, vir, orbs, and cabs orbitals in that order
//     vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
//     vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
//     obtain_orbitals(spincase, v_orbs1, v_orbs2);

     const Ref<OrbitalSpace> occ_act = r12eval()->occ_act(spin);
     const Ref<OrbitalSpace> vir = r12eval()->vir(spin);
     const Ref<OrbitalSpace> cabs = r12eval()->cabs_space_hcanonical(spin);

     const int nocc_act = occ_act->rank();

     // skip spincase if no electron of this kind
     if (nocc_act == 0)
       continue;

     const int nvir = vir->rank();
     const int ncabs = cabs->rank();
     if (debug_ >= DefaultPrintThresholds::N2)
       ExEnv::out0() << endl << spinletters << " D^b_c" << endl
                     << "number of occupied orbital: " << nocc_act << endl
                     << "number of virtual orbital: " << nvir << endl
                     << "number of cabs orbital: " << ncabs << endl;

     // calculate AlphaAlpha/BetaBeta part for D^c_b:
     // R^a'b_ij R^ij_a'c & R^ba'_ij R^ij_a'c contracted over i,j,a'
     const int nvir12 = nvir * nvir;
     double* RRapb_apc = new double[nvir12];
     double* RRbap_apc = new double[nvir12];
     fill_n(RRapb_apc, nvir12, 0.0);
     fill_n(RRbap_apc, nvir12, 0.0);

     // activate integrals
     Ref<DistArray4> apaii_f12_ints;
     activate_ints(cabs->id(), vir->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, apaii_f12_ints);

//     // test ints
//     // f12^b1b2_A'M or f12^A'M_k1k2
//     {
//       Ref<DistArray4> f12_ints;
//       activate_ints(occ_act->id(), occ_act->id(), cabs->id(), cabs->id(),
//                     descr_f12_key, moints4_rtime, f12_ints);
//       const int size_idx1 = f12_ints->ni();
//       const int size_idx2 = f12_ints->nj();
//       const int size_idx3 = f12_ints->nx();
//       const int size_idx4 = f12_ints->ny();
//       for (int idx1 = 0; idx1 < size_idx1; ++idx1) {
//         for(int idx2 = 0; idx2 < size_idx2; ++idx2) {
//           ExEnv::out0() << indent << "index1: " << idx1 <<"  index2:"<< idx2 << endl;
//
//           const double* f12_ij_blk = f12_ints->retrieve_pair_block(idx1, idx2, f12_idx);
//           for (int idx3 = 0; idx3 < size_idx3; ++idx3) {
//             for(int idx4 = 0; idx4 < size_idx4; ++idx4) {
//
//               ExEnv::out0() << indent << scprintf("%12.10f", *f12_ij_blk) << " ";
//               ++f12_ij_blk;
//             }
//             ExEnv::out0() << endl;
//           }
//           f12_ints->release_pair_block(idx1, idx2, f12_idx);
//           ExEnv::out0() << endl;
//         }
//         ExEnv::out0() << endl;
//       }
//     } // end of test ints


     Ref<DistArray4> aapii_f12_ints;
     activate_ints(vir->id(), cabs->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, aapii_f12_ints);

     // b is the 1st index, c the 2nd, a' the 3rd
     // R^a'b_ij R^ij_a'c
//     compute_RR2_sum_3idx(RR31_32, f12_idx,
//                          apaii_f12_ints, apaii_f12_ints,
//                          RRapb_apc);

     double* iter_array = RRapb_apc;
     const blasint one = 1;
     const blasint nocc12 = nocc_act * nocc_act;
     for (int c = 0; c < nvir; ++c) {
       for(int b = 0; b < nvir; ++b) {

         double f12f12_sum_ap = 0;
         for(int ap = 0; ap < ncabs; ++ap) {
           const double* f12_ij_blk1 = apaii_f12_ints->retrieve_pair_block(ap, b, f12_idx);
           const double* f12_ij_blk2 = apaii_f12_ints->retrieve_pair_block(ap, c, f12_idx);

           const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

           f12f12_sum_ap += f12f12_sum_ij;

           apaii_f12_ints->release_pair_block(ap, b, f12_idx);
           apaii_f12_ints->release_pair_block(ap, c, f12_idx);
         }
         *iter_array = f12f12_sum_ap;
         ++iter_array;
       }
     }

     // R^ba'_ij R^ij_a'c
//     compute_RR2_sum_3idx(RR13_32, f12_idx,
//                          aapii_f12_ints, apaii_f12_ints,
//                          RRbap_apc);

     iter_array = RRbap_apc;
     for (int c = 0; c < nvir; ++c) {
       for(int b = 0; b < nvir; ++b) {

         double f12f12_sum_ap = 0;
         for(int ap = 0; ap < ncabs; ++ap) {
           const double* f12_ij_blk1 = aapii_f12_ints->retrieve_pair_block(b, ap, f12_idx);
           const double* f12_ij_blk2 = apaii_f12_ints->retrieve_pair_block(ap, c, f12_idx);

           const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

           f12f12_sum_ap += f12f12_sum_ij;

           aapii_f12_ints->release_pair_block(b, ap, f12_idx);
           apaii_f12_ints->release_pair_block(ap, c, f12_idx);
         }
         *iter_array = f12f12_sum_ap;
         ++iter_array;
       }
     }

     apaii_f12_ints->deactivate();
     aapii_f12_ints->deactivate();

     if (debug_ >= DefaultPrintThresholds::mostN2) {
       print_intermediate(spinletters, "R^a'b_ij R^ij_a'c", RRapb_apc, nvir, nvir);
       print_intermediate(spinletters, "R^ba'_ij R^ij_a'c", RRbap_apc, nvir, nvir);
     }

     // calculate AlphaBeta part for D^c_b:
     //            alpha                      beta
     // RR1:   R^ba'_ij R^ij_ca'   or   R^a'b_ij R^ij_a'c
     // RR2:   R^a'b_ij R^ij_ca'   or   R^ba'_ij R^ij_a'c
     // RR3:   R^ba'_ij R^ij_a'c   or   R^a'b_ij R^ij_ca'
     // RR4:   R^a'b_ij R^ij_a'c   or   R^ba'_ij R^ij_ca'
     double* RR1 = NULL;
     double* RR2 = NULL;
     double* RR3 = NULL;
     double* RR4 = NULL;

     if (nspincases2 == 3) {

       RR1 = new double[nvir12];
       RR2 = new double[nvir12];
       RR3 = new double[nvir12];
       RR4 = new double[nvir12];
       fill_n(RR1, nvir12, 0.0);
       fill_n(RR2, nvir12, 0.0);
       fill_n(RR3, nvir12, 0.0);
       fill_n(RR4, nvir12, 0.0);

       if (spin == Alpha) {
         Ref<DistArray4> a1ap2i1i2_f12_ints;
         Ref<DistArray4> ap2a1i1i2_f12_ints;
         activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, a1ap2i1i2_f12_ints);
         activate_ints(cabs2->id(), vir1->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, ap2a1i1i2_f12_ints);

//         compute_RR2_sum_3idx(RR13_23, f12_idx,
//                              a1ap2i1i2_f12_ints, a1ap2i1i2_f12_ints,
//                              RR1); // R^ba'_ij R^ij_ca'
//         compute_RR2_sum_3idx(RR31_23, f12_idx,
//                              ap2a1i1i2_f12_ints, a1ap2i1i2_f12_ints,
//                              RR2); // R^a'b_ij R^ij_ca'
//         compute_RR2_sum_3idx(RR13_32, f12_idx,
//                              a1ap2i1i2_f12_ints, ap2a1i1i2_f12_ints,
//                              RR3); // R^ba'_ij R^ij_a'c
//         compute_RR2_sum_3idx(RR31_32, f12_idx,
//                              ap2a1i1i2_f12_ints, ap2a1i1i2_f12_ints,
//                              RR4); // R^a'b_ij R^ij_a'c

         const blasint nocc12 = nocc_alpha * nocc_beta;

         iter_array = RR1;
         for (int c = 0; c < nvir; ++c) {
           for(int b = 0; b < nvir; ++b) {

             double f12f12_sum_ap = 0;
             for(int ap = 0; ap < ncabs2; ++ap) {
               const double* f12_ij_blk1 = a1ap2i1i2_f12_ints->retrieve_pair_block(b, ap, f12_idx);
               const double* f12_ij_blk2 = a1ap2i1i2_f12_ints->retrieve_pair_block(c, ap, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_ap += f12f12_sum_ij;

               a1ap2i1i2_f12_ints->release_pair_block(b, ap, f12_idx);
               a1ap2i1i2_f12_ints->release_pair_block(c, ap, f12_idx);
             }
             *iter_array = f12f12_sum_ap;
             ++iter_array;
           }
         }

         iter_array = RR2;
         for (int c = 0; c < nvir; ++c) {
           for(int b = 0; b < nvir; ++b) {

             double f12f12_sum_ap = 0;
             for(int ap = 0; ap < ncabs2; ++ap) {
               const double* f12_ij_blk1 = ap2a1i1i2_f12_ints->retrieve_pair_block(ap, b, f12_idx);
               const double* f12_ij_blk2 = a1ap2i1i2_f12_ints->retrieve_pair_block(c, ap, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_ap += f12f12_sum_ij;

               ap2a1i1i2_f12_ints->release_pair_block(ap, b, f12_idx);
               a1ap2i1i2_f12_ints->release_pair_block(c, ap, f12_idx);
             }
             *iter_array = f12f12_sum_ap;
             ++iter_array;
           }
         }

         iter_array = RR3;
         for (int c = 0; c < nvir; ++c) {
           for(int b = 0; b < nvir; ++b) {

             double f12f12_sum_ap = 0;
             for(int ap = 0; ap < ncabs2; ++ap) {
               const double* f12_ij_blk1 = a1ap2i1i2_f12_ints->retrieve_pair_block(b, ap, f12_idx);
               const double* f12_ij_blk2 = ap2a1i1i2_f12_ints->retrieve_pair_block(ap, c, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_ap += f12f12_sum_ij;

               a1ap2i1i2_f12_ints->release_pair_block(b, ap, f12_idx);
               ap2a1i1i2_f12_ints->release_pair_block(ap, c, f12_idx);
             }
             *iter_array = f12f12_sum_ap;
             ++iter_array;
           }
         }

         iter_array = RR4;
         for (int c = 0; c < nvir; ++c) {
           for(int b = 0; b < nvir; ++b) {

             double f12f12_sum_ap = 0;
             for(int ap = 0; ap < ncabs2; ++ap) {
               const double* f12_ij_blk1 = ap2a1i1i2_f12_ints->retrieve_pair_block(ap, b, f12_idx);
               const double* f12_ij_blk2 = ap2a1i1i2_f12_ints->retrieve_pair_block(ap, c, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_ap += f12f12_sum_ij;

               ap2a1i1i2_f12_ints->release_pair_block(ap, b, f12_idx);
               ap2a1i1i2_f12_ints->release_pair_block(ap, c, f12_idx);
             }
             *iter_array = f12f12_sum_ap;
             ++iter_array;
           }
         }
         a1ap2i1i2_f12_ints->deactivate();
         ap2a1i1i2_f12_ints->deactivate();
       } else {
           Ref<DistArray4> a2ap1i1i2_f12_ints;
           Ref<DistArray4> ap1a2i1i2_f12_ints;
           activate_ints(vir2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, a2ap1i1i2_f12_ints);
           activate_ints(cabs1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, ap1a2i1i2_f12_ints);

           compute_RR2_sum_3idx(RR31_32, f12_idx,
                                ap1a2i1i2_f12_ints, ap1a2i1i2_f12_ints,
                                RR1); // R^a'b_ij R^ij_a'c
           compute_RR2_sum_3idx(RR13_32, f12_idx,
                                a2ap1i1i2_f12_ints, ap1a2i1i2_f12_ints,
                                RR2); // R^ba'_ij R^ij_a'c
           compute_RR2_sum_3idx(RR31_23, f12_idx,
                                ap1a2i1i2_f12_ints, a2ap1i1i2_f12_ints,
                                RR3); // R^a'b_ij R^ij_ca'
           compute_RR2_sum_3idx(RR13_23, f12_idx,
                                a2ap1i1i2_f12_ints, a2ap1i1i2_f12_ints,
                                RR4); // R^ba'_ij R^ij_ca'
           a2ap1i1i2_f12_ints->deactivate();
           ap1a2i1i2_f12_ints->deactivate();
       }

       if (debug_ >= DefaultPrintThresholds::mostN2) {
         const string spin_label = (spin == Alpha? "Alpha" : "Beta");
         const string RR1_label = (spin == Alpha? "R^ba'_ij R^ij_ca'" : "R^a'b_ij R^ij_a'c");
         const string RR2_label = (spin == Alpha? "R^a'b_ij R^ij_ca'" : "R^ba'_ij R^ij_a'c");
         const string RR3_label = (spin == Alpha? "R^ba'_ij R^ij_a'c" : "R^a'b_ij R^ij_ca'");
         const string RR4_label = (spin == Alpha? "R^a'b_ij R^ij_a'c" : "R^ba'_ij R^ij_ca'");

         print_intermediate(spin_label, RR1_label, RR1, nvir, nvir);
         print_intermediate(spin_label, RR2_label, RR2, nvir, nvir);
         print_intermediate(spin_label, RR3_label, RR3, nvir, nvir);
         print_intermediate(spin_label, RR4_label, RR4, nvir, nvir);
       }

     } else if (nspincases2 == 2) {
         RR1 = RRapb_apc;
         RR2 = RRbap_apc;
         RR3 = RRbap_apc;
         RR4 = RRapb_apc;
    }// end of AlphaBeta part for D^c_b

     // calculate D^c_b[s]
     double* Dcb = new double[nvir12];
     for (int c = 0;  c < nvir; ++c) {
       for (int b = 0; b < nvir; ++b) {

         const int idx = c * nvir + b;
         // AlphaAlpha/BetaBeta part
         double dbc_12 = C_1 * C_1 * (RRapb_apc[idx] - RRbap_apc[idx]);

         // AlphaBeta part
         if (nocc_alpha != 0 && nocc_beta != 0) {

           dbc_12 += pow(0.5 * (C_0 + C_1), 2) * RR1[idx]
                  + 0.25 * (C_0 * C_0 - C_1 * C_1) * (RR2[idx] + RR3[idx])
                  + pow(0.5 * (C_0 - C_1), 2) * RR4[idx];
         } // end of AlphaBeta part

         Dcb[idx] = dbc_12;
       } // end of looping over b
     } // end of calculating D^c_b

      delete[] RRapb_apc;
      delete[] RRbap_apc;
      if (nspincases2 == 3) {
        delete[] RR1;
        delete[] RR2;
        delete[] RR3;
        delete[] RR4;
      }

      print_intermediate(spinletters, " test D^c_b", Dcb, nvir, nvir);
      delete[] Dcb;
   } // end of spincase1 loop
}
// end of computation of D^c_b

// test function for D^c'_b' sum over a
void MP2R12Energy_Diag::compute_Dcpbp_a(const int nspincases1, const int nspincases2,
                                       const double C_0, const double C_1)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);

  // obtain occ, vir, orbs, and cabs orbitals for AlphaBeta case
  const SpinCase1 spin1 = case1(AlphaBeta);
  const SpinCase1 spin2 = case2(AlphaBeta);

  const Ref<OrbitalSpace> occ1_act = r12eval()->occ_act(spin1);
  const Ref<OrbitalSpace> vir1 = r12eval()->vir(spin1);
  const Ref<OrbitalSpace> orbs1 = r12eval()->orbs(spin1);
  const Ref<OrbitalSpace> cabs1 = r12eval()->cabs_space_hcanonical(spin1);

  const Ref<OrbitalSpace> occ2_act = r12eval()->occ_act(spin2);
  const Ref<OrbitalSpace> vir2 = r12eval()->vir(spin2);
  const Ref<OrbitalSpace> orbs2 = r12eval()->orbs(spin2);
  const Ref<OrbitalSpace> cabs2 = r12eval()->cabs_space_hcanonical(spin2);

  const int nvir1 = vir1->rank();
  const int nvir2 = vir2->rank();
  const int ncabs1 = cabs1->rank();
  const int ncabs2 = cabs2->rank();

  const int nocc_alpha = occ1_act->rank();
  const int nocc_beta = occ2_act->rank();

  if (debug_ >= DefaultPrintThresholds::N2)
    ExEnv::out0() << endl << "Number of alpha virtual orbital: " << nvir1
                  << endl << "Number of beta virtual orbital: " << nvir2
                  << endl << "number of alpha cabs orbital: " << ncabs1
                  << endl << "number of beta cabs orbital: " << ncabs2 << endl;

  for(int s = 0; s < nspincases1; ++s) {
     const SpinCase1 spin = static_cast<SpinCase1>(s);

     // spin=0 (Alpha) => AlphaAlpha case (1)
     // spin=1 (Beta) => BetaBeta case (2)
     const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
     string spinletters = to_string(spincase);

     // obtain orbitals
     const Ref<OrbitalSpace> occ_act = r12eval()->occ_act(spin);
     const Ref<OrbitalSpace> vir = r12eval()->vir(spin);
     const Ref<OrbitalSpace> cabs = r12eval()->cabs_space_hcanonical(spin);

     const int nocc_act = occ_act->rank();

     // skip spincase if no electron of this kind
     if (nocc_act == 0)
       continue;

     const int nvir = vir->rank();
     const int ncabs = cabs->rank();
     if (debug_ >= DefaultPrintThresholds::N2)
       ExEnv::out0() << endl << spinletters << " D^c'_b' sum over a" << endl
                     << "number of occupied orbital: " << nocc_act << endl
                     << "number of virtual orbital: " << nvir << endl
                     << "number of cabs orbital: " << ncabs << endl;

     // calculate AlphaAlpha/BetaBeta part for D^c'_b':
     // R^ab'_ij R^ij_ac' & R^b'a_ij R^ij_ac' contracted over i,j,a
     const int ncabs12 = ncabs * ncabs;
     double* RRabp_acp = new double[ncabs12];
     double* RRbpa_acp = new double[ncabs12];
     fill_n(RRabp_acp, ncabs12, 0.0);
     fill_n(RRbpa_acp, ncabs12, 0.0);

     // activate integrals
     Ref<DistArray4> aapii_f12_ints;
     activate_ints(vir->id(), cabs->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, aapii_f12_ints);

     Ref<DistArray4> apaii_f12_ints;
     activate_ints(cabs->id(), vir->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, apaii_f12_ints);

     // R^ab'_ij R^ij_ac'
     double* iter_array = RRabp_acp;
     const blasint one = 1;
     const blasint nocc12 = nocc_act * nocc_act;
     for (int cp = 0; cp < ncabs; ++cp) {
       for(int bp = 0; bp < ncabs; ++bp) {

         double f12f12_sum_a = 0;
         for(int a = 0; a < nvir; ++a) {
           const double* f12_ij_blk1 = aapii_f12_ints->retrieve_pair_block(a, bp, f12_idx);
           const double* f12_ij_blk2 = aapii_f12_ints->retrieve_pair_block(a, cp, f12_idx);

           const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

           f12f12_sum_a += f12f12_sum_ij;

           aapii_f12_ints->release_pair_block(a, bp, f12_idx);
           aapii_f12_ints->release_pair_block(a, cp, f12_idx);
         }
         *iter_array = f12f12_sum_a;
         ++iter_array;
       }
     }

     // R^b'a_ij R^ij_ac'
     iter_array = RRbpa_acp;
     for (int cp = 0; cp < ncabs; ++cp) {
       for(int bp = 0; bp < ncabs; ++bp) {

         double f12f12_sum_a = 0;
         for(int a = 0; a < nvir; ++a) {
           const double* f12_ij_blk1 = apaii_f12_ints->retrieve_pair_block(bp, a, f12_idx);
           const double* f12_ij_blk2 = aapii_f12_ints->retrieve_pair_block(a, cp, f12_idx);

           const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

           f12f12_sum_a += f12f12_sum_ij;

           apaii_f12_ints->release_pair_block(bp, a, f12_idx);
           aapii_f12_ints->release_pair_block(a, cp, f12_idx);
         }
         *iter_array = f12f12_sum_a;
         ++iter_array;
       }
     }

     aapii_f12_ints->deactivate();
     apaii_f12_ints->deactivate();

     if (debug_ >= DefaultPrintThresholds::mostN2) {
       print_intermediate(spinletters, "R^ab'_ij R^ij_ac'", RRabp_acp, ncabs, ncabs);
       print_intermediate(spinletters, "R^b'a_ij R^ij_ac'", RRbpa_acp, ncabs, ncabs);
     }

     // calculate AlphaBeta part for D^c'_b':
     //            alpha                      beta
     // RR1:   R^b'a_ij R^ij_c'a   or   R^a'b_ij R^ij_a'c
     // RR2:   R^ab'_ij R^ij_c'a   or   R^ba'_ij R^ij_a'c
     // RR3:   R^b'a_ij R^ij_ac'   or   R^a'b_ij R^ij_ca'
     // RR4:   R^ab'_ij R^ij_ac'   or   R^ba'_ij R^ij_ca'
     double* RR1 = NULL;
     double* RR2 = NULL;
     double* RR3 = NULL;
     double* RR4 = NULL;

     if (nspincases2 == 3) {

       RR1 = new double[ncabs12];
       RR2 = new double[ncabs12];
       RR3 = new double[ncabs12];
       RR4 = new double[ncabs12];
       fill_n(RR1, ncabs12, 0.0);
       fill_n(RR2, ncabs12, 0.0);
       fill_n(RR3, ncabs12, 0.0);
       fill_n(RR4, ncabs12, 0.0);

       const blasint nocc12 = nocc_alpha * nocc_beta;

       if (spin == Alpha) {
         Ref<DistArray4> ap1a2i1i2_f12_ints;
         Ref<DistArray4> a2ap1i1i2_f12_ints;
         activate_ints(cabs1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, ap1a2i1i2_f12_ints);
         activate_ints(vir2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, a2ap1i1i2_f12_ints);

         iter_array = RR1;  // R^b'a_ij R^ij_c'a
         for (int cp = 0; cp < ncabs; ++cp) {
           for(int bp = 0; bp < ncabs; ++bp) {

             double f12f12_sum_a = 0;
             for(int a = 0; a < nvir2; ++a) {
               const double* f12_ij_blk1 = ap1a2i1i2_f12_ints->retrieve_pair_block(bp, a, f12_idx);
               const double* f12_ij_blk2 = ap1a2i1i2_f12_ints->retrieve_pair_block(cp, a, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_a += f12f12_sum_ij;

               ap1a2i1i2_f12_ints->release_pair_block(bp, a, f12_idx);
               ap1a2i1i2_f12_ints->release_pair_block(cp, a, f12_idx);
             }
             *iter_array = f12f12_sum_a;
             ++iter_array;
           }
         }

         iter_array = RR2; // R^ab'_ij R^ij_c'a
         for (int cp = 0; cp < ncabs; ++cp) {
           for(int bp = 0; bp < ncabs; ++bp) {

             double f12f12_sum_a = 0;
             for(int a = 0; a < nvir2; ++a) {
               const double* f12_ij_blk1 = a2ap1i1i2_f12_ints->retrieve_pair_block(a, bp, f12_idx);
               const double* f12_ij_blk2 = ap1a2i1i2_f12_ints->retrieve_pair_block(cp, a, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_a += f12f12_sum_ij;

               a2ap1i1i2_f12_ints->release_pair_block(a, bp, f12_idx);
               ap1a2i1i2_f12_ints->release_pair_block(cp, a, f12_idx);
             }
             *iter_array = f12f12_sum_a;
             ++iter_array;
           }
         }

         iter_array = RR3;  // R^b'a_ij R^ij_ac'
         for (int cp = 0; cp < ncabs; ++cp) {
           for(int bp = 0; bp < ncabs; ++bp) {

             double f12f12_sum_a = 0;
             for(int a = 0; a < nvir2; ++a) {
               const double* f12_ij_blk1 = ap1a2i1i2_f12_ints->retrieve_pair_block(bp, a, f12_idx);
               const double* f12_ij_blk2 = a2ap1i1i2_f12_ints->retrieve_pair_block(a, cp, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_a += f12f12_sum_ij;

               ap1a2i1i2_f12_ints->release_pair_block(bp, a, f12_idx);
               a2ap1i1i2_f12_ints->release_pair_block(a, cp, f12_idx);
             }
             *iter_array = f12f12_sum_a;
             ++iter_array;
           }
         }

         iter_array = RR4; // R^ab'_ij R^ij_ac'
         for (int cp = 0; cp < ncabs; ++cp) {
           for(int bp = 0; bp < ncabs; ++bp) {

             double f12f12_sum_a = 0;
             for(int a = 0; a < nvir2; ++a) {
               const double* f12_ij_blk1 = a2ap1i1i2_f12_ints->retrieve_pair_block(a, bp, f12_idx);
               const double* f12_ij_blk2 = a2ap1i1i2_f12_ints->retrieve_pair_block(a, cp, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_a += f12f12_sum_ij;

               a2ap1i1i2_f12_ints->release_pair_block(a, bp, f12_idx);
               a2ap1i1i2_f12_ints->release_pair_block(a, cp, f12_idx);
             }
             *iter_array = f12f12_sum_a;
             ++iter_array;
           }
         }
         ap1a2i1i2_f12_ints->deactivate();
         a2ap1i1i2_f12_ints->deactivate();
       } else {
           Ref<DistArray4> a1ap2i1i2_f12_ints;
           Ref<DistArray4> ap2a1i1i2_f12_ints;
           activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, a1ap2i1i2_f12_ints);
           activate_ints(cabs2->id(), vir1->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, ap2a1i1i2_f12_ints);

           iter_array = RR4;  // R^b'a_ij R^ij_c'a
           for (int cp = 0; cp < ncabs; ++cp) {
             for(int bp = 0; bp < ncabs; ++bp) {

               double f12f12_sum_a = 0;
               for(int a = 0; a < nvir1; ++a) {
                 const double* f12_ij_blk1 = ap2a1i1i2_f12_ints->retrieve_pair_block(bp, a, f12_idx);
                 const double* f12_ij_blk2 = ap2a1i1i2_f12_ints->retrieve_pair_block(cp, a, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_a += f12f12_sum_ij;

                 ap2a1i1i2_f12_ints->release_pair_block(bp, a, f12_idx);
                 ap2a1i1i2_f12_ints->release_pair_block(cp, a, f12_idx);
               }
               *iter_array = f12f12_sum_a;
               ++iter_array;
             }
           }

           iter_array = RR3;  // R^ab'_ij R^ij_c'a
           for (int cp = 0; cp < ncabs; ++cp) {
             for(int bp = 0; bp < ncabs; ++bp) {

               double f12f12_sum_a = 0;
               for(int a = 0; a < nvir1; ++a) {
                 const double* f12_ij_blk1 = a1ap2i1i2_f12_ints->retrieve_pair_block(a, bp, f12_idx);
                 const double* f12_ij_blk2 = ap2a1i1i2_f12_ints->retrieve_pair_block(cp, a, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_a += f12f12_sum_ij;

                 a1ap2i1i2_f12_ints->release_pair_block(a, bp, f12_idx);
                 ap2a1i1i2_f12_ints->release_pair_block(cp, a, f12_idx);
               }
               *iter_array = f12f12_sum_a;
               ++iter_array;
             }
           }

           iter_array = RR2;  // R^b'a_ij R^ij_ac'
           for (int cp = 0; cp < ncabs; ++cp) {
             for(int bp = 0; bp < ncabs; ++bp) {

               double f12f12_sum_a = 0;
               for(int a = 0; a < nvir1; ++a) {
                 const double* f12_ij_blk1 = ap2a1i1i2_f12_ints->retrieve_pair_block(bp, a, f12_idx);
                 const double* f12_ij_blk2 = a1ap2i1i2_f12_ints->retrieve_pair_block(a, cp, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_a += f12f12_sum_ij;

                 ap2a1i1i2_f12_ints->release_pair_block(bp, a, f12_idx);
                 a1ap2i1i2_f12_ints->release_pair_block(a, cp, f12_idx);
               }
               *iter_array = f12f12_sum_a;
               ++iter_array;
             }
           }

           iter_array = RR1;  // R^ab'_ij R^ij_ac'
           for (int cp = 0; cp < ncabs; ++cp) {
             for(int bp = 0; bp < ncabs; ++bp) {

               double f12f12_sum_a = 0;
               for(int a = 0; a < nvir1; ++a) {
                 const double* f12_ij_blk1 = a1ap2i1i2_f12_ints->retrieve_pair_block(a, bp, f12_idx);
                 const double* f12_ij_blk2 = a1ap2i1i2_f12_ints->retrieve_pair_block(a, cp, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_a += f12f12_sum_ij;

                 a1ap2i1i2_f12_ints->release_pair_block(a, bp, f12_idx);
                 a1ap2i1i2_f12_ints->release_pair_block(a, cp, f12_idx);
               }
               *iter_array = f12f12_sum_a;
               ++iter_array;
             }
           }

           a1ap2i1i2_f12_ints->deactivate();
           ap2a1i1i2_f12_ints->deactivate();
       }

       if (debug_ >= DefaultPrintThresholds::mostN2) {
         const string spin_label = (spin == Alpha? "Alpha" : "Beta");
         const string RR1_label = (spin == Alpha? "R^b'a_ij R^ij_c'a" : "R^ab'_ij R^ij_ac'");
         const string RR2_label = (spin == Alpha? "R^ab'_ij R^ij_c'a" : "R^b'a_ij R^ij_ac'");
         const string RR3_label = (spin == Alpha? "R^b'a_ij R^ij_ac'" : "R^ab'_ij R^ij_c'a");
         const string RR4_label = (spin == Alpha? "R^ab'_ij R^ij_ac'" : "R^b'a_ij R^ij_c'a");

         print_intermediate(spin_label, RR1_label, RR1, ncabs, ncabs);
         print_intermediate(spin_label, RR2_label, RR2, ncabs, ncabs);
         print_intermediate(spin_label, RR3_label, RR3, ncabs, ncabs);
         print_intermediate(spin_label, RR4_label, RR4, ncabs, ncabs);
       }

     } else if (nspincases2 == 2) {
         RR1 = RRabp_acp;
         RR2 = RRbpa_acp;
         RR3 = RRbpa_acp;
         RR4 = RRabp_acp;
    }// end of AlphaBeta part for D^c'_b'

     // calculate D^c'_b'
     double* Dcpbp = new double[ncabs12];
     fill_n(Dcpbp, ncabs12, 0.0);

     for (int cp = 0;  cp < ncabs; ++cp) {
       for (int bp = 0; bp < ncabs; ++bp) {

         const int idx = cp * ncabs + bp;
         // AlphaAlpha/BetaBeta part
         double dbpcp_12 = C_1 * C_1 * (RRabp_acp[idx] - RRbpa_acp[idx]);

         // AlphaBeta part
         if (nocc_alpha != 0 && nocc_beta != 0) {

           dbpcp_12 += pow(0.5 * (C_0 + C_1), 2) * RR1[idx]
                  + 0.25 * (C_0 * C_0 - C_1 * C_1) * (RR2[idx] + RR3[idx])
                  + pow(0.5 * (C_0 - C_1), 2) * RR4[idx];
         } // end of AlphaBeta part

         Dcpbp[idx] = dbpcp_12;
       } // end of looping over b'
     } // end of calculating D^c'_b'

      delete[] RRabp_acp;
      delete[] RRbpa_acp;
      if (nspincases2 == 3) {
        delete[] RR1;
        delete[] RR2;
        delete[] RR3;
        delete[] RR4;
      }

      print_intermediate(spinletters, " test D^c'_b' sum over a", Dcpbp, ncabs, ncabs);
      delete[] Dcpbp;
      Dcpbp = NULL;
   } // end of spincase1 loop
}
// end of computation of D^c'_b' sum over a

// test function for D^c'_b' sum over a'
void MP2R12Energy_Diag::compute_Dcpbp_ap(const int nspincases1, const int nspincases2,
                                        const double C_0, const double C_1)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);

  // obtain occ, vir, orbs, and cabs orbitals for AlphaBeta case
  const SpinCase1 spin1 = case1(AlphaBeta);
  const SpinCase1 spin2 = case2(AlphaBeta);

  const Ref<OrbitalSpace> occ1_act = r12eval()->occ_act(spin1);
  const Ref<OrbitalSpace> vir1 = r12eval()->vir(spin1);
  const Ref<OrbitalSpace> orbs1 = r12eval()->orbs(spin1);
  // cabs_space_hcanonical only used for test propose
  const Ref<OrbitalSpace> cabs1 = r12eval()->cabs_space_hcanonical(spin1);

  const Ref<OrbitalSpace> occ2_act = r12eval()->occ_act(spin2);
  const Ref<OrbitalSpace> vir2 = r12eval()->vir(spin2);
  const Ref<OrbitalSpace> orbs2 = r12eval()->orbs(spin2);
  // cabs_space_hcanonical only used for test propose
  const Ref<OrbitalSpace> cabs2 = r12eval()->cabs_space_hcanonical(spin2);

  const int nvir1 = vir1->rank();
  const int nvir2 = vir2->rank();
  const int ncabs1 = cabs1->rank();
  const int ncabs2 = cabs2->rank();

  const int nocc_alpha = occ1_act->rank();
  const int nocc_beta = occ2_act->rank();

  if (debug_ >= DefaultPrintThresholds::N2)
    ExEnv::out0() << endl << "Number of alpha virtual orbital: " << nvir1
                  << endl << "Number of beta virtual orbital: " << nvir2
                  << endl << "number of alpha cabs orbital: " << ncabs1
                  << endl << "number of beta cabs orbital: " << ncabs2 << endl;

  for(int s = 0; s < nspincases1; ++s) {
     const SpinCase1 spin = static_cast<SpinCase1>(s);

     // spin=0 (Alpha) => AlphaAlpha case (1)
     // spin=1 (Beta) => BetaBeta case (2)
     const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
     string spinletters = to_string(spincase);

     // obtain orbitals
     const Ref<OrbitalSpace> occ_act = r12eval()->occ_act(spin);
     const Ref<OrbitalSpace> vir = r12eval()->vir(spin);
     const Ref<OrbitalSpace> cabs = r12eval()->cabs_space_hcanonical(spin);

     const int nocc_act = occ_act->rank();

     // skip spincase if no electron of this kind
     if (nocc_act == 0)
       continue;

     const int nvir = vir->rank();
     const int ncabs = cabs->rank();
     if (debug_ >= DefaultPrintThresholds::N2)
       ExEnv::out0() << endl << spinletters << " D^c'_b' sum over a'" << endl
                     << "number of occupied orbital: " << nocc_act << endl
                     << "number of virtual orbital: " << nvir << endl
                     << "number of cabs orbital: " << ncabs << endl;

     // calculate AlphaAlpha/BetaBeta part for D^c'_b':
     // R^a'b'_ij R^ij_a'c' & R^b'a'_ij R^ij_a'c' contracted over i,j,a'
     const int ncabs12 = ncabs * ncabs;
     double* RRapbp_apcp = new double[ncabs12];
     double* RRbpap_apcp = new double[ncabs12];
     fill_n(RRapbp_apcp, ncabs12, 0.0);
     fill_n(RRbpap_apcp, ncabs12, 0.0);

     // activate integrals
     Ref<DistArray4> apapii_f12_ints;
     activate_ints(cabs->id(), cabs->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, apapii_f12_ints);

     // R^a'b'_ij R^ij_a'c'
     double* iter_array = RRapbp_apcp;
     const blasint one = 1;
     const blasint nocc12 = nocc_act * nocc_act;
     for (int cp = 0; cp < ncabs; ++cp) {
       for(int bp = 0; bp < ncabs; ++bp) {

         double f12f12_sum_ap = 0;
         for(int ap = 0; ap < ncabs; ++ap) {
           const double* f12_ij_blk1 = apapii_f12_ints->retrieve_pair_block(ap, bp, f12_idx);
           const double* f12_ij_blk2 = apapii_f12_ints->retrieve_pair_block(ap, cp, f12_idx);

           const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

           f12f12_sum_ap += f12f12_sum_ij;

           apapii_f12_ints->release_pair_block(ap, bp, f12_idx);
           apapii_f12_ints->release_pair_block(ap, cp, f12_idx);
         }
         *iter_array = f12f12_sum_ap;
         ++iter_array;
       }
     }

     // R^b'a'_ij R^ij_a'c'
     iter_array = RRbpap_apcp;
     for (int cp = 0; cp < ncabs; ++cp) {
       for(int bp = 0; bp < ncabs; ++bp) {

         double f12f12_sum_ap = 0;
         for(int ap = 0; ap < ncabs; ++ap) {
           const double* f12_ij_blk1 = apapii_f12_ints->retrieve_pair_block(bp, ap, f12_idx);
           const double* f12_ij_blk2 = apapii_f12_ints->retrieve_pair_block(ap, cp, f12_idx);

           const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

           f12f12_sum_ap += f12f12_sum_ij;

           apapii_f12_ints->release_pair_block(bp, ap, f12_idx);
           apapii_f12_ints->release_pair_block(ap, cp, f12_idx);
         }
         *iter_array = f12f12_sum_ap;
         ++iter_array;
       }
     }

     apapii_f12_ints->deactivate();

     if (debug_ >= DefaultPrintThresholds::mostN2) {
       print_intermediate(spinletters, "R^a'b'_ij R^ij_a'c'", RRapbp_apcp, ncabs, ncabs);
       print_intermediate(spinletters, "R^b'a'_ij R^ij_a'c'", RRbpap_apcp, ncabs, ncabs);
     }

     // calculate AlphaBeta part for D^c'_b':
     //            alpha                      beta
     // RR1:   R^b'a'_ij R^ij_c'a'   or   R^a'b_ij R^ij_a'c
     // RR2:   R^a'b'_ij R^ij_c'a'   or   R^ba'_ij R^ij_a'c
     // RR3:   R^b'a'_ij R^ij_a'c'   or   R^a'b_ij R^ij_ca'
     // RR4:   R^a'b'_ij R^ij_a'c'   or   R^ba'_ij R^ij_ca'
     double* RR1 = NULL;
     double* RR2 = NULL;
     double* RR3 = NULL;
     double* RR4 = NULL;

     if (nspincases2 == 3) {

       RR1 = new double[ncabs12];
       RR2 = new double[ncabs12];
       RR3 = new double[ncabs12];
       RR4 = new double[ncabs12];
       fill_n(RR1, ncabs12, 0.0);
       fill_n(RR2, ncabs12, 0.0);
       fill_n(RR3, ncabs12, 0.0);
       fill_n(RR4, ncabs12, 0.0);

       const blasint nocc12 = nocc_alpha * nocc_beta;
       Ref<DistArray4> ap1ap2i1i2_f12_ints;
       Ref<DistArray4> ap2ap1i1i2_f12_ints;
       activate_ints(cabs1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, ap1ap2i1i2_f12_ints);
       activate_ints(cabs2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, ap2ap1i1i2_f12_ints);

       if (spin == Alpha) {
         iter_array = RR1;  // R^b'a_'ij R^ij_c'a'
         for (int cp = 0; cp < ncabs; ++cp) {
           for(int bp = 0; bp < ncabs; ++bp) {

             double f12f12_sum_ap = 0;
             for(int ap = 0; ap < ncabs2; ++ap) {
               const double* f12_ij_blk1 = ap1ap2i1i2_f12_ints->retrieve_pair_block(bp, ap, f12_idx);
               const double* f12_ij_blk2 = ap1ap2i1i2_f12_ints->retrieve_pair_block(cp, ap, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_ap += f12f12_sum_ij;

               ap1ap2i1i2_f12_ints->release_pair_block(bp, ap, f12_idx);
               ap1ap2i1i2_f12_ints->release_pair_block(cp, ap, f12_idx);
             }
             *iter_array = f12f12_sum_ap;
             ++iter_array;
           }
         }

         iter_array = RR2; // R^a'b'_ij R^ij_c'a'
         for (int cp = 0; cp < ncabs; ++cp) {
           for(int bp = 0; bp < ncabs; ++bp) {

             double f12f12_sum_ap = 0;
             for(int ap = 0; ap < ncabs2; ++ap) {
               const double* f12_ij_blk1 = ap2ap1i1i2_f12_ints->retrieve_pair_block(ap, bp, f12_idx);
               const double* f12_ij_blk2 = ap1ap2i1i2_f12_ints->retrieve_pair_block(cp, ap, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_ap += f12f12_sum_ij;

               ap2ap1i1i2_f12_ints->release_pair_block(ap, bp, f12_idx);
               ap1ap2i1i2_f12_ints->release_pair_block(cp, ap, f12_idx);
             }
             *iter_array = f12f12_sum_ap;
             ++iter_array;
           }
         }

         iter_array = RR3;  // R^b'a'_ij R^ij_a'c'
         for (int cp = 0; cp < ncabs; ++cp) {
           for(int bp = 0; bp < ncabs; ++bp) {

             double f12f12_sum_ap = 0;
             for(int ap = 0; ap < ncabs2; ++ap) {
               const double* f12_ij_blk1 = ap1ap2i1i2_f12_ints->retrieve_pair_block(bp, ap, f12_idx);
               const double* f12_ij_blk2 = ap2ap1i1i2_f12_ints->retrieve_pair_block(ap, cp, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_ap += f12f12_sum_ij;

               ap1ap2i1i2_f12_ints->release_pair_block(bp, ap, f12_idx);
               ap2ap1i1i2_f12_ints->release_pair_block(ap, cp, f12_idx);
             }
             *iter_array = f12f12_sum_ap;
             ++iter_array;
           }
         }

         iter_array = RR4; // R^a'b'_ij R^ij_a'c'
         for (int cp = 0; cp < ncabs; ++cp) {
           for(int bp = 0; bp < ncabs; ++bp) {

             double f12f12_sum_ap = 0;
             for(int ap = 0; ap < ncabs2; ++ap) {
               const double* f12_ij_blk1 = ap2ap1i1i2_f12_ints->retrieve_pair_block(ap, bp, f12_idx);
               const double* f12_ij_blk2 = ap2ap1i1i2_f12_ints->retrieve_pair_block(ap, cp, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_ap += f12f12_sum_ij;

               ap2ap1i1i2_f12_ints->release_pair_block(ap, bp, f12_idx);
               ap2ap1i1i2_f12_ints->release_pair_block(ap, cp, f12_idx);
             }
             *iter_array = f12f12_sum_ap;
             ++iter_array;
           }
         }

       } else {
           iter_array = RR4;  // R^b'a'_ij R^ij_c'a'
           for (int cp = 0; cp < ncabs; ++cp) {
             for(int bp = 0; bp < ncabs; ++bp) {

               double f12f12_sum_ap = 0;
               for(int ap = 0; ap < ncabs1; ++ap) {
                 const double* f12_ij_blk1 = ap2ap1i1i2_f12_ints->retrieve_pair_block(bp, ap, f12_idx);
                 const double* f12_ij_blk2 = ap2ap1i1i2_f12_ints->retrieve_pair_block(cp, ap, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_ap += f12f12_sum_ij;

                 ap2ap1i1i2_f12_ints->release_pair_block(bp, ap, f12_idx);
                 ap2ap1i1i2_f12_ints->release_pair_block(cp, ap, f12_idx);
               }
               *iter_array = f12f12_sum_ap;
               ++iter_array;
             }
           }

           iter_array = RR3;  // R^a'b'_ij R^ij_c'a'
           for (int cp = 0; cp < ncabs; ++cp) {
             for(int bp = 0; bp < ncabs; ++bp) {

               double f12f12_sum_ap = 0;
               for(int ap = 0; ap < ncabs1; ++ap) {
                 const double* f12_ij_blk1 = ap1ap2i1i2_f12_ints->retrieve_pair_block(ap, bp, f12_idx);
                 const double* f12_ij_blk2 = ap2ap1i1i2_f12_ints->retrieve_pair_block(cp, ap, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_ap += f12f12_sum_ij;

                 ap1ap2i1i2_f12_ints->release_pair_block(ap, bp, f12_idx);
                 ap2ap1i1i2_f12_ints->release_pair_block(cp, ap, f12_idx);
               }
               *iter_array = f12f12_sum_ap;
               ++iter_array;
             }
           }

           iter_array = RR2;  // R^b'a'_ij R^ij_a'c'
           for (int cp = 0; cp < ncabs; ++cp) {
             for(int bp = 0; bp < ncabs; ++bp) {

               double f12f12_sum_ap = 0;
               for(int ap = 0; ap < ncabs1; ++ap) {
                 const double* f12_ij_blk1 = ap2ap1i1i2_f12_ints->retrieve_pair_block(bp, ap, f12_idx);
                 const double* f12_ij_blk2 = ap1ap2i1i2_f12_ints->retrieve_pair_block(ap, cp, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_ap += f12f12_sum_ij;

                 ap2ap1i1i2_f12_ints->release_pair_block(bp, ap, f12_idx);
                 ap1ap2i1i2_f12_ints->release_pair_block(ap, cp, f12_idx);
               }
               *iter_array = f12f12_sum_ap;
               ++iter_array;
             }
           }

           iter_array = RR1;  // R^a'b'_ij R^ij_a'c'
           for (int cp = 0; cp < ncabs; ++cp) {
             for(int bp = 0; bp < ncabs; ++bp) {

               double f12f12_sum_ap = 0;
               for(int ap = 0; ap < ncabs1; ++ap) {
                 const double* f12_ij_blk1 = ap1ap2i1i2_f12_ints->retrieve_pair_block(ap, bp, f12_idx);
                 const double* f12_ij_blk2 = ap1ap2i1i2_f12_ints->retrieve_pair_block(ap, cp, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_ap += f12f12_sum_ij;

                 ap1ap2i1i2_f12_ints->release_pair_block(ap, bp, f12_idx);
                 ap1ap2i1i2_f12_ints->release_pair_block(ap, cp, f12_idx);
               }
               *iter_array = f12f12_sum_ap;
               ++iter_array;
             }
           }

       }
       ap1ap2i1i2_f12_ints->deactivate();
       ap2ap1i1i2_f12_ints->deactivate();

       if (debug_ >= DefaultPrintThresholds::mostN2) {
         const string spin_label = (spin == Alpha? "Alpha" : "Beta");
         const string RR1_label = (spin == Alpha? "R^b'a'_ij R^ij_c'a'" : "R^a'b'_ij R^ij_a'c'");
         const string RR2_label = (spin == Alpha? "R^a'b'_ij R^ij_c'a'" : "R^b'a'_ij R^ij_a'c'");
         const string RR3_label = (spin == Alpha? "R^b'a'_ij R^ij_a'c'" : "R^a'b'_ij R^ij_c'a'");
         const string RR4_label = (spin == Alpha? "R^a'b'_ij R^ij_a'c'" : "R^b'a'_ij R^ij_c'a'");

         print_intermediate(spin_label, RR1_label, RR1, ncabs, ncabs);
         print_intermediate(spin_label, RR2_label, RR2, ncabs, ncabs);
         print_intermediate(spin_label, RR3_label, RR3, ncabs, ncabs);
         print_intermediate(spin_label, RR4_label, RR4, ncabs, ncabs);
       }

     } else if (nspincases2 == 2) {
         RR1 = RRapbp_apcp;
         RR2 = RRbpap_apcp;
         RR3 = RRbpap_apcp;
         RR4 = RRapbp_apcp;
    }// end of AlphaBeta part for D^c'_b' sum over a'

     // calculate D^c'_b'sum over a'
     double* Dcpbp = new double[ncabs12];
     fill_n(Dcpbp, ncabs12, 0.0);

     for (int cp = 0;  cp < ncabs; ++cp) {
       for (int bp = 0; bp < ncabs; ++bp) {

         const int idx = cp * ncabs + bp;
         // AlphaAlpha/BetaBeta part
         double dbpcp_12 = C_1 * C_1 * (RRapbp_apcp[idx] - RRbpap_apcp[idx]);

         // AlphaBeta part
         if (nocc_alpha != 0 && nocc_beta != 0) {

           dbpcp_12 += pow(0.5 * (C_0 + C_1), 2) * RR1[idx]
                  + 0.25 * (C_0 * C_0 - C_1 * C_1) * (RR2[idx] + RR3[idx])
                  + pow(0.5 * (C_0 - C_1), 2) * RR4[idx];
         } // end of AlphaBeta part

         Dcpbp[idx] = dbpcp_12;
       } // end of looping over b'
     } // end of calculating D^c'_b' sum over a'

      delete[] RRapbp_apcp;
      delete[] RRbpap_apcp;
      if (nspincases2 == 3) {
        delete[] RR1;
        delete[] RR2;
        delete[] RR3;
        delete[] RR4;
      }

      print_intermediate(spinletters, " test D^c'_b' sum over a'", Dcpbp, ncabs, ncabs);
      delete[] Dcpbp;
      Dcpbp = NULL;
   } // end of spincase1 loop
}
// end of computation of D^c'_b' sum over a'

// test function for D^a'_a RR part
void MP2R12Energy_Diag::compute_Dapa_RR(const int nspincases1, const int nspincases2,
                                        const double C_0, const double C_1)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);

  // obtain occ, vir, orbs, and cabs orbitals for AlphaBeta case
  const SpinCase1 spin1 = case1(AlphaBeta);
  const SpinCase1 spin2 = case2(AlphaBeta);

  const Ref<OrbitalSpace> occ1_act = r12eval()->occ_act(spin1);
  const Ref<OrbitalSpace> vir1 = r12eval()->vir(spin1);
  const Ref<OrbitalSpace> orbs1 = r12eval()->orbs(spin1);
  //const Ref<OrbitalSpace> cabs1 = r12eval()->cabs_space_hcanonical(spin1);
  const Ref<OrbitalSpace> cabs1 = r12world->cabs_space(spin1);

  const Ref<OrbitalSpace> occ2_act = r12eval()->occ_act(spin2);
  const Ref<OrbitalSpace> vir2 = r12eval()->vir(spin2);
  const Ref<OrbitalSpace> orbs2 = r12eval()->orbs(spin2);
  //const Ref<OrbitalSpace> cabs2 = r12eval()->cabs_space_hcanonical(spin2);
  const Ref<OrbitalSpace> cabs2 = r12world->cabs_space(spin2);

  const int nvir1 = vir1->rank();
  const int nvir2 = vir2->rank();
  const int ncabs1 = cabs1->rank();
  const int ncabs2 = cabs2->rank();

  const int nocc_alpha = occ1_act->rank();
  const int nocc_beta = occ2_act->rank();

  if (debug_ >= DefaultPrintThresholds::N2)
    ExEnv::out0() << endl << "Number of alpha virtual orbital: " << nvir1
                  << endl << "Number of beta virtual orbital: " << nvir2
                  << endl << "number of alpha cabs orbital: " << ncabs1
                  << endl << "number of beta cabs orbital: " << ncabs2 << endl;

  for(int s = 0; s < nspincases1; ++s) {
     const SpinCase1 spin = static_cast<SpinCase1>(s);

     // spin=0 (Alpha) => AlphaAlpha case (1)
     // spin=1 (Beta) => BetaBeta case (2)
     const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
     string spinletters = to_string(spincase);

     // obtain orbitals
     const Ref<OrbitalSpace> occ_act = r12eval()->occ_act(spin);
     const Ref<OrbitalSpace> vir = r12eval()->vir(spin);
     //const Ref<OrbitalSpace> cabs = r12eval()->cabs_space_hcanonical(spin);
     const Ref<OrbitalSpace> cabs = r12world->cabs_space(spin);

     const int nocc_act = occ_act->rank();

     // skip spincase if no electron of this kind
     if (nocc_act == 0)
       continue;

     const int nvir = vir->rank();
     const int ncabs = cabs->rank();
     if (debug_ >= DefaultPrintThresholds::N2)
       ExEnv::out0() << endl << spinletters << " D^a'_a RR part" << endl
                     << "number of occupied orbital: " << nocc_act << endl
                     << "number of virtual orbital: " << nvir << endl
                     << "number of cabs orbital: " << ncabs << endl;

     // calculate AlphaAlpha/BetaBeta part for D^a'_a:
     // R^ab'_ij R^ij_a'b' & R^b'a_ij R^ij_a'b' contracted over i,j,b'
     const int ncabsvir =  ncabs * nvir;
     double* RRabp_apbp = new double[ncabsvir];
     double* RRbpa_apbp = new double[ncabsvir];
     fill_n(RRabp_apbp, ncabsvir, 0.0);
     fill_n(RRbpa_apbp, ncabsvir, 0.0);

     // activate integrals
     Ref<DistArray4> aapii_f12_ints;
     activate_ints(vir->id(), cabs->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, aapii_f12_ints);
     Ref<DistArray4> apaii_f12_ints;
     activate_ints(cabs->id(), vir->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, apaii_f12_ints);
     Ref<DistArray4> apapii_f12_ints;
     activate_ints(cabs->id(), cabs->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, apapii_f12_ints);

     // R^ab'_ij R^ij_a'b'
     double* iter_array = RRabp_apbp;
     const blasint one = 1;
     const blasint nocc12 = nocc_act * nocc_act;
     for (int ap = 0; ap < ncabs; ++ap) {
       for(int a = 0; a < nvir; ++a) {

         double f12f12_sum_bp = 0;
         for(int bp = 0; bp < ncabs; ++bp) {
           const double* f12_ij_blk1 = aapii_f12_ints->retrieve_pair_block(a, bp, f12_idx);
           const double* f12_ij_blk2 = apapii_f12_ints->retrieve_pair_block(ap, bp, f12_idx);

           const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

           f12f12_sum_bp += f12f12_sum_ij;

           aapii_f12_ints->release_pair_block(a, bp, f12_idx);
           apapii_f12_ints->release_pair_block(ap, bp, f12_idx);
         }
         *iter_array = f12f12_sum_bp;
         ++iter_array;
       }
     }

     // test R^ab'_ij R^ij_a'b' ?= R^b'a'_ij R^ij_b'a
//     // R^ab'_ij R^ij_a'b'
//     double* RRbpap_bpa = new double[ncabsvir];
//     fill_n(RRbpap_bpa, ncabsvir, 0.0);
//     iter_array = RRbpap_bpa;
//     for (int ap = 0; ap < ncabs; ++ap) {
//       for(int a = 0; a < nvir; ++a) {
//
//         double f12f12_sum_bp = 0;
//         for(int bp = 0; bp < ncabs; ++bp) {
//           const double* f12_ij_blk1 = apaii_f12_ints->retrieve_pair_block(bp, a, f12_idx);
//           const double* f12_ij_blk2 = apapii_f12_ints->retrieve_pair_block(bp, ap, f12_idx);
//
//           const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);
//
//           f12f12_sum_bp += f12f12_sum_ij;
//
//           apaii_f12_ints->release_pair_block(bp, a, f12_idx);
//           apapii_f12_ints->release_pair_block(bp, ap, f12_idx);
//         }
//         *iter_array = f12f12_sum_bp;
//         ++iter_array;
//       }
//     }
//     print_intermediate(spinletters, "R^b'a'_ij R^ij_b'a", RRbpap_bpa, ncabs, nvir);

     // R^b'a_ij R^ij_a'b'
     iter_array = RRbpa_apbp;
     for (int ap = 0; ap < ncabs; ++ap) {
       for(int a = 0; a < nvir; ++a) {

         double f12f12_sum_bp = 0;
         for(int bp = 0; bp < ncabs; ++bp) {
           const double* f12_ij_blk1 = apaii_f12_ints->retrieve_pair_block(bp, a, f12_idx);
           const double* f12_ij_blk2 = apapii_f12_ints->retrieve_pair_block(ap, bp, f12_idx);

           const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

           f12f12_sum_bp += f12f12_sum_ij;

           apaii_f12_ints->release_pair_block(bp, a, f12_idx);
           apapii_f12_ints->release_pair_block(ap, bp, f12_idx);
         }
         *iter_array = f12f12_sum_bp;
         ++iter_array;
       }
     }

     aapii_f12_ints->deactivate();
     apaii_f12_ints->deactivate();
     apapii_f12_ints->deactivate();

     if (debug_ >= DefaultPrintThresholds::mostN2) {
       print_intermediate(spinletters, "R^ab'_ij R^ij_a'b'", RRabp_apbp, ncabs, nvir);
       print_intermediate(spinletters, "R^b'a_ij R^ij_a'b'", RRbpa_apbp, ncabs, nvir);
     }

     // calculate AlphaBeta part for D^a'_a RR part:
     //            alpha                      beta
     // RR1:   R^ab'_ij R^ij_a'b'   or   R^b'a_ij R^ij_b'a'
     // RR2:   R^b'a_ij R^ij_a'b'   or   R^ab'_ij R^ij_b'a'
     // RR3:   R^ab'_ij R^ij_b'a'   or   R^b'a_ij R^ij_a'b'
     // RR4:   R^b'a_ij R^ij_b'a'   or   R^ab'_ij R^ij_a'b'
     double* RR1 = NULL;
     double* RR2 = NULL;
     double* RR3 = NULL;
     double* RR4 = NULL;

     if (nspincases2 == 3) {

       RR1 = new double[ncabsvir];
       RR2 = new double[ncabsvir];
       RR3 = new double[ncabsvir];
       RR4 = new double[ncabsvir];
       fill_n(RR1, ncabsvir, 0.0);
       fill_n(RR2, ncabsvir, 0.0);
       fill_n(RR3, ncabsvir, 0.0);
       fill_n(RR4, ncabsvir, 0.0);

       const blasint nocc12 = nocc_alpha * nocc_beta;
       Ref<DistArray4> ap1ap2i1i2_f12_ints;
       Ref<DistArray4> ap2ap1i1i2_f12_ints;
       activate_ints(cabs1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, ap1ap2i1i2_f12_ints);
       activate_ints(cabs2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, ap2ap1i1i2_f12_ints);

       if (spin == Alpha) {
         Ref<DistArray4> a1ap2i1i2_f12_ints;
         Ref<DistArray4> ap2a1i1i2_f12_ints;
         activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, a1ap2i1i2_f12_ints);
         activate_ints(cabs2->id(), vir1->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, ap2a1i1i2_f12_ints);

         iter_array = RR1;  // R^ab'_ij R^ij_a'b'
         for (int ap = 0; ap < ncabs; ++ap) {
           for(int a = 0; a < nvir; ++a) {

             double f12f12_sum_bp = 0;
             for(int bp = 0; bp < ncabs2; ++bp) {
               const double* f12_ij_blk1 = a1ap2i1i2_f12_ints->retrieve_pair_block(a, bp, f12_idx);
               const double* f12_ij_blk2 = ap1ap2i1i2_f12_ints->retrieve_pair_block(ap, bp, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_bp += f12f12_sum_ij;

               a1ap2i1i2_f12_ints->release_pair_block(a, bp, f12_idx);
               ap1ap2i1i2_f12_ints->release_pair_block(ap, bp, f12_idx);
             }
             *iter_array = f12f12_sum_bp;
             ++iter_array;
           }
         }

         iter_array = RR2; // R^b'a_ij R^ij_a'b'
         for (int ap = 0; ap < ncabs; ++ap) {
           for(int a = 0; a < nvir; ++a) {

             double f12f12_sum_bp = 0;
             for(int bp = 0; bp < ncabs2; ++bp) {
               const double* f12_ij_blk1 = ap2a1i1i2_f12_ints->retrieve_pair_block(bp, a, f12_idx);
               const double* f12_ij_blk2 = ap1ap2i1i2_f12_ints->retrieve_pair_block(ap, bp, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_bp += f12f12_sum_ij;

               ap2a1i1i2_f12_ints->release_pair_block(bp, a, f12_idx);
               ap1ap2i1i2_f12_ints->release_pair_block(ap, bp, f12_idx);
             }
             *iter_array = f12f12_sum_bp;
             ++iter_array;
           }
         }

         iter_array = RR3;  // R^ab'_ij R^ij_b'a'
         for (int ap = 0; ap < ncabs; ++ap) {
           for(int a = 0; a < nvir; ++a) {

             double f12f12_sum_bp = 0;
             for(int bp = 0; bp < ncabs2; ++bp) {
               const double* f12_ij_blk1 = a1ap2i1i2_f12_ints->retrieve_pair_block(a, bp, f12_idx);
               const double* f12_ij_blk2 = ap2ap1i1i2_f12_ints->retrieve_pair_block(bp, ap, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_bp += f12f12_sum_ij;

               a1ap2i1i2_f12_ints->release_pair_block(a, bp, f12_idx);
               ap2ap1i1i2_f12_ints->release_pair_block(bp, ap, f12_idx);
             }
             *iter_array = f12f12_sum_bp;
             ++iter_array;
           }
         }

         iter_array = RR4; // R^b'a_ij R^ij_b'a'
         for (int ap = 0; ap < ncabs; ++ap) {
           for(int a = 0; a < nvir; ++a) {

             double f12f12_sum_bp = 0;
             for(int bp = 0; bp < ncabs2; ++bp) {
               const double* f12_ij_blk1 = ap2a1i1i2_f12_ints->retrieve_pair_block(bp, a, f12_idx);
               const double* f12_ij_blk2 = ap2ap1i1i2_f12_ints->retrieve_pair_block(bp, ap, f12_idx);

               const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

               f12f12_sum_bp += f12f12_sum_ij;

               ap2a1i1i2_f12_ints->release_pair_block(bp, a, f12_idx);
               ap2ap1i1i2_f12_ints->release_pair_block(bp, ap, f12_idx);
             }
             *iter_array = f12f12_sum_bp;
             ++iter_array;
           }
         }

         a1ap2i1i2_f12_ints->deactivate();
         ap2a1i1i2_f12_ints->deactivate();

       } else {
           Ref<DistArray4> ap1a2i1i2_f12_ints;
           Ref<DistArray4> a2ap1i1i2_f12_ints;
           activate_ints(cabs1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, ap1a2i1i2_f12_ints);
           activate_ints(vir2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, a2ap1i1i2_f12_ints);

           iter_array = RR4;  // R^ab'_ij R^ij_a'b'
           for (int ap = 0; ap < ncabs; ++ap) {
             for(int a = 0; a < nvir; ++a) {

               double f12f12_sum_bp = 0;
               for(int bp = 0; bp < ncabs1; ++bp) {
                 const double* f12_ij_blk1 = a2ap1i1i2_f12_ints->retrieve_pair_block(a, bp, f12_idx);
                 const double* f12_ij_blk2 = ap2ap1i1i2_f12_ints->retrieve_pair_block(ap, bp, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_bp += f12f12_sum_ij;

                 a2ap1i1i2_f12_ints->release_pair_block(a, bp, f12_idx);
                 ap2ap1i1i2_f12_ints->release_pair_block(ap, bp, f12_idx);
               }
               *iter_array = f12f12_sum_bp;
               ++iter_array;
             }
           }

           iter_array = RR3;  // R^b'a_ij R^ij_a'b'
           for (int ap = 0; ap < ncabs; ++ap) {
             for(int a = 0; a < nvir; ++a) {

               double f12f12_sum_bp = 0;
               for(int bp = 0; bp < ncabs1; ++bp) {
                 const double* f12_ij_blk1 = ap1a2i1i2_f12_ints->retrieve_pair_block(bp, a, f12_idx);
                 const double* f12_ij_blk2 = ap2ap1i1i2_f12_ints->retrieve_pair_block(ap, bp, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_bp += f12f12_sum_ij;

                 ap1a2i1i2_f12_ints->release_pair_block(bp, a, f12_idx);
                 ap2ap1i1i2_f12_ints->release_pair_block(ap, bp, f12_idx);
               }
               *iter_array = f12f12_sum_bp;
               ++iter_array;
             }
           }

           iter_array = RR2;  // R^ab'_ij R^ij_b'a'
           for (int ap = 0; ap < ncabs; ++ap) {
             for(int a = 0; a < nvir; ++a) {

               double f12f12_sum_bp = 0;
               for(int bp = 0; bp < ncabs1; ++bp) {
                 const double* f12_ij_blk1 = a2ap1i1i2_f12_ints->retrieve_pair_block(a, bp, f12_idx);
                 const double* f12_ij_blk2 = ap1ap2i1i2_f12_ints->retrieve_pair_block(bp, ap, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_bp += f12f12_sum_ij;

                 a2ap1i1i2_f12_ints->release_pair_block(a, bp, f12_idx);
                 ap1ap2i1i2_f12_ints->release_pair_block(bp, ap, f12_idx);
               }
               *iter_array = f12f12_sum_bp;
               ++iter_array;
             }
           }

           iter_array = RR1;  // R^b'a_ij R^ij_b'a'
           for (int ap = 0; ap < ncabs; ++ap) {
             for(int a = 0; a < nvir; ++a) {

               double f12f12_sum_bp = 0;
               for(int bp = 0; bp < ncabs1; ++bp) {
                 const double* f12_ij_blk1 = ap1a2i1i2_f12_ints->retrieve_pair_block(bp, a, f12_idx);
                 const double* f12_ij_blk2 = ap1ap2i1i2_f12_ints->retrieve_pair_block(bp, ap, f12_idx);

                 const double f12f12_sum_ij = F77_DDOT(&nocc12, f12_ij_blk1, &one, f12_ij_blk2, &one);

                 f12f12_sum_bp += f12f12_sum_ij;

                 ap1a2i1i2_f12_ints->release_pair_block(bp, a, f12_idx);
                 ap1ap2i1i2_f12_ints->release_pair_block(bp, ap, f12_idx);
               }
               *iter_array = f12f12_sum_bp;
               ++iter_array;
             }
           }
           ap1a2i1i2_f12_ints->deactivate();
           a2ap1i1i2_f12_ints->deactivate();
       }
       ap1ap2i1i2_f12_ints->deactivate();
       ap2ap1i1i2_f12_ints->deactivate();

       if (debug_ >= DefaultPrintThresholds::mostN2) {
         const string spin_label = (spin == Alpha? "Alpha" : "Beta");
         const string RR1_label = (spin == Alpha? "R^ab'_ij R^ij_a'b'" : "R^b'a_ij R^ij_b'a'");
         const string RR2_label = (spin == Alpha? "R^b'a_ij R^ij_a'b'" : "R^ab'_ij R^ij_b'a'");
         const string RR3_label = (spin == Alpha? "R^ab'_ij R^ij_b'a'" : "R^b'a_ij R^ij_a'b'");
         const string RR4_label = (spin == Alpha? "R^b'a_ij R^ij_b'a'" : "R^ab'_ij R^ij_a'b'");

         print_intermediate(spin_label, RR1_label, RR1, ncabs, nvir);
         print_intermediate(spin_label, RR2_label, RR2, ncabs, nvir);
         print_intermediate(spin_label, RR3_label, RR3, ncabs, nvir);
         print_intermediate(spin_label, RR4_label, RR4, ncabs, nvir);
       }

     } else if (nspincases2 == 2) {
         RR1 = RRabp_apbp;
         RR2 = RRbpa_apbp;
         RR3 = RRbpa_apbp;
         RR4 = RRabp_apbp;
    }// end of AlphaBeta part for D^a'_a RR part

     // calculate D^a'_a RR part
     double* Dapa_RR = new double[ncabsvir];
     fill_n(Dapa_RR, ncabsvir, 0.0);

     for (int ap = 0;  ap < ncabs; ++ap) {
       for (int a = 0; a < nvir; ++a) {

         const int idx = ap * nvir + a;
         // AlphaAlpha/BetaBeta part
         double dapa_12 = C_1 * C_1 * (RRabp_apbp[idx] - RRbpa_apbp[idx]);

         // AlphaBeta part
         if (nocc_alpha != 0 && nocc_beta != 0) {

           dapa_12 += pow(0.5 * (C_0 + C_1), 2) * RR1[idx]
                  + 0.25 * (C_0 * C_0 - C_1 * C_1) * (RR2[idx] + RR3[idx])
                  + pow(0.5 * (C_0 - C_1), 2) * RR4[idx];
         } // end of AlphaBeta part

         Dapa_RR[idx] = dapa_12;
       } // end of looping over b'
     } // end of calculating D^c'_b' sum over a'

      delete[] RRabp_apbp;
      delete[] RRbpa_apbp;
      if (nspincases2 == 3) {
        delete[] RR1;
        delete[] RR2;
        delete[] RR3;
        delete[] RR4;
      }

      print_intermediate(spinletters, " test D^a'_a RR part", Dapa_RR, ncabs, nvir);
      delete[] Dapa_RR;
      Dapa_RR = NULL;
   } // end of spincase1 loop
}
// end of computation of D^a'_a RR part

// compute R * T2 (CC amplitudes) in D^a'_a
void MP2R12Energy_Diag::compute_RT2_apa(const int nspincases1, const int nspincases2,
                                        const double C_0, const double C_1,
                                        const vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                                        const vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                                        double* D_alpha, double* D_beta)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);

  // Obtain T2 amplitudes
  Ref<DistArray4> T2[NSpinCases2];
  for (int s = 0; s < nspincases2; ++s) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);

    if (r12eval()->dim_oo(spincase2).n() == 0)
      continue;
    T2[spincase2] = r12intermediates_->get_T2_cc(spincase2);

    if (debug_ >= DefaultPrintThresholds::allN2)
      _print(spincase2, T2[spincase2], prepend_spincase(spincase2,"CCSD T2 amplitudes:").c_str());
  } // end of T2

  const Ref<OrbitalSpace> occ1_act = v_orbs1_ab[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2_ab[0];
  const Ref<OrbitalSpace> vir1 = v_orbs1_ab[1];
  const Ref<OrbitalSpace> vir2 = v_orbs2_ab[1];
  const Ref<OrbitalSpace> cabs1 = v_orbs1_ab[3];
  const Ref<OrbitalSpace> cabs2 = v_orbs2_ab[3];

  const int nocc_alpha = occ1_act->rank();
  const int nocc_beta = occ2_act->rank();

  for(int s = 0; s < nspincases1; ++s) {
     const SpinCase1 spin = static_cast<SpinCase1>(s);

     // spin=0 (Alpha) => AlphaAlpha case (1)
     // spin=1 (Beta) => BetaBeta case (2)
     const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
     string spinletters = to_string(spincase);

     // obtain occ_act, vir, orbs, cabs, and occ orbitals in that order
     vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
     vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
     obtain_orbitals(spincase, v_orbs1, v_orbs2);

     const Ref<OrbitalSpace> occ_act = v_orbs1[0];
     const int nocc_act = occ_act->rank();

     // skip spincase if no electron of this kind
     if (nocc_act == 0)
       continue;

     const Ref<OrbitalSpace> vir = v_orbs1[1];
     const Ref<OrbitalSpace> cabs = v_orbs1[3];
     const int nvir = vir->rank();
     const int ncabs = cabs->rank();

     // test propose
     if (debug_ >= DefaultPrintThresholds::mostN2)
       ExEnv::out0() << endl << spinletters << " number of occupied orbital: " << nocc_act << endl
                           << "number of virtual orbital: " << nvir << endl
                           << "number of cabs: " << ncabs << endl;

     // calculate AlphaAlpha/BetaBeta part: RT2^a'b_ij T^ij_ab
     //                                   & RT2^ba'_ij T^ij_ab

     const int ncabs_vir= ncabs * nvir;
     double* RT2apb_ab= new double[ncabs_vir];
     fill_n(RT2apb_ab, ncabs_vir, 0.0);

     // RT2^ba'_ij T^ij_ab = -RT2^a'b_ij T^ij_ab
//     double* RT2bap_ab = new double[ncabs_vir];
//     fill_n(RT2bap_ab, ncabs_vir, 0.0);

     // activate integrals
     Ref<DistArray4> Rij_apb_ints;
     activate_ints(occ_act->id(), occ_act->id(), cabs->id(), vir->id(),
                   descr_f12_key, moints4_rtime, Rij_apb_ints);
     if (!r12eval()->dim_oo(spincase).n() == 0) {
       T2[s+1]->activate();

       // test code
#if 0
     {
       double* const RT2 = new double[ncabs_vir];
       fill_n(RT2, ncabs_vir, 0.0);

       const int nvir12_nocc12 = nvir * nvir * nocc_act * nocc_act;
       double* const T2ab_ij = new double[nvir12_nocc12];
       fill_n(T2ab_ij, nvir12_nocc12, 0.0);
       transform_T2ijab_T2abij(T2[s+1], T2ab_ij);

       Ref<DistArray4> Rapb_ij_ints;
       activate_ints(cabs->id(), vir->id(), occ_act->id(), occ_act->id(),
                     descr_f12_key, moints4_rtime, Rapb_ij_ints);


       // R^ij_a'b * T^ab_ij
       compute_RTmp2_sum_3idx(RT13_23, f12_idx, nvir,
                              Rapb_ij_ints, T2ab_ij,
                              RT2);

       Rapb_ij_ints->deactivate();
       delete[] T2ab_ij;

       print_intermediate(spinletters, "testing R^a'b_ij T2^ij_ab", RT2, ncabs, nvir);

//       fill_n(RT2, ncabs_vir, 0.0);
//       for (int i = 0; i < nocc_act; ++i) {
//         for(int j = 0; j < nocc_act; ++j) {
//           const double* const F12_blk = Rij_apb_ints->retrieve_pair_block(i, j, f12_idx);
//           const double* const T2_blk = T2[s+1]->retrieve_pair_block(i, j, 0);
//
//           for (int b = 0; b < nvir; ++b) {
//             double* iter_RT2 = RT2;
//             for(int ap = 0; ap < ncabs; ++ap) {
//
//               const double* iter_T2_blk = T2_blk + b;
//               for (int a = 0; a < nvir; ++a) {
//
//                 const int F12_blk_idx = ap * nvir + b;
//                 const int T2_blk_idx = a * nvir + b;
//                 const int RT2_idx = ap *nvir + a;
//                 const double RTapa = F12_blk[F12_blk_idx] * (*iter_T2_blk);
////                 RT2[RT2_idx] += RTapa;
//                 *iter_RT2 += RTapa;
//                 ++iter_RT2;
//
//                 iter_T2_blk += nvir;
//               } // end of a loop
//             } // end of a' loop
//           } // end of b loop
//
//           Rij_apb_ints->release_pair_block(i, j, f12_idx);
//           T2[s+1]->release_pair_block(i, j, 0);
//         } // end of j loop
//       } // end of i loop
//
//       print_intermediate(spinletters, "second testing R^a'b_ij T2^ij_ab (T2)", RT2, ncabs, nvir);

       delete[] RT2;
     }
#endif
       // end of test code

       // R^a'b_ij T2^ij_ab (T2: antisymmetrized):
       compute_RT2_sum_3idx(RT13_23, f12_idx,
                            Rij_apb_ints, T2[s+1],
                            RT2apb_ab);
       T2[s+1]->deactivate();
     }

     if (debug_ >= DefaultPrintThresholds::mostN2) {
       print_intermediate(spinletters, "R^a'b_ij T2^ij_ab", RT2apb_ab, ncabs, nvir);
     }

     // calculate AlphaBeta part:
     // RT2^a'b_ij T^ij_ab & RT2^ba'_ij T^ij_ab
     double* RT2_1 = NULL;
     double* RT2_2 = NULL;

     if (nocc_alpha != 0 && nocc_beta != 0) {

       RT2_1 = new double[ncabs_vir];
       RT2_2 = new double[ncabs_vir];
       fill_n(RT2_1, ncabs_vir, 0.0);
       fill_n(RT2_2, ncabs_vir, 0.0);
       T2[AlphaBeta]->activate();

       if (spin == Alpha) {
         // activate integrals
         Ref<DistArray4> Rapb_12_ints;  // R^i(1)j(2)_a'(1)b(2)
         Ref<DistArray4> Rbap_21_ints;  // R^i(1)j(2)_b(2)a'(1)

         if (nspincases2 == 3) {
           activate_ints(occ1_act->id(), occ2_act->id(), cabs1->id(), vir2->id(),
                         descr_f12_key, moints4_rtime, Rapb_12_ints);
         } else {
             Rapb_12_ints = Rij_apb_ints;
         }
         activate_ints(occ1_act->id(), occ2_act->id(), vir2->id(), cabs1->id(),
                       descr_f12_key, moints4_rtime, Rbap_21_ints);

         // R^i(1)j(2)_a'(1)b(2) T2^a(1)b(2)_i(1)j(2):
         compute_RT2_sum_3idx(RT13_23, f12_idx,
                              Rapb_12_ints, T2[AlphaBeta],
                              RT2_1);
         // R^i(1)j(2)_b(2)a'(1) T2^a(1)b(2)_i(1)j(2):
         compute_RT2_sum_3idx(RT31_23, f12_idx,
                              Rbap_21_ints, T2[AlphaBeta],
                              RT2_2);
         // test code
#if 0
         {
           double* const RT2 = new double[ncabs_vir];
           fill_n(RT2, ncabs_vir, 0.0);

           const int nvir1 = vir1->rank();
           const int nvir2 = vir2->rank();
           const int nvir12_nocc12 = nvir1 * nvir2 * nocc_alpha * nocc_beta;
           double* const T2ab_ij = new double[nvir12_nocc12];
           fill_n(T2ab_ij, nvir12_nocc12, 0.0);
           transform_T2ijab_T2abij(T2[AlphaBeta], T2ab_ij);

           Ref<DistArray4> Rapb_ij_ints;
           Ref<DistArray4> Rbap_ij_ints;
           activate_ints(cabs1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, Rapb_ij_ints);
           activate_ints(vir2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, Rbap_ij_ints);


           // R^ij_a'b * T^ab_ij
           compute_RTmp2_sum_3idx(RT13_23, f12_idx, nvir1,
                                  Rapb_ij_ints, T2ab_ij,
                                  RT2);
           print_intermediate("alpha", "testing R^a'b_ij T2^ij_ab", RT2, ncabs, nvir);

           // test compute_RTmp2_sum_3idx function explicitly
#if 0
         {
           const blasint nij = nocc_alpha * nocc_beta;
           const int na = nvir;
           const int nb = nvir2;
           double* const F12T2_1 = new double[ncabs_vir];
//           double* const F12T2_2 = new double[ncabs_vir];
           double* iter_F12T2_1 = F12T2_1;
//           double* iter_F12T2_2 = F12T2_2;
           const blasint one = 1;
           for (int ap = 0; ap < ncabs; ++ap) {
             for(int a = 0; a < na; ++a) {

               double F12T2_sum_idx3_1 = 0;
//               double F12T2_sum_idx3_2 = 0;
               for(int b = 0; b < nb; ++b) {
                 const double* f12_ij_blk1 = Rapb_ij_ints->retrieve_pair_block(ap, b, f12_idx);
//                 const double* f12_ij_blk2 = Rbap_ij_ints->retrieve_pair_block(b, ap, f12_idx);
                 const double* iter_T2 = T2ab_ij + (a * nb + b)* nij;

                 const double F12T2_sum_ij_1 = F77_DDOT(&nij, f12_ij_blk1, &one, iter_T2, &one);
//                 const double F12T2_sum_ij_2 = F77_DDOT(&nij, f12_ij_blk2, &one, iter_T2, &one);

                 F12T2_sum_idx3_1 += F12T2_sum_ij_1;
//                 F12T2_sum_idx3_2 += F12T2_sum_ij_2;

                 Rapb_ij_ints->release_pair_block(ap, b, f12_idx);
//                 Rbap_ij_ints->release_pair_block(b, ap, f12_idx);
               }
               *iter_F12T2_1 = F12T2_sum_idx3_1;
//               *iter_F12T2_2 = F12T2_sum_idx3_2;
               ++iter_F12T2_1;
//               ++iter_F12T2_2;
             }
           }

           print_intermediate("alpha", "testing R^ij_a'b * T^ab_ij for compute_RTmp2_sum_3idx", F12T2_1, ncabs, nvir);
//           print_intermediate(spinletters, "testing R^ij_ba' * T^ab_ij (alpha)", F12T2_2, ncabs, nvir);
           delete[] F12T2_1;
//           delete[] F12T2_2;
         }
#endif

           // test compute_RT2_sum_3idx explicitly
#if 0
           fill_n(RT2, ncabs_vir, 0.0);
           for (int i = 0; i < nocc_alpha; ++i) {
             for(int j = 0; j < nocc_beta; ++j) {
               const double* const F12_blk = Rapb_12_ints->retrieve_pair_block(i, j, f12_idx);
               const double* const T2_blk = T2[AlphaBeta]->retrieve_pair_block(i, j, 0);

               for (int b = 0; b < nvir2; ++b) {

                 for(int ap = 0; ap < ncabs; ++ap) {
                   for (int a = 0; a < nvir; ++a) {

                     const int F12_blk_idx = ap * nvir2 + b;
                     const int T2_blk_idx = a * nvir2 + b;
                     const int RT2_idx = ap *nvir + a;

                     const double RTapa = F12_blk[F12_blk_idx] * T2_blk[T2_blk_idx];
                     RT2[RT2_idx] += RTapa;

                   } // end of a loop
                 } // end of a' loop
               } // end of b loop

               Rapb_12_ints->release_pair_block(i, j, f12_idx);
               T2[AlphaBeta]->release_pair_block(i, j, 0);
             } // end of j loop
           } // end of i loop
           print_intermediate("alpha", "testing R^a'b_ij T2^ij_ab for compute_RT2_sum_3idx", RT2, ncabs, nvir);
#endif

           // R^ij_ba' * T^ab_ij
           compute_RTmp2_sum_3idx(RT31_23, f12_idx, nvir1,
                                  Rbap_ij_ints, T2ab_ij,
                                  RT2);
           print_intermediate("alpha", "testing R^ba'_ij T2^ij_ab", RT2, ncabs, nvir);

           Rapb_ij_ints->deactivate();
           Rbap_ij_ints->deactivate();
           delete[] T2ab_ij;
           delete[] RT2;
         }
#endif
         // end of test code

         if (nspincases2 == 3) Rapb_12_ints->deactivate();
         Rbap_21_ints->deactivate();

       } else {
           // activate integrals
           Ref<DistArray4> Rbap_12_ints;
           Ref<DistArray4> Rapb_21_ints;

           activate_ints(occ1_act->id(), occ2_act->id(), vir1->id(), cabs2->id(),
                         descr_f12_key, moints4_rtime, Rbap_12_ints);

           if (nspincases2 == 3) {
             activate_ints(occ1_act->id(), occ2_act->id(), cabs2->id(), vir1->id(),
                           descr_f12_key, moints4_rtime, Rapb_21_ints);
           } else {
               Rapb_21_ints = Rij_apb_ints;
           }

           // R^ba'_ij T2^ij_ba:
           compute_RT2_sum_3idx(RT31_32, f12_idx,
                                Rbap_12_ints, T2[AlphaBeta],
                                RT2_1);
           // R^a'b_ij T2^ij_ba:
           compute_RT2_sum_3idx(RT13_32, f12_idx,
                                Rapb_21_ints, T2[AlphaBeta],
                                RT2_2);

           // test code
#if 0
           {
             double* const RT2 = new double[ncabs_vir];
             fill_n(RT2, ncabs_vir, 0.0);

             const int nvir1 = vir1->rank();
             const int nvir2 = vir2->rank();
             const int nvir12_nocc12 = nvir1 * nvir2 * nocc_alpha * nocc_beta;
             double* const T2ab_ij = new double[nvir12_nocc12];
             fill_n(T2ab_ij, nvir12_nocc12, 0.0);
             transform_T2ijab_T2abij(T2[AlphaBeta], T2ab_ij);

             Ref<DistArray4> Rapb_ij_ints;
             Ref<DistArray4> Rbap_ij_ints;
             activate_ints(cabs2->id(), vir1->id(), occ1_act->id(), occ2_act->id(),
                           descr_f12_key, moints4_rtime, Rapb_ij_ints);
             activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                           descr_f12_key, moints4_rtime, Rbap_ij_ints);


             // R^ij_ba' * T^ba_ij
             compute_RTmp2_sum_3idx(RT31_32, f12_idx, nvir2,
                                    Rbap_ij_ints, T2ab_ij,
                                    RT2);
             print_intermediate("beta", "testing R^ba'_ij T2^ij_ba", RT2, ncabs, nvir);

             // R^ij_a'b * T^ba_ij
             compute_RTmp2_sum_3idx(RT13_32, f12_idx, nvir2,
                                    Rapb_ij_ints, T2ab_ij,
                                    RT2);
             print_intermediate("beta", "testing R^a'b_ij T2^ij_ba", RT2, ncabs, nvir);

             Rapb_ij_ints->deactivate();
             Rbap_ij_ints->deactivate();
             delete[] T2ab_ij;
             delete[] RT2;
           }
#endif
           // end of test code

           Rbap_12_ints->deactivate();
           if (nspincases2 == 3) Rapb_21_ints->deactivate();
       }
       T2[AlphaBeta]->deactivate();
       Rij_apb_ints->deactivate();

       if (debug_ >= DefaultPrintThresholds::mostN2) {
         const string spin_label = (spin == Alpha? "alpha" : "beta");
         const string RT2_1_label = (spin == Alpha? "RT2^a'b_ij T^ij_ab" : "RT2^ba'_ij T^ij_ba");
         const string RT2_2_label = (spin == Alpha? "RT2^ba'_ij T^ij_ab" : "RT2^a'b_ij T^ij_ba");

         print_intermediate(spin_label, RT2_1_label, RT2_1, ncabs, nvir);
         print_intermediate(spin_label, RT2_2_label, RT2_2, ncabs, nvir);
       }

     } // end of AlphaBeta part for D

     // calculate D^a'_a RT2 part
     const double* iter_RT2apb_ab = RT2apb_ab;
     const double* iter_RT2_1 = RT2_1;
     const double* iter_RT2_2 = RT2_2;

     double* iter_D = (spin == Alpha? D_alpha : D_beta);
     for (int idx1 = 0;  idx1 < ncabs; ++idx1) {

       double d_12 = 0.0;
       for (int idx2 = 0; idx2 < nvir; ++idx2) {

         // AlphaAlpha/BetaBeta part
         d_12 = C_1 * (*iter_RT2apb_ab);
//         ExEnv::out0() << spinletters << " part d^" << idx1 << "_" << idx2 << " = "
//                       << scprintf("%12.10f", d_12) << endl;
         ++iter_RT2apb_ab;

         // AlphaBeta part
         if (nocc_alpha != 0 && nocc_beta != 0) {
   //        ExEnv::out0() << "RR1: "  << scprintf("%12.10f", RR1)
   //                      << "  RR2: " << scprintf("%12.10f", RR2)
   //                      << "  RR3: "  << scprintf("%12.10f", RR3)
   //                      << "  RR4: " << scprintf("%12.10f", RR4)
   //                      << endl;

           d_12 += 0.5 * (C_0 + C_1) * (*iter_RT2_1) + 0.5 * (C_0 - C_1) * (*iter_RT2_2);
//           ExEnv::out0() << "AlphaBeta part: d^"  << idx1 << "_" << idx2 << " = "
//                         << scprintf("%12.10f", d_12) << endl;

           ++iter_RT2_1;
           ++iter_RT2_2;
         } // end of AlphaBeta part

         *iter_D = d_12;
         ++iter_D;
       } // end of looping over c
     } // end of calculating D^b_c[s]

     delete[] RT2apb_ab;
     if (nocc_alpha != 0 && nocc_beta != 0) {
       delete[] RT2_1;
       delete[] RT2_2;
     }
  } // end of spincase1 loop

}
// end of compute_RT2_apa

// compute MP2 T2 amplitude (non-antisymmetrized), stored in array (a,b,i,j)
//  T^i(spin1)j(spin2)_a(spin1)b(spin2) 4-dimension matrix
//= g^ij_ab / (e_i + e_j - e_a - e_b)
void MP2R12Energy_Diag::compute_T2_mp2(const vector< Ref<OrbitalSpace> >& v_orbs1,
                                       const vector< Ref<OrbitalSpace> >& v_orbs2,
                                       double* const T2ab_ij)
{
  // get moints4_rtime, descr_f12_key, and eri_idx
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type eri_type = r12world->r12tech()->corrfactor()->tbint_type_eri();
  const unsigned int eri_idx = descr_f12->intset(eri_type);

  // obtain orbitals
  const Ref<OrbitalSpace>& occ1_act = v_orbs1[0];
  const Ref<OrbitalSpace>& occ2_act = v_orbs2[0];
  const Ref<OrbitalSpace>& vir1 = v_orbs1[1];
  const Ref<OrbitalSpace>& vir2 = v_orbs2[1];
  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();
  const int nvir1 = vir1->rank();
  const int nvir2 = vir2->rank();

  // g^ab_ij
  Ref<DistArray4> a1a2i1i2_ints;
  activate_ints(vir1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                descr_f12_key, moints4_rtime,
                a1a2i1i2_ints);

  // get eigenvalues of Fock matrix
  const RefDiagSCMatrix evals_i1 = occ1_act->evals();
  const RefDiagSCMatrix evals_i2 = occ2_act->evals();
  const RefDiagSCMatrix evals_a1 = vir1->evals();
  const RefDiagSCMatrix evals_a2 = vir2->evals();

  double* iter_T2 = T2ab_ij;
//  if (spincase == AlphaBeta) {
    for (int a = 0; a < nvir1; ++a) {
      for (int b = 0; b < nvir2; ++b) {
        const double* gab_ij = a1a2i1i2_ints->retrieve_pair_block(a, b, eri_idx);

        for (int i = 0; i < nocc1_act; ++i) {
          for (int j = 0; j < nocc2_act; ++j) {
            *iter_T2 =  (*gab_ij)
                       / (evals_i1(i) + evals_i2(j) - evals_a1(a) - evals_a2(b));

            ++gab_ij;
            ++iter_T2;
          }
        }

        a1a2i1i2_ints->release_pair_block(a, b, eri_idx);
      }
    }
//  } else {
//      // AlphaAlpha or BetaBeta case
//      for (int a = 0; a < nvir1; ++a) {
//        for (int b = 0; b < nvir2; ++b) {
//          const double* const gab_ij = a1a2i1i2_ints->retrieve_pair_block(a, b, eri_idx);
//
//          const double* iter_gab_ij = gab_ij;
//          for (int i = 0; i < nocc1_act; ++i) {
//            const double* iter_gab_ji = gab_ij + i;
//
//            for (int j = 0; j < nocc2_act; ++j) {
//              *iter_T2 = (*iter_gab_ij - *iter_gab_ji)
//                        / (evals_i1(i) + evals_i2(j) - evals_a1(a) - evals_a2(b));
//
//               ++iter_gab_ij;
//               iter_gab_ji += nocc1_act;
//               ++iter_T2;
//             }
//           }
//           a1a2i1i2_ints->release_pair_block(a, b, eri_idx);
//         }
//       }
//    }

    a1a2i1i2_ints->deactivate();

    // test code for T2
#if 0
    string spinletters = to_string(spincase);
    Ref<DistArray4> i1i2a1a2_ints;
    activate_ints(occ1_act->id(), occ2_act->id(), vir1->id(),
                  vir2->id(), descr_f12_key, moints4_rtime,
                  i1i2a1a2_ints);

    //   T^ij_ab 4-dimension matrix
    // = g^ij_ab / (e_i + e_j - e_a -e_b)
    const int nocc12 = nocc1_act * nocc2_act;
    const int nvir12 = nvir1 * nvir2;
    double* Tij_ab = new double[nvir12 * nocc12];
    fill_n(Tij_ab, nvir12 * nocc12, 0.0);

    double* iter_Tij_ab = Tij_ab;
    if (spincase == AlphaBeta) {
      for (int i1 = 0; i1 < nocc1_act; ++i1) {
        for (int i2 = 0; i2 < nocc2_act; ++i2) {
          const double* gij_ab = i1i2a1a2_ints->retrieve_pair_block(i1, i2, eri_idx);

          for (int a1 = 0; a1 < nvir1; ++a1) {
            for (int a2 = 0; a2 < nvir2; ++a2, ++iter_Tij_ab, ++gij_ab) {
              *iter_Tij_ab =
                  (*gij_ab)
                      / (evals_i1(i1) + evals_i2(i2) - evals_a1(a1)
                          - evals_a2(a2));
            }
          }
          i1i2a1a2_ints->release_pair_block(i1, i2, eri_idx);
        }
      }
    } else {
        for (int i1 = 0; i1 < nocc1_act; ++i1) {
          for (int i2 = 0; i2 < nocc2_act; ++i2) {
            const double* gij_ab = i1i2a1a2_ints->retrieve_pair_block(i1, i2, eri_idx);
            const double* gji_ab = i1i2a1a2_ints->retrieve_pair_block(i2, i1, eri_idx);

            for (int a1 = 0; a1 < nvir1; ++a1) {
              for (int a2 = 0; a2 < nvir2; ++a2, ++iter_Tij_ab, ++gij_ab, ++gji_ab) {
                *iter_Tij_ab =
                    (*gij_ab - *gji_ab)
                        / (evals_i1(i1) + evals_i2(i2) - evals_a1(a1)
                            - evals_a2(a2));
              }
            }
            i1i2a1a2_ints->release_pair_block(i1, i2, eri_idx);
            i1i2a1a2_ints->release_pair_block(i2, i1, eri_idx);

          }
        }
    }
    i1i2a1a2_ints->deactivate();

    ExEnv::out0() << endl << spinletters << endl
                  << "number of occupied orbital: " << nocc1_act << " " << nocc2_act << endl
                  << "number of virtual orbital: " << nvir1 << " " << nvir2 << endl;
    print_T2ijab_mp2(spinletters, "T2^ij_ab",
                     nocc1_act, nocc2_act, nvir1, nvir2,
                     Tij_ab);
    delete[] Tij_ab;
#endif

}
// end of function: compute_T2_mp2

// compute R * T2  (mp2 amplitudes) in D^a'_a
void MP2R12Energy_Diag::compute_RTmp2_apa(const int nspincases1, const int nspincases2,
                                          const double C_0, const double C_1,
                                          const vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                                          const vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                                          double* const RTmp2_alpha, double* const RTmp2_beta)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);

  // compute T2^ab_ij amplitudes in MP2R12 method
  vector<double*> T2ab_ij(NSpinCases2);
  for(int s = 0; s < nspincases2; ++s) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);

    // obtain occ, vir, orbs, and cabs orbitals in that order
    vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
    vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
    obtain_orbitals(spincase2, v_orbs1, v_orbs2);


    const int nocc1_act = v_orbs1[0]->rank();
    const int nocc2_act = v_orbs2[0]->rank();
    const int nvir1 = v_orbs1[1]->rank();
    const int nvir2 = v_orbs2[1]->rank();
    const int nocc12 = nocc1_act * nocc2_act;
    const int nvir12 = nvir1 * nvir2;

    T2ab_ij[s] = new double[nvir12 * nocc12];
    fill_n(T2ab_ij[s], nvir12 * nocc12, 0.0);

    if (nocc1_act == 0 || nocc2_act == 0)
      continue;
     compute_T2_mp2(v_orbs1, v_orbs2, T2ab_ij[s]);

     // print T2^ab_ij
 #if 0
     string spinletters = to_string(spincase2);
     ExEnv::out0() << endl << spinletters << " MP2 T2^ab_ij" << endl
                   << "number of occupied orbital: " << nocc1_act << " " << nocc2_act << endl
                   << "number of virtual orbital: " << nvir1 << " " << nvir2 << endl;
     print_T2abij_mp2(spinletters, "T2^ab_ij",
                      nocc1_act, nocc2_act, nvir1, nvir2,
                      T2ab_ij[s]);
 #endif
  }

  if (nspincases1 == 1) {
    T2ab_ij[BetaBeta] = T2ab_ij[AlphaAlpha];
  }

  // preliminaries for AlphaBeta
  const Ref<OrbitalSpace> occ1_act = v_orbs1_ab[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2_ab[0];
  const Ref<OrbitalSpace> vir1 = v_orbs1_ab[1];
  const Ref<OrbitalSpace> vir2 = v_orbs2_ab[1];
  const Ref<OrbitalSpace> cabs1 = v_orbs1_ab[3];
  const Ref<OrbitalSpace> cabs2 = v_orbs2_ab[3];

  const int nocc_alpha = occ1_act->rank();
  const int nocc_beta = occ2_act->rank();
  const int nvir_alpha = vir1->rank();
  const int nvir_beta = vir2->rank();

  // Alpha: C1 * (R^IJ_A'B * T^AB_IJ - R^JI_A'B * T^AB_IJ)  (AlphaAlpha part)
  //     + [(C0+C1)/2 * R^IJ_A'B + (C0-C1)/2 * R^JI_A'B] * T^AB_IJ  (AlphaBeta part)
  //
  // Beta:  C1 * (R^IJ_A'B * T^AB_IJ- R^JI_A'B * T^AB_IJ)   (BetaBeta part)
  //     + [(C0+C1)/2 * R^IJ_BA'+ (C0-C1)/2 * R^JI_BA'] * T^BA_IJ  (AlphaBeta part)
  // where T^AB_IJ is not antisymmetrized

  for(int s=0; s<nspincases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);

    // spin=0 (Alpha) => AlphaAlpha case (1)
    // spin=1 (Beta) => BetaBeta case (2)
    const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
    string spinletters = to_string(spincase);

    // obtain occ, vir, orbs, and cabs orbitals in that order
    vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
    vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
    obtain_orbitals(spincase, v_orbs1, v_orbs2);

    Ref<OrbitalSpace> occ_act = v_orbs1[0];
    Ref<OrbitalSpace> vir = v_orbs1[1];
    Ref<OrbitalSpace> cabs = v_orbs1[3];

    const int nocc_act = occ_act->rank();
    // skip spincase if no electron of this kind
    if (nocc_act == 0)
      continue;

     const int nvir = vir->rank();
     const int ncabs = cabs->rank();
     const int ncabs_vir = ncabs * nvir;

     if (debug_ >= DefaultPrintThresholds::N2)
       ExEnv::out0() << endl << spinletters << " D^a'_a RT2 (mp2)" << endl
                     << "number of occupied orbital: " << nocc_act << endl
                     << "number of virtual orbital: " << nvir << endl
                     << "number of cabs orbital: " << ncabs << endl;

     // calculate AlphaAlpha/BetaBeta part for D^a'_a RT2 part:
     // R^ij_a'b * T2^ab_ij & R^ij_ba' * T2^ab_ij
     double* RT2apb_ab= new double[ncabs_vir];
     double* RT2bap_ab = new double[ncabs_vir];
     fill_n(RT2apb_ab, ncabs_vir, 0.0);
     fill_n(RT2bap_ab, ncabs_vir, 0.0);

     // activate integrals
     Ref<DistArray4> Rapb_ij_ints;
     Ref<DistArray4> Rbap_ij_ints;

     activate_ints(cabs->id(), vir->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, Rapb_ij_ints);
     activate_ints(vir->id(), cabs->id(), occ_act->id(), occ_act->id(),
                   descr_f12_key, moints4_rtime, Rbap_ij_ints);

     // R^ij_a'b * T^ab_ij
     compute_RTmp2_sum_3idx(RT13_23, f12_idx, nvir,
                            Rapb_ij_ints, T2ab_ij[spincase],
                            RT2apb_ab);
     // R^ij_ba' * T^ab_ij
     compute_RTmp2_sum_3idx(RT31_23, f12_idx, nvir,
                            Rbap_ij_ints, T2ab_ij[spincase],
                            RT2bap_ab);

     // test R^ij_a'b * T^ab_ij & R^ij_ba' * T^ab_ij
#if 0
     {
       const blasint nij = nocc_act * nocc_act;
       const int na = nvir;
       const int nb = nvir;
       double* const F12T2_1 = new double[ncabs_vir];
       double* const F12T2_2 = new double[ncabs_vir];
       double* iter_F12T2_1 = F12T2_1;
       double* iter_F12T2_2 = F12T2_2;

       const blasint one = 1;
       for (int ap = 0; ap < ncabs; ++ap) {
         for(int a = 0; a < na; ++a) {

           double F12T2_sum_idx3_1 = 0;
           double F12T2_sum_idx3_2 = 0;
           for(int b = 0; b < nb; ++b) {
             const double* f12_ij_blk1 = Rapb_ij_ints->retrieve_pair_block(ap, b, f12_idx);
             const double* f12_ij_blk2 = Rbap_ij_ints->retrieve_pair_block(b, ap, f12_idx);
             const double* iter_T2 = T2ab_ij[spincase] + (a * nb + b)* nij;

             const double F12T2_sum_ij_1 = F77_DDOT(&nij, f12_ij_blk1, &one, iter_T2, &one);
             const double F12T2_sum_ij_2 = F77_DDOT(&nij, f12_ij_blk2, &one, iter_T2, &one);

             F12T2_sum_idx3_1 += F12T2_sum_ij_1;
             F12T2_sum_idx3_2 += F12T2_sum_ij_2;

             Rapb_ij_ints->release_pair_block(ap, b, f12_idx);
             Rbap_ij_ints->release_pair_block(b, ap, f12_idx);
           }
           *iter_F12T2_1 = F12T2_sum_idx3_1;
           *iter_F12T2_2 = F12T2_sum_idx3_2;
           ++iter_F12T2_1;
           ++iter_F12T2_2;
         }
       }

       print_intermediate(spinletters, "testing R^ij_a'b * T^ab_ij", F12T2_1, ncabs, nvir);
       print_intermediate(spinletters, "testing R^ij_ba' * T^ab_ij", F12T2_2, ncabs, nvir);
       delete[] F12T2_1;
       delete[] F12T2_2;
     }
#endif
     // end of testing code

     Rapb_ij_ints->deactivate();
     Rbap_ij_ints->deactivate();

     if (debug_ >= DefaultPrintThresholds::N2){
       print_intermediate(spinletters, "R^ij_a'b * T^ab_ij", RT2apb_ab, ncabs, nvir);
       print_intermediate(spinletters, "R^ij_ba' * T^ab_ij", RT2bap_ab, ncabs, nvir);
     }

     // calculate AlphaBeta for D^a'_a:
     //           alpha case              beta case
     // RT2_1: R^ij_a'b * T^ab_ij     R^ij_ba' * T2^ba_ij
     // RT2_2: R^ji_a'b * T^ab_ij     R^ji_ba' * T2^ba_ij
     double* RT2_1 = NULL;
     double* RT2_2 = NULL;
     if (nspincases2 == 3) {

       RT2_1 = new double[ncabs_vir];
       RT2_2 = new double[ncabs_vir];
       fill_n(RT2_1, ncabs_vir, 0.0);
       fill_n(RT2_2, ncabs_vir, 0.0);

       if (spin == Alpha) {

         // activate integrals
         Ref<DistArray4> Rapb_12_ints;
         Ref<DistArray4> Rbap_21_ints;
         activate_ints(cabs1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, Rapb_12_ints);
         activate_ints(vir2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, Rbap_21_ints);

         // R^i(1)j(2)_a'(1)b(2) * T^a(1)b(2)_i(1)j(2)
         compute_RTmp2_sum_3idx(RT13_23, f12_idx, nvir_alpha,
                                Rapb_12_ints, T2ab_ij[AlphaBeta],
                                RT2_1);
         // R^i(1)j(2)_b(2)a'(1) * T^a(1)b(2)_i(1)j(2)
         compute_RTmp2_sum_3idx(RT31_23, f12_idx, nvir_alpha,
                                Rbap_21_ints, T2ab_ij[AlphaBeta],
                                RT2_2);

         // test R^ij_a'b * T^ab_ij & R^ij_ba' * T^ab_ij
#if 0
         {
           const blasint nij = nocc_alpha * nocc_beta;
           const int na = nvir;
           const int nb = nvir_beta;
           double* const F12T2_1 = new double[ncabs_vir];
           double* const F12T2_2 = new double[ncabs_vir];
           double* iter_F12T2_1 = F12T2_1;
           double* iter_F12T2_2 = F12T2_2;
           const blasint one = 1;
           for (int ap = 0; ap < ncabs; ++ap) {
             for(int a = 0; a < na; ++a) {

               double F12T2_sum_idx3_1 = 0;
               double F12T2_sum_idx3_2 = 0;
               for(int b = 0; b < nb; ++b) {
                 const double* f12_ij_blk1 = Rapb_12_ints->retrieve_pair_block(ap, b, f12_idx);
                 const double* f12_ij_blk2 = Rbap_21_ints->retrieve_pair_block(b, ap, f12_idx);
                 const double* iter_T2 = T2ab_ij[spincase] + (a * nvir_beta + b)* nij;

                 const double F12T2_sum_ij_1 = F77_DDOT(&nij, f12_ij_blk1, &one, iter_T2, &one);
                 const double F12T2_sum_ij_2 = F77_DDOT(&nij, f12_ij_blk2, &one, iter_T2, &one);

                 F12T2_sum_idx3_1 += F12T2_sum_ij_1;
                 F12T2_sum_idx3_2 += F12T2_sum_ij_2;

                 Rapb_12_ints->release_pair_block(ap, b, f12_idx);
                 Rbap_21_ints->release_pair_block(b, ap, f12_idx);
               }
               *iter_F12T2_1 = F12T2_sum_idx3_1;
               *iter_F12T2_2 = F12T2_sum_idx3_2;
               ++iter_F12T2_1;
               ++iter_F12T2_2;
             }
           }

           print_intermediate(spinletters, "testing R^ij_a'b * T^ab_ij (alpha)", F12T2_1, ncabs, nvir);
           print_intermediate(spinletters, "testing R^ij_ba' * T^ab_ij (alpha)", F12T2_2, ncabs, nvir);
           delete[] F12T2_1;
           delete[] F12T2_2;
         }
#endif
         // end of testing code

         Rapb_12_ints->deactivate();
         Rbap_21_ints->deactivate();

       } else {
           // activate integrals
           Ref<DistArray4> Rbap_12_ints;
           Ref<DistArray4> Rapb_21_ints;
           activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, Rbap_12_ints);
           activate_ints(cabs2->id(), vir1->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, Rapb_21_ints);

           // R^i(1)j(2)_b(1)a'(2) * T^b(1)a(2)_i(1)j(2)
           compute_RTmp2_sum_3idx(RT31_32, f12_idx, nvir_beta,
                                  Rbap_12_ints, T2ab_ij[AlphaBeta],
                                  RT2_1);
           // R^i(1)j(2)_a'(2)b(1) * T^b(1)a(2)_i(1)j(2)
           compute_RTmp2_sum_3idx(RT13_32, f12_idx, nvir_beta,
                                  Rapb_21_ints, T2ab_ij[AlphaBeta],
                                  RT2_2);

           // test R^ij_ba' * T^ba_ij & R^ij_a'b * T^ba_ij
#if 0
           {
             const blasint nij = nocc_alpha * nocc_beta;
             const int na = nvir;
             const int nb = nvir_alpha;
             double* const F12T2_1 = new double[ncabs_vir];
             double* const F12T2_2 = new double[ncabs_vir];
             double* iter_F12T2_1 = F12T2_1;
             double* iter_F12T2_2 = F12T2_2;

             const blasint one = 1;
             for (int ap = 0; ap < ncabs; ++ap) {
               for(int a = 0; a < na; ++a) {

                 double F12T2_sum_idx3_1 = 0;
                 double F12T2_sum_idx3_2 = 0;
                 for(int b = 0; b < nb; ++b) {
                   const double* f12_ij_blk1 = Rapb_21_ints->retrieve_pair_block(ap, b, f12_idx);   // R^ij_a'b
                   const double* f12_ij_blk2 = Rbap_12_ints->retrieve_pair_block(b, ap, f12_idx);   // R^ij_ba'
                   const double* iter_T2 = T2ab_ij[AlphaBeta] + (b * nvir_beta + a)* nij;                  // T^ba_ij

                   const double F12T2_sum_ij_1 = F77_DDOT(&nij, f12_ij_blk1, &one, iter_T2, &one);
                   const double F12T2_sum_ij_2 = F77_DDOT(&nij, f12_ij_blk2, &one, iter_T2, &one);

                   F12T2_sum_idx3_1 += F12T2_sum_ij_1;
                   F12T2_sum_idx3_2 += F12T2_sum_ij_2;

                   Rapb_21_ints->release_pair_block(ap, b, f12_idx);
                   Rbap_12_ints->release_pair_block(b, ap, f12_idx);
                 }
                 *iter_F12T2_1 = F12T2_sum_idx3_1;
                 *iter_F12T2_2 = F12T2_sum_idx3_2;
                 ++iter_F12T2_1;
                 ++iter_F12T2_2;
               }
             }

             print_intermediate(spinletters, "testing R^ij_ba' * T^ba_ij (beta)", F12T2_2, ncabs, nvir);
             print_intermediate(spinletters, "testing R^ij_a'b * T^ba_ij (beta)", F12T2_1, ncabs, nvir);
             delete[] F12T2_1;
             delete[] F12T2_2;
           }
#endif
           // end of testing code

           Rbap_12_ints->deactivate();
           Rapb_21_ints->deactivate();
       }

       if (debug_ >= DefaultPrintThresholds::N2) {
         const string spin_label = (spin == Alpha? "alpha" : "beta");
         const string RT2_1_label = (spin == Alpha? "R^ij_a'b * T^ab_ij" : "R^ij_ba' * T^ba_ij");
         const string RT2_2_label = (spin == Alpha? "R^ij_ba' * T^ab_ij" : "R^ij_a'b * T^ba_ij");

         print_intermediate(spin_label, RT2_1_label, RT2_1, ncabs, nvir);
         print_intermediate(spin_label, RT2_2_label, RT2_2, ncabs, nvir);
       }

     } else {
         RT2_1 = RT2apb_ab;
         RT2_2 = RT2bap_ab;
    } // end of AlphaBeta part for D^a'_a

     delete[] T2ab_ij[spincase];
     T2ab_ij[spincase] = NULL;

     // calculate D^a'_a RT2 part
     const double* iter_RT2apb_ab = RT2apb_ab;
     const double* iter_RT2bap_ab = RT2bap_ab;
     const double* iter_RT2_1 = RT2_1;
     const double* iter_RT2_2 = RT2_2;
     double* iter_RTmp2 = (spin == Alpha? RTmp2_alpha : RTmp2_beta);

     for (int idx1 = 0;  idx1 < ncabs; ++idx1) {
       for (int idx2 = 0; idx2 < nvir; ++idx2) {

         // AlphaAlpha/BetaBeta part
         double d_12 = C_1 * (*iter_RT2apb_ab - *iter_RT2bap_ab);
//         ExEnv::out0() << spinletters << " part d^" << idx1 << "_" << idx2 << " = "
//                       << scprintf("%12.10f", d_12) << endl;

         ++iter_RT2apb_ab;
         ++iter_RT2bap_ab;

         // AlphaBeta part
         if (nocc_alpha != 0 && nocc_beta != 0) {
   //        ExEnv::out0() << "RR1: "  << scprintf("%12.10f", RR1)
   //                      << "  RR2: " << scprintf("%12.10f", RR2)
   //                      << "  RR3: "  << scprintf("%12.10f", RR3)
   //                      << "  RR4: " << scprintf("%12.10f", RR4)
   //                      << endl;

           d_12 += 0.5 * (C_0 + C_1) * (*iter_RT2_1) + 0.5 * (C_0 - C_1) * (*iter_RT2_2);
//           ExEnv::out0() << "AlphaBeta part: d^"  << idx1 << "_" << idx2 << " = "
//                         << scprintf("%12.10f", d_12) << endl;

           ++iter_RT2_1;
           ++iter_RT2_2;
         } // end of AlphaBeta part

         *iter_RTmp2 = d_12;
         ++iter_RTmp2;
       } // end of looping over c
     } // end of calculating D^a'_a element

     //
     const double* const iter_D = (spin == Alpha? RTmp2_alpha : RTmp2_beta);
     print_intermediate(spinletters, "D^A'_A RT2 (mp2) part", iter_D, ncabs, nvir);

     delete[] RT2apb_ab;
     delete[] RT2bap_ab;
     if (nspincases2 == 3) {
       delete[] RT2_1;
       delete[] RT2_2;
     }

   } // end of loop spincase1

  // delete T2ab_ij
  if (nspincases2 > 1 ) {
    delete[] T2ab_ij[AlphaBeta];
    T2ab_ij[AlphaBeta] = NULL;
  }
}
// end of compute_RTmp2_apa

void MP2R12Energy_Diag::compute_density_diag()
{

  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
  Ref<MessageGrp> msg = r12world->world()->msg();
  int me = msg->me();
  int ntasks = msg->n();

  // Only diagonal ansatz is supported
  const bool diag = r12world->r12tech()->ansatz()->diag();
  if (diag == false)
    throw ProgrammingError("only diagonal ansatz supported", __FILE__,
                           __LINE__);
  // geminal coefficients
  const double C_0 = 1.0 / 2.0;
  const double C_1 = 1.0 / 4.0;

  // get the number of the unique spincase2 and spincase1
  const int nspincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  if (debug_ >= DefaultPrintThresholds::N2)
    ExEnv::out0() <<endl << "Number of unique spincases = " << nspincases2 << endl;
  const int nspincases1 = r12eval()->nspincases1();

  // obtain AlphaBeta occ, vir, orbs, cabs, and occ orbitals
  vector< Ref<OrbitalSpace> > v_orbs1_ab;
  vector< Ref<OrbitalSpace> > v_orbs2_ab;
  obtain_orbitals(AlphaBeta, v_orbs1_ab, v_orbs2_ab);

  // get the numbers of the occupied alpha and beta orbitals
  const Ref<OrbitalSpace> occ1_act = v_orbs1_ab[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2_ab[0];
  const Ref<OrbitalSpace> vir1 = v_orbs1_ab[1];
  const Ref<OrbitalSpace> vir2 = v_orbs2_ab[1];
  const Ref<OrbitalSpace> orbs1 = v_orbs1_ab[2];
  const Ref<OrbitalSpace> orbs2 = v_orbs2_ab[2];
  const Ref<OrbitalSpace> cabs1 = v_orbs1_ab[3];
  const Ref<OrbitalSpace> cabs2 = v_orbs2_ab[3];
  const Ref<OrbitalSpace> occ1 = v_orbs1_ab[4];
  const Ref<OrbitalSpace> occ2 = v_orbs2_ab[4];
  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();
  const int nvir1 = vir1->rank();
  const int nvir2 = vir2->rank();
  const int ncabs1 = cabs1->rank();
  const int ncabs2 = cabs2->rank();

  const int norb = nocc1_act + nvir1;
  const int norb_com = nocc1_act + nvir1 + ncabs1;

  ExEnv::out0() << endl << "**********************************************" << endl;
  ExEnv::out0() << endl << "Computing CCSD_F12 one-particle density" << endl<< endl;

  // test propose
//  if (debug_ >= DefaultPrintThresholds::N2)
    ExEnv::out0() << endl << "number of alpha active occupied orbital: " << nocc1_act << endl
                          << "number of beta active occupied orbital: " << nocc2_act << endl
                          << "number of alpha virtual orbital: " << nvir1 << endl
                          << "number of beta virtual orbital: " << nvir2 << endl
                          << "number of alpha cabs: " << ncabs1 << endl
                          << "number of beta cabs: " << ncabs2 << endl;

  // initialize arrays
  const int nocc11 = nocc1_act * nocc1_act;
  const int nvir11 = nvir1 * nvir1;
  const int ncabs11 = ncabs1 * ncabs1;
  const int ncabs1_vir1 = ncabs1 * nvir1;

  double* Dm_i_alpha = new double[nocc11];
  double* Dc_b_alpha = new double[nvir11];
  double* Dcp_bp_alpha_A = new double[ncabs11];
  double* Dcp_bp_alpha_Ap = new double[ncabs11];
  double* Dap_a_alpha_RR = new double[ncabs1_vir1];

  fill_n(Dm_i_alpha, nocc11, 0.0);
  fill_n(Dc_b_alpha, nvir11, 0.0);
  fill_n(Dcp_bp_alpha_A, ncabs11, 0.0);
  fill_n(Dcp_bp_alpha_Ap, ncabs11, 0.0);
  fill_n(Dap_a_alpha_RR, ncabs1_vir1, 0.0);

  double* Dm_i_beta;
  double* Dc_b_beta;
  double* Dcp_bp_beta_A;
  double* Dcp_bp_beta_Ap;
  double* Dap_a_beta_RR;
  if (nspincases1 == 2) {
    const int nocc22 = nocc2_act * nocc2_act;
    const int nvir22 = nvir2 * nvir2;
    const int ncabs22 = ncabs2 * ncabs2;
    const int ncabs2_vir2 = ncabs2 * nvir2;

    Dm_i_beta = new double[nocc22];
    Dc_b_beta = new double[nvir22];
    Dcp_bp_beta_A = new double[ncabs22];
    Dcp_bp_beta_Ap = new double[ncabs22];
    Dap_a_beta_RR = new double[ncabs2_vir2];

    fill_n(Dm_i_beta, nocc22, 0.0);
    fill_n(Dc_b_beta, nvir22, 0.0);
    fill_n(Dcp_bp_beta_A, ncabs22, 0.0);
    fill_n(Dcp_bp_beta_Ap, ncabs22, 0.0);
    fill_n(Dap_a_beta_RR, ncabs2_vir2, 0.0);
  } else if (nspincases1 == 1) {
      Dm_i_beta =  Dm_i_alpha;
      Dc_b_beta = Dc_b_alpha;
      Dcp_bp_beta_A = Dcp_bp_alpha_A;
      Dcp_bp_beta_Ap = Dcp_bp_alpha_Ap;
      Dap_a_beta_RR = Dap_a_alpha_RR;
  }

  //
  // D^m_i
  // Alpha d^M_I = C1^2 * (R^IJ_AB R^AB_MJ - R^JI_AB R^AB_MJ)  (I, J, M, A, B in alpha orbitals)
  //
  //             + [(C0+C1)/2]^2 * R^IJ_AB R^AB_MJ   (I, M, A in alpha orbitals, J, B inbeta orbitals)
  //             + (C0+C1)/2*(C0-C1)/2 * (R^JI_AB R^AB_MJ + R^IJ_AB R^AB_JM)
  //             + [(C0-C1)/2]^2 * R^JI_AB R^AB_JM
  //
  // Beta  d^M_I = C1^2 * (R^IJ_AB R^AB_MJ - R^JI_AB R^AB_MJ)  (I, J, A, B in Beta orbitals)
  //
  //             + [(C0+C1)/2]^2 * R^JI_AB R^AB_JM    (J, A in alpha orbitals, I, B in beta orbitals)
  //             + (C0+C1)/2*(C0-C1)/2 * (R^IJ_AB R^AB_JM + R^IJ_AB R^AB_MJ)
  //             + [(C0-C1)/2]^2 * R^IJ_AB R^AB_MJ

#if 0
  // test for compute_Dmi
  ExEnv::out0() << endl << "testing: D^i_i" << endl;
  compute_Dii_test(nspincases1, nspincases2, nocc1_act, nocc2_act, C_0, C_1);
#endif

  // compute D^m_i through R^ij_a'b' R_mj^a'b' + 2 R^ij_ab' R_mj^ab'
  //
  // R^ij_a'b' R_mj^a'b'
  // Alpha:  C1^2 * (R^IJ_A'B' R^A'B'_MJ - R^JI_A'B' R^A'B'_MJ) (I, J, M, A', B' in alpha orbitals)
  //         + [(C0+C1)/2]^2 * R^IJ_A'B' R^A'B_MJ  (I, M, A' in alpha orbitals, J, B' in beta orbitals)
  //         + (C0+C1)/2*(C0-C1)/2 * (R^JI_A'B' R^A'B'_MJ + R^IJ_A'B' R^A'B_JM)
  //         + [(C0-C1)/2]^2 * R^JI_A'B' R^A'B'_JM
  //
  // Beta    C1^2 * (R^IJ_A'B' R^A'B'_MJ - R^JI_A'B' R^A'B'_MJ) (I, J, M, A', B' in beta orbitals)
  //         + [(C0+C1)/2]^2 * R^IJ_A'B' R^A'B_MJ  (I, M, B' in beta orbitals, J, A' in alpha orbitals)
  //         + (C0+C1)/2*(C0-C1)/2 * (R^JI_A'B' R^A'B'_MJ + R^IJ_A'B' R^A'B_JM)
  //         + [(C0-C1)/2]^2 * R^JI_A'B' R^A'B'_JM
  //
  // 2 R^ij_ab' R_mj^ab'
  // Alpha: C1^2 * (R^IJ_AB' R^AB'_MJ - R^JI_AB' R^AB'_MJ
  //               - R^IJ_AB' R^AB'_JM - R^JI_AB' R^AB'_JM) (I, J, M, A, B' in alpha orbitals)
  //        + [(C0+C1)/2]^2 * R^IJ_AB' R^AB_MJ           (I, M, A in alpha orbitals, J, B' in beta orbitals)
  //        + (C0+C1)/2*(C0-C1)/2 * (R^JI_AB' R^AB'_MJ + R^IJ_AB' R^AB_JM)
  //        + [(C0-C1)/2]^2 * R^JI_AB' R^AB'_JM
  //        + [(C0+C1)/2]^2 * R^IJ_B'A R^B'A_MJ           (I, M, B' in alpha orbitals, J, A in beta orbitals)
  //        + (C0+C1)/2*(C0-C1)/2 * (R^JI_B'A R^B'A_MJ + R^IJ_B'A R^B'A_JM)
  //        + [(C0-C1)/2]^2 * R^JI_B'A R^B'A_JM
  //
  // Beta:  C1^2 * (R^IJ_AB' R^AB'_MJ - R^JI_AB' R^AB'_MJ
  //               - R^IJ_AB' R^AB'_JM - R^JI_AB' R^AB'_JM) (I, J, M, A, B' in beta orbitals)
  //        + [(C0+C1)/2]^2 * R^JI_AB' R^AB_JM           (I, M, B' in beta orbitals, J, A in alpha orbitals)
  //        + (C0+C1)/2*(C0-C1)/2 * (R^IJ_AB' R^AB'_JM + R^JI_AB' R^AB'_MJ)
  //        + [(C0-C1)/2]^2 * R^IJ_AB' R^AB'_MJ
  //        + [(C0+C1)/2]^2 * R^JI_B'A R^B'A_JM           (I, M, A in beta orbitals, J, B' in alpha orbitals)
  //        + (C0+C1)/2*(C0-C1)/2 * (R^IJ_B'A R^B'A_JM + R^JI_B'A R^B'A_MJ)
  //        + [(C0-C1)/2]^2 * R^IJ_B'A R^B'A_MJ

  compute_Dmi_2(nspincases1, nspincases2, C_0, C_1,
                v_orbs1_ab, v_orbs2_ab,
                Dm_i_alpha, Dm_i_beta);

  if (debug_ >= DefaultPrintThresholds::N2) {
    print_intermediate("Alpha", "D^m_i", Dm_i_alpha, nocc1_act, nocc1_act);
    if (nspincases1 == 2)
      print_intermediate("Beta", "D^m_i", Dm_i_beta, nocc2_act, nocc2_act);
  }

  // D^c_b:
  // Alpha d^C_B = C1^2 * (R^A'B_IJ R^IJ_A'C - R^BA'_IJ R^IJ_A'C) (I, J, A', B, C in alpha orbitals)
  //             + [(C0+C1)/2]^2 * R^BA'_IJ R^IJ_CA'  (I, B, C in alpha orbitals, J, A' in beta orbitals)
  //             + (C0+C1)/2*(C0-C1)/2 * (R^A'B_IJ R^IJ_CA' + R^BA'_IJ R^IJ_A'C)
  //             + [(C0-C1)/2]^2 * R^A'B_IJ R^IJ_A'C
  //
  // Beta  d^C_B = C1^2 * (R^A'B_IJ R^IJ_A'C - R^BA'_IJ R^IJ_A'C) (I, J, A', B, C in beta orbitals)
  //             + [(C0+C1)/2]^2 * R^A'B_IJ R^IJ_A'C  (I, A' in alpha orbitals, J, B, C in beta orbitals)
  //             + (C0+C1)/2*(C0-C1)/2 * (R^BA'_IJ R^IJ_A'C + R^A'B_IJ R^IJ_CA')
  //             + [(C0-C1)/2]^2 * R^BA'_IJ R^IJ_CA'
  compute_RR31_32_spin(vir_vir_cabs,
                       nspincases1, nspincases2, C_0, C_1,
                       v_orbs1_ab, v_orbs2_ab,
                       Dc_b_alpha, Dc_b_beta);

  if (debug_ >= DefaultPrintThresholds::N2) {
    print_intermediate("Alpha", "D^c_b", Dc_b_alpha, nvir1, nvir1);
    if (nspincases1 == 2)
      print_intermediate("Beta", "D^c_b", Dc_b_beta, nvir2, nvir2);
  }

 // test for D^c_b
#if 0
  ExEnv::out0() << endl << "testing D^c_b : " << endl;
  compute_Dcb(nspincases1, nspincases2, C_0, C_1);
#endif

  // D^c'_b':
  // Alpha d^C'_B' =
  // 1st sum over A:  C1^2 * (R^AB'_IJ R^IJ_AC' - R^B'A_IJ R^IJ_AC') (I, J, A, B', C' in alpha orbitals)
  //                + [(C0+C1)/2]^2 * R^B'A_IJ R^IJ_C'A  (I, B', C' in alpha orbitals, J, A in beta orbitals)
  //                + (C0+C1)/2*(C0-C1)/2 * (R^AB'_IJ R^IJ_C'A + R^B'A_IJ R^IJ_AC')
  //                + [(C0-C1)/2]^2 * R^AB'_IJ R^IJ_AC'
  //
  // 2nd sum over A': C1^2 * (R^A'B'_IJ R^IJ_A'C' - R^B'A'_IJ R^IJ_A'C') (I, J, A', B', C' in alpha orbitals)
  //                + [(C0+C1)/2]^2 * R^B'A'_IJ R^IJ_C'A'  (I, B', C' in alpha orbitals, J, A' in beta orbitals)
  //                + (C0+C1)/2*(C0-C1)/2 * (R^A'B'_IJ R^IJ_C'A' + R^B'A'_IJ R^IJ_A'C')
  //                + [(C0-C1)/2]^2 * R^A'B'_IJ R^IJ_A'C'
  //
  // Beta  d^C'_B' =
  // 1st sum over A:  C1^2 * (R^AB'_IJ R^IJ_AC' - R^B'A_IJ R^IJ_AC') (I, J, A, B', C' in beta orbitals)
  //                + [(C0+C1)/2]^2 * R^AB'_IJ R^IJ_AC'  (I, A in alpha orbitals, J, B', C' in beta orbitals)
  //                + (C0+C1)/2*(C0-C1)/2 * (R^B'A_IJ R^IJ_AC' + R^AB'_IJ R^IJ_C'A)
  //                + [(C0-C1)/2]^2 * R^B'A_IJ R^IJ_C'A

  // 2nd sum over A': C1^2 * (R^A'B'_IJ R^IJ_A'C' - R^B'A'_IJ R^IJ_A'C') (I, J, A', B, C in beta orbitals)
  //                + [(C0+C1)/2]^2 * R^A'B'_IJ R^IJ_A'C'  (I, A in alpha orbitals, J, B', C' in beta orbitals)
  //                + (C0+C1)/2*(C0-C1)/2 * (R^B'A'_IJ R^IJ_A'C' + R^A'B'_IJ R^IJ_C'A')
  //                + [(C0-C1)/2]^2 * R^B'A'_IJ R^IJ_C'A'
  compute_RR31_32_spin(cabs_cabs_vir,
                       nspincases1, nspincases2, C_0, C_1,
                       v_orbs1_ab, v_orbs2_ab,
                       Dcp_bp_alpha_A, Dcp_bp_beta_A);

  if (debug_ >= DefaultPrintThresholds::N2) {
    print_intermediate("Alpha", "D^c'_b' sum over a", Dcp_bp_alpha_A, ncabs1, ncabs1);
    if (nspincases1 == 2)
      print_intermediate("Beta", "D^c'_b' sum over a", Dcp_bp_beta_A, ncabs2, ncabs2);
  }

  // test for D^c'_b' sum over a
#if 0
    ExEnv::out0() << endl << "testing D^c'_b' sum over a : " << endl;
    compute_Dcpbp_a(nspincases1, nspincases2, C_0, C_1);
#endif

  compute_RR31_32_spin(cabs_cabs_cabs,
                       nspincases1, nspincases2, C_0, C_1,
                       v_orbs1_ab, v_orbs2_ab,
                       Dcp_bp_alpha_Ap, Dcp_bp_beta_Ap);

  if (debug_ >= DefaultPrintThresholds::N2) {
    print_intermediate("Alpha", "D^c'_b' sum over a'", Dcp_bp_alpha_Ap, ncabs1, ncabs1);
    if (nspincases1 == 2)
      print_intermediate("Beta", "D^c'_b' sum over a'", Dcp_bp_beta_Ap, ncabs2, ncabs2);
  }

  // test for D^c'_b' sum over a'
  #if 0
    ExEnv::out0() << endl << "testing D^c'_b' sum over a' : " << endl;
    compute_Dcpbp_ap(nspincases1, nspincases2, C_0, C_1);
  #endif

  // test code: trace of D^m_i = - trace of (D^b_c + D^b'_c')
#if 0
    // D^m_i computed through R^ij_a'b' R_ij^a'b' + 2 R^ij_ab' R_ij^ab'
    double Tr_Dmi = 0;
    for (int i = 0; i != nocc1_act; ++i) {
      const int idx = i * nocc1_act + i;
      Tr_Dmi += Dm_i_alpha[idx];
    }
    ExEnv::out0() << endl << "Trace of alpha D^m_i: " << scprintf("%12.10f", Tr_Dmi) << endl;

    Tr_Dmi = 0;
    for (int i = 0; i != nocc2_act; ++i) {
      const int idx = i * nocc2_act + i;
      Tr_Dmi += Dm_i_beta[idx];
    }
    ExEnv::out0() << endl << "Trace of beta D^m_i: " << scprintf("%12.10f", Tr_Dmi) << endl;


    double Tr_alpha = 0;
    double Tr_beta = 0;

    double Tr_Dcb = 0;
    for (int i = 0; i != nvir1; ++i) {
      const int idx = i * nvir1 + i;
      Tr_Dcb += Dc_b_alpha[idx];
    }
    ExEnv::out0() << endl << "Trace of alpha D^c_b: " << scprintf("%12.10f", Tr_Dcb) << endl;
    Tr_alpha += Tr_Dcb;

    Tr_Dcb = 0;
    for (int i = 0; i != nvir2; ++i) {
      const int idx = i * nvir2 + i;
      Tr_Dcb += Dc_b_beta[idx];
    }
    ExEnv::out0() << endl << "Trace of beta D^c_b: " << scprintf("%12.10f", Tr_Dcb) << endl;
    Tr_beta += Tr_Dcb;

    double Tr_Dcpbp = 0;
    for (int i = 0; i != ncabs1; ++i) {
      const int idx = i * ncabs1 + i;
      Tr_Dcpbp += Dcp_bp_alpha_A[idx];
    }
    ExEnv::out0() << endl << "Trace of alpha D^c'_b' sum over a part: " << scprintf("%12.10f", Tr_Dcpbp) << endl;
    Tr_alpha += Tr_Dcpbp;

    Tr_Dcpbp = 0;
    for (int i = 0; i != ncabs2; ++i) {
      const int idx = i * ncabs2 + i;
      Tr_Dcpbp += Dcp_bp_beta_A[idx];
    }
    ExEnv::out0() << endl << "Trace of beta D^c'_b' sum over a part: " << scprintf("%12.10f", Tr_Dcpbp) << endl;
    Tr_beta += Tr_Dcpbp;

    Tr_Dcpbp = 0;
    for (int i = 0; i != ncabs1; ++i) {
      const int idx = i * ncabs1 + i;
      Tr_Dcpbp += Dcp_bp_alpha_Ap[idx];
    }
    ExEnv::out0() << endl << "Trace of alpha D^c'_b' sum over a' part: " << scprintf("%12.10f", Tr_Dcpbp) << endl;
    Tr_alpha += Tr_Dcpbp;

    Tr_Dcpbp = 0;
    for (int i = 0; i != ncabs2; ++i) {
      const int idx = i * ncabs2 + i;
      Tr_Dcpbp += Dcp_bp_beta_Ap[idx];
    }
    ExEnv::out0() << endl << "Trace of beta D^c'_b' sum over a' part: " << scprintf("%12.10f", Tr_Dcpbp) << endl;
    Tr_beta += Tr_Dcpbp;

    // test: trace of D^m_i = - trace of (D^c_b + D^c'_b')
    ExEnv::out0() << endl << "Trace of alpha D^c_b plus D^c'_b': " << scprintf("%12.10f", Tr_alpha) << endl;
    ExEnv::out0() << endl << "Trace of beta D^c_b plus D^c'_b': " << scprintf("%12.10f", Tr_beta) << endl;
#endif
  //
  // D^a'_a
  // Alpha d^A'_A = 1/2 * C1 * (R^IJ_A'B - R^JI_A'B) * T^AB_IJ  (I, J in alpha orbitals)
  //               + [(C0+C1)/2 * R^IJ_A'B + (C0-C1)/2 * R^JI_A'B] * T^AB_IJ  (I in alpha orbital, J in beta orbital)
  //
  //               + C1^2 * (R^AB'_IJ R^IJ_A'B' - R^B'A_IJ R^IJ_A'B')  (I, J in alpha orbitals)
  //               + [(C0+C1)/2]^2 * R^AB'_IJ R^IJ_A'B'  (I, A', A in alpha orbitals, J, B' in beta orbitals)
  //               + (C0+C1)/2*(C0-C1)/2 * (R^B'A_IJ R^IJ_A'B' + R^AB'_IJ R^IJ_B'A')
  //               + [(C0-C1)/2]^2 * R^B'A_IJ R^IJ_B'A'
  //
  // Beta  d^A'_A = 1/2 * C1 * (R^IJ_BA'- R^JI_BA') * T^BA_IJ  (I, J in beta orbitals)
  //               + [(C0+C1)/2 * R^IJ_BA'+ (C0-C1)/2 * R^JI_BA'] * T^BA_IJ  (I in alpha orbital, J in beta orbital)
  //
  //               + C1^2 * (R^AB'_IJ R^IJ_A'B' - R^B'A_IJ R^IJ_A'B') (I, J in beta orbitals)
  //               + [(C0+C1)/2]^2 * R^B'A_IJ R^IJ_B'A'  (I, B' in alpha orbitals, J, A', A in beta orbitals)
  //               + (C0+C1)/2*(C0-C1)/2 * (R^AB'_IJ R^IJ_B'A' + R^B'A_IJ R^IJ_A'B')
  //               + [(C0-C1)/2]^2 * R^AB'_IJ R^IJ_A'B'
  compute_RR31_32_spin(cabs_vir_cabs,
                       nspincases1, nspincases2, C_0, C_1,
                       v_orbs1_ab, v_orbs2_ab,
                       Dap_a_alpha_RR, Dap_a_beta_RR);

  if (debug_ >= DefaultPrintThresholds::N2) {
    print_intermediate("Alpha", "D^a'_a (RR part)", Dap_a_alpha_RR, ncabs1, nvir1);
    if (nspincases1 == 2)
      print_intermediate("Beta", "D^a'_a (RR part)", Dap_a_beta_RR, ncabs2, nvir2);
  }
  // test for D^a'_a RR part
#if 0
//  print_intermediate("Alpha", "D^a'_a RR part", Dap_a_alpha_RR, ncabs1, nvir1);
//  print_intermediate("Beta", "D^a'_a RR part", Dap_a_beta_RR, ncabs2, nvir2);
  ExEnv::out0() << endl << "testing D^a'_a RR part: " << endl;
  compute_Dapa_RR(nspincases1, nspincases2, C_0, C_1);
#endif

  double* Dap_a_alpha_RT;
  double* Dap_a_beta_RT;
  if (this->r12eval()->ebc() == false
      || this->r12eval()->coupling() == true) {

    Dap_a_alpha_RT = new double[ncabs1_vir1];
    fill_n(Dap_a_alpha_RT, ncabs1_vir1, 0.0);
    if (nspincases1 == 2) {
      const int ncabs2_vir2 = ncabs2 * nvir2;
      Dap_a_beta_RT = new double[ncabs2_vir2];
      fill_n(Dap_a_beta_RT, ncabs2_vir2, 0.0);
    } else if (nspincases1 == 1) {
        Dap_a_beta_RT = Dap_a_alpha_RT;
    }

    if (r12intermediates_->T2_cc_computed()) {
      compute_RT2_apa(nspincases1, nspincases2, C_0, C_1,
                      v_orbs1_ab, v_orbs2_ab,
                      Dap_a_alpha_RT, Dap_a_beta_RT);

      if (debug_ >= DefaultPrintThresholds::N2) {
        print_intermediate("Alpha", "D^a'_a (RT2 part)", Dap_a_alpha_RT, ncabs1, nvir1);
        if (nspincases1 == 2)
          print_intermediate("Beta", "D^a'_a (RT2 part)", Dap_a_beta_RT, ncabs2, nvir2);
      }
    } else { // in case both ebc and coupling are specified
        compute_RTmp2_apa(nspincases1, nspincases2, C_0, C_1,
                          v_orbs1_ab,v_orbs2_ab,
                          Dap_a_alpha_RT, Dap_a_beta_RT);
    }
  }

  RefSCMatrix D[NSpinCases1];
  RefSCMatrix D_cc[NSpinCases1];
  RefSCDimension rowdim = new SCDimension(norb_com);
  RefSCDimension coldim = new SCDimension(norb_com);
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

  RefSCVector dipoles = localkit->vector(RefSCDimension(new SCDimension(3)));
  dipoles.assign(0.0);

  for (int s = 0; s < nspincases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);

    D[spin] = localkit->matrix(rowdim, coldim);
    D[spin].assign(0.0);

    const int nocc_act = (spin == Alpha? nocc1_act : nocc2_act);
    const double* iter_Dmi = (spin == Alpha? Dm_i_alpha : Dm_i_beta);
    for (int i = 0; i < nocc_act; ++i) {
        for (int j = 0; j < nocc_act; ++j, ++iter_Dmi){
          D[spin].set_element(i, j, *iter_Dmi);
        }
    }
    //const string spin_label = (spin == Alpha? "Alpha" : "Beta");
    //ExEnv::out0() << endl << spin_label << " trace of Dij: " << scprintf("%12.10f", D[spin].trace()) << endl;

    const double* iter_Dcb = (spin == Alpha? Dc_b_alpha : Dc_b_beta);
    for (int a = nocc_act; a < norb; ++a) {
        for (int b = nocc_act; b < norb; ++b, ++iter_Dcb){
          D[spin].set_element(a, b, *iter_Dcb);
        }
    }
    //D[spin].print(prepend_spincase(spin,"F12 one-particle density Dij+Dab:").c_str());
    //ExEnv::out0() << endl << spin_label << " trace of Dij + Dab: " << scprintf("%12.10f", D[spin].trace()) << endl;

    const double* iter_Dcpbp_a = (spin == Alpha? Dcp_bp_alpha_A : Dcp_bp_beta_A);
    const double* iter_Dcpbp_ap = (spin == Alpha? Dcp_bp_alpha_Ap : Dcp_bp_beta_Ap);
    for (int ap = norb; ap < norb_com; ++ap) {
        for (int bp = norb; bp < norb_com; ++bp, ++iter_Dcpbp_a, ++iter_Dcpbp_ap) {
          const double Dapbp =  *iter_Dcpbp_a + *iter_Dcpbp_ap;
          //const double Dapbp =  *iter_Dcpbp_ap ;
          D[spin].set_element(ap, bp, Dapbp);
        }
    }
    //D[spin].print(prepend_spincase(spin,"F12 one-particle density Dij+Dab+Da'b':").c_str());

    if (this->r12eval()->ebc() == false
        || this->r12eval()->coupling() == true) {

//      bool ebc = this->r12eval()->ebc() == false;
//      bool coupling = this->r12eval()->coupling() == true;
//      ExEnv::out0() << endl << ebc << "  " << coupling << endl;

      const double* iter_Dap_a_RR = (spin == Alpha? Dap_a_alpha_RR : Dap_a_beta_RR);
      const double* iter_Dap_a_RT = (spin == Alpha? Dap_a_alpha_RT : Dap_a_beta_RT);
      for (int ap = norb; ap < norb_com; ++ap) {
        for (int a = nocc_act; a < norb; ++a, ++iter_Dap_a_RR, ++iter_Dap_a_RT) {
            const double Dapa = *iter_Dap_a_RR + *iter_Dap_a_RT;
            D[spin].set_element(ap, a, Dapa);
          }
      }

    } else {
        const double* iter_Dap_a = (spin == Alpha? Dap_a_alpha_RR : Dap_a_beta_RR);
        for (int ap = norb; ap < norb_com; ++ap) {
          for (int a = nocc_act; a < norb; ++a, ++iter_Dap_a) {
              const double Dapa = *iter_Dap_a;
              D[spin].set_element(ap, a, Dapa);
            }
        }
    }
//    if (debug_ >= DefaultPrintThresholds::allN2)
      D[spin].print(prepend_spincase(spin,"F12 one-particle density:").c_str());


    //
    // Obtain CCSD one-particle density from Psi
    D_cc[spin] = r12intermediates_->get_1rdm_cc(spin);
    //D_cc[spin].print(prepend_spincase(spin,"CCSD one-particle density:").c_str());

#if 1
    // test code for D and D_cc:
    // trace of D = 0, as trace of Dij = - trace of (Dab + Da'b')
    const string spin_label = (spin == Alpha? "Alpha" : "Beta");
    ExEnv::out0() << endl << spin_label << " Trace of D_f12: " << scprintf("%12.10f", D[spin].trace())<< endl;
    // trace of D_cc = # of electrons
    ExEnv::out0() << endl << spin_label << " Trace of D_cc: " << scprintf("%12.10f", D_cc[spin].trace()) << endl;
#endif
    // add psi ccsd and r12 one-particle density together
    D[spin].accumulate_subblock(D_cc[spin], 0, norb-1, 0, norb-1, 0, 0);
    if (debug_ >= DefaultPrintThresholds::allN2)
      D[spin].print(prepend_spincase(spin,"CCSD_F12 one-particle density:").c_str());

    // transform D_cc to MPQC ordering
    // not RefSymmSCMatrix ?
    RefSCMatrix opdm = onepdm_transformed(spin, D[spin]);
   // opdm.print("CCSD_F12 one-particle density matrix MPQC ordering");

    // computing dipole integrals in MO basis
    RefSCMatrix MXX, MYY, MZZ, MXY, MXZ, MYZ;

//    const Ref<OrbitalSpace> space = r12world->ribs_space();
//    compute_multipole_ints(space, space,
//                           MX, MY, MZ,
//                           MXX, MYY, MZZ,
//                           MXY, MXZ, MYZ);
//
//    // ?? why need?
//    const RefSCDimension M_dim(new SCDimension(norb_com));
//    RefSCMatrix MX_nb = localkit->matrix(M_dim,M_dim);
//    RefSCMatrix MY_nb = localkit->matrix(M_dim,M_dim);
//    RefSCMatrix MZ_nb = localkit->matrix(M_dim,M_dim);
//    for(int i = 0; i < M_dim.n(); i++) {
//      for(int j = 0; j < M_dim.n(); j++) {
//        MX_nb.set_element(i,j,MX.get_element(i,j));
//        MY_nb.set_element(i,j,MY.get_element(i,j));
//        MZ_nb.set_element(i,j,MZ.get_element(i,j));
//      }
//    }
//    //MZ.print(prepend_spincase(spin,"mu(Z)_nb").c_str());

    RefSCMatrix MX_orbs, MY_orbs, MZ_orbs;
    const Ref<OrbitalSpace>& space_orbs = (spin == Alpha? orbs1 : orbs2);
    compute_multipole_ints(space_orbs, space_orbs,
                           MX_orbs, MY_orbs, MZ_orbs,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);
    //MZ_orbs.print(prepend_spincase(spin,"mu(Z)_nb in orbs").c_str());

    // experiment dipole ints of ribs ?= dipole ints of orbs + cabs
    RefSCMatrix MX_cabs, MY_cabs, MZ_cabs;
    const Ref<OrbitalSpace>& space_cabs = (spin == Alpha? cabs1 : cabs2);
    compute_multipole_ints(space_cabs, space_cabs,
                           MX_cabs, MY_cabs, MZ_cabs,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);
    //MZ_cabs.print(prepend_spincase(spin,"mu(Z)_nb in cabs").c_str());

    RefSCMatrix MX_cabsorbs, MY_cabsorbs, MZ_cabsorbs;
    compute_multipole_ints(space_cabs, space_orbs,
                           MX_cabsorbs, MY_cabsorbs, MZ_cabsorbs,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);
    //MZ_cabsorbs.print(prepend_spincase(spin,"mu(Z)_nb in cabs orbs").c_str());

    // clean
    MXX = 0;
    MYY = 0;
    MZZ = 0;
    MXY = 0;
    MXZ = 0;
    MYZ = 0;

    // obtain the dipole ints in orbs+cabs space
    const int nirreps = space_cabs->nblocks();
    const std::vector<unsigned int>& orbspi = space_orbs->block_sizes();
    const std::vector<unsigned int>& cabspi = space_cabs->block_sizes();
    std::vector<unsigned int> orbsoff(nirreps);
    std::vector<unsigned int> cabsoff(nirreps);
    std::vector<unsigned int> ribsoff(nirreps);

    orbsoff[0] = 0;
    cabsoff[0] = 0;
    ribsoff[0] = 0;

    for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
      orbsoff[irrep] = orbsoff[irrep-1] + orbspi[irrep-1];
      cabsoff[irrep] = cabsoff[irrep-1] + cabspi[irrep-1];
      ribsoff[irrep] = orbsoff[irrep] + cabsoff[irrep];
    }

    const RefSCDimension M_dim_ribs(new SCDimension(norb_com));
    RefSCMatrix MX_nb_ribs = localkit->matrix(M_dim_ribs,M_dim_ribs);
    RefSCMatrix MY_nb_ribs = localkit->matrix(M_dim_ribs,M_dim_ribs);
    RefSCMatrix MZ_nb_ribs = localkit->matrix(M_dim_ribs,M_dim_ribs);
    MZ_nb_ribs.assign(0.0);

    for (unsigned int h = 0; h < nirreps; ++h) {

      for (int p = 0; p < orbspi[h]; ++p) {
          for (int q = 0; q <= p; ++q) {
            const int idxp1 = p + orbsoff[h];
            const int idxq1 = q + orbsoff[h];
            const int idxp2 = p + ribsoff[h];
            const int idxq2 = q + ribsoff[h];

            MX_nb_ribs.set_element(idxp2,idxq2,MX_orbs.get_element(idxp1,idxq1));
            MY_nb_ribs.set_element(idxp2,idxq2,MY_orbs.get_element(idxp1,idxq1));
            MZ_nb_ribs.set_element(idxp2,idxq2,MZ_orbs.get_element(idxp1,idxq1));
          }
      }

      for (int ap = 0; ap < cabspi[h]; ++ap) {
        for (int bp = 0; bp <= ap; ++bp) {
          const int idx_ap1 = cabsoff[h] + ap;
          const int idx_bp1 = cabsoff[h] + bp;
          const int idx_ap2 = ribsoff[h] + orbspi[h] + ap;
          const int idx_bp2 = ribsoff[h] + orbspi[h] + bp;

          MX_nb_ribs.set_element(idx_ap2,idx_bp2,MX_cabs.get_element(idx_ap1,idx_bp1));
          MY_nb_ribs.set_element(idx_ap2,idx_bp2,MY_cabs.get_element(idx_ap1,idx_bp1));
          MZ_nb_ribs.set_element(idx_ap2,idx_bp2,MZ_cabs.get_element(idx_ap1,idx_bp1));
        }
      }

      for (int ap = 0; ap < cabspi[h]; ++ap) {
        for (int p = 0; p < orbspi[h]; ++p) {
          const int idx_ap1 = cabsoff[h] + ap;
          const int idx_p1 = orbsoff[h] + p;
          const int idx_ap2 = ribsoff[h] + orbspi[h] + ap;
          const int idx_p2 = ribsoff[h] + p;

          MX_nb_ribs.set_element(idx_ap2,idx_p2,MX_cabsorbs.get_element(idx_ap1,idx_p1));
          MY_nb_ribs.set_element(idx_ap2,idx_p2,MY_cabsorbs.get_element(idx_ap1,idx_p1));
          MZ_nb_ribs.set_element(idx_ap2,idx_p2,MZ_cabsorbs.get_element(idx_ap1,idx_p1));
        }
      }
    }
    // Symmetrize matrices
    for(int p = 0; p < norb_com; p++) {
      for(int q = 0; q < p; q++) {
        MX_nb_ribs.set_element(q,p,MX_nb_ribs.get_element(p,q));
        MY_nb_ribs.set_element(q,p,MY_nb_ribs.get_element(p,q));
        MZ_nb_ribs.set_element(q,p,MZ_nb_ribs.get_element(p,q));
      }
    }

    //MZ_nb_ribs.print(prepend_spincase(spin,"mu(Z)_nb in ribs").c_str());


    RefSCMatrix opdm_x_ribs = MX_nb_ribs * opdm;
    RefSCMatrix opdm_y_ribs = MY_nb_ribs * opdm;
    RefSCMatrix opdm_z_ribs = MZ_nb_ribs * opdm;

    const double dx_ribs = opdm_x_ribs.trace();
    const double dy_ribs = opdm_y_ribs.trace();
    const double dz_ribs = opdm_z_ribs.trace();
    ExEnv::out0() << endl << spin_label <<  " CCSD_F12 correlation dipole moment x y z: " << scprintf("%12.10f",dx_ribs)
                          << " " << scprintf("%12.10f",dy_ribs) << " " << scprintf("%12.10f",dz_ribs) << endl;

    dipoles[0] = dipoles[0] + dx_ribs;
    dipoles[1] = dipoles[1] + dy_ribs;
    dipoles[2] = dipoles[2] + dz_ribs;

    // clean
    MX_cabs = 0;
    MY_cabs = 0;
    MZ_cabs = 0;
    MX_cabsorbs = 0;
    MY_cabsorbs = 0;
    MZ_cabsorbs = 0;
    MX_nb_ribs = 0;
    MY_nb_ribs = 0;
    MZ_nb_ribs = 0;

    opdm_x_ribs = 0;
    opdm_y_ribs = 0;
    opdm_z_ribs = 0;

    // test code: compute the ccsd dipole moment
#if 1
    {
      const RefSCDimension M_orbs_dim(new SCDimension(norb));
      RefSCMatrix MX2_nb = localkit->matrix(M_orbs_dim,M_orbs_dim);
      RefSCMatrix MY2_nb = localkit->matrix(M_orbs_dim,M_orbs_dim);
      RefSCMatrix MZ2_nb = localkit->matrix(M_orbs_dim,M_orbs_dim);
      for(int i = 0; i < M_orbs_dim.n(); i++) {
        for(int j = 0; j < M_orbs_dim.n(); j++) {
          MX2_nb.set_element(i,j,MX_orbs.get_element(i,j));
          MY2_nb.set_element(i,j,MY_orbs.get_element(i,j));
          MZ2_nb.set_element(i,j,MZ_orbs.get_element(i,j));
        }
      }

      // copy psi ccsd density in MPQC ordering back
      RefSCMatrix opdm_cc = onepdm_transformed2(spin, D_cc[spin]);
      //opdm_cc.print("CCSD one-particle density matrix in MPQC ordering");


      RefSCMatrix opdm_x_cc = MX2_nb * opdm_cc;
      RefSCMatrix opdm_y_cc = MY2_nb * opdm_cc;
      RefSCMatrix opdm_z_cc = MZ2_nb * opdm_cc;

      const double dx_cc = opdm_x_cc.trace();
      const double dy_cc = opdm_y_cc.trace();
      const double dz_cc = opdm_z_cc.trace();
      ExEnv::out0() << endl << spin_label <<  " ccsd correlation dipole x y z: " << scprintf("%12.10f",dx_cc)
                            << " " << scprintf("%12.10f",dy_cc) << " " << scprintf("%12.10f",dz_cc) << endl;
      // clean
      opdm_cc = 0;
      MX2_nb = 0;
      MY2_nb = 0;
      MZ2_nb = 0;
      opdm_x_cc = 0;
      opdm_y_cc = 0;
      opdm_z_cc = 0;
    }
#endif
    MX_orbs = 0;
    MY_orbs = 0;
    MZ_orbs = 0;

  } // end of loop over nspincase1

  Ref<Molecule> molecule = orbs1->basis()->molecule();
  int natom = molecule->natom();
  // test
  molecule->print();
  double dx_nuc = 0.0, dy_nuc = 0.0, dz_nuc = 0.0;
  for(int i = 0; i < natom; i++) {
    dx_nuc += molecule->r(i,0) * molecule->Z(i);
    dy_nuc += molecule->r(i,1) * molecule->Z(i);
    dz_nuc += molecule->r(i,2) * molecule->Z(i);
  }
  ExEnv::out0() << endl <<  "nuclear dipole x y z: " << scprintf("%12.10f",dx_nuc)
                        << " " << scprintf("%12.10f",dy_nuc)
                        << " " << scprintf("%12.10f",dz_nuc) << endl;

  delete[] Dm_i_alpha;
  delete[] Dc_b_alpha;
  delete[] Dcp_bp_alpha_A;
  delete[] Dcp_bp_alpha_Ap;
  delete[] Dap_a_alpha_RR;

  if (this->r12eval()->ebc() == false
      || this->r12eval()->coupling() == true) {
    delete[] Dap_a_alpha_RT;
    if (nspincases1 == 2)
      delete[] Dap_a_beta_RT;
  }

  if (nspincases1 == 2) {
       delete[] Dm_i_beta;
       delete[] Dc_b_beta;
       delete[] Dcp_bp_beta_A;
       delete[] Dcp_bp_beta_Ap;
       delete[] Dap_a_beta_RR;
   }

  if (nspincases1 == 1) {
    D[Beta] = D[Alpha];

    dipoles[0] = - dipoles[0] * 2 + dx_nuc;
    dipoles[1] = - dipoles[1] * 2 + dy_nuc;
    dipoles[2] = - dipoles[2] * 2 + dz_nuc;
  } else {
      dipoles[0] = - dipoles[0] + dx_nuc;
      dipoles[1] = - dipoles[1] + dy_nuc;
      dipoles[2] = - dipoles[2] + dz_nuc;
  }
  //ExEnv::out0() << endl <<  "dipole x y z: " << dipoles[0] << " " << dipoles[1] << " " << dipoles[2] << endl;
  dipoles.print("dipole moment: mu(X) mu(Y) mu(Z)");


//  D_ccsdf12_[Alpha] = D[Alpha];
//  D_ccsdf12_[Beta] = D[Beta];

  // test code
#if 1
  // the trace of D should be the same as D_cc (# of electrons)
  ExEnv::out0() << endl << "Trace of alpha CCSD_F12 one-particle density: " << scprintf("%12.10f", D[Alpha].trace()) << endl;
  ExEnv::out0() << endl << "Trace of beta CCSD_F12 one-particle density: " << scprintf("%12.10f", D[Beta].trace()) << endl;
#endif

  ExEnv::out0() << endl << "End of the computation for CCSD_F12 one-particle density" << endl;
  ExEnv::out0() << endl << "**********************************************" << endl;
  return;
}
// end of compute_density_diag


//void MP2R12Energy_Diag::compute_dipole_ints(const SpinCase1 spin,
//                                               RefSCMatrix& MX, RefSCMatrix& MY,
//                                               RefSCMatrix& MZ) {
////  Ref<OrbitalSpace> space = r12eval_->orbs(spin);
//  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
//  const Ref<OrbitalSpace> space = r12world->ribs_space();
//
//  const Ref<GaussianBasisSet> bs = space->basis();
//  Ref<Integral> localints = space->integral()->clone();
//  localints->set_basis(bs,bs);
//
//  Ref<OneBodyInt> m_ints = localints->dipole(0);
//
//  // form AO moment matrices
//  RefSCMatrix vec_t = space->coefs().t();
//  RefSCMatrix vec = space->coefs();
//  RefSCDimension aodim1 = vec_t.coldim();
//  RefSCDimension aodim2 = vec.rowdim();
//  Ref<SCMatrixKit> aokit = bs->so_matrixkit();
//
//  RefSCMatrix mx_ao(aodim1, aodim2, aokit);
//  RefSCMatrix my_ao(aodim1, aodim2, aokit);
//  RefSCMatrix mz_ao(aodim1, aodim2, aokit);
//
//  mx_ao.assign(0.0);
//  my_ao.assign(0.0);
//  mz_ao.assign(0.0);
//
//  const int nshell = bs->nshell();
//  ExEnv::out0() << endl << "number of shells: " << nshell << endl;
//  for(int sh1 = 0; sh1 < nshell; sh1++) {
//    int bf1_offset = bs->shell_to_function(sh1);
//    int nbf1 = bs->shell(sh1).nfunction();
//
//    for(int sh2 = 0; sh2 <= sh1; sh2++) {
//      int bf2_offset = bs->shell_to_function(sh2);
//      int nbf2 = bs->shell(sh2).nfunction();
//
//      m_ints->compute_shell(sh1,sh2);
//      const double *m_ints_ptr = m_ints->buffer();
//
//      int bf1_index = bf1_offset;
//      for(int bf1 = 0; bf1 < nbf1; bf1++, bf1_index++, m_ints_ptr+=3*nbf2) {
//        int bf2_index = bf2_offset;
//
//        const double *ptr = m_ints_ptr;
//        for(int bf2 = 0; bf2 < nbf2; bf2++, bf2_index++) {
//
//          mx_ao.set_element(bf1_index, bf2_index, *(ptr++));
//          my_ao.set_element(bf1_index, bf2_index, *(ptr++));
//          mz_ao.set_element(bf1_index, bf2_index, *(ptr++));
//
//        }
//      }
//    }
//    // end of loop over sh2
//  }
//  // end of loop over sh1
//
//  // and clean up a bit
//  m_ints = 0;
//
//  // Symmetrize matrices
//  const int nbasis = bs->nbasis();
//  for(int bf1 = 0; bf1 < nbasis; bf1++) {
//    for(int bf2 = 0; bf2 <= bf1; bf2++) {
//      mx_ao(bf2,bf1) = mx_ao(bf1,bf2);
//      my_ao(bf2,bf1) = my_ao(bf1,bf2);
//      mz_ao(bf2,bf1) = mz_ao(bf1,bf2);
//    }
//  }
//
//  //mz_ao.print("mu(Z) in AO basis");
//
//  // finally, transform
////  RefSCMatrix mx_mo, my_mo, mz_mo;
//  MX = vec_t * mx_ao * vec;
//  MY = vec_t * my_ao * vec;
//  MZ = vec_t * mz_ao * vec;
//
//#if 0
//  // transform MZ from CC to QT ordering
//  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
//
//  const Ref<OrbitalSpace>& occ = r12eval()->occ_act(spin);
//  const Ref<OrbitalSpace>& vir = r12eval()->vir(spin);
//  const Ref<OrbitalSpace>& orbs = r12eval()->orbs(spin);
//  const Ref<OrbitalSpace>& cabs = r12world->cabs_space(spin);
//
//  const std::vector<unsigned int>& occpi = occ->block_sizes();
//  const std::vector<unsigned int>& uoccpi = vir->block_sizes();
//  const std::vector<unsigned int>& orbspi = orbs->block_sizes();
//  const std::vector<unsigned int>& cabspi = cabs->block_sizes();
//  const int nirreps = occ->nblocks();
//
//  const unsigned int nocc = occ->rank();
//  const unsigned int norbs = orbs->rank();
//  ExEnv::out0() << endl << "In dipole function: nirreps = " << nirreps << " norbs = " << norbs << endl;
//
//  std::vector<unsigned int> occioff(nirreps);
//  std::vector<unsigned int> uoccioff(nirreps);
//  std::vector<unsigned int> orbsoff(nirreps);
// // std::vector<unsigned int> cabsoff(nirreps);
//
//  occioff[0] = 0;
//  uoccioff[0] = 0;
//  orbsoff[0] = 0;
// // cabsoff[0] = 0;
//
//  for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
//    occioff[irrep] = occioff[irrep-1] + occpi[irrep-1];
//    uoccioff[irrep] = uoccioff[irrep-1] + uoccpi[irrep-1];
//    orbsoff[irrep] = orbsoff[irrep-1] + orbspi[irrep-1];
//    //cabsoff[irrep] = cabsoff[irrep-1] + cabsoff[irrep-1];
////    ExEnv::out0() << endl << "irrep: " << irrep << "  occpi: " << occpi[irrep]
////                                                << "  uoccpi: " << uoccpi[irrep]
////                                                << "  orbspi: " << orbspi[irrep] << endl;
//  }
//
//  const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
//  const RefSCDimension dim(new SCDimension(norbs));
//  RefSCMatrix MZ_mpqc = localkit->matrix(dim, dim);
//  MZ_mpqc.assign(0.0);
//  MZ_mpqc = vec_t * mz_ao * vec;
//  //MZ_mpqc.print(prepend_spincase(spin,"mpqc ordering mu(Z)").c_str());
//  MZ = localkit->matrix(dim, dim);
//  MZ.assign(0.0);
//
//  for (unsigned int h = 0; h < nirreps; ++h) {
//    const unsigned int irrep_offset = orbsoff[h];
//    const unsigned int i_offset = occioff[h];
//    const unsigned int j_offset = i_offset;
//
//    for (int i = 0; i < occpi[h]; ++i) {
//      for (int j = 0; j <= i; ++j) {
//        //ExEnv::out0() << endl << h << " " << i+irrep_offset << " " << j+irrep_offset << endl;
//        const double MZ_ij = MZ_mpqc.get_element(i+irrep_offset,j+irrep_offset);
//        MZ.set_element(i+i_offset, j+j_offset, MZ_ij);
//        //ExEnv::out0() << " MZ_ij = " << MZ_ij << endl;
//      }
//    }
//  }
//  //opdm_mat.print(prepend_spincase(spin,"MPQC ccsd opdm ij").c_str());
//
//  int irrep_offset = 0;
//  int nmopi = 0;
//  for (unsigned int h = 0; h < nirreps; ++h) {
//    //const unsigned int irrep_offset = orbsoff[h]+occioff[h+1];
//    irrep_offset = nmopi + occpi[h];
//    const unsigned int a_offset = nocc + uoccioff[h];
//    const unsigned int b_offset = a_offset;
//
//    for (int a = 0; a < uoccpi[h]; ++a) {
//      for (int b = 0; b <= a; ++b) {
//        //ExEnv::out0() << endl << h << " " << a+a_offset << " " << b+b_offset  << endl;
//        const double MZ_ab = MZ_mpqc.get_element(a+irrep_offset,b+irrep_offset);
//        MZ.set_element(a+a_offset, b+b_offset, MZ_ab);
//        //ExEnv::out0() << " irrep_offset " << irrep_offset<< " MZ_ab = " << MZ_ab << endl;
//      }
//    }
//    nmopi += orbspi[h];
//    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsd opdm ab").c_str());
//  }
//
//  for (unsigned int h = 0; h < nirreps; ++h) {
//    const unsigned int irrep_offset_a = orbsoff[h]+occpi[h];
//    const unsigned int irrep_offset_i = orbsoff[h];
//    const unsigned int a_offset = nocc + uoccioff[h];
//    const unsigned int i_offset = occioff[h];
//
//    for (int a = 0; a < uoccpi[h]; ++a) {
//      for (int i = 0; i < occpi[h]; ++i) {
//        const double MZ_ai = MZ_mpqc.get_element(a+irrep_offset_a,i+irrep_offset_i);
//        MZ.set_element(a+a_offset, i+i_offset, MZ_ai);
//        //ExEnv::out0() << endl << h << " " << a+mpqc_ab_offset << " " << i+mpqc_ij_offset << " Dai = " << Dai << endl;
//      }
//    }
//    //opdm_mat.print(prepend_spincase(spin,"h ordm transfromed").c_str());
//
//  }
//
//
//  // Symmetrize matrices
//  for(int p = 0; p < norbs; p++) {
//    for(int q = 0; q < p; q++) {
//      MZ.set_element(q,p,MZ_mpqc.get_element(p,q));
//      //ExEnv::out0() << endl << "p q Dpq: "<< p << " " << q << " " << opdm_mat.get_element(q,p) << endl;
//    }
//  }
//
//  MZ_mpqc = 0;
//#endif
//
//  // and clean up a bit
//  mx_ao = 0;
//  my_ao = 0;
//  mz_ao = 0;
//
////  if (debug_ > 1) {
////    MX.print(prepend_spincase(spin,"mu(X)").c_str());
////    MY.print(prepend_spincase(spin,"mu(Y)").c_str());
////    MZ.print(prepend_spincase(spin,"mu(Z)").c_str());
////  }
//}
//// end of compute_dipole_ints

// Transform the density matrix in ribs space (ordered in occ, uocc, cabs)
//                        to the MPQC ordering (ordered in irreps)

RefSCMatrix MP2R12Energy_Diag::onepdm_transformed(const SpinCase1& spin,
                                                  const RefSCMatrix& D) {

  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();

  const Ref<OrbitalSpace>& occ = r12eval()->occ(spin);
  const Ref<OrbitalSpace>& vir = r12eval()->vir(spin);
  const Ref<OrbitalSpace>& orbs = r12eval()->orbs(spin);
  const Ref<OrbitalSpace>& cabs = r12world->cabs_space(spin);
  const Ref<OrbitalSpace> ribs = r12world->ribs_space();

  const std::vector<unsigned int>& occpi = occ->block_sizes();
  const std::vector<unsigned int>& uoccpi = vir->block_sizes();
  const std::vector<unsigned int>& orbspi = orbs->block_sizes();
  const std::vector<unsigned int>& cabspi = cabs->block_sizes();
  const std::vector<unsigned int>& ribspi = ribs->block_sizes();

  const int nirreps = occ->nblocks();
  const int nocc = occ->rank();
  const int norbs = orbs->rank();
  const int ncabs = cabs->rank();
  const int nribs = ribs->rank();
//  ExEnv::out0() << endl << "nocc: " << nocc << endl
//                        << "norbs: " << norbs << endl
//                        << "ncabs: " << ncabs << endl;

  std::vector<unsigned int> occoff(nirreps);
  std::vector<unsigned int> uoccoff(nirreps);
  //std::vector<unsigned int> orbsoff(nirreps);
  std::vector<unsigned int> cabsoff(nirreps);
  std::vector<unsigned int> ribsoff(nirreps);

  occoff[0] = 0;
  uoccoff[0] = 0;
  //orbsoff[0] = 0;
  cabsoff[0] = 0;
  ribsoff[0] = 0;

//  ExEnv::out0() << endl << "  occpi uoccpi cabspi  ribspi:" << endl
//                        << 0 << "  " << occpi[0] << "  " << uoccpi[0] << "  " << cabspi[0] << "   "<< ribspi[0] <<  endl;
  for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
    occoff[irrep] = occoff[irrep-1] + occpi[irrep-1];
    uoccoff[irrep] = uoccoff[irrep-1] + uoccpi[irrep-1];
    //orbsoff[irrep] = orbsoff[irrep-1] + orbspi[irrep-1];
    cabsoff[irrep] = cabsoff[irrep-1] + cabspi[irrep-1];
    ribsoff[irrep] = ribsoff[irrep-1] + ribspi[irrep-1];
//    ExEnv::out0() << irrep << "   " << occpi[irrep] << "  " << uoccpi[irrep] << "  " << cabspi[irrep]
//                                                                  << "   " << ribspi[irrep]  << endl;
  }

  const RefSCDimension opdm_dim(new SCDimension(nribs));
  const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix opdm_mat = localkit->matrix(opdm_dim, opdm_dim);
  opdm_mat.assign(0.0);

  for (unsigned int h = 0; h < nirreps; ++h) {
    // offset for D matrix in occ, uocc, cabs ordering
    const unsigned int i_offset = occoff[h];
    const unsigned int j_offset = i_offset;
    const unsigned int a_offset = nocc + uoccoff[h];
    const unsigned int b_offset = a_offset;
    const unsigned int ap_offset = norbs + cabsoff[h];
    const unsigned int bp_offset = ap_offset;

    const int mpqc_ij_offset = ribsoff[h];
    const int mpqc_ab_offset = ribsoff[h] + occpi[h];
    const int mpqc_apbp_offset = ribsoff[h] + occpi[h] + uoccpi[h];

    for (int i = 0; i < occpi[h]; ++i) {
      for (int j = 0; j <= i; ++j) {
        const double Dij = D.get_element(i+i_offset,j+j_offset);
        opdm_mat.set_element(i+mpqc_ij_offset, j+mpqc_ij_offset, Dij);
      }
    }
    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsdf12 opdm ij").c_str());


    for (int a = 0; a < uoccpi[h]; ++a) {
      for (int b = 0; b <= a; ++b) {
        const double Dab = D.get_element(a+a_offset,b+b_offset);
        opdm_mat.set_element(a+mpqc_ab_offset, b+mpqc_ab_offset, Dab);
      }
    }
    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsdf12 opdm ab").c_str());

    for (int a = 0; a < uoccpi[h]; ++a) {
      for (int i = 0; i < occpi[h]; ++i) {
        const double Dai = D.get_element(a+a_offset,i+i_offset);
        opdm_mat.set_element(a+mpqc_ab_offset, i+mpqc_ij_offset, Dai);
      }
    }
    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsdf12 opdm ai").c_str());

    for (int ap = 0; ap < cabspi[h]; ++ap) {
      for (int bp = 0; bp <= ap; ++bp) {
        const double Dapbp = D.get_element(ap+ap_offset,bp+bp_offset);
        //ExEnv::out0() << endl << h << " " << ap+ap_offset << " " << bp+bp_offset << " D Da'b' = " << Dapbp << endl;
        opdm_mat.set_element(ap+mpqc_apbp_offset, bp+mpqc_apbp_offset, Dapbp);
        //ExEnv::out0() << endl << h << " " << ap+mpqc_apbp_offset << " " << bp+mpqc_apbp_offset << " Da'b' = " << Dapbp << endl;
      }
    }
    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsd opdm a'b'").c_str());

    for (int ap = 0; ap < cabspi[h]; ++ap) {
      for (int a = 0; a < uoccpi[h]; ++a) {
        const double Dapa = D.get_element(ap+ap_offset,a+a_offset);
        opdm_mat.set_element(ap+mpqc_apbp_offset, a+mpqc_ab_offset, Dapa);
        //ExEnv::out0() << endl << h << " " << ap+mpqc_apbp_offset << " " << a+mpqc_ab_offset << " Da'a = " << Dapa << endl;
      }
    }
    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsd opdm a'a").c_str());

  }

  // Symmetrize matrices
  for(int p = 0; p < nribs; p++) {
    for(int q = 0; q < p; q++) {
      opdm_mat.set_element(q,p,opdm_mat.get_element(p,q));
    }
  }


  if(debug_>=DefaultPrintThresholds::mostN2) {
    opdm_mat.print(prepend_spincase(spin,"one-particle density matrix in MPQC ordering").c_str());
  }

  return(opdm_mat);
}
// end of onepdm_transformed


RefSCMatrix MP2R12Energy_Diag::onepdm_transformed2(const SpinCase1& spin,
                                                   const RefSCMatrix& D) {

  const Ref<OrbitalSpace>& occ = r12eval()->occ(spin);
  const Ref<OrbitalSpace>& vir = r12eval()->vir(spin);
  const Ref<OrbitalSpace>& orbs = r12eval()->orbs(spin);

  const std::vector<unsigned int>& occpi = occ->block_sizes();
  const std::vector<unsigned int>& uoccpi = vir->block_sizes();
  const std::vector<unsigned int>& orbspi = orbs->block_sizes();

  const int nirreps = occ->nblocks();
  const int nocc = occ->rank();
  const int norbs = orbs->rank();

  std::vector<unsigned int> occoff(nirreps);
  std::vector<unsigned int> uoccoff(nirreps);
  std::vector<unsigned int> orbsoff(nirreps);

  occoff[0] = 0;
  uoccoff[0] = 0;
  orbsoff[0] = 0;

  for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
    occoff[irrep] = occoff[irrep-1] + occpi[irrep-1];
    uoccoff[irrep] = uoccoff[irrep-1] + uoccpi[irrep-1];
    orbsoff[irrep] = orbsoff[irrep-1] + orbspi[irrep-1];
  }

  const RefSCDimension opdm_dim(new SCDimension(norbs));
  const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix opdm_mat = localkit->matrix(opdm_dim, opdm_dim);
  opdm_mat.assign(0.0);

  for (unsigned int h = 0; h < nirreps; ++h) {
    // offset for D matrix in occ, uocc, cabs ordering
    const unsigned int i_offset = occoff[h];
    const unsigned int j_offset = i_offset;
    const unsigned int a_offset = nocc + uoccoff[h];
    const unsigned int b_offset = a_offset;

    const int mpqc_ij_offset = orbsoff[h];
    const int mpqc_ab_offset = orbsoff[h] + occpi[h];;

    for (int i = 0; i < occpi[h]; ++i) {
      for (int j = 0; j <= i; ++j) {
        const double Dij = D.get_element(i+i_offset,j+j_offset);
        opdm_mat.set_element(i+mpqc_ij_offset, j+mpqc_ij_offset, Dij);
      }
    }
    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsdf12 opdm ij").c_str());


    for (int a = 0; a < uoccpi[h]; ++a) {
      for (int b = 0; b <= a; ++b) {
        const double Dab = D.get_element(a+a_offset,b+b_offset);
        opdm_mat.set_element(a+mpqc_ab_offset, b+mpqc_ab_offset, Dab);
      }
    }
    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsdf12 opdm ab").c_str());

    for (int a = 0; a < uoccpi[h]; ++a) {
      for (int i = 0; i < occpi[h]; ++i) {
        const double Dai = D.get_element(a+a_offset,i+i_offset);
        opdm_mat.set_element(a+mpqc_ab_offset, i+mpqc_ij_offset, Dai);
      }
    }
    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsdf12 opdm ai").c_str());

  }

  // Symmetrize matrices
  for(int p = 0; p < norbs; p++) {
    for(int q = 0; q < p; q++) {
      opdm_mat.set_element(q,p,opdm_mat.get_element(p,q));
    }
  }

  return(opdm_mat);
}
