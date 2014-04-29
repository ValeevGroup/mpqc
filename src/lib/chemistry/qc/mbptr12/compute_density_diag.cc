/*
 * compute_density_diag.cc
 *
 *  Created on: Nov 8, 2011
 *      Author: jinmei
 */

#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>
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

    //ExEnv::out0() << "C_0: "  << C_0 << " C_1: " << C_1 << endl;
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
//      ExEnv::out0() << spinletters << " part d^" << i << "_" << i << " = "
//                    << scprintf("%12.10f", dii_12) << endl;

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
//        ExEnv::out0() << "AlphaBeta part: d^"  << i << "_" << i << " = "
//                      << scprintf("%12.10f", dii_12) << endl;
      } // end of summing AlphaBeta X

      ExEnv::out0() << indent << spinletters << " d^"  << i << "_" << i << " = "
                    << scprintf("%12.10f", dii_12) << endl;
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

  // transform CC T2(i,j,a,b) to T2(a,b,i,j)
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

  // form MX, MY, and MZ in ribs from occ, vir, and cabs space
  void form_DipoleInts_inRibs(const SpinCase1 spin,
                              const vector< Ref<OrbitalSpace> >& v_orbs1_ab,
                              const vector< Ref<OrbitalSpace> >& v_orbs2_ab,
                              RefSCMatrix& MX_nb_ribs, RefSCMatrix& MY_nb_ribs,
                              RefSCMatrix& MZ_nb_ribs,
                              RefSCMatrix& MX_ribs_nfc, RefSCMatrix& MY_ribs_nfc,
                              RefSCMatrix& MZ_ribs_nfc)
  {
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
    const int nocc1 = occ1->rank();
    const int nocc2 = occ2->rank();
    const int nvir1 = vir1->rank();
    const int nvir2 = vir2->rank();
    const int ncabs1 = cabs1->rank();
    const int ncabs2 = cabs2->rank();

    const int nribs_fc = nocc1_act + nvir1 + ncabs1;
    const int nribs_nfc = nocc1 + nvir1 + ncabs1;;

    // dipole integrals in ribs: dipole integrals in occ_act, vir, and CABS
    RefSCMatrix MXX, MYY, MZZ, MXY, MXZ, MYZ;
    RefSCMatrix MX_occ_act, MY_occ_act, MZ_occ_act;
    const Ref<OrbitalSpace>& space_occ_act = (spin == Alpha? occ1_act : occ2_act);
    compute_multipole_ints(space_occ_act, space_occ_act,
                           MX_occ_act, MY_occ_act, MZ_occ_act,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);

    RefSCMatrix MX_vir, MY_vir, MZ_vir;
    const Ref<OrbitalSpace>& space_vir = (spin == Alpha? vir1 : vir2);
    compute_multipole_ints(space_vir, space_vir,
                           MX_vir, MY_vir, MZ_vir,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);

    RefSCMatrix MX_virocc_act, MY_virocc_act, MZ_virocc_act;
    compute_multipole_ints(space_vir, space_occ_act,
                           MX_virocc_act, MY_virocc_act, MZ_virocc_act,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);

    RefSCMatrix MX_cabs, MY_cabs, MZ_cabs;
    //RefSCMatrix MXX_cabs, MYY_cabs, MZZ_cabs, MXY_cabs, MXZ_cabs, MYZ_cabs;
    const Ref<OrbitalSpace>& space_cabs = (spin == Alpha? cabs1 : cabs2);
    compute_multipole_ints(space_cabs, space_cabs,
                           MX_cabs, MY_cabs, MZ_cabs,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);

    RefSCMatrix MX_cabsocc_act, MY_cabsocc_act, MZ_cabsocc_act;
    compute_multipole_ints(space_cabs, space_occ_act,
                           MX_cabsocc_act, MY_cabsocc_act, MZ_cabsocc_act,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);

    RefSCMatrix MX_cabsvir, MY_cabsvir, MZ_cabsvir;
    compute_multipole_ints(space_cabs, space_vir,
                           MX_cabsvir, MY_cabsvir, MZ_cabsvir,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);

    const int nirreps = space_cabs->nblocks();
    const std::vector<unsigned int>& occ_actpi = space_occ_act->block_sizes();
    const std::vector<unsigned int>& virpi = space_vir->block_sizes();
    const std::vector<unsigned int>& cabspi = space_cabs->block_sizes();

    std::vector<unsigned int> occ_actoff(nirreps);
    std::vector<unsigned int> viroff(nirreps);
    std::vector<unsigned int> cabsoff(nirreps);
    // orbs = occ_act + vir
    std::vector<unsigned int> orbsoff(nirreps);
    // orbs = occ_act + vir + cabs
    std::vector<unsigned int> ribsoff(nirreps);

    occ_actoff[0] = 0;
    viroff[0] = 0;
    cabsoff[0] = 0;
    orbsoff[0] = 0;
    ribsoff[0] = 0;

    std::vector<unsigned int> orbspi(nirreps);
    orbspi[0] = occ_actpi[0] + virpi[0];

    for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
      occ_actoff[irrep] = occ_actoff[irrep-1] + occ_actpi[irrep-1];
      viroff[irrep] = viroff[irrep-1] + virpi[irrep-1];
      cabsoff[irrep] = cabsoff[irrep-1] + cabspi[irrep-1];

      orbsoff[irrep] = occ_actoff[irrep] + viroff[irrep];
      ribsoff[irrep] = orbsoff[irrep] + cabsoff[irrep];

      orbspi[irrep] = occ_actpi[irrep] + virpi[irrep];
    }

    // obtain the dipole ints in ribs
    for (unsigned int h = 0; h < nirreps; ++h) {

      for (int i = 0; i < occ_actpi[h]; ++i) {
        for (int j = 0; j <= i; ++j) {
          const int idxi1 = i + occ_actoff[h];
          const int idxj1 = j + occ_actoff[h];
          const int idxi2 = i + ribsoff[h];
          const int idxj2 = j + ribsoff[h];

          MX_nb_ribs.set_element(idxi2,idxj2,MX_occ_act.get_element(idxi1,idxj1));
          MY_nb_ribs.set_element(idxi2,idxj2,MY_occ_act.get_element(idxi1,idxj1));
          MZ_nb_ribs.set_element(idxi2,idxj2,MZ_occ_act.get_element(idxi1,idxj1));
        }
      }

      for (int a = 0; a < virpi[h]; ++a) {
        for (int b = 0; b <= a; ++b) {
          const int idxa1 = a + viroff[h];
          const int idxb1 = b + viroff[h];
          const int idxa2 = a + ribsoff[h] + occ_actpi[h];
          const int idxb2 = b + ribsoff[h] + occ_actpi[h];

          MX_nb_ribs.set_element(idxa2,idxb2,MX_vir.get_element(idxa1,idxb1));
          MY_nb_ribs.set_element(idxa2,idxb2,MY_vir.get_element(idxa1,idxb1));
          MZ_nb_ribs.set_element(idxa2,idxb2,MZ_vir.get_element(idxa1,idxb1));
        }
      }

      for (int a = 0; a < virpi[h]; ++a) {
        for (int i = 0; i < occ_actpi[h]; ++i) {
          const int idx_a1 = viroff[h] + a;
          const int idx_i1 = occ_actoff[h] + i;
          const int idx_a2 = ribsoff[h] + occ_actpi[h] + a;
          const int idx_i2 = ribsoff[h] + i;

          MX_nb_ribs.set_element(idx_a2,idx_i2,MX_virocc_act.get_element(idx_a1,idx_i1));
          MY_nb_ribs.set_element(idx_a2,idx_i2,MY_virocc_act.get_element(idx_a1,idx_i1));
          MZ_nb_ribs.set_element(idx_a2,idx_i2,MZ_virocc_act.get_element(idx_a1,idx_i1));
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
        for (int i = 0; i < occ_actpi[h]; ++i) {
          const int idx_ap1 = cabsoff[h] + ap;
          const int idx_i1 = occ_actoff[h] + i;
          const int idx_ap2 = ribsoff[h] + orbspi[h] + ap;
          const int idx_i2 = ribsoff[h] + i;

          MX_nb_ribs.set_element(idx_ap2,idx_i2,MX_cabsocc_act.get_element(idx_ap1,idx_i1));
          MY_nb_ribs.set_element(idx_ap2,idx_i2,MY_cabsocc_act.get_element(idx_ap1,idx_i1));
          MZ_nb_ribs.set_element(idx_ap2,idx_i2,MZ_cabsocc_act.get_element(idx_ap1,idx_i1));
        }
      }

      for (int ap = 0; ap < cabspi[h]; ++ap) {
        for (int a = 0; a < virpi[h]; ++a) {
          const int idx_ap1 = cabsoff[h] + ap;
          const int idx_a1 = viroff[h] + a;
          const int idx_ap2 = ribsoff[h] + orbspi[h] + ap;
          const int idx_a2 = ribsoff[h] + occ_actpi[h] + a;

          MX_nb_ribs.set_element(idx_ap2,idx_a2,MX_cabsvir.get_element(idx_ap1,idx_a1));
          MY_nb_ribs.set_element(idx_ap2,idx_a2,MY_cabsvir.get_element(idx_ap1,idx_a1));
          MZ_nb_ribs.set_element(idx_ap2,idx_a2,MZ_cabsvir.get_element(idx_ap1,idx_a1));
        }
      }

    }
    // clean
    MX_occ_act = 0;
    MY_occ_act = 0;
    MZ_occ_act = 0;
    MX_vir = 0;
    MY_vir = 0;
    MZ_vir = 0;
    MX_virocc_act = 0;
    MY_virocc_act = 0;
    MZ_virocc_act = 0;
    MX_cabsocc_act = 0;
    MY_cabsocc_act = 0;
    MZ_cabsocc_act = 0;
    MX_cabsvir = 0;
    MY_cabsvir = 0;
    MZ_cabsvir = 0;

    // Symmetrize matrices
    for(int p = 0; p < nribs_fc; p++) {
      for(int q = 0; q < p; q++) {
        MX_nb_ribs.set_element(q,p,MX_nb_ribs.get_element(p,q));
        MY_nb_ribs.set_element(q,p,MY_nb_ribs.get_element(p,q));
        MZ_nb_ribs.set_element(q,p,MZ_nb_ribs.get_element(p,q));
      }
    }

    const Ref<OrbitalSpace>& space_orbs = (spin == Alpha? orbs1 : orbs2);
    RefSCMatrix MX_orbs, MY_orbs, MZ_orbs;
    compute_multipole_ints(space_orbs, space_orbs,
                           MX_orbs, MY_orbs, MZ_orbs,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);

    RefSCMatrix MX_cabsorbs, MY_cabsorbs, MZ_cabsorbs;
    compute_multipole_ints(space_cabs, space_orbs,
                           MX_cabsorbs, MY_cabsorbs, MZ_cabsorbs,
                           MXX, MYY, MZZ,
                           MXY, MXZ, MYZ);
    //MX_cabsorbs.print(prepend_spincase(spin,"mu(X) in cabs_orbs").c_str());
    MXX = 0;
    MYY = 0;
    MZZ = 0;
    MXY = 0;
    MXZ = 0;
    MYZ = 0;

    const std::vector<unsigned int>& orbspi_nfc = space_orbs->block_sizes();

    std::vector<unsigned int> orbsoff_nfc(nirreps);
    std::vector<unsigned int> ribsoff_nfc(nirreps);
    orbsoff_nfc[0] = 0;
    ribsoff_nfc[0] = 0;
    for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
      orbsoff_nfc[irrep] = orbsoff_nfc[irrep-1] + orbspi_nfc[irrep-1];
      ribsoff_nfc[irrep] = orbsoff_nfc[irrep] + cabsoff[irrep];
    }

    for (unsigned int h = 0; h < nirreps; ++h) {

      for (int p = 0; p < orbspi_nfc[h]; ++p) {
        for (int q = 0; q <= p; ++q) {
          const int idxp1 = p + orbsoff_nfc[h];
          const int idxq1 = q + orbsoff_nfc[h];
          const int idxp2 = p + ribsoff_nfc[h];
          const int idxq2 = q + ribsoff_nfc[h];

          MX_ribs_nfc.set_element(idxp2,idxq2,MX_orbs.get_element(idxp1,idxq1));
          MY_ribs_nfc.set_element(idxp2,idxq2,MY_orbs.get_element(idxp1,idxq1));
          MZ_ribs_nfc.set_element(idxp2,idxq2,MZ_orbs.get_element(idxp1,idxq1));
        }
      }

      for (int ap = 0; ap < cabspi[h]; ++ap) {
        for (int bp = 0; bp <= ap; ++bp) {
          const int idx_ap1 = cabsoff[h] + ap;
          const int idx_bp1 = cabsoff[h] + bp;
           const int idx_ap2 = ribsoff_nfc[h] + orbspi_nfc[h] + ap;
           const int idx_bp2 = ribsoff_nfc[h] + orbspi_nfc[h] + bp;

           MX_ribs_nfc.set_element(idx_ap2,idx_bp2,MX_cabs.get_element(idx_ap1,idx_bp1));
           MY_ribs_nfc.set_element(idx_ap2,idx_bp2,MY_cabs.get_element(idx_ap1,idx_bp1));
           MZ_ribs_nfc.set_element(idx_ap2,idx_bp2,MZ_cabs.get_element(idx_ap1,idx_bp1));
         }
       }

      for (int ap = 0; ap < cabspi[h]; ++ap) {
        for (int p = 0; p < orbspi_nfc[h]; ++p) {
          const int idx_ap1 = cabsoff[h] + ap;
          const int idx_p1 = orbsoff_nfc[h] + p;
          const int idx_ap2 = ribsoff_nfc[h] + orbspi_nfc[h] + ap;
          const int idx_p2 = ribsoff_nfc[h] + p;

          MX_ribs_nfc.set_element(idx_ap2,idx_p2,MX_cabsorbs.get_element(idx_ap1,idx_p1));
          MY_ribs_nfc.set_element(idx_ap2,idx_p2,MY_cabsorbs.get_element(idx_ap1,idx_p1));
          MZ_ribs_nfc.set_element(idx_ap2,idx_p2,MZ_cabsorbs.get_element(idx_ap1,idx_p1));
        }
      }
    }

    MX_orbs = 0;
    MY_orbs = 0;
    MZ_orbs = 0;
    MX_cabs = 0;
    MY_cabs = 0;
    MZ_cabs = 0;
    MX_cabsorbs = 0;
    MY_cabsorbs = 0;
    MZ_cabsorbs = 0;

    // Symmetrize matrices
    for(int p = 0; p < nribs_nfc; p++) {
      for(int q = 0; q < p; q++) {
        MX_ribs_nfc.set_element(q,p,MX_ribs_nfc.get_element(p,q));
        MY_ribs_nfc.set_element(q,p,MY_ribs_nfc.get_element(p,q));
        MZ_ribs_nfc.set_element(q,p,MZ_ribs_nfc.get_element(p,q));
      }
    }
    //MZ_ribs_nfc.print(prepend_spincase(spin,"mu(Z)_ribs_nfc").c_str());

  }
  // end of form_DipoleInts_inRibs

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
  //const Ref<OrbitalSpace> cabs1 = r12eval()->cabs_space_canonical(spin1);

  const Ref<OrbitalSpace> occ1 = r12eval()->occ(spin1);

  const Ref<OrbitalSpace> occ2_act = r12eval()->occ_act(spin2);
  const Ref<OrbitalSpace> vir2 = r12eval()->vir(spin2);
  const Ref<OrbitalSpace> orbs2 = r12eval()->orbs(spin2);
  const Ref<OrbitalSpace> cabs2 = r12world->cabs_space(spin2);
  // use canonical cabs for test propose
  //const Ref<OrbitalSpace> cabs2 = r12eval()->cabs_space_canonical(spin2);

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

// compute MP2 one-electron density matrix contribution
RefSCMatrix MP2R12Energy_Diag::compute_1rdm_mp2_test(const SpinCase1 spin)
{
  const int nspincase2 = (r12eval()->spin_polarized() ? 3 : 2);

  // obtain occ, vir, orbs, cabs, and occ orbitals
  // spin=0 (Alpha) => AlphaAlpha case (1)
  // spin=1 (Beta) => BetaBeta case (2)
  const SpinCase2 spincase = static_cast<SpinCase2>(spin+1);
  vector< Ref<OrbitalSpace> > v_orbs1;
  vector< Ref<OrbitalSpace> > v_orbs2;
  obtain_orbitals(spincase, v_orbs1, v_orbs2);

  // get the numbers of the occupied and vir orbitals
  const Ref<OrbitalSpace> occ_act = v_orbs1[0];
  const Ref<OrbitalSpace> vir = v_orbs1[1];
  const int nocc_act = occ_act->rank();
  const int nvir = vir->rank();
  const int norbs = nocc_act + nvir;

  // get eigenvalues of Fock matrix
  const RefDiagSCMatrix evals_i = occ_act->evals();
  const RefDiagSCMatrix evals_a = vir->evals();

  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
  Ref<TwoBodyIntDescr> descr_f12 =
      r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(), 0);
  const TwoBodyOper::type eri_type =
      r12world->r12tech()->corrfactor()->tbint_type_eri();

  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const unsigned int eri_idx = descr_f12->intset(eri_type);

  // D^i_j =
  // same spin part: (g^ik_ab g^ab_jk - g^ik_ab g^ab_jk)
  //                / (F^i_i + F^k_k - F^a_a - F^b_b)
  //                / (F^j_j + F^k_k - F^a_a - F^b_b)
  // different spin part: g^ik_ab g^ab_jk
  //                    / (F^i_i + F^k_k - F^a_a - F^b_b)
  //                    / (F^j_j + F^k_k - F^a_a - F^b_b)

  // AlphaAlpha or BetaBeta ints
  Ref<DistArray4> iiaa_ints;
  activate_ints(occ_act->id(), occ_act->id(), vir->id(),
                vir->id(), descr_f12_key, moints4_rtime,
                iiaa_ints);

  // AlphaBeta or BetaAlpha ints
  Ref<OrbitalSpace> occ1_act;
  Ref<OrbitalSpace> occ2_act;
  Ref<OrbitalSpace> vir1;
  Ref<OrbitalSpace> vir2;
  int nocc1_act = 0;
  int nocc2_act = 0;
  int nvir1 = 0;
  int nvir2 = 0;
  Ref<DistArray4> i1i2a1a2_ints;
  Ref<DistArray4> a1a2i1i2_ints;
  RefDiagSCMatrix evals_i1;
  RefDiagSCMatrix evals_i2;
  RefDiagSCMatrix evals_a1;
  RefDiagSCMatrix evals_a2;
  if (nspincase2 == 3) {
      // obtain AlphaBeta/BetaAlpha occ, vir, orbs, cabs, and occ orbitals
      // for alpha case: v_orbs1_ab has the orbital in alpha space
      // for beta case: v_orbs1_ab has the orbital in beta space
      vector< Ref<OrbitalSpace> > v_orbs1_ab;
      vector< Ref<OrbitalSpace> > v_orbs2_ab;
      if (spin == Alpha) {
        obtain_orbitals(AlphaBeta, v_orbs1_ab, v_orbs2_ab);
      } else {
          obtain_orbitals(AlphaBeta, v_orbs2_ab, v_orbs1_ab);
      }
      // get the numbers of the occupied alpha and beta orbitals
      occ1_act = v_orbs1_ab[0];
      occ2_act = v_orbs2_ab[0];
      vir1 = v_orbs1_ab[1];
      vir2 = v_orbs2_ab[1];
      nocc1_act = occ1_act->rank();
      nocc2_act = occ2_act->rank();
      nvir1 = vir1->rank();
      nvir2 = vir2->rank();

      activate_ints(occ1_act->id(), occ2_act->id(), vir1->id(),
                    vir2->id(), descr_f12_key, moints4_rtime,
                    i1i2a1a2_ints);
      activate_ints(vir1->id(), vir2->id(), occ1_act->id(),
                    occ2_act->id(), descr_f12_key, moints4_rtime,
                    a1a2i1i2_ints);

      evals_i1 = occ1_act->evals();
      evals_i2 = occ2_act->evals();
      evals_a1 = vir1->evals();
      evals_a2 = vir2->evals();
  }

  RefSCDimension rowdim = new SCDimension(norbs);
  RefSCDimension coldim = new SCDimension(norbs);
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

  RefSCMatrix D_mp2 = localkit->matrix(rowdim, coldim);
  D_mp2.assign(0.0);

  for (int i = 0; i < nocc_act; ++i) {
    const double Fii = evals_i(i);

    for (int j = 0; j <= i; ++j) {
      const double Fjj = evals_i(j);
      double Dmp2_ij_1 = 0;
      double Dmp2_ij_2 = 0;
      double Dmp2_ij_ab = 0;

      // AlphaAlpha or BetaBeta part
      for (int k = 0; k < nocc_act; ++k) {
        const double* gik_ab = iiaa_ints->retrieve_pair_block(i, k, eri_idx);
        const double* gki_ab = iiaa_ints->retrieve_pair_block(k, i, eri_idx);
        const double* gjk_ab = iiaa_ints->retrieve_pair_block(j, k, eri_idx);
        const double Fkk = evals_i(k);

        for (int a = 0; a < nvir; ++a) {
          const double Faa = evals_a(a);

          for (int b = 0; b < nvir; ++b) {
            const double Fbb = evals_a(b);
            const double denom = 1.0 / ((Fii + Fkk - Faa - Fbb)
                                     * (Fjj + Fkk - Faa - Fbb));
            //Dmp2_ij += ((*gik_ab - *gki_ab)*denom) * (*gjk_ab)*denom;
            Dmp2_ij_1 += (*gik_ab) * (*gjk_ab) * denom;
            Dmp2_ij_2 += (*gki_ab) * (*gjk_ab) * denom;
            ++gik_ab;
            ++gki_ab;
            ++gjk_ab;
          }
        }
        iiaa_ints->release_pair_block(i, k, eri_idx);
        iiaa_ints->release_pair_block(k, i, eri_idx);
        iiaa_ints->release_pair_block(j, k, eri_idx);
      }

      // AlphaBeta part
      if (nspincase2 == 3) {
        const double Fii = evals_i1(i);
        const double Fjj = evals_i1(j);

        for (int k = 0; k < nocc2_act; ++k) {
          const double* gik_ab = i1i2a1a2_ints->retrieve_pair_block(i, k, eri_idx);
          const double* gjk_ab = i1i2a1a2_ints->retrieve_pair_block(j, k, eri_idx);
          const double Fkk = evals_i2(k);

          for (int a = 0; a < nvir1; ++a) {
            const double Faa = evals_a1(a);

            for (int b = 0; b < nvir2; ++b) {
              const double Fbb = evals_a2(b);
              const double denom = 1.0 / ((Fii + Fkk - Faa - Fbb)
                                        * (Fjj + Fkk - Faa - Fbb));
              Dmp2_ij_ab += (*gik_ab) * (*gjk_ab) * denom;
              ++gik_ab;
              ++gjk_ab;
            }
          }
          i1i2a1a2_ints->release_pair_block(i, k, eri_idx);
          i1i2a1a2_ints->release_pair_block(j, k, eri_idx);
        }

      } else {
          Dmp2_ij_ab = Dmp2_ij_1;
      }

      const double Dmp2_ij = -(Dmp2_ij_1 - Dmp2_ij_2 + Dmp2_ij_ab);
      D_mp2.set_element(i, j, Dmp2_ij);
    }
  }
  iiaa_ints->deactivate();
  if (nspincase2 == 3) i1i2a1a2_ints->deactivate();

  // D^a_b
  Ref<DistArray4> aaii_ints;
  activate_ints(vir->id(), vir->id(), occ_act->id(),
                occ_act->id(), descr_f12_key, moints4_rtime,
                aaii_ints);

  for (int a = 0; a < nvir; ++a) {
    for (int b = 0; b <= a; ++b) {
      double Dmp2_ab_1 = 0;
      double Dmp2_ab_2 = 0;
      double Dmp2_ab_ab = 0;

      // AlphaAlpha or BetaBeta part
      for (int c = 0; c < nvir; ++c) {
        const double* gac_ij = aaii_ints->retrieve_pair_block(a, c, eri_idx);
        const double* gca_ij = aaii_ints->retrieve_pair_block(c, a, eri_idx);
        const double* gbc_ij = aaii_ints->retrieve_pair_block(b, c, eri_idx);

        const double Faa = evals_a(a);
        const double Fbb = evals_a(b);
        const double Fcc = evals_a(c);

//        double* T2 = new double[nocc_act * nocc_act];
//        double* iter_T2 = T2;
//        fill_n(T2, nocc_act * nocc_act, 0.0);
//        ExEnv::out0() << endl << a << " " << c << endl;

        for (int i = 0; i < nocc_act; ++i) {
          for (int j = 0; j < nocc_act; ++j) {
            const double Fii = evals_i(i);
            const double Fjj = evals_i(j);
            const double denom = 1.0 / ((Fii + Fjj - Faa - Fcc)
                                      * (Fii + Fjj - Fbb - Fcc));
            Dmp2_ab_1 += (*gac_ij) * (*gbc_ij) * denom;
            Dmp2_ab_2 += (*gca_ij) * (*gbc_ij) * denom;
//            *iter_T2 = (*gac_ij)
//                         / (Fii + Fjj - Faa - Fcc);
//            ++iter_T2;

            ++gac_ij;
            ++gca_ij;
            ++gbc_ij;
          }
        }
//        print_intermediate("AlphAlpha", "T test", T2, nocc_act, nocc_act);
//        delete[] T2;
        aaii_ints->release_pair_block(a, c, eri_idx);
        aaii_ints->release_pair_block(c, a, eri_idx);
        aaii_ints->release_pair_block(b, c, eri_idx);
      }

      // alpha beta part
      if (nspincase2 == 3) {
        const double Faa = evals_a1(a);
        const double Fbb = evals_a1(b);

        for (int c = 0; c < nvir2; ++c) {
          const double* gac_ab = a1a2i1i2_ints->retrieve_pair_block(a, c, eri_idx);
          const double* gbc_ab = a1a2i1i2_ints->retrieve_pair_block(b, c, eri_idx);
          const double Fcc = evals_a2(c);

          for (int i = 0; i < nocc1_act; ++i) {
            for (int j = 0; j < nocc2_act; ++j) {
              const double Fii = evals_i1(i);
              const double Fjj = evals_i2(j);

              const double denom = 1.0 / ((Fii + Fjj - Faa - Fcc)
                                        * (Fii + Fjj - Fbb - Fcc));
              Dmp2_ab_ab += (*gac_ab) * (*gbc_ab) * denom;
              ++gac_ab;
              ++gbc_ab;
            }
          }
          a1a2i1i2_ints->release_pair_block(a, c, eri_idx);
          a1a2i1i2_ints->release_pair_block(b, c, eri_idx);
        }

        } else {
            Dmp2_ab_ab = Dmp2_ab_1;
        }

      const double Dmp2_ab = Dmp2_ab_1 - Dmp2_ab_2 + Dmp2_ab_ab;
      D_mp2.set_element(a+nocc_act, b+nocc_act, Dmp2_ab);
    }
  }
  aaii_ints->deactivate();
  if (nspincase2 == 3) a1a2i1i2_ints->deactivate();

  // symmetrize matrix
  for (int p = 0; p < norbs; ++p) {
    for (int q = 0; q < p; ++q) {
      D_mp2.set_element(q, p, D_mp2.get_element(p,q));
    }
  }

  return D_mp2;
}
// end of function: compute_1rdm_mp2

// compute MP2 one-electron density matrix contribution
RefSCMatrix MP2R12Energy_Diag::compute_1rdm_mp2(const SpinCase1 spin)
{
  const int nspincase2 = (r12eval()->spin_polarized() ? 3 : 2);

  // obtain occ, vir, orbs, cabs, and occ orbitals
  // spin=0 (Alpha) => AlphaAlpha case (1)
  // spin=1 (Beta) => BetaBeta case (2)
  const SpinCase2 spincase = static_cast<SpinCase2>(spin+1);
  vector< Ref<OrbitalSpace> > v_orbs1;
  vector< Ref<OrbitalSpace> > v_orbs2;
  obtain_orbitals(spincase, v_orbs1, v_orbs2);

  // get the numbers of the occupied and vir orbitals
  const Ref<OrbitalSpace> occ_act = v_orbs1[0];
  const Ref<OrbitalSpace> vir = v_orbs1[1];
  const int nocc_act = occ_act->rank();
//  if (nocc_act == 0)
//    continue;
  const int nvir = vir->rank();

  // D^i_j =
  // same spin part: T^ik_ab T^ab_jk - T^ik_ab T^ab_jk
  // different spin part: T^ik_ab T^ab_jk
  const int norbs = nocc_act + nvir;
  const blasint noccocc = nocc_act * nocc_act;
  const int nvirvir = nvir * nvir;
  const int noccoccvir = nvir * nocc_act * nocc_act;

  double* T2ab_ij = new double[nvirvir * noccocc];
  fill_n(T2ab_ij, nvirvir * noccocc, 0.0);
  compute_T2_mp2(v_orbs1, v_orbs1, T2ab_ij);

//  // test for f12 corrected MP2 1e density matrix
//  const double C_0 = 1.0 / 2.0;
//  const double C_1 = 1.0 / 4.0;
//  compute_T2abij_mp2f12(spincase, C_0, C_1, T2ab_ij);

  // AlphaBeta or BetaAlpha ints
  double* T2ab_ij_ab;
  Ref<OrbitalSpace> occ1_act;
  Ref<OrbitalSpace> occ2_act;
  Ref<OrbitalSpace> vir1;
  Ref<OrbitalSpace> vir2;
  int nocc1_act = 0;
  int nocc2_act = 0;
  int nvir1 = 0;
  int nvir2 = 0;
  blasint nocc12 = 0;
  int nocc12vir = 0;
  if (nspincase2 == 3) {
    // obtain AlphaBeta/BetaAlpha occ, vir, orbs, cabs, and occ orbitals
    // for alpha case: v_orbs1_ab has the orbital in alpha space
    // for beta case: v_orbs1_ab has the orbital in beta space
    vector< Ref<OrbitalSpace> > v_orbs1_ab;
    vector< Ref<OrbitalSpace> > v_orbs2_ab;
    if (spin == Alpha) {
      obtain_orbitals(AlphaBeta, v_orbs1_ab, v_orbs2_ab);
    } else {
        obtain_orbitals(AlphaBeta, v_orbs2_ab, v_orbs1_ab);
    }
    // get the numbers of the occupied alpha and beta orbitals
    occ1_act = v_orbs1_ab[0];
    occ2_act = v_orbs2_ab[0];
    vir1 = v_orbs1_ab[1];
    vir2 = v_orbs2_ab[1];
    nocc1_act = occ1_act->rank();
    nocc2_act = occ2_act->rank();
    nvir1 = vir1->rank();
    nvir2 = vir2->rank();

    nocc12 = nocc1_act * nocc2_act;
    const int nvir12 = nvir1 * nvir2;
    nocc12vir = nvir2 * nocc12;
    T2ab_ij_ab = new double[nvir12 * nocc12];
    fill_n(T2ab_ij_ab, nvir12 * nocc12, 0.0);
    compute_T2_mp2(v_orbs1_ab, v_orbs2_ab, T2ab_ij_ab);
//    compute_T2abij_mp2f12(AlphaBeta, C_0, C_1, T2ab_ij_ab);
  }

  RefSCDimension rowdim = new SCDimension(norbs);
  RefSCDimension coldim = new SCDimension(norbs);
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

  RefSCMatrix D_mp2 = localkit->matrix(rowdim, coldim);
  D_mp2.assign(0.0);

  // D^i_j
  // same spin part: T^ik_ab T^ab_jk - T^ki_ab T^ab_jk
  // different spin part: T^ik_ab T^ab_jk
  for (int i = 0; i < nocc_act; ++i) {
    for (int j = 0; j <= i; ++j) {
      // AlphaAlpha or BetaBeta part
      double Dij_1 = 0;
      double Dij_2 = 0;

      for (int k = 0; k < nocc_act; ++k) {
        const double* iter_T2_ikab = T2ab_ij + i * nocc_act + k;
        const double* iter_T2_kiab = T2ab_ij + k * nocc_act + i;
        const double* iter_T2_jkab = T2ab_ij + j * nocc_act + k;

        for (int a = 0; a < nvir; ++a) {
          for (int b = 0; b < nvir; ++b) {
            Dij_1 += (*iter_T2_ikab) * (*iter_T2_jkab);
            Dij_2 += (*iter_T2_kiab) * (*iter_T2_jkab);

            iter_T2_ikab += noccocc;
            iter_T2_kiab += noccocc;
            iter_T2_jkab += noccocc;
          }
        }
      }

      // AlphaBeta part
      double Dij_3 = 0;
      if (nspincase2 == 3) {
        for (int k = 0; k < nocc2_act; ++k) {
          const double* iter_T2_ikab = T2ab_ij_ab + i * nocc_act + k;
          const double* iter_T2_jkab = T2ab_ij_ab + j * nocc_act + k;

          for (int a = 0; a < nvir1; ++a) {
            for (int b = 0; b < nvir2; ++b) {
              Dij_3 += (*iter_T2_ikab) * (*iter_T2_jkab);

              iter_T2_ikab += nocc12;
              iter_T2_jkab += nocc12;
            }
          }
        }
      } else {
          Dij_3 = Dij_1;
      }

      const double Dmp2_ij = -(Dij_1 - Dij_2 + Dij_3);
      D_mp2.set_element(i, j, Dmp2_ij);
    }
  }

  // D^a_b
  // same spin part: T^ac_ij T^ij_bc - T^ca_ij T^ij_bc
  // different spin part: T^ac_ij T^ij_bc
  const blasint one = 1; // for F77_DDOT
  for (int a = 0; a < nvir; ++a) {

    const double* iter_T2_bcij = T2ab_ij;
    for (int b = 0; b <= a; ++b) {
     // AlphaAlpha or BetaBeta part
      double Dab_1 = 0;
      double Dab_2 = 0;

      const double* iter_T2_acij = T2ab_ij + a * noccoccvir;
      const double* iter_T2_caij = T2ab_ij + a * noccocc;
      for (int c = 0; c < nvir; ++c) {
        const double TTac_bc = F77_DDOT(&noccocc, iter_T2_acij, &one, iter_T2_bcij, &one);
        const double TTca_bc = F77_DDOT(&noccocc, iter_T2_caij, &one, iter_T2_bcij, &one);
        Dab_1 += TTac_bc;
        Dab_2 += TTca_bc;

        iter_T2_acij += noccocc;
        iter_T2_bcij += noccocc;
        iter_T2_caij += noccoccvir;
      }

      // alpha beta part
      double Dab_3 = 0;
      if (nspincase2 == 3) {
        const double* iter_T2_acij_ab = T2ab_ij_ab + a * nocc12vir;
        const double* iter_T2_bcij_ab = T2ab_ij_ab + b * nocc12vir;
        for (int c = 0; c < nvir2; ++c) {
          const double TTac_bc = F77_DDOT(&nocc12, iter_T2_acij_ab, &one, iter_T2_bcij_ab, &one);
          Dab_3 += TTac_bc;

          iter_T2_acij_ab += nocc12;
          iter_T2_bcij_ab += nocc12;
        }
      } else {
          Dab_3 = Dab_1;
      }

      const double Dmp2_ab = Dab_1 - Dab_2 + Dab_3;
      D_mp2.set_element(a+nocc_act, b+nocc_act, Dmp2_ab);
    }
  }
  delete[] T2ab_ij;
  if (nspincase2 == 3) delete[] T2ab_ij_ab;

  // symmetrize matrix
  for (int p = 0; p < norbs; ++p) {
    for (int q = 0; q < p; ++q) {
      D_mp2.set_element(q, p, D_mp2.get_element(p,q));
    }
  }

  return D_mp2;
}
// end of function: compute_1rdm_mp2

// compute contribution for MP2F12 one-electron density matrix
// D_MP2F12 = T(MP2)T(MP2) + T(MP2)T(F12) + T(F12)T(MP2) + T(F12)T(F12)
// compute T(MP2)T(MP2) or T(F12)T(F12), whose left and right sides are the same
RefSCMatrix MP2R12Energy_Diag::compute_1rdm_mp2part(const SpinCase1 spin,
                                                    const int nocc1_act, const int nocc2_act,
                                                    const int nvir1, const int nvir2,
                                                    const double* const T2,
                                                    const double* const T2_ab)
{
  // get the numbers of the occupied and vir orbitals
  const int nocc_act = (spin == Alpha? nocc1_act : nocc2_act);
  const int nvir = (spin == Alpha? nvir1 : nvir2);
//  if (nocc_act == 0)
//    continue;
  const int norbs = nocc_act + nvir;
  const blasint noccocc = nocc_act * nocc_act;
  const int nvirvir = nvir * nvir;
  const int noccoccvir = nocc_act * nocc_act * nvir;

  const blasint nocc12 = nocc1_act * nocc2_act;
  const int nocc12vir2 = nocc12 * nvir2;

  RefSCDimension rowdim = new SCDimension(norbs);
  RefSCDimension coldim = new SCDimension(norbs);
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

  RefSCMatrix D = localkit->matrix(rowdim, coldim);
  D.assign(0.0);

  // D^i_j
  // same spin part: 1/2 T^ik_ab T^ab_jk (T is antisymmetrized)
  // different spin part: T^i1k2_a1b2 T^a1b2_j1k2
  //                   or T^k1i2_a1b2 T^a1b2_k1j2

  // preliminaries for different spin part
  const int ij_offset = (spin == Alpha? nocc2_act : 1);
  const int k_offset = (spin == Alpha? 1 : nocc2_act);
  const int size_k = (spin == Alpha? nocc2_act : nocc1_act);

  const int ab_offset = (spin == Alpha? nocc12vir2 : nocc12);
  const int c_offset = (spin == Alpha? nocc12 : nocc12vir2);
  const int size_c = (spin == Alpha? nvir2 : nvir1);

  for (int i = 0; i < nocc_act; ++i) {
    for (int j = 0; j <= i; ++j) {

      // AlphaAlpha or BetaBeta part: T^ik_ab T^ab_jk
      double Dij_1 = 0;
      for (int k = 0; k < nocc_act; ++k) {
        const double* iter_T2_ikab = T2 + i * nocc_act + k;
        const double* iter_T2_jkab = T2 + j * nocc_act + k;

        for (int a = 0; a < nvir; ++a) {
          for (int b = 0; b < nvir; ++b) {
            Dij_1 += (*iter_T2_ikab) * (*iter_T2_jkab);

            iter_T2_ikab += noccocc;
            iter_T2_jkab += noccocc;
          }
        }
      }

      // AlphaBeta part
      double Dij_2 = 0;
      // Alpha: T^i1k2_a1b2 T^a1b2_j1k2
      // Beta:  T^k1i2_a1b2 T^a1b2_k1j2
      for (int k = 0; k < size_k; ++k) {
        const double* iter_T2_ik = T2_ab + i * ij_offset + k * k_offset;   // T2^i1k2_a1b2 or T2^k1i2_a1b2
        const double* iter_T2_jk = T2_ab + j * ij_offset + k * k_offset;  // T2^j1k2_a1b2 or T2^k1j2_a1b2

        for (int a = 0; a < nvir1; ++a) {
          for (int b = 0; b < nvir2; ++b) {
            Dij_2 += (*iter_T2_ik) * (*iter_T2_jk);

            iter_T2_ik += nocc12;
            iter_T2_jk += nocc12;
          }
        }
      }

      const double Dij = -(Dij_1 * 0.5 + Dij_2);
      D.set_element(i, j, Dij);
    }
  }

  // D^a_b
  // same spin part: 1/2 T^ac_ij T^ij_bc (T is antisymmetrized)
  // different spin part: T^a1c2_i1j2 T^i1j2_b1c2
  //                   or T^c1a2_i1j2 T^i1j2_c1b2
  const blasint one = 1; // for F77_DDOT
  for (int a = 0; a < nvir; ++a) {
    const double* iter_T2_bcij = T2;

    for (int b = 0; b <= a; ++b) {
      const double* iter_T2_acij = T2 + a * noccoccvir;

     // AlphaAlpha or BetaBeta part: T^ac_ij T^ij_bc
      double Dab_1 = 0;
      for (int c = 0; c < nvir; ++c) {
        const double TTac_bc = F77_DDOT(&noccocc, iter_T2_acij, &one, iter_T2_bcij, &one);
        Dab_1 += TTac_bc;

        iter_T2_acij += noccocc;
        iter_T2_bcij += noccocc;
      }

      // AlphaBeta part
      double Dab_2 = 0;
      //Alpha: T^a1c2_i1j2 T^i1j2_b1c2
      //Beta: T^c1a2_i1j2 T^i1j2_c1b2

      const double* iter_T2_ac = T2_ab + a * ab_offset; // T^a1c2_i1j2 or T^c1a2_i1j2
      const double* iter_T2_bc = T2_ab + b * ab_offset;  // T^i1j2_b1c2 or T^i1j2_c1b2
      for (int c = 0; c < size_c; ++c) {
        const double TTac_bc = F77_DDOT(&nocc12, iter_T2_ac, &one, iter_T2_bc, &one);
        Dab_2 += TTac_bc;

        iter_T2_ac += c_offset;
        iter_T2_bc += c_offset;
      }
      const double Dab = Dab_1 * 0.5 + Dab_2;
      D.set_element(a+nocc_act, b+nocc_act, Dab);
    }
  }

  // symmetrize matrix
  for (int p = 0; p < norbs; ++p) {
    for (int q = 0; q < p; ++q) {
      D.set_element(q, p, D.get_element(p,q));
    }
  }

  return D;
}
// end of function: compute_1rdm_mp2part

// compute T(MP2)T(F12) or T(MP2)T(F12)
RefSCMatrix MP2R12Energy_Diag::compute_1rdm_mp2part(const SpinCase1 spin,
                                                    const int nocc1_act, const int nocc2_act,
                                                    const int nvir1, const int nvir2,
                                                    const double* const T2_left,
                                                    const double* const T2_right,
                                                    const double* const T2_ab_left,
                                                    const double* const T2_ab_right)
{
  // get the numbers of the occupied and vir orbitals
  const int nocc_act = (spin == Alpha? nocc1_act : nocc2_act);
  const int nvir = (spin == Alpha? nvir1 : nvir2);
//  if (nocc_act == 0)
//    continue;
  const int norbs = nocc_act + nvir;
  const blasint noccocc = nocc_act * nocc_act;
  const int nvirvir = nvir * nvir;
  const int noccoccvir = nocc_act * nocc_act * nvir;

  const blasint nocc12 = nocc1_act * nocc2_act;
  const int nocc12vir2 = nocc12 * nvir2;

  RefSCDimension rowdim = new SCDimension(norbs);
  RefSCDimension coldim = new SCDimension(norbs);
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

  RefSCMatrix D = localkit->matrix(rowdim, coldim);
  D.assign(0.0);

  // D^i_j
  // same spin part: 1/2 T^ik_ab T^ab_jk (T is antisymmetrized)
  // different spin part: T^i1k2_a1b2 T^a1b2_j1k2
  //                   or T^k1i2_a1b2 T^a1b2_k1j2

  // preliminaries for different spin part
  const int ij_offset = (spin == Alpha? nocc2_act : 1);
  const int k_offset = (spin == Alpha? 1 : nocc2_act);
  const int size_k = (spin == Alpha? nocc2_act : nocc1_act);

  const int ab_offset = (spin == Alpha? nocc12vir2 : nocc12);
  const int c_offset = (spin == Alpha? nocc12 : nocc12vir2);
  const int size_c = (spin == Alpha? nvir2 : nvir1);

  for (int i = 0; i < nocc_act; ++i) {
    for (int j = 0; j < nocc_act; ++j) {

      // AlphaAlpha or BetaBeta part: T^ik_ab T^ab_jk
      double Dij_1 = 0;
      for (int k = 0; k < nocc_act; ++k) {
        const double* iter_T2_ikab = T2_left + i * nocc_act + k;
        const double* iter_T2_jkab = T2_right + j * nocc_act + k;

        for (int a = 0; a < nvir; ++a) {
          for (int b = 0; b < nvir; ++b) {
            Dij_1 += (*iter_T2_ikab) * (*iter_T2_jkab);

            iter_T2_ikab += noccocc;
            iter_T2_jkab += noccocc;
          }
        }
      }

      // AlphaBeta part
      double Dij_2 = 0;
      // Alpha: T^i1k2_a1b2 T^a1b2_j1k2
      // Beta:  T^k1i2_a1b2 T^a1b2_k1j2
      for (int k = 0; k < size_k; ++k) {
        const double* iter_T2_ik = T2_ab_left + i * ij_offset + k * k_offset;   // T2^i1k2_a1b2 or T2^k1i2_a1b2
        const double* iter_T2_jk = T2_ab_right + j * ij_offset + k * k_offset;  // T2^j1k2_a1b2 or T2^k1j2_a1b2

        for (int a = 0; a < nvir1; ++a) {
          for (int b = 0; b < nvir2; ++b) {
            Dij_2 += (*iter_T2_ik) * (*iter_T2_jk);

            iter_T2_ik += nocc12;
            iter_T2_jk += nocc12;
          }
        }
      }

      const double Dij = -(Dij_1 * 0.5 + Dij_2);
      D.set_element(i, j, Dij);
    }
  }

  // D^a_b
  // same spin part: 1/2 T^ac_ij T^ij_bc (T is antisymmetrized)
  // different spin part: T^a1c2_i1j2 T^i1j2_b1c2
  //                   or T^c1a2_i1j2 T^i1j2_c1b2
  const blasint one = 1; // for F77_DDOT
  for (int a = 0; a < nvir; ++a) {
    const double* iter_T2_bcij = T2_right;

    for (int b = 0; b < nvir; ++b) {
      const double* iter_T2_acij = T2_left + a * noccoccvir;

     // AlphaAlpha or BetaBeta part: T^ac_ij T^ij_bc
      double Dab_1 = 0;
      for (int c = 0; c < nvir; ++c) {
        const double TTac_bc = F77_DDOT(&noccocc, iter_T2_acij, &one, iter_T2_bcij, &one);
        Dab_1 += TTac_bc;

        iter_T2_acij += noccocc;
        iter_T2_bcij += noccocc;
      }

      // AlphaBeta part
      double Dab_2 = 0;
      //Alpha: T^a1c2_i1j2 T^i1j2_b1c2
      //Beta: T^c1a2_i1j2 T^i1j2_c1b2

      const double* iter_T2_ac = T2_ab_left + a * ab_offset; // T^a1c2_i1j2 or T^c1a2_i1j2
      const double* iter_T2_bc = T2_ab_right + b * ab_offset;  // T^i1j2_b1c2 or T^i1j2_c1b2
      for (int c = 0; c < size_c; ++c) {
        const double TTac_bc = F77_DDOT(&nocc12, iter_T2_ac, &one, iter_T2_bc, &one);
        Dab_2 += TTac_bc;

        iter_T2_ac += c_offset;
        iter_T2_bc += c_offset;
      }
      const double Dab = Dab_1 * 0.5 + Dab_2;
      D.set_element(a+nocc_act, b+nocc_act, Dab);
    }
  }

  return D;
}
// end of function: compute_1rdm_mp2part

// compute D_MP2F12 = T(MP2)T(MP2) + 2 T(MP2)T(F12) + T(F12)T(F12)
//void MP2R12Energy_Diag::compute_1rdm_mp2f12(const int nspincases1,const int nspincases2,
//                                            const int C_0, const int C_1,
//                                            RefSCMatrix Dmp2f12[NSpinCases1])
//{
//  int norbs = 0;
//  vector<double*> T2_mp2(NSpinCases2);
//  vector<double*> T2_f12corr(NSpinCases2);
//  for(int s = 0; s < nspincases2; ++s) {
//    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
//
//    const SpinCase1 spin1 = case1(spincase2);
//    const SpinCase1 spin2 = case2(spincase2);
//
//    const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
//    const Ref<OrbitalSpace>& vir1 = r12eval()->vir(spin1);
//    const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
//    const Ref<OrbitalSpace>& vir2 = r12eval()->vir(spin2);
//
//    const int nocc1_act = occ1_act->rank();
//    const int nocc2_act = occ2_act->rank();
//    const int nvir1 = vir1->rank();
//    const int nvir2 = vir2->rank();
//
//    norbs = nocc1_act + nvir1;
//    const int nocc12 = nocc1_act * nocc2_act;
//    const int nvir12 = nvir1 * nvir2;
//    const int nocc12vir12 = nocc12 * nvir12;
//
//    if (nocc1_act == 0 || nocc2_act == 0)
//      continue;
//
//    T2_mp2[s] = new double[nocc12vir12];
//    T2_f12corr[s] = new double[nocc12vir12];
//    fill_n(T2_mp2[s], nocc12vir12, 0.0);
//    fill_n(T2_f12corr[s], nocc12vir12, 0.0);
//
//    compute_T2_mp2(spincase2, T2_mp2[s]);
//    compute_T2abij_f12corr(spincase2, C_0, C_1, T2_f12corr[s]);
//  }
//
//  if (nspincases1 == 1) {
//    T2_mp2[BetaBeta] = T2_mp2[AlphaAlpha];
//    T2_f12corr[BetaBeta] = T2_f12corr[AlphaAlpha];
//  }
//
////  RefSCMatrix Dmp2f12[NSpinCases1];
//  for (int s = 0; s < nspincases1; ++s) {
//    const SpinCase1 spin = static_cast<SpinCase1>(s);
//    const SpinCase2 spincase2 = static_cast<SpinCase2>(s+1);
//
//    // T(MP2)T(MP2)
//    Dmp2f12[spin] = compute_1rdm_mp2part(spin, nocc1_act, nocc2_act, nvir1, nvir2,
//                                         T2_mp2[AlphaAlpha], T2_mp2[AlphaBeta]);
////    const int norbs = Dmp2_mp2.nrow();
////    RefSCDimension rowdim = new SCDimension(norbs);
////    RefSCDimension coldim = new SCDimension(norbs);
////    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
////
////    Dmp2f12[s] = localkit->matrix(rowdim, coldim);
////    Dmp2f12[s].assign(0.0);
////
////    Dmp2f12[spin].accumulate_subblock(Dmp2_mp2, 0, norb-1, 0, norb-1, 0, 0);
////    Dmp2_mp2 = 0;
//
//    // + 2 * T(MP2)T(F12)
//    RefSCMatrix Dmp2_f12 = compute_1rdm_mp2part(spin, nocc1_act, nocc2_act, nvir1, nvir2,
//                                                T2_mp2[spincase2], T2_f12corr[spincase2],
//                                                T2_mp2[AlphaBeta], T2_f12corr[AlphaBeta]);
//    Dmp2_f12.scale(2);
//    Dmp2f12[spin].accumulate_subblock(Dmp2_f12, 0, norbs-1, 0, norbs-1, 0, 0);
//    Dmp2_f12 = 0;
//
//    // + T(F12)T(F12)
//    RefSCMatrix Df12_f12 = compute_1rdm_mp2part(spin, nocc1_act, nocc2_act, nvir1, nvir2,
//                                                T2_f12corr[spincase2], T2_f12corr[AlphaBeta]);
//    Dmp2f12[spin].accumulate_subblock(Df12_f12, 0, norbs-1, 0, norbs-1, 0, 0);
//    Df12_f12 = 0;
//  }
//  if (nspincases1 == 1) {
//    Dmp2f12[Beta] = Dmp2f12[Alpha];
//  }
//
//  for(int s = 0; s < nspincases2; ++s) {
//    delete[] T2_mp2[s];
//    delete[] T2_f12corr[s];
//  }
//}
// end of compute_1rdm_mp2f12

//// compute MP2 one-electron density matrix contribution due to
//// the f12 correction
//RefSCMatrix MP2R12Energy_Diag::compute_1rdm_mp2_f12corr(const SpinCase1 spin)
//{
//  const int nspincase2 = (r12eval()->spin_polarized() ? 3 : 2);
//
//  // obtain occ, vir, orbs, cabs, and occ orbitals
//  // spin=0 (Alpha) => AlphaAlpha case (1)
//  // spin=1 (Beta) => BetaBeta case (2)
//  const SpinCase2 spincase = static_cast<SpinCase2>(spin+1);
//  vector< Ref<OrbitalSpace> > v_orbs1;
//  vector< Ref<OrbitalSpace> > v_orbs2;
//  obtain_orbitals(spincase, v_orbs1, v_orbs2);
//
//  // get the numbers of the occupied and vir orbitals
//  const Ref<OrbitalSpace> occ_act = v_orbs1[0];
//  const Ref<OrbitalSpace> vir = v_orbs1[1];
//  const int nocc_act = occ_act->rank();
////  if (nocc_act == 0)
////    continue;
//  const int nvir = vir->rank();
//
//  // D^i_j =
//  // same spin part: T^ik_ab T^ab_jk - T^ik_ab T^ab_jk
//  // different spin part: T^ik_ab T^ab_jk
//  const int norbs = nocc_act + nvir;
//  const int noccocc = nocc_act * nocc_act;
//  const int nvirvir = nvir * nvir;
//  const int noccoccvir = noccocc * nvir;
//
//  double* T2ab_ij = new double[nvirvir * noccocc];
//  fill_n(T2ab_ij, nvirvir * noccocc, 0.0);
//  compute_T2_mp2(v_orbs1, v_orbs1, T2ab_ij);
//
//  double* T2_f12corr = new double[nvirvir * noccocc];
//  fill_n(T2_f12corr, nvirvir * noccocc, 0.0);
//  compute_T2abij_f12corr(spincase2, C_0, C_1, T2_f12corr);
//
//  // AlphaBeta or BetaAlpha ints
//  double* T2ab_ij_ab;
//  double* T2_f12corr_ab;
//  Ref<OrbitalSpace> occ1_act;
//  Ref<OrbitalSpace> occ2_act;
//  Ref<OrbitalSpace> vir1;
//  Ref<OrbitalSpace> vir2;
//  int nocc1_act = 0;
//  int nocc2_act = 0;
//  int nvir1 = 0;
//  int nvir2 = 0;
//  int nocc12 = 0;
//  int nocc12vir = 0;
//  if (nspincase2 == 3) {
//    // obtain AlphaBeta/BetaAlpha occ, vir, orbs, cabs, and occ orbitals
//    // for alpha case: v_orbs1_ab has the orbital in alpha space
//    // for beta case: v_orbs1_ab has the orbital in beta space
//    vector< Ref<OrbitalSpace> > v_orbs1_ab;
//    vector< Ref<OrbitalSpace> > v_orbs2_ab;
//    if (spin == Alpha) {
//      obtain_orbitals(AlphaBeta, v_orbs1_ab, v_orbs2_ab);
//    } else {
//        obtain_orbitals(AlphaBeta, v_orbs2_ab, v_orbs1_ab);
//    }
//    // get the numbers of the occupied alpha and beta orbitals
//    occ1_act = v_orbs1_ab[0];
//    occ2_act = v_orbs2_ab[0];
//    vir1 = v_orbs1_ab[1];
//    vir2 = v_orbs2_ab[1];
//    nocc1_act = occ1_act->rank();
//    nocc2_act = occ2_act->rank();
//    nvir1 = vir1->rank();
//    nvir2 = vir2->rank();
//
//    nocc12 = nocc1_act * nocc2_act;
//    const int nvir12 = nvir1 * nvir2;
//    nocc12vir = nvir2 * nocc12;
//    T2ab_ij_ab = new double[nvir12 * nocc12];
//    fill_n(T2ab_ij_ab, nvir12 * nocc12, 0.0);
//    compute_T2_mp2(v_orbs1_ab, v_orbs2_ab, T2ab_ij_ab);
//  }
//
//  RefSCDimension rowdim = new SCDimension(norbs);
//  RefSCDimension coldim = new SCDimension(norbs);
//  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
//
//  RefSCMatrix D_mp2 = localkit->matrix(rowdim, coldim);
//  D_mp2.assign(0.0);
//
//  // D^i_j
//  // same spin part: T^ik_ab T^ab_jk - T^ki_ab T^ab_jk
//  // different spin part: T^ik_ab T^ab_jk
//  for (int i = 0; i < nocc_act; ++i) {
//    for (int j = 0; j <= i; ++j) {
//      // AlphaAlpha or BetaBeta part
//      double Dij_1 = 0;
//      double Dij_2 = 0;
//
//      for (int k = 0; k < nocc_act; ++k) {
//        const double* iter_T2_ikab = T2ab_ij + i * nocc_act + k;
//        const double* iter_T2_kiab = T2ab_ij + k * nocc_act + i;
//        const double* iter_T2_jkab = T2ab_ij + j * nocc_act + k;
//
//        for (int a = 0; a < nvir; ++a) {
//          for (int b = 0; b < nvir; ++b) {
//            Dij_1 += (*iter_T2_ikab) * (*iter_T2_jkab);
//            Dij_2 += (*iter_T2_kiab) * (*iter_T2_jkab);
//
//            iter_T2_ikab += noccocc;
//            iter_T2_kiab += noccocc;
//            iter_T2_jkab += noccocc;
//          }
//        }
//      }
//
//      // AlphaBeta part
//      double Dij_3 = 0;
//      if (nspincase2 == 3) {
//        for (int k = 0; k < nocc2_act; ++k) {
//          const double* iter_T2_ikab = T2ab_ij_ab + i * nocc_act + k;
//          const double* iter_T2_jkab = T2ab_ij_ab + j * nocc_act + k;
//
//          for (int a = 0; a < nvir1; ++a) {
//            for (int b = 0; b < nvir2; ++b) {
//              Dij_3 += (*iter_T2_ikab) * (*iter_T2_jkab);
//
//              iter_T2_ikab += noccocc;
//              iter_T2_jkab += noccocc;
//            }
//          }
//        }
//      } else {
//          Dij_3 = Dij_1;
//      }
//
//      const double Dmp2_ij = -(Dij_1 - Dij_2 + Dij_3);
//      D_mp2.set_element(i, j, Dmp2_ij);
//    }
//  }
//
//  // D^a_b
//  // same spin part: T^ac_ij T^ij_bc - T^ca_ij T^ij_bc
//  // different spin part: T^ac_ij T^ij_bc
//  const blasint one = 1; // for F77_DDOT
//  for (int a = 0; a < nvir; ++a) {
//
//    const double* iter_T2_bcij = T2ab_ij;
//    for (int b = 0; b <= a; ++b) {
//     // AlphaAlpha or BetaBeta part
//      double Dab_1 = 0;
//      double Dab_2 = 0;
//
//      const double* iter_T2_acij = T2ab_ij + a * noccoccvir;
//      const double* iter_T2_caij = T2ab_ij + a * noccocc;
//      for (int c = 0; c < nvir; ++c) {
//        const double TTac_bc = F77_DDOT(&noccocc, iter_T2_acij, &one, iter_T2_bcij, &one);
//        const double TTca_bc = F77_DDOT(&noccocc, iter_T2_caij, &one, iter_T2_bcij, &one);
//        Dab_1 += TTac_bc;
//        Dab_2 += TTca_bc;
//
//        iter_T2_acij += noccocc;
//        iter_T2_bcij += noccocc;
//        iter_T2_caij += noccoccvir;
//      }
//
//      // alpha beta part
//      double Dab_3 = 0;
//      if (nspincase2 == 3) {
//        const double* iter_T2_acij_ab = T2ab_ij_ab + a * nocc12vir;
//        const double* iter_T2_bcij_ab = T2ab_ij_ab + b * nocc12vir;
//        for (int c = 0; c < nvir2; ++c) {
//          const double TTac_bc = F77_DDOT(&nocc12, iter_T2_acij_ab, &one, iter_T2_bcij_ab, &one);
//          Dab_3 += TTac_bc;
//
//          iter_T2_acij_ab += nocc12;
//          iter_T2_bcij_ab += nocc12;
//        }
//      } else {
//          Dab_3 = Dab_1;
//      }
//
//      const double Dmp2_ab = Dab_1 - Dab_2 + Dab_3;
//      D_mp2.set_element(a+nocc_act, b+nocc_act, Dmp2_ab);
//    }
//  }
//  delete[] T2ab_ij;
//  if (nspincase2 == 3) delete[] T2ab_ij_ab;
//
//  // symmetrize matrix
//  for (int p = 0; p < norbs; ++p) {
//    for (int q = 0; q < p; ++q) {
//      D_mp2.set_element(q, p, D_mp2.get_element(p,q));
//    }
//  }
//
//  return D_mp2;
//}
//// end of function: compute_1rdm_mp2_f12corr

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
      ExEnv::out0() << endl << indent << spinletters << " D^i_i (test):" << endl
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

//    const double* iter_Dmi = (spin == Alpha? Dm_i_alpha : Dm_i_beta);
//    print_intermediate(spinletters, "D^m_i:", iter_Dmi, nocc_act, nocc_act);

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
     if (r12eval()->dim_oo(spincase).n() != 0) {
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

void MP2R12Energy_Diag::compute_RT1_api(const int nspincases1, const int nspincases2,
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

  const blasint one = 1;

  // Obtain T1 amplitudes
  RefSCMatrix T1[NSpinCases1];
  for (int s = 0; s < nspincases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);
    T1[spin] = r12intermediates_->get_T1_cc(spin);

  }
  if (nspincases1 == 1) {
    T1[Beta] = T1[Alpha];
  }

  const Ref<OrbitalSpace> occ1_act = v_orbs1_ab[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2_ab[0];
  const Ref<OrbitalSpace> vir1 = v_orbs1_ab[1];
  const Ref<OrbitalSpace> vir2 = v_orbs2_ab[1];
  const Ref<OrbitalSpace> cabs1 = v_orbs1_ab[3];
  const Ref<OrbitalSpace> cabs2 = v_orbs2_ab[3];

  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();
  const blasint nvir1 = vir1->rank();
  const blasint nvir2 = vir2->rank();

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
     const blasint nvir = vir->rank();
     const int ncabs = cabs->rank();

     // test propose
//       ExEnv::out0() << endl << spinletters << " number of occupied orbital: " << nocc_act << endl
//                           << "number of virtual orbital: " << nvir << endl
//                           << "number of cabs: " << ncabs << endl;

     // calculate AlphaAlpha/BetaBeta part: R^a'a_ij T^j_a
     //                                   & R^a'a_ji T^j_a

     const int ncabs_occ= ncabs * nocc_act;
     double* const RT1= new double[ncabs_occ];   // R^a'a_ij T^j_a
     fill_n(RT1, ncabs_occ, 0.0);

     double* const RT2 = new double[ncabs_occ];  // R^a'a_ji T^j_a
     fill_n(RT2, ncabs_occ, 0.0);

     // activate integrals
     Ref<DistArray4> Rii_apa_ints;
     activate_ints(occ_act->id(), occ_act->id(), cabs->id(), vir->id(),
                   descr_f12_key, moints4_rtime, Rii_apa_ints);

     double* iter_RT1 = RT1;
     double* iter_RT2 = RT2;
     double* const raw_tja = new double[nvir];

     for (int ap = 0; ap < ncabs; ++ap) {
       for (int i = 0; i < nocc_act; ++i) {

         for (int j = 0; j < nocc_act; ++j) {
           const double* const f12_ijap_blk = Rii_apa_ints->retrieve_pair_block(i, j, f12_idx) + ap * nvir;
           const double* const f12_jiap_blk = Rii_apa_ints->retrieve_pair_block(j, i, f12_idx) + ap * nvir;

           RefSCVector tja = T1[spin].get_row(j);
           fill_n(raw_tja, nvir, 0.0);
           tja.convert(raw_tja);
//           double* iter_tja = raw_tja;
//           for (int a = 0; a < nvir; ++a) {
//             *iter_tja = tja.get_element(a);
//             ++iter_tja;
//           }

           // R^a'a_ij T^j_a & R^a'a_ji T^j_a sum over a
           const double f12t1_1 = F77_DDOT(&nvir, f12_ijap_blk, &one, raw_tja, &one);
           const double f12t1_2 = F77_DDOT(&nvir, f12_jiap_blk, &one, raw_tja, &one);

           *iter_RT1 += f12t1_1;
           *iter_RT2 += f12t1_2;
           Rii_apa_ints->release_pair_block(i, j, f12_idx);
           Rii_apa_ints->release_pair_block(j, i, f12_idx);
         }

         ++iter_RT1;
         ++iter_RT2;
       }
     }
     delete[] raw_tja;
     Rii_apa_ints->deactivate();

     if (debug_ >= DefaultPrintThresholds::mostN2) {
       print_intermediate(spinletters, "R^a'a_ij T^j_a", RT1, ncabs, nocc_act);
       print_intermediate(spinletters, "R^a'a_ji T^j_a", RT2, ncabs, nocc_act);
     }

     // calculate AlphaBeta part:
     // RT^a'a_ij T^j_a & RT^a'a_ji T^j_a
     double* RT1_ab = NULL;
     double* RT2_ab = NULL;

     if (nspincases2 == 3) {

       RT1_ab = new double[ncabs_occ];
       RT2_ab = new double[ncabs_occ];
       fill_n(RT1_ab, ncabs_occ, 0.0);
       fill_n(RT2_ab, ncabs_occ, 0.0);

       if (spin == Alpha) {
         // activate integrals
         Ref<DistArray4> Rii_12_ints;  // R_i(1)j(2)^a'(1)a(2)
         activate_ints(occ1_act->id(), occ2_act->id(), cabs1->id(), vir2->id(),
                       descr_f12_key, moints4_rtime, Rii_12_ints);
         Ref<DistArray4> Rii_21_ints;  // R_j(2)i(1)^a'(1)a(2)
         activate_ints(occ2_act->id(), occ1_act->id(), cabs1->id(), vir2->id(),
                       descr_f12_key, moints4_rtime, Rii_21_ints);

         // R_i(1)j(2)^a'(1)a(2) T_a(2)^j(2) &
         // R_j(2)i(1)^a'(1)a(2) T_a(2)^j(2)
         double* iter_RT1_ab = RT1_ab;
         double* iter_RT2_ab = RT2_ab;

         double* const raw_tja = new double[nvir2];

         for (int ap = 0; ap < ncabs; ++ap) {
           for (int i = 0; i < nocc_act; ++i) {

             for (int j = 0; j < nocc2_act; ++j) {
               const double* const f12_ijap_blk = Rii_12_ints->retrieve_pair_block(i, j, f12_idx) + ap * nvir2;
               const double* const f12_jiap_blk = Rii_21_ints->retrieve_pair_block(j, i, f12_idx) + ap * nvir2;

               RefSCVector tja = T1[Beta].get_row(j);
               fill_n(raw_tja, nvir2, 0.0);
               tja.convert(raw_tja);

               // sum over a(2)
               const double f12t1_1 = F77_DDOT(&nvir2, f12_ijap_blk, &one, raw_tja, &one);
               const double f12t1_2 = F77_DDOT(&nvir2, f12_jiap_blk, &one, raw_tja, &one);

               *iter_RT1_ab += f12t1_1;
               *iter_RT2_ab += f12t1_2;

               Rii_12_ints->release_pair_block(i, j, f12_idx);
               Rii_21_ints->release_pair_block(j, i, f12_idx);
             }

             ++iter_RT1_ab;
             ++iter_RT2_ab;
           }
         }
         delete[] raw_tja;
         Rii_12_ints->deactivate();
         Rii_21_ints->deactivate();

       } else {
           // activate integrals
           Ref<DistArray4> Rii_21_ints;  // R_i(2)j(1)^a'(2)a(1)
           activate_ints(occ2_act->id(), occ1_act->id(), cabs2->id(), vir1->id(),
                         descr_f12_key, moints4_rtime, Rii_21_ints);
           Ref<DistArray4> Rii_12_ints;  // R_j(1)i(2)^a'(2)a(1)
           activate_ints(occ1_act->id(), occ2_act->id(), cabs2->id(), vir1->id(),
                         descr_f12_key, moints4_rtime, Rii_12_ints);

           // R_i(2)j(1)^a'(2)a(1) T_a(1)^j(1) &
           // R_j(1)i(2)^a'(2)a(1) T_a(1)^j(1)
           double* iter_RT1_ab = RT1_ab;
           double* iter_RT2_ab = RT2_ab;

           double* const raw_tja = new double[nvir1];

           for (int ap = 0; ap < ncabs; ++ap) {
             for (int i = 0; i < nocc_act; ++i) {

               for (int j = 0; j < nocc1_act; ++j) {
                 const double* const f12_ijap_blk = Rii_21_ints->retrieve_pair_block(i, j, f12_idx) + ap * nvir1;
                 const double* const f12_jiap_blk = Rii_12_ints->retrieve_pair_block(j, i, f12_idx) + ap * nvir1;

                 RefSCVector tja = T1[Alpha].get_row(j);
                 fill_n(raw_tja, nvir1, 0.0);
                 tja.convert(raw_tja);

                 // sum over a(1)
                 const double f12t1_1 = F77_DDOT(&nvir1, f12_ijap_blk, &one, raw_tja, &one);
                 const double f12t1_2 = F77_DDOT(&nvir1, f12_jiap_blk, &one, raw_tja, &one);;

                 *iter_RT1_ab += f12t1_1;
                 *iter_RT2_ab += f12t1_2;

                 Rii_21_ints->release_pair_block(i, j, f12_idx);
                 Rii_12_ints->release_pair_block(j, i, f12_idx);
               }

               ++iter_RT1_ab;
               ++iter_RT2_ab;
             }
           }
           delete[] raw_tja;
           Rii_12_ints->deactivate();
           Rii_21_ints->deactivate();
       } // end of Beta

       if (debug_ >= DefaultPrintThresholds::mostN2) {
           print_intermediate(spinletters, "R^a'a_ij T^j_a", RT1_ab, ncabs, nocc_act);
           print_intermediate(spinletters, "R^a'a_ji T^j_a", RT2_ab, ncabs, nocc_act);
       }
     } else {
         RT1_ab = RT1;
         RT2_ab = RT2;
     }

     // calculate D^a'_i
     const double* iter_RT1_sp = RT1;
     const double* iter_RT2_sp = RT2;
     const double* iter_RT1_ab = RT1_ab;
     const double* iter_RT2_ab = RT2_ab;

     double* iter_D = (spin == Alpha? D_alpha : D_beta);

     for (int ap = 0;  ap < ncabs; ++ap) {
       for (int i = 0; i < nocc_act; ++i) {

         // AlphaAlpha/BetaBeta part
         double d_12 = C_1 * (*iter_RT1_sp - *iter_RT2_sp);
//         ExEnv::out0() << spinletters << " part d^" << idx1 << "_" << idx2 << " = "
//                       << scprintf("%12.10f", d_12) << endl;
         ++iter_RT1_sp;
         ++iter_RT2_sp;

         // AlphaBeta part
         if (nocc1_act != 0 && nocc2_act != 0) {
   //        ExEnv::out0() << "RR1: "  << scprintf("%12.10f", RR1)
   //                      << "  RR2: " << scprintf("%12.10f", RR2)
   //                      << "  RR3: "  << scprintf("%12.10f", RR3)
   //                      << "  RR4: " << scprintf("%12.10f", RR4)
   //                      << endl;

           d_12 += 0.5 * (C_0 + C_1) * (*iter_RT1_ab) + 0.5 * (C_0 - C_1) * (*iter_RT2_ab);
//           ExEnv::out0() << "AlphaBeta part: d^"  << idx1 << "_" << idx2 << " = "
//                         << scprintf("%12.10f", d_12) << endl;

           ++iter_RT1_ab;
           ++iter_RT2_ab;
         } // end of AlphaBeta part

         *iter_D = d_12;
         ++iter_D;
       } // end of looping over i
     } // end of calculating D^a'_i[s]

     delete[] RT1;
     delete[] RT2;
     if (nspincases2 == 3) {
       delete[] RT1_ab;
       delete[] RT2_ab;
     }
  } // end of spincase1 loop

}
// end of compute_RT1_api function

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

//        double* T2_test = new double[nocc1_act * nocc2_act];
//        double* iter_T2_test = T2_test;
//        fill_n(T2_test, nocc1_act * nocc2_act, 0.0);
//        ExEnv::out0() << endl << a << " " << b << endl;

        for (int i = 0; i < nocc1_act; ++i) {
          for (int j = 0; j < nocc2_act; ++j) {
            *iter_T2 =  (*gab_ij)
                       / (evals_i1(i) + evals_i2(j) - evals_a1(a) - evals_a2(b));

//            *iter_T2_test = *iter_T2;
//            ++iter_T2_test;

            ++gab_ij;
            ++iter_T2;
          }
        }
//        print_intermediate("AlphAlpha", "T2 test in compute_T2_mp2", T2_test, nocc1_act, nocc2_act);
//        delete[] T2_test;

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

// compute MP2 T2 amplitude (antisymmetrized), stored in array (a,b,i,j)
//  T^i(spin1)j(spin2)_a(spin1)b(spin2) 4-dimension matrix
//= g^ij_ab / (e_i + e_j - e_a - e_b)
void MP2R12Energy_Diag::compute_T2_mp2(const SpinCase2 spincase,
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
  const SpinCase1 spin1 = case1(spincase);
  const SpinCase1 spin2 = case2(spincase);

  const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
  const Ref<OrbitalSpace>& vir1 = r12eval()->vir(spin1);
  const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
  const Ref<OrbitalSpace>& vir2 = r12eval()->vir(spin2);

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
  if (spincase == AlphaBeta) {
    for (int a = 0; a < nvir1; ++a) {
      for (int b = 0; b < nvir2; ++b) {
        const double* gab_ij = a1a2i1i2_ints->retrieve_pair_block(a, b, eri_idx);

//        double* T2_test = new double[nocc1_act * nocc2_act];
//        double* iter_T2_test = T2_test;
//        fill_n(T2_test, nocc1_act * nocc2_act, 0.0);
//        ExEnv::out0() << endl << a << " " << b << endl;

        for (int i = 0; i < nocc1_act; ++i) {
          for (int j = 0; j < nocc2_act; ++j) {
            *iter_T2 =  (*gab_ij)
                       / (evals_i1(i) + evals_i2(j) - evals_a1(a) - evals_a2(b));

//            *iter_T2_test = *iter_T2;
//            ++iter_T2_test;

            ++gab_ij;
            ++iter_T2;
          }
        }
//        print_intermediate("AlphAlpha", "T2 test in compute_T2_mp2", T2_test, nocc1_act, nocc2_act);
//        delete[] T2_test;

        a1a2i1i2_ints->release_pair_block(a, b, eri_idx);
      }
    }
  } else {
      // AlphaAlpha or BetaBeta case
      for (int a = 0; a < nvir1; ++a) {
        for (int b = 0; b < nvir2; ++b) {
          const double* const gab_ij = a1a2i1i2_ints->retrieve_pair_block(a, b, eri_idx);

          const double* iter_gab_ij = gab_ij;
          for (int i = 0; i < nocc1_act; ++i) {
            const double* iter_gab_ji = gab_ij + i;

            for (int j = 0; j < nocc2_act; ++j) {
              *iter_T2 = (*iter_gab_ij - *iter_gab_ji)
                        / (evals_i1(i) + evals_i2(j) - evals_a1(a) - evals_a2(b));

               ++iter_gab_ij;
               iter_gab_ji += nocc1_act;
               ++iter_T2;
             }
           }
           a1a2i1i2_ints->release_pair_block(a, b, eri_idx);
         }
       }
    }

    a1a2i1i2_ints->deactivate();
}
// end of function: compute_T2_mp2

// compute F12 corrected T2 amplitude (antisymmetrized), stored in array (a,b,i,j)
// \tilde{T}^ij_ab =  - \bar{C}^ij_ab / (e_i + e_j - e_a - e_b)
// C^ij_ab = R^ij_a'b F^a'_a + R^ij_aa' F^a'_b
void MP2R12Energy_Diag::compute_T2abij_f12corr(const SpinCase2 spincase,
                                               const double C_0, const double C_1,
                                               double* const T2ab_ij_f12corr)
{
  const int nspincase2 = (r12eval()->spin_polarized() ? 3 : 2);

  // get moints4_rtime, descr_f12_key, and eri_idx
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type =
      r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);
  const TwoBodyOper::type eri_type =
      r12world->r12tech()->corrfactor()->tbint_type_eri();
  const unsigned int eri_idx = descr_f12->intset(eri_type);

  // obtain orbitals
  const SpinCase1 spin1 = case1(spincase);
  const SpinCase1 spin2 = case2(spincase);

  const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
  const Ref<OrbitalSpace>& vir1 = r12eval()->vir(spin1);
  const Ref<OrbitalSpace>& fvir1 = r12eval()->F_a_A(spin1);

  const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
  const Ref<OrbitalSpace>& vir2 = r12eval()->vir(spin2);
  const Ref<OrbitalSpace>& fvir2 = r12eval()->F_a_A(spin2);

  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();
  const int nvir1 = vir1->rank();
  const int nvir2 = vir2->rank();

  // R^ij_a'b f^a'_a
  Ref<DistArray4> AF1a2i1i2_ints = NULL;
  activate_ints(fvir1->id(), vir2->id(), occ1_act->id(),
                occ2_act->id(), descr_f12_key, moints4_rtime,
                AF1a2i1i2_ints);
  // R^ij_aa' f^a'_b
  Ref<DistArray4> a1AF2i1i2_ints = NULL;
  activate_ints(vir1->id(), fvir2->id(), occ1_act->id(),
                occ2_act->id(), descr_f12_key, moints4_rtime,
                a1AF2i1i2_ints);

  // get eigenvalues of Fock matrix
  const RefDiagSCMatrix evals_i1 = occ1_act->evals();
  const RefDiagSCMatrix evals_a1 = vir1->evals();

  double* iter_T2 = T2ab_ij_f12corr;
  if (spincase == AlphaBeta){

    if (nspincase2 == 3) {
      Ref<DistArray4> AF1a2i2i1_ints = NULL;
      Ref<DistArray4> a1AF2i2i1_ints = NULL;
      activate_ints(fvir1->id(), vir2->id(), occ2_act->id(),
                    occ1_act->id(), descr_f12_key, moints4_rtime,
                    AF1a2i2i1_ints);
      activate_ints(vir1->id(), fvir2->id(), occ2_act->id(),
                    occ1_act->id(), descr_f12_key, moints4_rtime,
                    a1AF2i2i1_ints);

      const RefDiagSCMatrix evals_i2 = occ2_act->evals();
      const RefDiagSCMatrix evals_a2 = vir2->evals();

      for (int a = 0; a < nvir1; ++a) {
        const double Faa = evals_a1(a); // F^a_a

        for (int b = 0; b < nvir2; ++b) {
          const double Fbb = evals_a2(b);        // F^b_b

          const double* const AF1a2_i1i2 = AF1a2i1i2_ints->retrieve_pair_block(a, b, f12_idx);
          const double* const AF1a2_i2i1 = AF1a2i2i1_ints->retrieve_pair_block(a, b, f12_idx);

          const double* const a1AF2_i1i2 = a1AF2i1i2_ints->retrieve_pair_block(a, b, f12_idx);
          const double* const a1AF2_i2i1 = a1AF2i2i1_ints->retrieve_pair_block(a, b, f12_idx);

          const double* R_afb_ij = AF1a2_i1i2;   // F^a'1_a1 R^i1j2_a'1b2
          const double* R_abf_ij = a1AF2_i1i2;   // F^a'2_b2 R^i1j2_a1a'2

          for (int i = 0; i < nocc1_act; ++i) {
            const double Fii = evals_i1(i);            // F^i_i
            const double* R_afb_ji = AF1a2_i2i1 + i;   // F^a'1_a1 R^j2i1_a'1b2
            const double* R_abf_ji = a1AF2_i2i1 + i;   // F^a'2_b2 R^j2i1_a1a'2

            for (int j = 0; j < nocc2_act; ++j) {
              const double Fjj = evals_i2(j);  // F^j_j

              const double denom = 1.0 / (Fii + Fjj - Faa - Fbb);
              *iter_T2 =  (0.5*(C_0 + C_1) * (*R_afb_ij + *R_abf_ij)
                          + 0.5*(C_0 - C_1) * (*R_afb_ji + *R_abf_ji)
                          ) * denom;

              ++iter_T2;
              ++R_afb_ij;
              ++R_abf_ij;
              R_afb_ji += nocc1_act;
              R_abf_ji += nocc1_act;
            }
          }
          AF1a2i1i2_ints->release_pair_block(a, b, f12_idx);
          AF1a2i2i1_ints->release_pair_block(a, b, f12_idx);
          a1AF2i1i2_ints->release_pair_block(a, b, f12_idx);
          a1AF2i2i1_ints->release_pair_block(a, b, f12_idx);
        }
      }

      AF1a2i2i1_ints->deactivate();
      a1AF2i2i1_ints->deactivate();
    } else {
        for (int a = 0; a < nvir1; ++a) {
          const double Faa = evals_a1(a); // F^a_a

          for (int b = 0; b < nvir2; ++b) {
            const double Fbb = evals_a1(b);        // F^b_b

            const double* const AFaii = AF1a2i1i2_ints->retrieve_pair_block(a, b, f12_idx);
            const double* const aAFii = a1AF2i1i2_ints->retrieve_pair_block(a, b, f12_idx);

            const double* R_afb_ij = AFaii;        // F^a'_a R^ij_a'b
            const double* R_abf_ij = aAFii;        // F^a'_b R^ij_aa'

            for (int i = 0; i < nocc1_act; ++i) {
              const double Fii = evals_i1(i);       // F^i_i
              const double* R_afb_ji = AFaii + i;   // F^a'_a R^ji_a'b
              const double* R_abf_ji = aAFii + i;   // F^a'_b R^ji_aa'

              for (int j = 0; j < nocc2_act; ++j) {
                const double Fjj = evals_i1(j);            // F^j_j
                const double denom = 1.0 / (Fii + Fjj - Faa - Fbb);

                *iter_T2 =  ( 0.5*(C_0 + C_1) * (*R_afb_ij + *R_abf_ij)
                            + 0.5*(C_0 - C_1) * (*R_afb_ji + *R_abf_ji)
                            ) * denom;

                ++iter_T2;
                ++R_afb_ij;
                ++R_abf_ij;
                R_afb_ji += nocc1_act;
                R_abf_ji += nocc1_act;

              }
            }
            AF1a2i1i2_ints->release_pair_block(a, b, f12_idx);
            a1AF2i1i2_ints->release_pair_block(a, b, f12_idx);
          }
        }

    }
    // end of AlphaBeta case
  } else {
      // AlphaAlpha or BetaBeta case
      for (int a = 0; a < nvir1; ++a) {
        const double Faa = evals_a1(a); // F^a_a

        for (int b = 0; b < nvir2; ++b) {
          const double Fbb = evals_a1(b);        // F^b_b

          const double* const AFaii = AF1a2i1i2_ints->retrieve_pair_block(a, b, f12_idx);
          const double* const aAFii = a1AF2i1i2_ints->retrieve_pair_block(a, b, f12_idx);

          const double* R_afb_ij = AFaii;        // F^a'_a R^ij_a'b
          const double* R_abf_ij = aAFii;        // F^a'_b R^ij_aa'

          for (int i = 0; i < nocc1_act; ++i) {
            const double Fii = evals_i1(i);       // F^i_i

            const double* R_afb_ji = AFaii + i;   // F^a'_a R^ji_a'b
            const double* R_abf_ji = aAFii + i;   // F^a'_b R^ji_aa'

            for (int j = 0; j < nocc2_act; ++j) {
              const double Fjj = evals_i1(j);       // F^j_j

              const double denom = 1.0 / (Fii+ Fjj - Faa - Fbb);
              *iter_T2 =  ( C_1 * (*R_afb_ij - *R_afb_ji
                                 + *R_abf_ij - *R_abf_ji)
                          ) * denom;

              ++iter_T2;
              ++R_afb_ij;
              ++R_abf_ij;
              R_afb_ji += nocc1_act;
              R_abf_ji += nocc1_act;
            }
          }
          AF1a2i1i2_ints->release_pair_block(a, b, f12_idx);
          a1AF2i1i2_ints->release_pair_block(a, b, f12_idx);
        }
      }
  }
  a1AF2i1i2_ints->deactivate();
}
// end of compute_T2abij_f12corr

// compute MP2F12 T2 amplitude = T2 (MP2) + T2 (F12 corrected)
void MP2R12Energy_Diag::compute_T2abij_mp2f12(const int nocc1_act, const int nocc2_act,
                                              const int nvir1, const int nvir2,
                                              const double* const T2ab_ij_mp2,
                                              const double* const T2ab_ij_f12corr,
                                              double* const T2abij_mp2f12)
{
  const double* iter_T2_mp2 = T2ab_ij_mp2;
  const double* iter_T2_f12corr = T2ab_ij_f12corr;
  double* iter_T2_mp2f12 = T2abij_mp2f12;

  for (int a = 0; a < nvir1; ++a) {
    for (int b = 0; b < nvir2; ++b) {
      for (int i = 0; i < nocc1_act; ++i) {
        for (int j = 0; j < nocc2_act; ++j) {
          *iter_T2_mp2f12 = *iter_T2_mp2 + *iter_T2_f12corr;

          ++iter_T2_mp2;
          ++iter_T2_f12corr;
          ++iter_T2_mp2f12;
        }
      }
    }
  }

}
// end of compute_T2abij_mp2f12

// compute MP2F12 T2 amplitude (antisymmetrized), stored in array (a,b,i,j)
// \tilde{T}^ij_ab = \bar{T}^ij_ab + \bar{C}^ij_ab / (e_i + e_j - e_a - e_b)
// C^ij_ab = R^ij_a'b F^a'_a + R^ij_aa' F^a'_b
void MP2R12Energy_Diag::compute_T2abij_mp2f12(const SpinCase2 spincase,
                                              const double C_0, const double C_1,
                                              double* const T2ab_ij_mp2f12)
{
  const int nspincase2 = (r12eval()->spin_polarized() ? 3 : 2);

  // get moints4_rtime, descr_f12_key, and eri_idx
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type =
      r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);
  const TwoBodyOper::type eri_type =
      r12world->r12tech()->corrfactor()->tbint_type_eri();
  const unsigned int eri_idx = descr_f12->intset(eri_type);

  // obtain orbitals
  const SpinCase1 spin1 = case1(spincase);
  const SpinCase1 spin2 = case2(spincase);

  const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
  const Ref<OrbitalSpace>& vir1 = r12eval()->vir(spin1);
  const Ref<OrbitalSpace>& fvir1 = r12eval()->F_a_A(spin1);

  const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
  const Ref<OrbitalSpace>& vir2 = r12eval()->vir(spin2);
  const Ref<OrbitalSpace>& fvir2 = r12eval()->F_a_A(spin2);

  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();
  const int nvir1 = vir1->rank();
  const int nvir2 = vir2->rank();

  // g^ab_ij
  Ref<DistArray4> a1a2i1i2_ints = NULL;
  activate_ints(vir1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                descr_f12_key, moints4_rtime,
                a1a2i1i2_ints);

  // R^ij_a'b f^a'_a
  Ref<DistArray4> AF1a2i1i2_ints = NULL;
  activate_ints(fvir1->id(), vir2->id(), occ1_act->id(),
                occ2_act->id(), descr_f12_key, moints4_rtime,
                AF1a2i1i2_ints);
  // R^ij_aa' f^a'_b
  Ref<DistArray4> a1AF2i1i2_ints = NULL;
  activate_ints(vir1->id(), fvir2->id(), occ1_act->id(),
                occ2_act->id(), descr_f12_key, moints4_rtime,
                a1AF2i1i2_ints);

  // get eigenvalues of Fock matrix
  const RefDiagSCMatrix evals_i1 = occ1_act->evals();
  const RefDiagSCMatrix evals_a1 = vir1->evals();

  double* iter_T2 = T2ab_ij_mp2f12;
  if (spincase == AlphaBeta){

    if (nspincase2 == 3) {
      Ref<DistArray4> AF1a2i2i1_ints = NULL;
      Ref<DistArray4> a1AF2i2i1_ints = NULL;
      activate_ints(fvir1->id(), vir2->id(), occ2_act->id(),
                    occ1_act->id(), descr_f12_key, moints4_rtime,
                    AF1a2i2i1_ints);
      activate_ints(vir1->id(), fvir2->id(), occ2_act->id(),
                    occ1_act->id(), descr_f12_key, moints4_rtime,
                    a1AF2i2i1_ints);

      const RefDiagSCMatrix evals_i2 = occ2_act->evals();
      const RefDiagSCMatrix evals_a2 = vir2->evals();

      for (int a = 0; a < nvir1; ++a) {
        const double Faa = evals_a1(a); // F^a_a

        for (int b = 0; b < nvir2; ++b) {
          const double Fbb = evals_a2(b);        // F^b_b

          const double* gab_ij = a1a2i1i2_ints->retrieve_pair_block(a, b, eri_idx);

          const double* const AF1a2_i1i2 = AF1a2i1i2_ints->retrieve_pair_block(a, b, f12_idx);
          const double* const AF1a2_i2i1 = AF1a2i2i1_ints->retrieve_pair_block(a, b, f12_idx);

          const double* const a1AF2_i1i2 = a1AF2i1i2_ints->retrieve_pair_block(a, b, f12_idx);
          const double* const a1AF2_i2i1 = a1AF2i2i1_ints->retrieve_pair_block(a, b, f12_idx);

          const double* R_afb_ij = AF1a2_i1i2;   // F^a'1_a1 R^i1j2_a'1b2
          const double* R_abf_ij = a1AF2_i1i2;   // F^a'2_b2 R^i1j2_a1a'2

          for (int i = 0; i < nocc1_act; ++i) {
            const double Fii = evals_i1(i);            // F^i_i
            const double* R_afb_ji = AF1a2_i2i1 + i;   // F^a'1_a1 R^j2i1_a'1b2
            const double* R_abf_ji = a1AF2_i2i1 + i;   // F^a'2_b2 R^j2i1_a1a'2

            for (int j = 0; j < nocc2_act; ++j) {
              const double Fjj = evals_i2(j);  // F^j_j

              const double denom = 1.0 / (Fii + Fjj - Faa - Fbb);
              *iter_T2 =  ((*gab_ij)
                          + 0.5*(C_0 + C_1) * (*R_afb_ij + *R_abf_ij)
                          + 0.5*(C_0 - C_1) * (*R_afb_ji + *R_abf_ji)
                          ) * denom;

              ++iter_T2;
              ++gab_ij;
              ++R_afb_ij;
              ++R_abf_ij;
              R_afb_ji += nocc1_act;
              R_abf_ji += nocc1_act;
            }
          }
          a1a2i1i2_ints->release_pair_block(a, b, eri_idx);
          AF1a2i1i2_ints->release_pair_block(a, b, f12_idx);
          AF1a2i2i1_ints->release_pair_block(a, b, f12_idx);
          a1AF2i1i2_ints->release_pair_block(a, b, f12_idx);
          a1AF2i2i1_ints->release_pair_block(a, b, f12_idx);
        }
      }

      AF1a2i2i1_ints->deactivate();
      a1AF2i2i1_ints->deactivate();
    } else {
        for (int a = 0; a < nvir1; ++a) {
          const double Faa = evals_a1(a); // F^a_a

          for (int b = 0; b < nvir2; ++b) {
            const double Fbb = evals_a1(b);        // F^b_b

            const double* gab_ij = a1a2i1i2_ints->retrieve_pair_block(a, b, eri_idx);
            const double* const AFaii = AF1a2i1i2_ints->retrieve_pair_block(a, b, f12_idx);
            const double* const aAFii = a1AF2i1i2_ints->retrieve_pair_block(a, b, f12_idx);

            const double* R_afb_ij = AFaii;        // F^a'_a R^ij_a'b
            const double* R_abf_ij = aAFii;        // F^a'_b R^ij_aa'

            for (int i = 0; i < nocc1_act; ++i) {
              const double Fii = evals_i1(i);       // F^i_i
              const double* R_afb_ji = AFaii + i;   // F^a'_a R^ji_a'b
              const double* R_abf_ji = aAFii + i;   // F^a'_b R^ji_aa'

              for (int j = 0; j < nocc2_act; ++j) {
                const double Fjj = evals_i1(j);            // F^j_j
                const double denom = 1.0 / (Fii + Fjj - Faa - Fbb);

                *iter_T2 =  ((*gab_ij)
                            + 0.5*(C_0 + C_1) * (*R_afb_ij + *R_abf_ij)
                            + 0.5*(C_0 - C_1) * (*R_afb_ji + *R_abf_ji)
                            ) * denom;

                ++iter_T2;
                ++gab_ij;
                ++R_afb_ij;
                ++R_abf_ij;
                R_afb_ji += nocc1_act;
                R_abf_ji += nocc1_act;

              }
            }
            a1a2i1i2_ints->release_pair_block(a, b, eri_idx);
            AF1a2i1i2_ints->release_pair_block(a, b, f12_idx);
            a1AF2i1i2_ints->release_pair_block(a, b, f12_idx);
            // print MP2F12 T2^ab_ij
            #if 0
            string spinletters = to_string(spincase);
            ExEnv::out0() <<endl << a << " " << b << endl;
            print_intermediate(spinletters, "Tabij", Tabij, nocc1_act, nocc2_act);
            #endif
          }
        }

    }
    // end of AlphaBeta case
  } else {
      // AlphaAlpha or BetaBeta case
      for (int a = 0; a < nvir1; ++a) {
        const double Faa = evals_a1(a); // F^a_a

        for (int b = 0; b < nvir2; ++b) {
          const double Fbb = evals_a1(b);        // F^b_b

          const double* const aaii = a1a2i1i2_ints->retrieve_pair_block(a, b, eri_idx);
          const double* const AFaii = AF1a2i1i2_ints->retrieve_pair_block(a, b, f12_idx);
          const double* const aAFii = a1AF2i1i2_ints->retrieve_pair_block(a, b, f12_idx);

          const double* gab_ij = aaii;           // g^ij_ab
          const double* R_afb_ij = AFaii;        // F^a'_a R^ij_a'b
          const double* R_abf_ij = aAFii;        // F^a'_b R^ij_aa'

          for (int i = 0; i < nocc1_act; ++i) {
            const double Fii = evals_i1(i);       // F^i_i

            const double* gab_ji = aaii + i;     // g^ji_ab
            const double* R_afb_ji = AFaii + i;   // F^a'_a R^ji_a'b
            const double* R_abf_ji = aAFii + i;   // F^a'_b R^ji_aa'

            for (int j = 0; j < nocc2_act; ++j) {
              const double Fjj = evals_i1(j);       // F^j_j

              const double denom = 1.0 / (Fii+ Fjj - Faa - Fbb);
              *iter_T2 =  ((*gab_ij - *gab_ji)
                          + C_1 * (*R_afb_ij - *R_abf_ji
                                 + *R_abf_ij - *R_abf_ji)
                          ) * denom;

              ++iter_T2;
              ++gab_ij;
              ++R_afb_ij;
              ++R_abf_ij;
              gab_ji += nocc1_act;
              R_afb_ji += nocc1_act;
              R_abf_ji += nocc1_act;
            }
          }
          a1a2i1i2_ints->release_pair_block(a, b, eri_idx);
          AF1a2i1i2_ints->release_pair_block(a, b, f12_idx);
          a1AF2i1i2_ints->release_pair_block(a, b, f12_idx);
        }
      }
  }
  a1a2i1i2_ints->deactivate();
  a1AF2i1i2_ints->deactivate();
}
// end of compute_mp2f12_T2abij

// compute CABS canonical T1 amplitude
//RefSCMatrix MP2R12Energy_Diag::compute_CABS_T1(SpinCase1 spin)
//{
//  // obtain occ, vir, orbs, and cabs orbitals in that order
//  const SpinCase2 spincase = static_cast<SpinCase2>(spin+1);
//  string spinletters = to_string(spincase);
//  vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
//  vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
//  obtain_orbitals(spincase, v_orbs1, v_orbs2);
//
//  // obtain orbitals
//  const Ref<OrbitalSpace> occ_act = v_orbs1[0];
//  const Ref<OrbitalSpace> vir = v_orbs1[1];
//  const Ref<OrbitalSpace> cabs = v_orbs1[3];
//
//  const int nocc_act = occ_act->rank();
//  const int nvir = vir->rank();
//  const int ncabs = cabs->rank();
//  const int nvir_com = nvir + ncabs;
//
////  ExEnv::out0() << endl << "number of active occupied orbital: " << nocc_act << endl
////                        << "number of virtual orbital: " << nvir << endl
////                        << "number of cabs: " << ncabs << endl;
//
//  // get the Fock matrix element
//  RefSCMatrix Fia = r12eval()->fock(occ_act, vir, spin);
////  Fia.print(prepend_spincase(spin,"Fia").c_str());
//
//  RefSCMatrix Fiap = r12eval()->fock(occ_act, cabs, spin);
////  Fiap.print(prepend_spincase(spin,"Fia'").c_str());
//
////  const RefDiagSCMatrix& i_evals = occ_act->evals();
////  const RefDiagSCMatrix& a_evals = vir->evals();
////  const RefDiagSCMatrix& ap_evals = cabs->evals();
////  i_evals.print(prepend_spincase(spin,"i_evals").c_str());
////  a_evals.print(prepend_spincase(spin,"a_evals").c_str());
////  ap_evals.print(prepend_spincase(spin,"ap_evals").c_str());
//
//  RefSCMatrix Fii = r12eval()->fock(occ_act, occ_act, spin);
//  RefSCMatrix Faa = r12eval()->fock(vir, vir, spin);
//  RefSCMatrix Fapap = r12eval()->fock(cabs, cabs, spin);
////  Fii.print(prepend_spincase(spin,"Fij").c_str());
////  Faa.print(prepend_spincase(spin,"Fab").c_str());
////  Fapap.print(prepend_spincase(spin,"Fa'b'").c_str());
//
//  RefSCDimension rowdim = new SCDimension(nocc_act);
//  RefSCDimension coldim = new SCDimension(nvir_com);
//  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
//  RefSCMatrix T = localkit->matrix(rowdim, coldim);
//  T.assign(0.0);
//
//  for (int i = 0; i < nocc_act; ++i) {
//
//    for (int a = 0; a < nvir; ++a) {
//      const double tia = Fia.get_element(i,a)
////                        /(i_evals(i) - a_evals(a));
//                         /(Fii.get_element(i,i) - Faa.get_element(a,a));
//      T.set_element(i, a, tia);
//    }
//
//    for (int ap = 0; ap < ncabs; ++ap) {
//      const double tiap = Fiap.get_element(i,ap)
////                         /(i_evals(i) - ap_evals(ap));
//                          /(Fii.get_element(i,i) - Fapap.get_element(ap,ap));
//      T.set_element(i, ap+nvir, tiap);
//    }
//  }
//  //T.print(prepend_spincase(spin,"CABS amplitude Tia").c_str());
//
//  return T;
//}
//// end of function: compute_cabs_t1

// compute one-electron density form CABS contribution
RefSCMatrix MP2R12Energy_Diag::compute_D_CABS(SpinCase1 spin) {

  // obtain occ, vir, orbs, and cabs orbitals in that order
  const SpinCase2 spincase = static_cast<SpinCase2>(spin+1);
  string spinletters = to_string(spincase);
  vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
  vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
  obtain_orbitals(spincase, v_orbs1, v_orbs2);

  // obtain orbitals
  const Ref<OrbitalSpace> occ = v_orbs1[4];
  const Ref<OrbitalSpace> vir = v_orbs1[1];
  const Ref<OrbitalSpace> cabs = v_orbs1[3];

  const int nocc = occ->rank();
  const int nvir = vir->rank();
  const int ncabs = cabs->rank();
  const int norbs = nocc + nvir;
  const int nvir_com = nvir + ncabs;
  const int norbs_com = nocc + nvir + ncabs;

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCDimension rowdim = new SCDimension(norbs_com);
  RefSCDimension coldim = new SCDimension(norbs_com);
  RefSCMatrix D = localkit->matrix(rowdim, coldim);
  D.assign(0.0);

  const RefSCMatrix T1_cabs = r12eval_->T1_cabs(spin);
  //T1_cabs.print(prepend_spincase(spin,"CABS T1 amplitude").c_str());

//  // test for T1_cabs for different symmetry
//  {
//  RefSCDimension rowdim = new SCDimension(T1_cabs.nrow()+T1_cabs.ncol());
//  RefSCDimension coldim = rowdim;
//
//  RefDiagSCMatrix evals = localkit->diagmatrix(rowdim);
//  RefSCMatrix evecs = localkit->matrix(rowdim, coldim);
//  evals.assign(0.0);
//
//  RefSymmSCMatrix T1 = localkit->symmmatrix(rowdim);
//  T1.assign(0.0);
//  //T1.copyRefSCMatrix(T1_cabs);
//  for(int i = 0; i < nocc; i++) {
//    for(int ap = 0; ap < ncabs; ap++) {
//      T1.set_element(i,ap+nocc,T1_cabs.get_element(i,ap));
//    }
//  }
//  //T1.print(prepend_spincase(spin,"CABS T1 amplitude (RefSymmSCMatrix)").c_str());
//
//  T1.diagonalize(evals, evecs);
//  T1 = 0;
//  evals.print("evals of T1_cabs");
//
//  evals = 0;
//  evecs = 0;
//  }

   RefSCDimension rowdim_occ = new SCDimension(nocc);
   RefSCDimension coldim_vir_com = new SCDimension(nvir_com);
   RefSCMatrix T = localkit->matrix(rowdim_occ, coldim_vir_com);
   T.assign(0.0);

   if (r12intermediates_->T2_cc_computed()) {
     RefSCMatrix T1_ccsd = r12intermediates_->get_T1_cc(spin);
     //T1_ccsd.print(prepend_spincase(spin,"CCSD T1 amplitude").c_str());
     const Ref<OrbitalSpace> occ_act = v_orbs1[0];
     const int nocc_act = occ_act->rank();
     const int nfzc = nocc - nocc_act;

//     // test for T1_ccsd for different symmetry
//     {
//     RefSCDimension rowdim = new SCDimension(T1_ccsd.nrow()+T1_ccsd.ncol());
//     RefSCDimension coldim = rowdim;
//
//     RefDiagSCMatrix evals = localkit->diagmatrix(rowdim);
//     RefSCMatrix evecs = localkit->matrix(rowdim, coldim);
//     evals.assign(0.0);
//
//     RefSymmSCMatrix T1 = localkit->symmmatrix(rowdim);
//     T1.assign(0.0);
//     //T1.copyRefSCMatrix(T1_ccsd);
//     for(int i = 0; i < nocc_act; i++) {
//       for(int a = 0; a < nvir; a++) {
//         T1.set_element(i,a+nocc_act,T1_ccsd.get_element(i,a));
//       }
//     }
//     //T1.print(prepend_spincase(spin,"CCSD T1 amplitude (RefSymmSCMatrix)").c_str());
//
//     T1.diagonalize(evals, evecs);
//     T1 = 0;
//     evals.print("evals of T1_ccsd");
//
//     evals = 0;
//     evecs = 0;
//     }

     T.accumulate_subblock(T1_ccsd, nfzc, nocc-1, 0, nvir-1, 0, 0);
   }

   const int nirreps = occ->nblocks();
   const std::vector<unsigned int>& occpi = occ->block_sizes();
   const std::vector<unsigned int>& virpi = vir->block_sizes();
   const std::vector<unsigned int>& cabspi = cabs->block_sizes();

   std::vector<unsigned int> occoff(nirreps);
   std::vector<unsigned int> viroff(nirreps);
   std::vector<unsigned int> cabsoff(nirreps);
   std::vector<unsigned int> vir_comoff(nirreps);

   occoff[0] = 0;
   viroff[0] = 0;
   cabsoff[0] = 0;
//   ExEnv::out0() << endl << "  occpi vir cabspi:" << endl
//                         << 0 << "  " << occpi[0] << "  " << virpi[0] << "  " << cabspi[0] <<  endl;

   for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
     occoff[irrep] = occoff[irrep-1] + occpi[irrep-1];
     viroff[irrep] = viroff[irrep-1] + virpi[irrep-1];
     cabsoff[irrep] = cabsoff[irrep-1] + cabspi[irrep-1];
     vir_comoff[irrep] = viroff[irrep] + cabsoff[irrep];
     //ExEnv::out0() << irrep << "   " << occpi[irrep] << "  " << virpi[irrep] << "  " << cabspi[irrep] << endl;
   }

   for (int h = 0; h < nirreps; ++h) {

     if (r12intermediates_->T2_cc_computed()) {

//       for (int i = 0; i < occpi[h]; ++i) {
//         for (int ap = 0; ap < cabspi[h]; ++ap) {
//           const int idxi1 = occoff[h] + i;
//           const int idxap1 = cabsoff[h] + ap;
//           const int idxi2 = idxi1;
//           const int idxap2 = nvir + cabsoff[h] + ap;
//           T.set_element(idxi2,idxap2,T1_cabs.get_element(idxi1,idxap1));
//         }
//       }
       for (int i = 0; i < nocc; ++i) {
         for (int ap = 0; ap < ncabs; ++ap) {
           T.set_element(i,nvir+ap,T1_cabs.get_element(i,ap));
         }
       }

     } else {
         for (int i = 0; i < occpi[h]; ++i) {
           for (int a = 0; a < virpi[h]; ++a) {
             const int idxi1 = occoff[h] + i;
             const int idxa1 = vir_comoff[h] + a;
             const int idxi2 = idxi1;
             const int idxa2 = viroff[h] + a;
             T.set_element(idxi2,idxa2,T1_cabs.get_element(idxi1,idxa1));
           }

           for (int ap = 0; ap < cabspi[h]; ++ap) {
             const int idxi1 = occoff[h] + i;
             const int idxap1 = vir_comoff[h] + virpi[h] + ap;
             const int idxi2 = idxi1;
             const int idxap2 = nvir + cabsoff[h] + ap;
             T.set_element(idxi2,idxap2,T1_cabs.get_element(idxi1,idxap1));
           }
         }
     }
   }

   if (debug_ >= DefaultPrintThresholds::mostN2)
     T.print(prepend_spincase(spin,"CABS amplitude Tia").c_str());

//   // test for T for different symmetry
//   {
//   RefSCDimension rowdim = new SCDimension(T.nrow()+T.ncol());
//   RefSCDimension coldim = rowdim;
//
//   RefDiagSCMatrix evals = localkit->diagmatrix(rowdim);
//   RefSCMatrix evecs = localkit->matrix(rowdim, coldim);
//   evals.assign(0.0);
//
//   RefSymmSCMatrix T1 = localkit->symmmatrix(rowdim);
//   T1.assign(0.0);
//   //T1.copyRefSCMatrix(T1_ccsd);
//   for(int i = 0; i < nocc; i++) {
//     for(int a = 0; a < nvir; a++) {
//       T1.set_element(i,a,T.get_element(i,a));
//     }
//   }
//
//   T1.diagonalize(evals, evecs);
//   T1 = 0;
//   evals.print("evals of CABS Singles T1 ");
//
//   evals = 0;
//   evecs = 0;
//   }

  if (!r12intermediates_->T2_cc_computed()) {

    // D^i_j = - t^i_a * t^a_j
    for (int i = 0; i < nocc; ++i) {
      RefSCVector tia = T.get_row(i);

      for (int j = 0; j <= i; ++j) {
        RefSCVector tja = T.get_row(j);
        const double Dij = - tia.dot(tja);
        D.set_element(i, j, Dij);
      }
    }

    // D^a_b = t^a_i * t^i_b
    for (int a = 0; a < nvir; ++a) {
      RefSCVector tai = T.get_column(a);

      for (int b = 0; b <= a; ++b) {
        RefSCVector tbi = T.get_column(b);
        const double Dab = tai.dot(tbi);
        D.set_element(a+nocc, b+nocc, Dab);
      }
    }

    // D^a_i = t^a_i
    for (int i = 0; i < nocc; ++i) {
      for (int a = 0; a < nvir; ++a) {
        const double Dai = T.get_element(i,a);
        D.set_element(a+nocc, i, Dai);
      }
    }

  }
  else {
      // D^i_j = - t^i_a' * t^a'_j
      for (int i = 0; i < nocc; ++i) {
        RefSCVector tia = T.get_row(i);

        for (int j = 0; j <= i; ++j) {
          RefSCVector tja = T.get_row(j);

          for (int a = 0; a < nvir; ++a) {
            tia.set_element(a, 0.0);
            tja.set_element(a, 0.0);
          }
          const double Dij = - tia.dot(tja);
          D.set_element(i, j, Dij);
        }
      }
  }
  // end of if (!r12intermediates_->T2_cc_computed())

  // D^a'_b= t^a'_i * t^i_b
  for (int ap = 0; ap < ncabs; ++ap) {
    RefSCVector tapi = T.get_column(ap+nvir);

    for (int b = 0; b < nvir; ++b) {
      RefSCVector tbi = T.get_column(b);;
      const double Dapb = tapi.dot(tbi);
      D.set_element(ap+norbs, b+nocc, Dapb);
    }
  }

  // D^a'_b' = t^a'_i * t^i_b'
  for (int ap = 0; ap < ncabs; ++ap) {
    RefSCVector tapi = T.get_column(ap+nvir);

    for (int bp = 0; bp <= ap; ++bp) {
      RefSCVector tbpi = T.get_column(bp+nvir);
      const double Dapbp = tapi.dot(tbpi);
      D.set_element(ap+norbs, bp+norbs, Dapbp);
    }
  }
//  D.print(prepend_spincase(spin,"Da'b' CABS").c_str());


  // D^a'_i = t^a'_i
  for (int i = 0; i < nocc; ++i) {
    for (int ap = 0; ap < ncabs; ++ap) {
      const double Dapi = T.get_element(i,ap+nvir);
      D.set_element(ap+norbs, i, Dapi);
    }
  }

  // Symmetrize matrices
  for(int p = 0; p < norbs_com; p++) {
    for(int q = 0; q < p; q++) {
      D.set_element(q,p,D.get_element(p,q));
    }
  }
  //D.print(prepend_spincase(spin,"D CABS").c_str());
  return D;
}
// end of compute_D_CABS

// test function for computing one-electron density form CABS contribution
RefSCMatrix MP2R12Energy_Diag::compute_D_CABS_test(SpinCase1 spin) {

  // obtain occ, vir, orbs, and cabs orbitals in that order
  const SpinCase2 spincase = static_cast<SpinCase2>(spin+1);
  string spinletters = to_string(spincase);
  vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
  vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
  obtain_orbitals(spincase, v_orbs1, v_orbs2);

  // obtain orbitals
  const Ref<OrbitalSpace> occ = v_orbs1[4];
  const Ref<OrbitalSpace> vir = v_orbs1[1];
  const Ref<OrbitalSpace> cabs = v_orbs1[3];

  const int nocc = occ->rank();
  const int nvir = vir->rank();
  const int ncabs = cabs->rank();
  const int norbs = nocc + nvir;
  const int nvir_com = nvir + ncabs;
  const int norbs_com = nocc + nvir + ncabs;

  const RefSCMatrix T1_cabs = r12eval_->T1_cabs(spin);

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCDimension rowdim_occ = new SCDimension(nocc);
  RefSCDimension coldim_vir_com = new SCDimension(nvir_com);
  RefSCDimension coldim_cabs = new SCDimension(ncabs);
  RefSCMatrix T;
  if (!r12intermediates_->T2_cc_computed()) {
    T = localkit->matrix(rowdim_occ, coldim_vir_com);
  } else {
      T = localkit->matrix(rowdim_occ, coldim_cabs);
  }
  T.assign(0.0);

  const int nirreps = cabs->nblocks();
  const std::vector<unsigned int>& occpi = occ->block_sizes();
  const std::vector<unsigned int>& virpi = vir->block_sizes();
  const std::vector<unsigned int>& cabspi = cabs->block_sizes();

  std::vector<unsigned int> occoff(nirreps);
  std::vector<unsigned int> viroff(nirreps);
  std::vector<unsigned int> cabsoff(nirreps);
  std::vector<unsigned int> vir_comoff(nirreps);

  occoff[0] = 0;
  viroff[0] = 0;
  cabsoff[0] = 0;

  for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
    occoff[irrep] = occoff[irrep-1] + occpi[irrep-1];
    viroff[irrep] = viroff[irrep-1] + virpi[irrep-1];
    cabsoff[irrep] = cabsoff[irrep-1] + cabspi[irrep-1];
    vir_comoff[irrep] = viroff[irrep] + cabsoff[irrep];
  }

  for (int h = 0; h < nirreps; ++h) {

    if (r12intermediates_->T2_cc_computed()) {

      for (int i = 0; i < occpi[h]; ++i) {
        for (int ap = 0; ap < cabspi[h]; ++ap) {
          const int idxi = occoff[h] + i;
          const int idxap = cabsoff[h] + ap;
          T.set_element(idxi,idxap,T1_cabs.get_element(idxi,idxap));
        }
      }
    } else {
        for (int i = 0; i < occpi[h]; ++i) {
          for (int a = 0; a < virpi[h]; ++a) {
            const int idxi1 = occoff[h] + i;
            const int idxa1 = vir_comoff[h] + a;
            const int idxi2 = idxi1;
            const int idxa2 = viroff[h] + a;
            T.set_element(idxi2,idxa2,T1_cabs.get_element(idxi1,idxa1));
          }

          for (int ap = 0; ap < cabspi[h]; ++ap) {
            const int idxi1 = occoff[h] + i;
            const int idxap1 = vir_comoff[h] + virpi[h] + ap;
            const int idxi2 = idxi1;
            const int idxap2 = nvir + cabsoff[h] + ap;
            T.set_element(idxi2,idxap2,T1_cabs.get_element(idxi1,idxap1));
          }
        }
    }
  }
  T.print(prepend_spincase(spin,"CABS amplitude Tia").c_str());

  RefSCDimension rowdim = new SCDimension(norbs_com);
  RefSCDimension coldim = new SCDimension(norbs_com);
  RefSCMatrix D = localkit->matrix(rowdim, coldim);
  D.assign(0.0);

  if (!r12intermediates_->T2_cc_computed()) {
    // D^i_j = - t^i_a * t^a_j
    for (int i = 0; i < nocc; ++i) {
      for (int j = 0; j <= i; ++j) {

        double Dij = 0.0;
        for (int a = 0; a < nvir+ncabs; ++a){
          Dij += T.get_element(i,a) * T.get_element(j,a);
        }

        D.set_element(i, j, -Dij);
      }
    }

    // D^a_b = t^a_i * t^i_b
    for (int a = 0; a < nvir; ++a) {
      for (int b = 0; b <= a; ++b) {

        double Dab = 0.0;
        for (int i = 0; i < nocc; ++i) {
          Dab += T.get_element(i,a) * T.get_element(i,b);
        }
        D.set_element(a+nocc, b+nocc, Dab);
      }
    }

    // D^a'_b' = t^a'_i * t^i_b'
    for (int ap = 0; ap < ncabs; ++ap) {
      for (int bp = 0; bp <= ap; ++bp) {

        double Dapbp = 0.0;
        for (int i = 0; i < nocc; ++i) {
          Dapbp += T.get_element(i,ap+nvir) * T.get_element(i,bp+nvir);
        }
        D.set_element(ap+norbs, bp+norbs, Dapbp);
      }
    }

    // D^a_i = t^a_i & D^a'_i = t^a'_i
    for (int i = 0; i < nocc; ++i) {

      for (int a = 0; a < nvir; ++a) {
        const double Dai = T.get_element(i,a);
        D.set_element(a+nocc, i, Dai);
      }
      for (int ap = 0; ap < ncabs; ++ap) {
        const double Dapi = T.get_element(i,ap+nvir);
        D.set_element(ap+norbs, i, Dapi);
      }
    }

    // D^a'_b= t^a'_i * t^i_b
    for (int ap = 0; ap < ncabs; ++ap) {
      for (int b = 0; b < nvir; ++b) {

        double Dapb = 0.0;
        for (int i = 0; i < nocc; ++i) {
          Dapb += T.get_element(i,ap+nvir) * T.get_element(i,b);
        }
        D.set_element(ap+norbs, b+nocc, Dapb);
      }
    }

  } else {
      // D^i_j = - t^i_a' * t^a'_j
      for (int i = 0; i < nocc; ++i) {
        for (int j = 0; j <= i; ++j) {

          double Dij = 0;
          for (int ap = 0; ap < ncabs; ++ap) {
            Dij += T.get_element(i,ap) * T.get_element(j,ap);
          }

          D.set_element(i, j, -Dij);
        }
      }

      // D^a'_b' = t^a'_i * t^i_b'
      for (int ap = 0; ap < ncabs; ++ap) {
        for (int bp = 0; bp <= ap; ++bp) {

          double Dapbp = 0.0;
          for (int i = 0; i < nocc; ++i) {
            Dapbp +=  T.get_element(i,ap) * T.get_element(i,bp);
          }
          D.set_element(ap+norbs, bp+norbs, Dapbp);
        }
      }

      // D^a'_i = t^a'_i
      for (int i = 0; i < nocc; ++i) {
        for (int ap = 0; ap < ncabs; ++ap) {
          const double Dapi = T.get_element(i,ap);
          D.set_element(ap+norbs, i, Dapi);
        }
      }
  }

  // Symmetrize matrices
  for(int p = 0; p < norbs_com; p++) {
    for(int q = 0; q < p; q++) {
      D.set_element(q,p,D.get_element(p,q));
    }
  }

  return D;
}
// end of compute_D_CABS_test

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
//    compute_T2_mp2(v_orbs1, v_orbs2, T2ab_ij[s]);
//    // test function: compute T2 = T2(mp2) + T2(F12 corrected) directly
//     compute_T2abij_mp2f12(spincase2, C_0, C_1, T2ab_ij[s]);

    double* T2_mp2 = new double[nvir12 * nocc12];
    double* T2_f12corr = new double[nvir12 * nocc12];
    fill_n(T2_mp2, nvir12 * nocc12, 0.0);
    fill_n(T2_f12corr, nvir12 * nocc12, 0.0);

    compute_T2_mp2(spincase2, T2_mp2);
    compute_T2abij_f12corr(spincase2, C_0, C_1, T2_f12corr);
    compute_T2abij_mp2f12(nocc1_act, nocc2_act, nvir1, nvir2,
                          T2_mp2, T2_f12corr, T2ab_ij[s]);

    delete[] T2_mp2;
    delete[] T2_f12corr;


     // print T2^ab_ij
 #if 0
     string spinletters = to_string(spincase2);
     ExEnv::out0() << endl << spinletters << " MP2 T2^ab_ij" << endl
                   << "number of occupied orbital: " << nocc1_act << " " << nocc2_act << endl
                   << "number of virtual orbital: " << nvir1 << " " << nvir2 << endl;
//     print_T2abij_mp2(spinletters, "T2^ab_ij",
//                      nocc1_act, nocc2_act, nvir1, nvir2,
//                      T2ab_ij[s]);

     const double* iter_T2_test = T2ab_ij[s];
     for (int a = 0; a < 2; ++a) {
       for (int b = 0; b < 2; ++b) {
         ExEnv::out0() << endl << a << " " << b << endl;
         print_intermediate(spinletters, " T2 test", iter_T2_test, nocc1_act, nocc2_act);

         iter_T2_test += nocc1_act * nocc1_act;
       }
     }

 #endif
  } // end of compute T2^ab_ij

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

  // Alpha: 1/2 * C1 * (R^IJ_A'B * T^AB_IJ - R^JI_A'B * T^AB_IJ)  (AlphaAlpha part)
  //     + [(C0+C1)/2 * R^IJ_A'B + (C0-C1)/2 * R^JI_A'B] * T^AB_IJ  (AlphaBeta part)
  //
  // Beta:  1/2 * C1 * (R^IJ_A'B * T^AB_IJ- R^JI_A'B * T^AB_IJ)   (BetaBeta part)
  //     + [(C0+C1)/2 * R^IJ_BA'+ (C0-C1)/2 * R^JI_BA'] * T^BA_IJ  (AlphaBeta part)
  // where T^AB_IJ is antisymmetrized

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

//     Rapb_ij_ints->deactivate();
//     Rbap_ij_ints->deactivate();
     delete[] T2ab_ij[spincase];
     T2ab_ij[spincase] = NULL;

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
     RT2_1 = new double[ncabs_vir];
     RT2_2 = new double[ncabs_vir];
     fill_n(RT2_1, ncabs_vir, 0.0);
     fill_n(RT2_2, ncabs_vir, 0.0);

     if (nspincases2 == 3) {

//       RT2_1 = new double[ncabs_vir];
//       RT2_2 = new double[ncabs_vir];
//       fill_n(RT2_1, ncabs_vir, 0.0);
//       fill_n(RT2_2, ncabs_vir, 0.0);

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
         // R^ij_a'b * T^ab_ij
         compute_RTmp2_sum_3idx(RT13_23, f12_idx, nvir_alpha,
                                Rapb_ij_ints, T2ab_ij[AlphaBeta],
                                RT2_1);
         // R^ij_ba' * T^ab_ij
         compute_RTmp2_sum_3idx(RT31_23, f12_idx, nvir_alpha,
                                Rbap_ij_ints, T2ab_ij[AlphaBeta],
                                RT2_2);

//         // R^ij_ba' * T^ba_ij
//         compute_RTmp2_sum_3idx(RT31_32, f12_idx, nvir_beta,
//                                Rbap_ij_ints, T2ab_ij[AlphaBeta],
//                                RT2_1);
//         // R^ij_a'b * T^ba_ij
//         compute_RTmp2_sum_3idx(RT13_32, f12_idx, nvir_beta,
//                                Rapb_ij_ints, T2ab_ij[AlphaBeta],
//                                RT2_2);

//         RT2_1 = RT2apb_ab;
//         RT2_2 = RT2bap_ab;
    } // end of AlphaBeta part for D^a'_a
     Rapb_ij_ints->deactivate();
     Rbap_ij_ints->deactivate();

     // calculate D^a'_a RT2 part
     const double* iter_RT2apb_ab = RT2apb_ab;
     const double* iter_RT2bap_ab = RT2bap_ab;
     const double* iter_RT2_1 = RT2_1;
     const double* iter_RT2_2 = RT2_2;
     double* iter_RTmp2 = (spin == Alpha? RTmp2_alpha : RTmp2_beta);

     for (int idx1 = 0;  idx1 < ncabs; ++idx1) {
       for (int idx2 = 0; idx2 < nvir; ++idx2) {

         // AlphaAlpha/BetaBeta part
         double d_12 = 0.5 * C_1 * (*iter_RT2apb_ab - *iter_RT2bap_ab);

         ++iter_RT2apb_ab;
         ++iter_RT2bap_ab;

         // AlphaBeta part
         if (nocc_alpha != 0 && nocc_beta != 0) {
           d_12 += 0.5 * (C_0 + C_1) * (*iter_RT2_1) + 0.5 * (C_0 - C_1) * (*iter_RT2_2);

           ++iter_RT2_1;
           ++iter_RT2_2;
         } // end of AlphaBeta part

         *iter_RTmp2 = d_12;
         ++iter_RTmp2;
       } // end of looping over c
     } // end of calculating D^a'_a element

     //
     const double* const iter_D = (spin == Alpha? RTmp2_alpha : RTmp2_beta);
     //print_intermediate(spinletters, "D^A'_A RT2 (mp2) part", iter_D, ncabs, nvir);

     delete[] RT2apb_ab;
     delete[] RT2bap_ab;
     delete[] RT2_1;
     delete[] RT2_2;
   } // end of loop spincase1

  delete[] T2ab_ij[AlphaBeta];
  T2ab_ij[AlphaBeta] = NULL;
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
  const int nocc1 = occ1->rank();
  const int nocc2 = occ2->rank();
  const int nvir1 = vir1->rank();
  const int nvir2 = vir2->rank();
  const int ncabs1 = cabs1->rank();
  const int ncabs2 = cabs2->rank();
  const int nfzc1 = nocc1 - nocc1_act;
  const int nfzc2 = nocc2 - nocc2_act;

  const int norb = nocc1_act + nvir1;
  const int norb_com = nocc1_act + nvir1 + ncabs1;
  const int nribs_nfc  = nocc1 + nvir1 + ncabs1;

  ExEnv::out0() << endl << "**********************************************" << endl;
  ExEnv::out0() << endl << "Computing F12 one-particle density" << endl<< endl;

  // test propose
//  if (debug_ >= DefaultPrintThresholds::N2)
    ExEnv::out0() << endl << "number of alpha active occupied orbital: " << nocc1_act << endl
                          << "number of beta active occupied orbital: " << nocc2_act << endl
                          << "number of alpha virtual orbital: " << nvir1 << endl
                          << "number of beta virtual orbital: " << nvir2 << endl
                          << "number of orbs: " << norb << endl
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
//  double* Dm_i_alpha_RI = new double[nocc11];
//  fill_n(Dm_i_alpha_RI, nocc11, 0.0);
//
//  double* Dm_i_beta_RI;
//  if (nspincases1 == 2) {
//    const int nocc22 = nocc2_act * nocc2_act;
//    Dm_i_beta_RI = new double[nocc11];
//    fill_n(Dm_i_beta_RI, nocc11, 0.0);
//  }
  compute_Dmi(nspincases1, nspincases2, C_0, C_1,
              v_orbs1_ab, v_orbs2_ab,
              Dm_i_alpha, Dm_i_beta);

  if (debug_ >= DefaultPrintThresholds::N2) {
    print_intermediate("Alpha", "D^m_i RI", Dm_i_alpha, nocc1_act, nocc1_act);
    if (nspincases1 == 2)
      print_intermediate("Beta", "D^m_i RI", Dm_i_beta, nocc2_act, nocc2_act);
  }

//  // test for compute_Dmi
//  ExEnv::out0() << indent << "testing: D^i_i" << endl;
//  compute_Dii_test(nspincases1, nspincases2, nocc1_act, nocc2_act, C_0, C_1);
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
#if 1
  compute_Dmi_2(nspincases1, nspincases2, C_0, C_1,
                v_orbs1_ab, v_orbs2_ab,
                Dm_i_alpha, Dm_i_beta);

  if (debug_ >= DefaultPrintThresholds::N2) {
    print_intermediate("Alpha", "D^m_i", Dm_i_alpha, nocc1_act, nocc1_act);
    if (nspincases1 == 2)
      print_intermediate("Beta", "D^m_i", Dm_i_beta, nocc2_act, nocc2_act);
  }
#endif

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

  //
  // compute MP2 & F12 corrected MP2 T2 amplitudes
  // T2(MP2) = \bar{g}^ab_ij / (F^i_i + F^j_j - F^a_a - F^b_b)
  // T2(F12 corrected) = C^ab_ij / (F^i_i + F^j_j - F^a_a - F^b_b)
  vector<double*> T2_mp2(NSpinCases2);
  vector<double*> T2_f12corr(NSpinCases2);
  if (!r12intermediates_->T2_cc_computed()) {
    for(int s = 0; s < nspincases2; ++s) {
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      const SpinCase1 spin1 = case1(spincase2);
      const SpinCase1 spin2 = case2(spincase2);

      const int nocc_act_1 = (spin1 == Alpha? nocc1_act : nocc2_act);
      const int nocc_act_2 = (spin2 == Alpha? nocc1_act : nocc2_act);
      const int nvir_1 =  (spin1 == Alpha? nvir1 : nvir2);
      const int nvir_2 =  (spin2 == Alpha? nvir1 : nvir2);

      // test
//      ExEnv::out0() << endl << spincase2 << endl << "number of nocc_act_1: " << nocc_act_1 << endl
//                            << "number of nocc_act_2: " << nocc_act_2 << endl
//                            << "number of nvir_1: " << nvir_1 << endl
//                            << "number of nvir_2: " << nvir_2 << endl;

      const int nocc12vir12 = nocc_act_1 * nocc_act_2 * nvir_1 * nvir_2;

      if (nocc_act_1 == 0 || nocc_act_2 == 0)
        continue;

      T2_mp2[s] = new double[nocc12vir12];
      T2_f12corr[s] = new double[nocc12vir12];
      fill_n(T2_mp2[s], nocc12vir12, 0.0);
      fill_n(T2_f12corr[s], nocc12vir12, 0.0);

      compute_T2_mp2(spincase2, T2_mp2[s]);
      compute_T2abij_f12corr(spincase2, C_0, C_1, T2_f12corr[s]);
    }

    if (nspincases1 == 1) {
      T2_mp2[BetaBeta] = T2_mp2[AlphaAlpha];
      T2_f12corr[BetaBeta] = T2_f12corr[AlphaAlpha];
    }
  }

  double* Dap_i_alpha;
  double* Dap_i_beta;
  if (r12intermediates_->T1_cc_computed()) {
    const int ncabs1_occ1 = ncabs1 * nocc1_act;

    Dap_i_alpha = new double[ncabs1_occ1];
    fill_n(Dap_i_alpha, ncabs1_occ1, 0.0);

    if (nspincases1 == 2) {
      const int ncabs2_occ2 = ncabs2 * nocc2_act;

      Dap_i_beta = new double[ncabs2_occ2];
      fill_n(Dap_i_beta, ncabs2_occ2, 0.0);
    } else {
        Dap_i_beta = Dap_i_alpha;
    }

    compute_RT1_api(nspincases1, nspincases2, C_0, C_1,
                    v_orbs1_ab, v_orbs2_ab,
                    Dap_i_alpha, Dap_i_beta);
  }
  //
  // compute MP2 or CCSD, F12 one electron densities & dipole moments
  RefSCDimension rowdim = new SCDimension(norb_com);
  RefSCDimension coldim = new SCDimension(norb_com);
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix D[NSpinCases1];

  RefSCVector dipoles_f12 = localkit->vector(RefSCDimension(new SCDimension(3)));
  RefSCVector dipoles_cabs = localkit->vector(RefSCDimension(new SCDimension(3)));
  dipoles_f12.assign(0.0);
  dipoles_cabs.assign(0.0);

  RefSCVector dipoles_mp2;
  RefSCVector dipoles_ccsd;
  RefSCVector dipoles_or_relax;
  if (!r12intermediates_->T2_cc_computed()) {
    dipoles_mp2 = localkit->vector(RefSCDimension(new SCDimension(3)));
    dipoles_mp2.assign(0.0);
  } else {
      dipoles_ccsd = localkit->vector(RefSCDimension(new SCDimension(3)));
      dipoles_ccsd.assign(0.0);
      if (r12intermediates_->Onerdm_relax_computed()) {
        dipoles_or_relax = localkit->vector(RefSCDimension(new SCDimension(3)));
        dipoles_or_relax.assign(0.0);
      }
  }

  for (int s = 0; s < nspincases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);
    const string spin_label = (spin == Alpha? "Alpha case:" : "Beta case:");
    ExEnv::out0() << endl << endl << spin_label << endl;

    const int nocc_act = (spin == Alpha? nocc1_act : nocc2_act);
    const int nvir = (spin == Alpha? nvir1 : nvir2);
    const int ncabs = (spin == Alpha? ncabs1 : ncabs2);

    D[spin] = localkit->matrix(rowdim, coldim);
    D[spin].assign(0.0);

    // Test: dipole momement from each contribution
    RefSCMatrix D_occ_act = localkit->matrix(rowdim, coldim);
    RefSCMatrix D_vir = localkit->matrix(rowdim, coldim);
    RefSCMatrix D_cabs = localkit->matrix(rowdim, coldim);
    RefSCMatrix D_cabsvir = localkit->matrix(rowdim, coldim);
    RefSCMatrix D_ebc = localkit->matrix(rowdim, coldim);;
    D_occ_act.assign(0.0);
    D_vir.assign(0.0);
    D_cabs.assign(0.0);
    D_cabsvir.assign(0.0);
    D_ebc.assign(0.0);

    const double* iter_Dmi = (spin == Alpha? Dm_i_alpha : Dm_i_beta);
    for (int i = 0; i < nocc_act; ++i) {
        for (int j = 0; j < nocc_act; ++j, ++iter_Dmi){
          D[spin].set_element(i, j, *iter_Dmi);

          D_occ_act.set_element(i, j, *iter_Dmi); // test propose
        }
    }
    //ExEnv::out0() << endl << " trace of Dij: " << scprintf("%12.10f", D[spin].trace()) << endl;
    ExEnv::out0() << endl << " trace of D_occ_act: " << scprintf("%12.10f", D_occ_act.trace()) << endl;

    const double* iter_Dcb = (spin == Alpha? Dc_b_alpha : Dc_b_beta);
    for (int a = nocc_act; a < norb; ++a) {
        for (int b = nocc_act; b < norb; ++b, ++iter_Dcb){
          D[spin].set_element(a, b, *iter_Dcb);

          D_vir.set_element(a, b, *iter_Dcb); // test
        }
    }
    //D[spin].print(prepend_spincase(spin,"F12 one-particle density Dij+Dab:").c_str());
    //ExEnv::out0() << endl << " trace of Dij + Dab: " << scprintf("%12.10f", D[spin].trace()) << endl;
    ExEnv::out0() << endl << " trace of D_vir: " << scprintf("%12.10f", D_vir.trace()) << endl;

    const double* iter_Dcpbp_a = (spin == Alpha? Dcp_bp_alpha_A : Dcp_bp_beta_A);
    const double* iter_Dcpbp_ap = (spin == Alpha? Dcp_bp_alpha_Ap : Dcp_bp_beta_Ap);
    for (int ap = norb; ap < norb_com; ++ap) {
        for (int bp = norb; bp < norb_com; ++bp, ++iter_Dcpbp_a, ++iter_Dcpbp_ap) {
          const double Dapbp =  *iter_Dcpbp_a + *iter_Dcpbp_ap;
          D[spin].set_element(ap, bp, Dapbp);

          D_cabs.set_element(ap, bp, Dapbp); // test
        }
    }
    //D[spin].print(prepend_spincase(spin,"F12 one-particle density Dij+Dab+Da'b':").c_str());
    //ExEnv::out0() << endl << " trace of Dij + Dab + Da'b': " << scprintf("%12.10f", D[spin].trace()) << endl;
    ExEnv::out0() << endl << " trace of D_cabs: " << scprintf("%12.10f", D_cabs.trace()) << endl;

    if (this->r12eval()->ebc() == false
        || this->r12eval()->coupling() == true) {

      const double* iter_Dap_a_RR = (spin == Alpha? Dap_a_alpha_RR : Dap_a_beta_RR);
      const double* iter_Dap_a_RT = (spin == Alpha? Dap_a_alpha_RT : Dap_a_beta_RT);
      for (int ap = norb; ap < norb_com; ++ap) {
        for (int a = nocc_act; a < norb; ++a, ++iter_Dap_a_RR, ++iter_Dap_a_RT) {
          const double Dapa = *iter_Dap_a_RR + *iter_Dap_a_RT;
          D[spin].set_element(ap, a, Dapa);
          D[spin].set_element(a, ap, Dapa);

          D_cabsvir.set_element(ap, a, *iter_Dap_a_RR); // test
          D_cabsvir.set_element(a, ap, *iter_Dap_a_RR); // test
          D_ebc.set_element(ap, a, *iter_Dap_a_RT); // test
          D_ebc.set_element(a, ap, *iter_Dap_a_RT); // test
        }
      }

    } else {
        const double* iter_Dap_a_RR = (spin == Alpha? Dap_a_alpha_RR : Dap_a_beta_RR);
        for (int ap = norb; ap < norb_com; ++ap) {
          for (int a = nocc_act; a < norb; ++a, ++iter_Dap_a_RR) {
            const double Dapa = *iter_Dap_a_RR;
            D[spin].set_element(ap, a, Dapa);
            D[spin].set_element(a, ap, Dapa);

            D_cabsvir.set_element(ap, a, Dapa); // test
            D_cabsvir.set_element(a, ap, Dapa); // test
          }
        }
    }

    if (debug_ >= DefaultPrintThresholds::allN2)
      D[spin].print(prepend_spincase(spin,"F12 one-particle density:").c_str());

    // Test code for D: trace of D = 0, as trace of Dij = - trace of (Dab + Da'b')
    //ExEnv::out0() << endl << "Trace of D_f12: " << scprintf("%12.10f", D[spin].trace())<< endl;
    ExEnv::out0() << endl << " trace of D_cabsvir: " << scprintf("%12.10f", D_cabsvir.trace()) << endl;

    //
    // Obtain dipole ints in ribs
    //
    // Dipole ints in occ_act, vir, cabs space
    const RefSCDimension M_dim_ribs(new SCDimension(norb_com));
    RefSCMatrix MX_nb_ribs = localkit->matrix(M_dim_ribs,M_dim_ribs);
    RefSCMatrix MY_nb_ribs = localkit->matrix(M_dim_ribs,M_dim_ribs);
    RefSCMatrix MZ_nb_ribs = localkit->matrix(M_dim_ribs,M_dim_ribs);
    MX_nb_ribs.assign(0.0);
    MY_nb_ribs.assign(0.0);
    MZ_nb_ribs.assign(0.0);

    // Dipole ints in occ, vir, cabs space
    const RefSCDimension dim_ribs_nfc(new SCDimension(nribs_nfc));
    RefSCMatrix MX_ribs_nfc = localkit->matrix(dim_ribs_nfc,dim_ribs_nfc);
    RefSCMatrix MY_ribs_nfc = localkit->matrix(dim_ribs_nfc,dim_ribs_nfc);
    RefSCMatrix MZ_ribs_nfc = localkit->matrix(dim_ribs_nfc,dim_ribs_nfc);
    MX_ribs_nfc.assign(0.0);
    MY_ribs_nfc.assign(0.0);
    MZ_ribs_nfc.assign(0.0);

    form_DipoleInts_inRibs(spin, v_orbs1_ab, v_orbs2_ab,
                           MX_nb_ribs, MY_nb_ribs, MZ_nb_ribs,
                           MX_ribs_nfc, MY_ribs_nfc, MZ_ribs_nfc);

    // test
    {
      const int nfzc = (spin == Alpha? nfzc1 : nfzc2);
      const int norbs_tot = norb + nfzc;
      const RefSCDimension M_dim_orbs(new SCDimension(norbs_tot));
      RefSCMatrix MX_nb_orbs = localkit->matrix(M_dim_orbs,M_dim_orbs);
      RefSCMatrix MY_nb_orbs = localkit->matrix(M_dim_orbs,M_dim_orbs);
      RefSCMatrix MZ_nb_orbs = localkit->matrix(M_dim_orbs,M_dim_orbs);
      MX_nb_orbs.assign(0.0);
      MY_nb_orbs.assign(0.0);
      MZ_nb_orbs.assign(0.0);

      const Ref<OrbitalSpace>& space_orbs = (spin == Alpha? orbs1 : orbs2);
      RefSCMatrix MX_orbs, MY_orbs, MZ_orbs;
      RefSCMatrix  MXX, MYY, MZZ, MXY, MXZ, MYZ;
      compute_multipole_ints(space_orbs, space_orbs,
                             MX_orbs, MY_orbs, MZ_orbs,
                             MXX, MYY, MZZ,
                             MXY, MXZ, MYZ);
      //MZ_orbs.print(prepend_spincase(spin,"MZ_orbs").c_str());
      MXX = 0;
      MYY = 0;
      MZZ = 0;
      MXY = 0;
      MXZ = 0;
      MYZ = 0;
    }

    //MZ_ribs_nfc.print(prepend_spincase(spin,"MZ_ribs_nfc").c_str());

    // Test codes: dipole moment from each contribution
#if 1
    RefSCMatrix opdm_occ_act = onepdm_transformed(spin, true, D_occ_act);
    D_occ_act = 0;

    RefSCMatrix opdm_x_occ_act = MX_nb_ribs * opdm_occ_act;
    RefSCMatrix opdm_y_occ_act = MY_nb_ribs * opdm_occ_act;
    RefSCMatrix opdm_z_occ_act = MZ_nb_ribs * opdm_occ_act;
    opdm_occ_act = 0;

    const double dx_occ_act = opdm_x_occ_act.trace();
    const double dy_occ_act = opdm_y_occ_act.trace();
    const double dz_occ_act = opdm_z_occ_act.trace();
    opdm_x_occ_act = 0;
    opdm_y_occ_act = 0;
    opdm_z_occ_act = 0;
    ExEnv::out0() << endl << "x y z dipole moments from occ_act: "
                          << scprintf("%12.10f", dx_occ_act) << "  "
                          << scprintf("%12.10f", dy_occ_act) << "  "
                          << scprintf("%12.10f", dz_occ_act) << endl;

    RefSCMatrix opdm_vir = onepdm_transformed(spin, true, D_vir);
    D_vir = 0;

    RefSCMatrix opdm_x_vir = MX_nb_ribs * opdm_vir;
    RefSCMatrix opdm_y_vir = MY_nb_ribs * opdm_vir;
    RefSCMatrix opdm_z_vir = MZ_nb_ribs * opdm_vir;
    opdm_vir = 0;

    const double dx_vir = opdm_x_vir.trace();
    const double dy_vir = opdm_y_vir.trace();
    const double dz_vir = opdm_z_vir.trace();
    opdm_x_vir = 0;
    opdm_y_vir = 0;
    opdm_z_vir = 0;
    ExEnv::out0() << endl << "x y z dipole moments from vir: "
                          << scprintf("%12.10f", dx_vir) << "  "
                          << scprintf("%12.10f", dy_vir) << "  "
                          << scprintf("%12.10f", dz_vir) << endl;

    RefSCMatrix opdm_orbs_cabs = onepdm_transformed(spin, true, D_cabs);
    D_cabs = 0;

    RefSCMatrix opdm_x_cabs = MX_nb_ribs * opdm_orbs_cabs;
    RefSCMatrix opdm_y_cabs = MY_nb_ribs * opdm_orbs_cabs;
    RefSCMatrix opdm_z_cabs = MZ_nb_ribs * opdm_orbs_cabs;
    opdm_orbs_cabs = 0;

    const double dx_cabs = opdm_x_cabs.trace();
    const double dy_cabs = opdm_y_cabs.trace();
    const double dz_cabs = opdm_z_cabs.trace();
    opdm_x_cabs = 0;
    opdm_y_cabs = 0;
    opdm_z_cabs = 0;
    ExEnv::out0() << endl << "x y z dipole moments from cabs: "
                          << scprintf("%12.10f", dx_cabs) << "  "
                          << scprintf("%12.10f", dy_cabs) << "  "
                          << scprintf("%12.10f", dz_cabs) << endl;

    RefSCMatrix opdm_cabsvir = onepdm_transformed(spin, true, D_cabsvir);
    D_cabsvir = 0;

    RefSCMatrix opdm_x_cabsvir = MX_nb_ribs * opdm_cabsvir;
    RefSCMatrix opdm_y_cabsvir = MY_nb_ribs * opdm_cabsvir;
    RefSCMatrix opdm_z_cabsvir = MZ_nb_ribs * opdm_cabsvir;
    opdm_cabsvir = 0;

    const double dx_cabsvir = opdm_x_cabsvir.trace();
    const double dy_cabsvir = opdm_y_cabsvir.trace();
    const double dz_cabsvir = opdm_z_cabsvir.trace();
    opdm_x_cabsvir = 0;
    opdm_y_cabsvir = 0;
    opdm_z_cabsvir = 0;
    ExEnv::out0() << endl << "x y z dipole moments from cabsvir: "
                          << scprintf("%12.10f", dx_cabsvir) << "  "
                          << scprintf("%12.10f", dy_cabsvir) << "  "
                          << scprintf("%12.10f", dz_cabsvir) << endl;

    if (this->r12eval()->ebc() == false
        || this->r12eval()->coupling() == true) {

      RefSCMatrix opdm_ebc = onepdm_transformed(spin, true, D_ebc);
      D_ebc = 0;

      RefSCMatrix opdm_x_ebc = MX_nb_ribs * opdm_ebc;
      RefSCMatrix opdm_y_ebc = MY_nb_ribs * opdm_ebc;
      RefSCMatrix opdm_z_ebc = MZ_nb_ribs * opdm_ebc;
      opdm_ebc = 0;

      const double dx_ebc = opdm_x_ebc.trace();
      const double dy_ebc = opdm_y_ebc.trace();
      const double dz_ebc = opdm_z_ebc.trace();
      opdm_x_ebc = 0;
      opdm_y_ebc = 0;
      opdm_z_ebc = 0;
      ExEnv::out0() << endl << "x y z dipole moments from ebc (false): "
                            << scprintf("%12.10f", dx_ebc) << "  "
                            << scprintf("%12.10f", dy_ebc) << "  "
                            << scprintf("%12.10f", dz_ebc) << endl;
    }


    if (r12intermediates_->T2_cc_computed()) {
      RefSCMatrix D_cabsocc = localkit->matrix(rowdim, coldim);;
      D_cabsocc.assign(0.0);

      const double* iter_Dap_i = (spin == Alpha? Dap_i_alpha : Dap_i_beta);
      for (int ap = norb; ap < norb_com; ++ap) {
        for (int i = 0; i < nocc_act; ++i, ++iter_Dap_i) {
          const double Dapi = *iter_Dap_i;
          D_cabsocc.set_element(ap, i, Dapi);
          D_cabsocc.set_element(i, ap, Dapi);
        }
      }
      //D_cabsocc.print("D_cabsocc");
      ExEnv::out0() << endl << "Trace of D^a'_i: " << scprintf("%12.10f", D_cabsocc.trace())<< endl;

      RefSCMatrix opdm_cabsocc = onepdm_transformed(spin, true, D_cabsocc);
      D_cabsocc = 0;

      RefSCMatrix opdm_x_cabsocc = MX_nb_ribs * opdm_cabsocc;
      RefSCMatrix opdm_y_cabsocc = MY_nb_ribs * opdm_cabsocc;
      RefSCMatrix opdm_z_cabsocc = MZ_nb_ribs * opdm_cabsocc;
      opdm_cabsocc = 0;

      const double dx_cabsocc = opdm_x_cabsocc.trace();
      const double dy_cabsocc = opdm_y_cabsocc.trace();
      const double dz_cabsocc = opdm_z_cabsocc.trace();
      opdm_x_cabsocc = 0;
      opdm_y_cabsocc = 0;
      opdm_z_cabsocc = 0;
      ExEnv::out0() << endl << "x y z dipole moments from cabsocc: "
                            << scprintf("%12.10f", dx_cabsocc) << "  "
                            << scprintf("%12.10f", dy_cabsocc) << "  "
                            << scprintf("%12.10f", dz_cabsocc) << endl;
    }
#endif

    // Transform D to MPQC ordering
    // RefsymmMatrix ??
    RefSCMatrix opdm_f12 = onepdm_transformed(spin, true, D[spin]);
    //opdm_f12.print("F12 one-particle density matrix MPQC ordering");

    RefSCMatrix opdm_x_f12 = MX_nb_ribs * opdm_f12;
    RefSCMatrix opdm_y_f12 = MY_nb_ribs * opdm_f12;
    RefSCMatrix opdm_z_f12 = MZ_nb_ribs * opdm_f12;
    opdm_f12 = 0;

    const double dx_f12 = opdm_x_f12.trace();
    const double dy_f12 = opdm_y_f12.trace();
    const double dz_f12 = opdm_z_f12.trace();
    opdm_x_f12 = 0;
    opdm_y_f12 = 0;
    opdm_z_f12 = 0;

    dipoles_f12[0] = dipoles_f12[0] + dx_f12;
    dipoles_f12[1] = dipoles_f12[1] + dy_f12;
    dipoles_f12[2] = dipoles_f12[2] + dz_f12;

    // Obtain one-particle density from CABS contribution
    RefSCMatrix D_CABS_single = compute_D_CABS(spin);
    //D_CABS_single.print(prepend_spincase(spin,"one-particle density from CABS contribution:").c_str());
    ExEnv::out0() << endl << endl
                  << "Trace of D_cabs_single: " << scprintf("%12.10f", D_CABS_single.trace())<< endl;

//    // add CABS contribution to D
//    D[spin].accumulate(D_cabs);
//    if (debug_ >= DefaultPrintThresholds::mostN2)
//      D[spin].print(prepend_spincase(spin,"F12 one-particle density with CABS contribution:").c_str());
#if 0
    // test for D_CABS
    const RefSCMatrix D_cabs_test = compute_D_CABS_test(spin);
    D_cabs_test.print(prepend_spincase(spin,"Test: one-particle density from CABS contribution:").c_str());
#endif

    RefSCMatrix opdm_cabs_single = onepdm_transformed(spin, false, D_CABS_single);
    //opdm_cabs_single.print(prepend_spincase(spin,"opdm_cabs_single").c_str());
    //print_Mathematica_form(opdm_cabs_single);
    D_CABS_single = 0;
    RefSCMatrix opdm_x_cabs_single = MX_ribs_nfc * opdm_cabs_single;
    RefSCMatrix opdm_y_cabs_single = MY_ribs_nfc * opdm_cabs_single;
    RefSCMatrix opdm_z_cabs_single = MZ_ribs_nfc * opdm_cabs_single;

//    // test for MX_ribs_nfc
//    RefSCMatrix opdm_x_cabs_single = MX_nb_ribs * opdm_cabs_single;
//    RefSCMatrix opdm_y_cabs_single = MY_nb_ribs * opdm_cabs_single;
//    RefSCMatrix opdm_z_cabs_single = MZ_nb_ribs * opdm_cabs_single;

//    // test for CABS Singles contribution
//    {
//    RefSCDimension rowdim = new SCDimension(opdm_cabs_single.nrow());
//    RefSCDimension coldim = new SCDimension(opdm_cabs_single.ncol());
//
//    RefDiagSCMatrix evals = localkit->diagmatrix(rowdim);
//    RefSCMatrix evecs = localkit->matrix(rowdim, coldim);
//    evals.assign(0.0);
//
//    RefSymmSCMatrix DCabsSingles = localkit->symmmatrix(rowdim);
//    DCabsSingles.assign(0.0);
//    DCabsSingles.copyRefSCMatrix(opdm_cabs_single);
//    //DCabsSingles.print(prepend_spincase(spin,"one-particle density from CABS contribution symmetric version:").c_str());
//    ExEnv::out0() << endl
//                  << "Trace of opdm_cabs_single(RefSymmSCMatrix): " << scprintf("%12.10f", opdm_cabs_single.trace())<< endl;
//
//    DCabsSingles.diagonalize(evals, evecs);
//    DCabsSingles = 0;
//    evals.print("evals of opdm_cabs_single");
//
//    RefSymmSCMatrix MZRibs = localkit->symmmatrix(rowdim);
//    MZRibs.assign(0.0);
//
//    ExEnv::out0() << endl << endl
//                  << "Trace of MZ_ribs_nfc: " << scprintf("%12.10f", MZ_ribs_nfc.trace())<< endl;
////    MXRibs.copyRefSCMatrix(MX_ribs_nfc);
//    MZRibs.copyRefSCMatrix(MZ_ribs_nfc);
//    ExEnv::out0() << endl
//                  << "Trace of MXRibs(RefSymmSCMatrix): " << scprintf("%12.10f", MZRibs.trace())<< endl;
//
//    evals.assign(0.0);
//    MZRibs.diagonalize(evals, evecs);
//    MZRibs = 0;
//    evals.print("evals of MZ");
//
//    RefSymmSCMatrix opdm_z = localkit->symmmatrix(rowdim);
//    opdm_z.assign(0.0);
//    opdm_z.copyRefSCMatrix(opdm_z_cabs_single);
//
//    evals.assign(0.0);
//    opdm_z.diagonalize(evals, evecs);
//    opdm_z = 0;
//    evals.print("evals of opdm_z_cabs_single");
//
//    evals = 0;
//    evecs = 0;
//    }
//    // end of test

    opdm_cabs_single = 0;
    MX_ribs_nfc = 0;
    MY_ribs_nfc = 0;
    MZ_ribs_nfc = 0;

    const double dx_cabs_single = opdm_x_cabs_single.trace();
    const double dy_cabs_single = opdm_y_cabs_single.trace();
    const double dz_cabs_single = opdm_z_cabs_single.trace();
    opdm_x_cabs_single = 0;
    opdm_y_cabs_single = 0;
    opdm_z_cabs_single = 0;
    ExEnv::out0() << endl << "x y z dipole moments from cabs_single: "
                          << scprintf("%12.10f", dx_cabs_single) << "  "
                          << scprintf("%12.10f", dy_cabs_single) << "  "
                          << scprintf("%12.10f", dz_cabs_single) << endl;

    dipoles_cabs[0] = dipoles_cabs[0] + dx_cabs_single;
    dipoles_cabs[1] = dipoles_cabs[1] + dy_cabs_single;
    dipoles_cabs[2] = dipoles_cabs[2] + dz_cabs_single;

    // Compute PSI CCSD (+ CCSD_F12 orbital relaxation)
    //      or MP2 density to r12 one-particle density
    if (r12intermediates_->T2_cc_computed()) {

      // Obtain CCSD one-particle density from Psi
      RefSCMatrix D_cc = r12intermediates_->get_1rdm_cc(spin);
      // test: trace of D_cc = # of electrons
      ExEnv::out0() << endl << "Trace of D_cc: " << scprintf("%12.10f", D_cc.trace()) << endl;
//      D[spin].accumulate_subblock(D_cc, 0, norb-1, 0, norb-1, 0, 0); // problem: D_cc include the frozen occ
//      if (debug_ >= DefaultPrintThresholds::allN2)
//      D[spin].print(prepend_spincase(spin,"CCSD_F12 one-particle density:").c_str());

      const int nfzc = (spin == Alpha? nfzc1 : nfzc2);
      const int nocc = nocc_act + nfzc;
      const int norbs_tot = norb + nfzc;

      // Test code: compute the ccsd dipole moment
  #if 1
      const RefSCDimension M_dim_orbs(new SCDimension(norbs_tot));
      RefSCMatrix MX_nb_orbs = localkit->matrix(M_dim_orbs,M_dim_orbs);
      RefSCMatrix MY_nb_orbs = localkit->matrix(M_dim_orbs,M_dim_orbs);
      RefSCMatrix MZ_nb_orbs = localkit->matrix(M_dim_orbs,M_dim_orbs);
      MX_nb_orbs.assign(0.0);
      MY_nb_orbs.assign(0.0);
      MZ_nb_orbs.assign(0.0);

      const Ref<OrbitalSpace>& space_orbs = (spin == Alpha? orbs1 : orbs2);
      RefSCMatrix MX_orbs, MY_orbs, MZ_orbs;
      RefSCMatrix  MXX, MYY, MZZ, MXY, MXZ, MYZ;
      compute_multipole_ints(space_orbs, space_orbs,
                             MX_orbs, MY_orbs, MZ_orbs,
                             MXX, MYY, MZZ,
                             MXY, MXZ, MYZ);
      MXX = 0;
      MYY = 0;
      MZZ = 0;
      MXY = 0;
      MXZ = 0;
      MYZ = 0;

      for(int p = 0; p < norbs_tot; p++) {
        for(int q = 0; q < norbs_tot; q++) {
          MX_nb_orbs.set_element(q,p,MX_orbs.get_element(p,q));
          MY_nb_orbs.set_element(q,p,MY_orbs.get_element(p,q));
          MZ_nb_orbs.set_element(q,p,MZ_orbs.get_element(p,q));
        }
      }
      MX_orbs = 0;
      MY_orbs = 0;
      MZ_orbs = 0;
      //MZ_nb_orbs.print(prepend_spincase(spin,"mu(Z)_nb in orbs").c_str());

//      // test for CABS Singles contribution
//      {
//      RefSCDimension rowdim = new SCDimension(MZ_nb_orbs.nrow());
//      RefSCDimension coldim = new SCDimension(MZ_nb_orbs.ncol());
//
//      RefDiagSCMatrix evals = localkit->diagmatrix(rowdim);
//      RefSCMatrix evecs = localkit->matrix(rowdim, coldim);
//      evals.assign(0.0);
//
//      RefSymmSCMatrix MZ = localkit->symmmatrix(rowdim);
//      MZ.assign(0.0);
//      MZ.copyRefSCMatrix(MZ_nb_orbs);
//
//      MZ.diagonalize(evals, evecs);
//      MZ = 0;
//      evals.print("evals of MZ_nb_orbs");
//
//      evals = 0;
//      evecs = 0;
//      }
//      // end of test

    // copy psi ccsd density in MPQC ordering back
    RefSCMatrix opdm_cc = onepdm_transformed2(spin, D_cc);
    //opdm_cc.print("CCSD one-particle density matrix in MPQC ordering");

    RefSCMatrix opdm_x_cc = MX_nb_orbs * opdm_cc;
    RefSCMatrix opdm_y_cc = MY_nb_orbs * opdm_cc;
    RefSCMatrix opdm_z_cc = MZ_nb_orbs * opdm_cc;
    opdm_cc = 0;

    const double dx_cc = opdm_x_cc.trace();
    const double dy_cc = opdm_y_cc.trace();
    const double dz_cc = opdm_z_cc.trace();
    // clean
    opdm_x_cc = 0;
    opdm_y_cc = 0;
    opdm_z_cc = 0;
    ExEnv::out0() << endl << "x y z dipole moments from CCSD: "
                          << scprintf("%12.10f", dx_cc) << "  "
                          << scprintf("%12.10f", dy_cc) << "  "
                          << scprintf("%12.10f", dz_cc) << endl;

    dipoles_ccsd[0] = dipoles_ccsd[0] + dx_cc;
    dipoles_ccsd[1] = dipoles_ccsd[1] + dy_cc;
    dipoles_ccsd[2] = dipoles_ccsd[2] + dz_cc;
#endif

    // Obtain CCSD_F12 orbital relaxation contribution to 1rdm
    if (r12intermediates_->Onerdm_relax_computed()) {
      RefSCMatrix D_ccsdf12_relax = r12intermediates_->get_1rdm_relax(spin);
      RefSCMatrix D_relax = localkit->matrix(D_cc.rowdim(), D_cc.coldim());
      D_relax.assign(0.0);

      for (int a = 0; a < nvir; ++a) {
        for (int i = 0; i < nocc; ++i){
          const double Dai = D_ccsdf12_relax.get_element(a,i);
          const int idx_a = a + nocc;
          D_relax.set_element(idx_a, i, Dai);
          D_relax.set_element(i, idx_a, Dai);
        }
      }
      //D_relax.print(prepend_spincase(spin,"one-particle density from orbital relax contribution:").c_str());

      // copy D_relax in MPQC ordering back
      RefSCMatrix opdm_relax = onepdm_transformed2(spin, D_relax);

      RefSCMatrix opdm_x_relax = MX_nb_orbs * opdm_relax;
      RefSCMatrix opdm_y_relax = MY_nb_orbs * opdm_relax;
      RefSCMatrix opdm_z_relax = MZ_nb_orbs * opdm_relax;
      opdm_relax = 0;

      const double dx_relax = opdm_x_relax.trace();
      const double dy_relax = opdm_y_relax.trace();
      const double dz_relax = opdm_z_relax.trace();
      // clean
      opdm_x_relax = 0;
      opdm_y_relax = 0;
      opdm_z_relax = 0;
      ExEnv::out0() << endl << "x y z dipole moments from CCSD_F12 orbital relaxtion: "
                            << scprintf("%12.10f", dx_relax) << "  "
                            << scprintf("%12.10f", dy_relax) << "  "
                            << scprintf("%12.10f", dz_relax) << endl;

      dipoles_or_relax[0] = dipoles_or_relax[0] + dx_relax;
      dipoles_or_relax[1] = dipoles_or_relax[1] + dy_relax;
      dipoles_or_relax[2] = dipoles_or_relax[2] + dz_relax;
    }
    MX_nb_orbs = 0;
    MY_nb_orbs = 0;
    MZ_nb_orbs = 0;
    } else {
//        RefSCMatrix Dmp2 = compute_1rdm_mp2(spin);
//        //Dmp2.print("MP2 one-electron density matrix");
//        RefSCMatrix D_mp2_test = compute_1rdm_mp2_test(spin);

        //
        // MP2F12 one-electron density and dipole moment
        // MP2 part: T(MP2)T(MP2)
        const SpinCase2 spincase2 = static_cast<SpinCase2>(s+1);
        RefSCMatrix Dmp2mp2 = compute_1rdm_mp2part(spin, nocc1_act, nocc2_act, nvir1, nvir2,
                                                   T2_mp2[spincase2],T2_mp2[AlphaBeta]);
        //Dmp2f12.print("Dmp2_mp2 density matrix");
        //ExEnv::out0() << endl << spin_label << " Trace of D_MP2F12 (MP2 part): " << scprintf("%12.10f", Dmp2f12.trace())<< endl;

        RefSCMatrix D_MP2 = localkit->matrix(rowdim, coldim);
        D_MP2.assign(0.0);
        D_MP2.accumulate_subblock(Dmp2mp2, 0, norb-1, 0, norb-1, 0, 0);
        Dmp2mp2 = 0;
        ExEnv::out0() << endl << endl
                      << "Trace of D_MP2: " << scprintf("%12.10f", D_MP2.trace())<< endl;

        RefSCMatrix opdm_mp2 = onepdm_transformed(spin, true, D_MP2);
        D_MP2 = 0;
        RefSCMatrix opdm_x_mp2 = MX_nb_ribs * opdm_mp2;
        RefSCMatrix opdm_y_mp2 = MY_nb_ribs * opdm_mp2;
        RefSCMatrix opdm_z_mp2 = MZ_nb_ribs * opdm_mp2;
        opdm_mp2 = 0;
        const double dx_mp2 = opdm_x_mp2.trace();
        const double dy_mp2 = opdm_y_mp2.trace();
        const double dz_mp2 = opdm_z_mp2.trace();
        opdm_x_mp2 = 0;
        opdm_y_mp2 = 0;
        opdm_z_mp2 = 0;
        ExEnv::out0() << endl << "x y z dipole moments of mp2: "
                               << scprintf("%12.10f", dx_mp2) << "  "
                               << scprintf("%12.10f", dy_mp2) << "  "
                               << scprintf("%12.10f", dz_mp2) << endl;

        // + T(MP2)T(F12)
        RefSCMatrix Dmp2_f12corr = compute_1rdm_mp2part(spin, nocc1_act, nocc2_act, nvir1, nvir2,
                                                   T2_mp2[spincase2], T2_f12corr[spincase2],
                                                   T2_mp2[AlphaBeta], T2_f12corr[AlphaBeta]);
        //Dmp2_f12.print("Dmp2_f12 density matrix");
        //Dmp2f12.accumulate(Dmp2_f12);
        //Dmp2_f12 = 0;

        // + T(F12)T(MP2)
        RefSCMatrix Df12_mp2 = compute_1rdm_mp2part(spin, nocc1_act, nocc2_act, nvir1, nvir2,
                                                    T2_f12corr[spincase2], T2_mp2[spincase2],
                                                    T2_f12corr[AlphaBeta], T2_mp2[AlphaBeta]);
        //Df12_mp2.print("Df12_mp2 density matrix");
        Dmp2_f12corr.accumulate(Df12_mp2);
        Df12_mp2 = 0;

        // + T(F12)T(F12)
        RefSCMatrix Df12_f12 = compute_1rdm_mp2part(spin, nocc1_act, nocc2_act, nvir1, nvir2,
                                                    T2_f12corr[spincase2], T2_f12corr[AlphaBeta]);
        Dmp2_f12corr.accumulate(Df12_f12);
        Df12_f12 = 0;

        //Dmp2f12.print("MP2F12 one-electron density matrix");
        ExEnv::out0() << endl << "Trace of D_MP2_F12corr: " << scprintf("%12.10f", Dmp2_f12corr.trace())<< endl;
//      D[spin].accumulate_subblock(Dmp2f12, 0, norb-1, 0, norb-1, 0, 0);
//      if (debug_ >= DefaultPrintThresholds::allN2)
//      D[spin].print(prepend_spincase(spin,"MP2F12 one-particle density:").c_str());

        // Test code for MP2F12 1rdm
#if 0
        double* T2_mp2f12 = new double[nvir *nvir *nocc_act *nocc_act];
        fill_n(T2_mp2f12, nvir *nvir *nocc_act *nocc_act, 0.0);
        compute_T2abij_mp2f12(spincase2, C_0, C_1, T2_mp2f12);

        double* T2_mp2f12_ab = new double[nvir1 *nvir2 *nocc1_act *nocc2_act];
        fill_n(T2_mp2f12_ab, nvir1 *nvir2 *nocc1_act *nocc2_act, 0.0);
        compute_T2abij_mp2f12(AlphaBeta, C_0, C_1, T2_mp2f12_ab);

        RefSCMatrix Dmp2f12_mp2f12 = compute_1rdm_mp2part(spin, nocc1_act, nocc2_act, nvir1, nvir2,
                                                          T2_mp2f12, T2_mp2f12_ab);
        delete[] T2_mp2f12;
        delete[] T2_mp2f12_ab;

        Dmp2f12_mp2f12.print("Test: MP2F12 one-electron density matrix");
        ExEnv::out0() << endl << spin_label << " Test: trace of D_MP2F12: " << scprintf("%12.10f", Dmp2f12_mp2f12.trace())<< endl;
        Dmp2f12_mp2f12 = 0;
#endif

        RefSCMatrix D_MP2F12 = localkit->matrix(rowdim, coldim);
        D_MP2F12.assign(0.0);
        D_MP2F12.accumulate_subblock(Dmp2_f12corr, 0, norb-1, 0, norb-1, 0, 0);
        Dmp2_f12corr = 0;

        RefSCMatrix opdm_mp2f12 = onepdm_transformed(spin, true, D_MP2F12);
        D_MP2F12 = 0;

        RefSCMatrix opdm_x_mp2f12 = MX_nb_ribs * opdm_mp2f12;
        RefSCMatrix opdm_y_mp2f12 = MY_nb_ribs * opdm_mp2f12;
        RefSCMatrix opdm_z_mp2f12 = MZ_nb_ribs * opdm_mp2f12;
        opdm_mp2f12 = 0;
        const double dx_mp2f12 = opdm_x_mp2f12.trace();
        const double dy_mp2f12 = opdm_y_mp2f12.trace();
        const double dz_mp2f12 = opdm_z_mp2f12.trace();
        opdm_x_mp2f12 = 0;
        opdm_y_mp2f12 = 0;
        opdm_z_mp2f12 = 0;
        ExEnv::out0() << endl << "x y z dipole moments from mp2 f12 coupling: "
                               << scprintf("%12.10f", dx_mp2f12) << "  "
                               << scprintf("%12.10f", dy_mp2f12) << "  "
                               << scprintf("%12.10f", dz_mp2f12) << endl;

        dipoles_mp2[0] = dipoles_mp2[0] + dx_mp2;
        dipoles_mp2[1] = dipoles_mp2[1] + dy_mp2;
        dipoles_mp2[2] = dipoles_mp2[2] + dz_mp2;

        dipoles_f12[0] = dipoles_f12[0] + dx_mp2f12;
        dipoles_f12[1] = dipoles_f12[1] + dy_mp2f12;
        dipoles_f12[2] = dipoles_f12[2] + dz_mp2f12;
    }

//    RefSCMatrix opdm = onepdm_transformed(spin, D[spin]);
//    // opdm.print("CCSD_F12 one-particle density matrix MPQC ordering");
//
//    RefSCMatrix opdm_x = MX_nb_ribs * opdm;
//    RefSCMatrix opdm_y = MY_nb_ribs * opdm;
//    RefSCMatrix opdm_z = MZ_nb_ribs * opdm;
//    opdm = 0;

//    const double dx = opdm_x.trace();
//    const double dy = opdm_y.trace();
//    const double dz = opdm_z.trace();
//    opdm_x = 0;
//    opdm_y = 0;
//    opdm_z = 0;

    // clean
    MX_nb_ribs = 0;
    MY_nb_ribs = 0;
    MZ_nb_ribs = 0;
  }
  // end of loop over nspincase1

  if (nspincases1 == 1) {
    D[Beta] = D[Alpha];
  }

  // clean
  delete[] Dm_i_alpha;
  delete[] Dc_b_alpha;
  delete[] Dcp_bp_alpha_A;
  delete[] Dcp_bp_alpha_Ap;
  delete[] Dap_a_alpha_RR;
  if (nspincases1 == 2) {
    delete[] Dm_i_beta;
    delete[] Dc_b_beta;
    delete[] Dcp_bp_beta_A;
    delete[] Dcp_bp_beta_Ap;
    delete[] Dap_a_beta_RR;
   }
  if (this->r12eval()->ebc() == false
      || this->r12eval()->coupling() == true) {
    delete[] Dap_a_alpha_RT;
    if (nspincases1 == 2)
      delete[] Dap_a_beta_RT;
  }
  if (!r12intermediates_->T2_cc_computed()) {
    for(int s = 0; s < nspincases2; ++s) {
      delete[] T2_mp2[s];
      delete[] T2_f12corr[s];
    }
  }

  // compute the nuclear dipole moment
  Ref<Molecule> molecule = orbs1->basis()->molecule();
  int natom = molecule->natom();
  // test
  //molecule->print();
  double dx_nuc = 0.0, dy_nuc = 0.0, dz_nuc = 0.0;
  for(int i = 0; i < natom; i++) {
    dx_nuc += molecule->r(i,0) * molecule->Z(i);
    dy_nuc += molecule->r(i,1) * molecule->Z(i);
    dz_nuc += molecule->r(i,2) * molecule->Z(i);
  }
  ExEnv::out0() << endl << endl <<  "Nuclear dipole x y z: "
                        << scprintf("%12.10f",dx_nuc) << " "
                        << scprintf("%12.10f",dy_nuc) << " "
                        << scprintf("%12.10f",dz_nuc) << endl;

  if (nspincases1 == 1) {
    dipoles_f12[0] = - dipoles_f12[0] * 2;
    dipoles_f12[1] = - dipoles_f12[1] * 2;
    dipoles_f12[2] = - dipoles_f12[2] * 2;

    dipoles_cabs[0] = - dipoles_cabs[0] * 2;
    dipoles_cabs[1] = - dipoles_cabs[1] * 2;
    dipoles_cabs[2] = - dipoles_cabs[2] * 2;

    if (r12intermediates_->T2_cc_computed()) {
      dipoles_ccsd[0] = - dipoles_ccsd[0] * 2;
      dipoles_ccsd[1] = - dipoles_ccsd[1] * 2;
      dipoles_ccsd[2] = - dipoles_ccsd[2] * 2;

      if (r12intermediates_->Onerdm_relax_computed()) {
        dipoles_or_relax[0] = - dipoles_or_relax[0] * 2;
        dipoles_or_relax[1] = - dipoles_or_relax[1] * 2;
        dipoles_or_relax[2] = - dipoles_or_relax[2] * 2;
      }
    } else{
        dipoles_mp2[0] = - dipoles_mp2[0] * 2;
        dipoles_mp2[1] = - dipoles_mp2[1] * 2;
        dipoles_mp2[2] = - dipoles_mp2[2] * 2;
    }
  } else {
      dipoles_f12[0] = - dipoles_f12[0];
      dipoles_f12[1] = - dipoles_f12[1];
      dipoles_f12[2] = - dipoles_f12[2];

      dipoles_cabs[0] = - dipoles_cabs[0];
      dipoles_cabs[1] = - dipoles_cabs[1];
      dipoles_cabs[2] = - dipoles_cabs[2];

      if (r12intermediates_->T2_cc_computed()) {
        dipoles_ccsd[0] = - dipoles_ccsd[0];
        dipoles_ccsd[1] = - dipoles_ccsd[1];
        dipoles_ccsd[2] = - dipoles_ccsd[2];
        if (r12intermediates_->Onerdm_relax_computed()) {
          dipoles_or_relax[0] = - dipoles_or_relax[0];
          dipoles_or_relax[1] = - dipoles_or_relax[1];
          dipoles_or_relax[2] = - dipoles_or_relax[2];
        }
      } else{
          dipoles_mp2[0] = - dipoles_mp2[0];
          dipoles_mp2[1] = - dipoles_mp2[1];
          dipoles_mp2[2] = - dipoles_mp2[2];
      }
  }

  RefSCVector dipoles = localkit->vector(RefSCDimension(new SCDimension(3)));
  dipoles.assign(0.0);

  dipoles[0] = dx_nuc + dipoles_f12[0] + dipoles_cabs[0];
  dipoles[1] = dy_nuc + dipoles_f12[1] + dipoles_cabs[1];
  dipoles[2] = dz_nuc + dipoles_f12[2] + dipoles_cabs[2];
  if (r12intermediates_->T2_cc_computed()) {
    dipoles[0] = dipoles[0] + dipoles_ccsd[0];
    dipoles[1] = dipoles[1] + dipoles_ccsd[1];
    dipoles[2] = dipoles[2] + dipoles_ccsd[2];
    if (r12intermediates_->Onerdm_relax_computed()) {
        dipoles[0] = dipoles[0] + dipoles_or_relax[0];
        dipoles[1] = dipoles[1] + dipoles_or_relax[1];
        dipoles[2] = dipoles[2] + dipoles_or_relax[2];
    }
  } else {
      dipoles[0] = dipoles[0] + dipoles_mp2[0];
      dipoles[1] = dipoles[1] + dipoles_mp2[1];
      dipoles[2] = dipoles[2] + dipoles_mp2[2];
  }

  dipoles.print("dipole moment: mu(X) mu(Y) mu(Z)");
  dipoles_f12.print("F12 contribution to dipole moment: mu(X) mu(Y) mu(Z)");
  dipoles_cabs.print("CABS Singles contribution to dipole moment: mu(X) mu(Y) mu(Z)");
  if (r12intermediates_->T2_cc_computed()) {
    dipoles_ccsd.print("CCSD contribution (nuclear excluded) to dipole moment: mu(X) mu(Y) mu(Z)");
    dipoles_or_relax.print("CCSD_F12 orbital relaxation contribution to dipole moment: mu(X) mu(Y) mu(Z)");
  } else {
      dipoles_mp2.print("MP2 contribution to dipole moment: mu(X) mu(Y) mu(Z)");
  }

  ExEnv::out0() << endl << "End of the computation for F12 one-particle density" << endl;
  ExEnv::out0() << endl << "**********************************************" << endl;
  return;
}
// end of compute_density_diag

// Transform the density matrix in ribs space (ordered in occ, uocc, cabs)
//                        to the MPQC ordering (ordered in irreps)
RefSCMatrix MP2R12Energy_Diag::onepdm_transformed(const SpinCase1& spin, const bool frozen_core,
                                                  const RefSCMatrix& D) {

  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();

  Ref<OrbitalSpace> occ;
  if (frozen_core == true) {
    occ = r12eval()->occ_act(spin);
  } else {
      occ = r12eval()->occ(spin);
  }
  const Ref<OrbitalSpace>& vir = r12eval()->vir(spin);
//  const Ref<OrbitalSpace>& orbs = r12eval()->orbs(spin);
  const Ref<OrbitalSpace>& cabs = r12world->cabs_space(spin);
  const Ref<OrbitalSpace> ribs = r12world->ribs_space();

  const std::vector<unsigned int>& occpi = occ->block_sizes();
  const std::vector<unsigned int>& virpi = vir->block_sizes();
  const std::vector<unsigned int>& cabspi = cabs->block_sizes();
//  const std::vector<unsigned int>& orbspi = orbs->block_sizes();
//  const std::vector<unsigned int>& ribspi = ribs->block_sizes();

  const int nirreps = occ->nblocks();
  const int nocc = occ->rank();
  const int nvir = vir->rank();
  const int ncabs = cabs->rank();

//  const int norbs = orbs->rank();
//  const int nribs = ribs->rank();
  const int norbs = nocc + nvir;
  const int nribs = norbs + ncabs;

//  ExEnv::out0() << endl << "nocc: " << nocc << endl
//                        << "norbs: " << norbs << endl
//                        << "ncabs: " << ncabs << endl;

  std::vector<unsigned int> occoff(nirreps);
  std::vector<unsigned int> viroff(nirreps);
  std::vector<unsigned int> orbsoff(nirreps);
  std::vector<unsigned int> cabsoff(nirreps);
  std::vector<unsigned int> ribsoff(nirreps);

  occoff[0] = 0;
  viroff[0] = 0;
  orbsoff[0] = 0;
  cabsoff[0] = 0;
  ribsoff[0] = 0;

//  ExEnv::out0() << endl << "irrep  occ_actpi virpi cabspi:" << endl
//                        << 0 << "  " << occ_actpi[0] << "  " << virpi[0] << "  " << cabspi[0] <<  endl;
  for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
    occoff[irrep] = occoff[irrep-1] + occpi[irrep-1];
    viroff[irrep] = viroff[irrep-1] + virpi[irrep-1];
    cabsoff[irrep] = cabsoff[irrep-1] + cabspi[irrep-1];
//    orbsoff[irrep] = orbsoff[irrep-1] + orbspi[irrep-1];
//    ribsoff[irrep] = ribsoff[irrep-1] + ribspi[irrep-1];
    ribsoff[irrep] = occoff[irrep] + viroff[irrep] + cabsoff[irrep];
//    ExEnv::out0() << irrep << "   " << occ_actpi[irrep] << "  " << virpi[irrep] << "  " << cabspi[irrep] << endl;
  }

  const RefSCDimension opdm_dim(new SCDimension(nribs));
  const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix opdm_mat = localkit->matrix(opdm_dim, opdm_dim);
  opdm_mat.assign(0.0);

  for (unsigned int h = 0; h < nirreps; ++h) {
    // offset for D matrix in occ, uocc, cabs ordering
    const unsigned int i_offset = occoff[h];
    const unsigned int j_offset = i_offset;
    const unsigned int a_offset = nocc + viroff[h];
    const unsigned int b_offset = a_offset;
    const unsigned int ap_offset = norbs + cabsoff[h];
    const unsigned int bp_offset = ap_offset;

    const int mpqc_ij_offset = ribsoff[h];
    const int mpqc_ab_offset = ribsoff[h] + occpi[h];
    const int mpqc_apbp_offset = ribsoff[h] + occpi[h] + virpi[h];

    for (int i = 0; i < occpi[h]; ++i) {
      for (int j = 0; j <= i; ++j) {
        const double Dij = D.get_element(i+i_offset,j+j_offset);
        opdm_mat.set_element(i+mpqc_ij_offset, j+mpqc_ij_offset, Dij);
      }
    }
    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsdf12 opdm ij").c_str());

    for (int a = 0; a < virpi[h]; ++a) {
      for (int b = 0; b <= a; ++b) {
        const double Dab = D.get_element(a+a_offset,b+b_offset);
        opdm_mat.set_element(a+mpqc_ab_offset, b+mpqc_ab_offset, Dab);
      }
    }
    //opdm_mat.print(prepend_spincase(spin,"MPQC ccsdf12 opdm ab").c_str());

    for (int a = 0; a < virpi[h]; ++a) {
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
      for (int i = 0; i < occpi[h]; ++i) {
        const double Dapi = D.get_element(ap+ap_offset,i+i_offset);
        opdm_mat.set_element(ap+mpqc_apbp_offset, i+mpqc_ij_offset, Dapi);
        //ExEnv::out0() << endl << h << " " << ap+mpqc_apbp_offset << " " << i+mpqc_ij_offset << " Da'i = " << Dapi << endl;
      }
    }

    for (int ap = 0; ap < cabspi[h]; ++ap) {
      for (int a = 0; a < virpi[h]; ++a) {
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

//  const Ref<OrbitalSpace>& occ = r12eval()->occ_act(spin);
  const Ref<OrbitalSpace>& occ = r12eval()->occ(spin);
  const Ref<OrbitalSpace>& vir = r12eval()->vir(spin);
  const Ref<OrbitalSpace>& orbs = r12eval()->orbs(spin);

  const std::vector<unsigned int>& occpi = occ->block_sizes();
  const std::vector<unsigned int>& uoccpi = vir->block_sizes();
//  const std::vector<unsigned int>& orbspi = orbs->block_sizes();

  const int nirreps = occ->nblocks();
  const int nocc = occ->rank();
  const int nvir = vir->rank();
  const int norbs = nocc + nvir;

  std::vector<unsigned int> occoff(nirreps);
  std::vector<unsigned int> uoccoff(nirreps);
  std::vector<unsigned int> orbsoff(nirreps);

  occoff[0] = 0;
  uoccoff[0] = 0;
  orbsoff[0] = 0;

  for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
    occoff[irrep] = occoff[irrep-1] + occpi[irrep-1];
    uoccoff[irrep] = uoccoff[irrep-1] + uoccpi[irrep-1];
    orbsoff[irrep] = occoff[irrep] + uoccoff[irrep];
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
