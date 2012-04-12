/*
 * compute_density_diag.cc
 *
 *  Created on: Nov 8, 2011
 *      Author: jinmei
 */

#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <util/misc/print.h>
#include <math/scmat/blas.h>

using namespace std;
using namespace sc;

// case1: R^idx1 idx3_ij * R^ij_idx2 idx3, eg: R^ba'_ij R^ij_ca'
// case2: R^idx3 idx1_ij * R^ij_idx2 idx3, eg: R^a'b_ij R^ij_ca'
// case3: R^idx1 idx3_ij * R^ij_idx3 idx2, eg: R^ba'_ij R^ij_a'c
// case4: R^idx3 idx1_ij * R^ij_idx3 idx2, eg: R^a'b_ij R^ij_a'c
enum idx_cases{RR13_23, RR31_23, RR13_32, RR31_32};

// orbs1_orbs2_orbs3: D^orbs1_orbs2 sum over orbs2
enum orbitals_cases{vir_vir_cabs, cabs_cabs_vir, cabs_cabs_cabs, cabs_vir_cabs};

namespace {
  // print out integrals
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

//  // antisymmetrize blk_ij
//  void antisymmetrize(double* blk_ij, const double* blk_ji, const int size)
//  {
//    double* iter_blk_ij = blk_ij;
//    double* iter_blk_ji = blk_ji;
//    for (int i = 0; i < size; ++i) {
//      for (int j = 0; j < size; ++j, ++iter_blk_ij, ++iter_blk_ji) {
//        *iter_blk_ij -= *iter_blk_ji;
//      }
//    }
//  } // end

  // print tension with 1 index
  void print_intermediate(const string& spinlabel, const string& label,
                          const double* const array, const int size_idx)
  {
     ExEnv::out0() << indent << spinlabel << " " << label << endl;

     const int col_print = 7;

     // print the column label
     const int size_col = (size_idx > col_print? col_print : size_idx);
     for (int i = 1; i <= size_col; ++i) {
         ExEnv::out0() << indent << scprintf("%12d", i) << "  ";
     }
     ExEnv::out0() << endl << indent << 1 << " " ;

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

//  // print tension with 2 indices
//  void print_intermediate(const string& spinlabel, const string& label,
//                          const double* const array,
//                          const int size_idx1, const int size_idx2)
//  {
//     ExEnv::out0() << indent << spinlabel << " " << label << endl;
//
//     const int size_col1 = (size_idx2 > 5? 5 : size_idx2);
//
//     const double* iter = array;
//     for (int idx1 = 1 ; idx1 <= size_idx1; ++idx1) {
//       // print the column label
//       for (int i = 1; i <= size_col1; ++i) {
//         ExEnv::out0() << indent << scprintf("%12d", i) << "  ";
//       }
//       ExEnv::out0() << endl << indent << idx1 << " ";
//
//       for (int idx2 = 1; idx2 <= size_idx2; ++idx2, ++iter) {
//         ExEnv::out0() << indent << scprintf("%12.10f", *iter) << "  ";
//
//         // print the column label for the next 5 elements
//         const int left_col = idx2 % 5;
//         if (left_col == 0) {
//           ExEnv::out0() << endl;
//
//           const int size_col2 = ((idx2 + 5) > size_idx2? size_idx2 : idx2 + 5);
//           for (int i = idx2 + 1; i <= size_col2; ++i) {
//             ExEnv::out0() << indent << scprintf("%12d", i) << "  ";
//           }
//           ExEnv::out0() << endl << indent << idx1 << " ";
//         }
//
//       } // end of looping idx2
//       ExEnv::out0() << endl << endl;
//     } // end of looping idx1
//     ExEnv::out0() << endl;
//   }

//  // test propose
//  // print tension with 2 indices without change line
//  void print_intermediate(const string& spinlabel, const string& label,
//                          const double* const array,
//                          const int size_idx1, const int size_idx2)
//  {
//     ExEnv::out0() << indent << spinlabel << " " << label << endl;
//
//     const double* iter = array;
//     for (int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
//       ExEnv::out0() << endl << indent << idx1 << " ";
//
//       for (int idx2 = 0; idx2 < size_idx2; ++idx2) {
//         ExEnv::out0() << indent << scprintf("%12.10f", *iter) << "  ";
//         ++iter;
//       }
//       ExEnv::out0() << endl;
//     }
//     ExEnv::out0() << endl;
//   } // end

  // print tension with 2 indices
  void print_intermediate(const string& spinlabel, const string& label,
                          const double* const array,
                          const int size_idx1, const int size_idx2)
  {
     ExEnv::out0() << indent << spinlabel << " " << label << endl;

     const int col_print = 7;
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
        ExEnv::out0() << "***";
        f12_ints->release_pair_block(*f12_blk_idx1, *f12_blk_idx2, f12_idx);
      }
      ExEnv::out0() << endl;
    }
  }

  void print_RT2(const string& spinlabel, const string& label, const double* const RT2,
                 const int nap, const int nb, const int na)
  {
    ExEnv::out0() << indent << spinlabel << " " << label << endl;

    const double* iter = RT2;
    for (int ap = 0 ; ap < nap; ++ap) {
      ExEnv::out0() << indent << "a' = " << ap << endl;

      for (int b = 0; b < nb; ++b) {
        for (int a = 0; a < na; ++a) {
          ExEnv::out0() << indent << scprintf("%12.10f", *iter) << "  ";
          ++iter;
        }
        ExEnv::out0() << endl;
      }
      ExEnv::out0() << endl;
    }
  }
  // print tension with 3 indices
  void print_intermediate(const SpinCase1 spin, const string& label, const double* const array,
                          const int nidx1, const int nidx2, const int nidx3)
  {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    ExEnv::out0() << indent << spinletter << " " << label << endl;
    const double* iter = array;
    for (int i1 = 0 ; i1 < nidx1; ++i1) {
      for (int i2 = 0; i2 < nidx2; ++i2) {
        for (int i3 = 0; i3 < nidx3; ++i3, ++iter) {
                ExEnv::out0() << indent << scprintf("%12.10f", *iter) << "  ";
        }
        ExEnv::out0() << endl;
      }
      ExEnv::out0() << endl;
    }
  }

  void get_conjugation(const double* array_ij, double* array_ji,
                       const int ni, const int nj)
  {
    const double* iter_ij = array_ij;

    for (int i = 0 ; i < ni; ++i) {

      double* iter_ji = array_ji + i;
      for (int j = 0; j < nj; ++j) {
        *iter_ji = *iter_ij;

        ++iter_ij;
        iter_ji += ni;
      }
    }
  }


  // compute (f12f12)^b1b2_k1k2 sum over index 3:
  // which is j for (f12f12)^ij_ij and its other permutaitions
  // and index 1 and 2 are i
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
  } // end of function: compute_F12F12_sum_idx3

  // compute R * R contracted over 3 indices
  // here they are: i, j, idx3
  void compute_RR2_sum_3idx(const int RRb1b2_k1k2, const unsigned int f12_idx,
                           Ref<DistArray4>& f12_ints1, Ref<DistArray4>& f12_ints2,
                           double* f12f12_array)
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
//          ExEnv::out0() << indent << scprintf("%12.10f", f12f12_sum_idx3) << " ";
          *iter_array = f12f12_sum_idx3;
          ++iter_array;
        }
//        ExEnv::out0() << endl;
      }
  } // end of compute_RR2_sum_3dix

  // compute R^b1b2_ab * R^ab_k1k2 contracted over a, b, and j
  // b1, b2, k1, k2 can be i or j (only two indices)
  // external index: j
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
           //        ExEnv::out0() << indent << scprintf("%12.10f", f12f12_sum_ij) << " ";

        f12f12_sum_idx3 += f12f12_sum_idx;

        f12_ints1->release_pair_block((*f12_blk1_idx1), (*f12_blk1_idx2), f12_idx);
        f12_ints2->release_pair_block((*f12_blk2_idx1), (*f12_blk2_idx2), f12_idx);
      }
      *iter_array = f12f12_sum_idx3;
      ++iter_array;
    }
  //    ExEnv::out0() << endl;
  } // end of compute_RR1_sum_3dix

  void compute_Dii_spin(const SpinCase1 spin, const double C_0, const double C_1,
                        int nocc_alpha, int nocc_beta,
                        const double* RRij_ij, const double* RRji_ij,
                        const double* RRi1i2_i1i2, const double* RRi1i2_i2i1,
                        const double* RRi2i1_i1i2, const double* RRi2i1_i2i1,
                        RefSCMatrix& Dii)
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

      Dii.set_element(i,i,dii_12);
    }

  } // end of function: compute_Dii_spin

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

  const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
  const Ref<OrbitalSpace>& vir1 = r12eval()->vir(spin1);
  const Ref<OrbitalSpace>& orbs1 = r12eval()->orbs(spin1);
  const Ref<OrbitalSpace>& cabs1 = r12world->cabs_space(spin1);
  const Ref<OrbitalSpace>& occ1 = r12eval()->occ(spin1);

  const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
  const Ref<OrbitalSpace>& vir2 = r12eval()->vir(spin2);
  const Ref<OrbitalSpace>& orbs2 = r12eval()->orbs(spin2);
  const Ref<OrbitalSpace>& cabs2 = r12world->cabs_space(spin2);
  const Ref<OrbitalSpace>& occ2 = r12eval()->occ(spin2);

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
} // end of obtain_orbitals

// activate f12f12_ints for computing X
void MP2R12Energy_Diag::activate_ints_f12f12(Ref<TwoBodyFourCenterMOIntsRuntime>& moints4_rtime, const int b1b2_k1k2,
                                             const vector< Ref<OrbitalSpace> >& v_orbs1,
                                             const vector< Ref<OrbitalSpace> >& v_orbs2,
                                             const string&  descr_f12f12_key, Ref<DistArray4>& f12f12_ints)
{
  // obtain orbitals
  const Ref<OrbitalSpace> occ1_act = v_orbs1[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2[0];

  // indices for (f12f12)^b1b2_k1k2
  const Ref<OrbitalSpace>* bra1 = &occ1_act;
  const Ref<OrbitalSpace>* bra2 = &occ2_act;
  const Ref<OrbitalSpace>* ket1 = &occ1_act;
  const Ref<OrbitalSpace>* ket2 = &occ2_act;

  switch (b1b2_k1k2) {
  case ij_ij:
    // default value
  break;

  case ij_ji:
    ket1 = &occ2_act;
    ket2 = &occ1_act;
  break;

  case ji_ij:
    bra1 = &occ2_act;
    bra2 = &occ1_act;
  break;

  case ji_ji:
    bra1 = &occ2_act;
    bra2 = &occ1_act;
    ket1 = &occ2_act;
    ket2 = &occ1_act;
  break;

  default:
    ExEnv::out0() << "There is no such index";
    break;
  }

  // (f12f12)^b1b2_k1k2
  activate_ints((*bra1)->id(), (*bra2)->id(), (*ket1)->id(), (*ket2)->id(),
                descr_f12f12_key, moints4_rtime, f12f12_ints);
} // end of function: activate_ints_f12f12

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
} // end of function: activate_ints_X_f12

// compute AlphaAlpha/BetaBeta:
// R^ij_ab R^ab_ij and R^ji_ab R^ab_ij
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
  activate_ints_f12f12(moints4_rtime, ij_ij,
                       v_orbs1, v_orbs2,
                       descr_f12f12_key, ijij_f12f12_ints);
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

} // end of function: compute_RRii_ii

// compute the openshell AlphaBeta:
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
  Ref<DistArray4> ijij_f12f12_ints;
  Ref<DistArray4> jiij_f12f12_ints;
  Ref<DistArray4> ijji_f12f12_ints;
  Ref<DistArray4> jiji_f12f12_ints;
  activate_ints_f12f12(moints4_rtime, ij_ij,
                       v_orbs1, v_orbs2,
                       descr_f12f12_key, ijij_f12f12_ints);
  activate_ints_f12f12(moints4_rtime, ji_ij,
                       v_orbs1, v_orbs2,
                       descr_f12f12_key, jiij_f12f12_ints);
  activate_ints_f12f12(moints4_rtime, ij_ji,
                       v_orbs1, v_orbs2,
                       descr_f12f12_key, ijji_f12f12_ints);
  activate_ints_f12f12(moints4_rtime, ji_ji,
                       v_orbs1, v_orbs2,
                       descr_f12f12_key, jiji_f12f12_ints);
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
} // end of function: compute_Rii_ii

// compute MP2 T2 amplitude (antisymmetrized), stored in array (a,b,i,j)
//  T^i(spin1)j(spin2)_a(spin1)b(spin2) 4-dimension matrix
//= g^ij_ab / (e_i + e_j - e_a - e_b)
void MP2R12Energy_Diag::compute_T2_mp2(const SpinCase2 spincase,
                                       const vector< Ref<OrbitalSpace> >& v_orbs1,
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
  activate_ints(vir1->id(), vir2->id(), occ1_act->id(),
                occ2_act->id(), descr_f12_key, moints4_rtime,
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

} // end of function: compute_T2_mp2

// compute F12T2(a',b,a) = R^ij_a'b * T^ab_ij   (alpha)
//                      or R^ij_ba' * T^ba_ij   (beta)
// contracted indices: ij
void MP2R12Energy_Diag::compute_F12T2_ij(const SpinCase1 spin, const int na,
                                         const unsigned int f12_idx,
                                         Ref<DistArray4>& f12_ints, const double* T2_ints,
                                         double* F12T2)
{
  // get the sizes of the external indices: a'b
  //                    & internal indices: ij
  int nap = f12_ints->ni();
  int nb = f12_ints->nj();
  const blasint nij =  f12_ints->nx() * f12_ints->ny();

  // get the right block of f12_ints: R^ij_a'b or R^ij_ba'
  int ap = 0;
  int b = 0;
  int* f12_blk_idx1 = &ap;
  int* f12_blk_idx2 = &b;

  // Get the right start iterator for T^ab_ij or T^ba_ij
  int iter_T2ab_start = 0;
  int iter_T2ba_start = 0;
  int* iter_T2_start = &iter_T2ab_start;

  // Get the right offset for T^ab_ij or T^ba_ij in the inner loop
  int offset_T2ab = nb * nij;
  int offset_T2ba = nij;
  int* offset_T2 = &offset_T2ab;

  if (spin == Beta) {
    nap = f12_ints->nj();
    nb = f12_ints->ni();

    f12_blk_idx1 = &b;
    f12_blk_idx2 = &ap;

    iter_T2_start = &iter_T2ba_start;
    offset_T2 = &offset_T2ba;
  }

  // indices of F12T2: : a', b, and a
  // which means its indices increases one each loop
  double* iter_F12T2 = F12T2;

#if 0
  // test code: print f12_ij_blk
  const SpinCase2 spincase = static_cast<SpinCase2>(spin + 1);
  string spinletters = to_string(spincase);
  const string f12_ij_label = (spin == Alpha ? "R^ij_a'b" : "R^ij_ba'");
  ExEnv::out0() << spinletters << " "<< f12_ij_label << " (ap' by b)"
                << " nap = " << nap << " nb = " << nb << "nij = " << nij << endl;
#endif

  const blasint one = 1; // for F77_DDOT
  for(ap = 0; ap < nap; ++ap) {
    for (b = 0; b < nb; ++b) {
      const double* f12_ij_blk = f12_ints->retrieve_pair_block(*f12_blk_idx1, *f12_blk_idx2, f12_idx);

      iter_T2ab_start = b * nij;
      iter_T2ba_start = b * na * nij;
      const double* iter_T2 = T2_ints + *iter_T2_start;
      for(int a = 0; a < na; ++a) {
        const double f12t2 = F77_DDOT(&nij, f12_ij_blk, &one, iter_T2, &one);
        *iter_F12T2 = f12t2;

        ++iter_F12T2;
        iter_T2 += *offset_T2;
      }

//      // test code: print f12_ij_blk
//      const double* iter = f12_ij_blk;
//      for (int i=0 ; i < nij; ++i) {
//        ExEnv::out0() << indent << scprintf("%12.10f", *iter) << " ";
//        ++iter;
//      }
//       ExEnv::out0() << "   ";

      f12_ints->release_pair_block(*f12_blk_idx1, *f12_blk_idx2, f12_idx);
    }
//    ExEnv::out0() << endl;
  }

} // end of compute_F12T2_ij

// compute F12T2(a',b,a) = R^ji_a'b * T^ab_ij   (alpha)
//                      or R^ji_ba' * T^ba_ij   (beta)
// contracted indices: i & j
void MP2R12Energy_Diag::compute_F12T2_ji(const SpinCase1 spin, const int na,
                                         const unsigned int f12_idx,
                                         Ref<DistArray4>& f12_ints, const double* T2_ints,
                                         double* F12T2)
{
  // get the sizes of the external indices: a'b
  //                    & internal indices: ij
  int nap = f12_ints->ni();
  int nb = f12_ints->nj();
  const int ni = f12_ints->ny();
  const int nj = f12_ints->nx();
  const blasint nij = ni * nj;

  // get the right block of f12_ints: R^ij_a'b or R^ij_ba'
  int ap = 0;
  int b = 0;
  int* f12_blk_idx1 = &ap;
  int* f12_blk_idx2 = &b;

  // Get the right start iterator for T^ab_ij or T^ba_ij
  int iter_T2ab_start = 0;
  int iter_T2ba_start = 0;
  int* iter_T2_start = &iter_T2ab_start;

  // Get the right offset for T^ab_ij or T^ba_ij in the inner loop
  int offset_T2ab = nb * nij;
  int offset_T2ba = nij;
  int* offset_T2 = &offset_T2ab;

  if (spin == Beta) {
    nap = f12_ints->nj();
    nb = f12_ints->ni();

    f12_blk_idx1 = &b;
    f12_blk_idx2 = &ap;

    iter_T2_start = &iter_T2ba_start;
    offset_T2 = &offset_T2ba;
  }

  // indices of F12T2: : a', b, and a
  // which means its indices increases one each loop
  double* iter_F12T2 = F12T2;

#if 0
  // test code: print f12_ij_blk
  const SpinCase2 spincase = static_cast<SpinCase2>(spin + 1);
  string spinletters = to_string(spincase);
  const string f12_ji_label = (spin == Alpha ? "R^ji_a'b" : "R^ji_ba'");
  ExEnv::out0() << spinletters << " "<< f12_ji_label << " (ap' by b)"
                << " nap = " << nap << " nb = " << nb << "nij = " << nij << endl;
#endif

  const blasint one = 1; // for F77_DDOT
  for(ap = 0; ap < nap; ++ap) {
    for (b = 0; b < nb; ++b) {
      const double* f12_ij_blk = f12_ints->retrieve_pair_block(*f12_blk_idx1, *f12_blk_idx2, f12_idx);
      double* f12_ji_blk = new double[nij];
      get_conjugation(f12_ij_blk, f12_ji_blk, ni, nj);

      iter_T2ab_start = b * nij;
      iter_T2ba_start = b * na * nij;
      const double* iter_T2 = T2_ints + *iter_T2_start;
      for(int a = 0; a < na; ++a) {
        const double f12t2 = F77_DDOT(&nij, f12_ji_blk, &one, iter_T2, &one);
        *iter_F12T2 = f12t2;

        ++iter_F12T2;
        iter_T2 += *offset_T2;
      }

//      // test code: print f12_ij_blk
//      const double* iter = f12_ji_blk;
//      for (int i = 0 ; i < nij; ++i) {
//        ExEnv::out0() << indent << scprintf("%12.10f", *iter) << " ";
//        ++iter;
//      }
//       ExEnv::out0() << "   ";

      delete[] f12_ji_blk;
      f12_ints->release_pair_block(*f12_blk_idx1, *f12_blk_idx2, f12_idx);
    }
//    ExEnv::out0() << endl;
  }
} // end of compute_F12T2_ji

// compute R^ij_a'b * T2^ab_ij & R^ji_a'b * T2^ab_ij
// T2: mp2 amplitude
void MP2R12Energy_Diag::compute_RT2_mp2(const SpinCase1 spincase1, const SpinCase2 spincase2,
                                        const vector< Ref<OrbitalSpace> >& v_orbs1,
                                        const vector< Ref<OrbitalSpace> >& v_orbs2,
                                        const double* T2ab_ij, double* RT2ij, double* RT2ji)
{
  // get moints4_rtime, descr_f12_key, and f12_idx
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);

  // obtain orbitals
  const Ref<OrbitalSpace>& occ1_act = v_orbs1[0];
  const Ref<OrbitalSpace>& occ2_act = v_orbs2[0];
  const Ref<OrbitalSpace>& vir1 = v_orbs1[1];
  const Ref<OrbitalSpace>& vir2 = v_orbs2[1];
  const Ref<OrbitalSpace>& cabs1 = v_orbs1[3];
  const Ref<OrbitalSpace>& cabs2 = v_orbs2[3];
  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();
  const int nvir1 = vir1->rank();
  const int nvir2 = vir2->rank();
  const int nocc12 = nocc1_act * nocc2_act;
  const int nvir12 = nvir1 * nvir2;

  const int nspincase2 = (r12eval()->spin_polarized() ? 3 : 2);

  Ref<DistArray4> i1i2_ints = NULL;    // ap1a2i1i2 or a1ap2i1i2
  Ref<DistArray4> i2i1_ints = NULL;    // ap1a2i2i1 or a1ap2i2i1
  if (spincase1 == Alpha){
    // AlphaAlpha, AlphaBeta for spin = Alpha
    // R^IJ_A'B  (A':alpha)
    activate_ints(cabs1->id(), vir2->id(), occ1_act->id(),
                  occ2_act->id(), descr_f12_key, moints4_rtime,
                  i1i2_ints);
    // R^JI_A'B  (A':alpha)
    if (nspincase2 == 3) {
      activate_ints(cabs1->id(), vir2->id(), occ2_act->id(),
                    occ1_act->id(), descr_f12_key, moints4_rtime,
                    i2i1_ints);
    } else {
        i2i1_ints = i1i2_ints;
    }

  } else {
      // BetaBeta, AlphaBeta for spin = Beta
      // R^IJ_BA' (A':beta)
      activate_ints(vir1->id(), cabs2->id(), occ1_act->id(),
                    occ2_act->id(), descr_f12_key, moints4_rtime,
                    i1i2_ints);
      //  R^JI_BA' (A':beta)
      if (nspincase2 == 3) {
        activate_ints(vir1->id(), cabs2->id(), occ2_act->id(),
                      occ1_act->id(), descr_f12_key, moints4_rtime,
                      i2i1_ints);
      } else {
          i2i1_ints = i1i2_ints;
      }
  }
  const int na = (spincase1 == Alpha? nvir1 : nvir2);

#if 1
//  Ref<DistArray4> i1i2i1i2_ints = NULL;
//  activate_ints(occ1_act->id(), occ1_act->id(), occ1_act->id(),
//                occ2_act->id(), descr_f12_key, moints4_rtime,
//                i1i2i1i2_ints);
//  ExEnv::out0() << " i1i2i1i2_ints:" << endl;
//  for(int i1 = 0; i1 < nocc1_act; ++i1) {
//    for (int i2 = 0; i2 < nocc2_act; ++i2) {
//      const double* f12_blk = i1i2i1i2_ints->retrieve_pair_block(i1, i2, f12_idx);
//
//      const double* iter = f12_blk;
//      for (int i=0 ; i < nocc12; ++i) {
//        ExEnv::out0() << indent << scprintf("%12.10f", *iter) << " ";
//      }
//      ExEnv::out0() << "   ";
//      i1i2i1i2_ints->release_pair_block(i1, i2, f12_idx);
//    }
//    ExEnv::out0() << endl;
//  }

  // test code: print f12_ij_blk
  const int s = spincase1;
  const SpinCase2 spincase = static_cast<SpinCase2>(s+1);
  string spinletters = to_string(spincase);
  const string f12_ij_label = (spincase1 == Alpha ?
                               "R^ij_a'b" : "R^ij_ba'");
  const string f12_ji_label = (spincase1 == Alpha ?
                               "R^ji_a'b" : "R^ji_ba'");

  const int ncabs1 = cabs1->rank();

  print_f12_ints(f12_ij_label, spincase1,
                 i1i2_ints, f12_idx,
                 ncabs1, nvir2, nocc1_act, nocc2_act);

//    // R^IJ_A'B  (A':alpha)
//    Ref<DistArray4> i1i2ap1a2_ints = NULL;
//    activate_ints(occ1_act->id(), occ2_act->id(), cabs1->id(),
//                  vir2->id(),descr_f12_key, moints4_rtime,
//                  i1i2ap1a2_ints);
//    ExEnv::out0() << "i1i2ap1a2_ints: " << endl;
//    for(int i1 = 0; i1 < nocc1_act; ++i1) {
//      for (int i2 = 0; i2 < nocc2_act; ++i2) {
//        ExEnv::out0() << i1 << " " << i2 << endl;
//        const double* f12_blk = i1i2ap1a2_ints->retrieve_pair_block(i1, i2, f12_idx);
//
//        const double* iter = f12_blk;
//        for (int ap = 0 ; ap < ncabs1; ++ap) {
//            for (int b = 0 ; b < nvir2; ++b) {
//              ExEnv::out0() << indent << scprintf("%12.10f", *iter) << " ";
//            }
//            ExEnv::out0() << endl;
//        }
//        i1i2ap1a2_ints->release_pair_block(i1, i2, f12_idx);
//        ExEnv::out0() << endl;
//      }
//    }

//    ExEnv::out0() << spinletters << " "<< f12_ji_label << " (b by a')" << endl;
//    for (int b = 0; b < nvir2; ++b) {
//      for(int ap = 0; ap < ncabs1; ++ap) {
//        const double* f12_blk = i2i1_ints->retrieve_pair_block(ap, b, f12_idx);
//
//        const double* iter = f12_blk;
//        for (int i=0 ; i < nocc12; ++i) {
//          ExEnv::out0() << indent << scprintf("%12.10f", *iter) << " ";
//
//          ++iter;
//        }
//        ExEnv::out0() << "   ";
//        i2i1_ints->release_pair_block(ap, b, f12_idx);
//      }
//      ExEnv::out0() << endl;
//    }

  print_f12_ints(f12_ji_label, spincase1,
                 i2i1_ints, f12_idx,
                 ncabs1, nvir2, nocc2_act, nocc1_act);
  #endif

  // R^ij_a'b * T^ab_ij (a':alpha) or  R^ij_ba' * T^bb_ij (a':beta)
  compute_F12T2_ij(spincase1, na, f12_idx, i1i2_ints, T2ab_ij, RT2ij);
  // R^ji_a'b * T^ab_ij (a':alpha)  or R^ji_ba' * T^ba_ij (a':beta)
  compute_F12T2_ji(spincase1, na, f12_idx, i2i1_ints, T2ab_ij, RT2ji);

  i1i2_ints->deactivate();
  if (nspincase2 == 3) i2i1_ints->deactivate();
} // end of function: compute_RT2_mp2


void MP2R12Energy_Diag::compute_Dii(const int nspincases1, const int nspincases2,
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
      ExEnv::out0() << endl << spinletters << " D^i_i" << endl
                    << "number of occupied orbital: " << nocc_act << endl;

    // skip spincase if no electron of this kind
    if (nocc_act == 0)
      continue;
    const int nocc12 = nocc_act * nocc_act;

    // initialize D^i_i[s]
    RefSCDimension rowdim_occ = new SCDimension(nocc_act);
    RefSCDimension coldim_occ = new SCDimension(nocc_act);
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    Dii_[s] = localkit->matrix(rowdim_occ, coldim_occ);
    Dii_[s].assign(0.0);

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
          ExEnv::out0() << endl << "AlphaBeta D^i_i"
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
                     RRi1i2_i1i2, RRi1i2_i2i1, RRi2i1_i1i2, RRi2i1_i2i1,
                     Dii_[s]);

     delete[] RRij_ij;
     delete[] RRji_ij;
  //    // Assign m to a subblock of this.
  //    virtual void assign_subblock(SCMatrix *m, int, int, int, int, int=0, int=0) =0;
  //    D[s].assign_subblock(Dii[s],);
  } // end of spincase1 loop

  if (nspincases2 == 3) {
    delete[] RRi1i2_i1i2;
    delete[] RRi1i2_i2i1;
    delete[] RRi2i1_i1i2;
    delete[] RRi2i1_i2i1;
  }
} // end of computation of Di_i

// compute R^ab_b1b2 R^k1k2 _ab (a, b represent complete virtual orbitals)
// sum over a, b, and other index which would be j here and labeled as idx3
void MP2R12Energy_Diag::compute_RR_sum_abj(const int RRb1b2_k1k2, const int nocc_act,
                                           const int f12f12_idx, const int f12_idx,
                                           Ref<DistArray4>& f12f12_ints,
                                           vector< Ref<DistArray4> >& v_f12_ints1, vector< Ref<DistArray4> >& v_f12_ints2,
                                           double* RR_result)
{
  // R^ab_b1b2 R^k1k2 _ab = (f12f12)^b1b2_k1k2
  //                       - f12^b1b2_pq f12^pq_k1k2
  //                       - f12^b1b2_pq_ma' f12^ma'_k1k2
  //                       - f12^b1b2_pq_a'm f12^a'm_k1k2

  // add (f12f12)^b1b2_k1k2
  double* RR_part1 = new double[nocc_act];
  fill_n(RR_part1, nocc_act, 0.0);
  compute_F12F12_sum_idx3(RRb1b2_k1k2,
                          f12f12_idx, f12f12_ints,
                          RR_part1);

  double* iter_RR_result = RR_result;
  const double* iter_RR_part1 = RR_part1;
  for (int i = 0; i != nocc_act; ++i) {
    *iter_RR_result += *iter_RR_part1;

    ++iter_RR_result;
    ++iter_RR_part1;
  }
  delete[] RR_part1;

  // test
  //print_intermediate("", "(f12f12) part", RR_result, nocc_act);

  if (v_f12_ints1.size() != 3)
    ExEnv::out0() << "Then number of integrals in computing R^ab_b1b2 R^k1k2 _ab are wrong"
                  << endl;

  // subtract f12^b1b2_pq f12^pq_k1k2, f12^b1b2_pq_ma' f12^ma'_k1k2,
  // and f12^b1b2_pq_a'm f12^a'm_k1k2
  const int nocc12 = nocc_act * nocc_act;
  double* RR_part2 = new double[nocc12];
  for (int i = 0; i != 3; ++i) {
    fill_n(RR_part2, nocc12, 0.0);

    Ref<DistArray4> f12_ints1 = v_f12_ints1[i];
    Ref<DistArray4> f12_ints2 = v_f12_ints2[i];
    compute_RR1_sum_3idx(RRb1b2_k1k2, f12_idx,
                        f12_ints1, f12_ints2,
                        RR_part2);

    iter_RR_result = RR_result;
    const double* iter_RR_part2 = RR_part2;
    for (int i = 0; i != nocc_act; ++i) {
      *iter_RR_result -= *iter_RR_part2;

      ++iter_RR_result;
      ++iter_RR_part2;
    }
    // test
    //print_intermediate("", "(f12) part", RR_result, nocc_act);
  }
  delete[] RR_part2;

} // end of compute_RR_sum_abj


void MP2R12Energy_Diag::compute_Dii_2(const int nspincases1, const int nspincases2,
                                    const int nocc_alpha, const int nocc_beta,
                                    const double C_0, const double C_1)
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

  // prelimiaries for AlphaBeta parts
  // activate integrals for R^i1i2_a1a2 R^a1a2_i1i2, R^i2i1_a1a2 R^a1a2_i1i2
  //                        R^i1i2_a1a2 R^a1a2_i2i1, R^i2i2_a1a2 R^a1a2_i2i1
  // i: occupied orbital, a: complete virtual orbitals
  // 1: alpha orbtial, 2: beta orbital

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
    // obtain occ, vir, orbs, and cabs orbitals
    vector< Ref<OrbitalSpace> > v_orbs1;
    vector< Ref<OrbitalSpace> > v_orbs2;
    obtain_orbitals(AlphaBeta, v_orbs1, v_orbs2);

    const int nocc1_act = v_orbs1[0]->rank();
    const int nocc2_act = v_orbs2[0]->rank();
    const int nocc12 = nocc1_act* nocc2_act;
    // test
      ExEnv::out0() << endl << "Number of alpha occupied orbital: " << nocc1_act
                    << endl << "Number of beta occupied orbital: " << nocc2_act << endl;

    // activate (f12f12) ints
     activate_ints_f12f12(moints4_rtime, ij_ij,
                          v_orbs1, v_orbs2,
                          descr_f12f12_key, i1i2_f12f12_i1i2);
     activate_ints_f12f12(moints4_rtime, ji_ij,
                          v_orbs1, v_orbs2,
                          descr_f12f12_key, i2i1_f12f12_i1i2);
     activate_ints_f12f12(moints4_rtime, ij_ji,
                          v_orbs1, v_orbs2,
                          descr_f12f12_key, i1i2_f12f12_i2i1);
     activate_ints_f12f12(moints4_rtime, ji_ji,
                          v_orbs1, v_orbs2,
                          descr_f12f12_key, i2i1_f12f12_i2i1);
     // activate f12 ints
     activate_ints_X_f12(moints4_rtime, "ij",
                         v_orbs1, v_orbs2,
                         descr_f12_key, f12_ints_i1i2);
     activate_ints_X_f12(moints4_rtime, "ji",
                         v_orbs1, v_orbs2,
                         descr_f12_key, f12_ints_i2i1);
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

    const int nocc_act = v_orbs1[0]->rank();
    if (debug_ >= DefaultPrintThresholds::N2)
      ExEnv::out0() << endl << spinletters << " D^i_i 2" << endl
                    << "number of occupied orbital: " << nocc_act << endl;

    // skip spincase if no electron of this kind
    if (nocc_act == 0)
      continue;

    // initialize D^i_i[s]
    RefSCDimension rowdim_occ = new SCDimension(nocc_act);
    RefSCDimension coldim_occ = new SCDimension(nocc_act);
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    Dii_[s] = localkit->matrix(rowdim_occ, coldim_occ);
    Dii_[s].assign(0.0);

    // calculate AlphaAlpha/BetaBeta part for D^i_i:
    // R^ij_ab R^ab_ij & R^ji_ab R^ab_ij
    // which are contracted over a, b, and j
    double* RRij_ij = new double[nocc_act];
    double* RRji_ij = new double[nocc_act];
    fill_n(RRij_ij, nocc_act, 0.0);
    fill_n(RRji_ij, nocc_act, 0.0);

    // R^ij_ab R^ab_ij = (f12f12)^ij_ij
    //                 - f12^ij_pq f12^pq_ij
    //                 - f12^ij_ma' f12^ma'_ij
    //                 - f12^ij_a'm f12^a'm_ij

    // activate_ints: (f12f12)^ij_ij and (f12)^ij_k1k2 ints
    Ref<DistArray4> f12f12_ints_ijij;
    vector<Ref<DistArray4> > f12_ints_ij;

    activate_ints_f12f12(moints4_rtime, ij_ij,
                         v_orbs1, v_orbs2,
                         descr_f12f12_key, f12f12_ints_ijij);
    activate_ints_X_f12(moints4_rtime, "ij",
                        v_orbs1, v_orbs2,
                        descr_f12_key, f12_ints_ij);

    // for RR^ij_ij:
    // i is  1st and 2st index (for label), j: 3rd index (sum over j)
    compute_RR_sum_abj(RR13_23, nocc_act,
                       f12f12_idx, f12_idx,
                       f12f12_ints_ijij, f12_ints_ij, f12_ints_ij,
                       RRij_ij);

    // R^ji_ab R^ab_ij:
    Ref<DistArray4> f12f12_ints_jiij;
    vector<Ref<DistArray4> > f12_ints_ji;

    activate_ints_f12f12(moints4_rtime, ji_ij,
                         v_orbs1, v_orbs2,
                         descr_f12f12_key, f12f12_ints_jiij);
    activate_ints_X_f12(moints4_rtime, "ji",
                        v_orbs1, v_orbs2,
                        descr_f12_key, f12_ints_ji);

    compute_RR_sum_abj(RR31_23, nocc_act,
                       f12f12_idx, f12_idx,
                       f12f12_ints_jiij, f12_ints_ji, f12_ints_ij,
                       RRji_ij);

    f12f12_ints_ijij->deactivate();
    f12f12_ints_jiij->deactivate();
    for (int i = 0; i != f12_ints_ij.size(); ++i) {
      f12_ints_ij[i]->deactivate();
      f12_ints_ji[i]->deactivate();
    }

    if (debug_ >= DefaultPrintThresholds::N2) {
      print_intermediate(spinletters, "R^ij_ab R^ab_ij", RRij_ij, nocc_act);
      print_intermediate(spinletters, "R^ji_ab R^ab_ij", RRji_ij, nocc_act);
    }

    // calculate AlphaBeta part for D^i_i:
    double* RR1 = NULL;
    double* RR2 = NULL;
    double* RR3 = NULL;
    double* RR4 = NULL;
    if (nspincases2 == 3) {

      RR1 = new double[nocc_act];
      RR2 = new double[nocc_act];
      RR3 = new double[nocc_act];
      RR4 = new double[nocc_act];
      fill_n(RR1, nocc_act, 0.0);
      fill_n(RR2, nocc_act, 0.0);
      fill_n(RR3, nocc_act, 0.0);
      fill_n(RR4, nocc_act, 0.0);

      if (spin == Alpha) {
        // R^i1j2_a1b2 R^a1b2_i1j2
        compute_RR_sum_abj(RR13_23, nocc_act,
                           f12f12_idx, f12_idx, i1i2_f12f12_i1i2,
                           f12_ints_i1i2, f12_ints_i1i2, RR1);
        // R^j2i1_a1b2 R^a1b2_i1j2
        compute_RR_sum_abj(RR31_23, nocc_act,
                           f12f12_idx, f12_idx, i2i1_f12f12_i1i2,
                           f12_ints_i2i1, f12_ints_i1i2, RR2);
        // R^i1j2_a1b2 R^a1b2_j2i1
        compute_RR_sum_abj(RR13_32, nocc_act,
                           f12f12_idx, f12_idx, i1i2_f12f12_i2i1,
                           f12_ints_i1i2, f12_ints_i2i1, RR3);
        // R^j2i1_a1b2 R^a1b2_j2i1
        compute_RR_sum_abj(RR31_32, nocc_act,
                           f12f12_idx, f12_idx, i2i1_f12f12_i2i1,
                           f12_ints_i2i1, f12_ints_i2i1, RR4);
      } else {
          // R^j1i2_a1b2 R^a1b2_j1i2
          compute_RR_sum_abj(RR31_32, nocc_act,
                             f12f12_idx, f12_idx, i1i2_f12f12_i1i2,
                             f12_ints_i1i2, f12_ints_i1i2, RR1);
          // R^i2j1_a1b2 R^a1b2_j1i2
          compute_RR_sum_abj(RR13_32, nocc_act,
                             f12f12_idx, f12_idx, i2i1_f12f12_i1i2,
                             f12_ints_i2i1, f12_ints_i1i2, RR2);
          // R^j1i2_a1b2 R^a1b2_i2j1
          compute_RR_sum_abj(RR31_23, nocc_act,
                             f12f12_idx, f12_idx, i1i2_f12f12_i2i1,
                             f12_ints_i1i2, f12_ints_i2i1, RR3);
          // R^i2j1_a1b2 R^a1b2_i2j1
          compute_RR_sum_abj(RR13_23, nocc_act,
                             f12f12_idx, f12_idx, i2i1_f12f12_i2i1,
                             f12_ints_i2i1, f12_ints_i2i1, RR4);
      }

      if (debug_ >= DefaultPrintThresholds::N2) {
        const string spinlabel = (spin == Alpha? "alpha" : "beta");
        const string RR1_label = (spin == Alpha? "R^i1j2_a1b2 R^a1b2_i1j2"
                                               : "R^j1i2_a1b2 R^a1b2_j1i2");
        const string RR2_label = (spin == Alpha? "R^j2i1_a1b2 R^a1b2_i1j2"
                                               : "R^i2j1_a1b2 R^a1b2_j1i2");
        const string RR3_label = (spin == Alpha? "R^i1j2_a1b2 R^a1b2_j2i1"
                                               : "R^j1i2_a1b2 R^a1b2_i2j1");
        const string RR4_label = (spin == Alpha? "R^j2i1_a1b2 R^a1b2_j2i1"
                                               : "R^i2j1_a1b2 R^a1b2_i2j1");
        print_intermediate(spinlabel, RR1_label, RR1, nocc_act);
        print_intermediate(spinlabel, RR2_label, RR2, nocc_act);
        print_intermediate(spinlabel, RR3_label, RR3, nocc_act);
        print_intermediate(spinlabel, RR4_label, RR4, nocc_act);
      }

    } else if (nspincases2 == 2) {
        RR1 = RRij_ij;
        RR2 = RRji_ij;
        RR3 = RRji_ij;
        RR4 = RRij_ij;
   }


    // calculate D^i_i
    const double* iter_RRij_ij = RRij_ij;
    const double* iter_RRji_ij = RRji_ij;
    const double* iter_RR1 = RR1;
    const double* iter_RR2 = RR2;
    const double* iter_RR3 = RR3;
    const double* iter_RR4 = RR4;

    for (int i = 0;  i < nocc_act; ++i) {

      // AlphaAlpha/BetaBeta part
      double dii_12 = C_1 * C_1 * (*iter_RRij_ij - *iter_RRji_ij);
      ExEnv::out0() << spinletters << " part d^" << i << "_" << i << " = "
                    << scprintf("%12.10f", dii_12) << endl;

      ++iter_RRij_ij;
      ++iter_RRji_ij;

      // AlphaBeta part
      if (nocc_alpha != 0 && nocc_beta != 0) {
  //        ExEnv::out0() << "RR1: "  << scprintf("%12.10f", RR1)
  //                      << "  RR2: " << scprintf("%12.10f", RR2)
  //                      << "  RR3: "  << scprintf("%12.10f", RR3)
  //                      << "  RR4: " << scprintf("%12.10f", RR4)
  //                      << endl;

         dii_12 += pow(0.5 * (C_0 + C_1), 2) * (*iter_RR1)
                + 0.25 * (C_0 * C_0 - C_1 * C_1) * (*iter_RR2 + *iter_RR3)
                + pow(0.5 * (C_0 - C_1), 2) * (*iter_RR4);
         ExEnv::out0() << "AlphaBeta part: d^"  << i << "_" << i << " = "
                        << scprintf("%12.10f", dii_12) << endl;

         ++iter_RR1;
         ++iter_RR2;
         ++iter_RR3;
         ++iter_RR4;
       } // end of AlphaBeta part

      Dii_[s].set_element(i, i, dii_12);
    }

     delete[] RRij_ij;
     delete[] RRji_ij;
     if (nspincases2 == 3) {
       delete[] RR1;
       delete[] RR2;
       delete[] RR3;
       delete[] RR4;
     }
  //    // Assign m to a subblock of this.
  //    virtual void assign_subblock(SCMatrix *m, int, int, int, int, int=0, int=0) =0;
  //    D[s].assign_subblock(Dii[s],);
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
} // end of computation of Di_i 2nd version


void MP2R12Energy_Diag::compute_Dbc(const int nspincases1, const int nspincases2,
                                    const double C_0, const double C_1)
{
  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);

  // obtain occ, vir, orbs, and cabs orbitals for AlphaBeta case
  vector< Ref<OrbitalSpace> > v_orbs1_ab;
  vector< Ref<OrbitalSpace> > v_orbs2_ab;
  obtain_orbitals(AlphaBeta, v_orbs1_ab, v_orbs2_ab);

  const Ref<OrbitalSpace> occ1_act = v_orbs1_ab[0];
  const Ref<OrbitalSpace> occ2_act = v_orbs2_ab[0];
  const Ref<OrbitalSpace> vir1 = v_orbs1_ab[1];
  const Ref<OrbitalSpace> vir2 = v_orbs2_ab[1];
  const Ref<OrbitalSpace> cabs1 = v_orbs1_ab[3];
  const Ref<OrbitalSpace> cabs2 = v_orbs2_ab[3];
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
     vector< Ref<OrbitalSpace> > v_orbs1;          // orbitals of spin1
     vector< Ref<OrbitalSpace> > v_orbs2;          // orbitals of spin2
     obtain_orbitals(spincase, v_orbs1, v_orbs2);

     const Ref<OrbitalSpace>& occ_act = v_orbs1[0];
     const int nocc_act = occ_act->rank();

     // skip spincase if no electron of this kind
     if (nocc_act == 0)
       continue;

     const Ref<OrbitalSpace> vir = v_orbs1[1];
     const Ref<OrbitalSpace> cabs = v_orbs1[3];
     const int nvir = vir->rank();
     const int ncabs = cabs->rank();
     if (debug_ >= DefaultPrintThresholds::N2)
       ExEnv::out0() << endl << spinletters << " D^b_c" << endl
                     << "number of occupied orbital: " << nocc_act << endl
                     << "number of virtual orbital: " << nvir << endl
                     << "number of cabs orbital: " << ncabs << endl;

     // initialize D^b_c[s]
     RefSCDimension rowdim = new SCDimension(nvir);
     RefSCDimension coldim = new SCDimension(nvir);
     Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
     Dbc_[s] = localkit->matrix(rowdim, coldim);
     Dbc_[s].assign(0.0);

     // calculate AlphaAlpha/BetaBeta part for D^b_c:
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
     compute_RR2_sum_3idx(RR31_32, f12_idx,
                          apaii_f12_ints, apaii_f12_ints,
                          RRapb_apc);
     // R^ba'_ij R^ij_a'c
     compute_RR2_sum_3idx(RR13_32, f12_idx,
                          aapii_f12_ints, apaii_f12_ints,
                          RRbap_apc);

     apaii_f12_ints->deactivate();
     aapii_f12_ints->deactivate();

     if (debug_ >= DefaultPrintThresholds::N2) {
       print_intermediate(spinletters, "R^a'b_ij R^ij_a'c", RRapb_apc, nvir, nvir);
       print_intermediate(spinletters, "R^ba'_ij R^ij_a'c", RRbap_apc, nvir, nvir);
     }

     // calculate AlphaBeta part for D^b_c:
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

         compute_RR2_sum_3idx(RR13_23, f12_idx,
                              a1ap2i1i2_f12_ints, a1ap2i1i2_f12_ints,
                              RR1); // R^ba'_ij R^ij_ca'
         compute_RR2_sum_3idx(RR31_23, f12_idx,
                              ap2a1i1i2_f12_ints, a1ap2i1i2_f12_ints,
                              RR2); // R^a'b_ij R^ij_ca'
         compute_RR2_sum_3idx(RR13_32, f12_idx,
                              a1ap2i1i2_f12_ints, ap2a1i1i2_f12_ints,
                              RR3); // R^ba'_ij R^ij_a'c
         compute_RR2_sum_3idx(RR31_32, f12_idx,
                              ap2a1i1i2_f12_ints, ap2a1i1i2_f12_ints,
                              RR4); // R^a'b_ij R^ij_a'c
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

       if (debug_ >= DefaultPrintThresholds::N2) {
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
    }// end of AlphaBeta part for D^b_c

     // calculate D^b_c[s]
     const double* iter_RRapb_apc = RRapb_apc;
     const double* iter_RRbap_apc = RRbap_apc;
     const double* iter_RR1 = RR1;
     const double* iter_RR2 = RR2;
     const double* iter_RR3 = RR3;
     const double* iter_RR4 = RR4;

     for (int b = 0;  b < nvir; ++b) {

       double dbc_12 = 0.0;
       for (int c = 0; c < nvir; ++c) {

         // AlphaAlpha/BetaBeta part
         dbc_12 = C_1 * C_1 * (*iter_RRapb_apc - *iter_RRbap_apc);
         ExEnv::out0() << spinletters << " part d^" << b << "_" << c << " = "
                       << scprintf("%12.10f", dbc_12) << endl;

         ++iter_RRapb_apc;
         ++iter_RRbap_apc;

         // AlphaBeta part
         if (nocc_alpha != 0 && nocc_beta != 0) {
   //        ExEnv::out0() << "RR1: "  << scprintf("%12.10f", RR1)
   //                      << "  RR2: " << scprintf("%12.10f", RR2)
   //                      << "  RR3: "  << scprintf("%12.10f", RR3)
   //                      << "  RR4: " << scprintf("%12.10f", RR4)
   //                      << endl;

           dbc_12 += pow(0.5 * (C_0 + C_1), 2) * (*iter_RR1)
                  + 0.25 * (C_0 * C_0 - C_1 * C_1) * (*iter_RR2 + *iter_RR3)
                  + pow(0.5 * (C_0 - C_1), 2) * (*iter_RR4);
           ExEnv::out0() << "AlphaBeta part: d^"  << b << "_" << c << " = "
                         << scprintf("%12.10f", dbc_12) << endl;

           ++iter_RR1;
           ++iter_RR2;
           ++iter_RR3;
           ++iter_RR4;
         } // end of AlphaBeta part

         Dbc_[s].set_element(b, c, dbc_12);
       } // end of looping over c

     } // end of calculating D^b_c[s]

      delete[] RRapb_apc;
      delete[] RRbap_apc;
      if (nspincases2 == 3) {
        delete[] RR1;
        delete[] RR2;
        delete[] RR3;
        delete[] RR4;
      }
   //    // Assign m to a subblock of this.
   //    virtual void assign_subblock(SCMatrix *m, int, int, int, int, int=0, int=0) =0;
   //    D[s].assign_subblock(Dbc[s],);
   } // end of spincase1 loop

} // end of computation of Dbc
void MP2R12Energy_Diag::compute_Dapa(const int nspincases1, const int nspincases2,
                                     const int nocc_alpha, const int nocc_beta,
                                     const double C_0, const double C_1,
                                     RefSCMatrix Dapa[NSpinCases1])
{
  // compute T2^ab_ij amplitudes in MP2R12 method
  vector<double*> T2ab_ij(NSpinCases2);
  for(int s=0; s<nspincases2; ++s) {
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
     compute_T2_mp2(spincase2, v_orbs1, v_orbs2, T2ab_ij[s]);

     // print T2^ab_ij
 #if 1
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

    const int nocc_act = v_orbs1[0]->rank();
    // skip spincase if no electron of this kind
    if (nocc_act == 0)
      continue;

     const int nvir = v_orbs1[1]->rank();
     const int ncabs = v_orbs1[3]->rank();
     const int ncabs_vir12 = ncabs * nvir * nvir;
     const int ncom_orbs = nocc_act + nvir + ncabs;

     if (debug_ >= DefaultPrintThresholds::N2)
       ExEnv::out0() << endl << spinletters << " D^a'_a" << endl
                     << "number of occupied orbital: " << nocc_act << endl
                     << "number of virtual orbital: " << nvir << endl
                     << "number of cabs orbital: " << ncabs << endl;
     // initialize D^ap_a[s]
     RefSCDimension rowdim_cabs = new SCDimension(ncabs);
     RefSCDimension coldim_vir = new SCDimension(nvir);
     Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
     Dapa[s] = localkit->matrix(rowdim_cabs, coldim_vir);
     Dapa[s].assign(0.0);

     // calculate AlphaAlpha/BetaBeta part for D^a'_a:
     // AlphaAlpha: R^ij_a'b * T2^ab_ij & R^ji_a'b * T2^ab_ij
     // BetaBeta: R^ij_ba' * T2^ba_ij & R^ji_ba' * T2^ba_ij
     double* RT2ij = new double[ncabs_vir12];
     double* RT2ji = new double[ncabs_vir12];
     fill_n(RT2ij, ncabs_vir12, 0.0);
     fill_n(RT2ji, ncabs_vir12, 0.0);

     compute_RT2_mp2(spin, spincase,
                     v_orbs1, v_orbs2,
                     T2ab_ij[spincase], RT2ij, RT2ji);
     delete[] T2ab_ij[spincase];

     if (debug_ >= DefaultPrintThresholds::N2){
       const string  RT2ij_label = (spin == Alpha ?
                                   "R^ij_a'b * T^ab_ij" : "R^ij_ba' * T2^ba_ij");
       const string  RT2ji_label = (spin == Alpha ?
                                   "R^ji_a'b * T^ab_ij" : "R^ji_ba' * T2^ba_ij");
       print_RT2(spinletters, RT2ij_label, RT2ij,
                 ncabs, nvir, nvir);
       print_RT2(spinletters, RT2ji_label, RT2ij,
                 ncabs, nvir, nvir);
     }

     // calculate AlphaBeta for D^a'_a:
     double* RT2ij_ab = NULL;
     double* RT2ji_ab = NULL;
     if (nspincases2 == 3 && nocc_alpha != 0 && nocc_beta != 0) {
       // obtain occ, vir, orbs, and cabs orbitals
       vector< Ref<OrbitalSpace> > v_orbs1;
       vector< Ref<OrbitalSpace> > v_orbs2;
       obtain_orbitals(AlphaBeta, v_orbs1, v_orbs2);

       // R^ij_a'b * T2^ab_ij & R^ji_a'b * T2^ab_ij
       // for R * T2, AlphaBeta part for spin alpha and beta is different
 //        const int nvir1 = v_orbs1[1]->rank();
 //        const int nvir2 = v_orbs2[1]->rank();
 //        const int nvir12 = nvir1 * nvir2;
 //        const int ncabs1 = v_orbs1[3]->rank();
 //        const int ncabs2 = v_orbs2[3]->rank();
 //        const int ncabs_vir12 = (spin == Alpha ? ncabs1 * nvir12 : ncabs2 * nvir12);
 //
 //        RT2ij_ab = new double[ncabs_vir12];
 //        RT2ji_ab = new double[ncabs_vir12];
 //        fill_n(RT2ij_ab, ncabs_vir12, 0.0);
 //        fill_n(RT2ji_ab, ncabs_vir12, 0.0);
 //
 //        compute_RT2_mp2(spin, AlphaBeta,
 //                        v_orbs1, v_orbs2,
 //                        T2ab_ij[AlphaBeta], RT2ij_ab, RT2ji_ab);
 //        delete[] T2ab_ij[AlphaBeta];

     } else if (nspincases2 == 2 && nocc_alpha != 0 && nocc_beta != 0) {
         RT2ij_ab = RT2ij;
         RT2ji_ab = RT2ji;
    } // end of AlphaBeta part for D^a'_a

   // calculate D^a'_a
 //    if (this->r12eval()->coupling() == true) {
 //      compute_Dapa_spin(spin, C_0, C_1,
 //                   nocc_alpha, nocc_beta, ncabs,
 //                   RT2ij, RT2ji, RT2ij_ab, RT2ji_ab,
 //                   Dapa[s]);
 //      delete[] RT2ij;
 //      delete[] RT2ji;
 //      if (nspincases2 == 3) {
 //        delete[] RT2ij_ab;
 //        delete[] RT2ji_ab;
 //      }
 //    }
   //    /// Assign m to a subblock of this.
   //    virtual void assign_subblock(SCMatrix *m, int, int, int, int, int=0, int=0) =0;
   // D[s].assign_subblock(Dii[s],);
   } // end of loop spincase1

} // end of compute_Dapa


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
                                             double* D_alpha, double* D_beta)
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
    ExEnv::out0() << endl << "number of alpha active occupied orbital: " << occ1_act->rank() << endl
                  <<"number of beta active occupied orbital: " << occ2_act->rank() << endl
                  << "number of alpha virtula orbital: " << vir1->rank() << endl
                  <<"number of beta virtual orbital: " << vir2->rank() << endl
                  << "number of alpha cabs: " << cabs1->rank() << endl
                  <<"number of beta cabs: " << cabs2->rank() << endl;

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
    idx3_orbs1 = cabs1;
    idx3_orbs2 = cabs2;
  break;

  case cabs_vir_cabs:
    idx1_orbs1 = cabs1;
    idx1_orbs2 = cabs2;
    idx2_orbs1 = vir1;
    idx2_orbs2 = vir2;
    idx3_orbs1 = cabs1;
    idx3_orbs2 = cabs2;
  break;

  default:
    ExEnv::out0() << "There is no such case";
    break;
  }

  // test
  if (debug_ >= DefaultPrintThresholds::N2) {
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
       idx3_orbs = cabs;
     break;

     case cabs_vir_cabs:
       size_idx1 = ncabs;
       size_idx2 = ncabs;
       idx1_orbs = cabs;
       idx2_orbs = vir;
       idx3_orbs = cabs;
     break;
     }

     if (debug_ >= DefaultPrintThresholds::N2) {
       int size_idx3 = ncabs;
       if (orbitals_label == cabs_cabs_vir) size_idx3 = nvir;
       ExEnv::out0() << endl << spinletters << endl
                     << "number of index1 orbital: " << idx1_orbs->rank() << endl
                     << "number of index2 orbital: " << idx2_orbs->rank() << endl
                     << "number of index3 orbital: " << idx3_orbs->rank() << endl;
     }

     // test
    {
       Ref<DistArray4> f12_ints1;
       activate_ints(cabs->id(), vir->id(), occ_act->id(), occ_act->id(),
                     descr_f12_key, moints4_rtime, f12_ints1);
//       print_f12_ints(spinletters, "R^A'B_IJ ints", f12_idx, f12_ints1);

       Ref<DistArray4> f12_ints2;
       activate_ints(vir->id(), cabs->id(), occ_act->id(), occ_act->id(),
                     descr_f12_key, moints4_rtime, f12_ints2);
//       print_f12_ints(spinletters, "R^B'A_IJ ints", f12_idx, f12_ints2);

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
             const double* blk2 = f12_ints1->retrieve_pair_block(idx3, idx2, f12_idx);
             const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

             RR_sum_idx3 += blk12;
             f12_ints1->release_pair_block(idx3, idx1, f12_idx);
             f12_ints1->release_pair_block(idx3, idx2, f12_idx);

           }
           *iter_RR = RR_sum_idx3;
           ++iter_RR;
         }
       }
       print_intermediate(spinletters, "testing result for R^A'B_IJ R^IJ_A'C", RR, size_idx1, size_idx2);

       iter_RR = RR;
       for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
         for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

           double RR_sum_idx3 = 0;
           for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
             const double* blk1 = f12_ints2->retrieve_pair_block(idx1, idx3, f12_idx);
             const double* blk2 = f12_ints1->retrieve_pair_block(idx3, idx2, f12_idx);
             const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

             RR_sum_idx3 += blk12;
             f12_ints2->release_pair_block(idx1, idx3, f12_idx);
             f12_ints1->release_pair_block(idx3, idx2, f12_idx);

           }
           *iter_RR = RR_sum_idx3;
           ++iter_RR;
         }
       }
       print_intermediate(spinletters, "testing result for R^BA'_IJ R^IJ_A'C", RR, size_idx1, size_idx2);
       delete[] RR;

       f12_ints1->deactivate();
       f12_ints2->deactivate();
     } // end of testing ints


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

     if (debug_ >= DefaultPrintThresholds::N2) {
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
        {
           Ref<DistArray4> f12_ints1;
           activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, f12_ints1);

           Ref<DistArray4> f12_ints2;
           activate_ints(cabs2->id(), vir1->id(), occ1_act->id(), occ2_act->id(),
                         descr_f12_key, moints4_rtime, f12_ints2);

           const blasint one = 1;
           const blasint nocc12 = nocc_alpha * nocc_beta;
           const int size_idx1 = vir1->rank();
           const int size_idx2 = vir1->rank();
           const int size_idx3 = cabs2->rank();
           double* RR = new double[size_idx1 * size_idx2];

           double* iter_RR = RR;
           for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
             for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

               double RR_sum_idx3 = 0;
               for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                 const double* blk1 = f12_ints1->retrieve_pair_block(idx1, idx3, f12_idx);
                 const double* blk2 = f12_ints1->retrieve_pair_block(idx2, idx3, f12_idx);
                 const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                 RR_sum_idx3 += blk12;
                 f12_ints1->release_pair_block(idx1, idx3, f12_idx);
                 f12_ints1->release_pair_block(idx2, idx3, f12_idx);

               }
               *iter_RR = RR_sum_idx3;
               ++iter_RR;
             }
           }
           print_intermediate(spinletters, "testing result for R^BA'_IJ R^IJ_CA'", RR, size_idx1, size_idx2);

           iter_RR = RR;
           for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
             for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

               double RR_sum_idx3 = 0;
               for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                 const double* blk1 = f12_ints2->retrieve_pair_block(idx3, idx1, f12_idx);
                 const double* blk2 = f12_ints1->retrieve_pair_block(idx2, idx3, f12_idx);
                 const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                 RR_sum_idx3 += blk12;
                 f12_ints2->release_pair_block(idx3, idx1, f12_idx);
                 f12_ints1->release_pair_block(idx2, idx3, f12_idx);

               }
               *iter_RR = RR_sum_idx3;
               ++iter_RR;
             }
           }
           print_intermediate(spinletters, "testing result for R^AB'_IJ R^IJ_C'A", RR, size_idx1, size_idx2);

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
           print_intermediate(spinletters, "testing result for R^B'A_IJ R^IJ_AC'", RR, size_idx1, size_idx2);

           iter_RR = RR;
            for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
              for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

                double RR_sum_idx3 = 0;
                for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                  const double* blk1 = f12_ints2->retrieve_pair_block(idx3, idx1, f12_idx);
                  const double* blk2 = f12_ints2->retrieve_pair_block(idx3, idx2, f12_idx);
                  const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                  RR_sum_idx3 += blk12;
                  f12_ints2->release_pair_block(idx3, idx1, f12_idx);
                  f12_ints2->release_pair_block(idx3, idx2, f12_idx);

                }
                *iter_RR = RR_sum_idx3;
                ++iter_RR;
              }
            }
            print_intermediate(spinletters, "testing result for R^AB'_IJ R^IJ_AC'", RR, size_idx1, size_idx2);

           delete[] RR;
           f12_ints1->deactivate();
           f12_ints2->deactivate();
         } // end of testing ints

       } else {
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

           // test for AlphaBeta of beta spin
          {
             Ref<DistArray4> f12_ints1;
             activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                           descr_f12_key, moints4_rtime, f12_ints1);

             Ref<DistArray4> f12_ints2;
             activate_ints(cabs2->id(), vir1->id(), occ1_act->id(), occ2_act->id(),
                           descr_f12_key, moints4_rtime, f12_ints2);

             const blasint one = 1;
             const blasint nocc12 = nocc_alpha * nocc_beta;
             const int size_idx1 = cabs2->rank();
             const int size_idx2 = cabs2->rank();
             const int size_idx3 = vir1->rank();
             double* RR = new double[size_idx1 * size_idx2];

             double* iter_RR = RR;
             for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
               for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

                 double RR_sum_idx3 = 0;
                 for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                   const double* blk1 = f12_ints1->retrieve_pair_block(idx3, idx1, f12_idx);
                   const double* blk2 = f12_ints1->retrieve_pair_block(idx3, idx2, f12_idx);
                   const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                   RR_sum_idx3 += blk12;
                   f12_ints1->release_pair_block(idx3, idx1, f12_idx);
                   f12_ints1->release_pair_block(idx3, idx2, f12_idx);

                 }
                 *iter_RR = RR_sum_idx3;
                 ++iter_RR;
               }
             }
             print_intermediate(spinletters, "testing result for R^AB'_IJ R^IJ_AC'", RR, size_idx1, size_idx2);

             iter_RR = RR;
             for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
               for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

                 double RR_sum_idx3 = 0;
                 for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                   const double* blk1 = f12_ints2->retrieve_pair_block(idx1, idx3, f12_idx);
                   const double* blk2 = f12_ints1->retrieve_pair_block(idx3, idx2, f12_idx);
                   const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                   RR_sum_idx3 += blk12;
                   f12_ints2->release_pair_block(idx1, idx3, f12_idx);
                   f12_ints1->release_pair_block(idx3, idx2, f12_idx);

                 }
                 *iter_RR = RR_sum_idx3;
                 ++iter_RR;
               }
             }
             print_intermediate(spinletters, "testing result for R^B'A_IJ R^IJ_AC'", RR, size_idx1, size_idx2);

             iter_RR = RR;
             for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
               for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

                 double RR_sum_idx3 = 0;
                 for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                   const double* blk1 = f12_ints1->retrieve_pair_block(idx3, idx1, f12_idx);
                   const double* blk2 = f12_ints2->retrieve_pair_block(idx2, idx3, f12_idx);
                   const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                   RR_sum_idx3 += blk12;
                   f12_ints1->release_pair_block(idx3, idx1, f12_idx);
                   f12_ints2->release_pair_block(idx2, idx3, f12_idx);
                 }
                 *iter_RR = RR_sum_idx3;
                 ++iter_RR;
               }
             }
             print_intermediate(spinletters, "testing result for R^AB'_IJ R^IJ_C'A", RR, size_idx1, size_idx2);

             iter_RR = RR;
              for(int idx1 = 0 ; idx1 < size_idx1; ++idx1) {
                for(int idx2 = 0; idx2 < size_idx2; ++idx2) {

                  double RR_sum_idx3 = 0;
                  for(int idx3 = 0; idx3 < size_idx3; ++idx3) {
                    const double* blk1 = f12_ints2->retrieve_pair_block(idx1, idx3, f12_idx);
                    const double* blk2 = f12_ints2->retrieve_pair_block(idx2, idx3, f12_idx);
                    const double blk12 = F77_DDOT(&nocc12, blk1, &one, blk2, &one);

                    RR_sum_idx3 += blk12;
                    f12_ints2->release_pair_block(idx1, idx3, f12_idx);
                    f12_ints2->release_pair_block(idx2, idx3, f12_idx);

                  }
                  *iter_RR = RR_sum_idx3;
                  ++iter_RR;
                }
              }
              print_intermediate(spinletters, "testing result for R^B'A_IJ R^IJ_C'A", RR, size_idx1, size_idx2);

             delete[] RR;
             f12_ints1->deactivate();
             f12_ints2->deactivate();
           } // end of testing ints
       }

       if (debug_ >= DefaultPrintThresholds::N2) {
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

       double d_12 = 0.0;
       for (int idx2 = 0; idx2 < size_idx2; ++idx2) {

         // AlphaAlpha/BetaBeta part
         d_12 = C_1 * C_1 * (*iter_RR31_32_array - *iter_RR13_32_array);
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

     //
     iter_D = (spin == Alpha? D_alpha : D_beta);
     print_intermediate(spinletters, "D^B_C", iter_D, size_idx1, size_idx2);

     delete[] RR31_32_array;
     delete[] RR13_32_array;
     if (nspincases2 == 3) {
       delete[] RR1;
       delete[] RR2;
       delete[] RR3;
       delete[] RR4;
     }
  } // end of spincase1 loop

  if (nspincases2 == 2)
    D_beta = D_alpha;
} // end of computation of RR_sum_ijidx2_spinform


void MP2R12Energy_Diag::compute_density_diag() {

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


  // test propose
  ExEnv::out0() << endl << "number of alpha active occupied orbital: " << nocc1_act << endl
                        <<"number of beta active occupied orbital: " << nocc2_act << endl
                        << "number of alpha virtula orbital: " << nvir1 << endl
                        <<"number of beta virtual orbital: " << nvir2 << endl
                        << "number of alpha cabs: " << ncabs1 << endl
                        <<"number of beta cabs: " << ncabs2 << endl;
  // test ints for open shell and close shell
 {
    Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
    Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

    Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
    const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
    const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
    const unsigned int f12_idx = descr_f12->intset(f12_type);

    const SpinCase1 spin1 = case1(AlphaBeta);
    const SpinCase1 spin2 = case2(AlphaBeta);

    const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
    const Ref<OrbitalSpace>& vir1 = r12eval()->vir(spin1);
    const Ref<OrbitalSpace>& orbs1 = r12eval()->orbs(spin1);
//    const Ref<OrbitalSpace>& cabs1 = r12world->cabs_space(spin1);
    const Ref<OrbitalSpace>& occ1 = r12eval()->occ(spin1);

    const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
    const Ref<OrbitalSpace>& vir2 = r12eval()->vir(spin2);
    const Ref<OrbitalSpace>& orbs2 = r12eval()->orbs(spin2);
//    const Ref<OrbitalSpace>& cabs2 = r12world->cabs_space(spin2);
    const Ref<OrbitalSpace>& occ2 = r12eval()->occ(spin2);

    // use CABS orbitals canonicalized by diagonalizing CABS/CABS H(core) (cheaper than the Fock matrix)
    Ref<OrbitalSpace> cabs1_canon = r12eval()->cabs_space_hcanonical(spin1);
    Ref<OrbitalSpace> cabs2_canon = r12eval()->cabs_space_hcanonical(spin2);

    Ref<DistArray4> f12_ints;
//    cabs1->print_detail();
    cabs1_canon->print_detail();
    cabs2_canon->print_detail();
    activate_ints(cabs1_canon->id(), cabs2_canon->id(), occ1_act->id(), occ2_act->id(),
                  descr_f12_key, moints4_rtime, f12_ints);
    print_f12_ints("AlphaBeta", "R^AB'_IJ ints", f12_idx, f12_ints);
 }


  //
  // D^i_i
  // Alpha d^I_I = C1^2 * (R^IJ_AB R^AB_IJ - R^JI_AB R^AB_IJ)  (I, J, A, B in alpha orbitals)
  //
  //             + [(C0+C1)/2]^2 * R^IJ_AB R^AB_IJ   (I, A in alpha orbitals, J, B inbeta orbitals)
  //             + (C0+C1)/2*(C0-C1)/2 * (R^IJ_AB R^AB_JI + R^JI_AB R^AB_IJ)
  //             + [(C0-C1)/2]^2 * R^JI_AB R^AB_JI
  //
  // Beta  d^I_I = C1^2 * (R^IJ_AB R^AB_IJ - R^JI_AB R^AB_IJ)  (I, J, A, B in Beta orbitals)
  //
  //             + [(C0+C1)/2]^2 * R^JI_AB R^AB_JI    (J, A in alpha orbitals, I, B in beta orbitals)
  //             + (C0+C1)/2*(C0-C1)/2 * (R^IJ_AB R^AB_JI + R^JI_AB R^AB_IJ)
  //             + [(C0-C1)/2]^2 * R^IJ_AB R^AB_IJ

  // test
//{
//    Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
//    Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
//
//    Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
//    const string descr_f12_key = moints4_rtime->descr_key(descr_f12);
//    const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
//    const unsigned int f12_idx = descr_f12->intset(f12_type);
//
//    Ref<TwoBodyIntDescr> descr_f12f12 =
//        r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(), 0, 0);
//    const std::string descr_f12f12_key = moints4_rtime->descr_key(descr_f12f12);
//    const TwoBodyOper::type f12f12_type =
//        r12world->r12tech()->corrfactor()->tbint_type_f12f12();
//    const unsigned int f12f12_idx = descr_f12f12->intset(f12f12_type);
//
//  const int nocc12 = nocc1_act * nocc2_act;
//  double* Xij_ij = new double[nocc12];
//  fill_n(Xij_ij, nocc12, 0.0);
//
//  // store all the ints
//  std::vector<Ref<DistArray4> > f12_ij_ints;
//  std::vector<std::string> VX_output;
//
//  // X^ij_ij -= f^ij_pq f^pq_ij
//  Ref<DistArray4> i1i2p1p2_ints;
//  activate_ints(occ1_act->id(), occ2_act->id(), orbs1->id(), orbs2->id(),
//                descr_f12_key, moints4_rtime, i1i2p1p2_ints);
//  f12_ij_ints.push_back(i1i2p1p2_ints);
//  VX_output.push_back("diag-pq contribution");
//
//  // X^ij_ij -= f^ij_ma' f^ma'_ij
//  Ref<DistArray4> i1i2i1a2_ints;
//  activate_ints(occ1_act->id(), occ2_act->id(), occ1->id(), cabs2->id(),
//                descr_f12_key, moints4_rtime, i1i2i1a2_ints);
//
//  f12_ij_ints.push_back(i1i2i1a2_ints);
//  VX_output.push_back("diag-pq-ma' contribution");
//
//  // X^ij_ij -= f^ij_a'm f^a'm_ij
//  // TODO: for RHF can simply scale previous contribution by 2 and "symmetrize" at the end V^ij_ij = 0.5 * (V^ij_ij + V^ji_ji)
//  Ref<DistArray4> i1i2a1i2_ints;
//  activate_ints(occ1_act->id(), occ2_act->id(), cabs1->id(), occ2->id(),
//                descr_f12_key, moints4_rtime, i1i2a1i2_ints);
//
//  f12_ij_ints.push_back(i1i2a1i2_ints);
//  VX_output.push_back("diag-pq-ma'-a'm contribution");
//
//  // X^ij_ij += (f12f12)^ij_ij
//  Ref<DistArray4> i1i2i1i2_f12f12_ints;
//  activate_ints(occ1_act->id(), occ2_act->id(), occ1_act->id(),
//                occ2_act->id(), descr_f12f12_key, moints4_rtime,
//                i1i2i1i2_f12f12_ints);
//
//  compute_VX(ij_ij, VX_output, f12f12_idx, i1i2i1i2_f12f12_ints, f12_idx,
//             f12_idx, f12_ij_ints, f12_ij_ints, Xij_ij);
//  delete[] Xij_ij;
//}


  compute_Dii(nspincases1, nspincases2, nocc1_act, nocc2_act, C_0, C_1);

  // test
  compute_Dii_2(nspincases1, nspincases2, nocc1_act, nocc2_act, C_0, C_1);

  // D^b_c:
  // Alpha d^B_C = C1^2 * (R^A'B_IJ R^IJ_A'C - R^BA'_IJ R^IJ_A'C) (I, J, A', B, C in alpha orbitals)
  //             + [(C0+C1)/2]^2 * R^BA'_IJ R^IJ_CA'  (I, B, C in alpha orbitals, J, A' in beta orbitals)
  //             + (C0+C1)/2*(C0-C1)/2 * (R^A'B_IJ R^IJ_CA' + R^BA'_IJ R^IJ_A'C)
  //             + [(C0-C1)/2]^2 * R^A'B_IJ R^IJ_A'C
  //
  // Beta  d^B_C = C1^2 * (R^A'B_IJ R^IJ_A'C - R^BA'_IJ R^IJ_A'C) (I, J, A', B, C in beta orbitals)
  //             + [(C0+C1)/2]^2 * R^A'B_IJ R^IJ_A'C  (I, A' in alpha orbitals, J, B, C in beta orbitals)
  //             + (C0+C1)/2*(C0-C1)/2 * (R^BA'_IJ R^IJ_A'C + R^A'B_IJ R^IJ_CA')
  //             + [(C0-C1)/2]^2 * R^BA'_IJ R^IJ_CA'
//  compute_Dbc(nspincases1, nspincases2, C_0, C_1);
  // test
//  const int nvir11 = nvir1 * nvir1;
//  const int nvir22 = nvir2 * nvir2;
//  double* Db_c_alpha = new double[nvir11];
//  double* Db_c_beta = new double[nvir22];
//  fill_n(Db_c_alpha, nvir11, 0.0);
//  fill_n(Db_c_beta, nvir22, 0.0);
//
//  ExEnv::out0() << endl << "D^b_c : " << endl;
//  compute_RR31_32_spin(vir_vir_cabs,
//                       nspincases1, nspincases2, C_0, C_1,
//                       v_orbs1_ab, v_orbs2_ab,
//                       Db_c_alpha, Db_c_beta);
//  delete[] Db_c_alpha;
//  delete[] Db_c_beta;

  // D^b'_c':
  // Alpha d^B'_C' =
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
  // Beta  d^B'_C' =
  // 1st sum over A:  C1^2 * (R^AB'_IJ R^IJ_AC' - R^B'A_IJ R^IJ_AC') (I, J, A, B', C' in beta orbitals)
  //                + [(C0+C1)/2]^2 * R^AB'_IJ R^IJ_AC'  (I, A in alpha orbitals, J, B', C' in beta orbitals)
  //                + (C0+C1)/2*(C0-C1)/2 * (R^B'A_IJ R^IJ_AC' + R^AB'_IJ R^IJ_C'A)
  //                + [(C0-C1)/2]^2 * R^B'A_IJ R^IJ_C'A

  // 2nd sum over A': C1^2 * (R^A'B'_IJ R^IJ_A'C' - R^B'A'_IJ R^IJ_A'C') (I, J, A', B, C in beta orbitals)
  //                + [(C0+C1)/2]^2 * R^A'B'_IJ R^IJ_A'C'  (I, A in alpha orbitals, J, B', C' in beta orbitals)
  //                + (C0+C1)/2*(C0-C1)/2 * (R^B'A'_IJ R^IJ_A'C' + R^A'B'_IJ R^IJ_C'A')
  //                + [(C0-C1)/2]^2 * R^B'A'_IJ R^IJ_C'A'
//  const int ncabs11 = ncabs1 * ncabs1;
//  const int ncabs22 = ncabs2 * ncabs2;
//  double* Dbp_cp_alpha_A = new double[ncabs11];
//  double* Dbp_cp_beta_A = new double[ncabs22];
//  fill_n(Dbp_cp_alpha_A, ncabs11, 0.0);
//  fill_n(Dbp_cp_beta_A, ncabs22, 0.0);
//
//  ExEnv::out0() << endl << "D^b'_c' part1 (sum over a)" << endl;
//  compute_RR31_32_spin(cabs_cabs_vir,
//                       nspincases1, nspincases2, C_0, C_1,
//                       v_orbs1_ab, v_orbs2_ab,
//                       Dbp_cp_alpha_A, Dbp_cp_beta_A);

//  double* Dbp_cp_alpha_Ap = new double[ncabs11];
//  double* Dbp_cp_beta_Ap = new double[ncabs22];
//  fill_n(Dbp_cp_alpha_Ap, ncabs11, 0.0);
//  fill_n(Dbp_cp_beta_Ap, ncabs22, 0.0);
//
//  ExEnv::out0() << endl << "D^b'_c' part2 (sum over a')" << endl;
//  RR_sum_ijidx3_w_spin(cabs_cabs_cabs,
//                       nspincases1, nspincases2, C_0, C_1,
//                       v_orbs1_ab, v_orbs2_ab,
//                       Dbp_cp_alpha_Ap, Dbp_cp_beta_Ap);

//  delete[] Dbp_cp_alpha_A;
//  delete[] Dbp_cp_beta_A;
//  delete[] Dbp_cp_alpha_Ap;
//  delete[] Dbp_cp_beta_Ap;

  //
  // D^a'_a
  // Alpha d^A'_A = 1/2 * C1 * (R^IJ_A'B - R^JI_A'B) * T^AB_IJ  (I, J in alpha orbitals)
  //               + [(C0+C1)/2 * R^IJ_A'B + (C0-C1)/2 * R^JI_A'B] * T^AB_IJ  (I in alpha orbital, J in beta orbital)
  //
  //               + 1/2 * C1 * (R^B'A'_IJ R^IJ_B'A - R^A'B'_IJ R^IJ_B'A  (I, J in alpha orbitals)
  //               + [(C0+C1)/2]^2 * R^A'B'_IJ R^IJ_AB'  (I, A', A in alpha orbitals, J, B' in beta orbitals)
  //               + (C0+C1)/2*(C0-C1)/2 * (R^B'A'_IJ R^IJ_AB' + R^A'B'_IJ R^IJ_B'A)
  //               + [(C0-C1)/2]^2 * R^B'A'_IJ R^IJ_B'A
  //
  // Beta  d^A'_A = 1/2 * C1 * (R^IJ_BA'- R^JI_BA') * T^BA_IJ  (I, J in beta orbitals)
  //               + [(C0+C1)/2 * R^IJ_BA'+ (C0-C1)/2 * R^JI_BA'] * T^BA_IJ  (I in alpha orbital, J in beta orbital)
  //
  //               + C1^2 * (R^B'A'_IJ R^IJ_B'A - R^A'B'_IJ R^IJ_B'A) (I, J in beta orbitals)
  //               + [(C0+C1)/2]^2 * R^B'A'_IJ R^IJ_B'A  (I, B' in alpha orbitals, J, A', A in beta orbitals)
  //               + (C0+C1)/2*(C0-C1)/2 * (R^A'B'_IJ R^IJ_B'A + R^B'A'_IJ R^IJ_AB')
  //               + [(C0-C1)/2]^2 * R^A'B'_IJ R^IJ_AB'
//  const int ncabs1_vir1 = ncabs1 * nvir1;
//  const int ncabs2_vir2 = ncabs2 * nvir2;
//  double* Dap_a_alpha_RR = new double[ncabs1_vir1];
//  double* Dap_a_beta_RR = new double[ncabs2_vir2];
//  fill_n(Dap_a_alpha_RR, ncabs1_vir1, 0.0);
//  fill_n(Dap_a_beta_RR, ncabs2_vir2, 0.0);
//
//  ExEnv::out0() << endl << "D^a'_a RR part" << endl;
//  RR_sum_ijidx3_w_spin(cabs_vir_cabs,
//                       nspincases1, nspincases2, C_0, C_1,
//                       v_orbs1_ab, v_orbs2_ab,
//                       Dap_a_alpha_RR, Dap_a_beta_RR);
//
//  delete[] Dap_a_alpha_RR;
//  delete[] Dap_a_beta_RR;

//  if (this->r12eval()->coupling() == true) {
//    compute_Dapa(nspincases1, nspincase2, nocc_alpha, nocc_beta, C_0, C_1, Dapa[NSpinCases1]);
//  }

  // D[Alpha].print(prepend_spincase(Alpha,"mp2r12 one electron density:").c_str());

//  if (nspincases1 == 1) {
//     D_[Beta] = D_[Alpha];
//   }

//  density_diag_[Alpha] = D[Alpha];
//  density_diag_[Beta] = D[Beta];

  ExEnv::out0() << endl << "end of the computation of one electron density" << endl;
  ExEnv::out0() << endl << "**********************************************" << endl;
  ExEnv::out0() << endl << "**********************************************" << endl;
  return;
}


