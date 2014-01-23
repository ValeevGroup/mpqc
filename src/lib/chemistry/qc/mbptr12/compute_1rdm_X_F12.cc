/*
 * compute_1rdm_Z.cc
 *
 *  Created on: Mar 4, 2013
 *      Author: jinmei
 */

#include <cassert>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/wfn/spin.h>
#include <math/scmat/blas.h>
#include <math/scmat/local.h>
#include <util/misc/print.h>
#include <chemistry/qc/lcao/utils.h>

using namespace std;
using namespace sc;

enum AIidx_cases{ax_ix, xa_ix, ax_xi, xa_xi};
enum XorV_ints_cases{occ1occ2, occ2occ1, vir1occ2, occ1vir2};

namespace {
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

  void activate_ints(const std::string& bra1_id, const std::string& bra2_id,
                     const std::string& ket1_id, const std::string& ket2_id,
                     const std::string& descr_key, Ref<TwoBodyFourCenterMOIntsRuntime>& moints4_rtime,
                     Ref<DistArray4>& ints)
  {
    const std::string ints_key = ParsedTwoBodyFourCenterIntKey::key(bra1_id, bra2_id,
                                                                    ket1_id, ket2_id,
                                                                    descr_key,
                                                                    TwoBodyIntLayout::b1b2_k1k2);
    Ref<TwoBodyMOIntsTransform> ints_tform = moints4_rtime->get(ints_key);
    ints_tform->compute();
    ints = ints_tform->ints_distarray4();
    ints->activate();
  }

  // activate the three f12_ints for V or X
  void activate_ints_VorX(const Ref<OrbitalSpace>& bk1, const Ref<OrbitalSpace>& bk2,
                          const SpinCase2 spincase2, const Ref<R12IntEval>& r12eval,
                          const string& descr_f12_key, vector<Ref<DistArray4> >& ints)
  {
    Ref<R12WavefunctionWorld> r12world = r12eval->r12world();
    Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    const Ref<OrbitalSpace> occ1 = r12eval->occ(spin1);
    const Ref<OrbitalSpace> vir1 = r12eval->vir(spin1);
    const Ref<OrbitalSpace> orbs1 = r12eval->orbs(spin1);
    const Ref<OrbitalSpace> cabs1 = r12world->cabs_space(spin1);

    const Ref<OrbitalSpace> occ2 = r12eval->occ(spin2);
    const Ref<OrbitalSpace> vir2 = r12eval->vir(spin2);
    const Ref<OrbitalSpace> orbs2 = r12eval->orbs(spin2);
    const Ref<OrbitalSpace> cabs2 = r12world->cabs_space(spin2);

    // (ints)^b1b2_PQ or (ints)^PQ_k1k2
    Ref<DistArray4> p1p2_ints;
    activate_ints(bk1->id(), bk2->id(), orbs1->id(), orbs2->id(),
                  descr_f12_key, moints4_rtime, p1p2_ints);
    ints.push_back(p1p2_ints);

    // (ints)^b1b2_MA' or (ints)^MA'_k1k2
    Ref<DistArray4> i1a2_ints;
    activate_ints(bk1->id(), bk2->id(), occ1->id(), cabs2->id(),
                  descr_f12_key, moints4_rtime, i1a2_ints);
    ints.push_back(i1a2_ints);

    // (ints)^b1b2_A'M or (ints)^A'M_k1k2
    Ref<DistArray4> a1i2_ints;
    activate_ints(bk1->id(), bk2->id(), cabs1->id(), occ2->id(),
                  descr_f12_key, moints4_rtime, a1i2_ints);
    ints.push_back(a1i2_ints);
  }
  // end of function: activate_ints_X_f12

  // compute (f12f12)^b1b2_k1k2 or (f12/r12)^b1b2_k1k2 sum over index 3
  void compute_ints1ints2_sum3idx(const int b1b2_k1k2, const unsigned int ints12_idx,
                                  const Ref<DistArray4>& ints12, double* const result)
  {
    // index1, 2: a, i, index3: x(here: o)
    int index1 = 0;
    int index2 = 0;
    int index3 = 0;
    int size_idx1 = ints12->ni();
    int size_idx2 = ints12->nx();
    int size_idx3 = ints12->nj();
    // size of ket2
    const int size_k2 = ints12->ny();

    int* ints_idx1 = &index1;
    int* ints_idx2 = &index3;

    int* blk_idx1 = &index2;
    int* blk_idx2 = &index3;

    switch (b1b2_k1k2) {
    case ax_ix:
      // default value
    break;

    case xa_ix:
      size_idx1 = ints12->nj();
      size_idx3 = ints12->ni();
      ints_idx1 = &index3;
      ints_idx2 = &index1;
    break;

    case ax_xi:
      size_idx2 = ints12->ny();
      blk_idx1 = &index3;
      blk_idx2 = &index2;
    break;

    case xa_xi:
      size_idx1 = ints12->nj();
      size_idx2 = ints12->ny();
      size_idx3 = ints12->ni();
      ints_idx1 = &index3;
      ints_idx2 = &index1;
      blk_idx1 = &index3;
      blk_idx2 = &index2;
    break;

    default:
      ExEnv::out0() << "There is no such index";
      break;
    }

    double* iter_result = result;
    for(index1 = 0; index1 < size_idx1; ++index1) {
      for(index2 = 0; index2 < size_idx2; ++index2) {

        double sum_idx3 = 0;
        for(index3 = 0; index3 < size_idx3; ++index3) {
          const double* blk = ints12->retrieve_pair_block(*ints_idx1, *ints_idx2, ints12_idx);
          const int element = (*blk_idx1) * size_k2 + *blk_idx2;
          sum_idx3 += blk[element];

          ints12->release_pair_block(*ints_idx1, *ints_idx2, ints12_idx);
        }
        *iter_result = sum_idx3;
        ++iter_result;

      }
    }
  }
  // end of function: compute_ints1ints2_sum_3idx

  // compute g^{a idx1}_{i idx2} * d^idx2_idx1 or g^{idx1 a}_{i idx2} * d^idx2_idx1
  // d is CABS Singles 1e density
   void compute_GDai(const string ai_label,
                     const Ref<DistArray4>& ints, const unsigned int ints_idx,
                     const RefSCMatrix& D_cabs, const int idx1_offset, const int idx2_offset,
                     double* const result)
   {
     int size_a = ints->ni();
     int size_i = ints->nx();
     int size_idx1 = ints->nj(); // contracted index
     const blasint size_idx2 = ints->ny(); // size of second contracted index

     int a = 0;
     int idx1 = 0;
     int* ints_idx1 = &a;
     int* ints_idx2 = &idx1;

     if (ai_label == "ia") {
       size_a = ints->nj();
       size_idx1 = ints->ni();

       ints_idx1 = &idx1;
       ints_idx2 = &a;
     }

     const blasint one = 1;

     const int size_Drow = D_cabs.ncol();
     double* const D_row = new double[size_Drow];

     fill_n(result, size_a * size_i, 0.0);

     for(a = 0; a < size_a; ++a) {
       for(idx1 = 0; idx1 < size_idx1; ++idx1) {
         // a idx1 block or idx1 a block
         const double* g_blk = ints->retrieve_pair_block(*ints_idx1, *ints_idx2, ints_idx);

         RefSCVector D_row_idx1 = D_cabs.get_row(idx1+idx1_offset);
         fill_n(D_row, size_Drow, 0.0);
         D_row_idx1.convert(D_row);
         const double* const D_blk_idx2 = D_row + idx2_offset;

         const double* iter_g_blk = g_blk;
         double* iter_result = result + a * size_i;
         for(int i = 0; i < size_i; ++i) {
           const double gD = F77_DDOT(&size_idx2, iter_g_blk, &one, D_blk_idx2, &one);
           *iter_result += gD;

           ++iter_result;
           iter_g_blk += size_idx2;
         }

         ints->release_pair_block(*ints_idx1, *ints_idx2, ints_idx);
       }
     }
   }
   // end of function: compute_GDai

   // compute F^a_ixd1 * d^idx1_i
   // d is CABS Singles 1e density
    void compute_FDai(const RefSCMatrix& F, const RefSCMatrix& D_cabs,
                      const int size_i, const int idx1_offset,
                      double* const result)
    {
      int size_a = F.nrow();
      int size_idx1 = F.ncol(); // contracted index
      ExEnv::out0() << endl << "size_a: " << size_a
                    << "  size_idx1: " << size_idx1 << endl;

      const blasint one = 1;
//      const int size_Drow = D_cabs.ncol();
//      double* const F_blk = new double[size_idx1];
//      double* const D_row_array = new double[size_Drow];

      double* iter_result = result;
      for(int a = 0; a < size_a; ++a) {
//        ExEnv::out0() << "a: " << a << endl;
//        RefSCVector F_row_a = F.get_row(a);
//        fill_n(F_blk, size_idx1, 0.0);
//        F_row_a.convert(F_blk);

        for(int i = 0; i < size_i; ++i) {
//          RefSCVector D_row_i = D_cabs.get_row(i);
//          fill_n(D_row_array, size_Drow, 0.0);
//          D_row_i.convert(D_row_array);
//          const double* const D_blk = D_row_array + idx1_offset;
//
//          const double FD = F77_DDOT(&size_idx1, F_blk, &one, D_blk, &one);
//            *iter_result = FD;

          double sum_idx = 0.0;
          for (int idx1 = 0; idx1 < size_idx1; ++idx1) {
            sum_idx += F(a,idx1) * D_cabs(i, idx1+idx1_offset);
          }

          *iter_result = sum_idx;
          ++iter_result;
        }
      }

//      delete[] F_blk;
//      delete[] D_row_array;
    }
    // end of function: compute_FDai

    // compute F^i_p' * d^p'_a
    // d is CABS Singles 1e density
     void compute_FDai(const RefSCMatrix& F, const RefSCMatrix& D_cabs,
                       const int size_a, const int a_offset, const int idx1_offset,
                       double* const result)
     {
       int size_i = F.nrow();
       int size_idx1 = F.ncol(); // contracted index
       ExEnv::out0() << endl << "size_i: " << size_i
                     << "  size_idx1: " << size_idx1 << endl;

       const blasint one = 1;

       double* iter_result = result;
       for(int a = 0; a < size_a; ++a) {
         for(int i = 0; i < size_i; ++i) {

           double sum_idx = 0.0;
           for (int idx1 = 0; idx1 < size_idx1; ++idx1) {
             sum_idx += F(i,idx1) * D_cabs(a+a_offset,idx1+idx1_offset);
           }

           *iter_result = sum_idx;
           ++iter_result;
         }
       }

     }
     // end of function: compute_FDai

     // compute F * D in general form
     // D is CABS Singles 1e density
      void compute_FDai(const int size_a, const int size_i,
                        const string F_label, const RefSCMatrix& F,
                        const RefSCMatrix& D_cabs,
                        const int idx1_offset,
                        double* const result)
      {
        // get the right block
        int a = 0;
        int i = 0;
        int idx1 = 0; // contracted index

        // Select block indices
        // F^a_a' * D^a'_i or F^a_m * D^m_i
        int* F_idx1 = &a;
        int* F_idx2 = &idx1;

        int* d_idx1 = &idx1;
        int* d_idx2 = &i;

        int d_idx1_offset = idx1_offset;
        int d_idx2_offset = 0;

        int size_idx1 = F.ncol();

        if (F_label == "i") {
          // F^a'_i * D^a_a'
          F_idx1 = &idx1;
          F_idx2 = &i;

          d_idx1 = &a;
          d_idx2 = &idx1;
          d_idx1_offset = F.ncol(); // # of occ
          d_idx2_offset = idx1_offset;
          size_idx1 = F.nrow();
        }

//        ExEnv::out0() << endl << "  size_idx1: " << size_idx1
//            << " d_idx1_offset: " << d_idx1_offset
//            <<" d_idx2_offset: " << d_idx2_offset << endl;

        const blasint one = 1;
        double* iter_result = result;
        for(a = 0; a < size_a; ++a) {
          for(i = 0; i < size_i; ++i) {

            double sum_idx = 0.0;
            for (idx1 = 0; idx1 < size_idx1; ++idx1) {
              sum_idx += F(*F_idx1, *F_idx2)
                         * D_cabs(*d_idx1+d_idx1_offset, *d_idx2+d_idx2_offset);
            }

            *iter_result = sum_idx;
            ++iter_result;
          }
        }
      }
      // end of function: compute_FDai

   // compute one-electron density form CABS contribution
   RefSCMatrix compute_D_CABS(SpinCase1 spin, const Ref<R12IntEval>& r12eval,
                              const Ref<R12EnergyIntermediates>& r12intermediates)
   {
     Ref<R12WavefunctionWorld> r12world = r12eval->r12world();

     const Ref<OrbitalSpace> occ = r12eval->occ(spin);
     const Ref<OrbitalSpace> vir = r12eval->vir(spin);
     const Ref<OrbitalSpace> cabs = r12world->cabs_space(spin);

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

     const RefSCMatrix T1_cabs = r12eval->T1_cabs(spin);

      RefSCDimension rowdim_occ = new SCDimension(nocc);
      RefSCDimension coldim_vir_com = new SCDimension(nvir_com);
      RefSCMatrix T = localkit->matrix(rowdim_occ, coldim_vir_com);
      T.assign(0.0);

#define TEST_CC_CABS  0
#if TEST_CC_CABS
      if (r12intermediates->T2_cc_computed()) {
        RefSCMatrix T1_ccsd = r12intermediates->get_T1_cc(spin);
        const Ref<OrbitalSpace> occ_act = r12eval->occ_act(spin);
        const int nocc_act = occ_act->rank();
        const int nfzc = nocc - nocc_act;

        T.accumulate_subblock(T1_ccsd, nfzc, nocc-1, 0, nvir-1, 0, 0);
      }
#endif
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

      for (unsigned int irrep = 1; irrep < nirreps; ++irrep) {
        occoff[irrep] = occoff[irrep-1] + occpi[irrep-1];
        viroff[irrep] = viroff[irrep-1] + virpi[irrep-1];
        cabsoff[irrep] = cabsoff[irrep-1] + cabspi[irrep-1];
        vir_comoff[irrep] = viroff[irrep] + cabsoff[irrep];
      }

      for (int h = 0; h < nirreps; ++h) {

#if TEST_CC_CABS
        if (r12intermediates->T2_cc_computed()) {
          for (int i = 0; i < nocc; ++i) {
            for (int ap = 0; ap < ncabs; ++ap) {
              T.set_element(i,nvir+ap,T1_cabs.get_element(i,ap));
            }
          }

        } else {
#endif
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
#if TEST_CC_CABS
        }
#endif
      }

#if TEST_CC_CABS
     if (!r12intermediates->T2_cc_computed()) {
#endif
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
#if TEST_CC_CABS
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
#endif
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

#undef TEST_CC_CABS
     return D;
   }
   // end of compute_D_CABS

} // end of namespace

// compute (rr)^a_i or (gr)^a_i contracted over 3 indices:
// o, w, and index 3 or c, d, and index 3
void compute_ai_sum_3idx(const int b1b2_k1k2,
                         const unsigned int ints1_idx, const unsigned int ints2_idx,
                         const Ref<DistArray4>& ints1, const Ref<DistArray4>& ints2,
                         double* const array_AI)
{
  // idx1:a, idx2:i, idx3 is the other contracted index
  // eg: R^ow_ac' g^ic'_ow: idx1 is a, idx2 is i, and idx3 is c'
  int idx1 = 0;
  int idx2 = 0;
  int idx3 = 0;

  int size_idx1 = ints1->ni();
  int size_idx2 = ints2->ni();
  int size_idx3 = ints1->nj();
  // get the number of the contracted indices
  // eg: R^ow_ac' g^ic'_ow: sum_idx is ow
  const blasint size_sum_idx = ints1->nx() * ints1->ny();

  // get block indices
  int* blk1_idx1 = &idx1;
  int* blk1_idx2 = &idx3;
  int* blk2_idx1 = &idx2;
  int* blk2_idx2 = &idx3;

  switch (b1b2_k1k2) {
  case ax_ix:
    // default value
  break;

  case xa_ix:
    size_idx1 = ints1->nj();
    size_idx3 = ints1->ni();
    blk1_idx1 = &idx3;
    blk1_idx2 = &idx1;
  break;

  case ax_xi:
    size_idx2 = ints2->nj();
    blk2_idx1 = &idx3;
    blk2_idx2 = &idx2;
  break;

  case xa_xi:
    size_idx1 = ints1->nj();
    size_idx2 = ints2->nj();
    size_idx3 = ints1->ni();
    blk1_idx1 = &idx3;
    blk1_idx2 = &idx1;
    blk2_idx1 = &idx3;
    blk2_idx2 = &idx2;
  break;

  default:
    ExEnv::out0() << "There is no such index case";
    break;
  }

  const blasint one = 1; // for F77_DDOT
  double* iter_array = array_AI;

    for (idx1 = 0; idx1 < size_idx1; ++idx1) {
      for(idx2 = 0; idx2 < size_idx2; ++idx2) {

        double sum_idx3 = 0;
        for(idx3 = 0; idx3 < size_idx3; ++idx3) {
          const double* blk1 = ints1->retrieve_pair_block((*blk1_idx1), (*blk1_idx2), ints1_idx);
          const double* blk2 = ints2->retrieve_pair_block((*blk2_idx1), (*blk2_idx2), ints2_idx);

          const double sum_idx = F77_DDOT(&size_sum_idx, blk1, &one, blk2, &one);

          sum_idx3 += sum_idx;

          ints1->release_pair_block((*blk1_idx1), (*blk1_idx2), ints1_idx);
          ints2->release_pair_block((*blk2_idx1), (*blk2_idx2), ints2_idx);
        }
        *iter_array = sum_idx3;
        ++iter_array;
      }
    }
}
// end of compute_ai_sum_3dix

// compute R^SigmaDelta_b1b2 R^k1k2_SigmaDelta
// or g^SigmaDelta_b1b2 R^SigmaDelta_ab
// sum over complete virtual orbitals a, b, and another index
// index 1,2: a, i; index 3: o    eg: g^SigmaDelta_ao * R^SigmaDelta_io
void compute_ai_sum_abx(const int b1b2_k1k2, const int ints1ints2_idx,
                        const int ints1_idx, const int ints2_idx,
                        const Ref<DistArray4>& ints1ints2,
                        const vector< Ref<DistArray4> >& v_ints1,
                        const vector< Ref<DistArray4> >& v_ints2,
                        double* const result)
{
  int nvir = ints1ints2->ni();
  int nocc = ints1ints2->nx();

  if (b1b2_k1k2 == xa_ix || b1b2_k1k2 == xa_xi) {
    nvir = ints1ints2->nj();
  }
  if (b1b2_k1k2 == ax_xi || b1b2_k1k2 == xa_xi) {
    nocc = ints1ints2->ny();
  }
  const int nvir_occ = nvir * nocc;

  // (ints1)^ab_b1b2 (ints2)^k1k2 _ab
  // =  (ints1ints2)^b1b2_k1k2
  //  - (ints1)^b1b2_pq (ints2)^pq_k1k2
  //  - (ints1)^b1b2_pq_ma' (ints2)^ma'_k1k2
  //  - (ints1)^b1b2_pq_a'm (ints2)^a'm_k1k2

  // add (ints1ints2)^b1b2_k1k2
  compute_ints1ints2_sum3idx(b1b2_k1k2, ints1ints2_idx,
                             ints1ints2, result);

  if (v_ints1.size() != 3)
    ExEnv::out0() << "Then number of integrals in computing (ints)^ab_b1b2 (ints)^k1k2 _ab are wrong"
                  << endl;

  // subtract (ints)^b1b2_pq (ints)^pq_k1k2,
  //          (ints)^b1b2_pq_ma' (ints)^ma'_k1k2,
  //      and (ints)^b1b2_pq_a'm (ints)^a'm_k1k2
  for (int i = 0; i != 3; ++i) {
    const Ref<DistArray4> ints1 = v_ints1[i];
    const Ref<DistArray4> ints2 = v_ints2[i];

    double* XR_part = new double[nvir_occ];
    fill_n(XR_part, nvir_occ, 0.0);
    compute_ai_sum_3idx(b1b2_k1k2, ints1_idx, ints2_idx,
                        ints1, ints2, XR_part);

    double* iter_result = result;
    const double* iter_XR_part = XR_part;
    for (int a = 0; a != nvir; ++a) {
      for (int i = 0; i != nocc; ++i) {
        *iter_result -= *iter_XR_part;

        ++iter_result;
        ++iter_XR_part;
      }
    }

    delete[] XR_part;
    XR_part = NULL;
  }
}
// end of compute_ai_sum_abj

// compute 1/2 \bar{\tilde{R}}_{ac' }^{ow}\bar{g}_{ow}^{ic'}
void compute_RGow(SpinCase1 spin, const double C_0, const double C_1,
                  const Ref<R12IntEval>& r12eval, int debug,
                  double* const Xai_RGow) {

  Ref<R12WavefunctionWorld> r12world = r12eval->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);

  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const TwoBodyOper::type eri_type = r12world->r12tech()->corrfactor()->tbint_type_eri();
  const unsigned int f12_idx = descr_f12->intset(f12_type);
  const unsigned int eri_idx = descr_f12->intset(eri_type);

  const int nspincases2 = (r12eval->spin_polarized() ? 3 : 2);

  const SpinCase1 spin1 = Alpha;
  const SpinCase1 spin2 = Beta;
  const Ref<OrbitalSpace> occ1_act = r12eval->occ_act(spin1);
  const Ref<OrbitalSpace> occ2_act = r12eval->occ_act(spin2);
  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();

  const Ref<OrbitalSpace> occ1 = r12eval->occ(spin1);
  const Ref<OrbitalSpace> occ2 = r12eval->occ(spin2);
  const int nocc1 = occ1->rank();
  const int nocc2 = occ2->rank();

  const Ref<OrbitalSpace> occ_act = (spin == Alpha? occ1_act : occ2_act);
  const int nocc_act = (spin == Alpha? nocc1_act : nocc2_act);

  const Ref<OrbitalSpace> occ = (spin == Alpha? occ1 : occ2);
  const int nocc = (spin == Alpha? nocc1 : nocc2);

  const Ref<OrbitalSpace> vir = r12eval->vir(spin);
  const Ref<OrbitalSpace> cabs = r12world->cabs_space(spin);
  const int nvir = vir->rank();

  const int nvir_occ = nvir * nocc;

  // calculate AlphaAlpha/BetaBeta part:
  // R^ow_ac' g^ic'_ow & R^ow_c'a g^ic'_ow
  double* Rg1 = new double[nvir_occ];
  double* Rg2 = new double[nvir_occ];
  fill_n(Rg1, nvir_occ, 0.0);
  fill_n(Rg2, nvir_occ, 0.0);

  // activate integrals
  Ref<DistArray4> acp_ii_ints;
  Ref<DistArray4> cpa_ii_ints;
  Ref<DistArray4> icp_ii_ints;

  activate_ints(vir->id(), cabs->id(), occ_act->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, acp_ii_ints);
  activate_ints(cabs->id(), vir->id(), occ_act->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, cpa_ii_ints);
  activate_ints(occ->id(), cabs->id(), occ_act->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, icp_ii_ints);

  // R^ow_ac' g^ic'_ow
  compute_ai_sum_3idx(ax_ix, f12_idx, eri_idx,
                      acp_ii_ints, icp_ii_ints,
                      Rg1);

  // R^ow_c'a g^ic'_ow
  compute_ai_sum_3idx(xa_ix, f12_idx, eri_idx,
                      cpa_ii_ints, icp_ii_ints,
                      Rg2);

  acp_ii_ints->deactivate();
  cpa_ii_ints->deactivate();
  icp_ii_ints->deactivate();

  if (debug >= DefaultPrintThresholds::allN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    print_intermediate(spinletter, "R^ow_ac' g^ic'_ow", Rg1, nvir, nocc);
    print_intermediate(spinletter, "R^ow_c'a g^ic'_ow", Rg2, nvir, nocc);
  }

  // calculate AlphaBeta part:
  //             alpha case                   beta case
  // Rg1_ab: R^o1w2_a1c'2 g^i1c'2_o1w2 or R^o1w2_c'2a1 g^i1c'2_o1w2
  // Rg2_ab: R^o1w2_c'1a2 g^c'1i2_o1w2 or R^o1w2_a2c'1 g^c'1i2_o1w2

  double* Rg1_ab = NULL;
  double* Rg2_ab = NULL;
  if (nspincases2 == 3) {

     const Ref<OrbitalSpace> vir1 = r12eval->vir(spin1);
     const Ref<OrbitalSpace> cabs1 = r12world->cabs_space(spin1);
     const Ref<OrbitalSpace> vir2 = r12eval->vir(spin2);
     const Ref<OrbitalSpace> cabs2 = r12world->cabs_space(spin2);

     Rg1_ab = new double[nvir_occ];
     Rg2_ab = new double[nvir_occ];
     fill_n(Rg1_ab, nvir_occ, 0.0);
     fill_n(Rg2_ab, nvir_occ, 0.0);

     if (spin == Alpha) {

       Ref<DistArray4> acpii_1212_ints;   // R^o1w2_a1c'2
       Ref<DistArray4> cpaii_2112_ints;   // R^o1w2_c'2a1
       Ref<DistArray4> icpii_1212_ints;   // g^i1c'2_o1w2

       activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, acpii_1212_ints);
       activate_ints(cabs2->id(), vir1->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, cpaii_2112_ints);
       activate_ints(occ1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, icpii_1212_ints);


       // Rg1_ab: R^o1w2_a1c'2 g^i1c'2_o1w2
       compute_ai_sum_3idx(ax_ix, f12_idx, eri_idx,
                           acpii_1212_ints, icpii_1212_ints,
                           Rg1_ab);
       // Rg2_ab: R^o1w2_c'2a1 g^i1c'2_o1w2
       compute_ai_sum_3idx(xa_ix, f12_idx, eri_idx,
                           cpaii_2112_ints, icpii_1212_ints,
                           Rg2_ab);

       acpii_1212_ints->deactivate();
       cpaii_2112_ints->deactivate();
       icpii_1212_ints->deactivate();

       if (debug >= DefaultPrintThresholds::allN2) {
         print_intermediate("Alpha", "R^o1w2_a1c'2 g^i1c'2_o1w2", Rg1_ab, nvir, nocc);
         print_intermediate("Alpha", "R^o1w2_c'2a1 g^i1c'2_o1w2", Rg2_ab, nvir, nocc);
       }
     } else {
         // R^o1w2_c'1a2 g^c'1i2_o1w2 or R^o1w2_a2c'1 g^c'1i2_o1w2
         Ref<DistArray4> cpaii_1212_ints;   // R^o1w2_c'1a2
         Ref<DistArray4> acpii_2112_ints;   // R^o1w2_a2c'1
         Ref<DistArray4> cpiii_1212_ints;   // g^c'1i2_o1w2

         activate_ints(cabs1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, cpaii_1212_ints);
         activate_ints(vir2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, acpii_2112_ints);
         activate_ints(cabs1->id(), occ2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, cpiii_1212_ints);


         // R^o1w2_c'1a2 g^c'1i2_o1w2
         compute_ai_sum_3idx(xa_xi, f12_idx, eri_idx,
                             cpaii_1212_ints, cpiii_1212_ints,
                             Rg1_ab);
         // R^o1w2_a2c'1 g^c'1i2_o1w2
         compute_ai_sum_3idx(ax_xi, f12_idx,eri_idx,
                             acpii_2112_ints, cpiii_1212_ints,
                             Rg2_ab);

         cpaii_1212_ints->deactivate();
         acpii_2112_ints->deactivate();
         cpiii_1212_ints->deactivate();

         if (debug >= DefaultPrintThresholds::allN2) {
           print_intermediate("Beta", "R^o1w2_c'1a2 g^c'1i2_o1w2", Rg1_ab, nvir, nocc);
           print_intermediate("Beta", "R^o1w2_a2c'1 g^c'1i2_o1w2", Rg2_ab, nvir, nocc);
         }
     }
  } else if (nspincases2 == 2) {
       Rg1_ab = Rg1;
       Rg2_ab = Rg2;
  }
  // end of AlphaBeta part

  // calculate Rg^a_i
  const double* iter_Rg1 = Rg1;
  const double* iter_Rg2 = Rg2;
  const double* iter_Rg1_ab = Rg1_ab;
  const double* iter_Rg2_ab = Rg2_ab;

  double* iter_Xai = Xai_RGow;
   for (int a = 0;  a < nvir; ++a) {
     for (int i = 0; i < nocc; ++i) {

       // AlphaAlpha/BetaBeta part
       double x_12 = C_1 * (*iter_Rg1 - *iter_Rg2);
       ++iter_Rg1;
       ++iter_Rg2;

       // AlphaBeta part
       if (nocc1 != 0 && nocc2 != 0) {
         x_12 += 0.5 * (C_0 + C_1) * (*iter_Rg1_ab)
                + 0.5 * (C_0 - C_1) * (*iter_Rg2_ab);
         ++iter_Rg1_ab;
         ++iter_Rg2_ab;
       } // end of AlphaBeta part

       *iter_Xai = x_12;
       ++iter_Xai;
     }
   } // end of calculating Xai

   delete[] Rg1;
   delete[] Rg2;
   if (nspincases2 == 3) {
     delete[] Rg1_ab;
     delete[] Rg2_ab;
   }

}
// end of compute_RGow

// compute 1/2 \bar{\tilde{R}}^oi}_{\sigma\delta} \bar{g}^{\sigma\delta}_{ao}
void compute_RG_SigmaDelta(SpinCase1 spin, const double C_0, const double C_1,
                  const Ref<R12IntEval>& r12eval, int debug,
                  double* const Xai_RG_SigmaDelta) {

  Ref<R12WavefunctionWorld> r12world = r12eval->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);

  const TwoBodyOper::type f12eri_type = r12world->r12tech()->corrfactor()->tbint_type_f12eri();
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const TwoBodyOper::type eri_type = r12world->r12tech()->corrfactor()->tbint_type_eri();
  const unsigned int f12eri_idx = descr_f12->intset(f12eri_type);
  const unsigned int f12_idx = descr_f12->intset(f12_type);
  const unsigned int eri_idx = descr_f12->intset(eri_type);

  const int nspincases2 = (r12eval->spin_polarized() ? 3 : 2);

  const SpinCase1 spin1 = Alpha;
  const SpinCase1 spin2 = Beta;
  const Ref<OrbitalSpace> occ1_act = r12eval->occ_act(spin1);
  const Ref<OrbitalSpace> occ2_act = r12eval->occ_act(spin2);
  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();

  const Ref<OrbitalSpace> occ_act = (spin == Alpha? occ1_act : occ2_act);
  const int nocc_act = (spin == Alpha? nocc1_act : nocc2_act);

  const Ref<OrbitalSpace> vir = r12eval->vir(spin);
  const int nvir = vir->rank();
  const int nvir_occ_act = nvir * nocc_act;

  // calculate AlphaAlpha/BetaBeta part:
  //   R^io_SigmaDelta g^SigmaDelta_ao
  // & R^oi_SigmaDelta g^SigmaDelta_ao

  double* Rg1 = new double[nvir_occ_act];
  double* Rg2 = new double[nvir_occ_act];
  fill_n(Rg1, nvir_occ_act, 0.0);
  fill_n(Rg2, nvir_occ_act, 0.0);

  // g^SigmaDelta_ao R^io_SigmaDelta
  //  = (f12/r12)_ao^io - g_ao^pq f12_pq^io
  //  - g_ao^ma' f12_ma'^io - g_ao^a'm f12_a'm_io

  // activate integrals
  Ref<DistArray4> aiii_ints;
  vector<Ref<DistArray4> > v_ai_ints;
  vector<Ref<DistArray4> > v_ii_ints;

  activate_ints(vir->id(), occ_act->id(), occ_act->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, aiii_ints);
  const SpinCase2 spincase2 = (spin == Alpha? AlphaAlpha: BetaBeta);
  activate_ints_VorX(vir, occ_act, spincase2,
                     r12eval, descr_f12_key, v_ai_ints);
  activate_ints_VorX(occ_act, occ_act, spincase2,
                     r12eval, descr_f12_key, v_ii_ints);

  // g^SigmaDelta_ao R^io_SigmaDelta
  compute_ai_sum_abx(ax_ix, f12eri_idx, eri_idx, f12_idx,
                     aiii_ints, v_ai_ints, v_ii_ints, Rg1);
  // g^SigmaDelta_ao R^oi_SigmaDelta
  compute_ai_sum_abx(ax_xi, f12eri_idx, eri_idx, f12_idx,
                     aiii_ints, v_ai_ints, v_ii_ints, Rg2);

  aiii_ints->deactivate();
  for (int i = 0; i != v_ai_ints.size(); ++i) {
    v_ai_ints[i]->deactivate();
    v_ii_ints[i]->deactivate();
  }


  if (debug >= DefaultPrintThresholds::allN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    print_intermediate(spinletter, "g^SigmaDelta_ao R^io_SigmaDelta", Rg1, nvir, nocc_act);
    print_intermediate(spinletter, "g^SigmaDelta_ao R^oi_SigmaDelta", Rg2, nvir, nocc_act);
  }

    // test for Xai_RG_SigmaDelta
    {
#if 0
      RefSCDimension rowdim = new SCDimension(nvir);
      RefSCDimension coldim = new SCDimension(nocc);
      Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
      RefSCMatrix Xai = localkit->matrix(rowdim, coldim);
      Xai.assign(0.0);

      const Ref<OrbitalSpace>& orbs = r12eval->orbs(spin);
      // Vpq_ij
      const SpinCase2 spincase2 = (spin == Alpha? AlphaAlpha : BetaBeta);
      std::vector< Ref<DistArray4> > Vpq_vec = r12eval->V_distarray4(spincase2, orbs, orbs);
      MPQC_ASSERT(Vpq_vec.size() == 1);
      Ref<DistArray4> Vpq_ij = Vpq_vec[0];
      Ref<DistArray4> Vaj_ij;
      sc::map(Vpq_ij, occ, occ, orbs, orbs, Vaj_ij, occ, occ, vir, occ);

      Vaj_ij->activate();

      for (int a = 0;  a < nvir; ++a) {
        for (int i = 0; i < nocc; ++i) {

          double Va_i = 0.0;
          for (int j = 0; j < nocc; ++j) {
            const double* const Vij_blk = Vaj_ij->retrieve_pair_block(i, j, 0);
            const int aj_ij = a * nocc + j;
            Va_i += Vij_blk[aj_ij];

            Vaj_ij->release_pair_block(i, j, 0);
          }

          Xai.set_element(a, i, Va_i);
        }
      }
      Vaj_ij->deactivate();
      Xai.print(prepend_spincase(spin,"Rg_SigmaDelta from Vai_ij:").c_str());

      double* Rg12 = new double[nvir_occ];
      fill_n(Rg12, nvir_occ, 0.0);

      double* iter_Rg12 = Rg12;
      double* iter_Rg1 = Rg1;
      double* iter_Rg2 = Rg2;
      for (int a = 0;  a < nvir; ++a) {
        for (int i = 0; i < nocc; ++i) {
          *iter_Rg12 = *iter_Rg1 - *iter_Rg2;

          ++iter_Rg12;
          ++iter_Rg1;
          ++iter_Rg2;
        }
      }
      const string spinletter = (spin == Alpha? "AlphaAlpha" : "BetaBeta");
      print_intermediate(spinletter, "antisymmetrized Rg_SigmaDelta", Rg12, nvir, nocc);
      delete[] Rg12;
#endif
    } // end of test

  // calculate AlphaBeta part:
  // Rg1_ab alpha: g^Sigma1Delta1_a1o2 R^i1o2_Sigma1Delta2
  //        beta : g^Sigma1Delta1_a1o2 R^o2i1_Sigma1Delta2
  // Rg2_ab alpha: g^Sigma1Delta1_o1a2 R^o1i2_Sigma1Delta2
  //        beta : g^Sigma1Delta1_o1a2 R^i2o1_Sigma1Delta2

  double* Rg1_ab = NULL;
  double* Rg2_ab = NULL;
  if (nspincases2 == 3) {

     const Ref<OrbitalSpace> vir1 = r12eval->vir(spin1);
     const Ref<OrbitalSpace> vir2 = r12eval->vir(spin2);

     vector<Ref<DistArray4> > v_i1i2_ints;
     vector<Ref<DistArray4> > v_i2i1_ints;

     activate_ints_VorX(occ1_act, occ2_act, AlphaBeta,
                        r12eval, descr_f12_key, v_i1i2_ints);
     activate_ints_VorX(occ2_act, occ1_act, AlphaBeta,
                        r12eval, descr_f12_key, v_i2i1_ints);

     Rg1_ab = new double[nvir_occ_act];
     Rg2_ab = new double[nvir_occ_act];
     fill_n(Rg1_ab, nvir_occ_act, 0.0);
     fill_n(Rg2_ab, nvir_occ_act, 0.0);

     if (spin == Alpha) {

       Ref<DistArray4> a1i2i1i2_ints;
       Ref<DistArray4> a1i2i2i1_ints;
       vector<Ref<DistArray4> > v_a1i2_ints;

       activate_ints(vir1->id(), occ2_act->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, a1i2i1i2_ints);
       activate_ints(vir1->id(), occ2_act->id(), occ2_act->id(), occ1_act->id(),
                     descr_f12_key, moints4_rtime, a1i2i2i1_ints);
       activate_ints_VorX(vir1, occ2_act, AlphaBeta,
                          r12eval, descr_f12_key, v_a1i2_ints);

       // g^Sigma1Delta2_a1o2 R^i1o2_Sigma1Delta2
       compute_ai_sum_abx(ax_ix, f12eri_idx, eri_idx, f12_idx,
                          a1i2i1i2_ints, v_a1i2_ints, v_i1i2_ints, Rg1_ab);
       // g^Sigma1Delta2_a1o2 R^o2i1_Sigma1Delta2
       compute_ai_sum_abx(ax_xi, f12eri_idx, eri_idx, f12_idx,
                          a1i2i2i1_ints, v_a1i2_ints, v_i2i1_ints, Rg2_ab);

       a1i2i1i2_ints->deactivate();
       a1i2i2i1_ints->deactivate();
       for (int i = 0; i != v_a1i2_ints.size(); ++i) {
         v_a1i2_ints[i]->deactivate();
       }

       if (debug >= DefaultPrintThresholds::allN2) {
         print_intermediate("Alpha", "g^Sigma1Delta1_a1o2 R^i1o2_Sigma1Delta2", Rg1_ab, nvir, nocc_act);
         print_intermediate("Alpha", "g^Sigma1Delta2_a1o2 R^o2i1_Sigma1Delta2", Rg2_ab, nvir, nocc_act);
       }

       // test for g^Sigma1Delta1_a1o2 R^i1o2_Sigma1Delta2
       {
   #if 0
         RefSCDimension rowdim = new SCDimension(nvir);
         RefSCDimension coldim = new SCDimension(nocc);
         Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
         RefSCMatrix Xai = localkit->matrix(rowdim, coldim);
         Xai.assign(0.0);

         const Ref<OrbitalSpace>& orbs1 = r12eval->orbs(spin1);
         const Ref<OrbitalSpace>& orbs2 = r12eval->orbs(spin2);
         // Vpq_ij
         //const SpinCase2 spincase2 = static_cast<SpinCase2>(spin+1);
         std::vector< Ref<DistArray4> > Vpq_vec = r12eval->V_distarray4(AlphaBeta, orbs1, orbs2);
         MPQC_ASSERT(Vpq_vec.size() == 1);
         Ref<DistArray4> Vpq_ij = Vpq_vec[0];
         Ref<DistArray4> Vaj_ij;
         sc::map(Vpq_ij, occ1, occ2, orbs1, orbs2, Vaj_ij, occ1, occ2, vir1, occ2);

         Vaj_ij->activate();

         for (int a = 0;  a < nvir; ++a) {
           for (int i = 0; i < nocc; ++i) {

             double Va_i = 0.0;
             for (int j = 0; j < nocc2; ++j) {
               const double* const Vij_blk = Vaj_ij->retrieve_pair_block(i, j, 0);
               const int aj_ij = a * nocc2 + j;
               Va_i += Vij_blk[aj_ij];

               Vaj_ij->release_pair_block(i, j, 0);
             }

             Xai.set_element(a, i, Va_i);
           }
         }
         Vaj_ij->deactivate();
         Xai.print(prepend_spincase(spin,"Rg1_ab_SigmaDelta from Vai_ij:").c_str());
   #endif
       } // end of test
     } else {
         Ref<DistArray4> i1a2i1i2_ints;
         Ref<DistArray4> i1a2i2i1_ints;
         vector<Ref<DistArray4> > v_i1a2_ints;

         activate_ints(occ1_act->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, i1a2i1i2_ints);
         activate_ints(occ1_act->id(), vir2->id(), occ2_act->id(), occ1_act->id(),
                       descr_f12_key, moints4_rtime, i1a2i2i1_ints);
         activate_ints_VorX(occ1_act, vir2, AlphaBeta,
                            r12eval, descr_f12_key, v_i1a2_ints);

         // g^Sigma1Delta1_o1a2 R^o1i2_Sigma1Delta2
         compute_ai_sum_abx(xa_xi, f12eri_idx, eri_idx, f12_idx,
                            i1a2i1i2_ints, v_i1a2_ints, v_i1i2_ints, Rg1_ab);
         // g^Sigma1Delta1_o1a2 R^i2o1_Sigma1Delta2
         compute_ai_sum_abx(xa_ix, f12eri_idx, eri_idx, f12_idx,
                            i1a2i2i1_ints, v_i1a2_ints, v_i2i1_ints, Rg2_ab);

         i1a2i1i2_ints->deactivate();
         i1a2i2i1_ints->deactivate();
         for (int i = 0; i != v_i1a2_ints.size(); ++i) {
           v_i1a2_ints[i]->deactivate();
         }

         if (debug >= DefaultPrintThresholds::allN2) {
           print_intermediate("Beta", "g^Sigma1Delta1_o1a2 R^o1i2_Sigma1Delta2", Rg1_ab, nvir, nocc_act);
           print_intermediate("Beta", "g^Sigma1Delta1_o1a2 R^i2o1_Sigma1Delta2", Rg2_ab, nvir, nocc_act);
         }
         // test for g^Sigma1Delta1_o1a2 R^o1i2_Sigma1Delta2
         {
     #if 0
           RefSCDimension rowdim = new SCDimension(nvir);
           RefSCDimension coldim = new SCDimension(nocc);
           Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
           RefSCMatrix Xai = localkit->matrix(rowdim, coldim);
           Xai.assign(0.0);

           const int nvir2 = vir2->rank();

           const Ref<OrbitalSpace>& orbs1 = r12eval->orbs(spin1);
           const Ref<OrbitalSpace>& orbs2 = r12eval->orbs(spin2);
           // Vpq_ij
           //const SpinCase2 spincase2 = static_cast<SpinCase2>(spin+1);
           std::vector< Ref<DistArray4> > Vpq_vec = r12eval->V_distarray4(AlphaBeta, orbs1, orbs2);
           MPQC_ASSERT(Vpq_vec.size() == 1);
           Ref<DistArray4> Vpq_ij = Vpq_vec[0];
           Ref<DistArray4> Vaj_ij;
           sc::map(Vpq_ij, occ1, occ2, orbs1, orbs2, Vaj_ij, occ1, occ2, occ1, vir2);

           Vaj_ij->activate();

           for (int a = 0;  a < nvir; ++a) {
             for (int i = 0; i < nocc; ++i) {

               double Va_i = 0.0;
               for (int j = 0; j < nocc1; ++j) {
                 const double* const Vij_blk = Vaj_ij->retrieve_pair_block(j, i, 0);
                 const int ja_ji = j * nvir2 + a;
                 Va_i += Vij_blk[ja_ji];

                 Vaj_ij->release_pair_block(j, i, 0);
               }

               Xai.set_element(a, i, Va_i);
             }
           }
           Vaj_ij->deactivate();
           Xai.print(prepend_spincase(spin,"Rg1_ab_SigmaDelta from Vai_ij:").c_str());
     #endif
         } // end of test
     }
     for (int i = 0; i != v_i1i2_ints.size(); ++i) {
       v_i1i2_ints[i]->deactivate();
       v_i2i1_ints[i]->deactivate();
     }
  } else if (nspincases2 == 2) {
       Rg1_ab = Rg1;
       Rg2_ab = Rg2;
  }
  // end of AlphaBeta part

  // calculate Rg^a_i
  const double* iter_Rg1 = Rg1;
  const double* iter_Rg2 = Rg2;
  const double* iter_Rg1_ab = Rg1_ab;
  const double* iter_Rg2_ab = Rg2_ab;

  double* iter_Xai = Xai_RG_SigmaDelta;
   for (int a = 0;  a < nvir; ++a) {
     for (int i = 0; i < nocc_act; ++i) {

       // AlphaAlpha/BetaBeta part
       double x_12 = C_1 * (*iter_Rg1 - *iter_Rg2);
       ++iter_Rg1;
       ++iter_Rg2;

       // AlphaBeta part
       if (nocc1_act != 0 && nocc2_act != 0) {
         x_12 += 0.5 * (C_0 + C_1) * (*iter_Rg1_ab)
                + 0.5 * (C_0 - C_1) * (*iter_Rg2_ab);
         ++iter_Rg1_ab;
         ++iter_Rg2_ab;
       } // end of AlphaBeta part

       *iter_Xai = x_12;
       ++iter_Xai;
     }
   } // end of calculating Xai

   delete[] Rg1;
   delete[] Rg2;
   if (nspincases2 == 3) {
     delete[] Rg1_ab;
     delete[] Rg2_ab;
   }

}
// end of compute_RG_SigmaDelta

// compute 1/2 \bar{\tilde{R}}^{ab'}_{ow}
//          *  F^{alpha}_i
//          *  \bar{\tilde{R}}^{ow}_{\alpha b'}
void compute_RFR_ow(SpinCase1 spin, const string orbs_focc,
                    const double C_0, const double C_1,
                    const Ref<R12IntEval>& r12eval, int debug,
                    double* const Xai_RFR_ow) {

  Ref<R12WavefunctionWorld> r12world = r12eval->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const unsigned int f12_idx = descr_f12->intset(f12_type);

  const int nspincases2 = (r12eval->spin_polarized() ? 3 : 2);

  const SpinCase1 spin1 = Alpha;
  const SpinCase1 spin2 = Beta;

  const Ref<OrbitalSpace> occ1_act = r12eval->occ_act(spin1);
  const Ref<OrbitalSpace> occ2_act = r12eval->occ_act(spin2);
  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();

  const Ref<OrbitalSpace> occ1 = r12eval->occ(spin1);
  const Ref<OrbitalSpace> occ2 = r12eval->occ(spin2);
  const int nocc1 = occ1->rank();
  const int nocc2 = occ2->rank();

  const Ref<OrbitalSpace> occ_act = (spin == Alpha? occ1_act : occ2_act);
  const Ref<OrbitalSpace> occ = (spin == Alpha? occ1 : occ2);
  const int nocc = (spin == Alpha? nocc1 : nocc2);

  const Ref<OrbitalSpace> vir = r12eval->vir(spin);
  const Ref<OrbitalSpace> cabs = r12world->cabs_space(spin);

  assert (orbs_focc == "vir" || orbs_focc == "cabs");
  Ref<OrbitalSpace> focc;
  if (orbs_focc == "vir") {
    focc = r12eval->F_m_a(spin);
  } else if (orbs_focc == "cabs") {
      focc = r12eval->F_m_A(spin);
  }

  const int nvir = vir->rank();
  const int nvir_occ = nvir * nocc;

  // calculate AlphaAlpha/BetaBeta part:
  // R^ab'_ow R^ow_fi_alpha b' & R^b'a_ow R^ow_fi_alpha b'
  double* RR1 = new double[nvir_occ];
  double* RR2 = new double[nvir_occ];
  fill_n(RR1, nvir_occ, 0.0);
  fill_n(RR2, nvir_occ, 0.0);

  // activate integrals
  Ref<DistArray4> abp_ii_ints;
  Ref<DistArray4> bpa_ii_ints;
  Ref<DistArray4> fibp_ii_ints;

  activate_ints(vir->id(), cabs->id(), occ_act->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, abp_ii_ints);
  activate_ints(cabs->id(), vir->id(), occ_act->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, bpa_ii_ints);
  activate_ints(focc->id(), cabs->id(), occ_act->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, fibp_ii_ints);

  // R^ab'_ow R^ow_fi b'
  compute_ai_sum_3idx(ax_ix, f12_idx, f12_idx,
                      abp_ii_ints, fibp_ii_ints,
                      RR1);

  // R^b'a_ow R^ow_fi b'
  compute_ai_sum_3idx(xa_ix, f12_idx, f12_idx,
                      bpa_ii_ints, fibp_ii_ints,
                      RR2);

  abp_ii_ints->deactivate();
  bpa_ii_ints->deactivate();
  fibp_ii_ints->deactivate();

  if (debug >= DefaultPrintThresholds::allN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    print_intermediate(spinletter, "R^ab'_ow R^ow_fi b'", RR1, nvir, nocc);
    print_intermediate(spinletter, "R^b'a_ow R^ow_fi b'", RR2, nvir, nocc);
  }

  // calculate AlphaBeta part:
  //             alpha case                   beta case
  // RR1_ab: R^a1b'2_o1w2 R^o1w2_fi1 b'2 or R^b'1a2_o1w2 R^o1w2_b'1 fi2
  // RR2_ab: R^b'2a1_o1w2 R^o1w2_fi1 b'2 or R^a2b'1_o1w2 R^o1w2_b'1 fi2
  // RR3_ab: R^a1b'2_o1w2 R^o1w2_b'2 fi1 or R^b'1a2_o1w2 R^o1w2_fi2 b'1
  // RR4_ab: R^b'2a1_o1w2 R^o1w2_b'2 fi1 or R^a2b'1_o1w2 R^o1w2_fi2 b'1

  double* RR1_ab = NULL;
  double* RR2_ab = NULL;
  double* RR3_ab = NULL;
  double* RR4_ab = NULL;
  if (nspincases2 == 3) {

     const Ref<OrbitalSpace> vir1 = r12eval->vir(spin1);
     const Ref<OrbitalSpace> cabs1 = r12world->cabs_space(spin1);
     const Ref<OrbitalSpace> vir2 = r12eval->vir(spin2);
     const Ref<OrbitalSpace> cabs2 = r12world->cabs_space(spin2);

     RR1_ab = new double[nvir_occ];
     RR2_ab = new double[nvir_occ];
     RR3_ab = new double[nvir_occ];
     RR4_ab = new double[nvir_occ];
     fill_n(RR1_ab, nvir_occ, 0.0);
     fill_n(RR2_ab, nvir_occ, 0.0);
     fill_n(RR3_ab, nvir_occ, 0.0);
     fill_n(RR4_ab, nvir_occ, 0.0);

     if (spin == Alpha) {
       Ref<OrbitalSpace> focc1;
       if (orbs_focc == "vir") {
         focc1 = r12eval->F_m_a(spin1);
       } else if (orbs_focc == "cabs") {
           focc1 = r12eval->F_m_A(spin1);
       }

       Ref<DistArray4> abpii_1212_ints;   // R^a1b'2_o1w2
       Ref<DistArray4> bpaii_2112_ints;   // R^b'2a1_o1w2
       Ref<DistArray4> fibpii_1212_ints;  // R^o1w2_fi1 b'2
       Ref<DistArray4> bpfiii_2112_ints;  // R^o1w2_b'2 fi1

       activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, abpii_1212_ints);
       activate_ints(cabs2->id(), vir1->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, bpaii_2112_ints);
       activate_ints(focc1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, fibpii_1212_ints);
       activate_ints(cabs2->id(), focc1->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, bpfiii_2112_ints);


       // R^a1b'2_o1w2 R^o1w2_fi1 b'2
       compute_ai_sum_3idx(ax_ix, f12_idx, f12_idx,
                           abpii_1212_ints, fibpii_1212_ints,
                           RR1_ab);
       // R^b'2a1_o1w2 R^o1w2_fi1 b'2
       compute_ai_sum_3idx(xa_ix, f12_idx, f12_idx,
                           bpaii_2112_ints, fibpii_1212_ints,
                           RR2_ab);
       // R^a1b'2_o1w2 R^o1w2_b'2 fi1
       compute_ai_sum_3idx(ax_xi, f12_idx, f12_idx,
                           abpii_1212_ints, bpfiii_2112_ints,
                           RR3_ab);
       // R^b'2a1_o1w2 R^o1w2_b'2 fi1
       compute_ai_sum_3idx(xa_xi, f12_idx, f12_idx,
                           bpaii_2112_ints, bpfiii_2112_ints,
                           RR4_ab);

       abpii_1212_ints->deactivate();
       bpaii_2112_ints->deactivate();
       fibpii_1212_ints->deactivate();
       bpfiii_2112_ints->deactivate();

       if (debug >= DefaultPrintThresholds::allN2) {
         print_intermediate("Alpha", "R^a1b'2_o1w2 R^o1w2_fi1 b'2", RR1_ab, nvir, nocc);
         print_intermediate("Alpha", "R^b'2a1_o1w2 R^o1w2_fi1 b'2", RR2_ab, nvir, nocc);
         print_intermediate("Alpha", "R^a1b'2_o1w2 R^o1w2_b'2 fi1", RR3_ab, nvir, nocc);
         print_intermediate("Alpha", "R^b'2a1_o1w2 R^o1w2_b'2 fi1", RR4_ab, nvir, nocc);
       }
     } else {
         Ref<OrbitalSpace> focc2;
         if (orbs_focc == "vir") {
           focc2 = r12eval->F_m_a(spin2);
         } else if (orbs_focc == "cabs") {
             focc2 = r12eval->F_m_A(spin2);
         }

         Ref<DistArray4> bpaii_1212_ints;   // R^b'1a2_o1w2
         Ref<DistArray4> abpii_2112_ints;   // R^a2b'1_o1w2
         Ref<DistArray4> bpfiii_1212_ints;  // R^o1w2_b'1 fi2
         Ref<DistArray4> fibpii_2112_ints;  // R^o1w2_fi2 b'1

         activate_ints(cabs1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, bpaii_1212_ints);
         activate_ints(vir2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, abpii_2112_ints);
         activate_ints(cabs1->id(), focc2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, bpfiii_1212_ints);
         activate_ints(focc2->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, fibpii_2112_ints);


         // R^b'1a2_o1w2 R^o1w2_b'1 fi2
         compute_ai_sum_3idx(xa_xi, f12_idx, f12_idx,
                             bpaii_1212_ints, bpfiii_1212_ints,
                             RR1_ab);
         // R^a2b'1_o1w2 R^o1w2_b'1 fi2
         compute_ai_sum_3idx(ax_xi, f12_idx, f12_idx,
                             abpii_2112_ints, bpfiii_1212_ints,
                             RR2_ab);
         // R^b'1a2_o1w2 R^o1w2_fi2 b'1
         compute_ai_sum_3idx(xa_ix, f12_idx, f12_idx,
                             bpaii_1212_ints, fibpii_2112_ints,
                             RR3_ab);
         // R^a2b'1_o1w2 R^o1w2_fi2 b'1
         compute_ai_sum_3idx(ax_ix, f12_idx, f12_idx,
                             abpii_2112_ints, fibpii_2112_ints,
                             RR4_ab);

         bpaii_1212_ints->deactivate();
         abpii_2112_ints->deactivate();
         bpfiii_1212_ints->deactivate();
         fibpii_2112_ints->deactivate();

         if (debug >= DefaultPrintThresholds::allN2) {
           print_intermediate("Beta", "R^b'1a2_o1w2 R^o1w2_b'1 fi2", RR1_ab, nvir, nocc);
           print_intermediate("Beta", "R^a2b'1_o1w2 R^o1w2_b'1 fi2", RR2_ab, nvir, nocc);
           print_intermediate("Beta", "R^b'1a2_o1w2 R^o1w2_fi2 b'1", RR3_ab, nvir, nocc);
           print_intermediate("Beta", "R^a2b'1_o1w2 R^o1w2_fi2 b'1", RR4_ab, nvir, nocc);
         }
     }
  } else if (nspincases2 == 2) {
       RR1_ab = RR1;
       RR2_ab = RR2;
       RR3_ab = RR2;
       RR4_ab = RR1;
  }
  // end of AlphaBeta part

  // calculate RR^a_i
  const double* iter_RR1 = RR1;
  const double* iter_RR2 = RR2;
  const double* iter_RR1_ab = RR1_ab;
  const double* iter_RR2_ab = RR2_ab;
  const double* iter_RR3_ab = RR3_ab;
  const double* iter_RR4_ab = RR4_ab;

  double* iter_Xai = Xai_RFR_ow;
   for (int a = 0;  a < nvir; ++a) {
     for (int i = 0; i < nocc; ++i) {

       // AlphaAlpha/BetaBeta part
       double x_12 = C_1 * C_1 * (*iter_RR1 - *iter_RR2);
       ++iter_RR1;
       ++iter_RR2;

       // AlphaBeta part
       if (nocc1 != 0 && nocc2 != 0) {
         x_12 += pow(0.5 * (C_0 + C_1), 2) * (*iter_RR1_ab)
               + 0.25 * (C_0 * C_0 - C_1 * C_1) * (*iter_RR2_ab + *iter_RR3_ab)
               + pow(0.5 * (C_0 - C_1), 2) * (*iter_RR4_ab);
         ++iter_RR1_ab;
         ++iter_RR2_ab;
         ++iter_RR3_ab;
         ++iter_RR4_ab;
       } // end of AlphaBeta part

       *iter_Xai = x_12;
       ++iter_Xai;
     }
   } // end of calculating Xai

   delete[] RR1;
   delete[] RR2;
   if (nspincases2 == 3) {
     delete[] RR1_ab;
     delete[] RR2_ab;
     delete[] RR3_ab;
     delete[] RR4_ab;
   }

}
// end of compute_RFR_ow

// compute 1/2 \bar{\tilde{R}}_{ic' }^{ow}\bar{g}_{ow}^{ac'}
void compute_RGow2(SpinCase1 spin, const double C_0, const double C_1,
                  const Ref<R12IntEval>& r12eval, int debug,
                  double* const Xai_RGow) {

  Ref<R12WavefunctionWorld> r12world = r12eval->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);

  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const TwoBodyOper::type eri_type = r12world->r12tech()->corrfactor()->tbint_type_eri();
  const unsigned int f12_idx = descr_f12->intset(f12_type);
  const unsigned int eri_idx = descr_f12->intset(eri_type);

  const int nspincases2 = (r12eval->spin_polarized() ? 3 : 2);

  const SpinCase1 spin1 = Alpha;
  const SpinCase1 spin2 = Beta;
  const Ref<OrbitalSpace> occ1_act = r12eval->occ_act(spin1);
  const Ref<OrbitalSpace> occ2_act = r12eval->occ_act(spin2);
  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();

  const Ref<OrbitalSpace> occ_act = (spin == Alpha? occ1_act : occ2_act);
  const int nocc_act = (spin == Alpha? nocc1_act : nocc2_act);

  const Ref<OrbitalSpace> vir = r12eval->vir(spin);
  const Ref<OrbitalSpace> cabs = r12world->cabs_space(spin);
  const int nvir = vir->rank();

  const int nvir_occ_act = nvir * nocc_act;

  // calculate AlphaAlpha/BetaBeta part:
  // g^ac'_ow R^ow_ic' & g^ac'_ow R^ow_c'i
  double* Rg1 = new double[nvir_occ_act];
  double* Rg2 = new double[nvir_occ_act];
  fill_n(Rg1, nvir_occ_act, 0.0);
  fill_n(Rg2, nvir_occ_act, 0.0);

  // activate integrals
  Ref<DistArray4> acp_ii_ints;
  Ref<DistArray4> icp_ii_ints;
  Ref<DistArray4> cpi_ii_ints;

  activate_ints(vir->id(), cabs->id(), occ_act->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, acp_ii_ints);
  activate_ints(occ_act->id(), cabs->id(), occ_act->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, icp_ii_ints);
  activate_ints(cabs->id(), occ_act->id(), occ_act->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, cpi_ii_ints);


  // g^ac'_ow R^ow_ic'
  compute_ai_sum_3idx(ax_ix, eri_idx, f12_idx,
                      acp_ii_ints, icp_ii_ints,
                      Rg1);

  // g^ac'_ow R^ow_c'i
  compute_ai_sum_3idx(ax_xi, eri_idx, f12_idx,
                      acp_ii_ints, cpi_ii_ints,
                      Rg2);

  acp_ii_ints->deactivate();
  icp_ii_ints->deactivate();
  cpi_ii_ints->deactivate();

  if (debug >= DefaultPrintThresholds::allN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    print_intermediate(spinletter, "R^ow_ic' g^ac'_ow", Rg1, nvir, nocc_act);
    print_intermediate(spinletter, "R^ow_c'i g^ac'_ow", Rg2, nvir, nocc_act);
  }

  // calculate AlphaBeta part:
  //             alpha case                   beta case
  // Rg1_ab: g^a1c'2_o1w2 R^o1w2_i1c'2  or g^c'1a2_o1w2 R^o1w2_c'1i2
  // Rg2_ab: g^a1c'2_o1w2 R^o1w2_c'2i1  or g^c'1a2_o1w2 R^o1w2_i2c'1

  double* Rg1_ab = NULL;
  double* Rg2_ab = NULL;
  if (nspincases2 == 3) {

     const Ref<OrbitalSpace> vir1 = r12eval->vir(spin1);
     const Ref<OrbitalSpace> cabs1 = r12world->cabs_space(spin1);
     const Ref<OrbitalSpace> vir2 = r12eval->vir(spin2);
     const Ref<OrbitalSpace> cabs2 = r12world->cabs_space(spin2);

     Rg1_ab = new double[nvir_occ_act];
     Rg2_ab = new double[nvir_occ_act];
     fill_n(Rg1_ab, nvir_occ_act, 0.0);
     fill_n(Rg2_ab, nvir_occ_act, 0.0);

     if (spin == Alpha) {

       Ref<DistArray4> icpii_1212_ints;   // R^o1w2_i1c'2
       Ref<DistArray4> cpiii_2112_ints;   // R^o1w2_c'2i1
       Ref<DistArray4> acpii_1212_ints;   // g^a1c'2_o1w2

       activate_ints(occ1_act->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, icpii_1212_ints);
       activate_ints(cabs2->id(), occ1_act->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, cpiii_2112_ints);
       activate_ints(vir1->id(), cabs2->id(), occ1_act->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, acpii_1212_ints);



       // Rg1_ab: g^a1c'2_o1w2 R^o1w2_i1c'2
       compute_ai_sum_3idx(ax_ix, eri_idx, f12_idx,
                           acpii_1212_ints, icpii_1212_ints,
                           Rg1_ab);
       // Rg2_ab: g^a1c'2_o1w2 R^o1w2_c'2i1
       compute_ai_sum_3idx(ax_xi, eri_idx, f12_idx,
                           acpii_1212_ints, cpiii_2112_ints,
                           Rg2_ab);

       icpii_1212_ints->deactivate();
       cpiii_2112_ints->deactivate();
       acpii_1212_ints->deactivate();


       if (debug >= DefaultPrintThresholds::allN2) {
         print_intermediate("Alpha", "R^o1w2_i1c'2 g^a1c'2_o1w2", Rg1_ab, nvir, nocc_act);
         print_intermediate("Alpha", "R^o1w2_c'2i1 g^a1c'2_o1w2", Rg2_ab, nvir, nocc_act);
       }
     } else {
         // g^c'1a2_o1w2 R^o1w2_c'1i2 or g^c'1a2_o1w2 R^o1w2_i2c'1
         Ref<DistArray4> cpiii_1212_ints;   // R^o1w2_c'1i2
         Ref<DistArray4> icpii_2112_ints;   // R^o1w2_i2c'1
         Ref<DistArray4> cpaii_1212_ints;   // g^c'1a2_o1w2

         activate_ints(cabs1->id(), occ2_act->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, cpiii_1212_ints);
         activate_ints(occ2_act->id(), cabs1->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, icpii_2112_ints);
         activate_ints(cabs1->id(), vir2->id(), occ1_act->id(), occ2_act->id(),
                       descr_f12_key, moints4_rtime, cpaii_1212_ints);



         // g^c'1a2_o1w2 R^o1w2_c'1i2
         compute_ai_sum_3idx(xa_xi, eri_idx, f12_idx,
                             cpaii_1212_ints, cpiii_1212_ints,
                             Rg1_ab);
         // g^c'1a2_o1w2 R^o1w2_i2c'1
         compute_ai_sum_3idx(xa_ix, eri_idx, f12_idx,
                             cpaii_1212_ints, icpii_2112_ints,
                             Rg2_ab);

         cpiii_1212_ints->deactivate();
         icpii_2112_ints->deactivate();
         cpaii_1212_ints->deactivate();

         if (debug >= DefaultPrintThresholds::allN2) {
           print_intermediate("Beta", "R^o1w2_c'1i2 g^c'1a2_o1w2", Rg1_ab, nvir, nocc_act);
           print_intermediate("Beta", "R^o1w2_i2c'1 g^c'1a2_o1w2", Rg2_ab, nvir, nocc_act);
         }
     }
  } else if (nspincases2 == 2) {
       Rg1_ab = Rg1;
       Rg2_ab = Rg2;
  }
  // end of AlphaBeta part

  // calculate Rg^a_i
  const double* iter_Rg1 = Rg1;
  const double* iter_Rg2 = Rg2;
  const double* iter_Rg1_ab = Rg1_ab;
  const double* iter_Rg2_ab = Rg2_ab;

  double* iter_Xai = Xai_RGow;
   for (int a = 0;  a < nvir; ++a) {
     for (int i = 0; i < nocc_act; ++i) {

       // AlphaAlpha/BetaBeta part
       double x_12 = C_1 * (*iter_Rg1 - *iter_Rg2);
       ++iter_Rg1;
       ++iter_Rg2;

       // AlphaBeta part
       if (nocc1_act != 0 && nocc2_act != 0) {
         x_12 += 0.5 * (C_0 + C_1) * (*iter_Rg1_ab)
                + 0.5 * (C_0 - C_1) * (*iter_Rg2_ab);
         ++iter_Rg1_ab;
         ++iter_Rg2_ab;
       } // end of AlphaBeta part

       *iter_Xai = x_12;
       ++iter_Xai;
     }
   } // end of calculating Xai

   delete[] Rg1;
   delete[] Rg2;
   if (nspincases2 == 3) {
     delete[] Rg1_ab;
     delete[] Rg2_ab;
   }

}
// end of compute_RGow2


// compute 1/2 \bar{\tilde{R}}^ao}_{\sigma\delta} \bar{g}^{\sigma\delta}_io}
void compute_RG_SigmaDelta2(SpinCase1 spin, const double C_0, const double C_1,
                  const Ref<R12IntEval>& r12eval, int debug,
                  double* const Xai_RG_SigmaDelta) {

  Ref<R12WavefunctionWorld> r12world = r12eval->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();

  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);

  const TwoBodyOper::type f12eri_type = r12world->r12tech()->corrfactor()->tbint_type_f12eri();
  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
  const TwoBodyOper::type eri_type = r12world->r12tech()->corrfactor()->tbint_type_eri();
  const unsigned int f12eri_idx = descr_f12->intset(f12eri_type);
  const unsigned int f12_idx = descr_f12->intset(f12_type);
  const unsigned int eri_idx = descr_f12->intset(eri_type);

  const int nspincases2 = (r12eval->spin_polarized() ? 3 : 2);

  const SpinCase1 spin1 = Alpha;
  const SpinCase1 spin2 = Beta;
  const Ref<OrbitalSpace> occ1_act = r12eval->occ_act(spin1);
  const Ref<OrbitalSpace> occ2_act = r12eval->occ_act(spin2);
  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();

  const Ref<OrbitalSpace> occ1 = r12eval->occ(spin1);
  const Ref<OrbitalSpace> occ2 = r12eval->occ(spin2);
  const int nocc1 = occ1->rank();
  const int nocc2 = occ2->rank();
  const int nocc = (spin == Alpha? nocc1 : nocc2);

  const Ref<OrbitalSpace> occ_act = (spin == Alpha? occ1_act : occ2_act);
  const Ref<OrbitalSpace> occ = (spin == Alpha? occ1 : occ2);

  const Ref<OrbitalSpace> vir = r12eval->vir(spin);
  const int nvir = vir->rank();
  const int nvir_occ = nvir * nocc;

  // calculate AlphaAlpha/BetaBeta part:
  //   R^ao_SigmaDelta g^SigmaDelta_io
  // & R^oa_SigmaDelta g^SigmaDelta_io

  double* Rg1 = new double[nvir_occ];
  double* Rg2 = new double[nvir_occ];
  fill_n(Rg1, nvir_occ, 0.0);
  fill_n(Rg2, nvir_occ, 0.0);

  // R^ao_SigmaDelta g^SigmaDelta_io
  //  = (f12/r12)^ao_io - f12_pq^ao g_io^pq
  //  - f12_ma'^ao g_io^ma' - f12_a'm_ao g_io^a'm

  // activate integrals
  Ref<DistArray4> aiii_ints;
  Ref<DistArray4> iaii_ints;
  vector<Ref<DistArray4> > v_ai_ints;
  vector<Ref<DistArray4> > v_ia_ints;
  vector<Ref<DistArray4> > v_ii_ints;

  activate_ints(vir->id(), occ_act->id(), occ->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, aiii_ints);
  activate_ints(occ_act->id(), vir->id(), occ->id(), occ_act->id(),
                descr_f12_key, moints4_rtime, iaii_ints);
  const SpinCase2 spincase2 = (spin == Alpha? AlphaAlpha: BetaBeta);
  activate_ints_VorX(vir, occ_act, spincase2,
                     r12eval, descr_f12_key, v_ai_ints);
  activate_ints_VorX(occ_act, vir, spincase2,
                     r12eval, descr_f12_key, v_ia_ints);
  activate_ints_VorX(occ, occ_act, spincase2,
                     r12eval, descr_f12_key, v_ii_ints);

  // R^ao_SigmaDelta g^SigmaDelta_io
  compute_ai_sum_abx(ax_ix, f12eri_idx, f12_idx, eri_idx,
                     aiii_ints, v_ai_ints, v_ii_ints, Rg1);
  // R^oa_SigmaDelta g^SigmaDelta_io
  compute_ai_sum_abx(xa_ix, f12eri_idx, f12_idx, eri_idx,
                     iaii_ints, v_ia_ints, v_ii_ints, Rg2);

  aiii_ints->deactivate();
  iaii_ints->deactivate();
  for (int i = 0; i != v_ai_ints.size(); ++i) {
    v_ai_ints[i]->deactivate();
    v_ia_ints[i]->deactivate();
    v_ii_ints[i]->deactivate();
  }


  if (debug >= DefaultPrintThresholds::allN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    print_intermediate(spinletter, "R^ao_SigmaDelta g^SigmaDelta_io", Rg1, nvir, nocc);
    print_intermediate(spinletter, "R^oa_SigmaDelta g^SigmaDelta_io", Rg2, nvir, nocc);
  }

  // calculate AlphaBeta part:
  // Rg1_ab alpha: R^a1o2_Sigma1Delta2 g^Sigma1Delta1_i1o2
  //        beta : R^o1a2_Sigma1Delta2 g^Sigma1Delta1_o1i2
  //
  // Rg2_ab alpha: R^o1a2_Sigma1Delta2 g^Sigma1Delta1_i1o2
  //        beta : R^a2o1_Sigma1Delta2 g^Sigma1Delta1_o1i2

  double* Rg1_ab = NULL;
  double* Rg2_ab = NULL;
  if (nspincases2 == 3) {

     const Ref<OrbitalSpace> vir1 = r12eval->vir(spin1);
     const Ref<OrbitalSpace> vir2 = r12eval->vir(spin2);

     Rg1_ab = new double[nvir_occ];
     Rg2_ab = new double[nvir_occ];
     fill_n(Rg1_ab, nvir_occ, 0.0);
     fill_n(Rg2_ab, nvir_occ, 0.0);

     if (spin == Alpha) {

       Ref<DistArray4> a1i2i1i2_ints;
       Ref<DistArray4> i2a1i1i2_ints;
       vector<Ref<DistArray4> > v_a1i2_ints;
       vector<Ref<DistArray4> > v_i2a1_ints;
       vector<Ref<DistArray4> > v_i1i2_ints;

       activate_ints(vir1->id(), occ2_act->id(), occ1->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, a1i2i1i2_ints);
       activate_ints(occ2_act->id(), vir1->id(), occ1->id(), occ2_act->id(),
                     descr_f12_key, moints4_rtime, i2a1i1i2_ints);
       activate_ints_VorX(vir1, occ2_act, AlphaBeta,
                          r12eval, descr_f12_key, v_a1i2_ints);
       activate_ints_VorX(occ2_act, vir1, AlphaBeta,
                          r12eval, descr_f12_key, v_i2a1_ints);
       activate_ints_VorX(occ1, occ2_act, AlphaBeta,
                          r12eval, descr_f12_key, v_i1i2_ints);

       // R^a1o2_Sigma1Delta2 g^Sigma1Delta1_i1o2
       compute_ai_sum_abx(ax_ix, f12eri_idx, f12_idx, eri_idx,
                          a1i2i1i2_ints, v_a1i2_ints, v_i1i2_ints, Rg1_ab);
       // R^o1a2_Sigma1Delta2 g^Sigma1Delta1_i1o2
       compute_ai_sum_abx(xa_ix, f12eri_idx,f12_idx, eri_idx,
                          i2a1i1i2_ints, v_i2a1_ints, v_i1i2_ints, Rg2_ab);

       a1i2i1i2_ints->deactivate();
       i2a1i1i2_ints->deactivate();
       for (int i = 0; i != v_a1i2_ints.size(); ++i) {
         v_a1i2_ints[i]->deactivate();
         v_i2a1_ints[i]->deactivate();
         v_i1i2_ints[i]->deactivate();
       }

       if (debug >= DefaultPrintThresholds::allN2) {
         print_intermediate("Alpha", "R^a1o2_Sigma1Delta2 g^Sigma1Delta1_i1o2", Rg1_ab, nvir, nocc);
         print_intermediate("Alpha", "R^o1a2_Sigma1Delta2 g^Sigma1Delta1_i1o2", Rg2_ab, nvir, nocc);
       }

     } else {
         Ref<DistArray4> i1a2i1i2_ints;
         Ref<DistArray4> a2i1i1i2_ints;
         vector<Ref<DistArray4> > v_i1a2_ints;
         vector<Ref<DistArray4> > v_a2i1_ints;
         vector<Ref<DistArray4> > v_i1i2_ints;

         activate_ints(occ1_act->id(), vir2->id(), occ1_act->id(), occ2->id(),
                       descr_f12_key, moints4_rtime, i1a2i1i2_ints);
         activate_ints(vir2->id(), occ1_act->id(), occ1_act->id(), occ2->id(),
                       descr_f12_key, moints4_rtime, a2i1i1i2_ints);
         activate_ints_VorX(occ1_act, vir2, AlphaBeta,
                            r12eval, descr_f12_key, v_i1a2_ints);
         activate_ints_VorX(vir2, occ1_act, AlphaBeta,
                            r12eval, descr_f12_key, v_a2i1_ints);
         activate_ints_VorX(occ1_act, occ2, AlphaBeta,
                            r12eval, descr_f12_key, v_i1i2_ints);


         // R^o1a2_Sigma1Delta2 g^Sigma1Delta1_o1i2
         compute_ai_sum_abx(xa_xi, f12eri_idx, f12_idx, eri_idx,
                            i1a2i1i2_ints, v_i1a2_ints, v_i1i2_ints, Rg1_ab);
         // R^a2o1_Sigma1Delta2 g^Sigma1Delta1_o1i2
         compute_ai_sum_abx(ax_xi, f12eri_idx, f12_idx, eri_idx,
                            a2i1i1i2_ints, v_a2i1_ints, v_i1i2_ints, Rg2_ab);

         i1a2i1i2_ints->deactivate();
         a2i1i1i2_ints->deactivate();
         for (int i = 0; i != v_i1a2_ints.size(); ++i) {
           v_i1a2_ints[i]->deactivate();
           v_a2i1_ints[i]->deactivate();
           v_i1i2_ints[i]->deactivate();
         }

         if (debug >= DefaultPrintThresholds::allN2) {
           print_intermediate("Beta", "R^o1a2_Sigma1Delta2 g^Sigma1Delta1_o1i2", Rg1_ab, nvir, nocc);
           print_intermediate("Beta", "R^a2o1_Sigma1Delta2 g^Sigma1Delta1_o1i2", Rg2_ab, nvir, nocc);
         }

     }
  } else if (nspincases2 == 2) {
       Rg1_ab = Rg1;
       Rg2_ab = Rg2;
  }
  // end of AlphaBeta part

  // calculate Rg^a_i
  const double* iter_Rg1 = Rg1;
  const double* iter_Rg2 = Rg2;
  const double* iter_Rg1_ab = Rg1_ab;
  const double* iter_Rg2_ab = Rg2_ab;

  double* iter_Xai = Xai_RG_SigmaDelta;
   for (int a = 0;  a < nvir; ++a) {
     for (int i = 0; i < nocc; ++i) {

       // AlphaAlpha/BetaBeta part
       double x_12 = C_1 * (*iter_Rg1 - *iter_Rg2);
       ++iter_Rg1;
       ++iter_Rg2;

       // AlphaBeta part
       if (nocc1_act != 0 && nocc2_act != 0) {
         x_12 += 0.5 * (C_0 + C_1) * (*iter_Rg1_ab)
                + 0.5 * (C_0 - C_1) * (*iter_Rg2_ab);
         ++iter_Rg1_ab;
         ++iter_Rg2_ab;
       } // end of AlphaBeta part

       *iter_Xai = x_12;
       ++iter_Xai;
     }
   } // end of calculating Xai

   delete[] Rg1;
   delete[] Rg2;
   if (nspincases2 == 3) {
     delete[] Rg1_ab;
     delete[] Rg2_ab;
   }

}
// end of compute_RG_SigmaDelta2

//// compute 1/2 \bar{\tilde{R}}^{ow}_{\alpha\beta}
////          *  F^a_o
////          *  \bar{\tilde{R}}^{\alpha\beta}_iw
//void compute_RFR_AlphaBeta(SpinCase1 spin, const double C_0, const double C_1,
//                           const Ref<R12IntEval>& r12eval, int debug,
//                           double* const RFR_AlphaBeta) {
//
//  Ref<R12WavefunctionWorld> r12world = r12eval->r12world();
//  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
//
//  Ref<TwoBodyIntDescr> descr_f12f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0,0);
//  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
//  const std::string descr_f12f12_key = moints4_rtime->descr_key(descr_f12f12);
//  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);
//
//  const TwoBodyOper::type f12f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12f12();
//  const TwoBodyOper::type f12_type = r12world->r12tech()->corrfactor()->tbint_type_f12();
//
//  const unsigned int f12f12_idx = descr_f12f12->intset(f12f12_type);
//  const unsigned int f12_idx = descr_f12->intset(f12_type);
//
//  const int nspincases2 = (r12eval->spin_polarized() ? 3 : 2);
//
//  const SpinCase1 spin1 = Alpha;
//  const SpinCase1 spin2 = Beta;
//  const Ref<OrbitalSpace> occ1 = r12eval->occ(spin1);
//  const Ref<OrbitalSpace> occ2 = r12eval->occ(spin2);
//  const int nocc1 = occ1->rank();
//  const int nocc2 = occ2->rank();
//
//  const Ref<OrbitalSpace> occ = (spin == Alpha? occ1 : occ2);
//  const int nocc = (spin == Alpha? nocc1 : nocc2);
//
//  const Ref<OrbitalSpace> vir = r12eval->vir(spin);
//  const Ref<OrbitalSpace> fvir = r12eval->F_m_a(spin);
//  ExEnv::out0() << endl << indent << "Rank of fvir (F_m_a): " << fvir->rank() << endl;
//  const int nvir = vir->rank();
//  const int nvir_occ = nvir * nocc;
//
//  // calculate AlphaAlpha/BetaBeta part:
//  //   R^fa_o w_AlphaBeta R^AlphaBeta_iw
//  // & R^w fa_o _AlphaBeta R^AlphaBeta_iw
//
//  double* RR1 = new double[nvir_occ];
//  double* RR2 = new double[nvir_occ];
//  fill_n(RR1, nvir_occ, 0.0);
//  fill_n(RR2, nvir_occ, 0.0);
//
//  // R^fa_o w_AlphaBeta R^AlphaBeta_iw
//  //  = (f12f12)^fa_o w_iw - f12^fa_o w_pq f12^pq_iw
//  //  - f12^fa_o w_ma' f12^ma'_iw - f12^fa_o w_a'm f12^a'm_iw
//
//  // activate integrals
//  Ref<DistArray4> faiii_ints;
//  vector<Ref<DistArray4> > v_fai_ints;
//  vector<Ref<DistArray4> > v_ii_ints;
//
//  activate_ints(fvir->id(), occ->id(), occ->id(), occ->id(),
//                descr_f12f12_key, moints4_rtime, faiii_ints);
//  const SpinCase2 spincase2 = (spin == Alpha? AlphaAlpha: BetaBeta);
//  activate_ints_VorX(fvir, occ, spincase2,
//                     r12eval, descr_f12_key, v_fai_ints);
//  activate_ints_VorX(occ, occ, spincase2,
//                     r12eval, descr_f12_key, v_ii_ints);
//
//  // R^fa_o w_AlphaBeta R^AlphaBeta_iw
//  compute_ai_sum_abx(ax_ix, f12f12_idx, f12_idx, f12_idx,
//                     faiii_ints, v_fai_ints, v_ii_ints, RR1);
//  // R^w fa_o _AlphaBeta R^AlphaBeta_iw
//  compute_ai_sum_abx(ax_xi, f12f12_idx, f12_idx, f12_idx,
//                     faiii_ints, v_fai_ints, v_ii_ints, RR2);
//
//  faiii_ints->deactivate();
//  for (int i = 0; i != v_fai_ints.size(); ++i) {
//    v_fai_ints[i]->deactivate();
//    v_ii_ints[i]->deactivate();
//  }
//
//
////  if (debug >= DefaultPrintThresholds::mostN2) {
//    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
//    print_intermediate(spinletter, "R^fa_o w_AlphaBeta R^AlphaBeta_iw", RR1, nvir, nocc);
//    print_intermediate(spinletter, "R^w fa_o _AlphaBeta R^AlphaBeta_iw", RR2, nvir, nocc);
////  }
//
//    // test for Xai_RG_AlphaBeta
//    {
//#if 0
//      RefSCDimension rowdim = new SCDimension(nvir);
//      RefSCDimension coldim = new SCDimension(nocc);
//      Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
//      RefSCMatrix Xai = localkit->matrix(rowdim, coldim);
//      Xai.assign(0.0);
//
//      const Ref<OrbitalSpace>& orbs = r12eval->orbs(spin);
//      // Vpq_ij
//      //const SpinCase2 spincase2 = static_cast<SpinCase2>(spin+1);
//      std::vector< Ref<DistArray4> > Vpq_vec = r12eval->V_distarray4(AlphaAlpha, orbs, orbs);
//      MPQC_ASSERT(Vpq_vec.size() == 1);
//      Ref<DistArray4> Vpq_ij = Vpq_vec[0];
//      Ref<DistArray4> Vaj_ij;
//      sc::map(Vpq_ij, occ, occ, orbs, orbs, Vaj_ij, occ, occ, vir, occ);
//
//      Vaj_ij->activate();
//
//      for (int a = 0;  a < nvir; ++a) {
//        for (int i = 0; i < nocc; ++i) {
//
//          double Va_i = 0.0;
//          for (int j = 0; j < nocc; ++j) {
//            //ExEnv::out0() << "a j i:" << a << " " << j << " " << i << endl;
//            const double* const Vij_blk = Vaj_ij->retrieve_pair_block(i, j, 0);
//            const int aj_ij = a * nocc + j;
//            Va_i += Vij_blk[aj_ij];
//
//            Vaj_ij->release_pair_block(i, j, 0);
//          }
//
//          Xai.set_element(a, i, Va_i);
//        }
//      }
//      Vaj_ij->deactivate();
//      Xai.print(prepend_spincase(spin,"Xai from Vai_ij:").c_str());
//
//      double* Rg12 = new double[nvir_occ];
//      fill_n(Rg12, nvir_occ, 0.0);
//
//      double* iter_Rg12 = Rg12;
//      double* iter_Rg1 = Rg1;
//      double* iter_Rg2 = Rg2;
//      for (int a = 0;  a < nvir; ++a) {
//        for (int i = 0; i < nocc; ++i) {
//          *iter_Rg12 = *iter_Rg1 - *iter_Rg2;
//
//          ++iter_Rg12;
//          ++iter_Rg1;
//          ++iter_Rg2;
//        }
//      }
//      print_intermediate(spinletter, "antisymmetrized G * R", Rg12, nvir, nocc);
//      delete[] Rg12;
//#endif
//    } // end of test
//
//  // calculate AlphaBeta part:
//  //                     Alpha
//  // Rg1_ab: R^fa1_o w2_Alpha1Beta2  R^Alpha1Beta2_i1w2
//  // Rg2_ab: R^w2 fa1_o _Alpha1Beta2 R^Alpha1Beta2_i1w2
//  // Rg3_ab: R^fa1_o w2_Alpha1Beta2  R^Alpha1Beta2_w2i1
//  // Rg4_ab: R^w2 fa1_o _Alpha1Beta2 R^Alpha1Beta2_w2i1
//  //                     Beta
//  // Rg1_ab: R^w1 fa2_o _Alpha1Beta2 R^Alpha1Beta2_w1i2
//  // Rg2_ab: R^fa2_o w1_Alpha1Beta2  R^Alpha1Beta2_w1i2
//  // Rg3_ab: R^w1 fa2_o _Alpha1Beta2 R^Alpha1Beta2_i2w1
//  // Rg4_ab: R^fa2_o w1_Alpha1Beta2  R^Alpha1Beta2_i2w1
//
//  double* RR1_ab = NULL;
//  double* RR2_ab = NULL;
//  double* RR3_ab = NULL;
//  double* RR4_ab = NULL;
//  if (nspincases2 == 3) {
//     const Ref<OrbitalSpace> vir1 = r12eval->vir(spin1);
//     const Ref<OrbitalSpace> vir2 = r12eval->vir(spin2);
//
//     vector<Ref<DistArray4> > v_i1i2_ints;
//     vector<Ref<DistArray4> > v_i2i1_ints;
//
//     activate_ints_VorX(occ1, occ2, AlphaBeta,
//                        r12eval, descr_f12_key, v_i1i2_ints);
//     activate_ints_VorX(occ2, occ1, AlphaBeta,
//                        r12eval, descr_f12_key, v_i2i1_ints);
//
//     RR1_ab = new double[nvir_occ];
//     RR2_ab = new double[nvir_occ];
//     RR3_ab = new double[nvir_occ];
//     RR4_ab = new double[nvir_occ];
//     fill_n(RR1_ab, nvir_occ, 0.0);
//     fill_n(RR2_ab, nvir_occ, 0.0);
//     fill_n(RR3_ab, nvir_occ, 0.0);
//     fill_n(RR4_ab, nvir_occ, 0.0);
//
//     if (spin == Alpha) {
//       const Ref<OrbitalSpace> fvir1 = r12eval->F_m_a(spin1);
//       Ref<DistArray4> fa1i2i1i2_ints;
//       Ref<DistArray4> i2fa1i1i2_ints;
//       Ref<DistArray4> fa1i2i2i1_ints;
//       Ref<DistArray4> i2fa1i2i1_ints;
//       vector<Ref<DistArray4> > v_fa1i2_ints;
//       vector<Ref<DistArray4> > v_i2fa1_ints;
//
//       activate_ints(fvir1->id(), occ2->id(), occ1->id(), occ2->id(),
//                     descr_f12_key, moints4_rtime, fa1i2i1i2_ints);
//       activate_ints(occ2->id(), fvir1->id(), occ1->id(), occ2->id(),
//                     descr_f12_key, moints4_rtime, i2fa1i1i2_ints);
//       activate_ints(fvir1->id(), occ2->id(), occ2->id(), occ1->id(),
//                     descr_f12_key, moints4_rtime, fa1i2i2i1_ints);
//       activate_ints(occ2->id(), fvir1->id(), occ2->id(), occ1->id(),
//                     descr_f12_key, moints4_rtime, i2fa1i2i1_ints);
//       activate_ints_VorX(fvir1, occ2, AlphaBeta,
//                          r12eval, descr_f12_key, v_fa1i2_ints);
//       activate_ints_VorX(occ2, fvir1, AlphaBeta,
//                          r12eval, descr_f12_key, v_i2fa1_ints);
//
//       // R^fa1_o w2_Alpha1Beta2  R^Alpha1Beta2_i1w2
//       compute_ai_sum_abx(ax_ix, f12f12_idx, f12_idx, f12_idx,
//                          fa1i2i1i2_ints, v_fa1i2_ints, v_i1i2_ints, RR1_ab);
//       // R^w2 fa1_o _Alpha1Beta2 R^Alpha1Beta2_i1w2
//       compute_ai_sum_abx(xa_ix, f12f12_idx, f12_idx, f12_idx,
//                          i2fa1i1i2_ints, v_i2fa1_ints, v_i1i2_ints, RR2_ab);
//       // R^fa1_o w2_Alpha1Beta2  R^Alpha1Beta2_w2i1
//       compute_ai_sum_abx(ax_xi, f12f12_idx, f12_idx, f12_idx,
//                          fa1i2i2i1_ints, v_fa1i2_ints, v_i2i1_ints, RR3_ab);
//       // R^w2 fa1_o _Alpha1Beta2 R^Alpha1Beta2_w2i1
//       compute_ai_sum_abx(xa_xi, f12f12_idx, f12_idx, f12_idx,
//                          i2fa1i2i1_ints, v_i2fa1_ints, v_i2i1_ints, RR4_ab);
//
//       fa1i2i1i2_ints->deactivate();
//       i2fa1i1i2_ints->deactivate();
//       fa1i2i2i1_ints->deactivate();
//       i2fa1i2i1_ints->deactivate();
//       for (int i = 0; i != v_fa1i2_ints.size(); ++i) {
//         v_fa1i2_ints[i]->deactivate();
//         v_i2fa1_ints[i]->deactivate();
//       }
//
//       //if (debug >= DefaultPrintThresholds::mostN2) {
//         print_intermediate("Alpha", "R^fa1_o w2_Alpha1Beta2  R^Alpha1Beta2_i1w2", RR1_ab, nvir, nocc);
//         print_intermediate("Alpha", "R^w2 fa1_o _Alpha1Beta2 R^Alpha1Beta2_i1w2", RR2_ab, nvir, nocc);
//         print_intermediate("Alpha", "R^fa1_o w2_Alpha1Beta2  R^Alpha1Beta2_w2i1", RR3_ab, nvir, nocc);
//         print_intermediate("Alpha", "R^w2 fa1_o _Alpha1Beta2 R^Alpha1Beta2_w2i1", RR4_ab, nvir, nocc);
//       //}
//     } else {
//         const Ref<OrbitalSpace> fvir2 = r12eval->F_m_a(spin2);
//         Ref<DistArray4> i1fa2i1i2_ints;
//         Ref<DistArray4> fa2i1i1i2_ints;
//         Ref<DistArray4> i1fa2i2i1_ints;
//         Ref<DistArray4> fa2i1i2i1_ints;
//         vector<Ref<DistArray4> > v_i1fa2_ints;
//         vector<Ref<DistArray4> > v_fa2i1_ints;
//
//         activate_ints(occ1->id(), fvir2->id(), occ1->id(), occ2->id(),
//                       descr_f12_key, moints4_rtime, i1fa2i1i2_ints);
//         activate_ints(fvir2->id(), occ1->id(), occ1->id(), occ2->id(),
//                       descr_f12_key, moints4_rtime, fa2i1i1i2_ints);
//         activate_ints(occ1->id(), fvir2->id(), occ2->id(), occ1->id(),
//                       descr_f12_key, moints4_rtime, i1fa2i2i1_ints);
//         activate_ints(fvir2->id(), occ1->id(), occ2->id(), occ1->id(),
//                       descr_f12_key, moints4_rtime, fa2i1i2i1_ints);
//         activate_ints_VorX(occ1, fvir2, AlphaBeta,
//                            r12eval, descr_f12_key, v_i1fa2_ints);
//         activate_ints_VorX(fvir2, occ1, AlphaBeta,
//                            r12eval, descr_f12_key, v_fa2i1_ints);
//
//         // R^w1 fa2_o _Alpha1Beta2 R^Alpha1Beta2_w1i2
//         compute_ai_sum_abx(xa_xi, f12f12_idx, f12_idx, f12_idx,
//                            i1fa2i1i2_ints, v_i1fa2_ints, v_i1i2_ints, RR1_ab);
//         // R^fa2_o w1_Alpha1Beta2  R^Alpha1Beta2_w1i2
//         compute_ai_sum_abx(ax_xi, f12f12_idx, f12_idx, f12_idx,
//                            fa2i1i1i2_ints, v_fa2i1_ints, v_i1i2_ints, RR2_ab);
//         // R^w1 fa2_o _Alpha1Beta2 R^Alpha1Beta2_i2w1
//         compute_ai_sum_abx(xa_ix, f12f12_idx, f12_idx, f12_idx,
//                            i1fa2i2i1_ints, v_i1fa2_ints, v_i2i1_ints, RR3_ab);
//         // R^fa2_o w1_Alpha1Beta2  R^Alpha1Beta2_i2w1
//         compute_ai_sum_abx(ax_ix, f12f12_idx, f12_idx, f12_idx,
//                            fa2i1i2i1_ints, v_fa2i1_ints, v_i2i1_ints, RR4_ab);
//
//         i1fa2i1i2_ints->deactivate();
//         fa2i1i1i2_ints->deactivate();
//         i1fa2i2i1_ints->deactivate();
//         fa2i1i2i1_ints->deactivate();
//         for (int i = 0; i != v_i1fa2_ints.size(); ++i) {
//           v_i1fa2_ints[i]->deactivate();
//           v_fa2i1_ints[i]->deactivate();
//         }
//
//         //if (debug >= DefaultPrintThresholds::mostN2) {
//           print_intermediate("Alpha", "R^w1 fa2_o _Alpha1Beta2 R^Alpha1Beta2_w1i2", RR1_ab, nvir, nocc);
//           print_intermediate("Alpha", "R^fa2_o w1_Alpha1Beta2  R^Alpha1Beta2_w1i2", RR2_ab, nvir, nocc);
//           print_intermediate("Alpha", "R^w1 fa2_o _Alpha1Beta2 R^Alpha1Beta2_i2w1", RR3_ab, nvir, nocc);
//           print_intermediate("Alpha", "R^fa2_o w1_Alpha1Beta2  R^Alpha1Beta2_i2w1", RR4_ab, nvir, nocc);
//         //}
//     }
//     for (int i = 0; i != v_i1i2_ints.size(); ++i) {
//       v_i1i2_ints[i]->deactivate();
//       v_i2i1_ints[i]->deactivate();
//     }
//  } else if (nspincases2 == 2) {
//       RR1_ab = RR1;
//       RR2_ab = RR2;
//       RR3_ab = RR2;
//       RR4_ab = RR1;
//  }
//  // end of AlphaBeta part
//
//  // calculate Rg^a_i
//  const double* iter_RR1 = RR1;
//  const double* iter_RR2 = RR2;
//  const double* iter_RR1_ab = RR1_ab;
//  const double* iter_RR2_ab = RR2_ab;
//  const double* iter_RR3_ab = RR3_ab;
//  const double* iter_RR4_ab = RR4_ab;
//
//  double* iter_RFR = RFR_AlphaBeta;
//   for (int a = 0;  a < nvir; ++a) {
//     for (int i = 0; i < nocc; ++i) {
//
//       // AlphaAlpha/BetaBeta part
//       double x_12 = C_1 * C_1 * ((*iter_RR1) - (*iter_RR2));
//       ++iter_RR1;
//       ++iter_RR2;
//
//       // AlphaBeta part
//       if (nocc1 != 0 && nocc2 != 0) {
//         x_12 += pow(0.5 * (C_0 + C_1), 2) * (*iter_RR1_ab)
//                   + 0.25 * (C_0 * C_0 - C_1 * C_1) * (*iter_RR2_ab + *iter_RR3_ab)
//                   + pow(0.5 * (C_0 - C_1), 2) * (*iter_RR4_ab);
//         ++iter_RR1_ab;
//         ++iter_RR2_ab;
//         ++iter_RR3_ab;
//         ++iter_RR4_ab;
//       } // end of AlphaBeta part
//
//       *iter_RFR = x_12;
//       ++iter_RFR;
//     }
//   } // end of calculating Xai
//
//   //print_intermediate(spinletter, "Xai_RFR_AlphaBeta", RFR_AlphaBeta, nvir, nocc);
//
//   delete[] RR1;
//   delete[] RR2;
//   if (nspincases2 == 3) {
//     delete[] RR1_ab;
//     delete[] RR2_ab;
//     delete[] RR3_ab;
//     delete[] RR4_ab;
//   }
//
//}
//// end of compute_RFR_AlphaBeta

RefSCMatrix sc::Onerdm_X_F12(SpinCase1 spin, const Ref<R12IntEval>& r12eval, int debug) {

  // geminal coefficients
  const double C_0 = 1.0/2.0;
  const double C_1 = 1.0/4.0;

  const Ref<OrbitalSpace> occ_act = r12eval->occ_act(spin);
  const Ref<OrbitalSpace> occ = r12eval->occ(spin);
  const Ref<OrbitalSpace> vir = r12eval->vir(spin);
  const int nocc_act = occ_act->rank();
  const int nocc = occ->rank();
  const int nvir = vir->rank();
  const int nvir_occ_act = nvir * nocc_act;
  const int nvir_occ = nvir * nocc;

  RefSCDimension rowdim = new SCDimension(nvir);
  RefSCDimension coldim = new SCDimension(nocc);
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix Xai = localkit->matrix(rowdim, coldim);
  Xai.assign(0.0);

  // compute contribution from
  // 1/2 \bar{\tilde{R}}^{ow}_{ac'} \bar{g}^{ic'}_{ow}
  double* RGow = new double[nvir_occ];
  fill_n(RGow, nvir_occ, 0.0);
  compute_RGow(spin, C_0, C_1, r12eval, debug, RGow);
  if (debug >= DefaultPrintThresholds::mostN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    print_intermediate(spinletter, "Xai_RGow", RGow, nvir, nocc);
  }
  // 1/2 \bar{\tilde{R}}^{ow}_{ic'} \bar{g}^{ac'}_{ow}
  double* RGow2 = new double[nvir_occ_act];
  fill_n(RGow2, nvir_occ_act, 0.0);
  compute_RGow2(spin, C_0, C_1, r12eval, debug, RGow2);
  if (debug >= DefaultPrintThresholds::mostN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    print_intermediate(spinletter, "Xai_RGow2", RGow2, nvir, nocc_act);
  }

  // compute contribution from
  // 1/2 \bar{\tilde{R}}^io}_{\sigma\delta} \bar{g}^{\sigma\delta}_{ao}
  double* RG_SigmaDelta = new double[nvir_occ_act];
  fill_n(RG_SigmaDelta, nvir_occ_act, 0.0);
  compute_RG_SigmaDelta(spin, C_0, C_1, r12eval, debug, RG_SigmaDelta);
  if (debug >= DefaultPrintThresholds::mostN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    print_intermediate(spinletter, "Xai_RG_SigmaDelta", RG_SigmaDelta, nvir, nocc_act);
  }
  // 1/2 \bar{\tilde{R}}^ao}_{\sigma\delta} \bar{g}^{\sigma\delta}_{io}
  double* RG_SigmaDelta2 = new double[nvir_occ];
  fill_n(RG_SigmaDelta2, nvir_occ, 0.0);
  compute_RG_SigmaDelta2(spin, C_0, C_1, r12eval, debug, RG_SigmaDelta2);
  if (debug >= DefaultPrintThresholds::mostN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    print_intermediate(spinletter, "Xai_RG_SigmaDelta2", RG_SigmaDelta2, nvir, nocc);
  }

  // compute contribution from
  // 1/2 \bar{\tilde{R}}^{ab'}_{ow} F^{alpha}_i \bar{\tilde{R}}^{ow}_{\alpha b'}
  // = 1/2 \bar{\tilde{R}}^{ab'}_{ow} \bar{\tilde{R}}^{ow}_{F_i_\alpha b'}
  // = 0
//  double* RFR1_ow = new double[nvir_occ];
//  fill_n(RFR1_ow, nvir_occ, 0.0);
//  compute_RFR_ow(spin, "vir", C_0, C_1, r12eval, debug, RFR1_ow);
//  print_intermediate(spinletter, "Xai_RFR1_ow", RFR1_ow, nvir, nocc);

  double* RFR_ow = new double[nvir_occ];
  fill_n(RFR_ow, nvir_occ, 0.0);
  compute_RFR_ow(spin, "cabs", C_0, C_1, r12eval, debug, RFR_ow);
  if (debug >= DefaultPrintThresholds::mostN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
    print_intermediate(spinletter, "Xai_RFR_ow", RFR_ow, nvir, nocc);
  }

  // compute contribution from
  // 1/2 \bar{\tilde{R}}^{ow}_{\alpha\beta} F^a_o \bar{\tilde{R}}^{\alpha\beta}_iw
  // = 1/2 \bar{\tilde{R}}^{F_a_o}_{\alpha\beta} \bar{\tilde{R}}^{\alpha\beta}_iw
  // = 0
//  double* Xai_RFR_AlphaBeta = new double[nvir_occ];
//  fill_n(Xai_RFR_AlphaBeta, nvir_occ, 0.0);
//  compute_RFR_AlphaBeta(spin, C_0, C_1, r12eval, debug, Xai_RFR_AlphaBeta);

  const double* iter_RG1 = RGow;
  const double* iter_RG2 = RG_SigmaDelta;
  const double* iter_RG12 = RGow2;
  const double* iter_RG22 = RG_SigmaDelta2;
  const double* iter_RFR_ow = RFR_ow;
//  const double* iter_RFR1_ow = RFR1_ow;
//  const double* iter_Xai_RFR_AlphaBeta = Xai_RFR_AlphaBeta;
  const int nfzc = nocc - nocc_act;
  for (int a = 0;  a < nvir; ++a) {
    for (int i = 0; i < nfzc; ++i) {
//      const double element_ai = *iter_RG1 - *iter_RG22;
//                              //+ *iter_RFR_ow * 2.0;
//      Xai.set_element(a,i,element_ai);

      ++iter_RG1;
      ++iter_RG22;
      ++iter_RFR_ow;
    }

    for (int i = nfzc; i < nocc; ++i) {
      const double element_ai = *iter_RG1 + *iter_RG12
                              - *iter_RG2 - *iter_RG22;
                              //+ *iter_RFR_ow * 2.0;
      Xai.set_element(a,i,element_ai);

      ++iter_RG1;
      ++iter_RG2;
      ++iter_RG12;
      ++iter_RG22;
      ++iter_RFR_ow;
    }
  }
  delete [] RGow;
  delete [] RG_SigmaDelta;
  delete [] RFR_ow;

  delete [] RGow2;
  delete [] RG_SigmaDelta2;

  if (debug >= DefaultPrintThresholds::mostN2) {
    Xai.print(prepend_spincase(spin,"F12 Xai:").c_str());
  }

 return Xai;

}


// compute different dg contraction contributions to Xai (CABS Singles)
// open-shell version
void compute_dg_contri(const Ref<OrbitalSpace>& vir, const Ref<OrbitalSpace>& occ,
                       const Ref<OrbitalSpace>& orbs1, const Ref<OrbitalSpace>& orbs2,
                       const Ref<OrbitalSpace>& orbs1_ab, const Ref<OrbitalSpace>& orbs2_ab,
                       const RefSCMatrix& D_CABS, const RefSCMatrix& D_CABS2,
                       const int offset_orbs1, const int offset_orbs2,
                       const int offset_orbs1_ab, const int offset_orbs2_ab,
                       const Ref<R12IntEval>& r12eval, double* const dg_contri) {

  Ref<R12WavefunctionWorld> r12world = r12eval->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);

  const TwoBodyOper::type eri_type = r12world->r12tech()->corrfactor()->tbint_type_eri();
  const unsigned int eri_idx = descr_f12->intset(eri_type);

  const int nocc = occ->rank();
  const int nvir = vir->rank();
  const int nvir_occ = nvir * nocc;

  // same spin part
  double* const dg1 = new double[nvir_occ];
  double* const dg2 = new double[nvir_occ];
  fill_n(dg1, nvir_occ, 0.0);
  fill_n(dg2, nvir_occ, 0.0);

  // activate integrals
  Ref<DistArray4> g_ints1;
  Ref<DistArray4> g_ints2;
  activate_ints(vir->id(), orbs1->id(), occ->id(), orbs2->id(),
                descr_f12_key, moints4_rtime, g_ints1);
  activate_ints(orbs1->id(), vir->id(), occ->id(), orbs2->id(),
                descr_f12_key, moints4_rtime, g_ints2);

  // g^a orbs1_i orbs2 * d^orbs2_orbs1
  compute_GDai("ai", g_ints1, eri_idx,
               D_CABS, offset_orbs1, offset_orbs2,
               dg1);
  //print_intermediate("same spin: ", "g^a orbs1_i orbs2 * d^orbs2_orbs1", dg1, nvir, nocc);
  // g^orbs1 a_i orbs2 * d^orbs2_orbs1
  compute_GDai("ia", g_ints2, eri_idx,
               D_CABS, offset_orbs1, offset_orbs2,
               dg2);
  //print_intermediate("same spin: ", "g^orbs1 a_i orbs2 * d^orbs2_orbs1", dg2, nvir, nocc);

  // test code: g^q'a_ip' * d^p'_q'
  #if 0
  {
  const int ncabs = orbs1->rank();
  const int norbs = nocc + nvir;

  double* const dg = new double[nvir_occ];
  fill_n(dg, nvir_occ, 0.0);

  double* iter = dg;
  for(int a = 0; a < nvir; ++a) {
    for(int i = 0; i < nocc; ++i) {

      double sum_pq = 0.0;
      for(int qp = 0;  qp < ncabs; ++qp) {
        const double* g_qpa_blk = g_ints2->retrieve_pair_block(qp, a, eri_idx);

        for(int pp = 0;  pp < ncabs; ++pp) {
            sum_pq += g_qpa_blk[i * ncabs + pp]
                    * D_CABS.get_element(pp + norbs, qp + norbs);
        }
        g_ints2->release_pair_block(qp, a, eri_idx);
      }

      *iter = sum_pq;
      ++iter;
    }
  }
  print_intermediate("same spin", "g^q'a_ip' * d^p'_q'", dg, nvir, nocc);
  delete[] dg;

  }
  #endif

  g_ints1->deactivate();
  g_ints2->deactivate();

  // AlphaBeta or BetaAlpha part
  double* const dg_ab = new double[nvir_occ];
  fill_n(dg_ab, nvir_occ, 0.0);

  Ref<DistArray4> g_ints_ab;
  activate_ints(vir->id(), orbs1_ab->id(), occ->id(), orbs2_ab->id(),
                descr_f12_key, moints4_rtime, g_ints_ab);

  // g^a orbs1_i orbs2 * d^orbs2_orbs1
  compute_GDai("ai", g_ints_ab, eri_idx,
               D_CABS2, offset_orbs1_ab, offset_orbs2_ab,
               dg_ab);
  g_ints_ab->deactivate();
  //print_intermediate("different spin: ", "g^a orbs1_i orbs2 * d^orbs2_orbs1", dg1, nvir, nocc);

  double* iter_dg = dg_contri;
  const double* iter_dg1 = dg1;
  const double* iter_dg2 = dg2;
  const double* iter_dg_ab = dg_ab;
  for(int a = 0; a < nvir; ++a) {
    for(int i = 0; i < nocc; ++i) {
      *iter_dg = *iter_dg1 - *iter_dg2 + *iter_dg_ab;

      ++iter_dg;
      ++iter_dg1;
      ++iter_dg2;
      ++iter_dg_ab;
    }
  }
  delete[] dg1;
  delete[] dg2;
  delete[] dg_ab;
}
// open-shell version


void compute_dg_contri(const Ref<OrbitalSpace>& vir, const Ref<OrbitalSpace>& occ,
                       const Ref<OrbitalSpace>& orbs1, const Ref<OrbitalSpace>& orbs2,
                       const RefSCMatrix D_CABS,
                       const int offset_orbs1, const int offset_orbs2,
                       const Ref<R12IntEval>& r12eval, double* const dg_contri) {

  Ref<R12WavefunctionWorld> r12world = r12eval->r12world();
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);

  const TwoBodyOper::type eri_type = r12world->r12tech()->corrfactor()->tbint_type_eri();
  const unsigned int eri_idx = descr_f12->intset(eri_type);

  const int nocc = occ->rank();
  const int nvir = vir->rank();
  const int nvir_occ = nvir * nocc;

  // same spin part
  double* dg1 = new double[nvir_occ];
  double* dg2 = new double[nvir_occ];
  fill_n(dg1, nvir_occ, 0.0);
  fill_n(dg2, nvir_occ, 0.0);

  // activate integrals
  Ref<DistArray4> g_ints1;
  Ref<DistArray4> g_ints2;
  activate_ints(vir->id(), orbs1->id(), occ->id(), orbs2->id(),
                descr_f12_key, moints4_rtime, g_ints1);
  activate_ints(orbs1->id(), vir->id(), occ->id(), orbs2->id(),
                descr_f12_key, moints4_rtime, g_ints2);

  // g^a orbs1_i orbs2 * d^orbs2_orbs1
  compute_GDai("ai", g_ints1, eri_idx,
               D_CABS, offset_orbs1, offset_orbs2,
               dg1);
  // g^orbs1 a_i orbs2 * d^orbs2_orbs1
  compute_GDai("ia", g_ints2, eri_idx,
               D_CABS, offset_orbs1, offset_orbs2,
               dg2);
  g_ints1->deactivate();
  g_ints2->deactivate();

  double* iter_dg = dg_contri;
  const double* iter_dg1 = dg1;
  const double* iter_dg2 = dg2;
  for(int a = 0; a < nvir; ++a) {
    for(int i = 0; i < nocc; ++i) {
      *iter_dg = (*iter_dg1) * 2 - *iter_dg2;

      ++iter_dg;
      ++iter_dg1;
      ++iter_dg2;
    }
  }
  delete[] dg1;
  delete[] dg2;
}
// close-shell version

// compute orbital relaxation contributions from CABS Singles
RefSCMatrix sc::Onerdm_X_CABS_Singles(SpinCase1 spin,
                                      const Ref<R12IntEval>& r12eval,
                                      const Ref<R12EnergyIntermediates>& r12intermediates,
                                      int debug) {

  Ref<R12WavefunctionWorld> r12world = r12eval->r12world();

  const Ref<OrbitalSpace> occ = r12eval->occ(spin);
  const Ref<OrbitalSpace> vir = r12eval->vir(spin);
  const Ref<OrbitalSpace> cabs = r12world->cabs_space(spin);
  const int nocc = occ->rank();
  const int nvir = vir->rank();
  const int norbs = nocc + nvir;
  const int nvir_occ = nvir * nocc;

  RefSCDimension rowdim = new SCDimension(nvir);
  RefSCDimension coldim = new SCDimension(nocc);
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix Xai = localkit->matrix(rowdim, coldim);
  Xai.assign(0.0);

  RefSCMatrix D_CABS = compute_D_CABS(spin, r12eval, r12intermediates);
  //D_CABS.print("CABS 1e density");

  // compute contribution from
  // d^p'_q' g^aq'_ip' & d^p'_q' g^q'a_ip'
  double* const dg_cabs = new double[nvir_occ];
  fill_n(dg_cabs, nvir_occ, 0.0);

  // d^p'_c g^ac_ip' & d^p'_c g^ca_ip'
  double* const dg_vir_cabs = new double[nvir_occ];
  fill_n(dg_vir_cabs, nvir_occ, 0.0);

  // d^c_q' g^aq'_ic & d^c_q' g^q'a_ic
  double* const dg_cabs_vir = new double[nvir_occ];
  fill_n(dg_cabs_vir, nvir_occ, 0.0);

  // d^k_l g^al_ik & d^k_l g^la_ik
  double* const dg_occ = new double[nvir_occ];
  fill_n(dg_occ, nvir_occ, 0.0);

  // d^p'_k g^ak_ip' & d^p'_k g^ka_ip'
  double* const dg_occ_cabs = new double[nvir_occ];
  fill_n(dg_occ_cabs, nvir_occ, 0.0);

  // d^k_p' g^ap'_ik & d^k_p' g^p'a_ik
  double* const dg_cabs_occ = new double[nvir_occ];
  fill_n(dg_cabs_occ, nvir_occ, 0.0);

  // F^a_a' * D^a'_i
  double* const dg_fd1 = new double[nvir_occ];
  fill_n(dg_fd1, nvir_occ, 0.0);
  RefSCMatrix FaA = r12eval->fock(vir,cabs,spin);
//  FaA.print("F^a_A'");
  compute_FDai(nvir, nocc,
               "a", FaA,
               D_CABS, norbs,
               dg_fd1);
//  // d^p'_i F^a_p' (test)
//  double* const dg_fd1 = new double[nvir_occ];
//  fill_n(dg_fd1, nvir_occ, 0.0);
//  compute_FDai(FaA, D_CABS, nocc, norbs, dg_fd1);

  // F^a'_i * D^a_a' (CC level)
  double* const dg_fd2 = new double[nvir_occ];
  fill_n(dg_fd2, nvir_occ, 0.0);
  RefSCMatrix FAi = r12eval->fock(cabs,occ,spin);
//  FAi.print("F^A_i");
  compute_FDai(nvir, nocc,
               "i", FAi,
               D_CABS, norbs,
               dg_fd2);
//  // d^a_p' F^p'_i (test)
//  double* const dg_fd2 = new double[nvir_occ];
//  fill_n(dg_fd2, nvir_occ, 0.0);
//  RefSCMatrix FiA = r12eval->fock(occ,cabs,spin);
//  compute_FDai(FiA, D_CABS, nvir, nocc, norbs, dg_fd2);

//  // F^a_m * D^m_i
//  double* const dg_fd3 = new double[nvir_occ];
//  fill_n(dg_fd3, nvir_occ, 0.0);
//  RefSCMatrix Fam = r12eval->fock(vir,occ,spin);
//  Fam.print("F^a_m");
//  compute_FDai(nvir, nocc,
//               "a", Fam,
//               D_CABS, 0,
//               dg_fd3);

  // F^m_i * D^a_m (MP2-R12 level)
  double* const dg_fd3 = new double[nvir_occ];
  fill_n(dg_fd3, nvir_occ, 0.0);
  RefSCMatrix Fmi = r12eval->fock(occ,occ,spin);
//  Fmi.print("Fmi");
//  D_CABS.print("D CABS Singles");
  compute_FDai(nvir, nocc,
               "i", Fmi,
               D_CABS, 0,
               dg_fd3);

  // F^a_c * D^c_i (MP2-R12 level)
  double* const dg_fd4 = new double[nvir_occ];
  fill_n(dg_fd4, nvir_occ, 0.0);
  RefSCMatrix Fac = r12eval->fock(vir,vir,spin);
  //Fac.print("Fac");
  compute_FDai(nvir, nocc,
               "a", Fac,
               D_CABS, nocc,
               dg_fd4);

  const int nspincases2 = (r12eval->spin_polarized() ? 3 : 2);
  if (nspincases2 == 3) {
    // obtain the opposite spin
    const SpinCase1 spin2 = (spin == Alpha? Beta : Alpha);

    const Ref<OrbitalSpace> occ2 = r12eval->occ(spin2);
    const Ref<OrbitalSpace> vir2 = r12eval->vir(spin2);
    const Ref<OrbitalSpace> cabs2 = r12world->cabs_space(spin2);
    const int nocc2 = occ2->rank();
    const int nvir2 = vir2->rank();
    const int norbs2 = nocc2 + nvir2;

    RefSCMatrix D_CABS2 = compute_D_CABS(spin2, r12eval, r12intermediates);

    compute_dg_contri(vir, occ, cabs, cabs,
                      cabs2, cabs2,
                      D_CABS, D_CABS2,
                      norbs, norbs, norbs2, norbs2,
                      r12eval, dg_cabs);

    compute_dg_contri(vir, occ, vir, cabs,
                      vir2, cabs2,
                      D_CABS, D_CABS2,
                      nocc, norbs, nocc2, norbs2,
                      r12eval, dg_vir_cabs);

    compute_dg_contri(vir, occ, cabs, vir,
                      cabs2, vir2,
                      D_CABS, D_CABS2,
                      norbs, nocc, norbs2, nocc2,
                      r12eval, dg_cabs_vir);

    const int zero = 0;
    compute_dg_contri(vir, occ, occ, occ,
                      occ2, occ2,
                      D_CABS, D_CABS2,
                      zero, zero, zero, zero,
                      r12eval, dg_occ);

    compute_dg_contri(vir, occ, occ, cabs,
                      occ2, cabs2,
                      D_CABS, D_CABS2,
                      zero, norbs, zero, norbs2,
                      r12eval, dg_occ_cabs);

    compute_dg_contri(vir, occ, cabs, occ,
                      cabs2, occ2,
                      D_CABS, D_CABS2,
                      norbs, zero, norbs2, zero,
                      r12eval, dg_cabs_occ);
  } else {
      compute_dg_contri(vir, occ, cabs, cabs,
                        D_CABS, norbs, norbs,
                        r12eval, dg_cabs);

      compute_dg_contri(vir, occ, vir, cabs,
                        D_CABS, nocc, norbs,
                        r12eval, dg_vir_cabs);

      compute_dg_contri(vir, occ, cabs, vir,
                        D_CABS, norbs, nocc,
                        r12eval, dg_cabs_vir);

      const int zero = 0;
      compute_dg_contri(vir, occ, occ, occ,
                        D_CABS, zero, zero,
                        r12eval, dg_occ);

      compute_dg_contri(vir, occ, occ, cabs,
                        D_CABS, zero, norbs,
                        r12eval, dg_occ_cabs);

      compute_dg_contri(vir, occ, cabs, occ,
                        D_CABS, norbs, zero,
                        r12eval, dg_cabs_occ);
  }

  if (debug >= DefaultPrintThresholds::mostN2) {
    const string spinletter = (spin == Alpha? "Alpha" : "Beta");
//    print_intermediate(spinletter, "CABS_Singles Xai_cabs", dg_cabs, nvir, nocc);
//    print_intermediate(spinletter, "CABS_Singles Xai_occ", dg_occ, nvir, nocc);
//    print_intermediate(spinletter, "CABS_Singles Xai_occ_cabs", dg_occ_cabs, nvir, nocc);
//    print_intermediate(spinletter, "CABS_Singles Xai_cabs_occ", dg_cabs_occ, nvir, nocc);
    print_intermediate(spinletter, "CABS_Singles F^a_a' * D^a'_i", dg_fd1, nvir, nocc);
    print_intermediate(spinletter, "CABS_Singles F^a'_i * D^a_a'", dg_fd2, nvir, nocc);
//    print_intermediate(spinletter, "CABS_Singles F^a_m * D^m_i", dg_fd3, nvir, nocc);
    print_intermediate(spinletter, "CABS_Singles F^m_i * D^a_m", dg_fd3, nvir, nocc);
    print_intermediate(spinletter, "CABS_Singles F^a_c * D^c_i", dg_fd4, nvir, nocc);
  }

  const double* iter_dg_cabs = dg_cabs;
  const double* iter_dg_vir_cabs = dg_vir_cabs;
  const double* iter_dg_cabs_vir = dg_cabs_vir;
  const double* iter_dg_occ = dg_occ;
  const double* iter_dg_occ_cabs = dg_occ_cabs;
  const double* iter_dg_cabs_occ = dg_cabs_occ;
  const double* iter_dg_fd1 = dg_fd1;
  const double* iter_dg_fd2 = dg_fd2;
  const double* iter_dg_fd3 = dg_fd3;
  const double* iter_dg_fd4 = dg_fd4;
  for (int a = 0;  a < nvir; ++a) {
    for (int i = 0; i < nocc; ++i) {
      const double element_ai =  *iter_dg_cabs
                               + *iter_dg_vir_cabs + *iter_dg_cabs_vir
//                               - *iter_dg_occ
//                               + *iter_dg_occ_cabs
//                               + *iter_dg_cabs_occ
//                               2.0 *(
//                                 *iter_dg_fd1
//                               + *iter_dg_fd4
//                               - *iter_dg_fd2
//                               - *iter_dg_fd3
//                               )
                               ;
      Xai.set_element(a,i,element_ai);

      ++iter_dg_cabs;
      ++iter_dg_vir_cabs;
      ++iter_dg_cabs_vir;
      ++iter_dg_occ;
      ++iter_dg_occ_cabs;
      ++iter_dg_cabs_occ;
//      ++iter_dg_fd1;
//      ++iter_dg_fd2;
//      ++iter_dg_fd3;
//      ++iter_dg_fd4;
    }
  }

  delete[] dg_cabs;
  delete[] dg_vir_cabs;
  delete[] dg_cabs_vir;
  delete[] dg_occ;
  delete[] dg_occ_cabs;
  delete[] dg_cabs_occ;
  delete[] dg_fd1;
  delete[] dg_fd2;
  delete[] dg_fd3;
  delete[] dg_fd4;

  return Xai;
}
