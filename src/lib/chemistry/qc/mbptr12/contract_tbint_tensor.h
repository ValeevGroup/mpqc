//
// contract_tbint_tensor.h
//
// Copyright (C) 2005 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_contracttbinttensor_h
#define _chemistry_qc_mbptr12_contracttbinttensor_h

#include <unistd.h>
#include <cmath>
#include <util/misc/regtime.h>
#include <util/misc/consumableresources.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/lcao/utils.h>
#include <chemistry/qc/lcao/utils.impl.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <util/misc/print.h>

namespace sc {

  template<typename DataProcess_Bra, typename DataProcess_Ket,
      typename DataProcess_BraKet, bool CorrFactorInBra, bool CorrFactorInKet,
      bool CorrFactorInInt>
  void R12IntEval::contract_tbint_tensor(
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
                                         const std::vector<std::string>& tformkeys_ket) {
    // are external spaces of particles 1 and 2 equivalent?
    const bool part1_strong_equiv_part2 = (space1_bra == space2_bra
        && space1_ket == space2_ket);
    // can external spaces of particles 1 and 2 be equivalent?
    const bool part1_weak_equiv_part2 = (space1_bra->rank()
        == space2_bra->rank() && space1_ket->rank() == space2_ket->rank());
    // Check correct semantics of this call : if antisymmetrize then particles must be equivalent
    bool correct_semantics = (antisymmetrize && (part1_weak_equiv_part2
        || part1_strong_equiv_part2)) || !antisymmetrize;
    // also
    correct_semantics
        = (correct_semantics && (space1_intb->rank() == space1_intk->rank())
            && (space2_intb->rank() == space2_intk->rank()));
    // also dimensions of tpcontract must match those of space1_int and space2_int
    correct_semantics = (correct_semantics && (tpcontract->nrow()
        == space1_intb->rank()) && (tpcontract->ncol() == space2_intb->rank()));
    if (!correct_semantics)
      throw ProgrammingError(
                             "R12IntEval::contract_tbint_tensor_() -- incorrect call semantics",
                             __FILE__, __LINE__);

    //
    // How is permutational symmetry implemented?
    //
    // 1) if need to antisymmetrize && internal spaces for p1 and p2 are same, then
    // can antisymmetrize each integral explicitly and compute antisymmetric tensor
    // 2) inf need to antisymmetrize but internal spaces for p1 and p2 do not match,
    // then compute as AlphaBeta and antisymmetrize at the end. I have to allocate temporary
    // result.
    //

    // are internal spaces of particles 1 and 2 equivalent?
    const bool part1_intequiv_part2 = (space1_intb == space2_intb
        && space1_intk == space2_intk);
#if 0
    // antisymmetrization for weakly equivalent particles and nonmatching internal spaces
    // is probably incorrect semantics
    if (!part1_intequiv_part2 && ! part1_strong_equiv_part2 && antisymmetrize)
    throw ProgrammingError("R12IntEval::contract_tbint_tensor_() -- dubious call semantics",
        __FILE__,__LINE__);
#endif
    // Will antisymmetrize each integral? If no, then result will be computed
    // as AlphaBeta and antisymmetrized at the end
    const bool alphabeta = !(antisymmetrize && part1_strong_equiv_part2
        && part1_intequiv_part2);

    //
    // NOTE! Even if computing in AlphaBeta, internal sums can be over AlphaAlpha!!!
    // Logic should not become much more complicated. Only need time to implement.
    //

    const bool CorrFactorInBraInt = CorrFactorInBra && CorrFactorInInt;
    const bool CorrFactorInKetInt = CorrFactorInKet && CorrFactorInInt;

    const unsigned int nbrasets = (CorrFactorInBra ? corrfactor()->nfunctions()
                                                   : 1);
    const unsigned int nketsets = (CorrFactorInKet ? corrfactor()->nfunctions()
                                                   : 1);
    const unsigned int nintsets = (CorrFactorInInt ? corrfactor()->nfunctions()
                                                   : 1);
    const unsigned int nbraintsets = (CorrFactorInBraInt ? UINT_MAX : nbrasets
        * nintsets);
    const unsigned int nketintsets = (CorrFactorInKetInt ? UINT_MAX : nketsets
        * nintsets);

    //
    // create transforms, if needed
    //
    typedef std::vector<Ref<TwoBodyMOIntsTransform> > tformvec;

    // bra transforms
    const size_t num_tforms_bra = tformkeys_bra.size();
    tformvec transforms_bra(num_tforms_bra);
    for (unsigned int t = 0; t < num_tforms_bra; ++t) {
      transforms_bra[t] = moints_runtime4()->get(tformkeys_bra[t]);
    }

    // ket transforms
    const size_t num_tforms_ket = tformkeys_ket.size();
    tformvec transforms_ket(num_tforms_ket);
    for (unsigned int t = 0; t < num_tforms_ket; ++t) {
      transforms_ket[t] = moints_runtime4()->get(tformkeys_ket[t]);
    }

    //
    // Generate contract label
    //
    Timer tim_gen_tensor_contract("Generic tensor contract");
    std::string label;
    {
      std::ostringstream oss_bra;
      oss_bra << "<" << space1_bra->id() << " " << space2_bra->id()
          << (antisymmetrize ? "||" : "|") << space1_intb->id() << " "
          << space2_intb->id() << ">";
      const std::string label_bra = oss_bra.str();
      std::ostringstream oss_ket;
      oss_ket << "<" << space1_ket->id() << " " << space2_ket->id()
          << (antisymmetrize ? "||" : "|") << space1_intk->id() << " "
          << space2_intk->id() << ">";
      const std::string label_ket = oss_ket.str();
      std::ostringstream oss;
      oss << "<" << space1_bra->id() << " " << space2_bra->id()
          << (antisymmetrize ? "||" : "|") << space1_ket->id() << " "
          << space2_ket->id() << "> = " << label_bra << " . " << label_ket
          << "^T";
      label = oss.str();
    }
    ExEnv::out0() << std::endl << indent << "Entered generic contraction (" << label
        << ")" << std::endl;
    ExEnv::out0() << incindent;

    //
    // Construct maps
    //
    // WARNING: Assuming all transforms are over same spaces!!!
    //
    Ref<OrbitalSpace> tspace1_bra = transforms_bra[0]->space1();
    Ref<OrbitalSpace> tspace2_bra = transforms_bra[0]->space3();
    Ref<OrbitalSpace> tspace1_intb = transforms_bra[0]->space2();
    Ref<OrbitalSpace> tspace2_intb = transforms_bra[0]->space4();
    Ref<OrbitalSpace> tspace1_ket = transforms_ket[0]->space1();
    Ref<OrbitalSpace> tspace2_ket = transforms_ket[0]->space3();
    Ref<OrbitalSpace> tspace1_intk = transforms_ket[0]->space2();
    Ref<OrbitalSpace> tspace2_intk = transforms_ket[0]->space4();
    // maps spaceX to spaceX of the transform
    std::vector<unsigned int> map1_bra, map2_bra, map1_ket, map2_ket,
        map1_intb, map2_intb, map1_intk, map2_intk;
    // maps space2_intb to space1_intb of transform
    std::vector<unsigned int> map12_intb;
    // maps space1_intb to space2_intb of transform
    std::vector<unsigned int> map21_intb;
    // maps space2_intk to space1_intk of transform
    std::vector<unsigned int> map12_intk;
    // maps space1_intk to space2_intk of transform
    std::vector<unsigned int> map21_intk;

    { // bra maps
      map1_bra = *tspace1_bra << *space1_bra;
      map2_bra = *tspace2_bra << *space2_bra;
      map1_intb = *tspace1_intb << *space1_intb;
      map2_intb = *tspace2_intb << *space2_intb;
      // Will antisymmetrize the integrals? Then need ijkl AND ijlk
      if (!alphabeta) {
        if (tspace1_intb == tspace2_intb) {
          map12_intb = map1_intb;
          map21_intb = map2_intb;
        } else {
          map12_intb = *tspace1_intb << *space2_intb;
          map21_intb = *tspace2_intb << *space1_intb;
        }
      }
    }
    { // ket maps
      map1_ket = *tspace1_ket << *space1_ket;
      map2_ket = *tspace2_ket << *space2_ket;
      map1_intk = *tspace1_intk << *space1_intk;
      map2_intk = *tspace2_intk << *space2_intk;
      // Will antisymmetrize the integrals? Then need ijkl AND ijlk
      if (!alphabeta) {
        if (tspace1_intk == tspace2_intk) {
          map12_intk = map1_intk;
          map21_intk = map2_intk;
        } else {
          map12_intk = *tspace1_intk << *space2_intk;
          map21_intk = *tspace2_intk << *space1_intk;
        }
      }
    }

    const unsigned int bratform_block_ncols = tspace2_intb->rank();
    const unsigned int kettform_block_ncols = tspace2_intk->rank();
    const RefDiagSCMatrix evals1_bra = space1_bra->evals();
    const RefDiagSCMatrix evals2_bra = space2_bra->evals();
    const RefDiagSCMatrix evals1_ket = space1_ket->evals();
    const RefDiagSCMatrix evals2_ket = space2_ket->evals();
    const RefDiagSCMatrix evals1_intb = space1_intb->evals();
    const RefDiagSCMatrix evals2_intb = space2_intb->evals();
    const RefDiagSCMatrix evals1_intk = space1_intk->evals();
    const RefDiagSCMatrix evals2_intk = space2_intk->evals();

    // Using spinorbital iterators means I don't take into account perm symmetry
    // More efficient algorithm will require generic code
    const SpinCase2 S = (alphabeta ? AlphaBeta : AlphaAlpha);
    SpinMOPairIter
        iterbra(space1_bra->rank(), (alphabeta ? space2_bra : space1_bra)->rank(), S);
    SpinMOPairIter
        iterket(space1_ket->rank(), (alphabeta ? space2_ket : space1_ket)->rank(), S);
    SpinMOPairIter iterint(space1_intb->rank(),
                           (alphabeta ? space2_intb : space1_intb)->rank(), S);
    // size of one block of <space1_bra space2_bra|
    const unsigned int nbra = iterbra.nij();
    // size of one block of <space1_ket space2_ket|
    const unsigned int nket = iterket.nij();
    // size of one block of |space1_int space2_int>
    const unsigned int nint = iterint.nij();

    RefSCMatrix Tcontr;
    // Allocate storage for the result, if need to antisymmetrize at the end; else accumulate directly to T
    if (antisymmetrize && alphabeta) {
      Tcontr = T.kit()->matrix(new SCDimension(nbra * nbrasets),
                               new SCDimension(nket * nketsets));
      Tcontr.assign(0.0);
    } else
      Tcontr = T;

    // size of integral blocks to contract
    const unsigned int blksize_int = space1_intb->rank() * space2_intb->rank();
    double* T_ij = new double[blksize_int];
    double* T_kl = new double[blksize_int];

    //
    // Contraction loops
    //

    // outermost loop over contraction blocks to minimize number of activate()/deactivate() calls
    for (unsigned int fint = 0; fint < nintsets; ++fint) {

      unsigned int fbraoffset = 0;
      for (unsigned int fbra = 0; fbra < nbrasets; ++fbra, fbraoffset += nbra) {
        const unsigned int fbraint = fbra * nintsets + fint;
        Ref<TwoBodyMOIntsTransform> tformb = transforms_bra[fbraint];
        const Ref<TwoBodyIntDescr>& intdescrb = tformb->intdescr();
        const unsigned int intsetidx_bra = intdescrb->intset(tbint_type_bra);

        tformb->compute();
        Ref<DistArray4> accumb = tformb->ints_distarray4();
        accumb->activate();

        unsigned int fketoffset = 0;
        for (unsigned int fket = 0; fket < nketsets; ++fket, fketoffset += nket) {
          const unsigned int fketint = fket * nintsets + fint;
          Ref<TwoBodyMOIntsTransform> tformk = transforms_ket[fketint];
          const Ref<TwoBodyIntDescr>& intdescrk = tformk->intdescr();
          const unsigned int intsetidx_ket = intdescrk->intset(tbint_type_ket);

          tformk->compute();
          Ref<DistArray4> accumk = tformk->ints_distarray4();
          accumk->activate();

          if (debug_ >= DefaultPrintThresholds::diagnostics) {
            ExEnv::out0() << indent << "Using transforms " << tformb->name()
                << " and " << tformk->name() << std::endl;
          }

          // split work over tasks which have access to integrals
          // WARNING: assuming same accessibility for both bra and ket transforms
          std::vector<int> proc_with_ints;
          const int nproc_with_ints = accumb->tasks_with_access(proc_with_ints);
          const int me = r12world()->world()->msg()->me();
          const bool nket_ge_nevals = (nket >= nproc_with_ints);

          if (accumb->has_access(me)) {

            unsigned int ijkl = 0;
            for (iterbra.start(); iterbra; iterbra.next()) {
              const int ij = iterbra.ij();

              bool fetch_ij_block = false;
              if (nket_ge_nevals) {
                fetch_ij_block = true;
              } else {
                const int first_ij_task = ijkl % nproc_with_ints;
                const int last_ij_task = (ijkl + nket - 1) % nproc_with_ints;
                // figure out if for this ij there is ijkl to be handled by me
                if (!(first_ij_task > proc_with_ints[me] && proc_with_ints[me]
                    > last_ij_task))
                  fetch_ij_block = true;
              }
              if (!fetch_ij_block)
                continue;

              const unsigned int i = iterbra.i();
              const unsigned int j = iterbra.j();
              const unsigned int ii = map1_bra[i];
              const unsigned int jj = map2_bra[j];

              Timer tim_intsretrieve("MO ints retrieve");
              const double *ij_buf = accumb->retrieve_pair_block(ii, jj,
                                                                 intsetidx_bra);
              tim_intsretrieve.exit();
              if (debug_ >= DefaultPrintThresholds::allO2N2)
                ExEnv::outn() << indent << "task " << me
                    << ": obtained ij blocks" << std::endl;

              for (iterket.start(); iterket; iterket.next(), ijkl++) {
                const int kl = iterket.ij();

                const int ijkl_proc = ijkl % nproc_with_ints;
                if (ijkl_proc != proc_with_ints[me])
                  continue;

                const unsigned int k = iterket.i();
                const unsigned int l = iterket.j();
                const unsigned int kk = map1_ket[k];
                const unsigned int ll = map2_ket[l];

                if (debug_ >= DefaultPrintThresholds::allO2N2)
                  ExEnv::outn() << indent << "task " << me
                      << ": working on (i,j | k,l) = (" << i << "," << j
                      << " | " << k << "," << l << ")" << std::endl;
                tim_intsretrieve.enter("MO ints retrieve");
                const double *kl_buf =
                    accumk->retrieve_pair_block(kk, ll, intsetidx_ket);
                tim_intsretrieve.exit();
                if (debug_ >= DefaultPrintThresholds::allO2N2)
                  ExEnv::outn() << indent << "task " << me
                      << ": obtained kl blocks" << std::endl;

                // zero out intblocks
                memset(T_ij, 0, blksize_int * sizeof(double));
                memset(T_kl, 0, blksize_int * sizeof(double));

                if (debug_ >= DefaultPrintThresholds::allO2N2) {
                  ExEnv::out0() << indent << "i = " << i << " j = " << j
                      << " k = " << k << " l = " << l << incindent << std::endl;
                }

                for (iterint.start(); iterint; iterint.next()) {
                  const unsigned int m = iterint.i();
                  const unsigned int n = iterint.j();
                  const int mn = iterint.ij();
                  int MN_bra, MN_ket;
                  {
                    const unsigned int mm = map1_intb[m];
                    const unsigned int nn = map2_intb[n];
                    MN_bra = mm * bratform_block_ncols + nn;
                  }
                  {
                    const unsigned int mm = map1_intk[m];
                    const unsigned int nn = map2_intk[n];
                    MN_ket = mm * kettform_block_ncols + nn;
                  }

                  const double I_ijmn = ij_buf[MN_bra];
                  const double I_klmn = kl_buf[MN_ket];

                  if (debug_ >= DefaultPrintThresholds::mostO6) {
                    ExEnv::out0() << indent << " m = " << m << " n = " << n
                        << std::endl;
                  }

                  if (alphabeta) {
                    if (debug_ >= DefaultPrintThresholds::mostO6) {
                      ExEnv::out0() << indent << " <ij|mn> = " << I_ijmn
                          << std::endl << indent << " <kl|mn> = " << I_klmn << std::endl;
                    }
                    const double T_ijmn = DataProcess_Bra::I2T(I_ijmn, i, j, m,
                                                               n, evals1_bra,
                                                               evals1_intb,
                                                               evals2_bra,
                                                               evals2_intb);
                    const double T_klmn = DataProcess_Ket::I2T(I_klmn, k, l, m,
                                                               n, evals1_ket,
                                                               evals1_intk,
                                                               evals2_ket,
                                                               evals2_intk);
                    if (debug_ >= DefaultPrintThresholds::mostO6) {
                      ExEnv::out0() << indent << " <ij|T|mn> = " << T_ijmn
                          << std::endl << indent << " <kl|T|mn> = " << T_klmn
                          << std::endl;
                    }

                    T_ij[mn] = T_ijmn;
                    T_kl[mn] = T_klmn;
                  } else {

                    int NM_bra, NM_ket;
                    {
                      const unsigned int mm = map21_intb[m];
                      const unsigned int nn = map12_intb[n];
                      NM_bra = nn * bratform_block_ncols + mm;
                    }
                    {
                      const unsigned int mm = map21_intk[m];
                      const unsigned int nn = map12_intk[n];
                      NM_ket = nn * kettform_block_ncols + mm;
                    }

                    const double I_ijnm = ij_buf[NM_bra];
                    const double I_klnm = kl_buf[NM_ket];

                    if (debug_ >= DefaultPrintThresholds::mostO6) {
                      ExEnv::out0() << " <ij|mn> = " << I_ijmn << std::endl
                          << " <ij|nm> = " << I_ijnm << std::endl << " <kl|mn> = "
                          << I_klmn << std::endl << " <kl|nm> = " << I_klnm << std::endl;
                    }

                    const double T_ijmn = DataProcess_Bra::I2T(I_ijmn - I_ijnm,
                                                               i, j, m, n,
                                                               evals1_bra,
                                                               evals1_intb,
                                                               evals2_bra,
                                                               evals2_intb);
                    const double T_klmn = DataProcess_Ket::I2T(I_klmn - I_klnm,
                                                               k, l, m, n,
                                                               evals1_ket,
                                                               evals1_intk,
                                                               evals2_ket,
                                                               evals2_intk);
                    if (debug_ >= DefaultPrintThresholds::mostO6) {
                      ExEnv::out0() << indent << " <ij|T|mn> = " << T_ijmn
                          << std::endl << indent << " <kl|T|mn> = " << T_klmn
                          << std::endl;
                    }
                    T_ij[mn] = T_ijmn;
                    T_kl[mn] = T_klmn;
                  }

                } // int loop

                // contract matrices
                double T_ijkl = tpcontract->contract(T_ij, T_kl);
                if (debug_ >= DefaultPrintThresholds::allO2N2) {
                  ExEnv::out0() << decindent << indent << " <ij|kl> = "
                      << T_ijkl << std::endl;
                }
                T_ijkl = DataProcess_BraKet::I2T(T_ijkl, i, j, k, l,
                                                 evals1_bra, evals1_ket,
                                                 evals2_bra, evals2_ket);
                if (debug_ >= DefaultPrintThresholds::allO2N2) {
                  ExEnv::out0() << indent << " <ij|T|kl>(ref) = "
                      << Tcontr.get_element(ij + fbraoffset, kl + fketoffset)
                      << std::endl;
                  ExEnv::out0() << indent << " +<ij|T|kl> = " << T_ijkl << std::endl;
                }
                Tcontr.accumulate_element(ij + fbraoffset, kl + fketoffset,
                                          T_ijkl);
                if (debug_ >= DefaultPrintThresholds::allO2N2) {
                  ExEnv::out0() << indent << " <ij|T|kl>(new) = "
                      << Tcontr.get_element(ij + fbraoffset, kl + fketoffset)
                      << std::endl;
                }

                accumk->release_pair_block(kk, ll, intsetidx_ket);

              } // ket loop

              if (fetch_ij_block)
                accumb->release_pair_block(ii, jj, intsetidx_bra);

            } // bra loop
          } // loop over tasks with access

          //ExEnv::out0() << indent << "Accumb = " << accumb.pointer() << std::endl;
          //ExEnv::out0() << indent << "Accumk = " << accumk.pointer() << std::endl;
          //ExEnv::out0() << indent << "Accumb == Accumk : " << (accumb==accumk) << std::endl;
          if (accumb != accumk && accumk->data_persistent()) accumk->deactivate();

        } // ket blocks

        if (accumb->data_persistent()) accumb->deactivate();

      } // bra blocks
    } // int blocks

    if (antisymmetrize && alphabeta) {
      // antisymmetrization implies equivalent particles -- hence symmetrize before antisymmetrize
      symmetrize<false> (Tcontr, Tcontr, space1_bra, space1_ket);
      sc::antisymmetrize(T, Tcontr, space1_bra, space1_ket, true);
      Tcontr = 0;
    }

    delete[] T_ij;
    delete[] T_kl;

    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited generic contraction (" << label << ")"
        << std::endl;
    tim_gen_tensor_contract.exit();
  }

  template<bool CorrFactorInBra, bool CorrFactorInKet>
  void R12IntEval::contract_tbint_tensor(
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
                                         const std::vector<std::string>& tformkeys_ket) {

    MPQC_ASSERT(tformkeys_bra.size() == 1);
    MPQC_ASSERT(tformkeys_ket.size() == 1);

    std::vector< Ref<DistArray4> > results;
    contract_tbint_tensor<CorrFactorInBra, CorrFactorInKet>(
                          results,
                          tbint_type_bra,
                          tbint_type_ket,
                          scale,
                          space1_bra,
                          space2_bra,
                          space1_intb,
                          space2_intb,
                          space1_ket,
                          space2_ket,
                          space1_intk,
                          space2_intk,
                          antisymmetrize,
                          tformkeys_bra,
                          tformkeys_ket);

    RefSCMatrix Tclone = T.clone();
    Tclone << results[0];
    T.accumulate(Tclone);
  }

  template<bool CorrFactorInBra, bool CorrFactorInKet>
  void
  R12IntEval::contract_tbint_tensor(std::vector< Ref<DistArray4> >& results,
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
                                    const std::vector<std::string>& tformkeys_ket) {

    // are external spaces of particles 1 and 2 equivalent?
    const bool part1_strong_equiv_part2 = (space1_bra == space2_bra
        && space1_ket == space2_ket);
    // can external spaces of particles 1 and 2 be equivalent?
    const bool part1_weak_equiv_part2 = (space1_bra->rank()
        == space2_bra->rank() && space1_ket->rank() == space2_ket->rank());
    // Check correct semantics of this call : if antisymmetrize then particles must be equivalent
    bool correct_semantics = (antisymmetrize && (part1_weak_equiv_part2
        || part1_strong_equiv_part2)) || !antisymmetrize;
    // also
    correct_semantics
        = (correct_semantics && (space1_intb->rank() == space1_intk->rank())
            && (space2_intb->rank() == space2_intk->rank()));
    if (!correct_semantics)
      throw ProgrammingError(
                             "R12IntEval::contract_tbint_tensor_() -- incorrect call semantics",
                             __FILE__, __LINE__);

    //
    // How is permutational symmetry implemented?
    //
    // 1) if need to antisymmetrize && internal spaces for p1 and p2 are same, then
    // can antisymmetrize each integral explicitly and compute antisymmetric tensor:
    // \f$
    //   T_{ij}^{kl} = \sum_{p<q} (A_{ij}^{pq} - A_{ij}^{qp}) (B^{kl}_{pq} - B^{kl}_{qp})
    // \f$
    // 2) if need to antisymmetrize but internal spaces for p1 and p2 do not match,
    // then compute as AlphaBeta and antisymmetrize at the end. Temporary storage
    // will be allocated.
    //

    // are internal spaces of particles 1 and 2 equivalent?
    const bool part1_intequiv_part2 = (space1_intb == space2_intb
        && space1_intk == space2_intk);
    // If antisymmetrize == true and particles 1 and 2 are equivalent, then the target is computed as:
    // T_{ij}^{kl} = \sum_{p<q} (A_{ij}^{pq} - A_{ij}^{qp}) (B^{kl}_{pq} - B^{kl}_{qp})
    // i.e. each block of A and B is antisymmetrized first (for each i<j, k<l).
    // In all other cases, even when antisymmetric target T tensor is requested,
    // compute alpha-beta T and then antisymmetrize the result, i.e.
    // V_{ij}^{kl} = \sum_{pq} A_{ij}^{pq} B^{kl}_{pq}
    // for all i,j,k,l, then
    // T_{ij}^{kl} = V_{ij}{kl} - V_{ij}^{lk}
    // for i<j, k<l
    // this flag determines whether to contract as alpha-beta. If not, each source integral will be antisymmetrized
    // first
    const bool alphabeta = !(antisymmetrize && part1_strong_equiv_part2
        && part1_intequiv_part2);

    //
    // NOTE! Even if computing in AlphaBeta, internal sums can be over AlphaAlpha!!!
    // Logic should not become much more complicated. Only need time to implement.
    //

    const unsigned int nbrasets = (CorrFactorInBra ? corrfactor()->nfunctions()
                                                   : 1);
    const unsigned int nketsets = (CorrFactorInKet ? corrfactor()->nfunctions()
                                                   : 1);

    //
    // get transforms
    //
    typedef std::vector<Ref<TwoBodyMOIntsTransform> > tformvec;

    // bra transforms
    const size_t num_tforms_bra = tformkeys_bra.size();
    tformvec transforms_bra(num_tforms_bra);
    for (unsigned int t = 0; t < num_tforms_bra; ++t) {
      transforms_bra[t] = moints_runtime4()->get(tformkeys_bra[t]);
    }

    // ket transforms
    const size_t num_tforms_ket = tformkeys_ket.size();
    tformvec transforms_ket(num_tforms_ket);
    for (unsigned int t = 0; t < num_tforms_ket; ++t) {
      transforms_ket[t] = moints_runtime4()->get(tformkeys_ket[t]);
    }

    // create result DistArray4 objects, if needed
    bool accumulate = false;   // need to accumulate result to a given DistArray4?
    std::vector< Ref<DistArray4> > results_clone;
    if (results.empty()) {
      std::vector<std::string> tformkeys_result;
      { // create result objects by cloning transforms_bra[0]
        transforms_bra[0]->compute();
        for (unsigned int tb = 0; tb < num_tforms_bra; ++tb) {
          const std::string key_bra = tformkeys_bra[tb];
          const ParsedTwoBodyFourCenterIntKey pkey_bra(key_bra);
          const std::string parb = pkey_bra.params();
          std::string params_bra;
          if (!parb.empty()) { // remove [ and ]
            const std::string::size_type lbpos = parb.find_first_of('[');
            const std::string::size_type rbpos = parb.find_last_of(']');
            MPQC_ASSERT(lbpos < rbpos);
            params_bra = parb.substr(lbpos + 1, (rbpos - lbpos - 1));
          }

          for (unsigned int tk = 0; tk < num_tforms_ket; ++tk) {
            const std::string key_ket = tformkeys_ket[tk];
            const ParsedTwoBodyFourCenterIntKey pkey_ket(key_ket);
            const std::string park = pkey_ket.params();
            std::string params_ket;
            if (!park.empty()) { // remove [ and ]
              const std::string::size_type lbpos = park.find_first_of('[');
              const std::string::size_type rbpos = park.find_last_of(']');
              MPQC_ASSERT(lbpos < rbpos);
              params_ket = park.substr(lbpos + 1, (rbpos - lbpos - 1));
            }

            std::string params_result;
            if (CorrFactorInBra && CorrFactorInKet)
              params_result = '[' + params_bra + ',' + params_ket + ']';
            else {
              if (CorrFactorInBra && !CorrFactorInKet)
                params_result = '[' + params_bra + ']';
              if (!CorrFactorInBra && CorrFactorInKet)
                params_result = '[' + params_ket + ']';
            }

            const std::string descr_key = pkey_bra.oper() + "@" + pkey_ket.oper() + params_result;

            const std::string key_result =
                ParsedTwoBodyFourCenterIntKey::key(space1_bra->id(),
                                                   space2_bra->id(),
                                                   space1_ket->id(),
                                                   space2_ket->id(), descr_key,
                                                   TwoBodyIntLayout::b1b2_k1k2);

            tformkeys_result.push_back(key_result);

            DistArray4Dimensions dims(1, space1_bra->rank(),
                                      space2_bra->rank(), space1_ket->rank(),
                                      space2_ket->rank());
            Ref<DistArray4> result = transforms_bra[0]->ints_distarray4()->clone(dims);
            results.push_back(result);
          }
        }
      }
    } else { // results were given -- must clone them, and accumulate after contraction (can' currently do accumulation in place)
      accumulate = true;
      // make sure they have expected shape
      const size_t nresults = results.size();
      MPQC_ASSERT(nresults == num_tforms_bra*num_tforms_ket);
      for(int r=0; r<nresults; ++r) {
        MPQC_ASSERT(results[r]->num_te_types() == 1);
        MPQC_ASSERT(results[r]->ni() == space1_bra->rank());
        MPQC_ASSERT(results[r]->nj() == space2_bra->rank());
        MPQC_ASSERT(results[r]->nx() == space1_ket->rank());
        MPQC_ASSERT(results[r]->ny() == space2_ket->rank());
        MPQC_ASSERT(results[r]->storage() == DistArray4Storage_XY);
      }

      for(int r=0; r<nresults; ++r) {
        results_clone.push_back( results[r]->clone() );
      }
    }
    std::vector< Ref<DistArray4> >& results_ref = accumulate ? results_clone : results;

    //
    // Generate contract label
    //
    Timer tim_gen_tensor_contract("Generic tensor contract");
    std::string label;
    {
      std::ostringstream oss_bra;
      oss_bra << "<" << space1_bra->id() << " " << space2_bra->id()
          << (antisymmetrize ? "||" : "|") << space1_intb->id() << " "
          << space2_intb->id() << ">";
      const std::string label_bra = oss_bra.str();
      std::ostringstream oss_ket;
      oss_ket << "<" << space1_ket->id() << " " << space2_ket->id()
          << (antisymmetrize ? "||" : "|") << space1_intk->id() << " "
          << space2_intk->id() << ">";
      const std::string label_ket = oss_ket.str();
      std::ostringstream oss;
      oss << "<" << space1_bra->id() << " " << space2_bra->id()
          << (antisymmetrize ? "||" : "|") << space1_ket->id() << " "
          << space2_ket->id() << "> = " << label_bra << " . " << label_ket
          << "^T";
      label = oss.str();
    }
    ExEnv::out0() << std::endl << indent << "Entered generic contraction (" << label
        << ")" << std::endl;
    ExEnv::out0() << incindent;

    //
    // make sure that these are the needed transforms
    //
    {
      Ref<OrbitalSpace> tspace1_bra = transforms_bra[0]->space1();
      Ref<OrbitalSpace> tspace2_bra = transforms_bra[0]->space3();
      Ref<OrbitalSpace> tspace1_intb = transforms_bra[0]->space2();
      Ref<OrbitalSpace> tspace2_intb = transforms_bra[0]->space4();
      Ref<OrbitalSpace> tspace1_ket = transforms_ket[0]->space1();
      Ref<OrbitalSpace> tspace2_ket = transforms_ket[0]->space3();
      Ref<OrbitalSpace> tspace1_intk = transforms_ket[0]->space2();
      Ref<OrbitalSpace> tspace2_intk = transforms_ket[0]->space4();
      MPQC_ASSERT(tspace1_bra == space1_bra);
      MPQC_ASSERT(tspace2_bra == space2_bra);
      MPQC_ASSERT(tspace1_ket == space1_ket);
      MPQC_ASSERT(tspace2_ket == space2_ket);
      MPQC_ASSERT(tspace1_intb == space1_intb);
      MPQC_ASSERT(tspace2_intb == space2_intb);
      MPQC_ASSERT(tspace1_intk == space1_intk);
      MPQC_ASSERT(tspace2_intk == space2_intk);
    }

    //
    // Contraction loops
    //

    unsigned int fbraoffset = 0;
    unsigned int fbraket = 0;
    for (unsigned int fbra = 0; fbra < nbrasets; ++fbra) {
      Ref<TwoBodyMOIntsTransform> tformb = transforms_bra[fbra];
      const Ref<TwoBodyIntDescr>& intdescrb = tformb->intdescr();
      const unsigned int intsetidx_bra = intdescrb->intset(tbint_type_bra);

      tformb->compute();
      Ref<DistArray4> accumb = tformb->ints_distarray4();

      unsigned int fketoffset = 0;
      for (unsigned int fket = 0; fket < nketsets; ++fket, ++fbraket) {
          Ref<TwoBodyMOIntsTransform> tformk = transforms_ket[fket];
          const Ref<TwoBodyIntDescr>& intdescrk = tformk->intdescr();
          const unsigned int intsetidx_ket = intdescrk->intset(tbint_type_ket);

          tformk->compute();
          Ref<DistArray4> accumk = tformk->ints_distarray4();


          if (debug_ >= DefaultPrintThresholds::diagnostics) {
            ExEnv::out0() << indent << "Using transforms " << tformb->name()
            << " and " << tformk->name() << std::endl;
          }

          sc::contract34(results_ref[fbraket],
                         scale,
                         accumb,
                         intsetidx_bra,
                         accumk,
                         intsetidx_ket,
                         debug_);

          if (antisymmetrize) {
            sc::antisymmetrize(results_ref[fbraket]);
          }

          if (accumulate)
            sc::axpy(results_ref[fbraket], 1.0, results[fbraket]);

          //ExEnv::out0() << indent << "Accumb = " << accumb.pointer() << std::endl;
          //ExEnv::out0() << indent << "Accumk = " << accumk.pointer() << std::endl;
          //ExEnv::out0() << indent << "Accumb == Accumk : " << (accumb==accumk) << std::endl;
      } // ket blocks

    } // bra blocks

    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited generic contraction (" << label << ")"
        << std::endl;
    tim_gen_tensor_contract.exit();

  }

}

#endif

