//
// compute_tbint_tensor.h
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

#ifdef __GNUG__
#pragma interface
#endif

#include <cmath>
#include <util/misc/regtime.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/utils.impl.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/print.h>

#ifndef _chemistry_qc_mbptr12_computetbinttensor_h
#define _chemistry_qc_mbptr12_computetbinttensor_h

namespace sc {

  template <typename DataProcess,
            bool CorrFactorInBra,
            bool CorrFactorInKet>
    void
    R12IntEval::compute_tbint_tensor(RefSCMatrix& T,
                                     TwoBodyOper::type tbint_type,
                                     const Ref<OrbitalSpace>& space1_bra,
                                     const Ref<OrbitalSpace>& space1_ket,
                                     const Ref<OrbitalSpace>& space2_bra,
                                     const Ref<OrbitalSpace>& space2_ket,
                                     bool antisymmetrize,
                                     const std::vector<std::string>& transform_keys)
    {
      // are particles 1 and 2 equivalent?
      const bool part1_strong_equiv_part2 = (space1_bra==space2_bra && space1_ket==space2_ket);
      const bool part1_weak_equiv_part2 = (space1_bra->rank()==space2_bra->rank() && space1_ket->rank()==space2_ket->rank());
      // Check correct semantics of this call : if antisymmetrize then particles must be equivalent
      const bool correct_semantics = (antisymmetrize && part1_weak_equiv_part2) ||
                                     !antisymmetrize;
      if (!correct_semantics)
        throw ProgrammingError("R12IntEval::compute_tbint_tensor_() -- incorrect call semantics",
                               __FILE__,__LINE__);

      // If need to antisymmetrize but particles are not truly equivalent, compute as AlphaBeta and antisymmetrize
      const bool alphabeta = !(antisymmetrize && part1_strong_equiv_part2);

      const bool CorrFactorInBraKet = CorrFactorInBra && CorrFactorInKet;

      const unsigned int nbrasets = (CorrFactorInBra ? corrfactor()->nfunctions() : 1);
      const unsigned int nketsets = (CorrFactorInKet ? corrfactor()->nfunctions() : 1);
      const unsigned int nsets = (CorrFactorInBraKet ? -1 : nbrasets*nketsets);

      // create transforms, if needed
      typedef std::vector< Ref<TwoBodyMOIntsTransform> > tformvec;
      const size_t num_tforms = transform_keys.size();
      tformvec transforms(num_tforms);
      for(unsigned int t=0; t<num_tforms; ++t) {
        transforms[t] = moints_runtime4()->get(transform_keys[t]);
      }

      Timer tim_generic_tensor("Generic tensor");
      std::ostringstream oss;
      oss << "<" << space1_bra->id() << " " << space2_bra->id() << (antisymmetrize ? "||" : "|")
      << space1_ket->id() << " " << space2_ket->id() << ">";
      const std::string label = oss.str();
      ExEnv::out0() << std::endl << indent
      << "Entered generic tensor (" << label << ") evaluator" << std::endl;
      ExEnv::out0() << incindent;

      //
      // WARNING: Assuming all transforms are over same spaces!!!
      //
      Ref<OrbitalSpace> tspace1_bra = transforms[0]->space1();
      Ref<OrbitalSpace> tspace1_ket = transforms[0]->space2();
      Ref<OrbitalSpace> tspace2_bra = transforms[0]->space3();
      Ref<OrbitalSpace> tspace2_ket = transforms[0]->space4();

      // maps spaceX to spaceX of the transform
      std::vector<unsigned int> map1_bra, map1_ket, map2_bra, map2_ket;
      // maps space2_ket to space1_ket of transform
      std::vector<unsigned int> map12_ket;
      // maps space1_ket to space2_ket of transform
      std::vector<unsigned int> map21_ket;
      map1_bra = *tspace1_bra<<*space1_bra;
      map1_ket = *tspace1_ket<<*space1_ket;
      map2_bra = *tspace2_bra<<*space2_bra;
      map2_ket = *tspace2_ket<<*space2_ket;
      if (!alphabeta) {
        if (tspace1_ket == tspace2_ket) {
          map12_ket = map1_ket;
          map21_ket = map2_ket;
        }
        else {
          map12_ket = *tspace1_ket<<*space2_ket;
          map21_ket = *tspace2_ket<<*space1_ket;
        }
      }

      const unsigned int tblock_ncols = tspace2_ket->rank();
      const RefDiagSCMatrix evals1_bra = space1_bra->evals();
      const RefDiagSCMatrix evals1_ket = space1_ket->evals();
      const RefDiagSCMatrix evals2_bra = space2_bra->evals();
      const RefDiagSCMatrix evals2_ket = space2_ket->evals();

      // Using spinorbital iterators means I don't take into account perm symmetry
      // More efficient algorithm will require generic code
      const SpinCase2 S = (alphabeta ? AlphaBeta : AlphaAlpha);
      SpinMOPairIter iterbra(space1_bra,space2_bra,S);
      SpinMOPairIter iterket(space1_ket,space2_ket,S);
      // size of one block of <space1_bra space2_bra|
      const unsigned int nbra = iterbra.nij();
      // size of one block of |space1_ket space2_ket>
      const unsigned int nket = iterket.nij();

      RefSCMatrix Tresult;
      // Allocate storage for the result, if need to antisymmetrize at the end; else accumulate directly to T
      if (antisymmetrize && alphabeta) {
        Tresult = T.kit()->matrix(new SCDimension(nbra*nbrasets),
                                  new SCDimension(nket*nketsets));
        Tresult.assign(0.0);
      }
      else
        Tresult = T;

      unsigned int fbraket = 0;
      unsigned int fbraoffset = 0;
      for(unsigned int fbra=0; fbra<nbrasets; ++fbra,fbraoffset+=nbra) {

        unsigned int fketoffset = 0;
        for(unsigned int fket=0; fket<nketsets; ++fket,fketoffset+=nket,++fbraket) {

          Ref<TwoBodyMOIntsTransform> tform = transforms[fbraket];
          const Ref<TwoBodyIntDescr>& intdescr = tform->intdescr();
          const unsigned int intsetidx = intdescr->intset(tbint_type);

          if (debug_ > DefaultPrintThresholds::diagnostics)
            ExEnv::out0() << indent << "Using transform " << tform->name() << std::endl;

          tform->compute();
          Ref<DistArray4> accum = tform->ints_acc();
          accum->activate();

          // split work over tasks which have access to integrals
          std::vector<int> proc_with_ints;
          const int nproc_with_ints = accum->tasks_with_access(proc_with_ints);
          const int me = r12world()->world()->msg()->me();

          if (accum->has_access(me)) {
            for(iterbra.start(); iterbra; iterbra.next()) {
              const int ij = iterbra.ij();

              const int ij_proc = ij%nproc_with_ints;
              if (ij_proc != proc_with_ints[me])
                continue;

              const unsigned int i = iterbra.i();
              const unsigned int j = iterbra.j();
              const unsigned int ii = map1_bra[i];
              const unsigned int jj = map2_bra[j];

	      if (debug_ > DefaultPrintThresholds::allO2N2)
                ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << std::endl;
              Timer tim_intsretrieve("MO ints retrieve");
              const double *ij_buf = accum->retrieve_pair_block(ii,jj,intsetidx);
              tim_intsretrieve.exit();
	      if (debug_ > DefaultPrintThresholds::allO2N2)
                ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << std::endl;

              for(iterket.start(); iterket; iterket.next()) {
                const unsigned int a = iterket.i();
                const unsigned int b = iterket.j();
                const unsigned int aa = map1_ket[a];
                const unsigned int bb = map2_ket[b];
                const int AB = aa*tblock_ncols+bb;
                const int ab = iterket.ij();

                const double I_ijab = ij_buf[AB];

                if (alphabeta) {
		  if (debug_ > DefaultPrintThresholds::allO2N2) {
                    ExEnv::out0() << "i = " << i << " j = " << j << " a = " << a << " b = " << b
                    << " <ij|ab> = " << I_ijab << std::endl;
                  }
                  const double T_ijab = DataProcess::I2T(I_ijab,i,j,a,b,evals1_bra,evals1_ket,evals2_bra,evals2_ket);
                  Tresult.accumulate_element(ij+fbraoffset,ab+fketoffset,T_ijab);
                }
                else {
                  const int aa = map21_ket[a];
                  const int bb = map12_ket[b];
                  const int BA = bb*tblock_ncols+aa;
                  const double I_ijba = ij_buf[BA];
		  if (debug_ > DefaultPrintThresholds::allO2N2) {
                    ExEnv::out0() << "i = " << i << " j = " << j << " a = " << a << " b = " << b
                    << " <ij|ab> = " << I_ijab << " <ij|ba> = " << I_ijba << std::endl;
                  }
                  const double I_anti = I_ijab - I_ijba;
                  const double T_ijab = DataProcess::I2T(I_anti,i,j,a,b,evals1_bra,evals1_ket,evals2_bra,evals2_ket);
                  Tresult.accumulate_element(ij+fbraoffset,ab+fketoffset,T_ijab);
                }

              } // ket loop
              accum->release_pair_block(ii,jj,intsetidx);

            } // bra loop
          } // loop over tasks with access

          if (accum->data_persistent()) accum->deactivate();

        } // ket blocks
      } // bra blocks

      if (antisymmetrize && alphabeta) {
        // antisymmetrization implies equivalent particles -- hence symmetrize before antisymmetrize
        symmetrize<false>(Tresult,Tresult,space1_bra,space1_ket);
        sc::antisymmetrize(T,Tresult,space1_bra,space1_ket,true);
        Tresult = 0;
      }

      ExEnv::out0() << decindent;
      ExEnv::out0() << indent << "Exited generic tensor (" << label << ") evaluator" << std::endl;
      tim_generic_tensor.exit();
    }

  /// Contains classes used to compute many-body tensors
  namespace ManyBodyTensors {

    enum Sign {
      Minus = -1,
      Plus = +1
    };

    /// Tensor elements are <pq||rs>
    template <Sign sign = Plus>
    class Apply_Identity {
    public:
      static double I2T(double I, int i1, int i3, int i2, int i4,
      const RefDiagSCMatrix& evals1,
      const RefDiagSCMatrix& evals2,
      const RefDiagSCMatrix& evals3,
      const RefDiagSCMatrix& evals4)
      {
        if (sign == Plus)
          return I;
        else
          return -I;
      }
    };

    /// Applies (H0 - E0)
    template <Sign sign = Plus>
    class Apply_H0minusE0 {
    public:
      static double I2T(double I, int i1, int i3, int i2, int i4,
      const RefDiagSCMatrix& evals1,
      const RefDiagSCMatrix& evals2,
      const RefDiagSCMatrix& evals3,
      const RefDiagSCMatrix& evals4)
      {
        const double ediff = (- evals1(i1) - evals3(i3) + evals2(i2) + evals4(i4));
        if (sign == Plus)
          return I*ediff;
        else
          return -I*ediff;
      }
    };

    /// Applies (H0 - E0)^{-1}, e.g. MP2 T2 tensor elements are <ij||ab> /(e_i + e_j - e_a - e_b)
    template <Sign sign = Plus>
    class Apply_Inverse_H0minusE0 {
    public:
      static double I2T(double I, int i1, int i3, int i2, int i4,
      const RefDiagSCMatrix& evals1,
      const RefDiagSCMatrix& evals2,
      const RefDiagSCMatrix& evals3,
      const RefDiagSCMatrix& evals4)
      {
        const double ediff = (- evals1(i1) - evals3(i3) + evals2(i2) + evals4(i4));
        if (sign == Plus)
          return I/ediff;
        else
          return -I/ediff;
      }
    };

    /// Applies 1.0/sqrt(H0-E0)
    /// MP2 pseudo-T2 (S2) tensor elements are <ij||ab> /sqrt(|e_i + e_j - e_a - e_b|) such
    /// that MP2 pair energies are the diagonal elements of S2 * S2.t()
    template <Sign sign = Plus>
    class Apply_Inverse_Sqrt_H0minusE0 {
    public:
      static double I2T(double I, int i1, int i3, int i2, int i4,
                 const RefDiagSCMatrix& evals1,
                 const RefDiagSCMatrix& evals2,
                 const RefDiagSCMatrix& evals3,
                 const RefDiagSCMatrix& evals4)
      {
        const double ediff = (- evals1(i1) - evals3(i3) + evals2(i2) + evals4(i4));
        if (sign == Plus)
          return I/std::sqrt(std::fabs(ediff));
        else
          return -I/std::sqrt(std::fabs(ediff));
      }
    };

    typedef Apply_Identity<Plus> I_to_T;
    typedef Apply_Identity<Minus> I_to_mT;
    typedef Apply_Inverse_Sqrt_H0minusE0<Plus> ERI_to_S2;
    typedef Apply_Inverse_H0minusE0<Minus> ERI_to_T2;
  }

}

#endif

