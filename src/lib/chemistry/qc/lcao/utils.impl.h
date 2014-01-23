//
// utils.impl.h
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

#include <util/misc/scexception.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <math/mmisc/pairiter.h>

#ifndef _chemistry_qc_lcao_utilsimpl_h
#define _chemistry_qc_lcao_utilsimpl_h

namespace sc {

  template <PureSpinCase2 spin>
  RefSCMatrix spinadapt(const RefSCMatrix &A,
                        const Ref<OrbitalSpace> &bra,
                        const Ref<OrbitalSpace> &ket){
    SpatialMOPairIter_eq ij_iter(bra->rank());
    SpatialMOPairIter_eq kl_iter(ket->rank());
    const unsigned int brablock_size_ab = ij_iter.nij_ab();
    const unsigned int ketblock_size_ab = kl_iter.nij_ab();
    if (A.rowdim().n()%brablock_size_ab)
      throw ProgrammingError("sc::spinadapt() -- row dimension is not integer multiple of bra-space rank",__FILE__,__LINE__);
    if (A.coldim().n()%ketblock_size_ab)
      throw ProgrammingError("sc::spinadapt() -- col dimension is not integer multiple of ket-space rank",__FILE__,__LINE__);
    const unsigned int nbra_blocks = A.rowdim().n() / brablock_size_ab;
    const unsigned int nket_blocks = A.coldim().n() / ketblock_size_ab;
    const unsigned int brablock_size_aa = ij_iter.nij_aa();
    const unsigned int ketblock_size_aa = kl_iter.nij_aa();
    RefSCMatrix Aspinadapted(A.rowdim(),A.coldim(),A->kit());
    const double one_divby_sqrt_two=1.0/sqrt(2.0);
    unsigned int bra_offset_ab = 0;
    unsigned int bra_offset_aa = 0;
    for(int brablock=0; brablock<nbra_blocks; brablock++, bra_offset_ab += brablock_size_ab, bra_offset_aa += brablock_size_aa) {
      for(ij_iter.start();int(ij_iter);ij_iter.next()) {

        const int ij_ab = ij_iter.ij_ab();

        unsigned int ket_offset_ab = 0;
        unsigned int ket_offset_aa = 0;
        for(int ketblock=0; ketblock<nket_blocks; ketblock++, ket_offset_ab += ketblock_size_ab, ket_offset_aa += ketblock_size_aa) {
          for(kl_iter.start();int(kl_iter);kl_iter.next()) {

            const int kl_ab = kl_iter.ij_ab();
            const int lk_ab = kl_iter.ij_ba();
            double Aspinadapted_element;

            if(spin==Singlet){
              double prefactor=1.0;
              if(ij_iter.i()==ij_iter.j()) {
                prefactor*=one_divby_sqrt_two;
              }
              if(kl_iter.i()==kl_iter.j()) {
                prefactor*=one_divby_sqrt_two;
              }

              Aspinadapted_element=prefactor*(A.get_element(ij_ab+bra_offset_ab,kl_ab+ket_offset_ab)+A.get_element(ij_ab+bra_offset_ab,lk_ab+ket_offset_ab));
            }
            else if(spin==Triplet){
              Aspinadapted_element=A.get_element(ij_ab+bra_offset_ab,kl_ab+ket_offset_ab)-A.get_element(ij_ab+bra_offset_ab,lk_ab+ket_offset_ab);
            }
            else {
              throw ProgrammingError("sc::spinadapt() -- PureSpinCase2 spin is not adeqate",__FILE__,__LINE__);
            }

            Aspinadapted.set_element(ij_ab+bra_offset_ab,kl_ab+ket_offset_ab,Aspinadapted_element);
          }
        }
      }
    }

    return(Aspinadapted);
  }

  template <bool accumulate>
  void
  antisymmetrize(RefSCMatrix& Aanti, const RefSCMatrix& A,
                 const Ref<OrbitalSpace>& bra1,
                 const Ref<OrbitalSpace>& bra2,
                 const Ref<OrbitalSpace>& ket1,
                 const Ref<OrbitalSpace>& ket2)
  {
    const bool bra1_eq_bra2 = (bra1 == bra2);
    const bool ket1_eq_ket2 = (ket1 == ket2);
    if (!bra1_eq_bra2 && !ket1_eq_ket2)
      throw ProgrammingError("sc::antisymmetrize() -- the operation does not make sense: bra1!=bra2, ket1!=ket2",__FILE__,__LINE__);

    SpatialMOPairIter* ij_iter;
    SpatialMOPairIter* kl_iter;
    if (!bra1_eq_bra2)
      ij_iter = new SpatialMOPairIter_neq(bra1->rank(),bra2->rank());
    else
      ij_iter = new SpatialMOPairIter_eq(bra1->rank());
    if (!ket1_eq_ket2)
      kl_iter = new SpatialMOPairIter_neq(ket1->rank(),ket2->rank());
    else
      kl_iter = new SpatialMOPairIter_eq(ket1->rank());

    const unsigned int brablock_size_ab = ij_iter->nij_ab();
    const unsigned int ketblock_size_ab = kl_iter->nij_ab();
    const unsigned int brablock_size_aa = ij_iter->nij_aa();
    const unsigned int ketblock_size_aa = kl_iter->nij_aa();
    if (brablock_size_ab==0 || ketblock_size_ab==0)
      return;
    if (A.rowdim().n()%brablock_size_ab)
      throw ProgrammingError("sc::antisymmetrize() -- row dimension of Source is not integer multiple of bra-space rank",__FILE__,__LINE__);
    if (A.coldim().n()%ketblock_size_ab)
      throw ProgrammingError("sc::antisymmetrize() -- col dimension of Source is not integer multiple of ket-space rank",__FILE__,__LINE__);
    if (Aanti.rowdim().n()%brablock_size_aa)
      throw ProgrammingError("sc::antisymmetrize() -- row dimension of Result is not integer multiple of bra-space rank",__FILE__,__LINE__);
    if (Aanti.coldim().n()%ketblock_size_aa)
      throw ProgrammingError("sc::antisymmetrize() -- col dimension of Result is not integer multiple of ket-space rank",__FILE__,__LINE__);
    const unsigned int nbra_blocks = A.rowdim().n() / brablock_size_ab;
    const unsigned int nket_blocks = A.coldim().n() / ketblock_size_ab;
    if (Aanti.rowdim().n() / brablock_size_aa != nbra_blocks)
      throw ProgrammingError("sc::antisymmetrize() -- bra dimensions of Source and Result do not match",__FILE__,__LINE__);
    if (Aanti.coldim().n() / ketblock_size_aa != nket_blocks)
      throw ProgrammingError("sc::antisymmetrize() -- ket dimensions of Source and Result do not match",__FILE__,__LINE__);

    unsigned int bra_offset_ab = 0;
    unsigned int bra_offset_aa = 0;
    for(int brablock=0; brablock<nbra_blocks; brablock++, bra_offset_ab += brablock_size_ab, bra_offset_aa += brablock_size_aa) {
      for(ij_iter->start();int(*ij_iter);ij_iter->next()) {

        const int ij_aa = ij_iter->ij_aa();
        if (ij_aa == -1)
          continue;
        const int ij_ab = ij_iter->ij_ab();
        const int ji_ab = ij_iter->ij_ba();

        unsigned int ket_offset_ab = 0;
        unsigned int ket_offset_aa = 0;
        for(int ketblock=0; ketblock<nket_blocks; ketblock++, ket_offset_ab += ketblock_size_ab, ket_offset_aa += ketblock_size_aa) {
          for(kl_iter->start();int(*kl_iter);kl_iter->next()) {

            const int kl_aa = kl_iter->ij_aa();
            if (kl_aa == -1)
              continue;
            const int kl_ab = kl_iter->ij_ab();
            const int lk_ab = kl_iter->ij_ba();

            double Aanti_ijkl = (ket1_eq_ket2) ?
                                                A.get_element(ij_ab+bra_offset_ab,kl_ab+ket_offset_ab) -
                                                A.get_element(ij_ab+bra_offset_ab,lk_ab+ket_offset_ab)
                                               :
                                               A.get_element(ij_ab+bra_offset_ab,kl_ab+ket_offset_ab) -
                                               A.get_element(ji_ab+bra_offset_ab,kl_ab+ket_offset_ab);
            if (accumulate)
              Aanti.accumulate_element(ij_aa+bra_offset_aa,kl_aa+ket_offset_aa,Aanti_ijkl);
            else
              Aanti.set_element(ij_aa+bra_offset_aa,kl_aa+ket_offset_aa,Aanti_ijkl);
          }
        }
      }
    }

    delete ij_iter;
    delete kl_iter;
  }

  template <bool accumulate>
  void
  antisymmetrize(RefSymmSCMatrix& Aanti, const RefSymmSCMatrix& A,
                 const Ref<OrbitalSpace>& bra1)
  {
    SpatialMOPairIter* ij_iter = new SpatialMOPairIter_eq(bra1->rank());
    SpatialMOPairIter* kl_iter = new SpatialMOPairIter_eq(bra1->rank());

    const unsigned int block_size_ab = ij_iter->nij_ab();
    const unsigned int block_size_aa = ij_iter->nij_aa();
    if (block_size_ab==0)
      return;
    if (A.dim().n()%block_size_ab)
      throw ProgrammingError("sc::antisymmetrize() -- dimension of Source is not integer multiple of (bra1->rank)^2",__FILE__,__LINE__);
    if (Aanti.dim().n()%block_size_aa)
      throw ProgrammingError("sc::antisymmetrize() -- dimension of Result is not integer multiple of (bra1->rank * (bra1->rank-1)/2)",__FILE__,__LINE__);
    const unsigned int nblocks = A.dim().n() / block_size_ab;
    if (Aanti.dim().n() / block_size_aa != nblocks)
      throw ProgrammingError("sc::antisymmetrize() -- dimensions of Source and Result do not match",__FILE__,__LINE__);

    unsigned int bra_offset_ab = 0;
    unsigned int bra_offset_aa = 0;
    for(unsigned int brablock=0; brablock<nblocks; brablock++, bra_offset_ab += block_size_ab, bra_offset_aa += block_size_aa) {
      for(ij_iter->start();int(*ij_iter);ij_iter->next()) {

        const int ij_aa = ij_iter->ij_aa();
        if (ij_aa == -1)
          continue;
        const int ij_ab = ij_iter->ij_ab();
        const int ji_ab = ij_iter->ij_ba();

        unsigned int ket_offset_ab = 0;
        unsigned int ket_offset_aa = 0;
        for(unsigned int ketblock=0; ketblock<nblocks; ketblock++, ket_offset_ab += block_size_ab, ket_offset_aa += block_size_aa) {
          for(kl_iter->start();int(*kl_iter);kl_iter->next()) {

            const int kl_aa = kl_iter->ij_aa();
            if (kl_aa == -1)
              continue;
            const int kl_ab = kl_iter->ij_ab();
            const int lk_ab = kl_iter->ij_ba();

            double Aanti_ijkl = A.get_element(ij_ab+bra_offset_ab,kl_ab+ket_offset_ab) -
                                A.get_element(ij_ab+bra_offset_ab,lk_ab+ket_offset_ab);
            if (accumulate)
              Aanti.accumulate_element(ij_aa+bra_offset_aa,kl_aa+ket_offset_aa,Aanti_ijkl);
            else
              Aanti.set_element(ij_aa+bra_offset_aa,kl_aa+ket_offset_aa,Aanti_ijkl);
          }
        }
      }
    }

    delete ij_iter;
    delete kl_iter;
  }

  template <bool accumulate>
  void antisymmetrize(double* Aanti, const double* A,
                      const int n) {
    double * result_ptr = Aanti;
    for(int r=0; r<n; ++r) {
      int RC_offset = r*n;
      int CR = r;
      for(int c=0; c<r; ++c, ++result_ptr, CR+=n) {
        const double v = A[RC_offset + c] - A[CR];
        if (accumulate)
          *result_ptr += v;
        else
          *result_ptr = v;
      }
    }
  }


  template <bool Accumulate>
    void symmetrize(RefSCMatrix& Asymm, const RefSCMatrix& A,
                    const Ref<OrbitalSpace>& bra,
                    const Ref<OrbitalSpace>& ket)
    {
      if (A.rowdim().n() != Asymm.rowdim().n())
        throw ProgrammingError("sc::symmetrize() -- source and target matrices have different row dimensions",__FILE__,__LINE__);
      if (A.coldim().n() != Asymm.coldim().n())
        throw ProgrammingError("sc::symmetrize() -- source and target matrices have different column dimensions",__FILE__,__LINE__);
      SpatialMOPairIter_eq ij_iter(bra->rank());
      SpatialMOPairIter_eq kl_iter(ket->rank());
      const unsigned int brablock_size_ab = ij_iter.nij_ab();
      const unsigned int ketblock_size_ab = kl_iter.nij_ab();
      if (A.rowdim().n()%brablock_size_ab)
        throw ProgrammingError("sc::symmetrize() -- row dimension is not integer multiple of bra-space rank",__FILE__,__LINE__);
      if (A.coldim().n()%ketblock_size_ab)
        throw ProgrammingError("sc::symmetrize() -- col dimension is not integer multiple of ket-space rank",__FILE__,__LINE__);
      const unsigned int nbra_blocks = A.rowdim().n() / brablock_size_ab;
      const unsigned int nket_blocks = A.coldim().n() / ketblock_size_ab;

      unsigned int bra_offset_ab = 0;
      for(int brablock=0; brablock<nbra_blocks; brablock++, bra_offset_ab += brablock_size_ab) {
        for(ij_iter.start();int(ij_iter);ij_iter.next()) {

          const unsigned int IJ_ab = ij_iter.ij_ab() + bra_offset_ab;
          const unsigned int JI_ab = ij_iter.ij_ba() + bra_offset_ab;

          unsigned int ket_offset_ab = 0;
          for(int ketblock=0; ketblock<nket_blocks; ketblock++, ket_offset_ab += ketblock_size_ab) {
            for(kl_iter.start();int(kl_iter);kl_iter.next()) {

              const unsigned int KL_ab = kl_iter.ij_ab() + ket_offset_ab;
              const unsigned int LK_ab = kl_iter.ij_ba() + ket_offset_ab;

              const double A_ijkl = A.get_element(IJ_ab,KL_ab);
              const double A_ijlk = A.get_element(IJ_ab,LK_ab);
              const double A_jikl = A.get_element(JI_ab,KL_ab);
              const double A_jilk = A.get_element(JI_ab,LK_ab);

              {
                const double Asymm_ijkl = 0.5 * (A_ijkl + A_jilk);
                if (Accumulate) {
                  Asymm.accumulate_element(IJ_ab,KL_ab,Asymm_ijkl);
                  Asymm.accumulate_element(JI_ab,LK_ab,Asymm_ijkl);
                }
                else {
                  Asymm.set_element(IJ_ab,KL_ab,Asymm_ijkl);
                  Asymm.set_element(JI_ab,LK_ab,Asymm_ijkl);
                }
              }
              {
                const double Asymm_ijlk = 0.5 * (A_ijlk + A_jikl);
                if (Accumulate) {
                  Asymm.accumulate_element(IJ_ab,LK_ab,Asymm_ijlk);
                  Asymm.accumulate_element(JI_ab,KL_ab,Asymm_ijlk);
                }
                else {
                  Asymm.set_element(IJ_ab,LK_ab,Asymm_ijlk);
                  Asymm.set_element(JI_ab,KL_ab,Asymm_ijlk);
                }
              }

            } // end of kl
          } // endo f ket blocks
        } // end of ij
      } // end of bra blocks

    }

  template <bool Accumulate, sc::fastpairiter::PairSymm BraSymm, sc::fastpairiter::PairSymm KetSymm>
    void symmetrize12(RefSCMatrix& Asymm, const RefSCMatrix& A,
                      const Ref<OrbitalSpace>& bra1,
                      const Ref<OrbitalSpace>& bra2,
                      const Ref<OrbitalSpace>& ket1,
                      const Ref<OrbitalSpace>& ket2)
    {
      // doesn't make sense if Asymm == A
      MPQC_ASSERT(&Asymm != &A);
      using namespace sc::fastpairiter;
      using sc::fastpairiter::MOPairIter;

      // Detect inputs that violate semantics of this function
      if ( (BraSymm == AntiSymm && KetSymm == AntiSymm) ||
           (BraSymm == Symm && KetSymm == Symm) )
        throw ProgrammingError("sc::symmetrize12() -- the matrix is already symmetric",__FILE__,__LINE__);
      if ( (BraSymm == AntiSymm || BraSymm == Symm) && bra1 != bra2)
        throw ProgrammingError("sc::symmetrize12() -- bra dimension is anti/symmetrized, but bra1!=bra2",__FILE__,__LINE__);
      if ( (KetSymm == AntiSymm || KetSymm == Symm) && ket1 != ket2)
        throw ProgrammingError("sc::symmetrize12() -- ket dimension is anti/symmetrized, but ket1!=ket2",__FILE__,__LINE__);

      if (A.rowdim().n() != Asymm.rowdim().n())
        throw ProgrammingError("sc::symmetrize() -- source and target matrices have different row dimensions",__FILE__,__LINE__);
      if (A.coldim().n() != Asymm.coldim().n())
        throw ProgrammingError("sc::symmetrize() -- source and target matrices have different column dimensions",__FILE__,__LINE__);
      MOPairIter<BraSymm> ij_iter(bra1->rank(),bra2->rank());
      MOPairIter<KetSymm> kl_iter(ket1->rank(),ket2->rank());
      const unsigned int brablock_size = ij_iter.nij();
      const unsigned int ketblock_size = kl_iter.nij();
      if (A.rowdim().n()%brablock_size)
        throw ProgrammingError("sc::symmetrize() -- row dimension is not integer multiple of bra-space rank",__FILE__,__LINE__);
      if (A.coldim().n()%ketblock_size)
        throw ProgrammingError("sc::symmetrize() -- col dimension is not integer multiple of ket-space rank",__FILE__,__LINE__);
      const unsigned int nbra_blocks = A.rowdim().n() / brablock_size;
      const unsigned int nket_blocks = A.coldim().n() / ketblock_size;

      unsigned int bra_offset = 0;
      for(int brablock=0; brablock<nbra_blocks; brablock++, bra_offset += brablock_size) {
        for(ij_iter.start();int(ij_iter);ij_iter.next()) {

          const unsigned int IJ = ij_iter.ij() + bra_offset;
          const unsigned int JI = (BraSymm == ASymm) ? ij_iter.ij(ij_iter.j(),ij_iter.i()) + bra_offset
                                                     : IJ;
          const double ij_permfac = (BraSymm == AntiSymm) ? -1.0 : 1.0;

          unsigned int ket_offset = 0;
          for(int ketblock=0; ketblock<nket_blocks; ketblock++, ket_offset += ketblock_size) {
            for(kl_iter.start();int(kl_iter);kl_iter.next()) {

              const unsigned int KL = kl_iter.ij() + ket_offset;
              const unsigned int LK = (KetSymm == ASymm) ? kl_iter.ij(kl_iter.j(),kl_iter.i()) + ket_offset
                                                         : KL;
              const double kl_permfac = (KetSymm == AntiSymm) ? -1.0 : 1.0;

              const double A_ijkl = A.get_element(IJ,KL);
              const double A_jilk = A.get_element(JI,LK);
              const double Asymm_ijkl = 0.5 * (A_ijkl + ij_permfac*kl_permfac*A_jilk);

              if (Accumulate)
                Asymm.accumulate_element(IJ,KL,Asymm_ijkl);
              else
                Asymm.set_element(IJ,KL,Asymm_ijkl);

              //ExEnv::out0() << IJ << " " << KL << " " << JI << " " << LK << " " << A_ijkl << " " << A_jilk << " " << Asymm_ijkl << std::endl;

            } // end of kl
          } // end of ket blocks
        } // end of ij
      } // end of bra blocks

    }

  template <bool Accumulate,
    sc::fastpairiter::PairSymm SrcBraSymm, sc::fastpairiter::PairSymm SrcKetSymm,
    sc::fastpairiter::PairSymm DstBraSymm, sc::fastpairiter::PairSymm DstKetSymm
    >
    void symmetrize(RefSCMatrix& Asymm, const RefSCMatrix& A,
                      const Ref<OrbitalSpace>& bra1,
                      const Ref<OrbitalSpace>& bra2,
                      const Ref<OrbitalSpace>& ket1,
                      const Ref<OrbitalSpace>& ket2)
    {
      using namespace sc::fastpairiter;
      using sc::fastpairiter::MOPairIter;

      // Detect inputs that violate semantics of this function
      if (SrcBraSymm == DstBraSymm && SrcKetSymm == DstKetSymm)
        throw ProgrammingError("sc::symmetrize() -- nothing to be done",__FILE__,__LINE__);
      if ( (SrcBraSymm != ASymm && SrcKetSymm != ASymm) )
        throw ProgrammingError("sc::symmetrize() -- either bra of ket must be asymmetric",__FILE__,__LINE__);
      if ( (SrcBraSymm != ASymm && SrcBraSymm != DstBraSymm) )
        throw ProgrammingError("sc::symmetrize() -- can only change bra that is asymmetric",__FILE__,__LINE__);
      if ( (SrcKetSymm != ASymm && SrcKetSymm != DstKetSymm) )
        throw ProgrammingError("sc::symmetrize() -- can only change ket that is asymmetric",__FILE__,__LINE__);
      if (DstBraSymm!=ASymm && bra1 != bra2)
        throw ProgrammingError("sc::symmetrize12() -- bra is to be anti/symmetrized, but bra1!=bra2",__FILE__,__LINE__);
      if (DstKetSymm!=ASymm && ket1 != ket2)
        throw ProgrammingError("sc::symmetrize12() -- ket is to be anti/symmetrized, but ket1!=ket2",__FILE__,__LINE__);
      const bool transform_bra = DstBraSymm != ASymm && DstBraSymm != SrcBraSymm;
      const bool transform_ket = DstKetSymm != ASymm && DstKetSymm != SrcKetSymm;

      MOPairIter<SrcBraSymm> ij_srciter(bra1->rank(),bra2->rank());
      MOPairIter<SrcKetSymm> kl_srciter(ket1->rank(),ket2->rank());
      MOPairIter<DstBraSymm> ij_dstiter(bra1->rank(),bra2->rank());
      MOPairIter<DstKetSymm> kl_dstiter(ket1->rank(),ket2->rank());
      const unsigned int brablock_size_src = ij_srciter.nij();
      const unsigned int ketblock_size_src = kl_srciter.nij();
      const unsigned int brablock_size_dst = ij_dstiter.nij();
      const unsigned int ketblock_size_dst = kl_dstiter.nij();
      if (A.rowdim().n()%brablock_size_src)
        throw ProgrammingError("sc::symmetrize() -- row dimension of Source is not integer multiple of bra-space rank",__FILE__,__LINE__);
      if (A.coldim().n()%ketblock_size_src)
        throw ProgrammingError("sc::symmetrize() -- col dimension of Source is not integer multiple of ket-space rank",__FILE__,__LINE__);
      if (Asymm.rowdim().n()%brablock_size_dst)
        throw ProgrammingError("sc::symmetrize() -- row dimension of Result is not integer multiple of bra-space rank",__FILE__,__LINE__);
      if (Asymm.coldim().n()%ketblock_size_dst)
        throw ProgrammingError("sc::symmetrize() -- col dimension of Result is not integer multiple of ket-space rank",__FILE__,__LINE__);
      const unsigned int nbra_blocks = A.rowdim().n() / brablock_size_src;
      const unsigned int nket_blocks = A.coldim().n() / ketblock_size_src;
      if (nbra_blocks != Asymm.rowdim().n() / brablock_size_dst)
        throw ProgrammingError("sc::symmetrize() -- # of bra blocks in Source and Result do not match",__FILE__,__LINE__);
      if (nket_blocks != Asymm.coldim().n() / ketblock_size_dst)
        throw ProgrammingError("sc::symmetrize() -- # of bra blocks in Source and Result do not match",__FILE__,__LINE__);

      unsigned int bra_offset_src = 0;
      unsigned int bra_offset_dst = 0;
      for(int brablock=0; brablock<nbra_blocks;
          brablock++, bra_offset_src += brablock_size_src, bra_offset_dst += brablock_size_dst) {
        for(ij_dstiter.start();int(ij_dstiter);ij_dstiter.next()) {

          const unsigned int IJ_dst = ij_dstiter.ij() + bra_offset_dst;

          const unsigned int I = ij_srciter.i();
          const unsigned int J = ij_srciter.j();
          const unsigned IJ = ij_srciter.ij(I,J) + bra_offset_src;
          unsigned int JI;
          double ij_permfac;
          if (transform_bra) {
            JI = ij_srciter.ij(J,I) + bra_offset_src;
            ij_permfac = (DstBraSymm == AntiSymm) ? -1.0 : 1.0;
          }

          unsigned int ket_offset_src = 0;
          unsigned int ket_offset_dst = 0;
          for(int ketblock=0; ketblock<nket_blocks;
              ketblock++, ket_offset_src += ketblock_size_src, ket_offset_dst += ketblock_size_dst) {
            for(kl_dstiter.start();int(kl_dstiter);kl_dstiter.next()) {

              const unsigned int KL_dst = kl_dstiter.ij() + ket_offset_dst;

              const unsigned int K = kl_srciter.i();
              const unsigned int L = kl_srciter.j();
              const unsigned int KL = kl_srciter.ij(K,L) + ket_offset_src;
              unsigned int LK;
              double kl_permfac;
              if (transform_ket) {
                LK = kl_srciter.ij(L,K) + ket_offset_src;
                kl_permfac = (DstKetSymm == AntiSymm) ? -1.0 : 1.0;
              }

              double Asymm_ijkl;
              if (transform_bra) {
                const double A_ijkl = A.get_element(IJ,KL);
                const double A_jikl = A.get_element(JI,KL);
                Asymm_ijkl = A_ijkl + ij_permfac*A_jikl;
              }
              if (transform_ket) {
                const double A_ijkl = A.get_element(IJ,KL);
                const double A_ijlk = A.get_element(IJ,LK);
                Asymm_ijkl = A_ijkl + kl_permfac*A_ijlk;
              }
              if (transform_bra && transform_ket) {
                const double A_ijkl = A.get_element(IJ,KL);
                const double A_jikl = A.get_element(JI,KL);
                const double A_ijlk = A.get_element(IJ,LK);
                const double A_jilk = A.get_element(JI,LK);
                Asymm_ijkl = A_ijkl + ij_permfac*A_jikl + kl_permfac*A_ijlk + ij_permfac*kl_permfac*A_jilk;
              }

              if (Accumulate)
                Asymm.accumulate_element(IJ_dst,KL_dst,Asymm_ijkl);
              else
                Asymm.set_element(IJ_dst,KL_dst,Asymm_ijkl);

            } // end of kl
          } // end of ket blocks
        } // end of ij
      } // end of bra blocks

    }

}

#endif

