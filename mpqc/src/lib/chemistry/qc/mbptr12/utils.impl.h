//
// utils.impl.h
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#ifndef _chemistry_qc_mbptr12_utilsimpl_h
#define _chemistry_qc_mbptr12_utilsimpl_h

#include <util/class/scexception.h>
#include <chemistry/qc/mbptr12/moindexspace.h>
#include <chemistry/qc/mbptr12/pairiter.h>

namespace sc {

  template <bool accumulate>
  void
  antisymmetrize(RefSCMatrix& Aanti, const RefSCMatrix& A,
                 const Ref<MOIndexSpace>& bra1,
                 const Ref<MOIndexSpace>& bra2,
                 const Ref<MOIndexSpace>& ket1,
                 const Ref<MOIndexSpace>& ket2)
  {
    const bool bra1_eq_bra2 = (bra1 == bra2);
    const bool ket1_eq_ket2 = (ket1 == ket2);
    if (!bra1_eq_bra2 && !ket1_eq_ket2)
      throw ProgrammingError("sc::antisymmetrize() -- the operation does not make sense: bra1!=bra2, ket1!=ket2",__FILE__,__LINE__);

    SpatialMOPairIter* ij_iter;
    SpatialMOPairIter* kl_iter;
    if (!bra1_eq_bra2)
      ij_iter = new SpatialMOPairIter_neq(bra1,bra2);
    else
      ij_iter = new SpatialMOPairIter_eq(bra1);
    if (!ket1_eq_ket2)
      kl_iter = new SpatialMOPairIter_neq(ket1,ket2);
    else
      kl_iter = new SpatialMOPairIter_eq(ket1);

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
  }

  template <bool Accumulate>
    void symmetrize(RefSCMatrix& Asymm, const RefSCMatrix& A,
                    const Ref<MOIndexSpace>& bra,
                    const Ref<MOIndexSpace>& ket)
    {
      if (A.rowdim().n() != Asymm.rowdim().n())
        throw ProgrammingError("sc::symmetrize() -- source and target matrices have different row dimensions",__FILE__,__LINE__);
      if (A.coldim().n() != Asymm.coldim().n())
        throw ProgrammingError("sc::symmetrize() -- source and target matrices have different column dimensions",__FILE__,__LINE__);
      SpatialMOPairIter_eq ij_iter(bra);
      SpatialMOPairIter_eq kl_iter(ket);
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
    
}

#endif

