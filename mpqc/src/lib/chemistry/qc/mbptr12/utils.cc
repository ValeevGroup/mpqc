//
// utils.cc
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

#include <util/misc/scexception.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/pairiter.h>

void
sc::antisymmetrize(RefSCMatrix& Aanti, const RefSCMatrix& A,
                   const Ref<MOIndexSpace>& bra,
                   const Ref<MOIndexSpace>& ket)
{
  SpatialMOPairIter_eq ij_iter(bra);
  SpatialMOPairIter_eq kl_iter(ket);
  const unsigned int brablock_size_ab = ij_iter.nij_ab();
  const unsigned int ketblock_size_ab = kl_iter.nij_ab();
  if (A.rowdim().n()%brablock_size_ab)
    throw ProgrammingError("sc::antisymmetrize() -- row dimension is not integer multiple of bra-space rank",__FILE__,__LINE__);
  if (A.coldim().n()%ketblock_size_ab)
    throw ProgrammingError("sc::antisymmetrize() -- col dimension is not integer multiple of ket-space rank",__FILE__,__LINE__);
  const unsigned int nbra_blocks = A.rowdim().n() / brablock_size_ab;
  const unsigned int nket_blocks = A.coldim().n() / ketblock_size_ab;
  const unsigned int brablock_size_aa = ij_iter.nij_aa();
  const unsigned int ketblock_size_aa = kl_iter.nij_aa();
  
  unsigned int bra_offset_ab = 0;
  unsigned int bra_offset_aa = 0;
  for(int brablock=0; brablock<nbra_blocks; brablock++, bra_offset_ab += brablock_size_ab, bra_offset_aa += brablock_size_aa) {
    for(ij_iter.start();int(ij_iter);ij_iter.next()) {
      
      const int ij_aa = ij_iter.ij_aa();
      if (ij_aa == -1)
        continue;
      const int ij_ab = ij_iter.ij_ab();
      
      unsigned int ket_offset_ab = 0;
      unsigned int ket_offset_aa = 0;
      for(int ketblock=0; ketblock<nket_blocks; ketblock++, ket_offset_ab += ketblock_size_ab, ket_offset_aa += ketblock_size_aa) {
        for(kl_iter.start();int(kl_iter);kl_iter.next()) {
          
          const int kl_aa = kl_iter.ij_aa();
          if (kl_aa == -1)
            continue;
          const int kl_ab = kl_iter.ij_ab();
          const int lk_ab = kl_iter.ij_ba();
          
          double Aanti_ijkl = A.get_element(ij_ab+bra_offset_ab,kl_ab+ket_offset_ab) - A.get_element(ij_ab+bra_offset_ab,lk_ab+ket_offset_ab);
          Aanti.set_element(ij_aa+bra_offset_aa,kl_aa+ket_offset_aa,Aanti_ijkl);
        }
      }
    }
  }
}

