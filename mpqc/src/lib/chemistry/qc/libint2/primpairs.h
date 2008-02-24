//
// primpairs.h
//
// Copyright (C) 2001 Edward Valeev
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
#pragma interface
#endif

#ifndef _chemistry_qc_libint2_primpairs_h
#define _chemistry_qc_libint2_primpairs_h

#include <util/ref/ref.h>
#include <chemistry/qc/basis/basis.h>

namespace sc {

typedef struct {
      double P[3];
      double gamma;
      double ovlp;
} prim_pair_t;

/** PrimPairsLibint2 contains primitive pair data */
class PrimPairsLibint2 : public RefCount {
  Ref<GaussianBasisSet> bs1_;
  Ref<GaussianBasisSet> bs2_;
  int nprim1_;
  int nprim2_;
  prim_pair_t *prim_pair_;
  
  public:
  PrimPairsLibint2(const Ref<GaussianBasisSet>&,
		 const Ref<GaussianBasisSet>&);
  ~PrimPairsLibint2();

  prim_pair_t* prim_pair(int p1, int p2) const { return prim_pair_ + p1*nprim2_ + p2; };
  double P(int p1, int p2, int xyz) const { return prim_pair_[p1*nprim2_ + p2].P[xyz]; };
  double gamma(int p1, int p2) const { return prim_pair_[p1*nprim2_ + p2].gamma; };
  double ovlp(int p1, int p2) const { return prim_pair_[p1*nprim2_ + p2].ovlp; };

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
