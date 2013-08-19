//
// primpairs.h
//
// Copyright (C) 2001 Edward Valeev
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

/** PrimPairsLibint2 contains primitive pair data. The data can be viewed as a matrix,
 * with row and column dimensions referring to the primitives of basis set 1 and 2, respectively.
 * The ordering of primitives is in the order of appearance of shells in the basis sets. */
class PrimPairsLibint2 : public RefCount {
  friend class ShellPairLibint2;
  Ref<GaussianBasisSet> bs1_;
  Ref<GaussianBasisSet> bs2_;
  unsigned int nprim1_;
  unsigned int nprim2_;
  prim_pair_t *prim_pair_;
  
  std::vector<unsigned int> shell_to_prim1_;
  std::vector<unsigned int> shell_to_prim2_;

  public:
  PrimPairsLibint2(const Ref<GaussianBasisSet>&,
                   const Ref<GaussianBasisSet>&);
  ~PrimPairsLibint2();

  prim_pair_t* prim_pair(unsigned int p1, unsigned int p2) const { return prim_pair_ + p1*nprim2_ + p2; };
  double P(unsigned int p1, unsigned int p2, unsigned int xyz) const { return prim_pair_[p1*nprim2_ + p2].P[xyz]; };
  double gamma(unsigned int p1, unsigned int p2) const { return prim_pair_[p1*nprim2_ + p2].gamma; };
  double ovlp(unsigned int p1, unsigned int p2) const { return prim_pair_[p1*nprim2_ + p2].ovlp; };

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
