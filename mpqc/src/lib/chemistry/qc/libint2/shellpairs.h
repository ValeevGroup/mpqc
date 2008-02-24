//
// shellpairs.h
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

#ifndef _chemistry_qc_libint2_shellpairs_h
#define _chemistry_qc_libint2_shellpairs_h

#include <util/ref/ref.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/libint2/primpairs.h>

namespace sc {

class ShellPairsLibint2;

/** ShellPairLibint2 is an interface to PrimPairsLibint2. It provides all primitive pair data for a given shell pair */
class ShellPairLibint2 {
  const PrimPairsLibint2& prim_pairs_;
  int prim1_offset_;
  int prim2_offset_;

  public:
  ShellPairLibint2(const PrimPairsLibint2&);
  ~ShellPairLibint2();

  void init(const int,
	    const int);

  prim_pair_t* prim_pair(int p1, int p2) const { return prim_pairs_.prim_pair(p1+prim1_offset_,p2+prim2_offset_); };
  double prim_pair_P(int p1, int p2, int xyz) const { return prim_pairs_.P(p1+prim1_offset_,p2+prim2_offset_,xyz); };
  double prim_pair_gamma(int p1, int p2) const { return prim_pairs_.gamma(p1+prim1_offset_,p2+prim2_offset_); };
  double prim_pair_ovlp(int p1, int p2) const { return prim_pairs_.ovlp(p1+prim1_offset_,p2+prim2_offset_); }
};


/** ShellPairsLibint2 contains primitive pair data for all shell pairs */
class ShellPairsLibint2: virtual public SavableState {
  Ref<GaussianBasisSet> bs1_;
  Ref<GaussianBasisSet> bs2_;
  Ref<PrimPairsLibint2> prim_pairs_;
  ShellPairLibint2* shell_pair_;
  
  public:
  ShellPairsLibint2(const Ref<GaussianBasisSet>&,
		    const Ref<GaussianBasisSet>&);
  ShellPairsLibint2(const Ref<ShellPairsLibint2>&);
  ShellPairsLibint2(const Ref<KeyVal>&);
  ShellPairsLibint2(StateIn&);

  ~ShellPairsLibint2();

  void save_data_state(StateOut&);

  ShellPairLibint2* shell_pair(int si, int sj) const {
    shell_pair_->init(bs1_->shell_to_primitive(si), bs2_->shell_to_primitive(sj));
    return shell_pair_;
  }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
