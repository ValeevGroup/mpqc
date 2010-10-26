//
// primpairs.cc
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <util/misc/math.h>

#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/cints/primpairs.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}
inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

PrimPairsCints::PrimPairsCints(const Ref<GaussianBasisSet>& bs1,
				 const Ref<GaussianBasisSet>& bs2)
{
  bs1_ = bs1;
  bs2_ = bs2;
  nprim1_ = bs1_->nprimitive();
  nprim2_ = bs2_->nprimitive();
  prim_pair_ = new prim_pair_t[nprim1_*nprim2_];

  int nshell1 = bs1_->nshell();
  int nshell2 = bs2_->nshell();

  double A[3], B[3];
  for(int s1=0; s1<nshell1; s1++) {
    GaussianShell& shell1 = bs1_->shell(s1);
    int np1 = shell1.nprimitive();
    int p1_offset = bs1_->shell_to_primitive(s1);

    for(int xyz=0; xyz<3; xyz++)
      A[xyz] = bs1_->molecule()->r(bs1_->shell_to_center(s1),xyz);

    for(int s2=0; s2<nshell2; s2++) {
      GaussianShell& shell2 = bs2_->shell(s2);
      int np2 = shell2.nprimitive();
      int p2_offset = bs2_->shell_to_primitive(s2);

      for(int xyz=0; xyz<3; xyz++)
	B[xyz] = bs2_->molecule()->r(bs2_->shell_to_center(s2),xyz);

      double AB2 = (A[0]-B[0])*(A[0]-B[0]) + 
	(A[1]-B[1])*(A[1]-B[1]) +
	(A[2]-B[2])*(A[2]-B[2]);

      /*--- compute primitive pair data ---*/
      int p12_offset = nprim2_*p1_offset + p2_offset;
      prim_pair_t *row_offset = &(prim_pair_[p12_offset]);
      for(int p1=0; p1<np1; p1++, row_offset+=nprim2_) {
	prim_pair_t *pair_ptr = row_offset;
	for(int p2=0; p2<np2; p2++, pair_ptr++) {
	  double exp1, exp2, gamma, oog, t;
	  exp1 = shell1.exponent(p1);
	  exp2 = shell2.exponent(p2);
	  gamma = exp1 + exp2;
	  oog = 1.0 / gamma;
	  t = M_PI*oog;
	  
	  pair_ptr->gamma = gamma;
	  pair_ptr->ovlp = t*sqrt(t)*exp(-exp1*exp2*AB2*oog);
	  for(int xyz=0; xyz<3; xyz++)
	    pair_ptr->P[xyz] = (exp1*A[xyz] + exp2*B[xyz])*oog;
	}
      }
    }
  }
}

PrimPairsCints::~PrimPairsCints()
{
  delete[] prim_pair_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
