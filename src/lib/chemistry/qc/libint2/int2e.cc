//
// int2e.cc
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

#include <util/misc/formio.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/libint2/int2e.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}
inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

Int2eLibint2::Int2eLibint2(Integral *integral,
                 const Ref<GaussianBasisSet>& b1,
                 const Ref<GaussianBasisSet>& b2,
                 const Ref<GaussianBasisSet>& b3,
		 const Ref<GaussianBasisSet>& b4,
		 size_t storage) :
  integral_(integral),
  grp_(integral->messagegrp()),
  permute_(0),
  redundant_(1)
{
  bs1_ = b1;
  bs2_ = b2;
  bs3_ = b3;
  bs4_ = b4;

  if (bs2_ == 0) bs2_ = bs1_;
  if (bs3_ == 0) bs3_ = bs2_;
  if (bs4_ == 0) bs4_ = bs3_;

  /*--- Initialize storage ---*/
  init_storage(storage);

  /*--- allocate scratch for transformation ---*/
  if (bs1_->has_pure() || bs2_->has_pure() || bs3_->has_pure() || bs4_->has_pure() ||
      bs1_->max_ncontraction() != 1 || bs2_->max_ncontraction() != 1 ||
      bs3_->max_ncontraction() != 1 || bs4_->max_ncontraction() != 1) {
    // compute how much space one contraction quartet may need
    int nshell1 = bs1_->nshell();
    int maxncart1 = 0;
    for(int sh1=0; sh1<nshell1;sh1++) {
      int maxncart = bs1_->shell(sh1).max_cartesian();
      if (maxncart > maxncart1) maxncart1 = maxncart;
    }
    int nshell2 = bs2_->nshell();
    int maxncart2 = 0;
    for(int sh2=0; sh2<nshell2;sh2++) {
      int maxncart = bs2_->shell(sh2).max_cartesian();
      if (maxncart > maxncart2) maxncart2 = maxncart;
    }
    int nshell3 = bs3_->nshell();
    int maxncart3 = 0;
    for(int sh3=0; sh3<nshell3;sh3++) {
      int maxncart = bs3_->shell(sh3).max_cartesian();
      if (maxncart > maxncart3) maxncart3 = maxncart;
    }
    int nshell4 = bs4_->nshell();
    int maxncart4 = 0;
    for(int sh4=0; sh4<nshell4;sh4++) {
      int maxncart = bs4_->shell(sh4).max_cartesian();
      if (maxncart > maxncart4) maxncart4 = maxncart;
    }
    tformbuf_ = new double[maxncart1*maxncart2*maxncart3*maxncart4];
  }
  else {
    tformbuf_ = 0;
  }
}


Int2eLibint2::~Int2eLibint2()
{ 
  if (tformbuf_)
    delete[] tformbuf_;
  done_storage();
}

void
Int2eLibint2::bounds(const Ref<Log2Bounds>& b)
{
    bounds_ = b;
}

size_t
Int2eLibint2::storage_required_(const Ref<GaussianBasisSet>& b1,
				const Ref<GaussianBasisSet>& b2,
				const Ref<GaussianBasisSet>& b3,
				const Ref<GaussianBasisSet>& b4)
{
  size_t storage_required = 0;
  
  Ref<GaussianBasisSet> bs1 = b1;
  Ref<GaussianBasisSet> bs2 = b2;
  Ref<GaussianBasisSet> bs3 = b3;
  Ref<GaussianBasisSet> bs4 = b4;

  if (bs2 == 0)
    bs2 = bs1;
  if (bs3 == 0)
    bs3 = bs1;
  if (bs4 == 0)
    bs4 = bs1;

  if (bs1->has_pure() || bs2->has_pure() || bs3->has_pure() || bs4->has_pure() ||
      bs1->max_ncontraction() != 1 || bs2->max_ncontraction() != 1 ||
      bs3->max_ncontraction() != 1 || bs4->max_ncontraction() != 1) {
    // compute how much space one contraction quartet may need
    int nshell1 = bs1->nshell();
    int maxncart1 = 0;
    for(int sh1=0; sh1<nshell1;sh1++) {
      int maxncart = bs1->shell(sh1).max_cartesian();
      if (maxncart > maxncart1) maxncart1 = maxncart;
    }
    int nshell2 = bs2->nshell();
    int maxncart2 = 0;
    for(int sh2=0; sh2<nshell2;sh2++) {
      int maxncart = bs2->shell(sh2).max_cartesian();
      if (maxncart > maxncart2) maxncart2 = maxncart;
    }
    int nshell3 = bs3->nshell();
    int maxncart3 = 0;
    for(int sh3=0; sh3<nshell3;sh3++) {
      int maxncart = bs3->shell(sh3).max_cartesian();
      if (maxncart > maxncart3) maxncart3 = maxncart;
    }
    int nshell4 = bs4->nshell();
    int maxncart4 = 0;
    for(int sh4=0; sh4<nshell4;sh4++) {
      int maxncart = bs4->shell(sh4).max_cartesian();
      if (maxncart > maxncart4) maxncart4 = maxncart;
    }
    storage_required = maxncart1*maxncart2*maxncart3*maxncart4*sizeof(double);
  }
  else {
    storage_required = 0;
  }

  return storage_required;
}

int
Int2eLibint2::log2_bound(int s1, int s2, int s3, int s4)
{
    if (bounds_)
	return bounds_->log2_bound(s1,s2,s3,s4);
    else
	// 2^256 ~ 10^26
	return 256;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
