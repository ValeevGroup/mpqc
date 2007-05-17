//
// tbint.cc
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

#ifdef __GNUC__
#pragma implementation "junk.h"
#endif

#include <limits.h>
#include <util/class/scexception.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/intdescr.h>
#include <chemistry/qc/cints/tbintcints.h>
#include <chemistry/qc/cints/eri.h>
#include <chemistry/qc/cints/grt.h>

using namespace std;
using namespace sc;

inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

TwoBodyIntCints::TwoBodyIntCints(Integral*integral,
				 const Ref<GaussianBasisSet>& b1,
				 const Ref<GaussianBasisSet>& b2,
				 const Ref<GaussianBasisSet>& b3,
				 const Ref<GaussianBasisSet>& b4,
				 size_t storage, tbinteval int2etype):
    TwoBodyInt(integral,b1,b2,b3,b4), int2etype_(int2etype)
{
  // Which evaluator to use
  switch (int2etype_) {
  case erieval:
    int2ecints_ = new EriCints(integral,b1,b2,b3,b4,storage);
    num_tbint_types_ = 1;
    break;
  case grteval:
    int2ecints_ = new GRTCints(integral,b1,b2,b3,b4,storage);
    num_tbint_types_ = 4;
    break;
  default:
    ExEnv::errn() << scprintf("Tried to construct a two-electron integral evaluator of unimplemented or unknown type") << endl;
    fail();
  }

  buffer_ = int2ecints_->buffer(TwoBodyInt::eri);
  integral_->adjust_storage(int2ecints_->storage_used());
}

TwoBodyIntCints::~TwoBodyIntCints()
{
  integral_->adjust_storage(-int2ecints_->storage_used());
}

void
TwoBodyIntCints::compute_shell(int is, int js, int ks, int ls)
{
  int2ecints_->set_redundant(redundant());
  int2ecints_->compute_quartet(&is,&js,&ks,&ls);
}

int
TwoBodyIntCints::log2_shell_bound(int is, int js, int ks, int ls)
{
  return SCHAR_MAX;//int2ecints_->erep_4bound(is,js,ks,ls);
}

void
TwoBodyIntCints::set_integral_storage(size_t storage)
{
  int2ecints_->init_storage(storage);
}

unsigned int
TwoBodyIntCints::inttype(TwoBodyInt::tbint_type type) const
{
    switch(int2etype_) {
    case erieval:
	return TwoBodyIntDescrERI::intSet(type);
    case grteval:
	return TwoBodyIntDescrR12::intSet(type);
    }
}

TwoBodyInt::tbint_type
TwoBodyIntCints::inttype(unsigned int type) const
{
    switch(int2etype_) {
    case erieval:
	return TwoBodyIntDescrERI::intSet(type);
    case grteval:
	return TwoBodyIntDescrR12::intSet(type);
    }
}

//////////////////////////////////////////////////////////////////////////

TwoBodyDerivIntCints::TwoBodyDerivIntCints(Integral*integral,
					   const Ref<GaussianBasisSet>& b1,
					   const Ref<GaussianBasisSet>& b2,
					   const Ref<GaussianBasisSet>& b3,
					   const Ref<GaussianBasisSet>& b4,
					   size_t storage, tbinteval int2etype):
  TwoBodyDerivInt(integral,b1,b2,b3,b4)
{
  // Which evaluator to use
  switch (int2etype) {
  default:
    ExEnv::errn() << scprintf("Tried to construct a two-electron derivative integral evaluator of unimplemented or unknown type") << endl;
    fail();
  }

  //  int2ecints_ = new EriCints(integral,b1,b2,b3,b4,1,storage);
  //buffer_ = int2ecints_->buffer();
  //integral_->adjust_storage(int2ecints_->storage_used());
}

TwoBodyDerivIntCints::~TwoBodyDerivIntCints()
{
  //  integral_->adjust_storage(-int2ecints_->used_storage());
}

void
TwoBodyDerivIntCints::compute_shell(int is, int js, int ks, int ls,
                                 DerivCenters&c)
{
  int center;
  int sh[4], sz[4];

  sh[0]=is; sh[1]=js; sh[2]=ks; sh[3]=ls;

}


int
TwoBodyDerivIntCints::log2_shell_bound(int is, int js, int ks, int ls)
{
  return SCHAR_MAX;//int2ecints_->erep_4bound_1der(is,js,ks,ls);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
