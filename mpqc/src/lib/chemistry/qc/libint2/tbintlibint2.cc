//
// tbintlibint2.cc
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

#include <util/misc/scexception.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/libint2/tbintlibint2.h>
#include <chemistry/qc/libint2/eri.h>
#include <chemistry/qc/libint2/g12.h>

using namespace std;
using namespace sc;

inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

TwoBodyIntLibint2::TwoBodyIntLibint2(Integral*integral,
				 const Ref<GaussianBasisSet>& b1,
				 const Ref<GaussianBasisSet>& b2,
				 const Ref<GaussianBasisSet>& b3,
				 const Ref<GaussianBasisSet>& b4,
				 size_t storage, tbinteval int2etype,
                                 const ContractedGeminal& gbra,
                                 const ContractedGeminal& gket):
  TwoBodyInt(integral,b1,b2,b3,b4)
{
  // Which evaluator to use
  switch (int2etype) {
  case erieval:
    int2elibint2_ = new EriLibint2(integral,b1,b2,b3,b4,storage);
    num_tbint_types_ = 1;
    break;
  case g12eval:
    int2elibint2_ = new G12Libint2(integral,b1,b2,b3,b4,storage,gbra,gket);
    num_tbint_types_ = 6;
    break;
  default:
    throw FeatureNotImplemented("Tried to construct a two-electron integral \
evaluator of unimplemented or unknown type",__FILE__,__LINE__);
  }

  buffer_ = int2elibint2_->buffer();
  integral_->adjust_storage(int2elibint2_->storage_used());
}

TwoBodyIntLibint2::~TwoBodyIntLibint2()
{
  integral_->adjust_storage(-int2elibint2_->storage_used());
}

void
TwoBodyIntLibint2::compute_shell(int is, int js, int ks, int ls)
{
  int2elibint2_->set_redundant(redundant());
  int2elibint2_->compute_quartet(&is,&js,&ks,&ls);
}

int
TwoBodyIntLibint2::log2_shell_bound(int is, int js, int ks, int ls)
{
  return 10000000;//int2elibint2_->erep_4bound(is,js,ks,ls);
}

void
TwoBodyIntLibint2::set_integral_storage(size_t storage)
{
  int2elibint2_->init_storage(storage);
}

//////////////////////////////////////////////////////////////////////////

TwoBodyDerivIntLibint2::TwoBodyDerivIntLibint2(Integral*integral,
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

  //  int2elibint2_ = new EriLibint2(integral,b1,b2,b3,b4,1,storage);
  //buffer_ = int2elibint2_->buffer();
  //integral_->adjust_storage(int2elibint2_->storage_used());
}

TwoBodyDerivIntLibint2::~TwoBodyDerivIntLibint2()
{
  //  integral_->adjust_storage(-int2elibint2_->used_storage());
}

void
TwoBodyDerivIntLibint2::compute_shell(int is, int js, int ks, int ls,
                                 DerivCenters&c)
{
  int center;
  int sh[4], sz[4];

  sh[0]=is; sh[1]=js; sh[2]=ks; sh[3]=ls;

}


int
TwoBodyDerivIntLibint2::log2_shell_bound(int is, int js, int ks, int ls)
{
  return 0;//int2elibint2_->erep_4bound_1der(is,js,ks,ls);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
