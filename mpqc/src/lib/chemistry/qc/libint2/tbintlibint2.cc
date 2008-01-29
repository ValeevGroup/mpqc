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

#include <libint2/libint2.h>

#include <util/class/class.h>
#include <util/class/scexception.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/intdescr.h>
#include <chemistry/qc/libint2/tbintlibint2.h>
#include <chemistry/qc/libint2/bounds.h>
#include <chemistry/qc/libint2/bounds.timpl.h>
#if LIBINT2_SUPPORT_ERI
#  include <chemistry/qc/libint2/eri.h>
#endif
#if LIBINT2_SUPPORT_G12
# if LIBINT2_SUPPORT_T1G12
#  include <chemistry/qc/libint2/g12.h>
# else
#  include <chemistry/qc/libint2/g12nc.h>
# endif
#endif
#if LIBINT2_SUPPORT_GENG12
#  include <chemistry/qc/libint2/geng12.h>
#endif
#if LIBINT2_SUPPORT_G12DKH
#  include <chemistry/qc/libint2/g12dkh.h>
#endif

using namespace std;
using namespace sc;

namespace util {
    template <class To, class From>
    Ref<To> require_dynamic_cast(const Ref<From>& fptr) {
	To* tptr = dynamic_cast<To*>(fptr.pointer());
	if (!tptr)
	    throw ProgrammingError("Dynamic cast failed",__FILE__,__LINE__);
	return tptr;
    }
}

//////////////////////////////////////////////////////////////////////////

namespace sc {
    namespace libint2 {

#if LIBINT2_SUPPORT_ERI
	template<>
	Ref<EriLibint2>
	create_int2e(Integral*integral,
		     const Ref<GaussianBasisSet>& b1,
		     const Ref<GaussianBasisSet>& b2,
		     const Ref<GaussianBasisSet>& b3,
		     const Ref<GaussianBasisSet>& b4,
		     size_t storage,
		     const Ref<IntParams>& params)
	{
	    return new EriLibint2(integral,b1,b2,b3,b4,storage);
	}
#endif
	
#if LIBINT2_SUPPORT_G12
# if LIBINT2_SUPPORT_T1G12
	template<>
	Ref<G12Libint2>
	create_int2e(Integral*integral,
		     const Ref<GaussianBasisSet>& b1,
		     const Ref<GaussianBasisSet>& b2,
		     const Ref<GaussianBasisSet>& b3,
		     const Ref<GaussianBasisSet>& b4,
		     size_t storage,
		     const Ref<IntParams>& params)
	{
	    typedef IntParamsG12 IPType;
	    Ref<IPType> params_cast = util::require_dynamic_cast<IPType,IntParams>(params);
	    return new G12Libint2(integral,b1,b2,b3,b4,storage,params_cast->bra(),params_cast->ket());
	}
# else
	template<>
	Ref<G12NCLibint2>
	create_int2e(Integral*integral,
		     const Ref<GaussianBasisSet>& b1,
		     const Ref<GaussianBasisSet>& b2,
		     const Ref<GaussianBasisSet>& b3,
		     const Ref<GaussianBasisSet>& b4,
		     size_t storage,
		     const Ref<IntParams>& params)
	{
	    typedef IntParamsG12 IPType;
	    Ref<IPType> params_cast = util::require_dynamic_cast<IPType,IntParams>(params);
	    return new G12NCLibint2(integral,b1,b2,b3,b4,storage,params_cast->bra(),params_cast->ket());
	}
# endif
#endif
#if LIBINT2_SUPPORT_GENG12
	template<>
	Ref<GenG12Libint2>
	create_int2e(Integral*integral,
		     const Ref<GaussianBasisSet>& b1,
		     const Ref<GaussianBasisSet>& b2,
		     const Ref<GaussianBasisSet>& b3,
		     const Ref<GaussianBasisSet>& b4,
		     size_t storage,
		     const Ref<IntParams>& params)
	{
	    typedef IntParamsGenG12 IPType;
	    Ref<IPType> params_cast = util::require_dynamic_cast<IPType,IntParams>(params);
	    return new GenG12Libint2(integral,b1,b2,b3,b4,storage,params_cast->bra(),params_cast->ket());
	}
#endif
#if LIBINT2_SUPPORT_G12DKH
	   template<>
	    Ref<G12DKHLibint2>
	    create_int2e(Integral*integral,
	             const Ref<GaussianBasisSet>& b1,
	             const Ref<GaussianBasisSet>& b2,
	             const Ref<GaussianBasisSet>& b3,
	             const Ref<GaussianBasisSet>& b4,
	             size_t storage,
	             const Ref<IntParams>& params)
	    {
	        typedef IntParamsG12 IPType;
	        Ref<IPType> params_cast = util::require_dynamic_cast<IPType,IntParams>(params);
	        if (params_cast->bra() != params_cast->ket())
	          throw FeatureNotImplemented("G12DKH integrals are currently available for 1 correlation factor only",__FILE__,__LINE__);
	        return new G12DKHLibint2(integral,b1,b2,b3,b4,storage,params_cast->bra());
	    }
#endif
    }
}

TwoBodyIntLibint2::TwoBodyIntLibint2(Integral*integral,
				     const Ref<GaussianBasisSet>& b1,
				     const Ref<GaussianBasisSet>& b2,
				     const Ref<GaussianBasisSet>& b3,
				     const Ref<GaussianBasisSet>& b4,
				     size_t storage, tbinteval int2etype,
				     const Ref<IntParams>& params):
    TwoBodyInt(integral,b1,b2,b3,b4), int2etype_(int2etype)
{
    using sc::libint2::create_int2e;
  // Which evaluator to use
  switch (int2etype_) {
#if LIBINT2_SUPPORT_ERI
  case erieval:
  {
    typedef EriLibint2 Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_ = create_int2e<EriLibint2>(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    num_tbint_types_ = TwoBodyIntDescrERI::num_intsets;
    break;
  }
#endif
#if LIBINT2_SUPPORT_G12
# if LIBINT2_SUPPORT_T1G12
  case g12eval:
  {
    int2elibint2_ = create_int2e<G12Libint2>(integral,b1,b2,b3,b4,storage,params);
    num_tbint_types_ = TwoBodyIntDescrG12::num_intsets;
    break;
  }
# else
  case g12nceval:
  {
    int2elibint2_ = create_int2e<G12NCLibint2>(integral,b1,b2,b3,b4,storage,params);
    num_tbint_types_ = TwoBodyIntDescrG12NC::num_intsets;
    break;
  }
# endif
#endif
#if LIBINT2_SUPPORT_GENG12
  case geng12eval:
  {
    int2elibint2_ = create_int2e<GenG12Libint2>(integral,b1,b2,b3,b4,storage,params);
    num_tbint_types_ = TwoBodyIntDescrGenG12::num_intsets;
    break;
  }
#endif

#if LIBINT2_SUPPORT_G12DKH
  case g12dkheval:
  {
    int2elibint2_ = create_int2e<G12DKHLibint2>(integral,b1,b2,b3,b4,storage,params);
    num_tbint_types_ = TwoBodyIntDescrG12DKH::num_intsets;
    break;
  }
#endif
  default:
    throw FeatureNotImplemented("Tried to construct a two-electron integral \
evaluator of unimplemented or unknown type",__FILE__,__LINE__);
  }

  buffer_ = int2elibint2_->buffer(TwoBodyInt::eri);
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
  const int bound = int2elibint2_->log2_bound(is,js,ks,ls);
  return bound;
}

void
TwoBodyIntLibint2::set_integral_storage(size_t storage)
{
  int2elibint2_->init_storage(storage);
}

unsigned int
TwoBodyIntLibint2::inttype(TwoBodyInt::tbint_type type) const
{
    switch (int2etype_) {

    case erieval:
	return TwoBodyIntDescrERI::intSet(type);

    case g12eval:
	return TwoBodyIntDescrG12::intSet(type);

    case g12nceval:
	return TwoBodyIntDescrG12NC::intSet(type);

    case g12dkheval:
    return TwoBodyIntDescrG12DKH::intSet(type);

    case geng12eval:
	return TwoBodyIntDescrGenG12::intSet(type);
    }
}

TwoBodyInt::tbint_type
TwoBodyIntLibint2::inttype(unsigned int type) const
{
    switch (int2etype_) {

    case erieval:
	return TwoBodyIntDescrERI::intSet(type);

    case g12eval:
	return TwoBodyIntDescrG12::intSet(type);

    case g12nceval:
	return TwoBodyIntDescrG12NC::intSet(type);

    case g12dkheval:
    return TwoBodyIntDescrG12DKH::intSet(type);

    case geng12eval:
	return TwoBodyIntDescrGenG12::intSet(type);
    }
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
      throw ProgrammingError("TwoBodyDerivIntLibint2: tried to construct a two-electron derivative integral evaluator of unimplemented or unknown type",__FILE__,__LINE__);
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
  return int2elibint2_->log2_bound(is,js,ks,ls);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
