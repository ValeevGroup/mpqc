//
// tbintlibint2.cc
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

#include <libint2.h>

#include <util/class/class.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/intdescr.h>
#include <chemistry/qc/libint2/tbintlibint2.h>
#include <chemistry/qc/libint2/bounds.h>
#include <chemistry/qc/libint2/bounds.timpl.h>
#if LIBINT2_SUPPORT_ERI
#  include <chemistry/qc/libint2/tbosar.h>
#  include <chemistry/qc/libint2/g12nc.h>
#endif
#if LIBINT2_SUPPORT_G12
# if LIBINT2_SUPPORT_T1G12
#  include <chemistry/qc/libint2/g12.h>
# endif
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
    struct Int2eCreator<G12NCLibint2> {
        Ref<G12NCLibint2> operator()(Integral*integral,
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
    };
#endif

#if LIBINT2_SUPPORT_G12
# if LIBINT2_SUPPORT_T1G12
	template<>
	struct Int2eCreator<G12Libint2> {
	Ref<G12Libint2> operator()(Integral*integral,
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
	};
# endif
#endif
#if LIBINT2_SUPPORT_G12DKH
	   template<>
	   struct Int2eCreator<G12DKHLibint2> {
	    Ref<G12DKHLibint2> operator()(Integral*integral,
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
	   };
#endif
    }
}

TwoBodyIntLibint2::TwoBodyIntLibint2(Integral*integral,
				     const Ref<GaussianBasisSet>& b1,
				     const Ref<GaussianBasisSet>& b2,
				     const Ref<GaussianBasisSet>& b3,
				     const Ref<GaussianBasisSet>& b4,
				     size_t storage, TwoBodyOperSet::type int2etype,
				     const Ref<IntParams>& params):
    TwoBodyInt(integral,b1,b2,b3,b4), int2etype_(int2etype),
    descr_(TwoBodyOperSetDescr::instance(int2etype)),
    params_(params)
{
    using sc::libint2::Int2eCreator;
  // Which evaluator to use
  switch (int2etype_) {
#if LIBINT2_SUPPORT_ERI
  case TwoBodyOperSet::ERI:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::eri> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::G12NC:
  {
    typedef G12NCLibint2 Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::R12_0_G12:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::r12_0_g12> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::R12_m1_G12:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::r12_m1_g12> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::G12_T1_G12:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::g12t1g12> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::DeltaFunction:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::delta> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
#endif
#if LIBINT2_SUPPORT_G12
# if LIBINT2_SUPPORT_T1G12
  case TwoBodyOperSet::G12:
  {
    typedef G12Libint2 Int2e;
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    break;
  }
# endif
#endif

#if LIBINT2_SUPPORT_G12DKH
  case TwoBodyOperSet::G12DKH:
  {
    typedef G12DKHLibint2 Int2e;
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    break;
  }
#endif
  default:
    throw FeatureNotImplemented("Tried to construct a two-electron integral \
evaluator of unimplemented or unknown type",__FILE__,__LINE__);
  }

  buffer_ = int2elibint2_->buffer(TwoBodyOper::eri);
  integral_->adjust_storage(int2elibint2_->storage_used());
}

TwoBodyIntLibint2::TwoBodyIntLibint2(const TwoBodyIntLibint2& other):
    TwoBodyInt(other.integral(),
               other.basis1(),
               other.basis2(),
               other.basis3(),
               other.basis4()),
    int2etype_(other.int2etype_),
    descr_(TwoBodyOperSetDescr::instance(other.int2etype_)),
    params_(other.params_),
    int2elibint2_(other.int2elibint2_->clone())
{
  buffer_ = int2elibint2_->buffer(TwoBodyOper::eri);
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

bool
TwoBodyIntLibint2::cloneable() const
{
  return true;
}

Ref<TwoBodyInt>
TwoBodyIntLibint2::clone()
{
  return new TwoBodyIntLibint2(*this);
}

//////////////////////////////////////////////////////////////////////////

TwoBodyThreeCenterIntLibint2::TwoBodyThreeCenterIntLibint2(Integral*integral,
                     const Ref<GaussianBasisSet>& b1,
                     const Ref<GaussianBasisSet>& b2,
                     const Ref<GaussianBasisSet>& b3,
                     size_t storage, TwoBodyOperSet::type int2etype,
                     const Ref<IntParams>& params):
    TwoBodyThreeCenterInt(integral,b1,b2,b3), int2etype_(int2etype),
    descr_(TwoBodyOperSetDescr::instance(int2etype)),
    params_(params)
{
  Ref<GaussianBasisSet> b4 = GaussianBasisSet::unit();
  using sc::libint2::Int2eCreator;
  // Which evaluator to use
  switch (int2etype_) {
#if LIBINT2_SUPPORT_ERI
  case TwoBodyOperSet::ERI:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::eri> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::G12NC:
  {
    typedef G12NCLibint2 Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::R12_0_G12:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::r12_0_g12> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::R12_m1_G12:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::r12_m1_g12> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::G12_T1_G12:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::g12t1g12> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::DeltaFunction:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::delta> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
#endif
#if LIBINT2_SUPPORT_G12
# if LIBINT2_SUPPORT_T1G12
  case TwoBodyOperSet::G12:
  {
    typedef G12Libint2 Int2e;
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    break;
  }
# endif
#endif

#if LIBINT2_SUPPORT_G12DKH
  case TwoBodyOperSet::G12DKH:
  {
    typedef G12DKHLibint2 Int2e;
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    break;
  }
#endif
  default:
    throw FeatureNotImplemented("Tried to construct a two-electron integral \
evaluator of unimplemented or unknown type",__FILE__,__LINE__);
  }

  buffer_ = int2elibint2_->buffer(TwoBodyOper::eri);
  integral_->adjust_storage(int2elibint2_->storage_used());
}

TwoBodyThreeCenterIntLibint2::TwoBodyThreeCenterIntLibint2(const TwoBodyThreeCenterIntLibint2& other) :
    TwoBodyThreeCenterInt(other.integral(),
                          other.basis1(),
                          other.basis2(),
                          other.basis3()),
    int2etype_(other.int2etype_),
    descr_(TwoBodyOperSetDescr::instance(other.int2etype_)),
    params_(other.params_),
    int2elibint2_(other.int2elibint2_->clone())
{
  buffer_ = int2elibint2_->buffer(TwoBodyOper::eri);
  integral_->adjust_storage(int2elibint2_->storage_used());
}

TwoBodyThreeCenterIntLibint2::~TwoBodyThreeCenterIntLibint2()
{
  integral_->adjust_storage(-int2elibint2_->storage_used());
}

void
TwoBodyThreeCenterIntLibint2::compute_shell(int is, int js, int ks)
{
  int ls = 0;
  int2elibint2_->set_redundant(redundant());
  int2elibint2_->compute_quartet(&is,&js,&ks,&ls);
}

int
TwoBodyThreeCenterIntLibint2::log2_shell_bound(int is, int js, int ks)
{
  const int bound = int2elibint2_->log2_bound(is,js,ks,0);
  return bound;
}

bool
TwoBodyThreeCenterIntLibint2::cloneable() const
{
  return true;
}

Ref<TwoBodyThreeCenterInt>
TwoBodyThreeCenterIntLibint2::clone()
{
  return new TwoBodyThreeCenterIntLibint2(*this);
}

//////////////////////////////////////////////////////////////////////////

TwoBodyTwoCenterIntLibint2::TwoBodyTwoCenterIntLibint2(Integral*integral,
                     const Ref<GaussianBasisSet>& b1,
                     const Ref<GaussianBasisSet>& b3,
                     size_t storage, TwoBodyOperSet::type int2etype,
                     const Ref<IntParams>& params):
    TwoBodyTwoCenterInt(integral,b1,b3), int2etype_(int2etype),
    descr_(TwoBodyOperSetDescr::instance(int2etype)),
    params_(params)
{
  Ref<GaussianBasisSet> b2 = GaussianBasisSet::unit();
  Ref<GaussianBasisSet> b4 = b2;
  using sc::libint2::Int2eCreator;
  // Which evaluator to use
  switch (int2etype_) {
#if LIBINT2_SUPPORT_ERI
  case TwoBodyOperSet::ERI:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::eri> Int2e;
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    break;
  }
  case TwoBodyOperSet::G12NC:
  {
    typedef G12NCLibint2 Int2e;
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    break;
  }
  case TwoBodyOperSet::R12_0_G12:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::r12_0_g12> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::R12_m1_G12:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::r12_m1_g12> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::G12_T1_G12:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::g12t1g12> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
  case TwoBodyOperSet::DeltaFunction:
  {
    typedef TwoBodyOSARLibint2<TwoBodyOper::delta> Int2e;
    typedef BoundsLibint2<Int2e> Bounds;
    Ref<Bounds> bounds = new Bounds(integral,b1,b2,b3,b4,storage,params);
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    int2elibint2_->bounds(bounds);
    break;
  }
#endif
#if LIBINT2_SUPPORT_G12
# if LIBINT2_SUPPORT_T1G12
  case TwoBodyOperSet::G12:
  {
    typedef G12Libint2 Int2e;
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    break;
  }
# endif
#endif

#if LIBINT2_SUPPORT_G12DKH
  case TwoBodyOperSet::G12DKH:
  {
    typedef G12DKHLibint2 Int2e;
    Int2eCreator<Int2e> creator;
    int2elibint2_ = creator(integral,b1,b2,b3,b4,storage,params);
    break;
  }
#endif
  default:
    throw FeatureNotImplemented("Tried to construct a two-electron integral \
evaluator of unimplemented or unknown type",__FILE__,__LINE__);
  }

  buffer_ = int2elibint2_->buffer(TwoBodyOper::eri);
  integral_->adjust_storage(int2elibint2_->storage_used());
}

TwoBodyTwoCenterIntLibint2::TwoBodyTwoCenterIntLibint2(const TwoBodyTwoCenterIntLibint2& other):
    TwoBodyTwoCenterInt(other.integral(),
                        other.basis1(),
                        other.basis2()),
    int2etype_(other.int2etype_),
    descr_(TwoBodyOperSetDescr::instance(other.int2etype_)),
    params_(other.params_),
    int2elibint2_(other.int2elibint2_->clone())
{
  buffer_ = int2elibint2_->buffer(TwoBodyOper::eri);
  integral_->adjust_storage(int2elibint2_->storage_used());
}

TwoBodyTwoCenterIntLibint2::~TwoBodyTwoCenterIntLibint2()
{
  integral_->adjust_storage(-int2elibint2_->storage_used());
}

void
TwoBodyTwoCenterIntLibint2::compute_shell(int is, int ks)
{
  int js = 0;
  int ls = 0;
  int2elibint2_->set_redundant(redundant());
  int2elibint2_->compute_quartet(&is,&js,&ks,&ls);
}

int
TwoBodyTwoCenterIntLibint2::log2_shell_bound(int is, int ks)
{
  // no computational savings from Cauchy bounds here!
  return 256; // 2^256 \approx 10^26
}

bool
TwoBodyTwoCenterIntLibint2::cloneable() const
{
  return true;
}

Ref<TwoBodyTwoCenterInt>
TwoBodyTwoCenterIntLibint2::clone()
{
  return new TwoBodyTwoCenterIntLibint2(*this);
}

//////////////////////////////////////////////////////////////////

TwoBodyDerivIntLibint2::TwoBodyDerivIntLibint2(Integral*integral,
					   const Ref<GaussianBasisSet>& b1,
					   const Ref<GaussianBasisSet>& b2,
					   const Ref<GaussianBasisSet>& b3,
					   const Ref<GaussianBasisSet>& b4,
					   size_t storage, TwoBodyOperSet::type int2etype):
  TwoBodyDerivInt(integral,b1,b2,b3,b4)
{
  // Which evaluator to use
  switch (int2etype) {
  default:
      throw FeatureNotImplemented("IntegralLibint2 does not yet implement geometrical derivatives. Try IntegralV3.",__FILE__,__LINE__);
  }
}

TwoBodyDerivIntLibint2::~TwoBodyDerivIntLibint2()
{
}

void
TwoBodyDerivIntLibint2::compute_shell(int is, int js, int ks, int ls,
                                 DerivCenters&c)
{
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
