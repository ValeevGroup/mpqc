//
// bounds.h
//
// Copyright (C) 2007 Edward Valeev
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
#pragma interface
#endif

#ifndef _chemistry_qc_libint2_boundstimpl_h
#define _chemistry_qc_libint2_boundstimpl_h

#include <vector>
#include <util/class/scexception.h>
#include <chemistry/qc/libint2/int2e.h>

namespace libint2 {
    template <typename T>
    struct abs : public std::unary_function<const T&,T> {
	T operator()(const T& x) {
	    if (x < 0) return std::negate<T>()(x);
	}
    };
}

namespace sc {

    template <class Int2e>
    BoundsLibint2<Int2e>::BoundsLibint2(Integral*integral,
					const Ref<GaussianBasisSet>& b1,
					const Ref<GaussianBasisSet>& b2,
					const Ref<GaussianBasisSet>& b3,
					const Ref<GaussianBasisSet>& b4,
					size_t storage,
					const Ref<IntParams>& params)
    {
	equiv_12_34_ = b1->equiv(b3) && b2->equiv(b4);
	equiv_12_43_ = b1->equiv(b4) && b2->equiv(b3);

	Ref<MessageGrp> msg = integral->messagegrp();
	const int ntasks = msg->n();
	const int me = msg->me();

	Ref<Int2e> int12 = libint2::create_int2e<Int2e>(integral,b1,b2,b1,b2,storage,params);
	const int nsh1 = b1->nshell();
	const int nsh2 = b2->nshell();
	const int n12 = nsh1*nsh2;
	Q12_.resize(n12);  for(int i=0; i<n12; ++i) Q12_[i] = 0.0;
        double* buf12 = int12->buffer(0);
	int f12 = 0;
	for(int s1=0; s1<nsh1; ++s1) {
	    const int nf1 = b1->shell(s1).nfunction();
	    for(int s2=0; s2<nsh2; ++s2, ++f12) {
		const int nf2 = b2->shell(s2).nfunction();
		const int nf12 = nf1*nf2;

		if (f12%ntasks != me)
		    continue;

		int12->compute_quartet(&s1,&s2,&s1,&s2);

		std::transform(buf12,buf12+nf12,buf12,::libint2::abs<double>());
		const double max_elem = *(std::max_element(buf12,buf12+nf12));
		Q12_[f12] = Log2Bounds::bound_cast(max_elem);
	    }
	}
	// propagate Q12_ among the nodes
	{
	    double* q12 = new double[n12];
	    std::copy(Q12_.begin(),Q12_.end(),q12);
	    msg->sum(q12,n12);
	    std::copy(q12,q12+n12,Q12_.begin());
	}

	if (!equiv_12_34_ && !equiv_12_43_) {
	    throw FeatureNotImplemented("general basis set combinations are not handled yet in bounds computation",__FILE__,__LINE__);
	}
	else {
	    if (equiv_12_34_) {
		Q34_ = Q12_;
	    }
	    if (equiv_12_43_) {
		Q34_.resize(n12);
		int f12 = 0;
		for(int s2=0; s2<nsh2; ++s2) {
		    for(int s1=0; s1<nsh1; ++s1, ++f12) {
			Q34_[f12] = Q12_[s1*nsh2 + s2];
		    }
		}
	    }
	}
    }

    template <class Int2e>
    BoundsLibint2<Int2e>::~BoundsLibint2()
    {
    }

    template <class Int2e>
    Log2Bounds::int_bound_t
    BoundsLibint2<Int2e>::log2_bound(int sh1, int sh2, int sh3, int sh4) const
    {
	return Q12_[nsh2_*sh1 + sh2] + Q34_[nsh4_*sh3 + sh4];
    }
}

#endif // header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
