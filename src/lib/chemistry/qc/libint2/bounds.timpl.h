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

#ifndef _chemistry_qc_libint2_boundstimpl_h
#define _chemistry_qc_libint2_boundstimpl_h

#include <vector>
#include <cmath>
#include <algorithm>
#include <util/misc/scexception.h>
#include <chemistry/qc/libint2/int2e.h>

namespace {
    template <typename T>
    struct abssqrt : public std::unary_function<const T&,T> {
	T operator()(const T& x) {
	    if (x < 0) return std::sqrt(std::negate<T>()(x));
	    else return std::sqrt(x);
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
					const Ref<IntParams>& params) :
    nsh2_(b2->nshell()), nsh4_(b4->nshell())
    {
	equiv_12_34_ = b1->equiv(b3) && b2->equiv(b4);
	equiv_12_43_ = b1->equiv(b4) && b2->equiv(b3);
        equiv_1_2_ = b1->equiv(b2);
        equiv_3_4_ = b3->equiv(b4);

	Ref<MessageGrp> msg = integral->messagegrp();
	const int ntasks = msg->n();
	const int me = msg->me();

	const int nsh1 = b1->nshell();
	const int nsh2 = b2->nshell();
	const int n12 = nsh1*nsh2;
	{
	  libint2::Int2eCreator<Int2e> creator;
	    Ref<Int2e> int12 = creator(integral,b1,b2,b1,b2,storage,params);
	    Q12_.resize(n12);  for(int i=0; i<n12; ++i) Q12_[i] = (int_bound_t)0;
	    double* buf = int12->buffer(0);
	    int f12 = 0;
	    for(int s1=0; s1<nsh1; ++s1) {
		const int nf1 = b1->shell(s1).nfunction();
		for(int s2=0; s2<nsh2; ++s2, ++f12) {
		    const int nf2 = b2->shell(s2).nfunction();
		    const int nf = nf1*nf2*nf1*nf2;
		    
		    if (f12%ntasks != me)
			continue;
		    
		    int12->compute_quartet(&s1,&s2,&s1,&s2);
		    
		    std::transform(buf,buf+nf,buf,abssqrt<double>());
		    const double max_elem = *(std::max_element(buf,buf+nf));
		    Q12_[f12] = Log2Bounds::bound_cast(max_elem);
		}
	    }
	}
	// propagate Q12_ among the nodes
	{
	    int_bound_t* q12 = new int_bound_t[n12];
	    std::copy(Q12_.begin(),Q12_.end(),q12);
	    msg->sum(q12,n12);
	    std::copy(q12,q12+n12,Q12_.begin());
	}

	if (!equiv_12_34_ && !equiv_12_43_) {

        const int nsh3 = b3->nshell();
        const int nsh4 = b4->nshell();
        const int n34 = nsh3*nsh4;
        {
          libint2::Int2eCreator<Int2e> creator;
            Ref<Int2e> int34 = creator(integral,b3,b4,b3,b4,storage,params);
            Q34_.resize(n34);  for(int i=0; i<n34; ++i) Q34_[i] = (int_bound_t)0;
            double* buf = int34->buffer(0);
            int f34 = 0;
            for(int s3=0; s3<nsh3; ++s3) {
                const int nf3 = b3->shell(s3).nfunction();
                for(int s4=0; s4<nsh4; ++s4, ++f34) {
                    const int nf4 = b4->shell(s4).nfunction();
                    const int nf = nf3*nf4*nf3*nf4;

                    if (f34%ntasks != me)
                        continue;

                    int34->compute_quartet(&s3,&s4,&s3,&s4);

                    std::transform(buf,buf+nf,buf,abssqrt<double>());
                    const double max_elem = *(std::max_element(buf,buf+nf));
                    Q34_[f34] = Log2Bounds::bound_cast(max_elem);
                }
            }
        }
        // propagate Q34_ among the nodes
        {
            int_bound_t* q34 = new int_bound_t[n34];
            std::copy(Q34_.begin(),Q34_.end(),q34);
            msg->sum(q34,n34);
            std::copy(q34,q34+n34,Q34_.begin());
        }

	}
	else {
	    if (equiv_12_34_) {
		Q34_ = Q12_;
	    }
	    else if (equiv_12_43_) {
		Q34_.resize(n12);
		int f21 = 0;
		for(int s2=0; s2<nsh2; ++s2) {
		    for(int s1=0; s1<nsh1; ++s1, ++f21) {
			Q34_[f21] = Q12_[s1*nsh2 + s2];
		    }
		}
	    }
	}

    if (debugclass_ > 0) {
        ExEnv::out0() << indent << "BoundsLibint2::Q12 :" << std::endl;
        ExEnv::out0() << indent << "bs1:" << std::endl;
        b1->print_brief(ExEnv::out0());
        ExEnv::out0() << indent << "bs2:" << std::endl;
        b2->print_brief(ExEnv::out0());
        for(int s1=0; s1<nsh1; ++s1) {
            const int s2max = equiv_1_2_ ? s1 : nsh2-1;
            for(int s2=0; s2<=s2max; ++s2) {
                const int f12 = s1*nsh2 + s2;
                ExEnv::out0() << indent << s1 << " " << s2 << "  "
                              << int(Q12_[f12]) << std::endl;
            }
        }
    }

    if (debugclass_ > 0 && not equiv_12_34_ && not equiv_12_43_) {
      const int nsh3 = b3->nshell();
      const int nsh4 = b4->nshell();
        ExEnv::out0() << indent << "BoundsLibint2::Q34 :" << std::endl;
        ExEnv::out0() << indent << "bs3:" << std::endl;
        b3->print_brief(ExEnv::out0());
        ExEnv::out0() << indent << "bs4:" << std::endl;
        b4->print_brief(ExEnv::out0());
        for(int s3=0; s3<nsh3; ++s3) {
            const int s4max = equiv_3_4_ ? s3 : nsh4-1;
            for(int s4=0; s4<=s4max; ++s4) {
                const int f34 = s3*nsh4 + s4;
                ExEnv::out0() << indent << s3 << " " << s4 << "  "
                              << int(Q34_[f34]) << std::endl;
            }
        }
    }

    }

    template <class Int2e>
    BoundsLibint2<Int2e>::~BoundsLibint2()
    {
    }

    template <class Int2e>
    int
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
