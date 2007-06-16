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

#ifndef _chemistry_qc_libint2_bounds_h
#define _chemistry_qc_libint2_bounds_h

#include <util/class/scexception.h>
#include <chemistry/qc/basis/intparams.h>
#include <chemistry/qc/libint2/int2e.h>

namespace sc {

    /// Computes log2 bounds
    class Log2Bounds : virtual public RefCount {
    protected:
        // Set to non-zero to debug this and derived classes
        static const int debugclass_ = 2;
    public:
	typedef signed char int_bound_t;
	enum { int_bound_min = SCHAR_MIN, int_bound_max = SCHAR_MAX };
	Log2Bounds() {}
	virtual ~Log2Bounds() {}

	virtual int_bound_t log2_bound(int sh1, int sh2, int sh3, int sh4) const =0;
	static int_bound_t bound_cast(double);
    };

    /// Computes log2 bounds for a particular Int2e evaluator
    template <class Int2e>
    class BoundsLibint2 : public Log2Bounds {
    public:
	typedef Log2Bounds::int_bound_t int_bound_t;

	BoundsLibint2(Integral*integral,
		     const Ref<GaussianBasisSet>& b1,
		     const Ref<GaussianBasisSet>& b2,
		     const Ref<GaussianBasisSet>& b3,
		     const Ref<GaussianBasisSet>& b4,
		     size_t storage,
		     const Ref<IntParams>& params);
	~BoundsLibint2();

	/// Implements Log2Bounds::log2_bound()
	int_bound_t log2_bound(int sh1, int sh2, int sh3, int sh4) const;

    private:
	std::vector<int_bound_t> Q12_;
	std::vector<int_bound_t> Q34_;
	int nsh2_;
	int nsh4_;
	bool equiv_12_34_;
	bool equiv_12_43_;

    };

}

#endif // header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
