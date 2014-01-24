//
// assert.h
//
// Copyright (C) 2014 Edward Valeev
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

#ifndef mpqc_util_misc_assert_h
#define mpqc_util_misc_assert_h

#include <mpqc_config.h>

// null assert
#if MPQC_ASSERT_MODE == 0
#  define MPQC_ASSERT( a )
#endif

// std assert
#if MPQC_ASSERT_MODE == 1
#  include <cassert>
#  define MPQC_ASSERT( a ) assert(a)
#endif

// throw
#if MPQC_ASSERT_MODE == 2
#  include <util/misc/exception.h>
#  define MPQC_ASSERT( a ) if (not (a) ) throw sc::Exception("assertion failed", __FILE__, __LINE__)
#endif

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
