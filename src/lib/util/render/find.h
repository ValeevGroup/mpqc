//
// find.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#ifndef _util_render_find_h
#define _util_render_find_h

#include <util/render/stack.h>
#include <util/render/parameter.h>

namespace sc {

// cannot be used with g++ 2.6-94q4 and has other bugs anyway
#if 0
template <class T1, class T2>
void
find_parameter_in_stack(Stack<T1>& stack,
                        Parameter<T2>& (T1::*access)(),
                        T2& result
                        )
{
  int have_result = 0;
  for (int i=stack.n()-1; i>=0; i--) {
      if ((stack[i]->*access)().is_set()) {
          if (!have_result || (stack[i]->*access)().overrides()) {
              result = (stack[i]->*access)().value();
              have_result = 1;
            }
        }
    }
}
#endif

inline void
find_int_parameter_in_appearance_stack(Stack<Ref<Appearance> >& stack,
                        Parameter<int>& (Appearance::*access)(),
                        int& result
                        )
{
  int have_result = 0;
  for (int i=stack.n()-1; i>=0; i--) {
      if ((stack[i].pointer()->*access)().is_set()) {
          if (!have_result || (stack[i].pointer()->*access)().overrides()) {
              result = (stack[i].pointer()->*access)().value();
              have_result = 1;
            }
        }
    }
}

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
