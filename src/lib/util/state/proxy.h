//
// proxy.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifndef _util_state_proxy_h
#define _util_state_proxy_h

#include <util/class/proxy.h>
#include <util/state/state.h>

namespace sc {

/** @ingroup CoreState
 *  @ingroup CoreKeyVal
 *  Create a proxy for a SavableState object.
    This can be used to include an object that has been saved with
    SavableState in an input file for ParsedKeyVal.
*/
class SavableStateProxy: public DescribedClassProxy {
  private:
    Ref<SavableState> object_;
  public:
    /// Return the object referred to by this proxy.
    Ref<DescribedClass> object();
    /// Create the proxy with a Ref to a KeyVal object.
    SavableStateProxy(const Ref<KeyVal> &);
};

}

#endif
