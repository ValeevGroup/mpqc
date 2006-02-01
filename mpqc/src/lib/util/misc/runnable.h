//
// runnable.h
//
// Author: Curtis Janssen <cljanss@sandia.gov>
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
#pragma interface
#endif

#ifndef _util_misc_runnable_h
#define _util_misc_runnable_h

#include <string>

#include <util/class/class.h>

namespace sc {

    /** The Runnable class is a DescribedClass with a pure virtual
        run member. */
    class Runnable: virtual public DescribedClass {
      public:
        /// Executates an action as specified the derived class.
        virtual void run() = 0;
    };

    class TestRunnable: public Runnable {
        std::string test_;
      public:
        TestRunnable(const Ref<KeyVal> &);
        void run();
    };

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
