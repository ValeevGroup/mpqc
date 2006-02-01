//
// compute.cc
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
#pragma implementation
#endif

#include <util/keyval/keyval.h>
#include <util/misc/runnable.h>

static sc::ClassDesc Runnable_cd(typeid(sc::Runnable),
                                 "Runnable",1,
                                 "virtual public DescribedClass");


static sc::ClassDesc TestRunnable_cd(typeid(sc::TestRunnable),
                                     "TestRunnable",1,
                                     "public Runnable",
                                     0,
                                     sc::create<sc::TestRunnable>);


sc::TestRunnable::TestRunnable(const Ref<sc::KeyVal> &keyval):
  test_("default")
{
  test_ = keyval->stringvalue("test");
}

void
sc::TestRunnable::run()
{
  sc::ExEnv::out0() << indent
                    << "TestRunnable: test = " << test_ << std::endl;
}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
