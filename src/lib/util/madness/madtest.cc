//
// madtest.cc
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
#pragma implementation
#endif

#include <util/madness/init.h>
#include <util/madness/world.h>
#include <util/misc/exenv.h>
#define BOOST_TEST_MODULE test_world
#include <boost/test/included/unit_test.hpp>

using namespace sc;
using namespace mpqc;
using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE(test_world_constructor){
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char** argv = boost::unit_test::framework::master_test_suite().argv;
  ExEnv::init(argc, argv);
  MADNESSRuntime::initialize();

  mpqc::World world0;

#if MPQC_HAS_ELEMENTAL
  elem::DistMatrix<double> mat(5, 5);
#endif // MPQC_HAS_ELEMENTAL


  Ref<AssignedKeyVal> akv1 = new AssignedKeyVal;
  mpqc::World world1(akv1);
  assert(world0.madworld() == world1.madworld());
  mpqc::World world2(akv1);
  assert(world1.madworld() == world2.madworld());

  Ref<AssignedKeyVal> akv3 = new AssignedKeyVal;
  akv3->assign("key", "this is my world");
  // mpqc::World world3(akv3); // ERROR for now only default world is allowed

  MADNESSRuntime::finalize();
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
