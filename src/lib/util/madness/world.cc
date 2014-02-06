//
// world.cc
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

// includes go here
#include <util/madness/world.h>

using namespace mpqc;
using namespace sc;

ClassDesc
World::class_desc_(typeid(World),
                   "World",
                   1,               // version
                   "virtual public DescribedClass", // must match parent
                   create<World>,   // class is DefaultConstructible
                   create<World>,   // class is not KeyValConstructible
                   0  // class is not StateInConstructible
);

World::World() : key_("default"), world_(&madness::World::get_default()) {}

World::World(const Ref<KeyVal>& kv) {
  key_ = kv->stringvalue("key", KeyValValuestring("default"));
  // for now only use default world
  MPQC_ASSERT(key_ == "default");
  world_ = &madness::World::get_default();
}

World::~World() {
  // this may be necessary if this is a templated class
  const bool make_sure_class_desc_initialized = (&class_desc_ != 0);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
