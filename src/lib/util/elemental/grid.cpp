//
// grid.cpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma implementation
#endif

#include <util/elemental/grid.hpp>

using namespace mpqc;
using namespace sc;

ClassDesc Grid::class_desc_(typeid(Grid),
                           "Grid",
                           1,
                           "virtual public DescribedClass",
                           create<Grid>,
                           create<Grid>,
                           0
);

mpqc::Grid::Grid() : grid_(new elem::Grid(elem::mpi::COMM_WORLD)) { }

mpqc::Grid::Grid(const sc::Ref<sc::KeyVal>& kv) {
  MPQC_ASSERT(false);
}

mpqc::Grid::~Grid() {
  // Copying the madness style world.cc
  const bool make_sure_class_desc_initialialized = (&class_desc_ != 0);
}
