//
// grid.hpp
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

#ifndef MPQC_SRC_LIB_UTIL_ELEMENTAL_GRID_HPP
#define MPQC_SRC_LIB_UTIL_ELEMENTAL_GRID_HPP

#include <elemental.hpp>
#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <util/misc/assert.h>

namespace mpqc {
  /// Grid is a wrapper around elem::Grid
  class Grid : virtual public sc::DescribedClass {
  public:
    Grid();
    /**
     * A KeyVal Constructor
     */
    Grid(const sc::Ref<sc::KeyVal> &kv);
    ~Grid();

    const std::shared_ptr<elem::Grid> elemGrid() const {return grid_;}
    std::shared_ptr<elem::Grid> elemGrid() {return grid_;}

  private:
    static sc::ClassDesc class_desc_;
    std::shared_ptr<elem::Grid> grid_;
  };
}


#endif /* MPQC_SRC_LIB_UTIL_ELEMENTAL_GRID_HPP */
