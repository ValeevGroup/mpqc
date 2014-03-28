//
// world.h
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

#ifndef _mpqc_src_lib_util_madness_world_h
#define _mpqc_src_lib_util_madness_world_h

#include <world/worldfwd.h>
#include <util/class/class.h>
#include <util/keyval/keyval.h>

#ifdef MPQC_HAS_ELEMENTAL
#include <elemental.hpp>
#endif // MPQC_HAS_ELEMENTAL

namespace mpqc {

  /// World is a wrapper around madness::World
  class World : virtual public sc::DescribedClass {
    public:
      World();
      /** A KeyVal constructor is used to generate a World
          object from the input. The full list of keywords
          that are accepted is below.

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

          <tr><td><tt>key</tt><td>string<td>"default"<td>World key used to associate World objects with
          parts of the computation.

          </table>
       */
      World(const sc::Ref<sc::KeyVal>& kv);
      ~World();

      const std::string& key() const { return key_; }
      const madness::World* madworld() const { return world_; }
      madness::World* madworld() { return world_; }

#ifdef MPQC_HAS_ELEMENTAL
      const elem::Grid* elemGrid() const {return grid_;}
#endif // MPQC_HAS_ELEMENTAL

    private:
      static sc::ClassDesc class_desc_;

      std::string key_;
      madness::World* world_;

#ifdef MPQC_HAS_ELEMENTAL
      const elem::Grid* grid_;
#endif // MPQC_HAS_ELEMENTAL

  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
