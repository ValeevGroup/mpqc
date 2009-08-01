//
// twobodygrid.h
//
// Copyright (C) 2004 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_twobodygrid_h
#define _chemistry_qc_mbptr12_twobodygrid_h

#include <vector>
#include <util/ref/ref.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/molecule.h>

namespace sc {

  /** Class TwoBodyGrid describes a set of coordinates of 2 particles.
      */

class TwoBodyGrid : virtual public SavableState {

  std::string name_;               // String identifier for the grid
  std::vector<SCVector3> r1_;           // Cartesian coordinates of particle 1
  std::vector<SCVector3> r2_;           // Cartesian coordinates of particle 2
  SCVector3 O_;                   // The origin with respect to which r1 and r2 are defined

public:
  TwoBodyGrid(StateIn&);
  /** The KeyVal constructor

      <dl>

      <dt><tt>name</tt><dd> This string specifies the name for this grid.

      <dt><tt>origin</tt><dd> The 3-element vector specifies the Cartesian coordinates
      of the origin of the grid. The default is the origin of the reference coordinate frame.

      <dt><tt>polar</tt><dd> This boolean specifies whether coordinates of particles 1 and 2
      are specified as polar or Cartesian vectors. The default is no (Cartesian).

      <dt><tt>r1</tt><dd> This vector of 3-element vectors specifies the positions
      of particle 1. There is no default.

      <dt><tt>r2</tt><dd> See <tt>r1</tt>.

      <dt><tt>ndim</tt><dd> This integer specifies the number of dimensions for this grid.
      Valid values are 1 and 2. The default is 2, i.e. for every position of particle 1
      every position of particle 2 is included. Value of 1 is only valid if the number of elements
      in <tt>r1</tt> and <tt>r2</tt> is the same -- then every r1[i] r2[i] pair of positions is included.

      </dl>
    */
  TwoBodyGrid(const Ref<KeyVal>&);
  ~TwoBodyGrid();

  void save_data_state(StateOut&);

  const std::string& name() const;
  int nelem() const;
  const SCVector3& origin() const;

  /// Cartesian coordinates of particle 1 at grid point i with respect to origin O
  SCVector3 xyz1(int i, const SCVector3& O=SCVector3(0.0)) const;
  /// Cartesian coordinates of particle 2 at grid point i with respect to origin O
  SCVector3 xyz2(int i, const SCVector3& O=SCVector3(0.0)) const;
  /// Spherical polar coordinates of particle 1 at grid point i with respect to origin O
  SCVector3 rtp1(int i, const SCVector3& O=SCVector3(0.0)) const;
  /// Spherical polar coordinates of particle 2 at grid point i with respect to origin O
  SCVector3 rtp2(int i, const SCVector3& O=SCVector3(0.0)) const;
  /// Prints out this
  void print(std::ostream&o=ExEnv::out0()) const;

};


}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


