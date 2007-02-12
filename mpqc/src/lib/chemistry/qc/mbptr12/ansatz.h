//
// ansatz.h
//
// Copyright (C) 2006 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#ifndef _chemistry_qc_mbptr12_ansatz_h
#define _chemistry_qc_mbptr12_ansatz_h

#include <util/keyval/keyval.h>
#include <util/state/state.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/mbptr12/linearr12.h>

#ifdef __GNUC__
#pragma interface
#endif

namespace sc {
  
  class LinearR12Ansatz : virtual public SavableState {
    public:
    /** The KeyVal constructor.
    <dl>
    
    <dt><tt>orbital_product</tt><dd> This specifies how the geminal space is generated.
    Geminal functions are products of the correlation factor and 2 orbitals.
    This keyword specifies which orbital products are allowed.
    Valid choices are:
      <dl>
        <dt><tt>ij</tt><dd> Biproducts of occupied orbitals. This is the default.
        <dt><tt>pq</tt><dd> Biproducts of any Hartree-Fock orbitals. This has not been implemented yet.
      </dl>

    <dt><tt>projector</tt><dd> This specifies the form of the orthogonal projector.
    Valid values are:
      <dl>
        <dt><tt>2</tt><dd> (1-O1)(1-O2)(1-V1V2). This is the default.
        <dt><tt>3</tt><dd> 1-P1P2. Should be used ONLY for testing.
      </dl>

    <dt><tt>diag</tt><dd> Setting this to <tt>true</tt> will only keep the diagonal terms,
    which is equivalent to the "old" (pre-1992) form of R12 theory. The default is <tt>false</tt>,
    which corresponds to the orbital invariant ansatz of Klopper.
    
    </dl>
    */
    LinearR12Ansatz(const Ref<KeyVal>&);
    /// The StateIn constructor
    LinearR12Ansatz(StateIn&);
    /// The default constructor creates orbital-invariant ansatz with projector 2
    LinearR12Ansatz();
    ~LinearR12Ansatz();
    
    void save_data_state(StateOut&);
    void print(std::ostream& o =ExEnv::out0()) const;
    
    LinearR12::Projector projector() const;
    bool diag() const;
    LinearR12::OrbitalProduct orbital_product() const;
    
    private:
    LinearR12::Projector projector_;
    bool diag_;
    LinearR12::OrbitalProduct orbital_product_;
  };
  
}

#endif // include guard

