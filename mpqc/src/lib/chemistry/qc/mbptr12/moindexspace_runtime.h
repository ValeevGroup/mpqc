//
// moindexspace_runtime.h
//
// Copyright (C) 2008 Edward Valeev
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_moindexspaceruntime_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_moindexspaceruntime_h

#include <chemistry/qc/mbptr12/moindexspace.h>
#include <chemistry/qc/mbptr12/registry.h>

namespace sc {

  /** Smart runtime that knows how to build MOIndexSpace objects and knows relationships
      between MOIndexSpace objects. The objects are held by the global MOIndexSpaceRegistry.
    */
  class MOIndexSpaceRuntime : virtual public SavableState {
    public:
      MOIndexSpaceRuntime();
      MOIndexSpaceRuntime(StateIn&);
      void save_data_state(StateOut&);

      /**
       Finds or constructs a space with the given key. To be able to
       construct a transformed space base spaces need to be available.
       */
      Ref<MOIndexSpace> get(const std::string& key) const;

      /// declares that space represents the AO space supported by space->basis()
      void declare_aospace(const Ref<MOIndexSpace>& space);
      /// which space is the AO space for basis?
      Ref<MOIndexSpace> aospace(const Ref<GaussianBasisSet>& basis);

    private:

      typedef Registry<Ref<GaussianBasisSet>, Ref<MOIndexSpace>, detail::NonsingletonCreationPolicy > BasisToAOSpaceMap;
      Ref<BasisToAOSpaceMap> basis_to_aospace_map_;

  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
