//
// transform.h
//
// Copyright (C) 2009 Edward Valeev
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_transform_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_transform_h


namespace sc {

  struct MOIntsTransform {

    /// Describes the method of storing transformed MO integrals.
    struct StoreMethod {
      enum type { mem_posix = 0, posix = 1, mem_mpi = 2, mpi = 3, mem_only = 4 };
    };
    /// How integrals are stored. Type_13 means (ix|jy) integrals are stored as (ij|xy)
    enum StorageType {StorageType_First=0, StorageType_Last=1,
                      StorageType_12=0, StorageType_13=1};
    /// enumerates all known 4-center transform types
    enum TwoBodyTransformType {
      TwoBodyTransformType_ixjy=0, TwoBodyTransformType_ikjy=1,
      TwoBodyTransformType_ijxy=2, TwoBodyTransformType_iRjS=3,
      TwoBodyTransformType_ixjy_df=4,
      TwoBodyTransformType_ijR=5,
      TwoBodyTransformType_First=TwoBodyTransformType_ixjy,
      TwoBodyTransformType_Last=TwoBodyTransformType_ijR
    };
  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
