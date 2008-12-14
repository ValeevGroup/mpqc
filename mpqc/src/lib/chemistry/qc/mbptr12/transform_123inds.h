//
// transform_123inds.h
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

#ifndef _chemistry_qc_mbpt_transform123inds_h
#define _chemistry_qc_mbpt_transform123inds_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/misc/regtime.h>
#include <util/group/memory.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/transform_tbint.h>

namespace sc {

#define PRINT_BIGGEST_INTS 0

class TwoBodyMOIntsTransform_123Inds: public Thread {

    Ref<TwoBodyMOIntsTransform> tform_;
    Ref<TwoBodyInt> tbint_;
    Ref<ThreadLock> lock_;
    Ref<RegionTimer> timer_;

    int mythread_;
    int nthread_;
    int ni_;        // Number of i-indices handled in each pass
    int i_offset_;  // first i-index handled in this pass

    double tol_;
    int debug_;

    int aoint_computed_;

  public:
    TwoBodyMOIntsTransform_123Inds(const Ref<TwoBodyMOIntsTransform>& tform,
    int mythread, int nthread, const Ref<ThreadLock>& lock, const Ref<TwoBodyInt> &tbint,
    double tol, int debug);
    ~TwoBodyMOIntsTransform_123Inds();

    void set_i_offset(const int ioff) { i_offset_ = ioff; }
    void set_ni(const int nivalue) { ni_ = nivalue; }

    void run();
};

}

#endif

// //////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
