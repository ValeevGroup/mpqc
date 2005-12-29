//
// csgrade12.h
// based on csgrad.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
// Maintainer: LPS
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

#ifndef _chemistry_qc_mbpt_csgrade12_h
#define _chemistry_qc_mbpt_csgrade12_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/misc/regtime.h>
#include <util/group/memory.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/distshpair.h>

namespace sc {

#define PRINT_BIGGEST_INTS 0

class CSGradErep12Qtr: public Thread {
  private:
    Ref<MessageGrp> msg;
    Ref<MemoryGrp> mem;
    Ref<TwoBodyInt> tbint;
    Ref<GaussianBasisSet> basis;
    Ref<ThreadLock> lock;
    Ref<RegionTimer> timer;
    int mythread;
    int nthread;
    int ni;
    int nocc;
    int i_offset;
    int aoint_computed;
    int me;
    int nproc;
    double tol;
    double **scf_vector;
    int debug;
    int dynamic_;
    double print_percent_;
    int usep4_;
    DistShellPair::SharedData *shellpair_shared_data_;
  public:
    CSGradErep12Qtr(int mythread_a, int nthread_a,
                    int me_a, int nproc_a,
                    const Ref<MemoryGrp> &mem_a,
                    const Ref<MessageGrp> &msg_a,
                    const Ref<ThreadLock> &lock_a,
                    const Ref<GaussianBasisSet> &basis_a,
                    const Ref<TwoBodyInt> &tbint_a,
                    int nocc_a,
                    double **scf_vector_a,
                    double tol_a, int debug_a,
                    int dynamic_a, double print_percent_a,
                    DistShellPair::SharedData *shellpair_shared_data,
                    int usep4);
    ~CSGradErep12Qtr();

    void set_i_offset(int ioff) { i_offset = ioff; }
    void set_ni(int nivalue) { ni = nivalue; }

    void run();
};

}

#endif

// //////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
