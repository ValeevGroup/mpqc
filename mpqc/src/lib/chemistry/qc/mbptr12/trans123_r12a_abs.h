//
// trans123_r12a_abs.h
//
// Copyright (C) 2003 Edward Valeev
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

#ifndef _chemistry_qc_mbpt_trans123r12aabs_h
#define _chemistry_qc_mbpt_trans123r12aabs_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/misc/regtime.h>
#include <util/group/memory.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/integral.h>

namespace sc {

#define PRINT_BIGGEST_INTS 0

class R12A_ABS_123Qtr: public Thread {
  private:
    Ref<MessageGrp> msg;
    Ref<MemoryGrp> mem;
    Ref<TwoBodyInt> tbint;
    Ref<GaussianBasisSet> basis;
    Ref<GaussianBasisSet> aux_basis;
    Ref<ThreadLock> lock;
    Ref<RegionTimer> timer;
    int mythread;
    int nthread;
    int ni;
    int nocc;
    int nocc_act;
    int i_offset;
    int aoint_computed;
    int me;
    int nproc;
    double tol;
    double **scf_vector;
    int debug;
    int dynamic_;
    double print_percent_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    Ref<GaussianBasisSet> bs3_;
    Ref<GaussianBasisSet> bs4_;

  public:
    R12A_ABS_123Qtr(int mythread, int nthread,
		    int me, int nproc,
		    const Ref<MemoryGrp> &mem,
		    const Ref<MessageGrp> &msg,
		    const Ref<ThreadLock> &lock,
		    const Ref<GaussianBasisSet> &basis,
		    const Ref<GaussianBasisSet> &aux_basis,
		    const Ref<TwoBodyInt> &tbint,
		    int nocc, int nocc_act,
		    double **scf_vector,
		    double tol, int debug,
		    int dynamic, double print_percent);
    ~R12A_ABS_123Qtr();

    void set_i_offset(int ioff) { i_offset = ioff; }
    void set_ni(int nivalue) { ni = nivalue; }

    /// Returns the amount of memory needed for this thread
  //static size_t storage_required(int nocc_act, int nocc, Ref<GaussianBasisSet>& bs, Ref<GaussianBasisSet>& bs_aux) const;

    void run();
};

}

#endif

// //////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
