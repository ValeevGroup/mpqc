//
// hsosv1e1.h
// based on csgrade12.h
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

#ifndef _chemistry_qc_mbpt_hsosv1e1_h
#define _chemistry_qc_mbpt_hsosv1e1_h

#include <util/misc/regtime.h>
#include <util/group/memory.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/integral.h>

namespace sc {

class HSOSV1Erep1Qtr: public Thread {
  private:
    int mythread;
    int nthread;
    int me;
    int nproc;
    Ref<RegionTimer> timer;
    Ref<ThreadLock> lock;
    Ref<GaussianBasisSet> basis;
    Ref<TwoBodyInt> tbint;
    int ni,i_offset;
    double **scf_vector;
    double tol;
    int debug;

    double *trans_int1;
    double aoint_computed_;
    int nfuncmax;
    int nbasis;
    int nshell;

    int R,S,nr,ns;
  public:
    HSOSV1Erep1Qtr(int mythread_a, int nthread_a,
                   int me_a, int nproc_a,
                   const Ref<ThreadLock> &lock_a,
                   const Ref<GaussianBasisSet> &basis_a,
                   const Ref<TwoBodyInt> &tbint_a,
                   int ni_a, double **scf_vector_a,
                   double tol_a, int debug_a);
    ~HSOSV1Erep1Qtr();

    void run();

    void accum_buffer(double *buffer);
    void set_data(int R_a,int nr_a,int S_a,int ns_a,int ni_a, int ioffset_a);
    double aoint_computed() { return aoint_computed_; }
};

}

#endif

// //////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
