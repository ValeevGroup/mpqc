//
// csgrads2pdm.h
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

#ifndef _chemistry_qc_mbpt_csgrads2pdm_h
#define _chemistry_qc_mbpt_csgrads2pdm_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/misc/regtime.h>
#include <util/group/memory.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/integral.h>

namespace sc {

#define PRINT_BIGGEST_INTS 0

class CSGradS2PDM: public Thread {
  private:
    int mythread;
    int nthread;
    int me;
    int nproc;
    Ref<ThreadLock> lock;
    Ref<GaussianBasisSet> basis;
    Ref<TwoBodyDerivInt> tbintder;
    const double *PHF;
    const double *P2AO;
    int tol;
    int debug;
    int dynamic;

    double **ginter;
    double **hf_ginter;

    void accum_contrib(double **sum, double **contribs);
  public:
    CSGradS2PDM(int mythread_a, int nthread_a,
                int me_a, int nproc_a,
                const Ref<ThreadLock> &lock_a,
                const Ref<GaussianBasisSet> &basis_a,
                const Ref<TwoBodyDerivInt> &tbintder_a,
                const double *PHF_a, const double *P2AO_a,
                int tol_a, int debug_a, int dynamic_a);

    ~CSGradS2PDM();
    void accum_mp2_contrib(double **ginter);
    void accum_hf_contrib(double **hf_ginter);
    void run();
};

}

#endif

// //////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
