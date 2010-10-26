//
// csgrad34qb.h
// based on csgrad.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ibniels@ca.sandia.gov>
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

#ifndef _chemistry_qc_mbpt_csgrad34qb_h
#define _chemistry_qc_mbpt_csgrad34qb_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/misc/regtime.h>
#include <util/group/memory.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/distshpair.h>

namespace sc {

class CSGrad34Qbtr: public Thread {
  private:
    Ref<MessageGrp> msg;
    Ref<MemoryGrp> mem;
    Ref<TwoBodyInt> tbint;
    Ref<TwoBodyDerivInt> tbintder;
    Ref<GaussianBasisSet> basis;
    Ref<ThreadLock> lock;
    Ref<RegionTimer> timer;
    int mythread;
    int nthread;
    int ni;
    int nocc;
    int nfzc;
    int i_offset;
    int aointder_computed;
    int me;
    int nproc;
    double tol;
    double **scf_vector;
    int debug;
    int dynamic_;
    double print_percent_;
    int dograd;
    int natom;
    double *Lpi;
    double **ginter;
    DistShellPair::SharedData *shellpair_shared_data_;
  public:
    CSGrad34Qbtr(int mythread_a, int nthread_a,
                 int me_a, int nproc_a,
                 const Ref<MemoryGrp> &mem_a,
                 const Ref<MessageGrp> &msg_a,
                 const Ref<ThreadLock> &lock_a,
                 const Ref<GaussianBasisSet> &basis_a,
                 const Ref<TwoBodyInt> &tbint_a,
                 const Ref<TwoBodyDerivInt> &tbintder_a,
                 int nocc_a, int nfzc_a,
                 double **scf_vector_a,
                 double tol_a, int debug_a,
                 int dynamic_a, double print_percent_a,
                 DistShellPair::SharedData *shellpair_shared_data,
                 int dograd_a, int natom_a);
    ~CSGrad34Qbtr();

    void set_i_offset(int ioff) { i_offset = ioff; }
    void set_ni(int nivalue) { ni = nivalue; }
    void run();
    double *get_Lpi() { return Lpi ;}
    double **get_ginter() { return ginter ;}
    int get_aointder_computed() { return aointder_computed ;}
};

}

#endif

// //////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
