//
// distshpair.h
// based on mbpt/distsh.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
// Updated: Edward Valeev <evaleev@vt.edu>
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

#ifndef _chemistry_qc_basis_distshpair_h
#define _chemistry_qc_basis_distshpair_h

#include <util/misc/regtime.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/basis.h>

namespace sc {

/// Distributes shell pairs either statically or dynamically.
class DistShellPair {
  public:
    /** This is used to store data that must be shared between all
     * cooperating shell pairs.
     */
     class SharedData {
       public:
         volatile long int shellpair_;
         /// Construct and initialize.
         SharedData() { init(); }
         /** If this will be used to iterate through the shells again,
          * then init must be called.
          */
         void init() { shellpair_ = 0; }
     };
  private:
    Ref<MessageGrp> msg_;
    int nthread_;
    Ref<ThreadLock> lock_;
    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    bool bs1_eq_bs2_;
    bool task_dynamic_;
    bool thread_dynamic_;
    int debug_;
    // How often updates are printed (i.e. every 10% of total work)
    double print_percent_;
    SharedData *shared_;

    // Number of tasks handled by thread 0 in task 0:
    // if dynamic == true : it will distribute all of them
    // if dynamic == false : it will handle its share
    long int ntask_;
    // Print period
    long int print_interval_;
    // Index of the next task to be served
    long int current_shellpair_;

    // for dynamic load balancing
    int req_type_;
    int ans_type_;
    int ncpu_less_0_;
    void serve_tasks();

    // for static load balancing
    int S_, R_;          // NOTE: S is in bs1, R is in bs2
    int ncpu_;
    int incS_, incR_;
    int mythread_;

    // sorted work for dynamic load balancing
    int *cost_;
    int *Svec_;
    int *Rvec_;
    int *Ivec_;

    void init_dynamic_work();
  public:
    /** The DistShellPair class is used to distribute shell pair indices among tasks.

        Both static (round-robin) and dynamic methods are supported. */
    DistShellPair(const Ref<MessageGrp> &, int nthread, int mythread,
                  const Ref<ThreadLock>& lock,
                  const Ref<GaussianBasisSet>& bs1, const Ref<GaussianBasisSet>& bs2,
		  bool dynamic, SharedData *shared = 0);
    ~DistShellPair();
    /// Resets to the first shell pair.
    void init();
    /// How much stuff to print out.
    void set_debug(int d) { debug_ = d; }
    /** How often to print status from node 0.  If p > 100.0, then
        no printing will be done. */
    void set_print_percent(double p);
    /** Puts the current PQ shell pair into P and Q and returns 1.
        When there are no more shell pairs to be processed by this processor,
        0 is returned.  Once we start doing get_tasks, we have to go to the
        end if dynamic load balancing is used.

        P belongs to bs1, and Q belongs to bs2. If (bs1 == bs2) then
        P is greater or equal to Q. */
    int get_task(int &P, int &Q);
};

}

#endif

// //////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
