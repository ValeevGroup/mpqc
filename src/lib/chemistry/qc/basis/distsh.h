//
// distsh.h
// modeled after distshpair.h
//
// Copyright (C) 2009 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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

#ifndef _chemistry_qc_basis_distsh_h
#define _chemistry_qc_basis_distsh_h

#include <util/misc/regtime.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/basis.h>

namespace sc {

/// Distributes sets of shells either statically or dynamically.
class DistShell {
  public:
    /** This is used to store data that must be shared between all
        cooperating shell sets.
     */
     class SharedData {
       public:
         volatile long int shell_;
         /// Construct and initialize.
         SharedData() { init(); }
         /** If this will be used to iterate through the shells again,
          * then init must be called.
          */
         void init() { shell_ = 0; }
     };

    /** The DistShell class is used to distribute sets of shells to compute tasks.

        Both static (round-robin) and dynamic methods are supported. */
    DistShell(const Ref<MessageGrp> &, int nthread, int mythread,
              const Ref<ThreadLock>& lock,
              const Ref<GaussianBasisSet>& bs,
              bool dynamic,
              int max_nfunctions = -1,
              int max_nshell = 1,
              SharedData *shared = 0);
    ~DistShell();
    /// Resets to the first shell.
    void init();
    /// How much stuff to print out.
    void set_debug(int d) { debug_ = d; }
    /** How often to print status from node 0.  If p > 100.0, then
        no printing will be done. */
    void set_print_percent(double p);
    /** Copies "payload" (I is the first shell of the shell set, in the increasing-size order,
        and N is the number of shells)  and returns 1.
        When there are no more shells to be processed by this processor,
        0 is returned. */
    int get_task(int& I, int& N);

    /** maps shell index from the work ordering (as reported to get_task() )
        to its index within the basis set */
    int shell_index(int i) const;

  private:
    Ref<MessageGrp> msg_;
    int nthread_;
    Ref<ThreadLock> lock_;
    Ref<GaussianBasisSet> bs_;
    int max_nfunction_;
    int max_nshell_;
    bool task_dynamic_;
    bool thread_dynamic_;      //< only true is 1 task, multiple threads, and shared data given
    int debug_;
    // How often updates are printed (i.e. every 10% of total work)
    double print_percent_;
    SharedData *shared_;

    typedef std::vector<int> ShellIndexMap;
    ShellIndexMap shell_map_;

    typedef std::pair<int,int> Task;
    typedef std::vector<Task> Tasks;
    Tasks tasks_;

    // Number of tasks handled by thread 0 in task 0:
    // if dynamic == true : it will distribute all of them
    // if dynamic == false : it will handle its share
    long int ntask_;
    // Print period
    long int print_interval_;
    // Index of the next task to be served
    int current_task_;

    // for dynamic load balancing
    int req_type_;
    int ans_type_;
    void serve_tasks();

    // for static load balancing
    int ncpu_;
    int mythread_;

    // computes work units subject to max_nfunction and max_nshell constraints
    void init_work();

  public:
};

}

#endif

// //////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
