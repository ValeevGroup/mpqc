//
// r12ia.h
//
// Copyright (C) 2002 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_r12ia_h
#define _chemistry_qc_mbptr12_r12ia_h

#ifdef __GNUC__
#pragma interface
#endif

#include <vector>
#include <util/ref/ref.h>
#include <util/state/state.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/group/memory.h>
#include <chemistry/qc/basis/tbint.h>

using namespace std;

namespace sc {

/////////////////////////////////////////////////////////////////
/**R12IntsAcc contains a set of one or more distributed dense 4-index arrays.
   These 4-index quantities are typically AO or MO
   integrals. Since integrals are typically computed in sets -- for example, R12 methods
   require integrals of operators 1/r12, r12, [t1,r12], and [t2,r12] at the same time -- several sets of
   integrals are held together. The data can then be described as (ijxy)^O
   where i, j, x, and, y are indices with ranges [0,ni), [0,nj), etc. and O labels
   the operator type and ranges from 0 to num_te_types. The data is stored and accessed as follows:
   each ij block is a set of num_te_types base-0 contiguous 2-dimensional array with dimensions
   nx and ny. How blocks are stored and accessed is determined in the derived class.

   Public interface of R12IntsAcc is designed to accomodate the needs of the TwoBodyMOIntsTransform
   objects. Parallel AO->MO integral transforms are performed in single or multiple passes.
   In the latter case all (jxy)^O integrals are produces for a particular subrange of i.
   These integrals are contained in a MemoryGrp object. The contents of the MemoryGrp object is
   "stored" in accumulator using <tt>store_memorygrp(Ref<MemoryGrp>& mem, int ni)</tt>, where
   ni is the size of the subrange of i produced in this pass.
   After all batches have been stored, the content of R12IntsAcc needs to be "committed"
   using <tt>commit</tt>. After that blocks of MO integrals can be accessed using
   <tt>retrieve_pair_block</tt>.
    */

class R12IntsAcc: virtual public SavableState {
  public:
    // mem will be used to fetch data
    R12IntsAcc(int num_te_types, int ni, int nj, int nx, int ny);
    R12IntsAcc(StateIn&);
    virtual ~R12IntsAcc();
    void save_data_state(StateOut&);

    /// Types of two-body operators that R12IntsAcc understands
    typedef unsigned int tbint_type;
    static const unsigned int max_num_te_types_ = TwoBodyInt::max_num_tbint_types;

    /// The number of types of integrals that are being handled together
    int num_te_types() const { return num_te_types_; };
    /// Rank of index space i
    int ni() const { return ni_; }
    /// Rank of index space j
    int nj() const { return nj_; }
    /// Rank of index space x
    int nx() const { return nx_; }
    /// Rank of index space y
    int ny() const { return ny_; }
    /// Size of each block of the integrals of one type, in double words
    size_t blocksize() const { return blocksize_; };
    // return blksize_
    size_t blksize() const { return blksize_; }

    /// call this before operations on this object can begin
    virtual void activate() { active_ = true; }
    /// call this after operations on this object are finished. May destroy data (see data_persistent()).
    virtual void deactivate() { active_ = false; }
    /// if this returns false, call to deactivate may destroy data
    virtual bool data_persistent() const =0;
    /// Retrieves an ij block of integrals
    virtual const double * retrieve_pair_block(int i, int j, tbint_type oper_type) const =0;
    /// Releases the buffer that holds ij block of integrals
    virtual void release_pair_block(int i, int j, tbint_type oper_type) const =0;
    /// Stores an ij pair block of integrals
    virtual void store_pair_block(int i, int j, tbint_type oper_type, const double* ints) =0;

    int ij_index(int i, int j) const { return i*nj_ + j; };

    /// Can this block be accessed via retrieve_pair_block from this task?
    virtual bool is_avail(int i, int j) const =0;
    /// Is this block stored locally? If true, this implies that it can be retrieved efficiently (compare to is_avail).
    virtual bool is_local(int i, int j) const =0;
    /// Does this task have access to all blocks?
    virtual bool has_access(int proc) const =0;
    /** Returns the total number of tasks with access to integrals.
        If task i has access to the integrals, then twa_map[i] is its index among
        the tasks with access, -1 otherwise. */
    int tasks_with_access(vector<int>& twa_map) const;

  private:
    /// Set to nonzero to debug this and derived classes
    static const int classdebug_ = 0;

    int num_te_types_;  // Number of types of integrals in a block
    Ref<MessageGrp> msg_;
    int ni_, nj_;
    int nx_, ny_;
    size_t nxy_;        // nx_ * ny_  - the number of integrals of one type in a block
    size_t blksize_;    // the same in bytes
    size_t blocksize_;  // hence the size of the block of num_te_types of integrals is blksize_ * num_te_types

    bool active_;

   protected:
    // return nxy_
    size_t nxy() const { return nxy_; }
    // return active_
    bool active() const { return active_; }
    /// total number of tasks
    int ntasks() const { return msg_->n(); }
    /// rank of this task
    int me() const { return msg_->me(); }
    /// return debug level for this class
    int classdebug() const { return classdebug_; }

};


namespace detail {

  /** Stores all pair block of integrals held in mem
   in a layout assumed throughout R12IntEval.
   Let's suppose the number of tasks is nproc, nj is the number of j indices,
   ni is the number of i indices of integrals held in
   mem at the moment. Then all integrals with a given i and j
   are stored on task (i*nj+j)/nproc and this ij block is
   (i*nj+j)%nproc -th block on this task. Each ij block contains
   num_te_types_ subblocks of integrals. Each subblock of integrals
   has blksize_memgrp bytes allocated for it. Note that
   blksize_memgrp may be larger than blksize_ because an ij-block of partially
   transformed integrals may be larger than the block of fully transformed integrals.
   */
  void store_memorygrp(Ref<R12IntsAcc>& acc, Ref<MemoryGrp>& mem, int i_offset,
                       int ni, const size_t blksize_memgrp = 0);

  /** Reverse of store_memorygrp()
   */
  void restore_memorygrp(Ref<R12IntsAcc>& acc, Ref<MemoryGrp>& mem, int i_offset,
                         int ni, const size_t blksize_memgrp = 0);

}

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
