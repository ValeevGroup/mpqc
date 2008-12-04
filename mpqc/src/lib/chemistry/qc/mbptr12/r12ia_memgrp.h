//
// r12ia_memgrp.h
//
// Copyright (C) 2002 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_r12ia_memgrp_h
#define _chemistry_qc_mbptr12_r12ia_memgrp_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/ref/ref.h>
#include <util/group/memory.h>
#include <chemistry/qc/mbptr12/r12ia.h>

namespace sc {

/////////////////////////////////////////////////////////////////////
/** R12IntsAcc_MemoryGrp handles transformed integrals held in memory
    by MemoryGrp

    The ordering of integrals in MemoryGrp buffers is not specified
    to avoid having to reorder integrals
    Each pair block has size of num_te_types*nbasis1*nbasis2
*/

class R12IntsAcc_MemoryGrp: public R12IntsAcc {

    Ref<MemoryGrp> mem_; /// The MemoryGrp used by this accumulator to store integrals
    int nproc_;
    size_t blksize_memgrp_;  // The size of the ij-block in held memory (may be larger than blksize_)

    struct PairBlkInfo {
      double *ints_[max_num_te_types_];         // blocks corresponding to each operator type
      int refcount_[max_num_te_types_];          // number of references
      distsize_t offset_;    // global Memgrp offset in bytes
    } *pairblk_;
    
    /// Initialization tasks common to all constructors
    void init();
    /// total number of tasks
    int ntasks() const { return mem_->n(); }
    /// ID of this task
    int taskid() const { return mem_->me(); }

    /// Stores an ij pair block of integrals (assumes the block resides locally)
    void store_pair_block(int i, int j, double *ints);
    
  public:
    R12IntsAcc_MemoryGrp(Ref<MemoryGrp>&, int num_te_types, int ni, int nj, int nx, int ny);
    R12IntsAcc_MemoryGrp(StateIn&);
    ~R12IntsAcc_MemoryGrp();
    void save_data_state(StateOut&);

    /** Implements R12IntsAcc::store_memorygrp().
     mem must be the same Memorygrp used to construct this.
     This is a collective operation.
     */
    void store_memorygrp(Ref<MemoryGrp>& mem, int ni, const size_t blksize = 0);
    /// Implements R12IntsAcc::restore_memorygrp(). mem must be the same MemoryGrp used to construct this.
    void restore_memorygrp(Ref<MemoryGrp>& mem, int ioffset, int ni, const size_t blksize = 0) const;
    /// Done reading content - call set_localsize(0) on the associated MemoryGrp
    /// This is a collective operation. This accumulator cannot be activated again.
    void deactivate();
    /// Retrieves an ij pair block of integrals
    double* retrieve_pair_block(int i, int j, tbint_type oper_type);
    /// Releases an ij pair block of integrals (if needed)
    void release_pair_block(int i, int j, tbint_type oper_type);
    /// Is this block stored locally?
    bool is_local(int i, int j) const { return (ij_proc(i,j) == mem_->me());};
    /// In this implementation all blocks are globally available
    bool is_avail(int i, int j) const { return true;};
    /// Does this task have access to all the integrals?
    bool has_access(int proc) const { return true;};
    /// Cannot restart MemoryGrp-based accumulator
    bool can_restart() const { return false; };

    // Utility functions
    int ij_index(int i, int j) const { return i*nj_ + j; };
    int ij_proc(int i, int j) const { return ij_index(i,j)%nproc_;};
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
