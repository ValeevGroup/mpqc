//
// r12ia_memgrp.h
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

    Ref<MemoryGrp> mem_;     /// The MemoryGrp that holds integrals
    /// The size of the ij-block in held memory, may be larger than that returned by blksize().
    size_t blksize_memgrp_;

    /** Holds the information and data for ij-blocks.
        Local blocks are held in MemoryGrp whereas remote blocks are held in local storage. */
    struct PairBlkInfo {
      const double *ints_[max_num_te_types_];   // blocks corresponding to each operator type
      mutable int refcount_[max_num_te_types_];         // counts references
      distsize_t offset_;                       // global Memgrp offset in bytes
    } *pairblk_;

    /// Initialization tasks common to all constructors
    void init();
    // Utility functions
    int ij_proc(int i, int j) const { return ij_index(i,j)%ntasks();};

  public:
    R12IntsAcc_MemoryGrp(const Ref<MemoryGrp>&, int num_te_types,
                         int ni, int nj, int nx, int ny,
                         size_t memorygrp_blksize);
    R12IntsAcc_MemoryGrp(StateIn&);
    ~R12IntsAcc_MemoryGrp();
    void save_data_state(StateOut&);

    Ref<R12IntsAcc> clone(const R12IntsAccDimensions& dim = R12IntsAccDimensions::default_dim());

    // Implementation of R12IntsAcc_MemoryGrp
    void deactivate();
    /// implementation of R12IntsAcc::data_persistent()
    bool data_persistent() const { return false; }
    /// Stores an ij pair block of integrals (assumes the block resides locally)
    void store_pair_block(int i, int j, tbint_type oper_type, const double *ints);
    /// Retrieves an ij pair block of integrals
    const double* retrieve_pair_block(int i, int j, tbint_type oper_type) const;
    /// Releases an ij pair block of integrals (if needed)
    void release_pair_block(int i, int j, tbint_type oper_type) const;

    /// Is this block stored locally?
    bool is_local(int i, int j) const { return (ij_proc(i,j) == mem_->me());};
    /// In this implementation all blocks are globally available
    bool is_avail(int i, int j) const { return true;};
    /// Does this task have access to all the integrals?
    bool has_access(int proc) const { return true;};

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
