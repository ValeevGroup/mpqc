//
// r12ia.h
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

#ifndef _chemistry_qc_mbptr12_r12ia_h
#define _chemistry_qc_mbptr12_r12ia_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/ref/ref.h>
#include <util/state/state.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/group/memory.h>

namespace sc {

/////////////////////////////////////////////////////////////////
/** R12IntsAcc accumulates transformed (MO) integrals of (ix|jy) type
   where i and j are in the active occupied MO space O and x and y are in any
   spaces X and Y, respectively. Transformed integrals are usually computed
   using a parallel MO integrals transformation procedure. In general, such
   transformations will require multiple passes through AO integrals. Each pass
   produces a batch of transformed integrals. A batch is a set of integrals {(ix|jy)}
   in which i indices are in a finite subrange of O
   and x, j, and y take any of their allowed values. For example, if batch I contains
   all integrals (ix|jy) with i greater than or equal m but less than n, then batch I+1
   contains integrals (ix|jy) with i greater than n. Integrals in batch 0 have indices
   i greater than or equal to 0.
   
   After each pass the MO integrals are contained in a MemoryGrp object. The object is
   "stored" in accumulator using <tt>store_memorygrp(Ref<MemoryGrp>& mem, int ni)</tt>.
   After all batches have been stored, the content of R12IntsAcc needs to be "committed"
   using <tt>commit()</tt>. After that blocks of MO integrals can be accessed using
   <tt>retrieve_pair_block</tt>.
    */

class R12IntsAcc: virtual public SavableState {

    int num_te_types_;  // Number of types of integrals in a block (in R12 theories -- usually 3)

   protected:
    int nocc_act_;
    int nbasis1_, nbasis2_;
    size_t nbasis__2_;  // nbasis1_ * nbasis2_  - the size of a block of integrals of one type
    size_t blksize_;    // the same in bytes
    size_t blocksize_;  // hence the size of the block of num_te_types of integrals is nbasis__2_ * num_te_types
    
    int next_orbital_;  // The first index of the next batch to be stored
    bool committed_;    // Whether all data has been written out and ready to be read

    /// The index of the first orbital in the next integrals batch to be stored
    void inc_next_orbital(int ni);

  public:
    R12IntsAcc(int num_te_types, int nbasis1, int nbasis2, int nocc_act);
    R12IntsAcc(StateIn&);
    ~R12IntsAcc();
    void save_data_state(StateOut&);

    /// Types of two-body operators that R12IntsAcc understands
    enum tbint_type { eri=0, r12=1, r12t1=2, r12t2=3};
    static const int max_num_te_types_ = 4;

    /// The number of types of integrals that are being handled together
    int num_te_types() const { return num_te_types_; };
    /// Size of each block of the integrals of the same type, in double words
    size_t blocksize() const { return nbasis__2_; };
    /// The index of the first orbital in the next integrals batch to be stored
    int next_orbital() const;

    /// Stores all pair block of integrals held in mem
    /// in a layout assumed throughout MBPT2_R12
    /// ni is the range of index i of integrals held in mem at the moment
    virtual void store_memorygrp(Ref<MemoryGrp>& mem, int ni)=0;
    /// All member functions of this class and its children
    /// indices i and j don't include frozen orbitals
    /// Stores an ij pair block of integrals (assumes the block resides locally)
    virtual void store_pair_block(int i, int j, double *ints)=0;
    /// Commit the content of the accumulator for reading
    virtual void commit() { committed_ = true; };
    /// Has the content of the accumulator been commited for reading?
    bool is_commited() { return committed_; };
    /// Call when done reading content
    virtual void deactivate() {};
    /// Retrieves an ij pair block of integrals
    virtual double* retrieve_pair_block(int i, int j, tbint_type oper_type) =0;
    /// Releases an ij pair block of integrals (if needed)
    virtual void release_pair_block(int i, int j, tbint_type oper_type) =0;
    /// Is this block stored locally?
    virtual bool is_local(int i, int j) const =0;
    /// Is this block available to this task?
    virtual bool is_avail(int i, int j) const =0;
    /// Does this task have access to all the integrals?
    virtual bool has_access(int proc) const =0;
    /// Can this specialization be used in restarts?
    virtual bool can_restart() const =0;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
