//
// r12intsacc.h
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

#ifndef _chemistry_qc_mbptr12_r12intsacc_h
#define _chemistry_qc_mbptr12_r12intsacc_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/ref/ref.h>
#include <util/group/memory.h>
#include <unistd.h>

#include <mpi.h>

namespace sc {

/////////////////////////////////////////////////////////////////
// R12IntsAcc accumulates transformed integrals of (og|og) type
//   o = active occupied MO
//   g = general MO

class R12IntsAcc: public RefCount {

   protected:
    int nocc_, nocc_act_, nfzc_;
    int nbasis_;
    size_t nbasis__2_;  // nbasis_ squared - the size of a block of integrals of one type
    size_t blksize_;    // the same in bytes
    size_t blocksize_;  // hence the size of the block of num_te_types of integrals is nbasis__2_ * num_te_types
    int noso_;
    bool committed_;    // Whether all data has been written out and ready to be read

#define num_te_types_ 3

  public:
    R12IntsAcc(int nbasis, int noso, int nocc, int nfzc);
    ~R12IntsAcc();

    /// Types of two-body operators that r12IntsAcc understands
    enum tbint_type { eri=0, r12=1, r12t1=2};
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
    /// Done reading content
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
};

/////////////////////////////////////////////////////////////////////
// R12IntsAcc_MemoryGrp handles transformed integrals held in memory
// by MemoryGrp
//
// The ordering of integrals in MemoryGrp buffers is not specified
// to avoid having to reorder integrals
// Each pair block has size of num_te_types*nbasis*nbasis

class R12IntsAcc_MemoryGrp: public R12IntsAcc {

    Ref<MemoryGrp> mem_; // The MemoryGrp used by this accumulator to store integrals
    int nproc_;

    struct PairBlkInfo {
      double *ints_[num_te_types_];         // blocks corresponding to each operator type
      distsize_t offset_;    // global Memgrp offset in bytes
    } *pairblk_;

  public:
    R12IntsAcc_MemoryGrp(Ref<MemoryGrp>&, int nbasis, int noso, int nocc, int nfzc);
    ~R12IntsAcc_MemoryGrp();

    /// Stores all pair block of integrals held in mem
    /// mem must be the same as mem_ used to construct this
    /// This is a collective operation
    void store_memorygrp(Ref<MemoryGrp>& mem, int ni);
    /// Stores an ij pair block of integrals (assumes the block resides locally)
    void store_pair_block(int i, int j, double *ints);
    /// Done reading content - call set_localsize(0) on the associated MemoryGrp
    /// This is a collective operation
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

    // Utility functions
    int ij_index(int i, int j) const { return i*nocc_act_ + j; };
    int ij_proc(int i, int j) const { return ij_index(i,j)%nproc_;};
};

/////////////////////////////////////////////////////////////////////
// R12IntsAcc_Node0File handles transformed integrals stored in file
// on node 0 (file is a usual POSIX binary file)
//
// Transfering integrals to the file from nodes is done via MemoryGrp
//   given as an argument to store_memorygrp
// Remote retrieval is not possible
//
// The ordering of integrals in blocks is not specified
// to avoid having to reorder integrals
// Each pair block has size of num_te_types*nbasis*nbasis

class R12IntsAcc_Node0File: public R12IntsAcc {

    Ref<MemoryGrp> mem_; // The MemoryGrp associated with this accumulator
    int icounter_;      // the number of i's for which blocks have been written to 
    char *filename_;
    int datafile_;

    struct PairBlkInfo {
      double* ints_[num_te_types_];      // blocks corresponding to each operator type
      off_t offset_;      // location in file (in bytes)
    } *pairblk_;
    
  public:
    R12IntsAcc_Node0File(Ref<MemoryGrp>& mem, char *filename, int nbasis, int noso, int nocc, int nfzc, bool restart);
    ~R12IntsAcc_Node0File();

    /// Stores all pair block of integrals held in mem.
    /// By default blocks are appended to the end of the same file, i.e.
    /// they are assumed to have come from consecutive passes of
    /// the same transformation
    /// This is a collective operation
    void store_memorygrp(Ref<MemoryGrp>& mem, int ni);
    /// Stores an ij pair block of integrals to the file
    void store_pair_block(int i, int j, double *ints);
    /// Commit the content of the accumulator for reading - deactivate the associated MemoryGrp
    /// This is a collective operation
    void commit();
    /// Done reading content - activate the associated MemoryGrp
    /// This is a collective operation
    void deactivate();
    /// Retrieves an ij pair block of integrals from the file
    double* retrieve_pair_block(int i, int j, tbint_type oper_type);
    /// Releases an ij pair block of integrals
    void release_pair_block(int i, int j, tbint_type oper_type);
    /// Is this block stored locally?
    bool is_local(int i, int j) const { return (mem_->me() == 0);};
    /// In this implementation blocks are available only on node 0
    bool is_avail(int i, int j) const { return (mem_->me() == 0);};
    /// Does this task have access to all the integrals?
    bool has_access(int proc) const { return (proc == 0);};

    // Utility functions
    int ij_index(int i, int j) const { return i*nocc_act_ + j; };
    int ij_proc(int i, int j) const { return 0;};
};

//////////////////////////////////////////////////////////////////////////
// R12IntsAcc_MPIIOFile handles transformed integrals stored in a binary
// file accessed through MPI-IO. This is an abstract base for MPIIO-based
// accumulators using individual and collective I/O.
//
// The ordering of integrals in blocks is not specified
// to avoid having to reorder integrals
// Each pair block has size of num_te_types*nbasis*nbasis

class R12IntsAcc_MPIIOFile: public R12IntsAcc {

  protected:
    Ref<MemoryGrp> mem_; // The MemoryGrp associated with this accumulator
    int nproc_;
    size_t nints_per_block_;  // number of integrals per block = num_te_types*nbasis__2_
    int icounter_;      // the number of i's for which blocks have been written to 
    char *filename_;
    MPI_File datafile_;

    struct PairBlkInfo {
      double* ints_[num_te_types_];      // blocks corresponding to each operator type
      MPI_Offset offset_;      // location in file (in bytes)
    } *pairblk_;
    
  public:
    R12IntsAcc_MPIIOFile(Ref<MemoryGrp>& mem, char *filename, int nbasis, int noso, int nocc, int nfzc, bool restart);
    ~R12IntsAcc_MPIIOFile();

    /// Stores an ij pair block of integrals to the file
    void store_pair_block(int i, int j, double *ints);
    /// Commit the content of the accumulator for reading - deactivate the associated MemoryGrp
    /// This is a collective operation
    void commit();
    /// Done reading content - activate the associated MemoryGrp
    /// This is a collective operation
    void deactivate();
    /// Releases an ij pair block of integrals
    void release_pair_block(int i, int j, tbint_type oper_type);
    /// Is this block stored locally?
    bool is_local(int i, int j) const { return true;};
    /// In this implementation blocks are available everywhere
    bool is_avail(int i, int j) const { return true;};
    /// Does this task have access to all the integrals?
    bool has_access(int proc) const { return true;};

    // Utility functions
    int ij_index(int i, int j) const { return i*nocc_act_ + j; };
};

//////////////////////////////////////////////////////////////////////////////
// R12IntsAcc_MPIIOFile_Ind handles transformed integrals stored in a binary
// file accessed through MPI-IO individual I/O routines.
//
// The ordering of integrals in blocks is not specified
// to avoid having to reorder integrals
// Each pair block has size of num_te_types*nbasis*nbasis

class R12IntsAcc_MPIIOFile_Ind: public R12IntsAcc_MPIIOFile {

  public:
    R12IntsAcc_MPIIOFile_Ind(Ref<MemoryGrp>& mem, char *filename, int nbasis, int noso, int nocc, int nfzc, bool restart);
    ~R12IntsAcc_MPIIOFile_Ind();

    /// Stores all pair block of integrals held in mem.
    /// By default blocks are appended to the end of the same file, i.e.
    /// they are assumed to have come from consecutive passes of
    /// the same transformation
    /// This is a collective operation
    void store_memorygrp(Ref<MemoryGrp>& mem, int ni);
    /// Retrieves an ij pair block of integrals from the file
    double* retrieve_pair_block(int i, int j, tbint_type oper_type);
};

///////////////////////////////////////////////////////////////////////////////
// R12IntsAcc_MPIIOFile_Coll handles transformed integrals stored in a binary
// file accessed through MPI-IO collective I/O routines.
//
// The ordering of integrals in blocks is not specified
// to avoid having to reorder integrals
// Each pair block has size of num_te_types*nbasis*nbasis

class R12IntsAcc_MPIIOFile_Coll: public R12IntsAcc_MPIIOFile {

    MPI_Datatype filetype_;
    
  public:
    R12IntsAcc_MPIIOFile_Coll(Ref<MemoryGrp>& mem, char *filename, int nbasis, int noso, int nocc, int nfzc, bool restart);
    ~R12IntsAcc_MPIIOFile_Coll();

    /// Stores all pair block of integrals held in mem.
    /// By default blocks are appended to the end of the same file, i.e.
    /// they are assumed to have come from consecutive passes of
    /// the same transformation
    /// This is a collective operation
    void store_memorygrp(Ref<MemoryGrp>& mem, int ni);
    /// Retrieves an ij pair block of integrals from the file
    double* retrieve_pair_block(int i, int j, tbint_type oper_type);
};

}

#endif
