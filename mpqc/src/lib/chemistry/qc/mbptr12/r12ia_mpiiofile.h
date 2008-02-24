//
// r12ia_mpiiofile.h
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

#ifndef _chemistry_qc_mbptr12_r12ia_mpiiofile_h
#define _chemistry_qc_mbptr12_r12ia_mpiiofile_h

#ifdef __GNUC__
#pragma interface
#endif

#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include <util/ref/ref.h>
#include <util/group/memory.h>
#include <chemistry/qc/mbptr12/r12ia.h>

namespace sc {

//////////////////////////////////////////////////////////////////////////
/** R12IntsAcc_MPIIOFile handles transformed integrals stored in a binary
    file accessed through MPI-IO. This is an abstract base for MPIIO-based
    accumulators using individual and collective I/O.

    The ordering of integrals in blocks is not specified
    to avoid having to reorder integrals
    Each pair block has size of num_te_types*nbasis1*nbasis2
*/

class R12IntsAcc_MPIIOFile: public R12IntsAcc {

  protected:
    Ref<MemoryGrp> mem_; /// The MemoryGrp associated with this accumulator
    int nproc_;
    size_t nints_per_block_;  // number of integrals per block = num_te_types*nbasis__2_
    char *filename_;
    MPI_File datafile_;

    struct PairBlkInfo {
      double* ints_[max_num_te_types_];      // blocks corresponding to each operator type
      int refcount_[max_num_te_types_];      // number of references
      MPI_Offset offset_;      // location in file (in bytes)
    } *pairblk_;

    /// Utility function to check MPI I/O error codes.
    void check_error_code_(int errcod) const;
    /// Initialization tasks common to all constructors
    void init(bool restart);
    
    /// total number of tasks
    int ntasks() const { return mem_->n(); }
    /// ID of this task
    int taskid() const { return mem_->me(); }
    
  public:
    R12IntsAcc_MPIIOFile(Ref<MemoryGrp>& mem, const char *filename, int num_te_types,
                         int ni, int nj, int nx, int ny);
    R12IntsAcc_MPIIOFile(StateIn&);
    ~R12IntsAcc_MPIIOFile();
    void save_data_state(StateOut&);

    /// Stores an ij pair block of integrals to the file
    void store_pair_block(int i, int j, double *ints);
    /** Commit the content of the accumulator for reading - deactivate the associated MemoryGrp
        and activate the accumulator. This is a collective operation. */
    void commit();
    /// Call when ready to read content
    void activate();
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
    /// Can restart MPI-IO-based accumulator
    bool can_restart() const { return true; };

    // Utility functions
    int ij_index(int i, int j) const { return i*nj_ + j; };
};

//////////////////////////////////////////////////////////////////////////////
/** R12IntsAcc_MPIIOFile_Ind handles transformed integrals stored in a binary
    file accessed through MPI-IO individual I/O routines.
    
    The ordering of integrals in blocks is not specified
    to avoid having to reorder integrals
    Each pair block has size of num_te_types*nbasis*nbasis
*/

class R12IntsAcc_MPIIOFile_Ind: public R12IntsAcc_MPIIOFile {

  public:
    R12IntsAcc_MPIIOFile_Ind(Ref<MemoryGrp>& mem, const char *filename, int num_te_types,
                             int ni, int nj, int nx, int ny);
    R12IntsAcc_MPIIOFile_Ind(StateIn&);
    ~R12IntsAcc_MPIIOFile_Ind();
    void save_data_state(StateOut&);

    /** Stores all pair block of integrals held in mem.
        By default blocks are appended to the end of the same file, i.e.
        they are assumed to have come from consecutive passes of
        the same transformation.
        This is a collective operation.
        See documentation for R12IntsAcc::store_memorygrp()
        for more info.
      */
    void store_memorygrp(Ref<MemoryGrp>& mem, int ni, const size_t blksize = 0);
    /// Retrieves an ij pair block of integrals from the file
    double* retrieve_pair_block(int i, int j, tbint_type oper_type);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
