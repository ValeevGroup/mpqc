//
// r12ia_node0file.h
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

#ifndef _chemistry_qc_mbptr12_r12ia_node0file_h
#define _chemistry_qc_mbptr12_r12ia_node0file_h

#ifdef __GNUC__
#pragma interface
#endif

#include <unistd.h>
#include <util/ref/ref.h>
#include <util/group/memory.h>
#include <chemistry/qc/mbptr12/r12ia.h>

namespace sc {

/////////////////////////////////////////////////////////////////////
/** R12IntsAcc_Node0File handles transformed integrals stored in file
    on node 0 (file is a usual POSIX binary file)

    Transfering integrals to the file from nodes is done via MemoryGrp
    given as an argument to store_memorygrp
    Remote retrieval is not possible

    The ordering of integrals in blocks is not specified
    to avoid having to reorder integrals
    Each pair block has size of num_te_types*nbasis1*nbasis2
*/

class R12IntsAcc_Node0File: public R12IntsAcc {

    char *filename_;
    int datafile_;

    struct PairBlkInfo {
      // mutable since this data is only cached. offset is the only real data.
      mutable double* ints_[max_num_te_types_];      // blocks corresponding to each operator type
      mutable int refcount_[max_num_te_types_];      // number of references
      off_t offset_;      // location in file (in bytes)
    };
    PairBlkInfo* pairblk_;
    
    /// Initialization tasks common to all constructors
    void init(bool restart);
    // Check if the file operation went OK
    void check_filedescr_();
    // Utility functions
    int ij_proc(int i, int j) const { return 0;};
    
  public:
    R12IntsAcc_Node0File(const char *filename, int num_te_types, int ni, int nj, int nx, int ny);
    R12IntsAcc_Node0File(StateIn&);
    ~R12IntsAcc_Node0File();
    void save_data_state(StateOut&);

    /// implementation of R12IntsAcc::activate()
    void activate();
    /// implementation of R12IntsAcc::deactivate()
    void deactivate();
    /// implementation of R12IntsAcc::data_persistent()
    bool data_persistent() const { return true; }
    /// Stores an ij pair block of integrals to the file
    void store_pair_block(int i, int j, tbint_type oper_type, const double *ints);
    /// Retrieves an ij pair block of integrals from the file
    double* retrieve_pair_block(int i, int j, tbint_type oper_type) const;
    /// Releases an ij pair block of integrals
    void release_pair_block(int i, int j, tbint_type oper_type) const;
    
    /// Is this block stored locally?
    bool is_local(int i, int j) const { return (me() == 0);};
    /// In this implementation blocks are available only on node 0
    bool is_avail(int i, int j) const { return (me() == 0);};
    /// Does this task have access to all the integrals?
    bool has_access(int proc) const { return (proc == 0);};
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
