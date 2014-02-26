//
// distarray4_mpiiofile.h
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

#ifndef _chemistry_qc_mbptr12_distarray4_mpiiofile_h
#define _chemistry_qc_mbptr12_distarray4_mpiiofile_h

#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include <util/ref/ref.h>
#include <util/group/memory.h>
#include <util/misc/registry.h>
#include <math/distarray4/distarray4.h>

/// xlc 11 seems to dislike presence of a template member function that has the same name as a virtual function
/// remove after xlc is fixed
#define AVOID_XLC_BUG 1

namespace sc {

//////////////////////////////////////////////////////////////////////////
/** DistArray4_MPIIOFile handles transformed integrals stored in a binary
    file accessed through MPI-IO. This is an abstract base for MPIIO-based
    accumulators using individual and collective I/O.

    The ordering of integrals in blocks is not specified
    to avoid having to reorder integrals
    Each pair block has size of num_te_types*nbasis1*nbasis2
*/

class DistArray4_MPIIOFile: public DistArray4 {

  protected:

    size_t nints_per_block_;  // number of integrals per block = num_te_types*nbasis__2_
    char *filename_;
    MPI_File datafile_;

    // keep track of clones of this object to be able to create unique names
    typedef Registry<std::string,int,detail::NonsingletonCreationPolicy> ListOfClones;
    Ref<ListOfClones> clonelist_;
    void set_clonelist(const Ref<ListOfClones>& cl);

    struct PairBlkInfo {
      mutable double* ints_[max_num_te_types];      // blocks corresponding to each operator type
      mutable bool manage_[max_num_te_types];       // am I managing memory of this ints_?
      mutable int refcount_[max_num_te_types];      // number of references
      MPI_Offset offset_;      // location in file (in bytes)
    } *pairblk_;

    /// Utility function to check MPI I/O error codes.
    void check_error_code_(int errcod) const;
    /// Initialization tasks common to all constructors
    void init(bool restart);
    // Utility functions
    int ij_proc(int i, int j) const { return 0;};

#ifndef AVOID_XLC_BUG
    /// helps to implement clone of Derived class
    template <typename Derived> Ref<DistArray4> clone(const DistArray4Dimensions& dim =
                                                      DistArray4Dimensions::default_dim());
#endif

  public:
    DistArray4_MPIIOFile(const char *filename, int num_te_types,
                         int ni, int nj, int nx, int ny,
                         DistArray4Storage storage = DistArray4Storage_XY);
    DistArray4_MPIIOFile(StateIn&);
    ~DistArray4_MPIIOFile();
    void save_data_state(StateOut&);

     /// implementation of DistArray4::activate()
    void activate();
    /// implementation of DistArray4::deactivate()
    void deactivate();
    /// implementation of DistArray4::data_persistent()
    bool data_persistent() const { return true; }
    /// Releases an ij pair block of integrals
    void release_pair_block(int i, int j, tbint_type oper_type) const;

    /// Is this block stored locally?
    bool is_local(int i, int j) const { return true; }
    /// In this implementation blocks are available everywhere
    bool is_avail(int i, int j) const { return true; }
    /// Does this task have access to all the integrals?
    bool has_access(int proc) const { return true; }
};

namespace detail {
  void clone_filename(std::string& result, const char* original, int id);
}

#ifndef AVOID_XLC_BUG
template <typename Derived>
  Ref<DistArray4> DistArray4_MPIIOFile::clone(const DistArray4Dimensions& dim) {

    int id = 0;
    std::string clonename;
    using detail::clone_filename;
    clone_filename(clonename, this->filename_, id);
    if (clonelist_) {
      while (clonelist_->key_exists(clonename)) {
        ++id;
        clone_filename(clonename, this->filename_, id);
      }
    } else {
      clonelist_ = ListOfClones::instance();
    }
    clonelist_->add(clonename, id);

    Ref<Derived> result;
    if (dim == DistArray4Dimensions::default_dim())
      result = new Derived(clonename.c_str(), num_te_types(),
                           ni(), nj(),
                           nx(), ny(),
                           storage());
    else
      result = new Derived(clonename.c_str(), dim.num_te_types(),
                           dim.n1(), dim.n2(),
                           dim.n3(), dim.n4(),
                           dim.storage());

    result->set_clonelist(clonelist_);
    return result;
  }
#endif

//////////////////////////////////////////////////////////////////////////////
/** DistArray4_MPIIOFile_Ind handles transformed integrals stored in a binary
    file accessed through MPI-IO individual I/O routines.

    The ordering of integrals in blocks is not specified
    to avoid having to reorder integrals
    Each pair block has size of num_te_types*nbasis*nbasis
*/

class DistArray4_MPIIOFile_Ind: public DistArray4_MPIIOFile {

  public:
    DistArray4_MPIIOFile_Ind(const char *filename, int num_te_types,
                             int ni, int nj, int nx, int ny,
                             DistArray4Storage storage = DistArray4Storage_XY);
    DistArray4_MPIIOFile_Ind(StateIn&);
    ~DistArray4_MPIIOFile_Ind();
    void save_data_state(StateOut&);

    Ref<DistArray4> clone(const DistArray4Dimensions& dim = DistArray4Dimensions::default_dim());

    void store_pair_block(int i, int j, tbint_type oper_type, const double *ints);
    void store_pair_subblock(int i, int j, tbint_type oper_type,
                             int xstart, int xfence, int ystart, int yfence,
                             const double* ints);
    const double* retrieve_pair_block(int i, int j, tbint_type oper_type,
                                      double* buf = 0) const;
    void retrieve_pair_subblock(int i, int j, tbint_type oper_type,
                                int xstart, int xfence, int ystart, int yfence,
                                double* buf) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
