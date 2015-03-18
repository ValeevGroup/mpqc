//
// distarray4.h
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

#ifndef _chemistry_qc_mbptr12_distarray4_h
#define _chemistry_qc_mbptr12_distarray4_h

#include <vector>
#include <util/ref/ref.h>
#include <util/state/state.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/group/memory.h>
#include <util/group/message.h>

namespace sc {

  enum DistArray4Storage {DistArray4Storage_XY, DistArray4Storage_YX};

  struct DistArray4Dimensions {
    public:
      static const DistArray4Dimensions& default_dim() { return default_dim_; }
      DistArray4Dimensions(int num_te_types, int n1, int n2, int n3, int n4,
                           DistArray4Storage storage = DistArray4Storage_XY) :
        num_te_types_(num_te_types), n1_(n1), n2_(n2), n3_(n3), n4_(n4),
        storage_(storage)
        {
        }
      int num_te_types() const { return num_te_types_; }
      int n1() const { return n1_; }
      int n2() const { return n2_; }
      int n3() const { return n3_; }
      int n4() const { return n4_; }
      DistArray4Storage storage() const { return storage_; }
    private:
      static DistArray4Dimensions default_dim_;
      int num_te_types_;
      int n1_;
      int n2_;
      int n3_;
      int n4_;
      DistArray4Storage storage_;
  };
  bool operator==(const DistArray4Dimensions& A,
                  const DistArray4Dimensions& B);

/////////////////////////////////////////////////////////////////
/**DistArray4 contains a set of one or more distributed dense 4-index arrays.
   These 4-index quantities are typically AO or MO
   integrals. Since integrals are typically computed in sets -- for example, R12 methods
   require integrals of operators 1/r12, r12, [t1,r12], and [t2,r12] at the same time -- several sets of
   integrals are held together. The data can then be described as (ijxy)^O
   where i, j, x, and, y are indices with ranges [0,ni), [0,nj), etc. and O labels
   the operator type and ranges from 0 to num_te_types. The data is stored and accessed as follows:
   each ij block is a set of num_te_types base-0 contiguous 2-dimensional array with dimensions
   nx and ny. How blocks are stored and accessed is determined in the derived class. It is also possible
   to have the physical storage to be xy or yx. \sa DistArray4::storage()

   Public interface of DistArray4 is designed to accomodate the needs of the TwoBodyMOIntsTransform
   objects. Parallel AO->MO integral transforms are performed in single or multiple passes.
   In the latter case all (jxy)^O integrals are produces for a particular subrange of i.
   These integrals are contained in a MemoryGrp object. The contents of the MemoryGrp object is
   "stored" in accumulator using <tt>store_memorygrp(Ref<MemoryGrp>& mem, int ni)</tt>, where
   ni is the size of the subrange of i produced in this pass.
   After all batches have been stored, the content of DistArray4 needs to be "committed"
   using <tt>commit</tt>. After that blocks of MO integrals can be accessed using
   <tt>retrieve_pair_block</tt>.

    */

class DistArray4: virtual public SavableState {
  public:
    // mem will be used to fetch data
    DistArray4(int num_te_types, int ni, int nj, int nx, int ny,
               DistArray4Storage storage = DistArray4Storage_XY);
    DistArray4(StateIn&);
    virtual ~DistArray4();
    void save_data_state(StateOut&);

    /** how to clone. optional dim allows to obtain an object of the same type but different size.
        the default is to obtain an object of the same size. */
    virtual Ref<DistArray4> clone(const DistArray4Dimensions& dim = DistArray4Dimensions::default_dim()) =0;

    /// Types of two-body operators that DistArray4 understands
    typedef unsigned int tbint_type;
    static const unsigned int max_num_te_types = 14;

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
    /// physical storage of the integrals. The default storage is XY. Storage is not mutable.
    const DistArray4Storage& storage() const { return storage_; }
    /// Size of each block of the integrals of one type, in double words
    size_t blocksize() const { return blocksize_; };
    // same as blocksize(), but in bytes
    size_t blksize() const { return blksize_; }

    /// call this before operations on this object can begin
    virtual void activate() { active_ = true; }
    /// call this after operations on this object are finished. May destroy data (see data_persistent()).
    virtual void deactivate() { active_ = false; }
    /// if this returns false, call to deactivate may destroy data
    virtual bool data_persistent() const =0;
    /** Retrieves an ij block of integrals. Note that it comes stored according to storage().
        No locking is performed.

        \param buf specifies the buffer in which to write the data (if not provided, will allocate dynamically).
        this buffer will be used by subsequent retrieve_pair_block() requests until release_pair_block() is called.
      */
    virtual const double * retrieve_pair_block(int i, int j, tbint_type oper_type,
                                               double* buf = 0) const =0;
    /** Retrieves a rectangular subblock of ij block of integrals. \sa retrieve_pair_block()
      */
    virtual void retrieve_pair_subblock(int i, int j, tbint_type oper_type,
                                        int xstart, int xfence, int ystart, int yfence,
                                        double* buf) const =0;
    /// Releases the buffer that holds ij block of integrals. If it was allocated by DistArray4, it will be freed.
    virtual void release_pair_block(int i, int j, tbint_type oper_type) const =0;
    /// Stores an ij pair block of integrals. It is assumed to be stored according to storage().
    /// It does not have to be "local", i.e. is_local(i,j) can be false, but:
    ///  - nonlocal writes may require additional copies and/or communication.
    ///  - and currently fail for some implementations :-(
    /// Locking is performed, hence should be thread safe
    virtual void store_pair_block(int i, int j, tbint_type oper_type, const double* ints) =0;
    /// Stores an rectangular subblock of ij block of integrals. Most efficient if the block is contiguous, i.e. yfence - ystart = ny().
    /// \sa store_pair_block()
    virtual void store_pair_subblock(int i, int j, tbint_type oper_type,
                                     int xstart, int xfence, int ystart, int yfence,
                                     const double* ints) =0;

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
    int tasks_with_access(std::vector<int>& twa_map) const;

    const Ref<MessageGrp>& msg() const { return msg_; }

    // return active_
    bool active() const { return active_; }

  private:
    /// Set to nonzero to debug this and derived classes
    static const int classdebug_ = 0;

    int num_te_types_;  // Number of types of integrals in a block
    Ref<MessageGrp> msg_;
    int ni_, nj_;
    int nx_, ny_;
    DistArray4Storage storage_;
    size_t nxy_;        // nx_ * ny_  - the number of integrals of one type in a block
    size_t blksize_;    // the same in bytes
    size_t blocksize_;  // hence the size of the block of num_te_types of integrals is blksize_ * num_te_types

    bool active_;

   protected:
    // return nxy_
    size_t nxy() const { return nxy_; }
    /// total number of tasks
    int ntasks() const { return msg_->n(); }
    /// rank of this task
    int me() const { return msg_->me(); }
    /// return debug level for this class
    int classdebug() const { return classdebug_; }

};

Ref<DistArray4> make_distarray4(int num_te_types, int ni, int nj, int nx, int ny,
                                DistArray4Storage storage = DistArray4Storage_XY);

/// extracts te_type from A
Ref<DistArray4>
extract(const Ref<DistArray4>& A,
        unsigned int te_type,
        double scale = 1.0);

/** creates an array in which indices 2 and 3 are permuted
*/
Ref<DistArray4> permute23(const Ref<DistArray4>& src);

/** creates an array in which indices 3 and 4 are permuted
*/
Ref<DistArray4> permute34(const Ref<DistArray4>& src);

/** creates an array in which indices 1 and 2 are permuted*/
Ref<DistArray4> permute12(const Ref<DistArray4>& src);


/// axpy followed by scaling: Y += a*X; Y *= scale.
void
axpy(const Ref<DistArray4>& X,
     double a,
     const Ref<DistArray4>& Y,
     double scale = 1.0);

/** antisymmetrizes the 4-index array in-place:
    <ij|xy> =  ( (ij|xy) + (ji|yx) - (ij|yx) - (ji|xy) ) / 2

    \param A input tensor. on output contains the antisymmetrized tensor. Valid A will obey these conditions: A->ni() == A->nj(), A->nx() == A->ny()
  */
void antisymmetrize(const Ref<DistArray4>& A);

/** symmetrizes the 4-index array in-place:
    <ij|xy> =  ( (ij|xy) + (ji|yx) ) / 2

    \param A input tensor. on output contains the symmetrized tensor. Valid A will obey these conditions: A->ni() == A->nj(), A->nx() == A->ny()
  */
void symmetrize(const Ref<DistArray4>& A);

/// contracts ijxy ("bra") with klxy ("ket") to produce ijkl ("braket")
void contract34(Ref<DistArray4>& braket,
                double scale,
                const Ref<DistArray4>& bra,
                unsigned int intsetidx_bra,
                const Ref<DistArray4>& ket,
                unsigned int intsetidx_ket,
                int debug = 0);

class RefSCMatrix;

/// contracts ijxy("bra") with klxy ("ket", a RefSCMatrix) to produce ijkl ("braket"); The last two arguments
/// definie the dimension for the resultant DistArray4 object, since a RefSCMatrix itself only gives the the product of dimensions
void contract34_DA4_RefMat(Ref<DistArray4>& braket,
                                              double scale,
                                              const Ref<DistArray4>& bra,
                                              unsigned int intsetidx_bra,
                                              const RefSCMatrix& ket,
                                              const int MatBra1Dim, const int MatBra2Dim);


/// Contract  X^Aj_ik = DA^Ax_iy * M^jk_xy; This is written for spin free[2]R12, where rdm matrices can be
/// easily permuted to assume the needed form; while write a general function involving, say, M^jy_kx
/// would need to permute k&y, which then needs the dimension information. So we prefer to pre-arrange M
/// so that we can directly use contract34_DA4_RefMat.
/// procedure: DA'^Ai_xy <- DA^Ax_iy from permute23
///            then C^Ai_jk <- DA'^Ai_xy ** M^jk_xy from contract34_DA4_RefMat
///            then X^Aj_ik <- C^Ai_jk from permute23 again.


void contract_DA4_RefMat_k2b2_34(Ref<DistArray4>& braket,
                       double scale,
                       const Ref<DistArray4>& bra,
                       unsigned int intsetidx_bra,
                       const RefSCMatrix& ket,
                       const int MatBra1Dim, const int MatBra2Dim);


/// Contract  X^Aj_ik = DA^Ax_yi * M^jk_xy;
/// procedure: DA'^Ax_iy <- DA^Ax_yi from permute34
///            DA''^Ai_xy <- DA'^Ax_iy from permute23
///            then C^Ai_jk <- DA'^Ai_xy ** M^jk_xy from contract34_DA4_RefMat
///            then X^Aj_ik <- C^Ai_jk from permute23 again.

void contract_DA4_RefMat_k1b2_34(Ref<DistArray4>& braket,
                       double scale,
                       const Ref<DistArray4>& bra,
                       unsigned int intsetidx_bra,
                       const RefSCMatrix& ket,
                       const int MatBra1Dim, const int MatBra2Dim);


/// contracts ijxy with T_xz to produce ijzy
void contract3(const Ref<DistArray4>& ijxy, const RefSCMatrix& T, Ref<DistArray4>& ijzy);

/// contracts ijxy with T_yz to produce ijxz
void contract4(const Ref<DistArray4>& ijxy, const RefSCMatrix& T, Ref<DistArray4>& ijxz);



/// copies contents of src into dst
RefSCMatrix &
operator<<(RefSCMatrix& dst,
                        const Ref<DistArray4>& src);

/// copy a specific tensor to RefSCMatrix
RefSCMatrix &
copy_to_RefSCMat(RefSCMatrix& dst,
                        const Ref<DistArray4>& src, const int tensor_type);

RefSCMatrix
copy_to_RefSCMat(const Ref<DistArray4>& src, int tensor_index);







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
  void store_memorygrp(Ref<DistArray4>& acc, Ref<MemoryGrp>& mem, int i_offset,
                       int ni, const size_t blksize_memgrp = 0);

  /** Reverse of store_memorygrp(). storage specifies the target storage of integrals in mem.
   * Give acc->storage() to maintain the same storage as in acc.
   */
  void restore_memorygrp(Ref<DistArray4>& acc, Ref<MemoryGrp>& mem, int i_offset,
                         int ni, DistArray4Storage storage, const size_t blksize_memgrp = 0);

}

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
