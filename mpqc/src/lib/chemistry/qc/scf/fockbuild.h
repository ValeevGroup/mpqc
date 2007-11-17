//
// fockbuild.h --- a generic Fock matrix builder
//
// Based on code gbuild.h:
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: SNL
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

#ifndef _chemistry_qc_scf_fockbuild_h
#define _chemistry_qc_scf_fockbuild_h

#ifdef __GNUC__
#pragma interface
#endif

#include <scconfig.h>
#include <util/misc/regtime.h>
#include <util/group/thread.h>
#include <util/group/message.h>
#include <chemistry/qc/basis/integral.h>

#include <util/group/actmsg.h>
#include <chemistry/qc/scf/fockdist.h>

namespace sc {

class FockBuildMatrix: public RefCount {
  private:
    void operator = (const FockBuildMatrix&) {}
    FockBuildMatrix() {}
  protected:
    Ref<MessageGrp> msg_;
    int nI_, nJ_;
    int nproc_;
    int me_;

    Ref<FockBlocks> blocks1_;
    Ref<FockBlocks> blocks2_;
  public:
    FockBuildMatrix(const FockBuildMatrix&);
    FockBuildMatrix(const Ref<MessageGrp> &msg);

    const Ref<MessageGrp> &messagegrp() const { return msg_; }

    inline int block_offset(int I, int J) const {
      if (!symmetric()) {
          return I * blocks2_->nblock() + J;
        }
      else {
          if (J > I) {
              std::cout << "shell_block_offset: noncanonical indices: "
                        << I << ", " << J << std::endl;
              abort();
            }
          return (I*(I+1))/2 + J;
        }
    }
    inline int n_block() const {
      if (!symmetric()) {
          return blocks1_->nblock() * blocks2_->nblock();
        }
      else {
          return (blocks1_->nblock()*(blocks1_->nblock()+1))/2;
        }
    }
    inline int block_owning_proc(int I, int J) const {
      return block_offset(I,J)%nproc_;
    }
    inline bool block_is_owner(int I, int J) const {
      return block_owning_proc(I,J) == me_;
    }

    inline int shell_block_offset(int I, int J) const {
      if (!symmetric()) {
          return I * nJ_ + J;
        }
      else {
          if (J > I) {
              std::cout << "shell_block_offset: noncanonical indices: "
                        << I << ", " << J << std::endl;
              abort();
            }
          return (I*(I+1))/2 + J;
        }
    }
    inline int n_shell_block() const {
      if (!symmetric()) {
          return nI_ * nJ_;
        }
      else {
          return (nI_*(nI_+1))/2;
        }
    }
    inline int shell_block_owning_proc(int I, int J) const {
      return block_owning_proc(blocks1_->shell_to_block(I),
                               blocks2_->shell_to_block(J));
    }
    inline bool shell_block_is_owner(int I, int J) const {
      return shell_block_owning_proc(I,J) == me_;
    }

    virtual bool symmetric() const = 0;
    // This will average the off diagonal elements in each
    // diagonal block, iff the matrix is symmetric.
    virtual void fix_diagonal_blocks() const = 0;
    virtual void clear() = 0;
    // These allocate data to hold the SCMatrix.  The data is only
    // copied if copy is true, otherwise it is zeroed.
    virtual void scmat_to_data(const Ref<SymmSCMatrix> &m,
                       const Ref<GaussianBasisSet> &b, bool copy) = 0;
    virtual void scmat_to_data(const Ref<SCMatrix> &m,
                       const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2, bool copy) = 0;
    // This copies the held data back into the SCMatrix.
    virtual void data_to_scmat() const = 0;
    virtual void prefetch_block(int I, int J, int ifetch, int nfetch) = 0;
    virtual void finish_prefetch_block() = 0;
    virtual double *shell_block(int Ish, int Jsh) const = 0;
    virtual double *block(int Ish, int Jsh) const = 0;
    virtual FockBuildMatrix *copy(int unique_id,
                                  bool copy_data = false) const = 0;
    /// Zero out the data.
    virtual void zero_data() = 0;
    /// Accumulate fbm into this.
    virtual void accum(const Ref<FockBuildMatrix> &fbm) = 0;
    /// Accumulate remote contributions.
    virtual void accum_remote(const Ref<MessageGrp> &) = 0;
    /// The default print member does nothing.
    virtual void print() const;
    /// Flush the buffer cache (if it exists).
    virtual void flush();

    virtual void set_fockblocks(const Ref<FockBlocks> &b1,
                                const Ref<FockBlocks> &b2);

    const Ref<FockBlocks> &blocks1() { return blocks1_; }
    const Ref<FockBlocks> &blocks2() { return blocks2_; }

};

class ReplFockBuildMatrix: public FockBuildMatrix {
    double **blockpointers_; /** Points to the blocks within the data
                              * array. */

    Ref<SCMatrix> rectmat_;
    Ref<SymmSCMatrix> symmmat_;
    int ndata_;
    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    bool owns_data_;

    void data_to_symmat() const;
    void data_to_rectmat() const;
  public:
    ReplFockBuildMatrix(const Ref<MessageGrp> &msg);
    ReplFockBuildMatrix(const ReplFockBuildMatrix &,
                        bool copy_data = false);
    ~ReplFockBuildMatrix();
    bool symmetric() const;
    void fix_diagonal_blocks() const;
    void clear();
    void scmat_to_data(const Ref<SymmSCMatrix> &m,
                       const Ref<GaussianBasisSet> &b, bool copy);
    void scmat_to_data(const Ref<SCMatrix> &m,
                       const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2, bool copy);
    void data_to_scmat() const;
    void prefetch_block(int I, int J, int ifetch, int nfetch);
    void finish_prefetch_block();
    double *shell_block(int Ish, int Jsh) const;
    double *block(int Ish, int Jsh) const;
    FockBuildMatrix *copy(int unique_id, bool copy_data = false) const;
    void zero_data();
    void accum(const Ref<FockBuildMatrix> &fbm);
    void accum_remote(const Ref<MessageGrp> &);
    virtual void print() const;
};

class FockContribution;

class FockBuildAMG: public ActiveMessageGrp {
    Ref<FockContribution> fc_;
    Ref<MessageGrp> return_msg_;
  public:
    FockBuildAMG(const Ref<MessageGrp>& msg,
                 const Ref<ThreadGrp>& thr,
                 const Ref<FockContribution>& fc);
    const Ref<FockContribution> &contrib() const {
      return fc_;
    }
    const Ref<MessageGrp> &return_messagegrp() const {
      return return_msg_;
    }
};

class FockBuildAM: public ActiveMessage {
  protected:
    int matrix_, I_, J_, nIJ_;
    double *data_;
    int use_shell_blocks_;
  public:
    FockBuildAM(int matrix);
    ~FockBuildAM();
    virtual FockBuildAM* copy() = 0;
    FockBuildAM(StateIn &s);
    void save_data_state(StateOut &s);
    void set_info(int I, int J, int nIJ, double *data,
                  bool use_shell_blocks) {
      I_ = I;
      J_ = J;
      nIJ_ = nIJ;
      data_ = data;
      use_shell_blocks_ = use_shell_blocks;
    }
};

class FockBuildOp: public RefCount {
  private:
    int msg_type_;
    bool need_read_;
    bool need_write_;
    bool use_shell_blocks_;
    Ref<FockBuildAMG> fbamg_;
    Ref<FockBuildAM> fbam_;
    Ref<StateSend> statesend_;
  public:
    FockBuildOp(int unique_id,
                bool need_read,
                bool need_write,
                bool use_shell_blocks,
                const Ref<FockBuildAMG> &fbamg,
                const Ref<FockBuildAM>& fbam);
    void process_write(int node, int I, int J, int nIJ, double *data);
    void process_read(int node, int I, int J, int nIJ, double *data);
    void begin_prefetch(int node, int I, int J, int nIJ, int ifetch, int nfetch,
                     double *data, MessageGrp::MessageHandle &mh);
    void finish_prefetch(MessageGrp::MessageHandle &mh);
    bool need_read() const { return need_read_; }
    bool need_write() const { return need_write_; }
    FockBuildOp *copy(int unique_id) const;
    Ref<MessageGrp> messagegrp() const {
      return fbamg_->messagegrp();
    }
    Ref<MessageGrp> return_messagegrp() const {
      return fbamg_->return_messagegrp();
    }
};

class DistFockBuildMatrix: public FockBuildMatrix {
    double **blockpointers_; /** Points to the blocks within the data
                              * array. */

    Ref<SCMatrix> rectmat_;
    Ref<SymmSCMatrix> symmmat_;
    int ndata_;
    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    bool owns_data_;

    bool prefetch_blocks_;

    Ref<FockBuildOp> fockbuildop_;

    class IJ_t {
        uint64_t ij_;
      public:
        IJ_t(const IJ_t&IJ) { ij_ = IJ.ij_; }
        IJ_t(int i, int j) { ij_ = uint64_t(i)<<32|j; }
        int i() const { return ij_>>32; }
        int j() const { return ij_&0xffffffff; }
        bool operator < (const IJ_t&IJ) const { return ij_ < IJ.ij_; }
        bool operator > (const IJ_t&IJ) const { return ij_ > IJ.ij_; }
        bool operator >= (const IJ_t&IJ) const { return ij_ >= IJ.ij_; }
        bool operator <= (const IJ_t&IJ) const { return ij_ <= IJ.ij_; }
        bool operator == (const IJ_t&IJ) const { return ij_ == IJ.ij_; }
        void operator = (const IJ_t&IJ) { ij_ = IJ.ij_; }
    };
    typedef std::pair<double*,MessageGrp::MessageHandle> FetchData_t;

    mutable std::map<IJ_t,FetchData_t> prefetched_block_cache_;
    mutable std::map<IJ_t,double*> block_cache_;
    mutable std::map<IJ_t,double*> shell_block_cache_;

    std::map<IJ_t,double*> local_blocks_;
    std::map<IJ_t,double*> local_shell_blocks_;

    void clear_cache();

    void data_to_symmat() const;
    void data_to_rectmat() const;

    double *fetch_block(int I, int J) const;
    int block_size(int iblock, int jblock) const;
    void insert_shell_block_pointers(int iblock,int jblock,
                                     double *data,
                                     std::map<IJ_t,double*> &) const;
    void blockpointers_to_local_blocks();
    void local_blocks_to_blockpointers() const;
  public:
    DistFockBuildMatrix(bool prefetch_blocks,
                        const Ref<FockBuildOp> &fockbuildop);
    DistFockBuildMatrix(const DistFockBuildMatrix &,
                        int unique_id,
                        bool copy_data = false);
    ~DistFockBuildMatrix();
    bool symmetric() const;
    void fix_diagonal_blocks() const;
    void clear();
    void scmat_to_data(const Ref<SymmSCMatrix> &m,
                       const Ref<GaussianBasisSet> &b, bool copy);
    void scmat_to_data(const Ref<SCMatrix> &m,
                       const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2, bool copy);
    void data_to_scmat() const;
    void prefetch_block(int I, int J, int ifetch, int nfetch);
    void finish_prefetch_block();
    double *shell_block(int Ish, int Jsh) const;
    double *block(int Ish, int Jsh) const;
    FockBuildMatrix *copy(int unique_id, bool copy_data = false) const;
    void zero_data();
    void accum(const Ref<FockBuildMatrix> &fbm);
    void accum_remote(const Ref<MessageGrp> &);
    void flush();
};

class FockContribution: public RefCount {
  protected:
    double nint_;
  public:
    FockContribution();
    FockContribution(const FockContribution&);
    virtual ~FockContribution();
    /** This routine does not permute any indices.  The matrix elements are
        contracted with the integrals are F_IJ and P_KL.  If F is
        symmetric, then J>I will be ignored.
        K>=L.
    */
    virtual void contrib_e_J(double factor,
                             int I, int J, int K, int L,
                             int nI, int nJ, int nK, int nL,
                             const double * restrictxx buf) = 0;
    /** This routine does not permute any indices.  The matrix elements
        are contracted with the integrals are F_IK and P_JL.  If F is
        symmetric, then K>I will be ignored.
    */
    virtual void contrib_e_K(double factor,
                             int I, int J, int K, int L,
                             int nI, int nJ, int nK, int nL,
                             const double * restrictxx buf) = 0;
    virtual void contrib_p12_p13p24_J(double factor,
                                      int I, int J, int K, int L,
                                      int nI, int nJ, int nK, int nL,
                                      const double * restrictxx buf) = 0;
    virtual void contrib_p12_p13p24_K(double factor,
                                      int I, int J, int K, int L,
                                      int nI, int nJ, int nK, int nL,
                                      const double * restrictxx buf) = 0;
    virtual void contrib_p34_p13p24_J(double factor,
                                      int I, int J, int K, int L,
                                      int nI, int nJ, int nK, int nL,
                                      const double * restrictxx buf) = 0;
    virtual void contrib_p34_p13p24_K(double factor,
                                      int I, int J, int K, int L,
                                      int nI, int nJ, int nK, int nL,
                                      const double * restrictxx buf) = 0;
    virtual void contrib_p12_p34_J(double factor,
                                   int I, int J, int K, int L,
                                   int nI, int nJ, int nK, int nL,
                                   const double * restrictxx buf) = 0;
    virtual void contrib_p12_p34_K(double factor,
                                   int I, int J, int K, int L,
                                   int nI, int nJ, int nK, int nL,
                                   const double * restrictxx buf) = 0;
    virtual void contrib_p34_J(double factor,
                               int I, int J, int K, int L,
                               int nI, int nJ, int nK, int nL,
                               const double * restrictxx buf) = 0;
    virtual void contrib_p34_K(double factor,
                               int I, int J, int K, int L,
                               int nI, int nJ, int nK, int nL,
                               const double * restrictxx buf) = 0;
    virtual void contrib_p13p24_J(double factor,
                                  int I, int J, int K, int L,
                                  int nI, int nJ, int nK, int nL,
                                  const double * restrictxx buf) = 0;
    virtual void contrib_p13p24_K(double factor,
                                  int I, int J, int K, int L,
                                  int nI, int nJ, int nK, int nL,
                                  const double * restrictxx buf) = 0;
    virtual void contrib_all_J(double factor,
                               int I, int J, int K, int L,
                               int nI, int nJ, int nK, int nL,
                               const double * restrictxx buf) = 0;
    virtual void contrib_all_K(double factor,
                               int I, int J, int K, int L,
                               int nI, int nJ, int nK, int nL,
                               const double * restrictxx buf) = 0;
    virtual Ref<FockContribution> clone() = 0;

    virtual void set_fmat(int i, const Ref<SCMatrix> &) = 0;
    virtual void set_fmat(int i, const Ref<SymmSCMatrix> &) = 0;

    virtual void set_jmat(int i, const Ref<SCMatrix> &) = 0;
    virtual void set_jmat(int i, const Ref<SymmSCMatrix> &) = 0;

    virtual void set_kmat(int i, const Ref<SCMatrix> &) = 0;
    virtual void set_kmat(int i, const Ref<SymmSCMatrix> &) = 0;

    virtual void set_pmat(int i, const Ref<SymmSCMatrix> &) = 0;

    virtual double *jmat_shell_block(int i, int Ish, int Jsh) = 0;
    virtual double *kmat_shell_block(int i, int Ish, int Jsh) = 0;
    virtual const double *pmat_shell_block(int i, int Ish, int Jsh) = 0;

    virtual double *jmat_block(int i, int Ish, int Jsh) = 0;
    virtual double *kmat_block(int i, int Ish, int Jsh) = 0;
    virtual const double *pmat_block(int i, int Ish, int Jsh) = 0;

    /** Compute the maximum of the density in each block.  The pmax vector
        holds only the unique elements. */
    virtual signed char *compute_pmax() const = 0;

    /// Copy matrices to allow multiple threads to coexist.
    virtual void copy_matrices(int unique_id) = 0;
    /** Sum the Fock matrix contributions from different threads.  The
     *  passed specialization type must be the same as the specialization
     *  of this.
     */
    virtual void accum(const Ref<FockContribution> &) = 0;
    /** Sum the Fock matrix contributions from different processors.
        This might be a no-op for distributed matrices.
     */
    virtual void accum_remote(const Ref<MessageGrp> &) = 0;
    /// Push the internal Fock matrix data back into the original object.
    virtual void update() = 0;

    double nint() const { return nint_; }
    double &nint() { return nint_; }

    virtual void activate() = 0;
    virtual void sync() = 0;
    virtual void deactivate() = 0;

    virtual void flush() = 0;

    virtual void prefetch_blocks(int I, int J, int ifetch, int nfetch) = 0;
    virtual void finish_prefetch_blocks() = 0;

    virtual void set_fockblocks(const Ref<FockBlocks> &blocks_f1,
                                const Ref<FockBlocks> &blocks_f2,
                                const Ref<FockBlocks> &blocks_p) = 0;

    virtual Ref<ThreadLock> &get_lock(int i, int I, int J) = 0;
};

/** The GenericFockContribution class provides much of the infrastructure
    needed by FockContribution specializations.

    In addition to the members that GenericFockContribution provides,
    FockContribution classes must have members that actually do the
    work and sum in the contributions.  Due to two electron integral
    permutation symmetry, one integral may make several contributions.
    A contrib member exists for each subgroup of the full permutation
    group.

    <dl>

    <dt>void f()<dd>The f member.

    </dl>

    A lightweight copy CTOR should also exist for Contribution
    specializations.
*/
class GenericFockContribution: public FockContribution {
  protected:
    int nfmat_;     /// the number of Fock matrices
    std::vector<Ref<FockBuildMatrix> > jmats_;
    std::vector<Ref<FockBuildMatrix> > kmats_;
    std::vector<bool> k_is_j_;
    int npmat_;     /// the number of density matrices
    std::vector<Ref<FockBuildMatrix> > pmats_;
    Ref<GaussianBasisSet> f_b1_, f_b2_, p_b_;
    bool f_b1_equiv_f_b2;
    int nlocks_;
    std::vector<Ref<ThreadLock> > locks_;
    std::string fockbuildmatrixtype_;
    bool use_shell_blocks_;

    Ref<FockBuildAMG> fbamg_;

    FockBuildMatrix *fockbuildmatrix(int matrix,
                                     const std::string &type,
                                     const Ref<MessageGrp> &msg,
                                     const Ref<FockBuildAMG> &);

    class JLocator {
      public:
        double *operator()(GenericFockContribution *owner,
                           int i, int I, int J) {
          return owner->jmat_shell_block(i,I,J);
        }
    };

    class KLocator {
      public:
        double *operator()(GenericFockContribution *owner,
                           int i, int I, int J) {
          return owner->kmat_shell_block(i,I,J);
        }
    };

    template <class Locator>
    class JKBlock {
        GenericFockContribution *owner_;
        double *data_;
        int i_, I_, J_, nIJ_;
      public:
        JKBlock(GenericFockContribution *owner,
                int i, int I, int J, int nI, int nJ) {
          i_ = i;
          I_ = I;
          J_ = J;
          nIJ_ = nI*nJ;
          owner_ = owner;
          data_ = owner_->alloc_scratch(nIJ_);
        }
        ~JKBlock() {
          Locator l;
          int ilock, jlock;
          if (owner_->use_shell_blocks()) {
              ilock = I_;
              jlock = J_;
            }
          else {
              ilock = owner_->jmat(i_)->blocks1()->shell_to_block(I_);
              jlock = owner_->jmat(i_)->blocks2()->shell_to_block(J_);
            }
          const Ref<ThreadLock> &lock(
              owner_->get_lock(i_,ilock,jlock));
          lock->lock();
          double *real_data = l(owner_,i_,I_,J_);
//           std::cout << MessageGrp::get_default_messagegrp()->me()
//                     << ": writing back F block " << I_ << " " << J_
//                     << ": contrib[0] = " << data_[0]
//                     << " @ " << data_
//                     << " oldvalue[0] = " << real_data[0]
//                     << " @ " << real_data
//                     << " summedvalue[0] = " << real_data[0] + data_[0]
//                     << std::endl;
          for (int i=0; i<nIJ_; i++) {
//               std::cout << "  data_[" << i << "] = "
//                         << scprintf("%12.9f", data_[i])
//                         << std::endl;
              real_data[i] += data_[i];
            }
//           for (int i=0; i<nIJ_; i++) {
//               std::cout << "  real_data[" << i << "] = "
//                         << scprintf("%12.9f", real_data[i])
//                         << std::endl;
//             }
          lock->unlock();
          owner_->free_scratch(data_);
        }
        double *data() { return data_; }
    };

    class PBlock {
        GenericFockContribution *owner_;
        const double *data_;
      public:
        PBlock(GenericFockContribution *owner,
               int i, int I, int J, int nI, int nJ) {
          owner_ = owner;
          data_ = owner_->pmat_shell_block(i,I,J);
        }
        ~PBlock() {}
        const double *data() { return data_; }
    };

    GenericFockContribution(int nfmat, int npmat,
                            const Ref<GaussianBasisSet> &f_b1,
                            const Ref<GaussianBasisSet> &f_b2,
                            const Ref<GaussianBasisSet> &p_b,
                            const std::string &fockbuildmatrixtype);

    void pmax_contrib(const Ref<FockBuildMatrix> &mat,
                      signed char *pmax) const;

  public:
    double *jmat_shell_block(int i, int I, int J) {
      return jmats_[i]->shell_block(I,J);
    }
    bool jmat_symmetric(int i) const { return jmats_[i]->symmetric(); }
    double *kmat_shell_block(int i, int I, int J) {
      return kmats_[i]->shell_block(I,J);
    }
    bool kmat_symmetric(int i) const { return kmats_[i]->symmetric(); }
    const double *pmat_shell_block(int i, int I, int J) {
      return pmats_[i]->shell_block(I,J);
    }

    double *jmat_block(int i, int I, int J) {
      return jmats_[i]->block(I,J);
    }
    double *kmat_block(int i, int I, int J) {
      return kmats_[i]->block(I,J);
    }
    const double *pmat_block(int i, int I, int J) {
      return pmats_[i]->block(I,J);
    }

    Ref<ThreadLock> &get_lock(int i, int Ish, int Jsh) {
      int hash = (i+(Ish+1)*(Jsh+1))%nlocks_;
      return locks_[hash];
    }

    double *alloc_scratch(int size) {
      double *data = new double[size];
      memset(data,0,sizeof(double)*size);
      return data;
    }

    void free_scratch(double *data) {
      delete[] data;
    }

    void set_fmat(int i, const Ref<SCMatrix> &);
    void set_fmat(int i, const Ref<SymmSCMatrix> &);

    void set_jmat(int i, const Ref<SCMatrix> &);
    void set_jmat(int i, const Ref<SymmSCMatrix> &);

    void set_kmat(int i, const Ref<SCMatrix> &);
    void set_kmat(int i, const Ref<SymmSCMatrix> &);

    void set_pmat(int i, const Ref<SymmSCMatrix> &);

    void copy_matrices(int unique_id);
    void accum(const Ref<FockContribution> &);
    void accum_remote(const Ref<MessageGrp> &);
    void update();

    signed char* compute_pmax() const;

    ~GenericFockContribution();

    void activate();
    void sync();
    void deactivate();

    void prefetch_blocks(int I, int J, int ifetch, int nfetch);
    void finish_prefetch_blocks();

    void set_fockblocks(const Ref<FockBlocks> &blocks_f1,
                        const Ref<FockBlocks> &blocks_f2,
                        const Ref<FockBlocks> &blocks_p);

    void flush();

    const Ref<FockBuildMatrix> &jmat(int i) { return jmats_[i]; }
    const Ref<FockBuildMatrix> &kmat(int i) { return kmats_[i]; }
    const Ref<FockBuildMatrix> &pmat(int i) { return pmats_[i]; }

    bool use_shell_blocks() const { return use_shell_blocks_; }
};

/** The FockBuildThread class is used to actually build the Fock matrix.
    It is used by the FockBuilder class.
 */
class FockBuildThread : public Thread {
  protected:
    Ref<FockDistribution> fockdist_;
    Ref<FockContribution> contrib_;
    Ref<ThreadLock> lock_;
    Ref<Integral> integral_;
    double accuracy_;
    Ref<MessageGrp> msg_;
    int nthread_;
    int threadnum_;
    const signed char *pmax_;
    Ref<RegionTimer> timer_;
    bool prefetch_blocks_;

    int can_sym_offset(int i, int j) { return (i*(i+1))/2 + j; }
    int gen_sym_offset(int i, int j) {
      if (i>=j) { return can_sym_offset(i,j); }
      else      { return can_sym_offset(j,i); }
    }
  public:
    /// Each thread must be given a unique contribution, c.
    FockBuildThread(const Ref<FockDistribution> &fockdist,
                    const Ref<MessageGrp> &msg,
                    int nthread,
                    int threadnum,
                    bool prefetch_blocks,
                    const Ref<ThreadLock> &lock,
                    const Ref<Integral> &integral);
    void set_contrib(const Ref<FockContribution>&c) { contrib_ = c; }
    void set_accuracy(double acc) { accuracy_ = acc; }
    void set_pmax(const signed char *pmax) { pmax_ = pmax; }
    const Ref<RegionTimer> get_timer() const { return timer_; }
};

/** The FockBuildThread class is used to actually build the Fock matrix.
    It is used by the FockBuilder class.
 */
class FockBuildThread_F11_P11 : public FockBuildThread {
    Ref<GaussianBasisSet> basis_;
    Ref<FockBlocks> blocks_;
    Ref<PetiteList> pl_;
    Ref<TwoBodyInt> eri_;
    void prefetch_blocks(const Ref<FockDist> &dist,
                         int iblock, int jblock, int kblock, int lblock);
  public:
    /// Each thread must be given a unique contribution, c.
    FockBuildThread_F11_P11(const Ref<FockDistribution> &fockdist,
                            const Ref<MessageGrp> &msg,
                            int nthread,
                            int threadnum,
                            bool prefetch_blocks,
                            const Ref<ThreadLock> &lock,
                            const Ref<Integral> &integral,
                            const Ref<PetiteList> &pl,
                            const Ref<GaussianBasisSet> &basis1,
                            const Ref<GaussianBasisSet> &basis2/*not used*/,
                            const Ref<GaussianBasisSet> &basis3/*not used*/,
                            const Ref<FockBlocks> &blocks1,
                            const Ref<FockBlocks> &blocks2/*not used*/,
                            const Ref<FockBlocks> &blocks3/*not used*/);
    void run();
    void new_run();
    void old_run();
};

/** This is used to build the Fock matrix when none of the
    basis sets are equivalent.
 */
class FockBuildThread_F12_P33 : public FockBuildThread {
    Ref<GaussianBasisSet> basis1_;
    Ref<GaussianBasisSet> basis2_;
    Ref<GaussianBasisSet> basis3_;
    Ref<PetiteList> pl_;

    void run_J();
    void run_K();

    Ref<TwoBodyInt> eri_J_;
    Ref<TwoBodyInt> eri_K_;

  public:
    /// Each thread must be given a unique contribution, c.
    FockBuildThread_F12_P33(const Ref<FockDistribution> &fockdist,
                            const Ref<MessageGrp> &msg,
                            int nthread,
                            int threadnum,
                            bool prefetch_blocks,
                            const Ref<ThreadLock> &lock,
                            const Ref<Integral> &integral,
                            const Ref<PetiteList> &pl,
                            const Ref<GaussianBasisSet> &basis1,
                            const Ref<GaussianBasisSet> &basis2,
                            const Ref<GaussianBasisSet> &basis3,
                            const Ref<FockBlocks> &blocks1,
                            const Ref<FockBlocks> &blocks2,
                            const Ref<FockBlocks> &blocks3);
    void run();
};

/** The FockBuild class works with the FockBuildThread class to generate
    Fock matrices for both closed shell and open shell methods.  It uses a
    helper class, FockContribution, to do the work of forming
    contributions from the density matrices and integrals and placing these
    into the partial Fock matrices (the G matrices).
*/
class FockBuild: public RefCount {
    Ref<FockDistribution> fockdist_;
    Ref<FockContribution> contrib_;
    Ref<GaussianBasisSet> b_f1_;
    Ref<GaussianBasisSet> b_f2_;
    Ref<GaussianBasisSet> b_p_;
    Ref<MessageGrp> msg_;
    Ref<ThreadGrp> thr_;
    Ref<Integral> integral_;
    double accuracy_;
    Ref<PetiteList> pl_;
    bool prefetch_blocks_;

    Ref<FockBlocks> fb_f1_;
    Ref<FockBlocks> fb_f2_;
    Ref<FockBlocks> fb_p_;

    typedef FockBuildThread* (*FBT_CTOR)(const Ref<FockDistribution> &fockdist,
                                         const Ref<MessageGrp> &msg,
                                         int nthread,
                                         int threadnum,
                                         bool prefetch_blocks,
                                         const Ref<ThreadLock> &lock,
                                         const Ref<Integral> &integral,
                                         const Ref<PetiteList> &pl,
                                         const Ref<GaussianBasisSet> &basis1,
                                         const Ref<GaussianBasisSet> &basis2,
                                         const Ref<GaussianBasisSet> &basis3,
                                         const Ref<FockBlocks> &blocks1,
                                         const Ref<FockBlocks> &blocks2,
                                         const Ref<FockBlocks> &blocks3);


    // Build for the any case.  The thread constructing function is passed in.
    FockBuildThread **thread_;
    void init_threads(FBT_CTOR);
    void init_threads();
    void done_threads();

  public:
    /** Create a FockBuild object using b_f1 as the Fock matrix row
        dimension basis, b_f2 as the Fock matrix column dimension basis,
        and b_p as the density matrix dimensions.  If b_f2 is not given,
        then b_f1 is used.  If b_p1 is not given, then b_f1 is used.  If
        b_p2 is not given, then b_p1 is used.  If the following parameters
        are not given, then the global defaults are used: The msg parameter
        specifies the MessageGrp, thr gives the ThreadGrp, and integral
        gives the Integral.  */
    FockBuild(const Ref<FockDistribution> &fockdist,
              const Ref<FockContribution> &contrib,
              bool prefetch_blocks,
              const Ref<GaussianBasisSet> &b_f1,
              const Ref<GaussianBasisSet> &b_f2 = 0,
              const Ref<GaussianBasisSet> &b_p = 0,
              const Ref<MessageGrp> &msg=MessageGrp::get_default_messagegrp(),
              const Ref<ThreadGrp> &thr=ThreadGrp::get_default_threadgrp(),
              const Ref<Integral> &integral=Integral::get_default_integral());
    virtual ~FockBuild();

    /** Contruct the Fock matrices. */
    void build();

    const Ref<FockContribution> &contrib() const { return contrib_; }
    void set_accuracy(double acc) { accuracy_ = acc; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
