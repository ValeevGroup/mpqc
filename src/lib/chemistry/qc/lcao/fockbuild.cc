//
// fockbuild.cc --- a generic Fock matrix builder
//
// Based on code:
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

#include <iomanip>
#include <memory>
#include <stdexcept>
#include <vector>

#include <math.h>
#include <float.h>

#include <util/misc/regtime.h>
#include <util/misc/exenv.h>
#include <util/misc/scint.h>
#include <util/misc/autovec.h>
#include <util/misc/scexception.h>
#include <math/scmat/elemop.h>
#include <math/scmat/blkiter.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/lcao/fockbuild.h>
#include <chemistry/qc/lcao/fockdist.h>

#undef DEBUG
#define DEBUG 0

#undef DEBUG_DIST
#define DEBUG_DIST 0

#undef DETAILED_TIMINGS
#define DETAILED_TIMINGS 0

#undef SCF_USE_BOUNDS
#define SCF_USE_BOUNDS 1

namespace sc {

/////////////////////////////////////////////////////////////////
// FockBuildMatrixRectElemOp

class FockBuildMatrixRectElemOp: public SCElementOp {
    bool symm_;        /// true if matrix is symmetric
    bool data_to_mat_; /// true if we should accumulate data to matrix
    double **blocks_;  /// pointers to the data
    int ndata_;        /// the total number of data
    bool defer_collect_;
    Ref<SCBlockInfo> rowbi_;
    Ref<SCBlockInfo> colbi_;
  public:
    // If data_to_mat is true the data in blocks is summed into the
    // matrix, otherwise the data in the matrix is copied into the
    // blocks.
    FockBuildMatrixRectElemOp(double **blocks,
                              bool symm, bool data_to_mat,
                              int ndata,
                              const Ref<SCBlockInfo> &rowbi,
                              const Ref<SCBlockInfo> &colbi):
      symm_(symm),
      data_to_mat_(data_to_mat),
      blocks_(blocks),
      ndata_(ndata),
      defer_collect_(false),
      rowbi_(rowbi),
      colbi_(colbi) {
      double *data = blocks_[0];
      if (!data_to_mat_) for (int i=0; i<ndata; i++) data[i] = 0;
#if DEBUG
      std::cout << "blockinfo:" << std::endl;
      rowbi->print();
#endif
    }
    int has_collect() { return 1; }
    void defer_collect(int defer) { defer_collect_ = defer; }
    void collect(const Ref<MessageGrp> &msg) {
      if (!defer_collect_) msg->sum(blocks_[0], ndata_);
    }
    void collect(const Ref<SCElementOp>& op) {
      SCElementOp::collect(op);
    }
    int has_side_effects() { return data_to_mat_; }
    bool threadsafe() const { return true; }
    bool cloneable() const { return false; }
    /** This assumes that the blocks in the matrix actually
     *  correspond to the blocks in the SCDimension's BlockInfo.
     *  This assertion is not tested, however. */
    void process(SCMatrixBlockIter&iter) {
      for (iter.reset(); iter; ++iter) {
          int I, J, blocki, blockj;
          rowbi_->elem_to_block(iter.i(), I, blocki);
          colbi_->elem_to_block(iter.j(), J, blockj);
          int IJ;
          if (symm_) {
              IJ = (I*(I+1))/2 + J;
            }
          else {
              IJ = colbi_->nblock() * I + J;
            }
          int ni = rowbi_->size(I);
          int nj = colbi_->size(J);
          int blockoffset = blocki * nj + blockj;
          double *data = blocks_[IJ];
          double value;
          if (data_to_mat_) {
              iter.set(value = (data[blockoffset]+iter.get()));
            }
          else {
              value = data[blockoffset] = iter.get();
              // diagonal blocks hold the full square array,
              // rather than just the lower triangle
              if (symm_ && I == J) {
                  data[blockj*ni+blocki] = value;
                }
            }
#if DEBUG
          std::cout << "Data_" << I << J
                    << "("  << blocki << blockj
                    << "@" << std::setw(2) << blockoffset << ")";
          if (data_to_mat_) std::cout << " -> ";
          else              std::cout << " <- ";
          std::cout << "SCMat_" << iter.i() << iter.j();
          std::cout << "; value = "
                    << std::setw(6) << value << std::endl;
#endif
        }
    }
};

/////////////////////////////////////////////////////////////////
// FockBuildMatrix

FockBuildMatrix::FockBuildMatrix(const FockBuildMatrix&fbm)
{
  nI_ = fbm.nI_;
  nJ_ = fbm.nJ_;
  me_ = fbm.me_;
  nproc_ = fbm.nproc_;
  msg_ = fbm.msg_;
  blocks1_ = fbm.blocks1_;
  blocks2_ = fbm.blocks2_;
}

FockBuildMatrix::FockBuildMatrix(const Ref<MessageGrp> &msg):
  msg_(msg),
  nI_(0),
  nJ_(0)
{
  me_ = msg->me();
  nproc_ = msg->n();
}

void
FockBuildMatrix::set_fockblocks(const Ref<FockBlocks> &b1,
                                const Ref<FockBlocks> &b2)
{
  blocks1_ = b1;
  blocks2_ = b2;
}

void
FockBuildMatrix::print() const
{
}

void
FockBuildMatrix::flush()
{
}

/////////////////////////////////////////////////////////////////
// ReplFockBuildMatrix

ReplFockBuildMatrix::ReplFockBuildMatrix(const Ref<MessageGrp> &msg):
  FockBuildMatrix(msg),
  blockpointers_(0),
  ndata_(0),
  bs1_(0),
  bs2_(0),
  owns_data_(false)
{
}

ReplFockBuildMatrix::ReplFockBuildMatrix(const ReplFockBuildMatrix &fbm,
                                         bool copy_data):
  FockBuildMatrix(fbm)
{
  blockpointers_ = fbm.blockpointers_;
  rectmat_ = fbm.rectmat_;
  symmmat_ = fbm.symmmat_;
  ndata_ = fbm.ndata_;
  bs1_ = fbm.bs1_;
  bs2_ = fbm.bs2_;
  owns_data_ = false;

  if (copy_data) {
      int nIJ = n_shell_block();
      auto_vec<double*> tmp_blockpointers(new double*[nIJ]);
      auto_vec<double> tmp_data(new double[ndata_]);
      double *base1 = blockpointers_[0];
      double *base2 = tmp_data.get();
      for (int ij=0; ij<nIJ; ij++) {
          tmp_blockpointers[ij] = base2 + (blockpointers_[ij] - base1);
        }
      blockpointers_ = tmp_blockpointers.release();
      tmp_data.release(); // this data is now referenced by blockpointers_
      memcpy(blockpointers_[0], fbm.blockpointers_[0], sizeof(double)*nIJ);
      owns_data_ = true;
    }
}

ReplFockBuildMatrix::~ReplFockBuildMatrix()
{
  clear();
}

void
ReplFockBuildMatrix::clear()
{
  if (owns_data_ && blockpointers_ != 0) {
      delete[] blockpointers_[0];
      delete[] blockpointers_;
    }
  nI_ = nJ_ = ndata_ = 0;
  blockpointers_ = 0;
  symmmat_ = 0;
  rectmat_ = 0;
  bs1_ = 0;
  bs2_ = 0;
}

FockBuildMatrix *
ReplFockBuildMatrix::copy(int unique_id, bool copy_data) const
{
  ReplFockBuildMatrix *cp = new ReplFockBuildMatrix(*this, copy_data);
  return cp;
}

void
ReplFockBuildMatrix::zero_data()
{
  double *tmp = blockpointers_[0];
  for (int i=0; i<ndata_; i++) tmp[i] = 0.0;
}

void
ReplFockBuildMatrix::fix_diagonal_blocks() const
{
  if (!symmetric()) return;
  Ref<SCBlockInfo> bi = bs1_->basisdim()->blocks();
  for (int I=0; I<nI_; I++) {
      double *dat = shell_block(I,I);
      int nI = bi->size(I);
      for (int i=0; i<nI; i++) {
          for (int j=0; j<i; j++) {
              int ij = i*nI+j;
              int ji = j*nI+i;
              double val = 0.5*(dat[ij] + dat[ji]);
              dat[ij] = val;
              dat[ji] = val;
            }
        }
    }
}

void
ReplFockBuildMatrix::scmat_to_data(const RefSymmSCMatrix &m,
                                   const Ref<GaussianBasisSet> &b,
                                   bool copy)
{
  clear();

  symmmat_ = m;
  bs1_ = b;
  bs2_ = b;

  // bdim contains the shell blocking information.
  // the matrix dimension is an petite list AO dimension
  // which does not contain the needed shell information
  Ref<SCDimension> bdim = b->basisdim();

  if (bdim->n() != m->dim()->n()) {
      std::cout << "b->basisdim()" << std::endl;
      b->basisdim()->print();
      std::cout << "m->dim()" << std::endl;
      m->dim()->print();
      throw std::invalid_argument("scmat_to_data: bad dimension (sym)");
    }

  Ref<SCBlockInfo> biI(bdim->blocks());
  Ref<SCBlockInfo> biJ(bdim->blocks());
  nI_ = biI->nblock();
  nJ_ = biJ->nblock();

  ndata_ = 0;
  for (int I=0; I<nI_; I++) {
      for (int J=0; J<=I; J++) {
          ndata_ +=  biI->size(I) * biJ->size(J);
        }
    }

  auto_vec<double*> tmp_blockpointers(new double*[n_shell_block()]);
  auto_vec<double> tmp_data(new double[ndata_]);

  blockpointers_ = tmp_blockpointers.release();
  blockpointers_[0] = tmp_data.release();

  double *current_data = blockpointers_[0];
  for (int I=0; I<nI_; I++) {
      for (int J=0; J<=I; J++) {
          blockpointers_[shell_block_offset(I,J)] = current_data;
          current_data +=  biI->size(I) * biJ->size(J);
        }
    }

  if (copy) {
      bool symm = true;
      bool data_to_mat = false;
      Ref<SCElementOp> op
          = new FockBuildMatrixRectElemOp(blockpointers_,
                                          symm, data_to_mat,
                                          ndata_,
                                          bs1_->basisdim()->blocks(),
                                          bs2_->basisdim()->blocks());
      symmmat_->element_op(op);
    }
  else {
      double *tmp = blockpointers_[0];
      for (int i=0; i<ndata_; i++) tmp[i] = 0.0;
    }
}

void
ReplFockBuildMatrix::scmat_to_data(const RefSCMatrix &m,
                                   const Ref<GaussianBasisSet> &b1,
                                   const Ref<GaussianBasisSet> &b2,
                                   bool copy)
{
  if (!b1->basisdim()->equiv(m->rowdim())
      || !b2->basisdim()->equiv(m->coldim())) {
      throw std::invalid_argument("scmat_to_data: bad dimension (rect)");
    }

  rectmat_ = m;
  bs1_ = b1;
  bs2_ = b2;

  // bdim contains the shell blocking information.
  // the matrix dimension is an petite list AO dimension
  // which does not contain the needed shell information
  Ref<SCDimension> bdim1 = b1->basisdim();
  Ref<SCDimension> bdim2 = b2->basisdim();

  Ref<SCBlockInfo> biI(bdim1->blocks());
  Ref<SCBlockInfo> biJ(bdim2->blocks());
  nI_ = biI->nblock();
  nJ_ = biJ->nblock();

  ndata_ = 0;
  for (int I=0; I<nI_; I++) {
      for (int J=0; J<nJ_; J++) {
          ndata_ +=  biI->size(I) * biJ->size(J);
        }
    }

  auto_vec<double*> tmp_blockpointers(new double*[n_shell_block()]);
  auto_vec<double> tmp_data(new double[ndata_]);

  blockpointers_ = tmp_blockpointers.release();
  blockpointers_[0] = tmp_data.release();

  double *current_data = blockpointers_[0];
  for (int I=0; I<nI_; I++) {
      for (int J=0; J<nJ_; J++) {
          blockpointers_[shell_block_offset(I,J)] = current_data;
          current_data +=  biI->size(I) * biJ->size(J);
        }
    }

  if (copy) {
      bool symm = false;
      bool data_to_mat = false;
      Ref<SCElementOp> op
          = new FockBuildMatrixRectElemOp(blockpointers_,
                                          symm, data_to_mat,
                                          ndata_,
                                          bs1_->basisdim()->blocks(),
                                          bs2_->basisdim()->blocks());
      rectmat_->element_op(op);
    }
  else {
      double *tmp = blockpointers_[0];
      for (int i=0; i<ndata_; i++) tmp[i] = 0.0;
    }
}

void
ReplFockBuildMatrix::data_to_rectmat() const
{
  bool symm = false;
  bool data_to_mat = true;
  Ref<SCElementOp> op
      = new FockBuildMatrixRectElemOp(blockpointers_,
                                      symm, data_to_mat,
                                      ndata_,
                                      bs1_->basisdim()->blocks(),
                                      bs2_->basisdim()->blocks());
  rectmat_->element_op(op);
}

void
ReplFockBuildMatrix::data_to_symmat() const
{
  bool symm = true;
  bool data_to_mat = true;
  Ref<SCElementOp> op
      = new FockBuildMatrixRectElemOp(blockpointers_,
                                      symm, data_to_mat,
                                      ndata_,
                                      bs1_->basisdim()->blocks(),
                                      bs2_->basisdim()->blocks());
  symmmat_->element_op(op);
}

void
ReplFockBuildMatrix::data_to_scmat() const
{
  if (blockpointers_ == 0) return;
  if (symmmat_) {
      if (rectmat_) {
          throw std::runtime_error("data_to_scmat: both mats nonnull");
        }
      data_to_symmat();
    }
  else if (rectmat_) {
      data_to_rectmat();
    }
  else {
      throw std::runtime_error("data_to_scmat: both mats null");
    }
}

void
ReplFockBuildMatrix::accum(const Ref<FockBuildMatrix> &fa)
{
  Ref<ReplFockBuildMatrix> f;
  f << fa;
  if (ndata_ != f->ndata_) {
      throw std::invalid_argument("incompatible ReplFockBuildMatrix in accum");
    }
  double *src = f->blockpointers_[0];
  double *dst = blockpointers_[0];
  for (int i=0; i<ndata_; i++) dst[i] += src[i];
}

void
ReplFockBuildMatrix::accum_remote(const Ref<MessageGrp> &msg)
{
  if (blockpointers_)
    msg->sum(blockpointers_[0], ndata_);
//   std::cout << "after accum_remote blockpointers_[0][ndata_-2] = "
//             << blockpointers_[0][ndata_-2]
//             << std::endl;
}

void
ReplFockBuildMatrix::prefetch_block(int I, int J, int ifetch, int nfetch)
{
  // all data is already local so this is a no-op
}

void
ReplFockBuildMatrix::finish_prefetch_block()
{
  // all data is already local so this is a no-op
}

double *
ReplFockBuildMatrix::shell_block(int I, int J) const
{
  return blockpointers_[shell_block_offset(I,J)];
}

double *
ReplFockBuildMatrix::block(int I, int J) const
{
  // this should not be called.
  throw ProgrammingError("ReplFockBuildMatrix::block: "
                             "should never be called",
                             __FILE__, __LINE__);
  return 0;
}

bool
ReplFockBuildMatrix::symmetric() const
{
  return symmmat_;
}

void
ReplFockBuildMatrix::print() const
{
  Ref<SCBlockInfo> bi = bs1_->basisdim()->blocks();
  Ref<SCBlockInfo> bj = bs2_->basisdim()->blocks();
  for (int I=0; I<nI_; I++) {
      int nI = bi->size(I);
      for (int J=0; J<nJ_; J++) {
          if (symmetric() && J > I) continue;
          int nJ = bj->size(J);
          double *dat = shell_block(I,J);
          ExEnv::out0() << indent
                            << "Block " << I << ", " << J
                            << std::endl;
          for (int i=0; i<nI; i++) {
              for (int j=0; j<nJ; j++) {
                  if (symmetric() && I == J && j > i) continue;
                  ExEnv::out0() << indent
                                    << "  " << i << ", " << j << " : "
                                    << dat[i*nJ+j]
                                    << std::endl;
                }
            }
        }
    }

}

/////////////////////////////////////////////////////////////////
// FockBuildAMG

FockBuildAMG::FockBuildAMG(const Ref<MessageGrp>& msg,
                           const Ref<ThreadGrp>& thr,
                           const Ref<FockContribution>& fc):
  ActiveMessageGrp(msg,thr),
  fc_(fc)
{
  return_msg_ = msg->clone();
}

/////////////////////////////////////////////////////////////////
// FockBuildAM

static ClassDesc FockBuildAM_cd(
  typeid(FockBuildAM),"FockBuildAM",1,"public ActiveMessage",
  0, 0, 0);

FockBuildAM::FockBuildAM(int matrix):
  matrix_(matrix),
  I_(0), J_(0), nIJ_(0), data_(0) {
}

FockBuildAM::~FockBuildAM() {
  delete[] data_;
}

FockBuildAM::FockBuildAM(StateIn &s):
  SavableState(s),
  ActiveMessage(s)
{
  s.get(matrix_);
  s.get(I_);
  s.get(J_);
  s.get(nIJ_);
  s.get(use_shell_blocks_);
  int have_data;
  s.get(have_data);
  if (have_data) {
      data_ = new double[nIJ_];
      s.get_array_double(data_,nIJ_);
    }
  else {
      data_ = 0;
    }
}

void
FockBuildAM::save_data_state(StateOut &s)
{
  ActiveMessage::save_data_state(s);
  s.put(matrix_);
  s.put(I_);
  s.put(J_);
  s.put(nIJ_);
  s.put(use_shell_blocks_);
  int have_data = (data_ != 0);
  s.put(have_data);
  if (have_data) {
      s.put_array_double(data_,nIJ_);
    }
}

/////////////////////////////////////////////////////////////////
// FockBuildAMGetP

class FockBuildAMGetP: public FockBuildAM {
  public:
    FockBuildAMGetP(StateIn &s):
      SavableState(s),
      FockBuildAM(s) {
    };
    FockBuildAMGetP(int matrix): FockBuildAM(matrix) {}
    void run(int sender, int type, ActiveMessageGrp *context) {
      Ref<FockBuildAMG> fbamg;
      fbamg << context;
      const double *data;
      if (use_shell_blocks_) {
          data = fbamg->contrib()->pmat_shell_block(matrix_, I_, J_);
        }
      else {
          data = fbamg->contrib()->pmat_block(matrix_, I_, J_);
        }
#if DEBUG_DIST
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ": pmat = " << matrix_
                << ", I = " << I_
                << ", J = " << J_
                << ", nIJ = " << nIJ_
                << ", sender = " << sender
                << ", type = " << type
                << ", data[0] = " << data[0]
                << std::endl;
#endif
      fbamg->return_messagegrp()->sendt(sender,type,data,nIJ_,true);
    }
    virtual FockBuildAM* copy() { return new FockBuildAMGetP(matrix_); }
};

static ClassDesc FockBuildAMGetP_cd(
  typeid(FockBuildAMGetP),"FockBuildAMGetP",1,"public FockBuildAM",
  0, 0, create<FockBuildAMGetP>);

/////////////////////////////////////////////////////////////////
// FockBuildAMAccumJ

class FockBuildAMAccumJ: public FockBuildAM {
  public:
    FockBuildAMAccumJ(StateIn &s):
      SavableState(s),
      FockBuildAM(s) {
    }
    FockBuildAMAccumJ(int matrix): FockBuildAM(matrix) {}
    void run(int sender, int type, ActiveMessageGrp *context) {
      Ref<FockBuildAMG> fbamg;
      fbamg << context;
      double *data;
      if (use_shell_blocks_) {
          data = fbamg->contrib()->jmat_shell_block(matrix_, I_, J_);
        }
      else {
          data = fbamg->contrib()->jmat_block(matrix_, I_, J_);
        }
#if DEBUG_DIST
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ": jmat = " << matrix_
                << ", I = " << I_
                << ", J = " << J_
                << ", nIJ = " << nIJ_
                << ", sender = " << sender
                << ", type = " << type
                << ", contrib[0] = " << data_[0]
                << " oldvalue[0] = " << data[0]
                << " summedvalue[0] = " << data[0] + data_[0]
                << std::endl;
#endif
      ThreadLockHolder lockholder(
          fbamg->contrib()->get_lock(matrix_, I_, J_));
      for (int i=0; i<nIJ_; i++) {
          data[i] += data_[i];
        }
    }
    virtual FockBuildAM* copy() { return new FockBuildAMAccumJ(matrix_); }
};

static ClassDesc FockBuildAMAccumJ_cd(
  typeid(FockBuildAMAccumJ),"FockBuildAMAccumJ",1,"public FockBuildAM",
  0, 0, create<FockBuildAMAccumJ>);

/////////////////////////////////////////////////////////////////
// FockBuildAMAccumK

class FockBuildAMAccumK: public FockBuildAM {
  public:
    FockBuildAMAccumK(StateIn &s):
      SavableState(s),
      FockBuildAM(s) {
    }
    FockBuildAMAccumK(int matrix): FockBuildAM(matrix) {}
    void run(int sender, int type, ActiveMessageGrp *context) {
      Ref<FockBuildAMG> fbamg;
      fbamg << context;
      double *data;
      if (use_shell_blocks_) {
          data = fbamg->contrib()->kmat_shell_block(matrix_, I_, J_);
        }
      else {
          data = fbamg->contrib()->kmat_block(matrix_, I_, J_);
        }
      ThreadLockHolder lockholder(
          fbamg->contrib()->get_lock(matrix_, I_, J_));
#if DEBUG_DIST
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ": kmat = " << matrix_
                << ", I = " << I_
                << ", J = " << J_
                << ", nIJ = " << nIJ_
                << ", sender = " << sender
                << ", type = " << type
                << ", data_[0] = " << data_[0]
                << std::endl;
#endif
      for (int i=0; i<nIJ_; i++) {
          data[i] += data_[i];
        }
    }
    virtual FockBuildAM* copy() { return new FockBuildAMAccumK(matrix_); }
};

static ClassDesc FockBuildAMAccumK_cd(
  typeid(FockBuildAMAccumK),"FockBuildAMAccumK",1,"public FockBuildAM",
  0, 0, create<FockBuildAMAccumK>);

/////////////////////////////////////////////////////////////////
// FockBuildOp

FockBuildOp::FockBuildOp(int unique_id,
                         bool need_read,
                         bool need_write,
                         bool use_shell_blocks,
                         const Ref<FockBuildAMG> &fbamg,
                         const Ref<FockBuildAM>& fbam):
  msg_type_(unique_id),
  need_read_(need_read),
  need_write_(need_write),
  use_shell_blocks_(use_shell_blocks),
  fbamg_(fbamg),
  fbam_(fbam)
{
  statesend_ = fbamg->get_statesend();
  statesend_->type(msg_type_);
}

FockBuildOp *
FockBuildOp::copy(int unique_id) const {
  return new FockBuildOp(unique_id,need_read_,need_write_,use_shell_blocks_,
                         fbamg_, fbam_->copy());
}

void
FockBuildOp::process_write(int node, int I, int J, int nIJ, double *data)
{
#if DEBUG_DIST
  std::cout << MessageGrp::get_default_messagegrp()->me()
            << ": process_write: I = " << I
            << ", J = " << J
            << ", nIJ = " << nIJ
            << ", data[0] = " << data[0]
            << std::endl;
#endif
  fbam_->set_info(I,J,nIJ,data,use_shell_blocks_);
  fbamg_->send(node, statesend_, fbam_.pointer());
  fbam_->set_info(0,0,0,0,false);
}

void
FockBuildOp::process_read(int node, int I, int J, int nIJ, double *data)
{
#if DEBUG_DIST
  std::cout << MessageGrp::get_default_messagegrp()->me()
            << ": process read: I = " << I
            << ", J = " << J
            << ", nIJ = " << nIJ
            << std::endl;
#endif
  fbam_->set_info(I,J,nIJ,0,use_shell_blocks_);
  MessageGrp::MessageHandle handle;
  fbamg_->return_messagegrp()->nb_recvt(node, msg_type_, data, nIJ, handle);
  fbamg_->send(node, statesend_, fbam_.pointer());
  fbamg_->return_messagegrp()->wait(handle);
  fbam_->set_info(0,0,0,0,false);
}

void
FockBuildOp::begin_prefetch(int node, int I, int J, int nIJ, int ifetch, int nfetch,
                         double *data, MessageGrp::MessageHandle &mh)
{
#if DEBUG_DIST
  std::cout << MessageGrp::get_default_messagegrp()->me()
            << ": process fetch: I = " << I
            << ", J = " << J
            << ", nIJ = " << nIJ
            << std::endl;
#endif
  fbam_->set_info(I,J,nIJ,0,use_shell_blocks_);
  int unique_msg_type = nfetch*msg_type_ + ifetch;
  fbamg_->return_messagegrp()->nb_recvt(node, unique_msg_type, data, nIJ, mh);
  statesend_->type(unique_msg_type);
  fbamg_->send(node, statesend_, fbam_.pointer());
  fbam_->set_info(0,0,0,0,false);
}

void
FockBuildOp::finish_prefetch(MessageGrp::MessageHandle &mh)
{
#if DEBUG_DIST
  std::cout << MessageGrp::get_default_messagegrp()->me()
            << ": finish_prefetch"
            << std::endl;
#endif
  fbamg_->return_messagegrp()->wait(mh);
}

/////////////////////////////////////////////////////////////////
// DistFockBuildMatrix

DistFockBuildMatrix::DistFockBuildMatrix(
    bool prefetch_blocks,
    const Ref<FockBuildOp> &fockbuildop):
  FockBuildMatrix(fockbuildop->messagegrp()),
  blockpointers_(0),
  ndata_(0),
  bs1_(0),
  bs2_(0),
  owns_data_(false),
  prefetch_blocks_(prefetch_blocks),
  fockbuildop_(fockbuildop)
{
}

DistFockBuildMatrix::DistFockBuildMatrix(const DistFockBuildMatrix &fbm,
                                         int unique_id,
                                         bool copy_data):
  FockBuildMatrix(fbm)
{
  blockpointers_ = fbm.blockpointers_;
  rectmat_ = fbm.rectmat_;
  symmmat_ = fbm.symmmat_;
  ndata_ = fbm.ndata_;
  bs1_ = fbm.bs1_;
  bs2_ = fbm.bs2_;
  prefetch_blocks_ = fbm.prefetch_blocks_;
  owns_data_ = false;

  me_ = fbm.me_;
  nproc_ = fbm.nproc_;
  fockbuildop_ = fbm.fockbuildop_->copy(unique_id);

  local_blocks_ = fbm.local_blocks_;
  local_shell_blocks_ = fbm.local_shell_blocks_;

  if (copy_data) {
      for (std::map<IJ_t,double*>::iterator iter = local_blocks_.begin();
           iter != local_blocks_.end();
           iter++) {
          int iblock = iter->first.i();
          int jblock = iter->first.j();
          int size = block_size(iblock,jblock);
          double *src = iter->second;
          double *dst = new double[size];
          for (int i=0; i<size; i++) {
              dst[i] = src[i];
            }
          iter->second = dst;
          insert_shell_block_pointers(iblock,jblock,dst,local_shell_blocks_);
        }

      owns_data_ = true;
    }

  if (copy_data) {
      int nIJ = n_shell_block();
      auto_vec<double*> tmp_blockpointers(new double*[nIJ]);
      auto_vec<double> tmp_data(new double[ndata_]);
      double *base1 = blockpointers_[0];
      double *base2 = tmp_data.get();
      for (int ij=0; ij<nIJ; ij++) {
          tmp_blockpointers[ij] = base2 + (blockpointers_[ij] - base1);
        }
      blockpointers_ = tmp_blockpointers.release();
      tmp_data.release(); // this data is now referenced by blockpointers_
      memcpy(blockpointers_[0], fbm.blockpointers_[0], sizeof(double)*nIJ);
      owns_data_ = true;
    }
}

DistFockBuildMatrix::~DistFockBuildMatrix()
{
  clear();
}

void
DistFockBuildMatrix::clear()
{
  if (owns_data_ && blockpointers_ != 0) {
      delete[] blockpointers_[0];
      delete[] blockpointers_;
      for (std::map<IJ_t,double*>::iterator iter = local_blocks_.begin();
           iter != local_blocks_.end();
           iter++) {
          delete[] iter->second;
        }
    }
  local_blocks_.clear();
  local_shell_blocks_.clear();
  nI_ = nJ_ = ndata_ = 0;
  blockpointers_ = 0;
  symmmat_ = 0;
  rectmat_ = 0;
  bs1_ = 0;
  bs2_ = 0;
  clear_cache();
}

int
DistFockBuildMatrix::block_size(int iblock, int jblock) const
{
  Ref<SCBlockInfo> bi1 = bs1_->basisdim()->blocks();
  Ref<SCBlockInfo> bi2 = bs2_->basisdim()->blocks();

  int size = 0;
  for (int i=blocks1_->begin(iblock); i<blocks1_->end(iblock); i++) {
      int i_shell_size = bi1->size(i);
      int jend;
      if (symmetric() && iblock == jblock) jend = i+1;
      else jend = blocks2_->end(jblock);
      for (int j=blocks2_->begin(jblock); j<jend; j++) {
          int j_shell_size = bi2->size(j);
          size += i_shell_size * j_shell_size;
        }
    }

  return size;
}

void
DistFockBuildMatrix::insert_shell_block_pointers(
    int iblock, int jblock,
    double *block,
    std::map<IJ_t,double*> &shell_blocks) const
{
  Ref<SCBlockInfo> bi1 = bs1_->basisdim()->blocks();
  Ref<SCBlockInfo> bi2 = bs2_->basisdim()->blocks();

  double *offset = block;
  for (int i=blocks1_->begin(iblock); i<blocks1_->end(iblock); i++) {
      int i_shell_size = bi1->size(i);
      int jend;
      if (symmetric() && iblock == jblock) jend = i+1;
      else jend = blocks2_->end(jblock);
      for (int j=blocks2_->begin(jblock); j<jend; j++) {
          int j_shell_size = bi2->size(j);
          IJ_t key(i,j);
          shell_blocks[key] = offset;
#if DEBUG_DIST
          std::cout << scprintf("%2d: inserting shell block pointer for "
                                    "block %2d %2d shell %2d %2d @ ",
                                    me_, iblock, jblock, i, j)
                    << offset
                    << std::endl;
#endif
          offset += i_shell_size * j_shell_size;
        }
    }
}

FockBuildMatrix *
DistFockBuildMatrix::copy(int unique_id, bool copy_data) const
{
  DistFockBuildMatrix *cp
      = new DistFockBuildMatrix(*this, unique_id, copy_data);
  return cp;
}

void
DistFockBuildMatrix::zero_data()
{
  double *tmp = blockpointers_[0];
  for (int i=0; i<ndata_; i++) tmp[i] = 0.0;

  for (std::map<IJ_t,double*>::iterator iter = local_blocks_.begin();
       iter != local_blocks_.end();
       iter++) {
      int iblock = iter->first.i();
      int jblock = iter->first.j();
      int size = block_size(iblock,jblock);
      double *data = iter->second;
      for (int i=0; i<size; i++) {
          data[i] = 0.0;
        }
    }
}

void
DistFockBuildMatrix::fix_diagonal_blocks() const
{
  if (!symmetric()) return;
  Ref<SCBlockInfo> bi = bs1_->basisdim()->blocks();
  for (int I=0; I<nI_; I++) {
      // this is called after accum_remote, so must be run
      // on all blocks, on all nodes.
      //if (!shell_block_is_owner(I,I)) continue;
      //double *dat = shell_block(I,I);
      double *dat = blockpointers_[shell_block_offset(I,I)];
      int nI = bi->size(I);
      for (int i=0; i<nI; i++) {
          for (int j=0; j<i; j++) {
              int ij = i*nI+j;
              int ji = j*nI+i;
              double val = 0.5*(dat[ij] + dat[ji]);
#if DEBUG_DIST
              std::cout
                  << scprintf("%d: fixing %d%d %d%d: 1/2(%11.8f + %11.8f) -> %11.8f owner=%d", me_, I, I, i, j, dat[ij], dat[ji], val, shell_block_owning_proc(I,I))
                  << std::endl;
#endif
              dat[ij] = val;
              dat[ji] = val;
            }
        }
    }

  for (std::map<IJ_t,double*>::const_iterator
           iter = local_shell_blocks_.begin();
       iter != local_shell_blocks_.end();
       iter++) {
      int ishell = iter->first.i();
      int jshell = iter->first.j();
      if (ishell != jshell) continue;
      int nishell = bi->size(ishell);
      double *dat = iter->second;
      for (int i=0; i<nishell; i++) {
          for (int j=0; j<i; j++) {
              int ij = i*nishell+j;
              int ji = j*nishell+i;
              double val = 0.5*(dat[ij] + dat[ji]);
              dat[ij] = val;
              dat[ji] = val;
            }
        }
    }
}

void
DistFockBuildMatrix::scmat_to_data(const RefSymmSCMatrix &m,
                                   const Ref<GaussianBasisSet> &b,
                                   bool copy)
{
  clear();

  symmmat_ = m;
  bs1_ = b;
  bs2_ = b;

  // bdim contains the shell blocking information.
  // the matrix dimension is an petite list AO dimension
  // which does not contain the needed shell information
  Ref<SCDimension> bdim = b->basisdim();

  if (bdim->n() != m->dim()->n()) {
      std::cout << "b->basisdim()" << std::endl;
      b->basisdim()->print();
      std::cout << "m->dim()" << std::endl;
      m->dim()->print();
      throw std::invalid_argument("scmat_to_data: bad dimension (sym)");
    }

  Ref<SCBlockInfo> biI(bdim->blocks());
  Ref<SCBlockInfo> biJ(bdim->blocks());
  nI_ = biI->nblock();
  nJ_ = biJ->nblock();

  ndata_ = 0;
  for (int I=0; I<nI_; I++) {
      for (int J=0; J<=I; J++) {
          ndata_ +=  biI->size(I) * biJ->size(J);
        }
    }

  auto_vec<double*> tmp_blockpointers(new double*[n_shell_block()]);
  auto_vec<double> tmp_data(new double[ndata_]);

  blockpointers_ = tmp_blockpointers.release();
  blockpointers_[0] = tmp_data.release();

  double *current_data = blockpointers_[0];
  for (int I=0; I<nI_; I++) {
      for (int J=0; J<=I; J++) {
          blockpointers_[shell_block_offset(I,J)] = current_data;
          current_data +=  biI->size(I) * biJ->size(J);
        }
    }

  if (copy) {
      bool symm = true;
      bool data_to_mat = false;
      Ref<SCElementOp> op
          = new FockBuildMatrixRectElemOp(blockpointers_,
                                          symm, data_to_mat,
                                          ndata_,
                                          bs1_->basisdim()->blocks(),
                                          bs2_->basisdim()->blocks());
      symmmat_->element_op(op);
    }
  else {
      double *tmp = blockpointers_[0];
      for (int i=0; i<ndata_; i++) tmp[i] = 0.0;
    }

  if (prefetch_blocks_) {
      blockpointers_to_local_blocks();
    }
}

void
DistFockBuildMatrix::scmat_to_data(const RefSCMatrix &m,
                                   const Ref<GaussianBasisSet> &b1,
                                   const Ref<GaussianBasisSet> &b2,
                                   bool copy)
{
  if (!b1->basisdim()->equiv(m->rowdim())
      || !b2->basisdim()->equiv(m->coldim())) {
      throw std::invalid_argument("scmat_to_data: bad dimension (rect)");
    }

  rectmat_ = m;
  bs1_ = b1;
  bs2_ = b2;

  throw std::runtime_error("scmat_to_data: not yet impl (rect)");
}

void
DistFockBuildMatrix::data_to_rectmat() const
{
  throw std::runtime_error("data_to_rectmat: not yet impl");
}

void
DistFockBuildMatrix::data_to_symmat() const
{
  if (prefetch_blocks_) {
      local_blocks_to_blockpointers();
      msg_->sum(blockpointers_[0], ndata_);
    }

  bool symm = true;
  bool data_to_mat = true;
  Ref<SCElementOp> op
      = new FockBuildMatrixRectElemOp(blockpointers_,
                                      symm, data_to_mat,
                                      ndata_,
                                      bs1_->basisdim()->blocks(),
                                      bs2_->basisdim()->blocks());
  symmmat_->element_op(op);
}

void
DistFockBuildMatrix::data_to_scmat() const
{
  if (symmmat_) {
      if (rectmat_) {
          throw std::runtime_error("data_to_scmat: both mats nonnull");
        }
      data_to_symmat();
    }
  else if (rectmat_) {
      data_to_rectmat();
    }
  else {
      throw std::runtime_error("data_to_scmat: both mats null");
    }
}

void
DistFockBuildMatrix::blockpointers_to_local_blocks()
{
  for (int iblock=0; iblock<blocks1_->nblock(); iblock++) {
      int jblockend;
      if (symmetric()) jblockend=iblock+1;
      else jblockend=blocks2_->nblock();
      for (int jblock=0; jblock<jblockend; jblock++) {
          IJ_t key(iblock,jblock);
          int ndata = block_size(iblock,jblock);
          double *data = new double[ndata];
          local_blocks_[key] = data;
        }
    }

  for (std::map<IJ_t,double*>::iterator iter = local_blocks_.begin();
       iter != local_blocks_.end();
       iter++) {
      int iblock = iter->first.i();
      int jblock = iter->first.j();
      double *data = iter->second;
      insert_shell_block_pointers(iblock,jblock,data,local_shell_blocks_);
    }

  Ref<SCDimension> b1dim = bs1_->basisdim();
  Ref<SCDimension> b2dim = bs2_->basisdim();
  Ref<SCBlockInfo> biI(b1dim->blocks());
  Ref<SCBlockInfo> biJ(b2dim->blocks());
  for (std::map<IJ_t,double*>::iterator iter = local_shell_blocks_.begin();
       iter != local_shell_blocks_.end();
       iter++) {
      int ishell = iter->first.i();
      int jshell = iter->first.j();
      double *dst = iter->second;
      const double *src = blockpointers_[shell_block_offset(ishell,jshell)];
      int shell_block_size = biI->size(ishell) * biJ->size(jshell);
      memcpy(dst,src,sizeof(*dst)*shell_block_size);
    }
}

void
DistFockBuildMatrix::local_blocks_to_blockpointers() const
{
#if DEBUG_DIST
  std::cout << scprintf("%2d: local_blocks_to_blockpointers:", me_)
            << std::endl;
#endif
  for (int i=0; i<ndata_; i++) blockpointers_[0][i] = 0.0;
  Ref<SCDimension> b1dim = bs1_->basisdim();
  Ref<SCDimension> b2dim = bs2_->basisdim();
  Ref<SCBlockInfo> biI(b1dim->blocks());
  Ref<SCBlockInfo> biJ(b2dim->blocks());
  for (std::map<IJ_t,double*>::const_iterator iter
           = local_shell_blocks_.begin();
       iter != local_shell_blocks_.end();
       iter++) {
      int ishell = iter->first.i();
      int jshell = iter->first.j();
      const double *src = iter->second;
      double *dst = blockpointers_[shell_block_offset(ishell,jshell)];
      int shell_block_size = biI->size(ishell) * biJ->size(jshell);
      memcpy(dst,src,sizeof(*dst)*shell_block_size);
#if DEBUG_DIST
      std::cout << scprintf("%2d: shell block %2d %2d: data[0] = %12.8f",
                                me_, ishell, jshell, src[0])
                << std::endl;
#endif
    }
}

void
DistFockBuildMatrix::accum(const Ref<FockBuildMatrix> &fa)
{
  Ref<DistFockBuildMatrix> f;
  f << fa;
  if (ndata_ != f->ndata_) {
      throw std::invalid_argument("incompatible DistFockBuildMatrix in accum");
    }
  double *src = f->blockpointers_[0];
  double *dst = blockpointers_[0];
  for (int i=0; i<ndata_; i++) dst[i] += src[i];

  // this and f have the same local blocks
  std::map<IJ_t,double *>::iterator iter1,iter2;
  for (iter1 = local_blocks_.begin(), iter2 = f->local_blocks_.begin();
       iter1 != local_blocks_.end();
       iter1++, iter2++) {
      int iblock = iter1->first.i();
      int jblock = iter1->first.j();
      int size = block_size(iblock,jblock);
      double *dst = iter1->second;
      double *src = iter2->second;
      for (int i=0; i<size; i++) {
          dst[i] += src[i];
        }
    }
}

void
DistFockBuildMatrix::accum_remote(const Ref<MessageGrp> &msg)
{
  if (prefetch_blocks_) {
      // This is a no-op, since data is distributed.
    }
  else {
    if (blockpointers_)
      msg_->sum(blockpointers_[0], ndata_);
    }
}

bool
DistFockBuildMatrix::symmetric() const
{
  return symmmat_;
}

double *
DistFockBuildMatrix::shell_block(int I, int J) const
{
  int proc = shell_block_owning_proc(I,J);

#if DEBUG_DIST
  std::cout << MessageGrp::get_default_messagegrp()->me()
            << ": shell_block(" << I << "," << J << ")"
            << ", owner = " << proc
            << std::endl;
#endif

  if (proc == me_) {
      if (prefetch_blocks_) {
          std::map<IJ_t,double*>::const_iterator iter
              = local_shell_blocks_.find(IJ_t(I,J));
          if (iter == local_shell_blocks_.end()) {
              ExEnv::outn()
                  << indent
                  << MessageGrp::get_default_messagegrp()->me()
                  << ": ERROR: shell block "
                  << I << ", " << J
                  << " should be have been local, but was not found"
                  << std::endl;
              abort();
            }
#if DEBUG_DIST
          std::cout << MessageGrp::get_default_messagegrp()->me()
                    << ": block was local: data[0] = "
                    << *iter->second
                    << std::endl;
#endif
          return iter->second;
        }
#if DEBUG_DIST
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ": block was local: data[0] = "
                << blockpointers_[shell_block_offset(I,J)][0]
                << std::endl;
#endif
      return blockpointers_[shell_block_offset(I,J)];
    }

  IJ_t key(I,J);
  std::map<IJ_t,double*>::const_iterator
      location = shell_block_cache_.find(key);
  if (location != shell_block_cache_.end()) {
#if DEBUG_DIST
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ": block was in-cache: data[0] = "
                << location->second[0]
                << std::endl;
#endif
      return location->second;
    }

  if (prefetch_blocks_ && fockbuildop_->need_read()) {
      ExEnv::outn()
          << indent
          << MessageGrp::get_default_messagegrp()->me()
          << ": ERROR: shell block "
          << I << ", " << J
          << " should be have been prefetched, but was not found"
          << std::endl;
      abort();
    }

  if (prefetch_blocks_) {
      // When prefetch_blocks_ is true, the large block sizes are
      // used to accumulate data, too, so that must be handled as
      // a special case.

      // Note that need_read() must be false here.

      int iblock = blocks1_->shell_to_block(I);
      int jblock = blocks2_->shell_to_block(J);
      int ndata = block_size(iblock,jblock);
      double *data = new double[ndata];

#if DEBUG_DIST
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ": block is remote but write-only.  Zeroing."
                << std::endl;
#endif

      for (int i=0; i<ndata; i++) data[i] = 0.0;

      block_cache_[IJ_t(iblock,jblock)] = data;

      insert_shell_block_pointers(iblock,jblock,data,shell_block_cache_);

      // shell_block_cache_[IJ_t(I,J)] was initialized by the above
      // call to insert_shell_block_pointers
      return shell_block_cache_[IJ_t(I,J)];
    }
  else {
      Ref<SCBlockInfo> bi1 = bs1_->basisdim()->blocks();
      Ref<SCBlockInfo> bi2 = bs2_->basisdim()->blocks();
      int nI = bi1->size(I);
      int nJ = bi2->size(J);
      int nIJ = nI*nJ;

      double *data = new double[nIJ];
      if (fockbuildop_->need_read()) {
          // get the data
#if DEBUG_DIST
          std::cout << MessageGrp::get_default_messagegrp()->me()
                    << ": doing remote read of block"
                    << std::endl;
#endif
          fockbuildop_->process_read(proc,I,J,nIJ,data);
#if DEBUG_DIST
          std::cout << MessageGrp::get_default_messagegrp()->me()
                    << ": did remote read of block: data[0] = "
                    << data[0]
                    << std::endl;
#endif
        }
      else {
          // otherwise, just zero the data
#if DEBUG_DIST
          std::cout << MessageGrp::get_default_messagegrp()->me()
                    << ": block is remote but write-only.  Zeroing."
                    << std::endl;
#endif
          for (int i=0; i<nIJ; i++) data[i] = 0.0;
        }

      shell_block_cache_[key] = data;

      return data;
    }
}

double *
DistFockBuildMatrix::block(int I, int J) const
{
  IJ_t key(I,J);
  std::map<IJ_t, double*>::const_iterator iter = local_blocks_.find(key);
  if (iter == local_blocks_.end()) {
      throw ProgrammingError("DistFockBuildMatrix::block: "
                                 "request for nonexistent block",
                                 __FILE__, __LINE__);
    }
  return iter->second;
}

void
DistFockBuildMatrix::prefetch_block(int iblock, int jblock,
                                    int ifetch, int nfetch)
{
  int proc = block_owning_proc(iblock,jblock);

  if (proc == me_) return;

  if (!symmetric() || iblock>=jblock) {

      IJ_t key(iblock,jblock);
      if (prefetched_block_cache_.find(key)
          != prefetched_block_cache_.end()) {
          // the block has already been prefetched
          return;
        }
      FetchData_t &data = prefetched_block_cache_[key];

#if DEBUG_DIST
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ": prefetch_block("
                << iblock << ", " << jblock
                << ")"
                << std::endl;
#endif

      int ndata = block_size(iblock,jblock);

      data.first = new double[ndata];
      fockbuildop_->begin_prefetch(proc,iblock,jblock,ndata,
                                   ifetch, nfetch,
                                   data.first,data.second);
    }
}

void
DistFockBuildMatrix::finish_prefetch_block()
{
  std::map<IJ_t,FetchData_t>::iterator iter;
  for (iter=prefetched_block_cache_.begin(); iter != prefetched_block_cache_.end();
       iter++) {
      fockbuildop_->finish_prefetch(iter->second.second);
#if DEBUG_DIST
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ": " << iter->first.i()
                << ", " << iter->first.j()
                << ": prefetch finished"
                << std::endl;
#endif
      int iblock = iter->first.i();
      int jblock = iter->first.j();
      double *data = iter->second.first;
      block_cache_[iter->first] = data;
      insert_shell_block_pointers(iblock,jblock,data,shell_block_cache_);
    }
  prefetched_block_cache_.clear();
}

void
DistFockBuildMatrix::flush()
{
#if DEBUG_DIST
  std::cout << MessageGrp::get_default_messagegrp()->me()
            << ": DistFockBuildMatrix::flush(): called"
            << std::endl;
#endif
  if (fockbuildop_->need_write()) {
      if (prefetch_blocks_) {
          std::map<IJ_t,double*>::iterator iter;
          for (iter=block_cache_.begin();
               iter != block_cache_.end();
               iter++) {
              int iblock = iter->first.i();
              int jblock = iter->first.j();
              int ndata = block_size(iblock,jblock);
              double *data = iter->second;
              int proc = block_owning_proc(iblock,jblock);
#if DEBUG_DIST
              std::cout << MessageGrp::get_default_messagegrp()->me()
                        << ": flushing " << iblock << ", " << jblock
                        << " to " << proc
                        << " data[0] " << data[0]
                        << std::endl;
#endif
              fockbuildop_->process_write(proc,iblock,jblock,ndata,data);
            }
        }
      else {
          Ref<SCBlockInfo> bi1 = bs1_->basisdim()->blocks();
          Ref<SCBlockInfo> bi2 = bs2_->basisdim()->blocks();

          std::map<IJ_t,double*>::iterator iter;
          for (iter=shell_block_cache_.begin();
               iter != shell_block_cache_.end();
               iter++) {
              int I = iter->first.i();
              int J = iter->first.j();
              int nI = bi1->size(I);
              int nJ = bi2->size(J);
              int nIJ = nI*nJ;
              double *data = iter->second;
              int proc = shell_block_owning_proc(I,J);
#if DEBUG_DIST
              std::cout << MessageGrp::get_default_messagegrp()->me()
                        << ": flushing " << I << ", " << J
                        << " to " << proc
                        << " data[0] " << data[0]
                        << std::endl;
#endif
              fockbuildop_->process_write(proc,I,J,nIJ,data);
            }
        }
    }

  clear_cache();
}

void
DistFockBuildMatrix::clear_cache()
{
  if (prefetch_blocks_) {
      std::map<IJ_t,double*>::iterator iter;
      for (iter=block_cache_.begin(); iter != block_cache_.end();
           iter++) {
          double *data = iter->second;
          delete[] data;
        }
    }
  else {
      std::map<IJ_t,double*>::iterator iter;
      for (iter=shell_block_cache_.begin(); iter != shell_block_cache_.end();
           iter++) {
          double *data = iter->second;
          delete[] data;
        }
    }
#if DEBUG_DIST
  std::cout << me_ << ": clearing cache" << std::endl;
#endif
  block_cache_.clear();
  shell_block_cache_.clear();
}

/////////////////////////////////////////////////////////////////
// FockContribution


FockContribution::FockContribution():
  nint_(0)
{
}

FockContribution::~FockContribution()
{
}

/////////////////////////////////////////////////////////////////
// GenericFockContribution

GenericFockContribution::GenericFockContribution(
    int nfmat, int npmat,
    const Ref<GaussianBasisSet> &f_b1,
    const Ref<GaussianBasisSet> &f_b2,
    const Ref<GaussianBasisSet> &p_b,
    const std::string &fockbuildmatrixtype):
  nfmat_(nfmat),
  npmat_(npmat),
  jmats_(nfmat),
  kmats_(nfmat),
  k_is_j_(nfmat),
  pmats_(npmat),
  f_b1_(f_b1),
  f_b2_(f_b2),
  p_b_(p_b),
  fockbuildmatrixtype_(fockbuildmatrixtype)
{
  f_b1_equiv_f_b2 = f_b1_->equiv(f_b2);

  use_shell_blocks_
      = (fockbuildmatrixtype_ == "prefetched_distributed"?false:true);

  nlocks_ = 23;
  locks_.resize(nlocks_);

  Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();
  Ref<ThreadGrp>  thr = ThreadGrp::get_default_threadgrp();

  // FIXME: fbamg is holding a reference to this, making
  // the reference circular.
  Ref<FockBuildAMG> fbamg;
  if (fockbuildmatrixtype_ == "distributed"
      ||fockbuildmatrixtype_ == "prefetched_distributed") {
      fbamg_ = new FockBuildAMG(msg, thr, this);
    }

  for (int i=0; i<nlocks_; i++) {
      locks_[i] = thr->new_lock();
    }

  for (int i=0; i<nfmat; i++) {
      jmats_[i] = fockbuildmatrix(i,"J",msg,fbamg_);
      kmats_[i] = fockbuildmatrix(i,"K",msg,fbamg_);
    }

  for (int i=0; i<npmat; i++) {
      pmats_[i] = fockbuildmatrix(i,"P",msg,fbamg_);
    }
}


FockBuildMatrix *
GenericFockContribution::fockbuildmatrix(int matrix,
                                         const std::string &type,
                                         const Ref<MessageGrp> &msg,
                                         const Ref<FockBuildAMG> &fbamg)
{
  if (fockbuildmatrixtype_ == "replicated") {
      return new ReplFockBuildMatrix(msg);
    }
  else if (fockbuildmatrixtype_ == "prefetched_distributed"
           ||fockbuildmatrixtype_ == "distributed") {
      bool prefetch_blocks = (fockbuildmatrixtype_ == "prefetched_distributed");
      Ref<FockBuildAM> fbam;
      bool need_read, need_write;
      if (type == "J") {
          need_read = false;
          need_write = true;
          fbam = new FockBuildAMAccumJ(matrix);
        }
      else if (type == "K") {
          need_read = false;
          need_write = true;
          fbam = new FockBuildAMAccumK(matrix);
        }
      else if (type == "P") {
          need_read = true;
          need_write = false;
          fbam = new FockBuildAMGetP(matrix);
        }
      else {
          std::cout << "bad type: " << type
                    << std::endl;
          abort();
        }
      Ref<FockBuildOp> fockbuildop
          = new FockBuildOp(0, need_read, need_write, use_shell_blocks_,
                            fbamg, fbam);
      return new DistFockBuildMatrix(prefetch_blocks, fockbuildop);
    }
  else {
      throw ProgrammingError("GenericFockContribution::fockbuildmatrix: "
                                 "bad fockbuildmatrixtype",
                                 __FILE__, __LINE__);
    }
}

void
GenericFockContribution::set_fmat(int i, const RefSCMatrix &m)
{
  if (f_b1_equiv_f_b2) {
      throw ProgrammingError("set_fmat: rect but bases equiv",
                                 __FILE__, __LINE__);
    }
  jmats_[i]->scmat_to_data(m, f_b1_, f_b2_, false);
  kmats_[i] = jmats_[i];
  k_is_j_[i] = true;
}

void
GenericFockContribution::set_fmat(int i, const RefSymmSCMatrix &m)
{
  if (!f_b1_equiv_f_b2) {
      throw ProgrammingError("set_fmat: symm but bases not equiv",
                                 __FILE__, __LINE__);
    }
  jmats_[i]->scmat_to_data(m, f_b1_, false);
  kmats_[i] = jmats_[i];
  k_is_j_[i] = true;
}

void
GenericFockContribution::set_jmat(int i, const RefSCMatrix &m)
{
  if (f_b1_equiv_f_b2) {
      throw ProgrammingError("set_jmat: rect but bases equiv",
                                 __FILE__, __LINE__);
    }
  jmats_[i]->scmat_to_data(m, f_b1_, f_b2_, false);
  k_is_j_[i] = false;
}

void
GenericFockContribution::set_jmat(int i, const RefSymmSCMatrix &m)
{
  if (!f_b1_equiv_f_b2) {
      throw ProgrammingError("set_jmat: symm but bases not equiv",
                                 __FILE__, __LINE__);
    }
  jmats_[i]->scmat_to_data(m, f_b1_, false);
  k_is_j_[i] = false;
}

void
GenericFockContribution::set_kmat(int i, const RefSCMatrix &m)
{
  if (f_b1_equiv_f_b2) {
      throw ProgrammingError("set_kmat: rect but bases equiv",
                                 __FILE__, __LINE__);
    }
  kmats_[i]->scmat_to_data(m, f_b1_, f_b2_, false);
  k_is_j_[i] = false;
}

void
GenericFockContribution::set_kmat(int i, const RefSymmSCMatrix &m)
{
  if (!f_b1_equiv_f_b2) {
      throw ProgrammingError("set_kmat: symm but bases not equiv",
                                 __FILE__, __LINE__);
    }
  kmats_[i]->scmat_to_data(m, f_b1_, false);
  k_is_j_[i] = false;
}

void
GenericFockContribution::set_pmat(int i, const RefSymmSCMatrix &m)
{
  pmats_[i]->scmat_to_data(m, p_b_, true);
}

void
GenericFockContribution::accum_remote(const Ref<MessageGrp> &msg)
{
  for (int i=0; i<nfmat_; i++) {
      jmats_[i]->accum_remote(msg);
      if (!k_is_j_[i]) {
          kmats_[i]->accum_remote(msg);
        }
    }
  msg->sum(nint_);
}

void
GenericFockContribution::update()
{
  for (int i=0; i<nfmat_; i++) {
      jmats_[i]->fix_diagonal_blocks();
      jmats_[i]->data_to_scmat();
      if (!k_is_j_[i]) {
          kmats_[i]->fix_diagonal_blocks();
          kmats_[i]->data_to_scmat();
        }
    }
}

void
GenericFockContribution::copy_matrices(int unique_id)
{
  for (int i=0; i<nfmat_; i++) {
      jmats_[i] = jmats_[i]->copy(unique_id);
      if (k_is_j_[i]) {
          kmats_[i] = jmats_[i];
        }
      else {
          kmats_[i] = kmats_[i]->copy(unique_id);
        }
    }

  for (int i=0; i<npmat_; i++) {
      pmats_[i] = pmats_[i]->copy(unique_id);
    }
}

void
GenericFockContribution::accum(const Ref<FockContribution> &c_abs)
{
  // This is not needed since the storage for J & K
  // is shared among threads and updates are locked.
//   GenericFockContribution *c
//       = dynamic_cast<GenericFockContribution*>(c_abs.pointer());
//   if (nfmat_ != c->nfmat_
//       || npmat_ != c->npmat_) {
//       throw std::invalid_argument("merging incompatible fock contributions");
//     }
//   for (int i=0; i<nfmat_; i++) {
//       jmats_[i]->accum(c->jmats_[i]);
//       if (!k_is_j_[i]) {
//           kmats_[i]->accum(c->kmats_[i]);
//         }
//     }
//   for (int i=0; i<npmat_; i++) {
//       pmats_[i]->accum(c->pmats_[i]);
//     }
  nint_ += c_abs->nint();
}

GenericFockContribution::~GenericFockContribution()
{
}

signed char *
GenericFockContribution::compute_pmax() const
{
  int nb1 = p_b_->basisdim()->blocks()->nblock();
  int n12 = (nb1*(nb1+1))/2;
  signed char *pmax = new signed char[n12];
  for (int i=0; i<n12; i++) pmax[i] = SCHAR_MIN;

  for (int i=0; i<pmats_.size(); i++) {
      pmax_contrib(pmats_[i],pmax);
    }

  pmats_[0]->messagegrp()->max(pmax,n12);

  return pmax;
}

void
GenericFockContribution::pmax_contrib(const Ref<FockBuildMatrix> &mat,
                                      signed char *pmax) const
{
  double l2inv = 1.0/log(2.0);
  double tol = pow(2.0,-126.0);

  Ref<SCBlockInfo> bi1 = p_b_->basisdim()->blocks();

  int ish, jsh, ij;
  for (ish=ij=0; ish < bi1->nblock(); ish++) {
    int ni = bi1->size(ish);

    int jsh_fence;
    for (jsh=0; jsh <= ish; jsh++,ij++) {

      if (!mat->shell_block_is_owner(ish,jsh)) continue;

      int nj = bi1->size(jsh);

      double maxp=0, tmp;
      double *pmat_block = mat->shell_block(ish,jsh);
      int nij = ni*nj;

      for (int k=0; k<nij; k++) {
          if ((tmp=fabs(pmat_block[k])) > maxp) {
              maxp=tmp;
            }
      }

      if (maxp <= tol) maxp=tol;

      long power = long(ceil(log(maxp)*l2inv));
      signed char pmaxij;
      if (power < SCHAR_MIN) pmaxij = SCHAR_MIN;
      else if (power > SCHAR_MAX) pmaxij = SCHAR_MAX;
      else pmaxij = (signed char) power;

      if (pmaxij > pmax[ij]) pmax[ij] = pmaxij;
    }
  }
}

void
GenericFockContribution::activate()
{
  if (fbamg_) {
//       std::cout << "GenericFockContribution::activate(): activating"
//                 << std::endl;
      fbamg_->activate();
    }
//   else {
//       std::cout << "GenericFockContribution::activate(): fbamg_ null"
//                 << std::endl;
//     }
}

void
GenericFockContribution::sync()
{
  if (fbamg_) fbamg_->sync();
}

void
GenericFockContribution::deactivate()
{
  if (fbamg_) fbamg_->deactivate();
}

void
GenericFockContribution::flush()
{
  for (int i=0; i<nfmat_; i++) {
      jmats_[i]->flush();
      kmats_[i]->flush();
    }
  for (int i=0; i<npmat_; i++) {
      pmats_[i]->flush();
    }
}

void
GenericFockContribution::prefetch_blocks(int I, int J,
                                         int ifetch, int nfetch)
{
  int newnfetch = nfetch * npmat_;
  for (int i=0; i<npmat_; i++) {
      pmats_[i]->prefetch_block(I,J,ifetch*npmat_+i,newnfetch);
    }
}

void
GenericFockContribution::finish_prefetch_blocks()
{
  for (int i=0; i<npmat_; i++) {
      pmats_[i]->finish_prefetch_block();
    }
}

void
GenericFockContribution::set_fockblocks(const Ref<FockBlocks> &blocks_f1,
                                        const Ref<FockBlocks> &blocks_f2,
                                        const Ref<FockBlocks> &blocks_p)
{
  for (int i=0; i<nfmat_; i++) {
      jmats_[i]->set_fockblocks(blocks_f1,blocks_f2);
      if (!k_is_j_[i]) {
          kmats_[i]->set_fockblocks(blocks_f1,blocks_f2);
        }
    }
  for (int i=0; i<npmat_; i++) {
      pmats_[i]->set_fockblocks(blocks_p,blocks_p);
    }
}

/////////////////////////////////////////////////////////////////
// FockBuildThread

FockBuildThread::FockBuildThread(const Ref<FockDistribution> &fockdist,
                                 const Ref<MessageGrp> &msg,
                                 int nthread,
                                 int threadnum,
                                 bool prefetch_blocks,
                                 const Ref<ThreadLock> &lock,
                                 const Ref<Integral> &integral,
                                 bool compute_J,
                                 bool compute_K,
                                 double coef_K):
  fockdist_(fockdist),
  msg_(msg),
  nthread_(nthread),
  threadnum_(threadnum),
  prefetch_blocks_(prefetch_blocks),
  lock_(lock),
  integral_(integral),
  accuracy_(DBL_EPSILON),
  pmax_(0),
  compute_J_(compute_J),
  compute_K_(compute_K),
  coef_K_(coef_K)
{
  timer_ = new RegionTimer("FockBuildThreads",1,1);
}

/////////////////////////////////////////////////////////////////
// FockBuildThread_F11_P11

FockBuildThread_F11_P11::FockBuildThread_F11_P11(
    const Ref<FockDistribution> &fockdist,
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
    const Ref<FockBlocks> &blocks3,
    bool compute_J,
    bool compute_K,
    double coef_K
    ):
  FockBuildThread(fockdist,msg,nthread,threadnum,prefetch_blocks,lock,integral,
                  compute_J,compute_K,coef_K),
  pl_(pl),
  basis_(basis1),
  blocks_(blocks1)
{
  integral_->set_basis(basis_);
  eri_ = integral_->electron_repulsion();
  eri_->set_redundant(1);
  eri_->set_integral_storage(integral_->storage_unused()/nthread);
}

void
FockBuildThread_F11_P11::prefetch_blocks(const Ref<FockDist> &dist,
                                         int iblock, int jblock, int kblock, int lblock)
{
  int ifetch = 0;
  int nfetch = 8;
  if (iblock>=jblock) contrib_->prefetch_blocks(iblock,jblock,ifetch++,nfetch);
  if (iblock>=kblock) contrib_->prefetch_blocks(iblock,kblock,ifetch++,nfetch);
  if (iblock>=lblock) contrib_->prefetch_blocks(iblock,lblock,ifetch++,nfetch);
  if (jblock>=kblock) contrib_->prefetch_blocks(jblock,kblock,ifetch++,nfetch);
  if (jblock>=lblock) contrib_->prefetch_blocks(jblock,lblock,ifetch++,nfetch);
  if (kblock>=jblock) contrib_->prefetch_blocks(kblock,jblock,ifetch++,nfetch);
  if (kblock>=lblock) contrib_->prefetch_blocks(kblock,lblock,ifetch++,nfetch);
  if (lblock>=jblock) contrib_->prefetch_blocks(lblock,jblock,ifetch++,nfetch);
}

void
FockBuildThread_F11_P11::run()
{
  int l2tol = (int) (log(accuracy_)/log(2.0));

  Ref<FockDist> dist
      = fockdist_->fockdist(basis_, blocks_, pl_, msg_,
                            nthread_, threadnum_,
                            pmax_, eri_, l2tol);

  const double *buf = eri_->buffer();

  if (!dist->fixed_integral_map()
      || !fockdist_->cache_integrals()) {
      eri_->set_integral_storage(0);
    }


  int iblock, jblock, kblock, lblock;
  int iblock_next, jblock_next, kblock_next, lblock_next;
  bool got_blocks;
  got_blocks = dist->get_blocks(iblock_next,jblock_next,kblock_next,lblock_next);
  if (prefetch_blocks_ && got_blocks) {
#if DEBUG_DIST
      lock_->lock();
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ":" << threadnum_
                << ": prefetch_blocks for "
                << iblock_next
                << ", " << jblock_next
                << ", " << kblock_next
                << ", " << lblock_next
                << std::endl;
      lock_->unlock();
#endif
      prefetch_blocks(dist,
                      iblock_next,jblock_next,kblock_next,lblock_next);
    }
  while (got_blocks) {
      iblock = iblock_next;
      jblock = jblock_next;
      kblock = kblock_next;
      lblock = lblock_next;

      got_blocks = dist->get_blocks(iblock_next,jblock_next,kblock_next,lblock_next);

      if (prefetch_blocks_) {
#if DEBUG_DIST
          lock_->lock();
          std::cout << MessageGrp::get_default_messagegrp()->me()
                    << ":" << threadnum_
                    << ": finish_prefetch_blocks"
                    << std::endl;
          lock_->unlock();
#endif
          contrib_->finish_prefetch_blocks();
          if (got_blocks) {
#if DEBUG_DIST
              lock_->lock();
              std::cout << MessageGrp::get_default_messagegrp()->me()
                        << ":" << threadnum_
                        << ": prefetch_blocks for "
                        << iblock_next
                        << ", " << jblock_next
                        << ", " << kblock_next
                        << ", " << lblock_next
                        << std::endl;
              lock_->unlock();
#endif
              prefetch_blocks(dist,
                              iblock_next,jblock_next,kblock_next,lblock_next);
            }
        }

      int ibegin = dist->begin(iblock);
      int jbegin = dist->begin(jblock);
      int kbegin = dist->begin(kblock);
      int lbegin = dist->begin(lblock);
      int iend = dist->end(iblock);
      int jend = dist->end(jblock);
      int kend = dist->end(kblock);
      int lend = dist->end(lblock);

#if DEBUG || DEBUG_DIST
      lock_->lock();
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ":" << threadnum_
                << ": working on blocks: "
                << std::setw(2) << iblock
                << ", " << std::setw(2) << jblock
                << ", " << std::setw(2) << kblock
                << ", " << std::setw(2) << lblock
                << std::endl;
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ":" << threadnum_
                << ": begin_blocks: "
                << std::setw(2) << ibegin
                << ", " << std::setw(2) << jbegin
                << ", " << std::setw(2) << kbegin
                << ", " << std::setw(2) << lbegin
                << std::endl;
      std::cout << MessageGrp::get_default_messagegrp()->me()
                << ":" << threadnum_
                << ": end_blocks: "
                << std::setw(2) << iend
                << ", " << std::setw(2) << jend
                << ", " << std::setw(2) << kend
                << ", " << std::setw(2) << lend
                << std::endl;
      lock_->unlock();
#endif

      for (int i=ibegin; i<iend; i++) {
          if (!pl_->in_p1(i)) continue;
          int ni=basis_->shell(i).nfunction();
          for (int j=jbegin; j<jend && j<=i; j++) {
              int oij = can_sym_offset(i,j);
              if (!pl_->in_p2(oij)) continue;
              int nj=basis_->shell(j).nfunction();
              int pmaxij = pmax_[oij];
              for (int k=kbegin; k<kend && k <= i; k++) {
                  int nk=basis_->shell(k).nfunction();
//                   int pmaxijk=pmaxij, ptmp;
//                   if ((ptmp=pmax_[can_sym_offset(i,k)]-1) > pmaxijk) pmaxijk=ptmp;
//                   if ((ptmp=pmax_[gen_sym_offset(j,k)]-1) > pmaxijk) pmaxijk=ptmp;

                  int pmaxik = pmax_[can_sym_offset(i,k)];
                  int pmaxjk = pmax_[gen_sym_offset(j,k)];
                  int okl = can_sym_offset(k,lbegin);
                  for (int l=lbegin; l < lend && l <= (k==i?j:k); l++,okl++) {

                      int qijkl = pl_->in_p4(oij,okl,i,j,k,l);
                      if (!qijkl) continue;

                      double jfac = qijkl;
                      double kfac = qijkl*coef_K_;

//                       int pmaxijkl = pmaxijk;
//                       if ((ptmp=pmax_[okl]) > pmaxijkl) pmaxijkl=ptmp;
//                       if ((ptmp=pmax_[can_sym_offset(i,l)]-1) > pmaxijkl) pmaxijkl=ptmp;
//                       if ((ptmp=pmax_[gen_sym_offset(j,l)]-1) > pmaxijkl) pmaxijkl=ptmp;

                      int pmaxkl = pmax_[okl];
                      int pmaxil = pmax_[can_sym_offset(i,l)];
                      int pmaxjl = pmax_[gen_sym_offset(j,l)];

                      int pmax_J = (pmaxij>pmaxkl?pmaxij:pmaxkl);
                      int pmax_K = pmaxik;
                      if (pmaxil>pmax_K) pmax_K = pmaxil;
                      if (pmaxjk>pmax_K) pmax_K = pmaxjk;
                      if (pmaxjl>pmax_K) pmax_K = pmaxjl;

#if DEBUG || DEBUG_DIST
                      lock_->lock();
                      std::cout << MessageGrp::get_default_messagegrp()->me()
                                << ":" << threadnum_
                                << ": shell: "
                                << std::setw(2) << i
                                << ", " << std::setw(2) << j
                                << ", " << std::setw(2) << k
                                << ", " << std::setw(2) << l
                                << std::endl;
                      lock_->unlock();
#endif

#if SCF_USE_BOUNDS
                      int intmax = eri_->log2_shell_bound(i,j,k,l);
                      int max_J = intmax + pmax_J;
                      int max_K = intmax + pmax_K - 1;
                      bool doJ = compute_J_ && (max_J >= l2tol);
                      bool doK = compute_K_ && (max_K >= l2tol);
                      if (!doJ && !doK) continue;
#else
                      bool doJ = compute_J_, doK = compute_K_;
#endif

#if DETAILED_TIMINGS
                      Timer tim(timer_,"compute_shell");
#endif
                      eri_->compute_shell(i,j,k,l);
#if DETAILED_TIMINGS
                      tim.exit();
#endif

                      int e12 = (i==j);
                      int e34 = (k==l);
                      int e13e24 = (i==k) && (j==l);
                      int e_any = e12||e34||e13e24;
                      int nl=basis_->shell(l).nfunction();

#if DETAILED_TIMINGS
                      tim.enter("contribs");
#endif
                      if (e12) {
                          if (e34) {
                              if (e13e24) {
                                  // e12 e34 e13e24
                                  if (doJ) contrib_->contrib_e_J(jfac, i, j, k, l,
                                                                 ni, nj, nk, nl, buf);
                                  if (doK) contrib_->contrib_e_K(kfac, i, j, k, l,
                                                                 ni, nj, nk, nl, buf);
                                }
                              else {
                                  // e12 e34
                                  if (doJ) contrib_->contrib_p13p24_J(jfac, i, j, k, l,
                                                                      ni, nj, nk, nl, buf);
                                  if (doK) contrib_->contrib_p13p24_K(kfac, i, j, k, l,
                                                                      ni, nj, nk, nl, buf);
                                }
                            }
                          else {
                              // e12
                              if (doJ) contrib_->contrib_p34_p13p24_J(jfac, i, j, k, l,
                                                                      ni, nj, nk, nl, buf);
                              if (doK) contrib_->contrib_p34_p13p24_K(kfac, i, j, k, l,
                                                                      ni, nj, nk, nl, buf);
                            }
                        }
                      else if (e34) {
                          // e34
                          if (doJ) contrib_->contrib_p12_p13p24_J(jfac, i, j, k, l,
                                                                  ni, nj, nk, nl, buf);
                          if (doK) contrib_->contrib_p12_p13p24_K(kfac, i, j, k, l,
                                                                  ni, nj, nk, nl, buf);
                        }
                      else if (e13e24) {
                          // e13e24
                          if (doJ) contrib_->contrib_p12_p34_J(jfac, i, j, k, l,
                                                               ni, nj, nk, nl, buf);
                          if (doK) contrib_->contrib_p12_p34_K(kfac, i, j, k, l,
                                                               ni, nj, nk, nl, buf);
                        }
                      else {
                          // no equivalent indices
                          if (doJ) contrib_->contrib_all_J(jfac, i, j, k, l,
                                                           ni, nj, nk, nl, buf);
                          if (doK) contrib_->contrib_all_K(kfac, i, j, k, l,
                                                           ni, nj, nk, nl, buf);
                        }
#if DETAILED_TIMINGS
                      tim.exit();
#endif

//                       std::cout << "contrib_->nint() before = " << contrib_->nint()
//                                 << std::endl;
                      contrib_->nint() += (double) ni*nj*nk*nl;
//                       std::cout << "contrib_->nint() after = " << contrib_->nint()
//                                 << std::endl;
                    }
                }
            }
        }
#if DETAILED_TIMINGS
      Timer flush_time(timer_,"flush");
#endif
      contrib_->flush();
    }
}

/////////////////////////////////////////////////////////////////
// FockBuildThread_F12_P33

FockBuildThread_F12_P33::FockBuildThread_F12_P33(
    const Ref<FockDistribution> &fockdist,
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
    const Ref<FockBlocks> &blocks3,
    bool compute_J,
    bool compute_K,
    double coef_K
    ):
  FockBuildThread(fockdist,msg,nthread,threadnum,prefetch_blocks,lock,integral,
                  compute_J,compute_K,coef_K),
  pl_(pl),
  basis1_(basis1),
  basis2_(basis2),
  basis3_(basis3)
{
  integral_->set_basis(basis1_,basis2_,basis3_,basis3_);
  eri_J_ = integral_->electron_repulsion();
  eri_J_->set_redundant(1);
  eri_J_->set_integral_storage(integral_->storage_unused()/nthread/2);

  integral_->set_basis(basis1_,basis3_,basis2_,basis3_);
  eri_K_ = integral_->electron_repulsion();
  eri_K_->set_redundant(1);
  eri_K_->set_integral_storage(integral_->storage_unused()/nthread/2);
}

void
FockBuildThread_F12_P33::run()
{
  if (compute_J_) run_J();
  if (compute_K_) run_K();
}

void
FockBuildThread_F12_P33::run_J()
{
  bool I_equiv_J = basis1_->equiv(basis2_);

  int me=msg_->me();
  int nproc = msg_->n();
  sc_int_least64_t thread_index=0;
  sc_int_least64_t task_index=0;

  const double *buf = eri_J_->buffer();
  int nshell1 = basis1_->nshell();
  int nshell2 = basis2_->nshell();
  int nshell3 = basis3_->nshell();
  for (int I=0; I<nshell1; I++) {
      int nI = basis1_->shell(I).nfunction();
      int Jfence = I_equiv_J?I+1:nshell2;
      for (int J=0; J<Jfence; J++) {
          int nJ = basis2_->shell(J).nfunction();
          for (int K=0; K<nshell3; K++, task_index++) {
              if (task_index%nproc != me)
                  continue;

              thread_index++;
              if (thread_index % nthread_ != threadnum_)
                  continue;

              int nK = basis3_->shell(K).nfunction();
              for (int L=0; L<=K; L++) {
                  int nL = basis3_->shell(L).nfunction();
                  eri_J_->compute_shell(I,J,K,L);
                  // I, J permutations are ignored since just
                  // the lower triangle of J is needed if I and
                  // J are from the same basis set
                  if (K != L) {
                      contrib_->contrib_p34_J(1.0, I, J, K, L,
                                              nI, nJ, nK, nL, buf);
                    }
                  else {
                      contrib_->contrib_e_J(1.0, I, J, K, L,
                                            nI, nJ, nK, nL, buf);
                    }

                   //                      std::cout << "contrib_->nint() before = " << contrib_->nint()
                   //                                << std::endl;
                                        contrib_->nint() += (double) nI*nJ*nK*nL;
                   //                      std::cout << "contrib_->nint() after = " << contrib_->nint()
                   //                                << std::endl;

                }
            }
        }
    }
}

void
FockBuildThread_F12_P33::run_K()
{
  bool I_equiv_K = basis1_->equiv(basis2_);

  int me=msg_->me();
  int nproc = msg_->n();
  sc_int_least64_t thread_index=0;
  sc_int_least64_t task_index=0;

  const double *buf = eri_K_->buffer();
  int nshell1 = basis1_->nshell();
  int nshell2 = basis2_->nshell();
  int nshell3 = basis3_->nshell();
  for (int I=0; I<nshell1; I++) {
      int nI = basis1_->shell(I).nfunction();
      for (int J=0; J<nshell3; J++) {
          int nJ = basis3_->shell(J).nfunction();
          int Kfence = I_equiv_K?I+1:nshell2;
          for (int K=0; K<Kfence; K++) {
              if (task_index%nproc != me)
                  continue;

              thread_index++;
              if (thread_index % nthread_ != threadnum_)
                  continue;

              int nK = basis2_->shell(K).nfunction();
              for (int L=0; L<nshell3; L++) {
                  int nL = basis3_->shell(L).nfunction();
                  eri_K_->compute_shell(I,J,K,L);
                  // I, K permutations are ignored since just
                  // the lower triangle of K is needed if I and
                  // K are from the same basis set
                  contrib_->contrib_e_K(coef_K_, I, J, K, L,
                                        nI, nJ, nK, nL, buf);

                  //                       std::cout << "contrib_->nint() before = " << contrib_->nint()
                  //                                 << std::endl;
                                        contrib_->nint() += (double) nI*nJ*nK*nL;
                  //                       std::cout << "contrib_->nint() after = " << contrib_->nint()
                  //                                 << std::endl;

                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////
// FockBuild

FockBuild::FockBuild(const Ref<FockDistribution> &fockdist,
                     const Ref<FockContribution> &contrib,
                     bool prefetch_blocks,
                     const Ref<GaussianBasisSet> &b_f1,
                     const Ref<GaussianBasisSet> &b_f2,
                     const Ref<GaussianBasisSet> &b_p,
                     const Ref<MessageGrp> &msg,
                     const Ref<ThreadGrp> &thr,
                     const Ref<Integral> &integral):
  fockdist_(fockdist),
  contrib_(contrib),
  prefetch_blocks_(prefetch_blocks),
  accuracy_(DBL_EPSILON),
  b_f1_(b_f1),
  b_f2_(b_f2),
  b_p_(b_p),
  msg_(msg),
  thr_(thr),
  integral_(integral),
  thread_(0),
  compute_J_(true),
  compute_K_(true),
  coef_K_(1.0)
{
  if (b_f2_.null()) b_f2_ = b_f1_;
  if (b_p_.null()) b_p_ = b_f1_;

  integral_->set_basis(b_p_);
  pl_ = integral_->petite_list();

//   for (int ishell=0; ishell<b_p_->nshell(); ishell++) {
//       int icenter = b_p_->shell_to_center(ishell);
//       std::cout << " in_p1(" << std::setw(2) << ishell
//                 << ") = " << pl_->in_p1(ishell)
//                 << " center " << std::setw(2) <<  icenter
//                 << std::endl;
//     }

//   for (int ishell=0; ishell<b_p_->nshell(); ishell++) {
//       int icenter = b_p_->shell_to_center(ishell);
//       std::cout << " in_p1(" << std::setw(2) << ishell
//                 << ") = " << pl_->in_p1(ishell)
//                 << " center " << icenter
//                 << std::endl;
//       if (pl_->in_p1(ishell)) {
//           for (int jshell=0; jshell<=ishell; jshell++) {
//               int jcenter = b_p_->shell_to_center(jshell);
//               std::cout << "  lambda(" << std::setw(2) << ishell
//                         << "," << std::setw(2) << jshell << ") = "
//                         << pl_->lambda(ishell,jshell)
//                         << " centers " << std::setw(2)
//                         << std::setw(2) << icenter << ", "
//                         << std::setw(2) << jcenter
//                         << std::endl;
//             }
//         }
//     }

  fb_f1_ = fockdist_->fockblocks(b_f1_);
  fb_f2_ = fockdist_->fockblocks(b_f2_);
  fb_p_  = fockdist_->fockblocks(b_p_);

  contrib_->set_fockblocks(fb_f1_,fb_f2_,fb_p_);

  init_threads();
}

FockContribution::FockContribution(const FockContribution&)
{
  nint_ = 0.0;
}

FockBuild::~FockBuild()
{
  done_threads();
}

template <class T>
FockBuildThread *
create_FockBuildThread(const Ref<FockDistribution> &fockdist,
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
                       const Ref<FockBlocks> &blocks3,
                       bool compute_J,
                       bool compute_K,
                       double coef_K)
{
  if (coef_K == 0.0) compute_K = false;
  return new T(fockdist,msg,nthread,threadnum,prefetch_blocks,lock,integral,
               pl,basis1,basis2,basis3,blocks1,blocks2,blocks3,
               compute_J,compute_K,coef_K);
}

void
FockBuild::build()
{
  int nthread = thr_->nthread();
  std::vector<Ref<FockContribution> > contribs(nthread);

  contrib_->activate();

  signed char *pmax = contrib_->compute_pmax();

  Timer tim_fockbuild_only("build_only");

  contrib_->nint() = 0;

  for (int i=0; i<nthread; i++) {
      if (i==0) contribs[i] = contrib_;
      else {
          contribs[i] = contrib_->clone();
          contribs[i]->copy_matrices(i);
        }
      thread_[i]->set_accuracy(accuracy_);
      thread_[i]->set_compute_J(compute_J_);
      thread_[i]->set_compute_K(compute_K_);
      thread_[i]->set_coef_K(coef_K_);
      thread_[i]->set_pmax(pmax);
      thread_[i]->set_contrib(contribs[i]);
      thr_->add_thread(i, thread_[i]);
    }
  thr_->start_threads();
  thr_->wait_threads();
  contrib_->deactivate();
  delete[] pmax;
  for (int i=1; i<nthread; i++) {
      contrib_->accum(contribs[i]);
    }
  for (int i=0; i<nthread; i++) {
      Ref<RegionTimer> deftimer = RegionTimer::default_regiontimer();
      if (deftimer) deftimer->merge(thread_[i]->get_timer());
      thread_[i]->get_timer()->reset();
    }
  contrib_->accum_remote(msg_);

  tim_fockbuild_only.exit();

  contrib_->update();
}

void
FockBuild::done_threads()
{
  int nthread = thr_->nthread();
  for (int i=0; i<nthread; i++) {
      delete thread_[i];
    }
  delete[] thread_;
}

void
FockBuild::init_threads(FBT_CTOR ctor)
{
  int nthread = thr_->nthread();
  thread_ = new FockBuildThread*[nthread];
  Ref<ThreadLock> lock = thr_->new_lock();
  for (int i=0; i<nthread; i++) {
      thread_[i] = (*ctor)(fockdist_,msg_,nthread,i,
                           prefetch_blocks_,
                           lock,integral_,pl_,
                           b_f1_,b_f2_,b_p_,
                           fb_f1_,fb_f2_,fb_p_,
                           compute_J_, compute_K_, coef_K_);
    }
}

void
FockBuild::init_threads()
{
  if (b_f1_->equiv(b_f2_)) {
      if (b_f1_->equiv(b_p_)) {
          init_threads(create_FockBuildThread<FockBuildThread_F11_P11>);
          return;
        }
      else {
        }
    }
  else if (b_f1_->equiv(b_p_)) {
    }
  else if (b_f2_->equiv(b_p_)) {
    }
  // use the general code if this case not handled above
  init_threads(create_FockBuildThread<FockBuildThread_F12_P33>);
  return;
}

}


// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
