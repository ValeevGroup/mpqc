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

#ifdef __GNUC__
#pragma implementation
#endif

#include <memory>
#include <stdexcept>
#include <vector>

#include <util/misc/scint.h>
#include <util/misc/autovec.h>
#include <util/class/scexception.h>
#include <math/scmat/elemop.h>
#include <math/scmat/blkiter.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/fockbuild.h>

#undef DEBUG
#define DEBUG 0

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
      rowbi_(rowbi),
      colbi_(colbi),
      blocks_(blocks),
      ndata_(ndata),
      defer_collect_(false) {
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
    int has_side_effects() { return data_to_mat_; }
    bool threadsafe() { return true; }
    bool cloneable() { return false; }
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

FockBuildMatrix::FockBuildMatrix():
  owns_data_(false),
  blockpointers_(0),
  nI_(0),
  nJ_(0),
  ndata_(0)
{
}
    
FockBuildMatrix::~FockBuildMatrix()
{
  clear();
}

void
FockBuildMatrix::make_reference(const FockBuildMatrix &fbm)
{
  owns_data_ = false;
  blockpointers_ = fbm.blockpointers_;
  rectmat_ = fbm.rectmat_;
  symmmat_ = fbm.symmmat_;
  nI_ = fbm.nI_;
  nJ_ = fbm.nJ_;
  ndata_ = fbm.ndata_;
}

void
FockBuildMatrix::clear()
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

void
FockBuildMatrix::copy_data()
{
  if (owns_data_) return;
  int nIJ = nblock();
  auto_vec<double*> tmp_blockpointers(new double*[nIJ]);
  auto_vec<double> tmp_data(new double[ndata_]);
  double *base1 = blockpointers_[0];
  double *base2 = tmp_data.get();
  for (int ij=0; ij<nIJ; ij++) {
      tmp_blockpointers[ij] = base2 + (blockpointers_[ij] - base1);
    }
  blockpointers_ = tmp_blockpointers.release();
  owns_data_ = true;
}

void
FockBuildMatrix::fix_diagonal_blocks() const
{
  if (!symmetric()) return;
  Ref<SCBlockInfo> bi = bs1_->basisdim()->blocks();
  for (int I=0; I<nI_; I++) {
      double *dat = block(I,I);
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
FockBuildMatrix::scmat_to_data(const Ref<SymmSCMatrix> &m,
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

  auto_vec<double*> tmp_blockpointers(new double*[nblock()]);
  auto_vec<double> tmp_data(new double[ndata_]);

  blockpointers_ = tmp_blockpointers.release();
  blockpointers_[0] = tmp_data.release();
  owns_data_ = true;

  double *current_data = blockpointers_[0];
  for (int I=0; I<nI_; I++) {
      for (int J=0; J<=I; J++) {
          blockpointers_[offset(I,J)] = current_data;
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
      memset(blockpointers_[0], '0', ndata_*sizeof(double));
    }
}

void
FockBuildMatrix::scmat_to_data(const Ref<SCMatrix> &m,
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
FockBuildMatrix::data_to_rectmat() const
{
  throw std::runtime_error("data_to_rectmat: not yet impl");
}

void
FockBuildMatrix::data_to_symmat() const
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
FockBuildMatrix::data_to_scmat() const
{
  if (symmmat_.nonnull()) {
      if (rectmat_.nonnull()) {
          throw std::runtime_error("data_to_scmat: both mats nonnull");
        }
      data_to_symmat();
    }
  else if (rectmat_.nonnull()) {
      data_to_rectmat();
    }
  else {
      throw std::runtime_error("data_to_scmat: both mats null");
    }
}

void
FockBuildMatrix::accum(const FockBuildMatrix &f)
{
  if (ndata_ != f.ndata_) {
      throw std::invalid_argument("incompatible FockBuildMatrix in accum");
    }
  double *src = f.blockpointers_[0];
  double *dst = blockpointers_[0];
  for (int i=0; i<ndata_; i++) src[i] += dst[i];
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
    Ref<GaussianBasisSet> &f_b1,
    Ref<GaussianBasisSet> &f_b2,
    Ref<GaussianBasisSet> &p_b1,
    Ref<GaussianBasisSet> &p_b2):
  nfmat_(nfmat),
  npmat_(npmat),
  jmats_(nfmat),
  kmats_(nfmat),
  k_is_j_(nfmat),
  pmats_(npmat),
  f_b1_(f_b1),
  f_b2_(f_b2),
  p_b1_(p_b1),
  p_b2_(p_b2)
{
  f_b1_equiv_f_b2 = f_b1_->equiv(f_b2);
  p_b1_equiv_p_b2 = p_b1_->equiv(p_b2);
}

void
GenericFockContribution::set_fmat(int i, const Ref<SCMatrix> &m)
{
  if (f_b1_equiv_f_b2) {
      throw sc::ProgrammingError("set_fmat: rect but bases equiv",
                                 __FILE__, __LINE__);
    }
  jmats_[i].scmat_to_data(m, f_b1_, f_b2_, false);
  kmats_[i] = jmats_[i];
  k_is_j_[i] = true;
}

void
GenericFockContribution::set_fmat(int i, const Ref<SymmSCMatrix> &m)
{
  if (!f_b1_equiv_f_b2) {
      throw sc::ProgrammingError("set_fmat: symm but bases not equiv",
                                 __FILE__, __LINE__);
    }
  jmats_[i].scmat_to_data(m, f_b1_, false);
  kmats_[i] = jmats_[i];
  k_is_j_[i] = true;
}

void
GenericFockContribution::set_jmat(int i, const Ref<SCMatrix> &m)
{
  if (f_b1_equiv_f_b2) {
      throw sc::ProgrammingError("set_jmat: rect but bases equiv",
                                 __FILE__, __LINE__);
    }
  jmats_[i].scmat_to_data(m, f_b1_, f_b2_, false);
  k_is_j_[i] = false;
}

void
GenericFockContribution::set_jmat(int i, const Ref<SymmSCMatrix> &m)
{
  if (!f_b1_equiv_f_b2) {
      throw sc::ProgrammingError("set_jmat: symm but bases not equiv",
                                 __FILE__, __LINE__);
    }
  jmats_[i].scmat_to_data(m, f_b1_, false);
  k_is_j_[i] = false;
}

void
GenericFockContribution::set_kmat(int i, const Ref<SCMatrix> &m)
{
  if (f_b1_equiv_f_b2) {
      throw sc::ProgrammingError("set_kmat: rect but bases equiv",
                                 __FILE__, __LINE__);
    }
  kmats_[i].scmat_to_data(m, f_b1_, f_b2_, false);
  k_is_j_[i] = false;
}

void
GenericFockContribution::set_kmat(int i, const Ref<SymmSCMatrix> &m)
{
  if (!f_b1_equiv_f_b2) {
      throw sc::ProgrammingError("set_kmat: symm but bases not equiv",
                                 __FILE__, __LINE__);
    }
  kmats_[i].scmat_to_data(m, f_b1_, false);
  k_is_j_[i] = false;
}

void
GenericFockContribution::set_pmat(int i, const Ref<SCMatrix> &m)
{
  if (p_b1_equiv_p_b2) {
      throw sc::ProgrammingError("set_pmat: rect but bases equiv",
                                 __FILE__, __LINE__);
    }
  pmats_[i].scmat_to_data(m, f_b1_, f_b2_, true);
}

void
GenericFockContribution::set_pmat(int i, const Ref<SymmSCMatrix> &m)
{
  if (!p_b1_equiv_p_b2) {
      throw sc::ProgrammingError("set_pmat: symm but bases not equiv",
                                 __FILE__, __LINE__);
    }
  pmats_[i].scmat_to_data(m, p_b1_, true);
}

void
GenericFockContribution::update()
{
  for (int i=0; i<nfmat_; i++) {
      jmats_[i].fix_diagonal_blocks();
      jmats_[i].data_to_scmat();
      if (!k_is_j_[i]) {
          kmats_[i].fix_diagonal_blocks();
          kmats_[i].data_to_scmat();
        }
    }
}

void
GenericFockContribution::copy()
{
  for (int i=0; i<nfmat_; i++) {
      jmats_[i].copy_data();
      if (!k_is_j_[i]) {
          kmats_[i].copy_data();
        }
    }
  for (int i=0; i<npmat_; i++) {
      pmats_[i].copy_data();
    }
}

void
GenericFockContribution::accum(const Ref<FockContribution> &c_abs)
{
  GenericFockContribution *c
      = dynamic_cast<GenericFockContribution*>(c_abs.pointer());
  if (nfmat_ != c->nfmat_
      || npmat_ != c->npmat_) {
      throw std::invalid_argument("merging incompatible fock contributions");
    }
  for (int i=0; i<nfmat_; i++) {
      jmats_[i].accum(c->jmats_[i]);
      if (!k_is_j_[i]) {
          kmats_[i].accum(c->kmats_[i]);
        }
    }
  for (int i=0; i<npmat_; i++) {
      pmats_[i].accum(c->pmats_[i]);
    }
  nint_ += c_abs->nint();
}

GenericFockContribution::~GenericFockContribution()
{
}

signed char *
GenericFockContribution::compute_pmax(bool symmetric) const
{
  int nb1 = p_b1_->basisdim()->blocks()->nblock();
  int n12;
  if (symmetric) {
      n12 = (nb1*(nb1+1))/2;
    }
  else {
      int nb2 = p_b2_->basisdim()->blocks()->nblock();
      n12 = nb1*nb2;
    }
  signed char *pmax = new signed char[n12];
  for (int i=0; i<n12; i++) pmax[i] = SCHAR_MIN;

  for (int i=0; i<pmats_.size(); i++) {
      pmax_contrib(pmats_[i],pmax,symmetric);
    }
  
  return pmax;
}

void
GenericFockContribution::pmax_contrib(const FockBuildMatrix &mat,
                                      signed char *pmax,
                                      bool symmetric) const
{
  double l2inv = 1.0/log(2.0);
  double tol = pow(2.0,-126.0);
  
  Ref<SCBlockInfo> bi1 = p_b1_->basisdim()->blocks();
  Ref<SCBlockInfo> bi2 = p_b2_->basisdim()->blocks();

  int ish, jsh, ij;
  for (ish=ij=0; ish < bi1->nblock(); ish++) {
    int ni = bi1->size(ish);

    int jsh_fence;
    if (symmetric) jsh_fence = ish+1;
    else           jsh_fence = bi2->nblock();
    for (jsh=0; jsh < jsh_fence; jsh++,ij++) {
      int nj = bi2->size(jsh);
      
      double maxp=0, tmp;
      double *pmat_block = mat.block(ish,jsh);
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

/////////////////////////////////////////////////////////////////
// FockBuildThread

FockBuildThread::FockBuildThread(const Ref<MessageGrp> &msg,
                                 int nthread,
                                 int threadnum,
                                 const Ref<FockContribution>&c,
                                 const Ref<ThreadLock> &lock,
                                 const Ref<Integral> &integral,
                                 double acc, const signed char *pmax):
  msg_(msg),
  nthread_(nthread),
  threadnum_(threadnum),
  contrib_(c),
  lock_(lock),
  integral_(integral),
  accuracy_(acc),
  pmax_(pmax)
{
}

/////////////////////////////////////////////////////////////////
// FockBuildThread_F11_P11

FockBuildThread_F11_P11::FockBuildThread_F11_P11(
    const Ref<MessageGrp> &msg,
    int nthread,
    int threadnum,
    const Ref<FockContribution>&c,
    const Ref<ThreadLock> &lock,
    const Ref<Integral> &integral,
    double acc, const signed char *pmax,
    const Ref<PetiteList> &pl,
    const Ref<GaussianBasisSet> &basis1,
    const Ref<GaussianBasisSet> &basis2,
    const Ref<GaussianBasisSet> &basis3,
    const Ref<GaussianBasisSet> &basis4
    ):
  FockBuildThread(msg,nthread,threadnum,c,lock,integral,acc,pmax),
  pl_(pl),
  basis_(basis1)
{
}

void
FockBuildThread_F11_P11::run()
{
  int tol = (int) (log(accuracy_)/log(2.0));
  int me=msg_->me();
  int nproc = msg_->n();
  
  GaussianBasisSet& gbs = *basis_;
  PetiteList& pl = *pl_;

  integral_->set_basis(basis_);
  Ref<TwoBodyInt> eri = integral_->electron_repulsion();

  eri->set_redundant(1);
  const double *buf = eri->buffer();

  sc_int_least64_t threadind=0;
  sc_int_least64_t ijklind=0;

  for (int i=0; i < gbs.nshell(); i++) {
      if (!pl.in_p1(i))
          continue;

      int ni=gbs(i).nfunction();
        
      for (int j=0; j <= i; j++) {
          int oij = can_sym_offset(i,j);

          if (!pl.in_p2(oij))
              continue;

          int nj=gbs(j).nfunction();
          int pmaxij = pmax_[oij];

          for (int k=0; k <= i; k++, ijklind++) {
              if (ijklind%nproc != me)
                  continue;

              threadind++;
              if (threadind % nthread_ != threadnum_)
                  continue;
            
              int nk=gbs(k).nfunction();

              int pmaxijk=pmaxij, ptmp;
              if ((ptmp=pmax_[can_sym_offset(i,k)]-1) > pmaxijk) pmaxijk=ptmp;
              if ((ptmp=pmax_[gen_sym_offset(j,k)]-1) > pmaxijk) pmaxijk=ptmp;
        
              int okl = can_sym_offset(k,0);
              for (int l=0; l <= (k==i?j:k); l++,okl++) {
                  int pmaxijkl = pmaxijk;
                  if ((ptmp=pmax_[okl]) > pmaxijkl) pmaxijkl=ptmp;
                  if ((ptmp=pmax_[can_sym_offset(i,l)]-1) > pmaxijkl) pmaxijkl=ptmp;
                  if ((ptmp=pmax_[gen_sym_offset(j,l)]-1) > pmaxijkl) pmaxijkl=ptmp;

                  int qijkl = pl.in_p4(oij,okl,i,j,k,l);
                  if (!qijkl)
                      continue;

#ifdef SCF_CHECK_BOUNDS
                  double intbound
                      = pow(2.0,double(eri->log2_shell_bound(i,j,k,l)));
                  double pbound
                      = pow(2.0,double(pmaxijkl));
                  intbound *= qijkl;
                  GBuild<T>::contribution.set_bound(intbound, pbound);
#else
#  ifndef SCF_DONT_USE_BOUNDS
                  if (eri->log2_shell_bound(i,j,k,l)+pmaxijkl < tol)
                      continue;
#  endif
#endif

                  eri->compute_shell(i,j,k,l);

                  int e12 = (i==j);
                  int e34 = (k==l);
                  int e13e24 = (i==k) && (j==l);
                  int e_any = e12||e34||e13e24;
                  int nl=gbs(l).nfunction();

                  if (e12) {
                      if (e34) {
                          if (e13e24) {
                              // e12 e34 e13e24
                              contrib_->contrib_e_J(qijkl, i, j, k, l,
                                                    ni, nj, nk, nl, buf);
                              contrib_->contrib_e_K(qijkl, i, j, k, l,
                                                    ni, nj, nk, nl, buf);
                            }
                          else {
                              // e12 e34
                              contrib_->contrib_p13p24_J(qijkl, i, j, k, l,
                                                         ni, nj, nk, nl, buf);
                              contrib_->contrib_p13p24_K(qijkl, i, j, k, l,
                                                         ni, nj, nk, nl, buf);
                            }
                        }
                      else {
                          // e12
                          contrib_->contrib_p34_p13p24_J(qijkl, i, j, k, l,
                                                         ni, nj, nk, nl, buf);
                          contrib_->contrib_p34_p13p24_K(qijkl, i, j, k, l,
                                                         ni, nj, nk, nl, buf);
                        }
                    }
                  else if (e34) {
                      // e34
                      contrib_->contrib_p12_p13p24_J(qijkl, i, j, k, l,
                                                     ni, nj, nk, nl, buf);
                      contrib_->contrib_p12_p13p24_K(qijkl, i, j, k, l,
                                                     ni, nj, nk, nl, buf);
                    }
                  else if (e13e24) {
                      // e13e24
                      contrib_->contrib_p12_p34_J(qijkl, i, j, k, l,
                                                ni, nj, nk, nl, buf);
                      contrib_->contrib_p12_p34_K(qijkl, i, j, k, l,
                                                 ni, nj, nk, nl, buf);
                    }
                  else {
                      // no equivalent indices
                      contrib_->contrib_all_J(qijkl, i, j, k, l,
                                              ni, nj, nk, nl, buf);
                      contrib_->contrib_all_K(qijkl, i, j, k, l,
                                              ni, nj, nk, nl, buf);
                    }

                  contrib_->nint() += (double) ni*nj*nk*nl;
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////
// FockBuildThread_F12_P34

FockBuildThread_F12_P34::FockBuildThread_F12_P34(
    const Ref<MessageGrp> &msg,
    int nthread,
    int threadnum,
    const Ref<FockContribution>&c,
    const Ref<ThreadLock> &lock,
    const Ref<Integral> &integral,
    double acc, const signed char *pmax,
    const Ref<PetiteList> &pl,
    const Ref<GaussianBasisSet> &basis1,
    const Ref<GaussianBasisSet> &basis2,
    const Ref<GaussianBasisSet> &basis3,
    const Ref<GaussianBasisSet> &basis4
    ):
  FockBuildThread(msg,nthread,threadnum,c,lock,integral,acc,pmax),
  pl_(pl),
  basis1_(basis1),
  basis2_(basis2),
  basis3_(basis3),
  basis4_(basis4)
{
}

void
FockBuildThread_F12_P34::run()
{
  run_J();
  run_K();
}

void
FockBuildThread_F12_P34::run_J()
{
  ThreadLockHolder eri_alloc_lock(lock_);
  integral_->set_basis(basis1_,basis2_,basis3_,basis4_);
  Ref<TwoBodyInt> eri = integral_->electron_repulsion();
  eri_alloc_lock.unlock();
  const double *buf = eri->buffer();
  int nshell1 = basis1_->nshell();
  int nshell2 = basis2_->nshell();
  int nshell3 = basis3_->nshell();
  int nshell4 = basis4_->nshell();
  for (int I=0; I<nshell1; I++) {
      int nI = basis1_->shell(I).nfunction();
      for (int J=0; J<nshell2; J++) {
          int nJ = basis2_->shell(J).nfunction();
          for (int K=0; K<nshell3; K++) {
              int nK = basis3_->shell(K).nfunction();
              for (int L=0; L<nshell4; L++) {
                  int nL = basis4_->shell(L).nfunction();
                  eri->compute_shell(I,J,K,L);
                  contrib_->contrib_e_J(1.0, I, J, K, L, nI, nJ, nK, nL, buf);
                }
            }
        }
    }
}

void
FockBuildThread_F12_P34::run_K()
{
  ThreadLockHolder eri_alloc_lock(lock_);
  integral_->set_basis(basis1_,basis3_,basis2_,basis4_);
  Ref<TwoBodyInt> eri = integral_->electron_repulsion();
  eri_alloc_lock.unlock();
  const double *buf = eri->buffer();
  int nshell1 = basis1_->nshell();
  int nshell2 = basis2_->nshell();
  int nshell3 = basis3_->nshell();
  int nshell4 = basis4_->nshell();
  for (int I=0; I<nshell1; I++) {
      int nI = basis1_->shell(I).nfunction();
      for (int J=0; J<nshell3; J++) {
          int nJ = basis3_->shell(J).nfunction();
          for (int K=0; K<nshell2; K++) {
              int nK = basis2_->shell(K).nfunction();
              for (int L=0; L<nshell4; L++) {
                  int nL = basis4_->shell(L).nfunction();
                  eri->compute_shell(I,J,K,L);
                  contrib_->contrib_e_K(1.0, I, J, K, L, nI, nJ, nK, nL, buf);
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////
// FockBuildThread_F11_P22

FockBuildThread_F11_P22::FockBuildThread_F11_P22(
    const Ref<MessageGrp> &msg,
    int nthread,
    int threadnum,
    const Ref<FockContribution>&c,
    const Ref<ThreadLock> &lock,
    const Ref<Integral> &integral,
    double acc, const signed char *pmax,
    const Ref<PetiteList> &pl,
    const Ref<GaussianBasisSet> &basis1,
    const Ref<GaussianBasisSet> &basis2,
    const Ref<GaussianBasisSet> &basis3,
    const Ref<GaussianBasisSet> &basis4
    ):
  FockBuildThread(msg,nthread,threadnum,c,lock,integral,acc,pmax),
  pl_(pl),
  basis1_(basis1),
  basis2_(basis2),
  basis3_(basis3),
  basis4_(basis4)
{
}

void
FockBuildThread_F11_P22::run()
{
  run_J();
  run_K();
}

void
FockBuildThread_F11_P22::run_J()
{
  ThreadLockHolder eri_alloc_lock(lock_);
  integral_->set_basis(basis1_,basis2_,basis3_,basis4_);
  Ref<TwoBodyInt> eri = integral_->electron_repulsion();
  eri_alloc_lock.unlock();
  const double *buf = eri->buffer();
  int nshell1 = basis1_->nshell();
  int nshell2 = basis3_->nshell();
  for (int I=0; I<nshell1; I++) {
      int nI = basis1_->shell(I).nfunction();
      for (int J=0; J<=I; J++) {
          int nJ = basis2_->shell(J).nfunction();
          for (int K=0; K<nshell2; K++) {
              int nK = basis3_->shell(K).nfunction();
              for (int L=0; L<=K; L++) {
                  int nL = basis4_->shell(L).nfunction();
                  eri->compute_shell(I,J,K,L);
                  if (I == J) {
                      if (K == L) {
                          contrib_->contrib_e_J(1.0, I, J, K, L,
                                                nI, nJ, nK, nL, buf);
                        }
                      else {
                          contrib_->contrib_p34_J(1.0, I, J, K, L,
                                                  nI, nJ, nK, nL, buf);
                        }
                    }
                  else if (K == L) {
                      contrib_->contrib_p12_J(1.0, I, J, K, L,
                                              nI, nJ, nK, nL, buf);
                    }
                  else {
                      contrib_->contrib_p12_p34_J(1.0, I, J, K, L,
                                                  nI, nJ, nK, nL, buf);
                    }
                }
            }
        }
    }
}

void
FockBuildThread_F11_P22::run_K()
{
  ThreadLockHolder eri_alloc_lock(lock_);
  integral_->set_basis(basis1_,basis3_,basis2_,basis4_);
  Ref<TwoBodyInt> eri = integral_->electron_repulsion();
  eri_alloc_lock.unlock();
  const double *buf = eri->buffer();
  int nshell1 = basis1_->nshell();
  int nshell2 = basis3_->nshell();
  for (int I=0; I<nshell1; I++) {
      int nI = basis1_->shell(I).nfunction();
      for (int J=0; J<nshell2; J++) {
          int nJ = basis3_->shell(J).nfunction();
          for (int K=0; K<=I; K++) {
              int nK = basis2_->shell(K).nfunction();
              int Lfence = (I==K)?J+1:nshell2;
              for (int L=0; L<Lfence; L++) {
                  int nL = basis4_->shell(L).nfunction();
                  eri->compute_shell(I,J,K,L);
                  if (I == K && J == L) {
                      contrib_->contrib_e_K(1.0, I, J, K, L,
                                            nI, nJ, nK, nL, buf);
                    }
                  else {
                      contrib_->contrib_p13p24_K(1.0, I, J, K, L,
                                                 nI, nJ, nK, nL, buf);
                    }
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////
// FockBuild

FockBuild::FockBuild(const Ref<FockContribution> &contrib,
                     double acc,
                     const Ref<GaussianBasisSet> &b_f1,
                     const Ref<GaussianBasisSet> &b_f2,
                     const Ref<GaussianBasisSet> &b_p1,
                     const Ref<GaussianBasisSet> &b_p2,
                     const Ref<MessageGrp> &msg,
                     const Ref<ThreadGrp> &thr,
                     const Ref<Integral> &integral):
  contrib_(contrib),
  accuracy_(acc),
  b_f1_(b_f1),
  b_f2_(b_f2),
  b_p1_(b_p1),
  b_p2_(b_p2),
  msg_(msg),
  thr_(thr),
  integral_(integral)
{
  if (b_f2_.null()) b_f2_ = b_f1_;
  if (b_p1_.null()) b_p1_ = b_f1_;
  if (b_p2_.null()) b_p2_ = b_p1_;
}

FockBuild::~FockBuild()
{
}

template <class T>
FockBuildThread *
create_FockBuildThread(const Ref<MessageGrp> &msg,
                       int nthread,
                       int threadnum,
                       const Ref<FockContribution>&c,
                       const Ref<ThreadLock> &lock,
                       const Ref<Integral> &integral,
                       double acc, const signed char *pmax,
                       const Ref<PetiteList> &pl,
                       const Ref<GaussianBasisSet> &basis1,
                       const Ref<GaussianBasisSet> &basis2 = 0,
                       const Ref<GaussianBasisSet> &basis3 = 0,
                       const Ref<GaussianBasisSet> &basis4 = 0)
{
  return new T(msg,nthread,threadnum,c,lock,integral,acc,pmax,
               pl,basis1,basis2,basis3,basis4);
}

void
FockBuild::build_generic(FBT_CTOR ctor)
{
  int nthread = thr_->nthread();
  std::vector<Ref<FockContribution> > contribs(nthread);
  Ref<ThreadLock> lock = thr_->new_lock();
  
  Ref<PetiteList> pl = integral_->petite_list(b_p1_);

  signed char *pmax = contrib_->compute_pmax(true);

  for (int i=0; i<nthread; i++) {
      if (i==0) contribs[i] = contrib_;
      else contribs[i] = contrib_->clone();
      FockBuildThread *thread = (*ctor)(
          msg_,
          nthread,
          i,
          contribs[i],
          lock,
          integral_,
          accuracy_,
          pmax,
          pl,
          b_f1_,
          b_f2_,
          b_p1_,
          b_p2_);
      contribs[i]->copy();
      thr_->add_thread(i, thread);
    }
  thr_->start_threads();
  thr_->wait_threads();
  thr_->delete_threads();
  delete[] pmax;
  for (int i=1; i<nthread; i++) {
      contrib_->accum(contribs[i]);
    }
  contrib_->update();
}

void
FockBuild::build()
{
  if (b_f1_->equiv(b_f2_)) {
      if (b_p1_->equiv(b_p2_)) {
          if (b_f1_->equiv(b_p1_)) {
              build_generic(create_FockBuildThread<FockBuildThread_F11_P11>);
              return;
            }
          else {
//               build_generic(create_FockBuildThread<FockBuildThread_F11_P22>);
//               return;
            }
        }
    }
  else if (b_p1_->equiv(b_p2_)) {
    }
  else {
      build_generic(create_FockBuildThread<FockBuildThread_F12_P34>);
      return;
    }
  // use the general code if this case not handled above
  build_generic(create_FockBuildThread<FockBuildThread_F12_P34>);
  return;
}

}


// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
