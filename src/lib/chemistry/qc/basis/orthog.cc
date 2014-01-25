//
// orthog.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: LPS
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

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <util/state/stateio.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/orthog.h>
#include <math/scmat/elemop.h>
#include <math/scmat/blocked.h>

using namespace std;
using namespace sc;

static ClassDesc OverlapOrthog_cd(typeid(OverlapOrthog),"OverlapOrthog",1,
                                  "virtual public SavableState",
                                  0, 0, create<OverlapOrthog>);

OverlapOrthog::OverlapOrthog(
  OrthogMethod method,
  const RefSymmSCMatrix &overlap,
  const Ref<SCMatrixKit> &result_kit,
  double lindep_tolerance,
  int debug) : nlindep_(0), min_orthog_res_(0.0), max_orthog_res_(0.0)
{
  reinit(method,overlap,result_kit,lindep_tolerance,debug);
}

OverlapOrthog::OverlapOrthog(StateIn& si):
  SavableState(si)
{
  Ref<SCMatrixKit> kit;
  kit << si.override()->describedclassvalue("matrixkit");

  if (kit == 0) {
    throw std::runtime_error("OverlapOrthog::OverlapOrthog(StateIn& si): requires that a matrixkit be set up in the override info");
    }

  si.get(debug_);
  dim_ << SavableState::restore_state(si);
  orthog_dim_ << SavableState::restore_state(si);
  si.get(lindep_tol_);
  int i_orthog_method;
  si.get(i_orthog_method);
  orthog_method_ = OrthogMethod(i_orthog_method);

  orthog_trans_ = kit->matrix(orthog_dim_, dim_);
  orthog_trans_.restore(si);
  orthog_trans_inverse_ = kit->matrix(dim_, orthog_dim_);
  orthog_trans_inverse_.restore(si);
  si.get(min_orthog_res_);
  si.get(max_orthog_res_);
  si.get(nlindep_);
}

OverlapOrthog::~OverlapOrthog()
{
}

void
OverlapOrthog::save_data_state(StateOut& so)
{
  so.put(debug_);
  SavableState::save_state(dim_.pointer(), so);
  SavableState::save_state(orthog_dim_.pointer(), so);
  so.put(lindep_tol_);
  so.put(int(orthog_method_));

  orthog_trans_.save(so);
  orthog_trans_inverse_.save(so);
  so.put(min_orthog_res_);
  so.put(max_orthog_res_);
  so.put(nlindep_);
  // The overlap_ member is not saved since should not be needed.
  // The result_kit_ member is not saved since it depends on the
  // runtime environment.  It is given to the StateIn CTOR by
  // an overriding KeyVal.
}

void
OverlapOrthog::reinit(
  OrthogMethod method,
  const RefSymmSCMatrix &overlap,
  const Ref<SCMatrixKit> &result_kit,
  double lindep_tolerance,
  int debug)
{
  orthog_method_ = method;
  overlap_ = overlap;
  lindep_tol_ = lindep_tolerance;
  debug_ = debug;
  dim_ = overlap_.dim();
  result_kit_ = result_kit;
}

Ref<OverlapOrthog>
OverlapOrthog::copy() const
{
  Ref<OverlapOrthog> orthog
    = new OverlapOrthog(orthog_method_,
                        overlap_,
                        result_kit_,
                        lindep_tol_,
                        debug_);

  orthog->orthog_trans_ = orthog_trans_.copy();
  orthog->orthog_trans_inverse_ = orthog_trans_inverse_.copy();
  orthog->orthog_dim_ = orthog_dim_;
  orthog->min_orthog_res_ = min_orthog_res_;
  orthog->max_orthog_res_ = max_orthog_res_;
  orthog->nlindep_ = nlindep_;
  return orthog;
}

// computes intermediates needed to form orthogonalization matrices
// and their inverses.
void
OverlapOrthog::compute_overlap_eig(RefSCMatrix& overlap_eigvec,
                    RefDiagSCMatrix& overlap_isqrt_eigval,
                    RefDiagSCMatrix& overlap_sqrt_eigval)
{
  // first calculate S
  RefSymmSCMatrix M = overlap_;

  // Diagonalize M to get m and U
  RefSCMatrix U(M.dim(), M.dim(), M.kit());
  RefDiagSCMatrix m(M.dim(), M.kit());
  M.diagonalize(m,U);
  M = 0;

  // if given the number of linear dependencies to eliminate
  // compute the needed tolerance
  if (lindep_tol_ < 0.0) {
    const int target_nlindep = static_cast<int>(std::floor(-lindep_tol_));
    const int n = m.n();
    std::vector<double> mvec(n);
    for(int i=0; i<n; ++i) mvec[i] = m(i);
    std::sort(mvec.begin(), mvec.end());
    // zero linear dependencies is easy
    if (target_nlindep == 0)
      lindep_tol_ =  mvec[0] / (2.0 * mvec[n-1] );
    else {
      // exactly degenerate metric eigenvalues at cutoff are trouble!
      if (mvec[target_nlindep-1] == mvec[target_nlindep])
        throw AlgorithmException("degenerate metric eigenvalues during orthogonalization -- cannot ensure that nlindep is obeyed",
                                 __FILE__, __LINE__, this->class_desc());
      lindep_tol_ = (mvec[target_nlindep-1] + mvec[target_nlindep]) / (2.0 * mvec[n-1]);
    }
  }

  Ref<SCElementMaxAbs> maxabsop = new SCElementMaxAbs;
  m.element_op(maxabsop.pointer());
  double maxabs = maxabsop->result();
  double s_tol = lindep_tol_ * maxabs;

  double minabs = maxabs;
  BlockedDiagSCMatrix *bm = dynamic_cast<BlockedDiagSCMatrix*>(m.pointer());
  bool blocked;
  int nblocks;
  if (bm == 0) {
    blocked = false;
    nblocks = 1;
    }
  else {
    blocked = true;
    nblocks = bm->nblocks();
    }
  int i, j;
  double *pm_sqrt = new double[m->dim()->n()];
  double *pm_isqrt = new double[m->dim()->n()];
  int *pm_index = new int[m->dim()->n()];
  int *nfunc = new int[nblocks];
  int nfunctot = 0;
  nlindep_ = 0;
  for (i=0; i<nblocks; i++) {
      nfunc[i] = 0;
      if (blocked && bm->block(i) == 0) continue;
      int n;
      if (blocked) n = bm->block(i)->dim()->n();
      else n = m->dim()->n();
      double *pm = new double[n];
      if (blocked) bm->block(i)->convert(pm);
      else m->convert(pm);
      for (j=0; j<n; j++) {
          if (pm[j] > s_tol) {
              if (pm[j] < minabs) { minabs = pm[j]; }
              pm_sqrt[nfunctot] = sqrt(pm[j]);
              pm_isqrt[nfunctot] = 1.0/pm_sqrt[nfunctot];
              pm_index[nfunctot] = j;
              nfunc[i]++;
              nfunctot++;
            }
          else if (orthog_method_ == Symmetric) {
              pm_sqrt[nfunctot] = 0.0;
              pm_isqrt[nfunctot] = 0.0;
              pm_index[nfunctot] = j;
              nfunc[i]++;
              nfunctot++;
              nlindep_++;
            }
          else {
              nlindep_++;
            }
        }
      delete[] pm;
    }

  if (nlindep_ > 0 && orthog_method_ == Symmetric) {
    ExEnv::out0() << indent
                 << "WARNING: " << nlindep_
                 << " basis function"
                 << (dim_.n()-orthog_dim_.n()>1?"s":"")
                 << " ignored in symmetric orthogonalization."
                 << endl;
  }

  // make sure all nodes end up with exactly the same data
  MessageGrp::get_default_messagegrp()->bcast(nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(nfunc, nblocks);
  MessageGrp::get_default_messagegrp()->bcast(pm_sqrt,nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(pm_isqrt,nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(pm_index,nfunctot);

  if (orthog_method_ == Symmetric) {
    orthog_dim_ = new SCDimension(m->dim()->blocks(),
                                  "ortho basis (symmetric)");
    }
  else {
      orthog_dim_ = new SCDimension(nfunctot, nblocks,
                                nfunc, "ortho basis (canonical)");
      for (i=0; i<nblocks; i++) {
        orthog_dim_->blocks()->set_subdim(i, new SCDimension(nfunc[i]));
      }
    }

  overlap_eigvec = result_kit_->matrix(dim_, orthog_dim_);
  if (orthog_method_ == Symmetric) {
      overlap_eigvec.assign(U);
    }
  else {
      if (blocked) {
	  BlockedSCMatrix *bev
	      = dynamic_cast<BlockedSCMatrix*>(overlap_eigvec.pointer());
	  BlockedSCMatrix *bU
	      = dynamic_cast<BlockedSCMatrix*>(U.pointer());
	  int ifunc = 0;
	  for (i=0; i<bev->nblocks(); i++) {
	      if (bev->block(i) == 0) continue;
	      for (j=0; j<nfunc[i]; j++) {
		  RefSCVector col = bU->block(i)->get_column(pm_index[ifunc]);
		  bev->block(i)->assign_column(col,j);
		  col = 0;
		  ifunc++;
	      }
	  }
      }
      else {
	  for (j=0; j<nfunc[0]; j++) {
	      RefSCVector col = U->get_column(pm_index[j]);
	      overlap_eigvec->assign_column(col,j);
	      col = 0;
	  }
      }
    }

  overlap_sqrt_eigval = result_kit_->diagmatrix(orthog_dim_);
  overlap_sqrt_eigval->assign(pm_sqrt);
  overlap_isqrt_eigval = result_kit_->diagmatrix(orthog_dim_);
  overlap_isqrt_eigval->assign(pm_isqrt);

  delete[] nfunc;
  delete[] pm_sqrt;
  delete[] pm_isqrt;
  delete[] pm_index;
  
  max_orthog_res_ = maxabs;
  min_orthog_res_ = minabs;

  if (debug_ > 1) {
    overlap_.print("S");
    overlap_eigvec.print("S eigvec");
    overlap_isqrt_eigval.print("s^(-1/2) eigval");
    overlap_sqrt_eigval.print("s^(1/2) eigval");
  }
}

void
OverlapOrthog::compute_symmetric_orthog()
{
  RefSCMatrix overlap_eigvec;
  RefDiagSCMatrix overlap_isqrt_eigval;
  RefDiagSCMatrix overlap_sqrt_eigval;
  compute_overlap_eig(overlap_eigvec,
                      overlap_isqrt_eigval,
                      overlap_sqrt_eigval);

  orthog_trans_ = overlap_eigvec
    * overlap_isqrt_eigval
    * overlap_eigvec.t();
  orthog_trans_inverse_ = overlap_eigvec
    * overlap_sqrt_eigval
    * overlap_eigvec.t();
}

void
OverlapOrthog::compute_canonical_orthog()
{
  RefSCMatrix overlap_eigvec;
  RefDiagSCMatrix overlap_isqrt_eigval;
  RefDiagSCMatrix overlap_sqrt_eigval;
  compute_overlap_eig(overlap_eigvec,
                      overlap_isqrt_eigval,
                      overlap_sqrt_eigval);

  orthog_trans_ = overlap_isqrt_eigval * overlap_eigvec.t();
  orthog_trans_inverse_ = overlap_eigvec * overlap_sqrt_eigval;
}

void
OverlapOrthog::compute_gs_orthog()
{
  // Orthogonalize each subblock of the overlap.
  max_orthog_res_ = 1.0;
  min_orthog_res_ = 1.0;
  nlindep_ = 0;
  BlockedSymmSCMatrix *S
    = dynamic_cast<BlockedSymmSCMatrix *>(overlap_.pointer());
  int nblock = S->nblocks();
  Ref<BlockedSCMatrixKit> kit
    = dynamic_cast<BlockedSCMatrixKit*>(S->kit().pointer());
  Ref<SCMatrixKit> subkit = kit->subkit();
  RefSCMatrix *blockorthogs = new RefSCMatrix[nblock];
  int *nblockorthogs = new int[nblock];
  int northog = 0;
  for (int i=0; i<nblock; i++) {
    RefSymmSCMatrix Sblock = S->block(i);
    if (Sblock == 0) {
      blockorthogs[i] = 0;
      nblockorthogs[i] = 0;
      continue;
      }
    RefSCDimension dim = Sblock->dim();
    RefSCMatrix blockorthog(dim,dim,subkit);
    blockorthog->unit();
    double res;
    int nblockorthog = blockorthog->schmidt_orthog_tol(Sblock, lindep_tol_,
                                                       &res);
    if (res < min_orthog_res_) min_orthog_res_ = res;
    blockorthogs[i] = blockorthog;
    nblockorthogs[i] = nblockorthog;
    northog += nblockorthog;
    nlindep_ += dim.n() - nblockorthog;
  }

  // Construct the orthog basis SCDimension object.
  Ref<SCBlockInfo> blockinfo
    = new SCBlockInfo(northog, nblock, nblockorthogs);
  for (int i=0; i<nblock; i++) {
    blockinfo->set_subdim(i, new SCDimension(nblockorthogs[i]));
  }
  orthog_dim_ = new SCDimension(blockinfo, "ortho (Gram-Schmidt)");

  // Replace each blockorthog by a matrix with only linear independent columns
  for (int i=0; i<nblock; i++) {
    if (nblockorthogs[i] == 0) continue;
    RefSCMatrix old_blockorthog = blockorthogs[i];
    blockorthogs[i] = subkit->matrix(dim_->blocks()->subdim(i),
                                     orthog_dim_->blocks()->subdim(i));
    blockorthogs[i].assign_subblock(old_blockorthog,
                                    0, dim_->blocks()->subdim(i).n()-1,
                                    0, orthog_dim_->blocks()->subdim(i).n()-1);
  }

  // Compute the inverse of each orthogonalization block.
  RefSCMatrix *inverse_blockorthogs = new RefSCMatrix[nblock];
  for (int i=0; i<nblock; i++) {
    if (nblockorthogs[i] == 0) {
      inverse_blockorthogs[i] = 0;
      }
    else {
      inverse_blockorthogs[i] = blockorthogs[i].gi();
      }
  }

  // Construct the complete transformation matrices
  orthog_trans_ = result_kit_->matrix(dim_, orthog_dim_);
  orthog_trans_inverse_ = result_kit_->matrix(orthog_dim_, dim_);
  orthog_trans_.assign(0.0);
  orthog_trans_inverse_.assign(0.0);
  BlockedSCMatrix *X
    = dynamic_cast<BlockedSCMatrix*>(orthog_trans_.pointer());
  BlockedSCMatrix *Xi
    = dynamic_cast<BlockedSCMatrix*>(orthog_trans_inverse_.pointer());
  for (int i=0; i<nblock; i++) {
    if (nblockorthogs[i] == 0) continue;
    int nrow = blockorthogs[i].rowdim().n();
    int ncol = blockorthogs[i].coldim().n();
    X->block(i).assign_subblock(blockorthogs[i],
                                0, nrow-1, 0, ncol-1,
                                0, 0);
    Xi->block(i).assign_subblock(inverse_blockorthogs[i],
                                 0, ncol-1, 0, nrow-1,
                                 0, 0);
  }
  orthog_trans_ = orthog_trans_.t();
  orthog_trans_inverse_ = orthog_trans_inverse_.t();

  delete[] blockorthogs;
  delete[] inverse_blockorthogs;
  delete[] nblockorthogs;
}

void
OverlapOrthog::compute_orthog_trans()
{
  switch(orthog_method_) {
  case GramSchmidt:
    ExEnv::out0() << indent
                 << "Using Gram-Schmidt orthogonalization."
                 << endl;
    compute_gs_orthog();
    break;
  case Symmetric:
    compute_symmetric_orthog();
    ExEnv::out0() << indent
                 << "Using symmetric orthogonalization."
                 << endl;
    break;
  case Canonical:
    compute_canonical_orthog();
    ExEnv::out0() << indent
                 << "Using canonical orthogonalization."
                 << endl;
    break;
  default:
    ExEnv::outn() << "OverlapOrthog::compute_orthog_trans(): bad orthog method"
                 << endl;
    abort();
  }

  ExEnv::out0() << indent
               << "n(basis):        ";
  for (int i=0; i<dim_->blocks()->nblock(); i++) {
    ExEnv::out0() << scprintf(" %5d", dim_->blocks()->size(i));
  }
  ExEnv::out0() << endl;

  if (dim_.n() != orthog_dim_.n()) {
    ExEnv::out0() << indent
                 << "n(orthog basis): ";
    for (int i=0; i<orthog_dim_->blocks()->nblock(); i++) {
      ExEnv::out0() << scprintf(" %5d", orthog_dim_->blocks()->size(i));
      }
    ExEnv::out0() << endl;

    ExEnv::out0() << indent
                 << "WARNING: " << dim_.n() - orthog_dim_.n()
                 << " basis function"
                 << (dim_.n()-orthog_dim_.n()>1?"s":"")
                 << " discarded."
                 << endl;
    }
  ExEnv::out0() << indent
               << "Maximum orthogonalization residual = "
               << max_orthog_res_ << endl
               << indent
               << "Minimum orthogonalization residual = "
               << min_orthog_res_ << endl;

  if (debug_ > 0) {
    dim_.print();
    orthog_dim_.print();
    if (debug_ > 1) {
      orthog_trans_.print("basis to orthog basis");
      orthog_trans_inverse_.print("basis to orthog basis inverse");
      (orthog_trans_*overlap_
       *orthog_trans_.t()).print("X*S*X'",ExEnv::out0(),14);
      (orthog_trans_inverse_.t()*overlap_.gi()
       *orthog_trans_inverse_).print("X'^(-1)*S^(-1)*X^(-1)",
                                     ExEnv::out0(),14);
      (orthog_trans_
       *orthog_trans_inverse_).print("X*X^(-1)",ExEnv::out0(),14);
    }
  }
}

// returns the orthogonalization matrix
RefSCMatrix
OverlapOrthog::basis_to_orthog_basis()
{
  if (orthog_trans_ == 0) {
    compute_orthog_trans();
  }
  return orthog_trans_;
}

RefSCMatrix
OverlapOrthog::basis_to_orthog_basis_inverse()
{
  if (orthog_trans_inverse_ == 0) {
    compute_orthog_trans();
  }
  return orthog_trans_inverse_;
}

RefSCDimension
OverlapOrthog::dim()
{
  return dim_;
}

RefSCDimension
OverlapOrthog::orthog_dim()
{
  if (orthog_dim_ == 0) compute_orthog_trans();
  return orthog_dim_;
}

int
OverlapOrthog::nlindep()
{
  if (orthog_dim_ == 0) compute_orthog_trans();
  return nlindep_; 
}

RefSymmSCMatrix
OverlapOrthog::overlap_inverse()
{
  RefSCMatrix X = basis_to_orthog_basis();
  RefSymmSCMatrix result(X.coldim(), X.kit());
  result.assign(0.0);
  result.accumulate_symmetric_product(X.t());
  return result;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
