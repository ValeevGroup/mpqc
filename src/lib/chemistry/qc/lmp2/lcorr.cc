
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

// References:
// [1] Schutz et al., JCP v. 111, p. 5691 (1999)
// [2] Hampel & Werner, JCP v. 104, p. 6286 (1996)

#include <math.h>

#include <algorithm>
#include <iterator>
#include <map>
#include <set>

#include <util/misc/regtime.h>
#include <util/group/message.h>
#include <util/group/mstate.h>
#include <util/group/pregtime.h>
#include <util/state/state_bin.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/molecule/formula.h>

#include <math/optimize/diis.h>
#include <math/scmat/repl.h>
#include <math/scmat/blocked.h>

#include <chemistry/qc/lmp2/sma.h>
#include <chemistry/qc/lmp2/parallel.h>
#include <chemistry/qc/lmp2/extrap.h>
#include <chemistry/qc/lmp2/domain.h>
#include <chemistry/qc/lmp2/pop_local.h>

#include <chemistry/qc/lmp2/lcorr.h>

namespace sc {

static ClassDesc LCorr_cd(typeid(LCorr), "LCorr", 1, "public Wavefunction",
                        0, 0, 0);

LCorr::LCorr(const Ref<KeyVal> &keyval):
  Wavefunction(keyval)
{
  KeyValValuedouble default_W_eigval_threshold(1e-6);
  W_eigval_threshold_ = keyval->doublevalue("W_eigval_threshold",
                                            default_W_eigval_threshold);
}

LCorr::LCorr(StateIn &statein):Wavefunction(statein)
{
  statein.get(W_eigval_threshold_);
}

LCorr::~LCorr()
{
}

void
LCorr::save_data_state(StateOut &stateout)
{
  stateout.put(W_eigval_threshold_);
}

void
LCorr::print_parameters() const
{
  ExEnv::out0() << indent << "W_eigval_threshold = " << W_eigval_threshold_
                << std::endl;
}

void
LCorr::init_virb_to_bfns(const sma2::Range &vir)
{
  virb_to_bfns_.resize(vir.nblock());
  for (int i=0; i<vir.nblock(); i++) {
      virb_to_bfns_[i].resize(vir.block_size(i));
      for (int j=0; j<vir.block_size(i); j++) {
          virb_to_bfns_[i][j] = j + vir.block_offset(i);
        }
    }
}

void
LCorr::clear()
{
  virb_to_bfns_.clear();
  unique_W_.clear();
  unique_eigvals_.clear();
  unique_F_tilde_.clear();
}

void
LCorr::compute_W(domainmapvirbs_t &virset,
                 const sc::Ref<sc::GaussianBasisSet> &basis,
                 sc::RefSCMatrix &F_vir_mat, sc::RefSCMatrix &S_mat,
                 const sma2::Range &vir, int nocc_act,
                 double bound)
{
  // See if unique_W has already been computed for this virset.
  if (unique_W_.find(virset) != unique_W_.end()) {
      return;
    }

  int domain_dim = 0;
  for (domainmapvirbs_t::iterator it = virset.begin();
       it != virset.end(); it++) {
      domain_dim += virb_to_bfns_[*it].size();
    }

  RefSCDimension domain_blockdim(new SCDimension(domain_dim));
  Ref<SCMatrixKit> kit = F_vir_mat->kit();

  // Assign elements of S_domain and F_vir_domain
  RefSymmSCMatrix S_domain(domain_blockdim,kit);
  RefSymmSCMatrix F_vir_domain(domain_blockdim,kit);
  int r_local = 0;
  for (domainmapvirbs_t::iterator vir_it1 = virset.begin();
       vir_it1 != virset.end(); vir_it1++) {
      for (int index1=0; index1<virb_to_bfns_[*vir_it1].size(); index1++) {
          int r_global = virb_to_bfns_[*vir_it1][index1];
          int s_local = 0;
          for (domainmapvirbs_t::iterator vir_it2 = virset.begin();
               vir_it2 != virset.end(); vir_it2++) {
              for (int index2=0; index2<virb_to_bfns_[*vir_it2].size(); index2++) {
                  int s_global = virb_to_bfns_[*vir_it2][index2];
                  S_domain(r_local,s_local) = S_mat->get_element(r_global,s_global);
                  F_vir_domain(r_local,s_local) = F_vir_mat->get_element(r_global,s_global);
                  ++s_local;
                }
            }
          ++r_local;
        }
    }
      
  // Diagonalize S_domain
  RefSCMatrix S_domain_eigvecs(domain_blockdim,domain_blockdim,kit);
  RefDiagSCMatrix S_domain_eigvals(domain_blockdim,kit);
  S_domain.diagonalize(S_domain_eigvals, S_domain_eigvecs);

  // From S_domain_eigvecs, create the matrix X_domain ( = X_tilde_domain from [2])
  int nzero_eigval = 0;
  for (int index=0; index<domain_dim; index++) {
      if (fabs(S_domain_eigvals.get_element(index)) <= W_eigval_threshold_) {
          nzero_eigval++;
        }
    }
  RefSCDimension domain_blockdim_nonred(new SCDimension(domain_dim - nzero_eigval));
  RefSCMatrix X_domain(domain_blockdim,domain_blockdim_nonred,kit);
  int col_index = 0;
  for (int index1=0; index1<domain_dim; index1++) {
      double eigenvalue = S_domain_eigvals.get_element(index1);
      if (fabs(eigenvalue) > W_eigval_threshold_) {
          for (int index2=0; index2<domain_dim; index2++) {
              X_domain(index2,col_index)
                  = S_domain_eigvecs(index2,index1)/sqrt(fabs(eigenvalue));
            }
          col_index++;
        }
    }
 
  // Compute the transformed matrix F_bar_domain = X_domain^T*F_vir_domain*X_domain
  RefSymmSCMatrix F_bar_domain(domain_blockdim_nonred,kit);
  F_bar_domain.assign(0.0);
  F_bar_domain.accumulate_transform(X_domain, F_vir_domain, SCMatrix::TransposeTransform);

  // Diagonalize F_bar_domain
  RefDiagSCMatrix Lambda(domain_blockdim_nonred,kit);
  RefSCMatrix U(domain_blockdim_nonred,domain_blockdim_nonred,kit);
  F_bar_domain.diagonalize(Lambda,U);

  // Construct W_domain matrix
  RefSCMatrix W_domain(domain_blockdim,domain_blockdim_nonred,kit);
  W_domain = X_domain*U;

  // Determine index mapping for converting SCMatrix etc. to Array objects.
  std::vector<int> index_map1, index_map2;
  compute_W_index_maps(domain_blockdim->n(), domain_blockdim_nonred->n(),
                       virset, index_map1, index_map2, vir);

  // Place F_tilde (==F_vir_domain) into the map holding all F_tilde's
  unique_F_tilde_[virset].init(vir,vir);
  pack_submatrix_into_empty_array(F_vir_domain, unique_F_tilde_[virset], bound,
                                  index_map1, index_map1);

  // likewise for W
  sma2::Range domain_vir_nonred(domain_blockdim_nonred,1);
  unique_W_[virset].init(vir,domain_vir_nonred);
  pack_submatrix_into_empty_array(W_domain, unique_W_[virset], bound,
                                  index_map1, index_map2);

  // and create a vector to hold the eigenvalues
  std::vector<double> &eigvec = unique_eigvals_[virset];
  eigvec.resize(Lambda.dim()->n());
  for (int i=0; i<eigvec.size(); i++) {
      eigvec[i] = Lambda(i);
    }
}


sma2::Array<2> &
LCorr::unique_W(const domainmapvirbs_t &i)
{
  std::map<domainmapvirbs_t, sma2::Array<2> >::iterator
      res = unique_W_.find(i);
  if (res == unique_W_.end()) {
      throw sc::ProgrammingError("missing unique_W",
                                 __FILE__, __LINE__,
                                 class_desc());
    }
  return res->second;
}

std::vector<double> &
LCorr::unique_eigvals(const domainmapvirbs_t &i)
{
  std::map<domainmapvirbs_t, std::vector<double> >::iterator
      res = unique_eigvals_.find(i);
  if (res == unique_eigvals_.end()) {
      throw sc::ProgrammingError("missing unique_eigvals",
                                 __FILE__, __LINE__,
                                 class_desc());
    }
  return res->second;
}

sma2::Array<2> &
LCorr::unique_F_tilde(const domainmapvirbs_t &i)
{
  std::map<domainmapvirbs_t, sma2::Array<2> >::iterator
      res = unique_F_tilde_.find(i);
  if (res == unique_F_tilde_.end()) {
      throw sc::ProgrammingError("missing unique_F_tilde",
                                 __FILE__, __LINE__,
                                 class_desc());
    }
  return res->second;
}

void
LCorr::compute_W_index_maps(int blockdim, int blockdim_nonred,
                            std::set<int> &virbs,
                            std::vector<int> &index_map1,
                            std::vector<int> &index_map2,
                            const sma2::Range &vir)
{
  // Determine index mapping
  index_map1.resize(blockdim);
  index_map2.resize(blockdim_nonred);
  for (int index=0; index<blockdim_nonred; index++)
      index_map2[index] = index;
  int index = 0;
  for (std::set<int>::iterator set_iter = virbs.begin();
       set_iter != virbs.end(); set_iter++) {
      int current_virb = *set_iter;
      int virb_offset = vir.block_offset(current_virb);
      // loop over virtual orbs in current virb and assign index_map1 values
      int nao = vir.block_size(current_virb);
      for (int iao=0; iao<nao; iao++) {
          index_map1[index] = iao + virb_offset;
          index++;
        }
    }
}

void
LCorr::transform_array(sma2::Array<2> &A, sma2::Array<2> &B,
                       sma2::Array<2> &C, sma2::Array<2> &D,
                       const sc::Ref<sc::MessageGrp> &msg)
{
  // Transform A to D using transformation matrices B,C: D = B^T*A*C

  // Parallelize this function by creating a local copy (C_local)
  // of C containing C(p,q) for all p and a subset of q and letting
  // each node process only its local copy of C

  int me = msg->me();
  int nproc = msg->n();

  // Create local, partial, copy of C
  sma2::Array<2> C_local(C.index(0), C.index(1), "C_local", C.tolerance());
  const sma2::Array<2>::blockmap_t &C_bmap = C.blockmap();
  sma2::BlockInfo<2> C_blockinfo;
  for (sma2::Array<2>::blockmap_t::const_iterator bmapiter = C_bmap.begin();
       bmapiter != C_bmap.end(); bmapiter++) {

      int index2 = (bmapiter->first).block(1);
      if (index2%nproc != me) continue; 
      int index1 = (bmapiter->first).block(0);

      C_blockinfo.block(0) = index1;
      C_blockinfo.block(1) = index2;
      C_local.add_unallocated_block(C_blockinfo);
    }
  C_local.zero();
  C_local("p","q") += C("p","q");

  sma2::Array<2> D_half(A.index(0), C_local.index(1),
                        "D_half", D.tolerance());
  sma2::Array<2> tmp_result(D.index(0), D.index(1),
                            "tmp_result", D.tolerance());
  
  D_half("p","q")     |= A("p","r") * C_local("r","q");
  D_half("p","q")      = A("p","r") *~C_local("r","q");

  tmp_result("p","q") |= B("r","p") * D_half("r","q");
  tmp_result("p","q")  = B("r","p") *~D_half("r","q");

  // tmp_result now contains a partial contribution to D;
  // assign this partial contribution to D on each node
  // (pruning out small blocks) and accumulate results

  D("p","q") |= tmp_result("p","q");
  D.parallel_union(msg); // need this to get the full D allocated on each node

  D.zero();
  D("p","q") += tmp_result("p","q");
  D.parallel_accumulate(msg);

}

}
