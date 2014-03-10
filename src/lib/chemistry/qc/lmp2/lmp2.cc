
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
#include <util/misc/scexception.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/molecule/formula.h>

#include <math/optimize/diis.h>
#include <math/scmat/repl.h>
#include <math/scmat/blocked.h>

#include "sma.h"
#include "parallel.h"
#include "extrap.h"
#include "domain.h"
#include "pop_local.h"
#include "analyze_array.h"

#include "lmp2.h"

// If 1, use W and eigvals for the full virtual space
// Should be 0 for normal operation.
#define COMPLETE_W 0

namespace sc {

static void
create_virb_united_pair_domains(int nvirb, domainmap_t &domainmap,
                                std::vector<std::set<int> >
                                &virb_united_pair_domains)
{
  virb_united_pair_domains.resize(nvirb);
  for (int ivirb=0; ivirb<nvirb; ivirb++) {
      std::set<int> &virb_united_pair_domain = virb_united_pair_domains[ivirb];
      virb_united_pair_domain.clear();
      // For the current virb, create the union of all pair domains
      // that contain this virb
      for (domainmap_t::iterator i_occpair = domainmap.begin();
           i_occpair != domainmap.end();
           i_occpair++) {
          // if current set contains ivirb, include virbs of
          // current set in virb_united_pair_domain
          const std::set<int> &current_set = i_occpair->second;
          if (current_set.find(ivirb) != current_set.end()) {
              virb_united_pair_domain.insert(current_set.begin(),
                                             current_set.end());
            }
        }
    }
}

void
LMP2::compute_doubles_W()
{
  // Compute unique_W, which is used to project the residual 
  // in each iteration; unique_eigvals is also assigned

  sc::Ref<sc::GaussianBasisSet> basis = this->basis();
  int nocc_act = occ_act_.nindex();

  Ref<SCMatrixKit> kit = basis->matrixkit();
  RefSCDimension ao_dim(new SCDimension(basis->nbasis()));
  RefSCDimension occ_act_dim(new SCDimension(nocc_act));

  // Create a matrixkit that is not parallelized so that compute_doubles_W
  // will not try to parallelize diagonalization of small matrices
  sc::Ref<sc::MessageGrp> saved_msg_grp = sc::MessageGrp::get_default_messagegrp();
  sc::MessageGrp::set_default_messagegrp(new ProcMessageGrp);
  sc::Ref<sc::SCMatrixKit> localkit = new ReplSCMatrixKit;
  sc::MessageGrp::set_default_messagegrp(saved_msg_grp);

  // Make matrix copies of F_vir and S
  RefSCMatrix F_vir_mat(ao_dim,ao_dim,localkit);
  RefSCMatrix S_mat(ao_dim,ao_dim,localkit); 
  pack_array_into_matrix(F_vir_,F_vir_mat);
  pack_array_into_matrix(S_,S_mat);
  
 for (my_occ_pairs_t::const_iterator iter = my_occ_pairs_.begin();
       iter != my_occ_pairs_.end(); iter++) {

     domainmapvirbs_t &virset = domainmap_[*iter];
     compute_W(virset, basis, F_vir_mat, S_mat, vir_, nocc_act, bound_);
   }

#if COMPLETE_W // This uses a W and eigvals for the full virtual space
      domainmapvirbs_t complete_virset;
      for (int i=0; i<vir_.nblock(); i++) {
          complete_virset.insert(i);
        }
      compute_W(complete_virset, basis, F_vir_mat, S_mat,
                vir_, nocc_act, bound_);
#endif

}

static ClassDesc LMP2_cd(typeid(LMP2), "LMP2", 1, "public LCorr",
                        0, create<LMP2>, create<LMP2>);

LMP2::LMP2(const Ref<KeyVal> &keyval):
  LCorr(keyval)
{
  KeyValValueboolean false_value(0);
  KeyValValueboolean true_value(1);

  msg_ = sc::MessageGrp::get_default_messagegrp();

  ref_ << keyval->describedclassvalue("reference");
  if(ref_.null()) {
      throw std::runtime_error("bad input: LMP2 requires reference");
    }

  copy_orthog_info(ref_);
  
  KeyValValueint defmaxiter(100);
  max_iter_ = keyval->intvalue("max_iter",defmaxiter);

  KeyValValuedouble defSthreshold(1.0e-6);
  S_threshold_ = keyval->doublevalue("S_threshold",defSthreshold);

  KeyValValuedouble defintegralthreshold(1.0e-8);
  integral_threshold_ = keyval->doublevalue("integral_threshold",defintegralthreshold);

  KeyValValuedouble defq1threshold(1.0e-8);
  q1_threshold_ = keyval->doublevalue("q1_threshold",defq1threshold);

  KeyValValuedouble defq2threshold(1.0e-8);
  q2_threshold_ = keyval->doublevalue("q2_threshold",defq2threshold);

  KeyValValuedouble defq3threshold(1.0e-8);
  q3_threshold_ = keyval->doublevalue("q3_threshold",defq3threshold);

  KeyValValuedouble defq4threshold(1.0e-8);
  q4_threshold_ = keyval->doublevalue("q4_threshold",defq4threshold);

  KeyValValuedouble def_thres_fact(10.0);
  threshold_factor_ = keyval->doublevalue("threshold_factor",def_thres_fact);
  
  KeyValValueboolean defcompletedomains(0);
  completedomains_ = keyval->booleanvalue("complete_domains",defcompletedomains);
  
  KeyValValueboolean defsinglevirb(0);
  singlevirb_ = keyval->booleanvalue("single_virtual_block",defsinglevirb);

  KeyValValueboolean def_i_ge_j(1);
  i_ge_j_ = keyval->booleanvalue("i_ge_j",def_i_ge_j);

  KeyValValueboolean def_m_ge_n(1);
  m_ge_n_ = keyval->booleanvalue("m_ge_n",def_m_ge_n);

  KeyValValueboolean def_minimize_q2(1);
  minimize_q2_ = keyval->booleanvalue("minimize_q2",def_minimize_q2);

  KeyValValueboolean def_always_use_dist_t(0);
  always_use_dist_t_
      = keyval->booleanvalue("always_use_dist_t",def_always_use_dist_t);

  KeyValValuedouble def_dist_thres(15.0);  // (use 15.0 bohr by default)
  distance_threshold_ = keyval->doublevalue("distance_threshold",def_dist_thres);

  KeyValValuedouble def_comp_thres(0.02);  // (use 0.02 by default)
  completeness_threshold_ = keyval->doublevalue("completeness_threshold",def_comp_thres);

  if (keyval->exists("extrap_T")) {
      extrap_T_ << keyval->describedclassvalue("extrap_T");
    }
  else {
      extrap_T_ = new DIIS(6,8,0.0,2,1);
    }

  // By default, no core orbitals are frozen; frozen core orbitals
  // can be specified either with integer value or "auto"
  nfzc_ = keyval->intvalue("nfzc");
  std::string nfzc_charval = keyval->stringvalue("nfzc");
  if (nfzc_charval == "auto") {
      if (molecule()->max_z() > 30) {
          ExEnv::err0()
               << "LMP2: cannot use \"nfzc = auto\" for Z > 30" << std::endl;
          abort();
        }
      nfzc_ = molecule()->n_core_electrons()/2;
      ExEnv::out0() << indent
           << "LMP2: auto-freezing " << nfzc_ << " core orbitals" << std::endl;
    }

  KeyValValuestring def_occ_orbitals("pipek-mezey");
  occ_orbitals_ = keyval->stringvalue("occ_orbitals", def_occ_orbitals);
  if (occ_orbitals_ != "pipek-mezey"
      && occ_orbitals_ != "canonical") {
      throw std::runtime_error("LMP2: invalid occ_orbitals input");
    }

  KeyValValuestring def_vir_orbitals("projected_atomic");
  vir_orbitals_ = keyval->stringvalue("vir_orbitals", def_vir_orbitals);
  if (vir_orbitals_ != "projected_atomic"
      && vir_orbitals_ != "old_projected_atomic"
      && vir_orbitals_ != "canonical") {
      throw std::runtime_error("LMP2: invalid vir_orbitals input");
    }

  analyze_occ_orbs_ = keyval->booleanvalue("analyze_occ_orbs", KeyValValueboolean(false));
  if (analyze_occ_orbs_) {
    const unsigned int nocc = nelectron() / 2;
    const unsigned int nocc_act = nocc - nfzc_;
    emp2_ij_.resize(nocc_act*nocc_act); std::fill(emp2_ij_.begin(), emp2_ij_.end(), 0.0);
    dist_ij_ = emp2_ij_;
    r_i_.resize(nocc_act);  std::fill(r_i_.begin(), r_i_.end(), SCVector3(0.0, 0.0, 0.0));
  }
}

LMP2::LMP2(StateIn &statein):LCorr(statein)
{
  throw sc::FeatureNotImplemented("LMP2 StateIn CTOR",
                                  __FILE__, __LINE__, class_desc());
  ref_ << SavableState::restore_state(statein);
  statein.get(S_threshold_);
  statein.get(integral_threshold_);
  statein.get(q1_threshold_);
  statein.get(q2_threshold_);
  statein.get(q3_threshold_);
  statein.get(q4_threshold_);
}

void
LMP2::save_data_state(StateOut &stateout)
{
  LCorr::save_data_state(stateout);

  SavableState::save_state(ref_.pointer(),stateout);
  
  stateout.put(S_threshold_);
  stateout.put(integral_threshold_);
  stateout.put(q1_threshold_);
  stateout.put(q2_threshold_);
  stateout.put(q3_threshold_);
  stateout.put(q4_threshold_);
}

void
LMP2::compute(void)
{
#if defined(HAVE_IEEE_SET_FP_CONTROL)
  ieee_set_fp_control(0);
#endif

  if(gradient_needed()) {
      throw std::logic_error("LMP2 cannot do gradients");
    }

  double extra_hf_acc = 10.;
  ref_->set_desired_value_accuracy(desired_value_accuracy()
                                   / extra_hf_acc);
  double refenergy = ref_->energy();
  double lmp2energy = compute_lmp2_energy();

  double totalenergy = refenergy + lmp2energy;

  set_value(totalenergy);

  ExEnv::out0() << std::endl;
  ExEnv::out0() << indent
                << scprintf("E[Reference]    = % 15.10f", refenergy)
                << std::endl;
  std::string basemethod("LMP2");
  ExEnv::out0() << indent
                << scprintf("E_corr[%s]%*s = % 15.10f",
                            basemethod.c_str(), 7-basemethod.size(), "",
                            lmp2energy)
                << std::endl;
  ExEnv::out0() << indent
                << scprintf("E[%s]%*s = % 15.10f",
                            basemethod.c_str(), 12-basemethod.size(), "",
                            refenergy+lmp2energy)
                << std::endl;

  clear();
}

int
LMP2::nelectron(void)
{
  return ref_->nelectron();
}

RefSymmSCMatrix
LMP2::density(void)
{
  throw std::runtime_error("LMP2 cannot compute density");
  return 0;
}

double
LMP2::magnetic_moment() const
{
  return 0;
}

int
LMP2::value_implemented(void) const
{
  return 1;
}

double
LMP2::compute_ecorr_lmp2()
{
  Timer tim("ecorr");

  const unsigned nocc = this->nelectron() / 2;
  const unsigned nocc_act = nocc - nfzc_;

  sma2::Index r("r"), s("s");
  sma2::Array<0> ecorr;
  double ecorr_lmp2 = 0.0;
  for (my_occ_pairs_t::const_iterator iter = my_occ_pairs_.begin();
       iter != my_occ_pairs_.end();
       iter++) {
      sma2::Index i(iter->first-nfzc_);
      sma2::Index j(iter->second-nfzc_);
      if (j.value() > i.value()) continue;
      double f;
      if (i.value() != j.value()) f = 2.0;
      else                        f = 1.0;
      ecorr.zero();
      ecorr() +=  f * 2.0 * K_2occ_(i,j,r,s) * T_local_(i,j,r,s);
      ecorr() -=  f *       K_2occ_(i,j,s,r) * T_local_(i,j,r,s);
      const double eij = ecorr.value();
      ecorr_lmp2 += eij;

      if (analyze_occ_orbs_)
        emp2_ij_[i.value() * nocc_act + j.value()] = eij;
    }

  msg_->sum(ecorr_lmp2);
  if (analyze_occ_orbs_)
    msg_->sum(&emp2_ij_[0], emp2_ij_.size());

  return ecorr_lmp2;
}

LMP2::~LMP2()
{
}

void
LMP2::clear()
{
  LCorr::clear();

  P_.clear();
  L_.clear();
  domainmap_.clear();
  K_2occ_.clear();
  T_jirs_.clear();
  T_.clear();
  F_occ_.clear();
  F_vir_.clear();
  S_.clear();
  F_diag_.clear();
}

void
LMP2::rearrange_q2_all_ij(sma2::Array<4> &q2_K_oo)
{
  // At this point q2_K_oo has only those blocks with M>=N.  When it is
  // redistributed to q2_K_oo_redist, then we want all M, N.  We will
  // modify q2_K_oo so that it has all M, N, i, j.  It will not
  // have the correct M, N distribution, but will be set up so that will be
  // taken care of in the redistribution to i, j below.
  // IMPORTANT: q2_K_oo's blockmap is copied, because we will change it below.
  sma2::Array<4>::blockmap_t q2_K_oo_bm(q2_K_oo.blockmap());
  for (sma2::Array<4>::blockmap_t::iterator iter = q2_K_oo_bm.begin();
       iter != q2_K_oo_bm.end(); iter++) {
      sma2::BlockInfo<4> &bi = iter->first;
      int M = bi.block(0);
      int N = bi.block(1);

      // We have q2(MNij).
      // Add q2(NMji).
      if (M != N) {
          int i = bi.block(2);
          int j = bi.block(3);
          const double *src = iter->second;
          sma2::BlockInfo<4> bi_perm;
          bi_perm.block(0) = N;
          bi_perm.block(1) = M;
          bi_perm.block(2) = j;
          bi_perm.block(3) = i;
          sma2::Array<4>::blockmap_t::value_type &block
              = q2_K_oo.add_allocated_block(bi_perm);
          // Instead of block iters, we'll do this more efficiently by
          // noting that the blocksize for indices 2 & 3 are always 1.
          double *dst = block.second;
          int nN = q2_K_oo.index(0).block_size(N);
          int nM = q2_K_oo.index(0).block_size(M);
          for (int i=0, ij=0; i<nM; i++) {
              for (int j=0, ji=i; j<nN; j++, ij++, ji+=nM) {
                  dst[ji] = src[ij];
                }
            }
        }
    }
  q2_K_oo_bm.clear();
}

void
LMP2::rearrange_q2_i_ge_j(sma2::Array<4> &q2_K_oo)
{
  // At this point q2_K_oo has only those blocks with M>=N.  When it is
  // redistributed to q2_K_oo_redist, then we want all M, N.  We will
  // modify q2_K_oo so that it has all M, N, i>=j.  It will not
  // have the correct M, N distribution, but will be set up so that will be
  // taken care of in the redistribution to i, j below.

  int max_block_size = q2_K_oo.index(0).max_block_size();
  sc::auto_vec<double> tmpsrc(new double[max_block_size*max_block_size]);

  // IMPORTANT: q2_K_oo's blockmap is copied, because we will change it below.
  sma2::Array<4>::blockmap_t q2_K_oo_bm(q2_K_oo.blockmap());
  for (sma2::Array<4>::blockmap_t::iterator iter = q2_K_oo_bm.begin();
       iter != q2_K_oo_bm.end(); iter++) {
      sma2::BlockInfo<4> &bi = iter->first;
      int M = bi.block(0);
      int N = bi.block(1);
      int j = bi.block(2);
      int i = bi.block(3);

      double *src = iter->second;
      sma2::BlockInfo<4> bi_perm;
      bi_perm.block(0) = N;
      bi_perm.block(1) = M;
      bi_perm.block(2) = i;
      bi_perm.block(3) = j;
      int nN = q2_K_oo.index(0).block_size(N);
      int nM = q2_K_oo.index(0).block_size(M);

      // We have q2(MNij).
      // If M == N && i==j, keep.
      // If M == N && i>j , keep.
      // If M == N && i<j , delete.
      // If M > N  && i==j, keep and remap a copy to q2(NMji)
      // If M > N  && i>j , keep.
      // If M > N  && i<j , Remap to q2(NMji).

      if (M == N) {
          if (i < j) {
              // Remove the block.  Depending on the
              // sma2::Array::blockmap_t implementation, this might not
              // actually deallocate the metadata for the block.
              q2_K_oo.remove_block(bi);
            }
          // Otherwise keep the block as it is.
        }
      else if (i == j) {
          // add q2(NMji)
          sma2::Array<4>::blockmap_t::value_type &block
              = q2_K_oo.add_allocated_block(bi_perm);
          double *dst = block.second;
          for (int i=0, ij=0; i<nM; i++) {
              for (int j=0, ji=i; j<nN; j++, ij++, ji+=nM) {
                  dst[ji] = src[ij];
                }
            }
        }
      else if (i > j) {
          // keep it unchanged
        }
      else {
          // replace with q2(NMji)
          double *dst = src;
          // src and dst are the same so copy src to a temporary
          double *tmpsrc_ptr = tmpsrc.get();
          memcpy(tmpsrc.get(),src,nM*nN*sizeof(double));
          q2_K_oo.relocate_block(bi,bi_perm);
          for (int i=0, ij=0; i<nM; i++) {
              for (int j=0, ji=i; j<nN; j++, ij++, ji+=nM) {
                  dst[ji] = tmpsrc_ptr[ij];
                }
            }
        }
    }
  q2_K_oo_bm.clear();
}

void
LMP2::compute_K_2occ(sc::RefSCMatrix &T_schwarz, sc::RefSCVector &T_schwarz_maxvec,
                     double T_schwarz_maxval, sc::RefSCMatrix &Pmax,
                     sc::RefSCVector &Pmaxvec,
                     std::vector<std::vector<double> > &Lmax,
                     sc::RefSCMatrix &Lshellmax,
                     std::vector<std::vector<double> > &PPmax,
                     const std::vector<std::vector<double> > &Dmax,
                     double &Dmax_maxval,
                     std::vector<double> &Lmaxvec,
                     std::vector<std::multimap<double,int,std::greater<double> > > &L_map)
{
  Timer tim1("Compute K_2occ");
  
  int me = msg_->me();
  int nproc = msg_->n();

  int nshell = basis()->nshell();
  int nocc_act = occ_act_.nindex();

  K_2occ_.init(occ_act_,occ_act_,vir_,vir_,"K_2occ",bound_);

  // Sort the T_schwarz array so we only need to loop through the
  // set of indices until we encounter the first value too small
  // to contribute.
  std::vector<std::multimap<double,int,std::greater<double> > >
      T_schwarz_sorted(nshell);
  for (int i=0; i<nshell; i++) {
      for (int j=0; j<nshell; j++) {
          double T_schwarz_ij = T_schwarz(i,j);
          T_schwarz_sorted[i].insert(std::make_pair(T_schwarz_ij,j));
        }
    }

  double Lmax_maxval = 0.0;
  for (int i=0; i<nshell; i++) {
      if (Lmax_maxval < Lmaxvec[i]) Lmax_maxval = Lmaxvec[i];
    }

  // Generate Dmax_maxvec[r] = max_s Dmax[r][s]
  std::vector<double> Dmax_maxvec(nshell);
  for (int i=0; i<nshell; i++) {
      double dmax_i_max = 0.0;
      const std::vector<double> &Dmax_i = Dmax[i];
      for (int j=0; j<nshell; j++) {
          double val = Dmax_i[j];
          if (val > dmax_i_max) dmax_i_max = val;
        }
      Dmax_maxvec[i] = dmax_i_max;
    }

  // Place Lshellmax in a structure with faster access.
  std::vector<std::vector<double> > Lshellmax_vv(nshell);
  for (int i=0; i<nshell; i++) {
      Lshellmax_vv[i].resize(nocc_act);
      for (int j=0; j<nocc_act; j++) {
          Lshellmax_vv[i][j] = Lshellmax(i,j);
        }
    }

  // Construct an integral evaluator
  sc::Ref<sc::TwoBodyInt> aoints = ref_->integral()->electron_repulsion();

  sma2::Array<4> q2_oo(ao_,ao_,occ_act_,occ_act_,"q2_oo",bound_);
  sma2::Array<4> q2_K_oo(ao_,ao_,occ_act_,occ_act_,"q2_K_oo",bound_);
  sma2::BlockInfo<4> q2_oo_blockinfo;
  sma2::BlockInfo<4> q2_K_oo_blockinfo;

  ///////////////////////////////////////////////////////////////
  // Do the first half of the integral transformation for K_2occ
  // (This is done directly, storing only parts of v and q1)
  //////////////////////////////////////////////////////////////
  Timer tim2("Half transf. for K_2occ");
  int counter = 0;
  for (int M=0; M<nshell; M++) {

      std::vector<double> &PPmax_M = PPmax[M];
      double Pmaxvec_M = Pmaxvec(M);
      double T_schwarz_maxvec_M = T_schwarz_maxvec(M);

      int Nfence = m_ge_n_?M+1:nshell;
      for (int N=0; N<Nfence; N++) {

          if (counter++%nproc != me) continue;

          double T_schwarz_maxvec_N = T_schwarz_maxvec(N);
          double PPmax_MN = PPmax_M[N];
          if (T_schwarz_maxvec_M*T_schwarz_maxvec_N*Dmax_maxval*PPmax_MN < integral_threshold_) continue;

          sma2::Array<4> v(ao_,ao_,ao_,ao_,"v");
          sma2::BlockInfo<4> v_blockinfo;
          v_blockinfo.block(0) = M;
          v_blockinfo.block(2) = N;

          q2_K_oo_blockinfo.block(0) = M;
          q2_K_oo_blockinfo.block(1) = N;

          Timer tim3("Allocate v");
          std::multimap<double, int, std::greater<double> >::iterator
              R_iter,
              R_iter_begin = T_schwarz_sorted[M].begin(),
              R_iter_end = T_schwarz_sorted[M].end();
          for (R_iter=R_iter_begin; R_iter!=R_iter_end; R_iter++) {
              double T_schwarz_MR = R_iter->first;
              double partial_bound_R = T_schwarz_MR*T_schwarz_maxvec_N*PPmax_MN;
              // If the the following is true, then we know we will not
              // find any larger contributions for this R, so break.
              int R = R_iter->second;
              double Dmax_R_max = Dmax_maxvec[R];
              if (partial_bound_R*Dmax_maxval < integral_threshold_) break;
              if (partial_bound_R*Dmax_R_max  < integral_threshold_) continue;
              v_blockinfo.block(1) = R;
              
              // Consider the integral array (ri|sj) and include all
              // integrals required for this array
              const std::vector<double> &Dmax_R = Dmax[R];
              std::multimap<double, int, std::greater<double> >::iterator
                  S_iter,
                  S_iter_begin = T_schwarz_sorted[N].begin(),
                  S_iter_end = T_schwarz_sorted[N].end();
              for (S_iter=S_iter_begin; S_iter!=S_iter_end; S_iter++) {
                  double T_schwarz_NS = S_iter->first;
                  double partial_bound_S = T_schwarz_MR*PPmax_MN*T_schwarz_NS;
                  // If the the following is true, then we know we will not
                  // find any larger contributions for this S, so break.
                  if (partial_bound_S*Dmax_R_max < integral_threshold_) break;
                  int S = S_iter->second;
                  if (partial_bound_S*Dmax_R[S] >= integral_threshold_) {
                      v_blockinfo.block(3) = S;
                      v.add_new_unallocated_block(v_blockinfo);
                    }
                }
            }

          // Compute the AO integrals for the current M,N shell pair
          tim3.change("Compute v");
          sma2::pack_2e_integrals_into_shell_blocked_array(v, aoints);

          if (v.max_abs_element()*Dmax_maxval*PPmax_MN < integral_threshold_) continue;
          
          // Allocate q1 for the current M,N pair
          tim3.change("Allocate q1");
          sma2::Array<4> q1(ao_,ao_,ao_,occ_act_,"q1",bound_);
          sma2::BlockInfo<4> q1_blockinfo;
          q1_blockinfo.block(0) = M;
          q1_blockinfo.block(1) = N;
          const sma2::Array<4>::blockmap_t &bmap = v.blockmap();
          for (sma2::Array<4>::blockmap_t::const_iterator bmapiter = bmap.begin();
               bmapiter != bmap.end(); bmapiter++) {
              const sma2::BlockInfo<4> &binfo = bmapiter->first;
              int R = binfo.block(1);
              int S = binfo.block(3);
              
              q1_blockinfo.block(2) = R;
              
#ifdef USE_BOUND
              double vmax      = binfo.bound();
#else
              double vmax      = v.block_max_abs(bmapiter);
#endif
              double lmaxval   = Lmaxvec[R];

              if (vmax*lmaxval*Lmaxvec[S]*PPmax_MN < q1_threshold_) continue; // imbn test screening

              const std::vector<double> &Lmax_R = Lmax[R];
              
              std::multimap<double,int,std::greater<double> >::iterator L_iter;
              for (L_iter = L_map[S].begin(); L_iter != L_map[S].end(); L_iter++) {
                  double Lshellmax_Si = L_iter->first; // this is Lshellmax(S,i)
                  double partial_bound = vmax*Lshellmax_Si*PPmax_MN;
                  if (partial_bound*lmaxval < q1_threshold_) break;
                  int i = L_iter->second;
                  if (partial_bound*Lmax_R[i] >= q1_threshold_) {
                      q1_blockinfo.block(3) = i;
                      // add block, if not already allocated (for another S value)
                      q1.add_unallocated_block(q1_blockinfo);
                    }
                }
            }

          // Compute q1 for the current M,N pair
          tim3.change("1st q.t.");
          q1("mu","nu","rho","j") = ~v("mu","rho","nu","sig") * L_("sig","j");

          // Allocate q2 array
          tim3.change("Allocate q2_K_oo");
          const sma2::Array<4>::blockmap_t &q1_bmap = q1.blockmap();
          for (sma2::Array<4>::blockmap_t::const_iterator bmapiter = q1_bmap.begin();
               bmapiter != q1_bmap.end(); bmapiter++) {  // Loop over nonzero blocks of q1
              const sma2::BlockInfo<4> &binfo = bmapiter->first;
              int R = binfo.block(2);
              int j = binfo.block(3);

              double *data = bmapiter->second;

              q2_K_oo_blockinfo.block(2) = j;
              
              double q1max = q1.block_max_abs(bmapiter);  // max abs. q1 value in current block
 
              // Loop over occupied orbitals l in the paired_occ[j] set,
              std::set<int>::const_iterator iter;
              const std::set<int> &current_set = paired_occ_[j];
              const std::vector<double> &Lshellmax_R = Lshellmax_vv[R];
              for (iter = current_set.begin(); iter != current_set.end(); iter++) {
                  int l = *iter;
                  double Lshellmax_Rl = Lshellmax_R[l];
                  if (q1max*Lshellmax_Rl*PPmax_MN >= q2_threshold_) {
                      q2_K_oo_blockinfo.block(3) = l;
                      q2_K_oo.add_zeroed_block(q2_K_oo_blockinfo);
                    }
                }
            }

          // Compute q2_K_oo contribution from current M,N pair
          tim3.change("2nd q.t.");
          sma2::Index imu("mu",M), inu("nu",N), irho("rho"), ij("j"), il("l");
          q2_K_oo(imu,inu,ij,il).skip_bounds_update() += ~q1(imu,inu,irho,ij) * L_(irho,il);
          
        } // exit N loop
    } // exit M loop

  analyze_array(q2_K_oo,"q2_K_oo(1)",msg_,true);

  if (m_ge_n_) {
      // At this point q2_K_oo has only those blocks with M>=N.  Make it
      // hold all M, N, and either all i, j or i>=j, depending on i_ge_j_.
      if (i_ge_j_) {
          if (minimize_q2_) {
              // This is the primary code path.  The others are for testing
              // purposes only.
              rearrange_q2_i_ge_j(q2_K_oo);
            }
          else {
              rearrange_q2_all_ij(q2_K_oo);
            }
        }
      else {
          rearrange_q2_all_ij(q2_K_oo);
        }
    }

  analyze_array(q2_K_oo,"q2_K_oo(2)",msg_,true);
  
  // Redistribute q2_oo and q2_K_oo (distrib. the two occupied indices) 
  tim2.change("Redistribute q2");
  sma2::Array<4> q2_K_oo_redist("q2_K_oo_redist");
  sma2::Array<4> *final_q2_K_oo;

  if (msg_->n() > 1) {
      sma2::PairBlockDistrib<4> d_K_oo(3,2,pair_mapping_);
      q2_K_oo_redist.distributed_from_distributed(msg_,d_K_oo,q2_K_oo,true,true);
      final_q2_K_oo = &q2_K_oo_redist;
    }
  else {
      final_q2_K_oo = &q2_K_oo;
    }
      
  ////////////////////////////////////////////////
  // Allocate three quarter transformed integrals
  // and do third quarter transformation
  ////////////////////////////////////////////////

  tim2.change("Allocate q3_ovo");
  sma2::Array<4> q3_ovo(ao_,occ_act_,vir_,occ_act_,"q3_ovo",bound_);
  sma2::BlockInfo<4> q3_ovo_blockinfo;
  const sma2::Array<4>::blockmap_t &q2_K_oo_bmap = final_q2_K_oo->blockmap();
  for (sma2::Array<4>::blockmap_t::const_iterator bmapiter = q2_K_oo_bmap.begin();
       bmapiter != q2_K_oo_bmap.end(); bmapiter++) {
      
      const sma2::BlockInfo<4> &binfo = bmapiter->first;
      int M = binfo.block(0);
      int i = binfo.block(3);
      int N = binfo.block(1);
      int j = binfo.block(2);

      double *data = bmapiter->second;
      double Pmaxvec_M = Pmaxvec(M);

      q3_ovo_blockinfo.block(0) = M;
      q3_ovo_blockinfo.block(1) = i;
      q3_ovo_blockinfo.block(3) = j;

      if (i_ge_j_ && j > i) continue;

      double q2max = final_q2_K_oo->block_max_abs(bmapiter);  // max abs. q2_K_oo value in MNji block

      // Loop over virbs s in the [ij] domain
      std::pair<int,int> current_pair = std::make_pair(i+nfzc_,j+nfzc_);
      std::set<int> k_2occ_domain = domainmap_[current_pair];
      std::set<int>::iterator tmp_iter;
      for (tmp_iter = k_2occ_domain.begin(); tmp_iter != k_2occ_domain.end();
           tmp_iter++) {
          int s = *tmp_iter;
          if (q2max*Pmax(N,s)*Pmaxvec_M  >= q3_threshold_) {
              q3_ovo_blockinfo.block(2) = s;
              q3_ovo.add_unallocated_block(q3_ovo_blockinfo);
            }
        }
    }

  analyze_array(q3_ovo,"q3_ovo",msg_,true);

  // Third quarter transformation
  tim2.change("3rd q.t.");
  q3_ovo("mu","i","s","j") = ~(*final_q2_K_oo)("mu","nu","j","i") * P_("nu","s");

  ////////////////////////////////////////
  // Allocate fully transformed integrals
  // and do fourth quarter transformation
  ////////////////////////////////////////

  tim2.change("Allocate K_2occ");
  sma2::BlockInfo<4> K_2occ_blockinfo;
  const sma2::Array<4>::blockmap_t &q3_ovo_bmap = q3_ovo.blockmap();
  // Allocate K_2occ
  for (sma2::Array<4>::blockmap_t::const_iterator bmapiter = q3_ovo_bmap.begin();
       bmapiter != q3_ovo_bmap.end(); bmapiter++) {
      
      const sma2::BlockInfo<4> &binfo = bmapiter->first;
      int M = binfo.block(0);
      int i = binfo.block(1);
      int s = binfo.block(2);
      int j = binfo.block(3);

      double *data = bmapiter->second;

      K_2occ_blockinfo.block(0) = i;
      K_2occ_blockinfo.block(1) = j;
      K_2occ_blockinfo.block(3) = s;

      double q3max = q3_ovo.block_max_abs(bmapiter);  // max abs. q3_ovo value in Misj block
          
      // Loop over virbs r in the [ij] domain
      std::pair<int,int> current_pair = std::make_pair(i+nfzc_,j+nfzc_);
      std::set<int> k_2occ_domain = domainmap_[current_pair];
      std::set<int>::iterator tmp_iter;
      for (tmp_iter = k_2occ_domain.begin(); tmp_iter != k_2occ_domain.end();
           tmp_iter++) {
          int r = *tmp_iter;
          if (q3max*Pmax(M,r) >= q4_threshold_) {
              K_2occ_blockinfo.block(2) = r;
              K_2occ_.add_unallocated_block(K_2occ_blockinfo);
            }
        }
    }

  analyze_array(K_2occ_,"K_2occ",msg_,true);

  // Fourth quarter transformation
  tim2.change("4th q.t.");
  K_2occ_("i","j","r","s") = ~q3_ovo("mu","i","s","j")*P_("mu","r");
}

void
LMP2::compute_Schwarz_screening_quantities(sc::Ref<TwoBodyInt> &tbint,
                                            sc::RefSCMatrix &T_schwarz,
                                            sc::RefSCVector &T_schwarz_maxvec,
                                            double &T_schwarz_maxval,
                                            int nshell)
{
  
  // Compute quantities used for Schwarz screening:
  // T_schwarz(M,N)      = max |(mu nu|mu nu)|^(1/2)
  //      (mu in M, nu in N; M, N are shells);
  // T_schwarz_maxvec(M) = max{T_schwarx(M,N)} over all N;
  // T_schwarz_maxval    = max{T_schwarz(M,N)} over all M,N

  Timer tim("compute_Schwarz_screening_quantities");
  
  T_schwarz_maxval = 0.0;
  for (int M=0; M<nshell; M++) {
    double maxval = 0.0;
    for (int N=0; N<nshell; N++) {
        double tmp = tbint->log2_shell_bound(M,N,M,N);
        // T_schwarz(M,N) = sqrt(2^tmp), where tmp = log2(max integral)
        double val = pow(2,tmp/2); 
        T_schwarz(M,N) = val;
        if (val > T_schwarz_maxval) T_schwarz_maxval = val;
        if (val > maxval) maxval = val;
      }
    T_schwarz_maxvec(M) = maxval;
  }

}

void
LMP2::compute_P_screening_quantities(sc::RefSCMatrix &Pmax,
                                     sc::RefSCMatrix &Pmatrix,
                                     sc::RefSCVector &Pmaxvec,
                                     double &Pmax_maxval)
{

  // Compute Pmax and Pmaxvec
  // Pmax(M,uvirb) = max|P(mu,u)|, mu in M, u running over virtual orbitals in uvirb
  // Pmaxvec(M) = max|P(mu,u)|, mu in M, all u

  Timer tim("compute_P_screening_quantities");
  
  sc::Ref<sc::GaussianBasisSet> basis = this->basis();

  Pmax_maxval = 0.0;
  int nshell = basis->nshell();
  for (int M=0; M<nshell; M++) {
      double maxval = 0.0;
      int size_M = basis->shell(M).nfunction();
      int M_offset = basis->shell_to_function(M);
      for (int uvirb=0; uvirb<vir_.nblock(); uvirb++) {
          int u_offset = vir_.block_offset(uvirb);
          double maxval2 = 0.0;
          for (int mu=0; mu<size_M; mu++) {
              for (int u=0; u<vir_.block_size(uvirb); u++) {
                  // Find max of |Pmatrix(mu+M_offset,u+u_offset)| where mu is an AO
                  // in shell M and put this into Pmax(mu+M_offset,uvirb) if greater
                  // than current value of Pmax(mu+M_offset,uvirb)
                  double val = fabs(Pmatrix(mu+M_offset,u+u_offset));
                  if (val > maxval2) maxval2 = val;
                  if (val > maxval)  maxval  = val;
                }
            }
          Pmax(M,uvirb) = maxval2;
        }
      Pmaxvec(M) = maxval;
      if (maxval > Pmax_maxval) Pmax_maxval = maxval;
    }

}

void
LMP2::compute_L_and_D_screening_quantities( std::vector<std::vector<double> > &Dmax,      // out
                                            double &Dmax_maxval,        // out
                                            std::vector<std::vector<double> > &Lmax, // out
                                            std::vector<double> &Lmaxvec,   // out
                                            sc::RefSCMatrix &Lshellmax, // out
                                            std::vector<std::vector<int> > &L_blocks, // out
                                            std::vector<std::multimap<double,int,std::greater<double> > > &L_map,  // out
                                            double &Lmax_maxval, // out
                                            double T_schwarz_maxval, // in
                                            double Pmax_maxval  // in
                                            )
{
  // Pre-compute various screening matrices and vectors
  // Lmax(M,j) = max of |L(m,i)| for all pairs ij (with j fixed) and m in shell M
  // Lmaxvec[M] = max of Lmax(M,j)| for all j (equal to max of Lshellmax(M,j) over all j)
  // Lshellmax(M,j) = max of |L(m,j)| for m in shell M
  // Dmax(M,N) = max over all valid pairs ij (j<=i) of |L_mu_i*L_nu_j|
  //             for mu in shell M and nu in shell N
  // Dmax_maxval = double holding the max of Dmax(M,N) over all M,N
  // L_blocks: vector of vectors; one vector for each shell, S, containing
  //          the i's for which the block L(S,i) exists
  // L_map: a vector of multimaps; there are nshell elements in the vector, and for
  //        each element (shell S) the value Lshellmax(S,i) and the corresponding i are kept;
  //        elements are sorted in descending order
  
  Timer routine_timer("compute_L_and_D_screening_quantities");

  sc::Ref<sc::GaussianBasisSet> basis = this->basis();

  int nocc_act = occ_act_.nindex();

  Ref<SCMatrixKit> kit = basis->matrixkit();
  RefSCDimension ao_dim(new SCDimension(basis->nbasis()));
  RefSCDimension occ_act_dim(new SCDimension(nocc_act));
  sc::RefSCMatrix Lmatrix(ao_dim, occ_act_dim, kit);
  pack_array_into_matrix(L_,Lmatrix); // Make matrix copy of L

  int nshell = basis->nshell();

  Timer tim("0: L_blocks");
  
  // Create L_blocks
  L_blocks.resize(nshell);
  for (int shell_index=0; shell_index<nshell; shell_index++) {
      sma2::BlockInfo<2> lb;
      lb.block(0) = shell_index;
      lb.block(1) = 0; 
      sma2::BlockInfo<2> ub;
      ub.block(0) = shell_index;
      ub.block(1) = L_.index(1).nblock(); 
      sma2::Array<2>::blockmap_t::const_iterator start = L_.blockmap().lower_bound(lb);
      sma2::Array<2>::blockmap_t::const_iterator end   = L_.blockmap().upper_bound(ub);
      for (sma2::Array<2>::blockmap_t::const_iterator iter = start;
           iter != end; iter++) {
          int occ_index = (iter->first).block(1);
          L_blocks[shell_index].push_back(occ_index);
        }
    }

  int me = msg_->me();
  int nproc = msg_->n();

  tim.change("1: Lshellmax");

  for (int M=0; M<nshell; M++) {
      int size_M = basis->shell(M).nfunction();
      int M_offset = basis->shell_to_function(M);
      for (int iocc=0; iocc<nocc_act; iocc++) {
          double max = 0;
          for (int mu=0; mu<size_M; mu++) {
//*** imbn: do we need Lmatrix? (can we access individual elements of L?)              
              double tmp = fabs(Lmatrix(mu+M_offset,iocc));
              if (tmp > max) max = tmp;
            }
          Lshellmax(M,iocc) = max;
        }
    }

  tim.change("2: L_map");
  
  // Create L_map, and its transpose, Lt_map
  L_map.resize(nshell);
  std::vector<std::multimap<double, int, std::greater<double> > > Lt_map(nocc_act);
  for (int shell_index=0; shell_index<nshell; shell_index++) {
      for (int occ_index=0; occ_index<nocc_act; occ_index++) {
          double lvalue = Lshellmax(shell_index,occ_index);
          L_map[shell_index].insert(std::make_pair(lvalue,occ_index));
          Lt_map[occ_index].insert(std::make_pair(lvalue,shell_index));
        }
    }

  tim.change("3: Lmax");

  for (domainmap_t::iterator i_occpair = domainmap_.begin();
       i_occpair != domainmap_.end();
       i_occpair++) {
      int i = i_occpair->first.first - nfzc_;
      int j = i_occpair->first.second - nfzc_;
      domainmapvirbs_t &virset = i_occpair->second;
      for (int M=0; M<nshell; M++) {
          int size_M = basis->shell(M).nfunction();
          int M_offset = basis->shell_to_function(M);
          double maxval;

          maxval = 0.0;
          for (int mu=0; mu<size_M; mu++) {
              double val = fabs(Lmatrix(mu+M_offset,j)); 
              if (val > maxval) maxval = val;
            }
          if (maxval > Lmax[M][i]) Lmax[M][i] = maxval;

          maxval = 0.0;
          for (int mu=0; mu<size_M; mu++) {
              double val = fabs(Lmatrix(mu+M_offset,i)); 
              if (val > maxval) maxval = val;
            }
          if (maxval > Lmax[M][j]) Lmax[M][j] = maxval;
        }
    }
  Lmax_maxval = 0.0;
  for (int M=0; M<nshell; M++) {
      double maxval = 0.0;
      for (int i=0; i<nocc_act; i++) {
          double val = Lmax[M][i];
          if (val > maxval) maxval = val;
        }
      Lmaxvec[M] = maxval;
      if (maxval > Lmax_maxval) Lmax_maxval = maxval;
    }
  // imbn debug print
  ExEnv::out0() << indent << "Lmax_maxval: " << Lmax_maxval << std::endl;
  // imbn end of debug print

  tim.change("4: Dmax");

  Dmax.resize(nshell);
  for (int i=0; i<nshell; i++) {
      Dmax[i].resize(nshell);
      std::fill(Dmax[i].begin(), Dmax[i].end(), 0.0);
    }

  double dmax_threshold = integral_threshold_
                         /(T_schwarz_maxval
                           *T_schwarz_maxval
                           *Pmax_maxval
                           *Pmax_maxval);
  // i loops over all occupied orbitals
  for (int i=0; i<nocc_act; i++) {
      // j loops over all occupied orbitals paired with i
      for (std::set<int>::iterator j_iter = paired_occ_[i].begin();
           j_iter != paired_occ_[i].end();
           j_iter++) {
          int j = *j_iter;
          // M loops over all shells with Lmax(M,i) large enough to contribute
          for (std::multimap<double, int, std::greater<double> >::iterator
                   M_iter = Lt_map[i].begin();
               M_iter != Lt_map[i].end();
               M_iter++) {
              double L_M_i = M_iter->first;
              if (L_M_i*Lmax_maxval < dmax_threshold) break;
              int M = M_iter->second;
              std::vector<double> &Dmax_M = Dmax[M];
              // N loops over all shells with Lmax(N,j) large enough to contribute
              for (std::multimap<double, int, std::greater<double> >::iterator
                       N_iter = Lt_map[j].begin();
                   N_iter != Lt_map[j].end();
                   N_iter++) {
                  double L_N_j = N_iter->first;
                  double bound = L_M_i*L_N_j;
                  if (bound < dmax_threshold) break;
                  int N = N_iter->second;
                  double &Dmax_MN = Dmax_M[N];
                  if (Dmax_MN < bound) Dmax_MN = bound;
                }
            }
        }
    }
  Dmax_maxval = 0.0;
  for (int i=0; i<nshell; i++) {
      const std::vector<double> &Dmax_i = Dmax[i];
      for (int j=0; j<nshell; j++) {
          double val = Dmax_i[j];
          if (val > Dmax_maxval) Dmax_maxval = val;
        }
    }

  // imbn debug print
  int tmpcount = 0;
  for (int i=0; i<nshell; i++) {
      for (int j=0; j<nshell; j++) {
          if (Dmax[i][j] > 0.00001) {
              tmpcount++;
            }
        }
    }
  ExEnv::out0() << "Dmax: Number of elements greater than 0.00001: " 
                << tmpcount << std::endl;
  ExEnv::out0() << scprintf("Dmax: Fraction of elements greater than"
                            " 0.00001 = % 14.8f", double(tmpcount)/(nshell*nshell))
                << std::endl;
  // imbn end of debug print
}

void
LMP2::compute_PPmax_screening_matrix(std::vector<std::vector<double> > &PPmax,
                                     sc::RefSCMatrix &Pmax,
                                     const std::vector<std::set<int> >&
                                     virb_united_pair_domains,
                                     double T_schwarz_maxval,
                                     double Pmax_maxval,
                                     double Lmax_maxval)
{
  // Construct screening matrix PPmax with elements
  // PPmax(M,R) = max of |P(mu,u)*P(rho,t)| for mu in M, rho in R,
  //           (= Pmax(M,u)*Pmax(R,t) )
  //              and u,t running over all virb pairs with the
  //              restriction that u and t belong to the same domain
  //              (i.e., there exists a pair domain containing both u and t)

  Timer tim("compute_PPmax_screening_matrix");

  sc::Ref<sc::GaussianBasisSet> basis = this->basis();

  int me = msg_->me();
  int nproc = msg_->n();

  int nvir =vir_.nblock();
  int nshell = basis->nshell();

  typedef std::multimap<double,int,std::greater<double> > sorted_vector_t;
  std::vector<sorted_vector_t> Pmaxt_sorted(nshell);

  for (int i=0; i<nvir; i++) {
      sorted_vector_t &Pmaxt_sorted_i = Pmaxt_sorted[i];
      for (int j=0; j<nshell; j++) {
          double Pmax_ji = Pmax(j,i);
          Pmaxt_sorted_i.insert(std::make_pair(Pmax_ji,j));
        }
    }

  double PPmax_threshold = integral_threshold_
                          /(T_schwarz_maxval
                            *T_schwarz_maxval
                            *Lmax_maxval
                            *Lmax_maxval);
    
  PPmax.resize(nshell);
  for (int i=0; i<nshell; i++) {
      PPmax[i].resize(nshell);
      std::fill(PPmax[i].begin(), PPmax[i].end(), 0.0);
    }
  for (int u=0; u<nvir; u++) {
      for (std::set<int>::const_iterator t_iter
               = virb_united_pair_domains[u].begin();
           t_iter != virb_united_pair_domains[u].end();
           t_iter++) {
          int t = *t_iter;
          for (sorted_vector_t::const_iterator M_iter = Pmaxt_sorted[u].begin();
               M_iter != Pmaxt_sorted[u].end();
               M_iter++) {
              double Pmax_Mu = M_iter->first;
              if (Pmax_Mu*Pmax_maxval < PPmax_threshold) break;
              int M = M_iter->second;
              std::vector<double> &PPmax_M = PPmax[M];
              for (sorted_vector_t::const_iterator R_iter = Pmaxt_sorted[t].begin();
                   R_iter != Pmaxt_sorted[t].end();
                   R_iter++) {
                  double Pmax_Rt = R_iter->first;
                  if (Pmax_Mu*Pmax_Rt < PPmax_threshold) break;
                  int R = R_iter->second;
                  double val = Pmax_Mu*Pmax_Rt;
                  double &PPmax_MR = PPmax_M[R];
                  if (PPmax_MR < val) PPmax_MR = val;
                }
            }
        }
    }
  
  // imbn debug print
  int tmpcount = 0;
  for (int i=0; i<nshell; i++) {
      for (int j=0; j<nshell; j++) {
          if (fabs(PPmax[i][j]) > 0.00001) {
              tmpcount++;
            }
        }
    }
  ExEnv::out0() << "PPmax: Number of elements greater than 0.00001: " 
                << tmpcount << std::endl;
  ExEnv::out0() << scprintf("PPmax: Fraction of elements greater than"
                            " 0.00001 = % 14.8f", double(tmpcount)/(nshell*nshell))
                << std::endl;
  // imbn end of debug print
  
}

static void
allocate_T(sma2::Array<4> &T, domainmap_t &domainmap, int nfzc)
{
  sma2::BlockInfo<4> T_blockinfo;
  for (domainmap_t::const_iterator i_occpair = domainmap.begin();
       i_occpair != domainmap.end(); i_occpair++) {

      int i = i_occpair->first.first - nfzc;
      int j = i_occpair->first.second - nfzc;
      T_blockinfo.block(0) = i;
      T_blockinfo.block(1) = j;
      std::set<int> domain_virbs = i_occpair->second;
      
      for (std::set<int>::const_iterator virb_iter = domain_virbs.begin();
           virb_iter != domain_virbs.end(); virb_iter++) {

          T_blockinfo.block(2) = *virb_iter;
          
          for (std::set<int>::const_iterator virb_iter_2 = domain_virbs.begin();
               virb_iter_2 != domain_virbs.end(); virb_iter_2++) {

              T_blockinfo.block(3) = *virb_iter_2;
          
              // allocate current block
              T.add_new_unallocated_block(T_blockinfo);
            }
        }
    }
}

static void
allocate_R(sma2::Array<4> &R, const domainmap_t &domainmap, int nfzc,
           const my_occ_pairs_t &my_occ_pairs)
{
  sma2::BlockInfo<4> R_blockinfo;
  for (my_occ_pairs_t::const_iterator i_occpair = my_occ_pairs.begin();
       i_occpair != my_occ_pairs.end(); i_occpair++) {

      int i = i_occpair->first - nfzc;
      int j = i_occpair->second - nfzc;
      R_blockinfo.block(0) = i;
      R_blockinfo.block(1) = j;

      domainmap_t::const_iterator domain_iter = domainmap.find(*i_occpair);
      if (domain_iter == domainmap.end()) {
          throw sc::ProgrammingError("allocate_R: no entry in domainmap",
                                     __FILE__, __LINE__);
        }
      const std::set<int> &domain_virbs = domain_iter->second;
      
      for (std::set<int>::const_iterator virb_iter = domain_virbs.begin();
           virb_iter != domain_virbs.end(); virb_iter++) {

          R_blockinfo.block(2) = *virb_iter;
          
          for (std::set<int>::const_iterator virb_iter_2 = domain_virbs.begin();
               virb_iter_2 != domain_virbs.end(); virb_iter_2++) {

              R_blockinfo.block(3) = *virb_iter_2;
          
              // allocate current block
              R.add_new_unallocated_block(R_blockinfo);
            }
        }
    }
}

double
LMP2::compute_lmp2_energy()
{
  if(molecule()->point_group()->char_table().order() != 1) {
      throw std::runtime_error("LMP2 supports C1 symmetry only");
    }

  sc::Timer tim;
  
  tim.enter("LMP2");

  int me = msg_->me();
  int nproc = msg_->n();

  double tolerance = desired_value_accuracy();
  set_actual_value_accuracy(1.0);

  int nbasis = basis()->nbasis();
  int nshell = basis()->nshell();
  int nocc = ref_->nelectron()/2;
  int nocc_act = nocc - nfzc_;
  int natom = molecule()->natom();

  sc::MolecularFormula mf(molecule());

  bound_ = tolerance/threshold_factor_;

  ExEnv::out0() << std::endl;
  ExEnv::out0() << indent << "Beginning local correlation calculation"
                << std::endl;
  LCorr::print_parameters();
  ExEnv::out0() << indent << "max_iter           = " << max_iter_ << std::endl;
  ExEnv::out0() << indent << "nfunction          = " << nbasis << std::endl;
  ExEnv::out0() << indent << "nshell             = " << nshell << std::endl;
  ExEnv::out0() << indent << "natom              = " << natom
                << std::endl;
  ExEnv::out0() << indent << "nocc               = " << nocc
                << std::endl;
  ExEnv::out0() << indent << "nfzc               = " << nfzc_
                << std::endl;
  ExEnv::out0() << indent << "S_threshold        = " << S_threshold_ << std::endl;
  ExEnv::out0() << indent << "integral_threshold = " << integral_threshold_ << std::endl;
  ExEnv::out0() << indent << "q1_threshold       = " << q1_threshold_ << std::endl;
  ExEnv::out0() << indent << "q2_threshold       = " << q2_threshold_ << std::endl;
  ExEnv::out0() << indent << "q3_threshold       = " << q3_threshold_ << std::endl;
  ExEnv::out0() << indent << "q4_threshold       = " << q4_threshold_ << std::endl;
  ExEnv::out0() << indent << "tolerance          = " << tolerance << std::endl;
  ExEnv::out0() << indent << "threshold_factor   = " << threshold_factor_
                << std::endl;
  ExEnv::out0() << indent << "extrap_T:" << std::endl;
  ExEnv::out0() << incindent;
  extrap_T_->print();
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "complete_domains       = "
                << (completedomains_?"yes":"no")
                << std::endl;
  ExEnv::out0() << indent << "single_virtual_block   = "
                << (singlevirb_?"yes":"no")
                << std::endl;
  ExEnv::out0() << indent << "i_ge_j                 = "
                << (i_ge_j_?"yes":"no")
                << std::endl;
  ExEnv::out0() << indent << "m_ge_n                 = "
                << (m_ge_n_?"yes":"no")
                << std::endl;
  ExEnv::out0() << indent << "minimize_q2            = "
                << (minimize_q2_?"yes":"no")
                << std::endl;
  ExEnv::out0() << indent << "always_use_dist_t      = "
                << (always_use_dist_t_?"yes":"no")
                << std::endl;
  ExEnv::out0() << indent << "occ_orbitals           = "
                << occ_orbitals_
                << std::endl;
  ExEnv::out0() << indent << "vir_orbitals           = "
                << vir_orbitals_
                << std::endl;
  ExEnv::out0() << indent << "distance_threshold = " << distance_threshold_
                << " bohr" << std::endl;
  ExEnv::out0() << indent << "completeness_threshold = " << completeness_threshold_
                << std::endl;
  ExEnv::out0() << indent << "molecule:          " << mf.formula() << std::endl;
  ExEnv::out0() << std::endl;
  
  // NB: must use ShellBlocking (required for the AO integrals)
  ao_.init(basis(), sma2::Range::ShellBlocking);

  if (vir_orbitals_ == "canonical") {
      vir_.init(nbasis-nocc, nbasis-nocc);
    }
  else if (singlevirb_) {
      vir_.init(nbasis, nbasis);
    }
  else {
      vir_.init(basis(), sma2::Range::AtomBlocking);
    }

  init_virb_to_bfns(vir_);

  occ_act_.init(nocc_act, 1);

  sc::RefSymmSCMatrix ao_density = ref_->ao_density();

  tim.enter("Overlap");
  sc::RefSymmSCMatrix ao_overlap
      = ref_->integral()->petite_list(basis())->to_AO_basis(ref_->overlap());
  tim.exit("Overlap");

  sc::RefSCMatrix scf_local;
  if (occ_orbitals_ == "pipek-mezey") {
      msg_->sync(); // imbn test
      tim.enter("Localization");
      scf_local = pop_local_mo(ref_, nfzc_, ao_overlap, msg_);
      tim.exit("Localization");
    }
  else if (occ_orbitals_ == "canonical") {
      sc::RefSCMatrix scf_local_all;
      scf_local_all
          = (ref_->so_to_mo()*(ref_->integral()->petite_list()->aotoso()).t()).t();
      scf_local = convert_complete_to_occupied_vector_nosymm(ref_, nfzc_, scf_local_all);
    }
  else {
      throw std::runtime_error("bad occ_orbitals value");
    }


  L_.init(ao_,occ_act_,"L",bound_);
  sma2::pack_matrix_into_empty_array(scf_local,L_,bound_);

  if (extrap_T_) extrap_T_->set_tolerance(tolerance);

  // Create the transformation matrix P (P = I - ao_density*ao_overlap)
  // For projected atomic virtual orbitals, P is the matrix transforming from the
  // AO to the projected AO basis (projecting out the occupied space)
  // This matrix is called R tilde in ref [2].
  // For canonical virtual orbitals, P is the the AO to MO transformation
  tim.enter("Create P");
  P_.init(ao_,vir_,"P",bound_);
  if (vir_orbitals_ == "old_projected_atomic") {
      sc::RefSCMatrix Pmat(ao_overlap.dim(),ao_overlap.dim(),ao_density.kit());
      tim.enter("compute Pmat");
      Pmat = -0.5*ao_density*ao_overlap;
      Pmat->shift_diagonal(1.0);
      tim.exit("compute Pmat");
      tim.enter("pack Pmat");
      sma2::pack_matrix_into_empty_array(Pmat,P_,bound_);
      tim.exit("pack Pmat");
    }
  else if (vir_orbitals_ == "projected_atomic") {
      sc::RefSymmSCMatrix hole_occupation(ref_->oso_dimension(),
                                          ref_->basis_matrixkit());
      hole_occupation.assign(0.0);
      for (int i=0; i<ref_->oso_dimension()->n(); i++) {
          if (fabs(ref_->occupation(i)) < 1e-14) {
              hole_occupation(i,i) = 1.0;
            }
        }
      sc::RefSCMatrix oso_vector = ref_->oso_eigenvectors();
      sc::RefSCMatrix vector = ref_->so_to_orthog_so().t() * oso_vector;
      sc::RefSymmSCMatrix hole_density(ref_->so_dimension(), ref_->basis_matrixkit());
      hole_density.assign(0.0);
      hole_density.accumulate_transform(vector, hole_occupation);
      sc::RefSCMatrix Pmat = hole_density * ao_overlap;
      sma2::pack_matrix_into_empty_array(Pmat,P_,bound_);
    }
  else if (vir_orbitals_ == "canonical") {
      sc::RefSCMatrix scf_vec
          = (ref_->so_to_mo()*(ref_->integral()->petite_list()->aotoso()).t()).t();
      // This assumes all orbitals are in block 0, hence C1 symmetry.
      sc::BlockedSCMatrix *scf_vec_blocked
          = dynamic_cast<sc::BlockedSCMatrix*>(scf_vec.pointer());
      sc::RefSCMatrix scf_vir
          = scf_vec_blocked->block(0)->get_subblock(0, nbasis-1, nocc, nbasis-1);
      sma2::pack_matrix_into_empty_array(scf_vir,P_,bound_);
    }
  else {
      throw std::runtime_error("LMP2: bad vir_orbitals");
    }
  tim.exit("Create P");

  tim.enter("Create Fock matrices ");

  // Create the Fock matrix (F_ao) in the atomic orbital basis
  sma2::Array<2> F_ao(ao_,"F_ao",bound_);
  sma2::pack_matrix_into_empty_array(ref_->fock(0),F_ao,bound_);

  // Transform the AO Fock matrix to the localized occupied MO basis (F_ao -> F_occ)
  F_occ_.init(occ_act_,"F_occ",bound_);
  transform_array(F_ao, L_, L_, F_occ_, msg_);

  F_diag_.resize(occ_act_.nindex());
  sma2::extract_diagonal(F_occ_, F_diag_);

  // Transform the AO Fock matrix to the projected AO basis (F_ao -> F_vir)
  F_vir_.init(vir_,"F_vir",bound_);
  transform_array(F_ao, P_, P_, F_vir_, msg_);
  
  tim.exit("Create Fock matrices ");

  tim.enter("Transform overlap ");
  // Create the overlap matrix (S) in the projected AO basis
  sma2::Array<2> S_tmp;
  S_tmp.init(vir_,"S_tmp",bound_);
  S_.init(vir_,"S",bound_);
  sma2::Array<2> S_ao(ao_,"S_ao",bound_);
  sma2::pack_matrix_into_empty_array(ao_overlap,S_ao,bound_);
  transform_array(S_ao, P_, P_, S_tmp, msg_);

  // imbn test code
  analyze_array(S_tmp,"S_tmp",msg_,false);
  S_.assign_tol(S_tmp, S_threshold_);
  analyze_array(S_,"S",msg_,false);
  // imbn end of test code
  tim.exit("Transform overlap ");

  //****IMBN: Need to figure out the proper dimensions for various matrices;
  //  the dimensions ao_density.dim() and ao_overlap.dim() should be the same (?)
  //  but the blocking methods for the matrices might be different?

  tim.enter("Domain map");
  std::vector<std::vector<int> > domains;
  domains.resize(nocc_act); // need nocc_act domain vectors
  domainmap_.clear();

  if (completedomains_) {
      domainmapvirbs_t virbs;
      for (int virb=0; virb<vir_.nblock(); virb++) virbs.insert(virb);
      for (int i=nfzc_; i<nocc; i++) {
          domainmapvirbs_t::iterator tmp_iter;
          int count = 0;
          for (tmp_iter = virbs.begin(); tmp_iter != virbs.end();
               tmp_iter++) {
              domains[i-nfzc_].push_back(*tmp_iter);
              count++;
            }
          int jfence = i_ge_j_?i+1:nocc;
          for (int j=nfzc_; j<jfence; j++) {
              domainmap_[std::make_pair(i,j)] = virbs;
            }
        }
      ExEnv::out0() << indent << "Using complete domains" << std::endl;
    }
  else {
      create_domains(ref_, nfzc_, scf_local, domains,
                     domainmap_, distance_threshold_, completeness_threshold_,
                     !i_ge_j_, dist_ij_, S_ao, L_, bound_, msg_);
    }
  tim.exit("Domain map");

  // Create the vector paired_occ containing for each
  // active occupied mo, i, the set of mo's j for which the ij pair
  // is included in T and K_2occ.
  // NB: No i<->j symmetry is used in occ_pairs, i.e., for each i,
  //     all j's for which ij or ji is a valid pair are included
  // NB: occ. indices run from 0 to nocc_act
  tim.enter("Create paired_occ");
  paired_occ_.resize(nocc_act);
  for (domainmap_t::iterator i_occpair = domainmap_.begin();
       i_occpair != domainmap_.end(); i_occpair++) {
      int mo1 = i_occpair->first.first - nfzc_;
      int mo2 = i_occpair->first.second - nfzc_;
      paired_occ_[mo1].insert(mo2);
      paired_occ_[mo2].insert(mo1);
    }
  tim.exit("Create paired_occ");

  // imbn debug print
  // Print out pair domains:
  ExEnv::out0() << indent << "Pair Domains:" << std::endl;
  domainmap_t::iterator tmp_iter;
  for (tmp_iter = domainmap_.begin(); tmp_iter != domainmap_.end();
       tmp_iter++) {
      ExEnv::out0() << indent
                    << "  " << tmp_iter->first.first
                    << ", " << tmp_iter->first.second << ":";
      std::set<int>::iterator tmp2_iter;
      for (tmp2_iter = tmp_iter->second.begin();
           tmp2_iter != tmp_iter->second.end();
           tmp2_iter++) {
          ExEnv::out0() << " " << *tmp2_iter;
        }
      ExEnv::out0() << std::endl;
    }
  // imbn end of debug print

  
  ref_->integral()->set_storage(150000000);  // increase this to avoid recomputation of intermediates
//  ref_->integral()->set_storage(50000000);
  sc::Ref<TwoBodyInt> tbint = ref_->integral()->electron_repulsion();
  tbint->set_integral_storage(0);

  // Pre-compute quantities to be used for Schwarz screening:
  RefSCDimension shell_dim(new SCDimension(nshell));
  Ref<SCMatrixKit> kit = basis()->matrixkit();
  sc::RefSCMatrix T_schwarz(shell_dim, shell_dim, kit);
  sc::RefSCVector T_schwarz_maxvec(shell_dim, kit);
  double T_schwarz_maxval;
  compute_Schwarz_screening_quantities(tbint, T_schwarz, T_schwarz_maxvec,
                                       T_schwarz_maxval, nshell);

  // Pre-compute Pmax and Pmaxvec to be used for screening
  RefSCDimension ao_dim(new SCDimension(nbasis));
  RefSCDimension vir_block_dim(new SCDimension(vir_.nblock()));
  sc::RefSCMatrix Pmax(shell_dim, vir_block_dim, kit);
  sc::RefSCMatrix Pmatrix(ao_dim, ao_dim, kit);
  sc::RefSCVector Pmaxvec(shell_dim, kit);
  pack_array_into_matrix(P_,Pmatrix); // Make matrix copy of P
  double Pmax_maxval;
  compute_P_screening_quantities(Pmax, Pmatrix, Pmaxvec, Pmax_maxval);

  // Create a map containing the virb_united_pair_domain for each virb
  tim.enter("create_virb_united_pair_domains");
  std::vector<std::set<int> > virb_united_pair_domains;
  create_virb_united_pair_domains(vir_.nblock(), domainmap_,
                                  virb_united_pair_domains);
  tim.exit("create_virb_united_pair_domains");
  
  // Pre-compute L and D quantities used for screening
  std::vector<std::vector<double> > Dmax;
  double Dmax_maxval;
//imbn  sc::RefSCMatrix Lmax(shell_dim, occ_act_dim, kit);
//imbn  Lmax.assign(0.0); // Lmax must be initialized to zero
  std::vector<std::vector<double> > Lmax;
  Lmax.resize(nshell);
  for (int shellindex=0; shellindex<nshell; shellindex++) {
      Lmax[shellindex].resize(nocc_act);
      for (int occindex=0; occindex<nocc_act; occindex++) {
          Lmax[shellindex][occindex] = 0.0;
        }
    }
  std::vector<double> Lmaxvec(nshell);
  std::vector<std::vector<int> > L_blocks;
  std::vector<std::multimap<double,int,std::greater<double> > > L_map;

  RefSCDimension occ_act_dim(new SCDimension(nocc_act));
  sc::RefSCMatrix Lshellmax(shell_dim, occ_act_dim, kit);
  double Lmax_maxval;
  compute_L_and_D_screening_quantities(Dmax, Dmax_maxval,
                                       Lmax, Lmaxvec, Lshellmax,
                                       L_blocks, L_map,
                                       Lmax_maxval,
                                       T_schwarz_maxval,
                                       Pmax_maxval);

  // Pre-compute PPmax screening matrix
  std::vector<std::vector<double> > PPmax;
  compute_PPmax_screening_matrix(PPmax, Pmax, virb_united_pair_domains,
                                 T_schwarz_maxval,
                                 Pmax_maxval,
                                 Lmax_maxval);

  // Create vector my_occ_pairs of local occ-occ pairs (used to
  // distribute work in the iterative procedure)
  std::set<std::pair<int,int> > included_occ_pairs;
  for (domainmap_t::const_iterator i_occpair = domainmap_.begin();
       i_occpair != domainmap_.end(); i_occpair++) {
      included_occ_pairs.insert(std::make_pair(i_occpair->first.first-nfzc_,
                                               i_occpair->first.second-nfzc_));
    }
  pair_mapping_ = new sma2::PairMapping(msg_,included_occ_pairs);
  std::set<std::pair<int,int> > my_occ_pairs_tmp;
  pair_mapping_->local_pairs(my_occ_pairs_tmp); // my_occ_pairs_tmp now contains the pairs to be held locally
  my_occ_pairs_.clear();
  // (convert to a vector)
  for (std::set<std::pair<int,int> >::iterator tmp_iter = my_occ_pairs_tmp.begin();
       tmp_iter != my_occ_pairs_tmp.end(); tmp_iter++) {
      int i=tmp_iter->first;
      int j=tmp_iter->second;
      my_occ_pairs_.push_back(std::make_pair(i+nfzc_, j+nfzc_));
    }

  compute_K_2occ(T_schwarz, T_schwarz_maxvec,
                 T_schwarz_maxval, Pmax, Pmaxvec,
                 Lmax, Lshellmax, PPmax, Dmax, Dmax_maxval, Lmaxvec, L_map);

  // Allocate the replicated T
  if (!always_use_dist_t_) {
      T_.init(occ_act_,occ_act_,vir_,vir_);
      allocate_T(T_, domainmap_, nfzc_);
      analyze_array(T_,"T",msg_,false);
      T_.allocate_blocks();
      if (i_ge_j_) {
          sma2::IndexList T_jirs_indexlist(1,0,2,3);
          remap(T_jirs_,T_,T_jirs_indexlist);
        }
    }

  // Allocate the distributed T. 
  T_local_.init(occ_act_,occ_act_,vir_,vir_);
  T_local_("i","j","r","s") |= K_2occ_("i","j","r","s");
  T_local_.zero();
  analyze_array(T_local_,"T_local",msg_,true);

  T_n_element_ = T_local_.n_element();
  msg_->sum(T_n_element_);

  // Compute W (and assign unique_eigvals)
  tim.enter("Compute W");
  compute_doubles_W();
  tim.exit("Compute W");
  ExEnv::out0() << indent
                << "Number of unique W's for doubles: "
                << n_unique_W()
                << std::endl;
  

  double ecorr = 0.0;

  ecorr = iterate_LMP2_equations(tolerance/10.0, tolerance);

  // need to check that it actually converged
  set_actual_value_accuracy(tolerance);

  tim.exit("LMP2");

  // analyze occupied orbitals
  if (analyze_occ_orbs_) {
    tim.enter("Analyze occ orbs");
    analyze_occ_orbs(scf_local);
    tim.exit("Analyze occ orbs");
  }

  return ecorr;

}

void
LMP2::compute_delta_T(sma2::Array<4> &Res,
                      sma2::Array<4> &delta_T)
{
  // Compute the updates delta_T to T.
  // The updates are computed as outlined in Ref. [2].
  // First Res is transformed to Res_bar, then delta_T_bar is computed
  // from Res_bar, and, finally, delta_T is computed by
  // back-transformation of delta_T_bar.

  // Parallelize computation of delta_T by looping over locally held ij pairs

  delta_T.zero();
  Ref<SCMatrixKit> kit = basis()->matrixkit();

  int nproc = msg_->n();
  int me = msg_->me();

  Timer occ_pairs_timer("occ pairs loop");
  for (my_occ_pairs_t::const_iterator i_occpair = my_occ_pairs_.begin();
       i_occpair != my_occ_pairs_.end(); i_occpair++) {
      Timer occ_pairs_iter_timer("0: setup");

      std::pair<int,int> current_pair = *i_occpair;
      int i = current_pair.first - nfzc_;
      int j = current_pair.second - nfzc_;

      sma2::Index iI(i), iJ(j); // LMO basis (fixed indices)
      sma2::Index iR("r"), iS("s"); // PAO virtual basis
      sma2::Index iX("x"), iY("Y"); // The nonredundant virtual basis

      // Extract W_ij for current pair
#if COMPLETE_W // This uses a W and eigvals for the full virtual space
      domainmapvirbs_t complete_virset;
      for (int i=0; i<vir_.nblock(); i++) {
          complete_virset.insert(i);
        }
      sma2::Array<2> &W_complete = unique_W(complete_virset);
      std::vector<double> &eigvals_complete
          = unique_eigvals(complete_virset);
      sma2::Array<2> &W_ij = W_complete;
      std::vector<double> &eigvals_ij = eigvals_complete;
#else  // This is the usual way of doing things:
      sma2::Array<2> &W_ij
          = unique_W(domainmap_[current_pair]);
      std::vector<double> &eigvals_ij
          = unique_eigvals(domainmap_[current_pair]);
#endif

      sma2::Range ij_vir_nonred(W_ij.index(1));

      // Transform Res_ij to pseudocanonical basis from redundant PAO basis
      occ_pairs_iter_timer.change("1: PAO->PC");
      sma2::Array<2> RT_tmp_1(vir_,ij_vir_nonred);
      RT_tmp_1(iR,iX) |= Res(iI,iJ,iR,iS) * W_ij(iS,iX);
      RT_tmp_1(iR,iX)  = Res(iI,iJ,iR,iS) * W_ij(iS,iX);
      sma2::Array<2> T_bar_ij(ij_vir_nonred,ij_vir_nonred);
      T_bar_ij(iX,iY) |= RT_tmp_1(iR,iY) * W_ij(iR,iX);
      T_bar_ij(iX,iY)  = RT_tmp_1(iR,iY) * W_ij(iR,iX);

      // Apply denominator
      occ_pairs_iter_timer.change("2: denom");
      double F_diag_ij = F_diag_[i] + F_diag_[j];
      std::vector<const std::vector<double>*> denoms_ij;
      denoms_ij.push_back(&eigvals_ij);
      denoms_ij.push_back(&eigvals_ij);
      apply_denominator(T_bar_ij,-F_diag_ij,denoms_ij);
      denoms_ij.clear();
      T_bar_ij *= -1.0;

      occ_pairs_iter_timer.change("3: PC->PAO");
      RT_tmp_1.clear();
      RT_tmp_1(iR,iX) |= W_ij(iR,iY) * T_bar_ij(iY,iX);
      RT_tmp_1(iR,iX)  = W_ij(iR,iY) * ~T_bar_ij(iY,iX);
      sma2::Array<2> newT_ij(vir_,vir_);
      newT_ij(iR,iS) |= W_ij(iS,iX) * RT_tmp_1(iR,iX);
      newT_ij(iR,iS)  = W_ij(iS,iX) * ~RT_tmp_1(iR,iX);

      // Assign newT to T
      occ_pairs_iter_timer.change("4: assign");
      delta_T(iI,iJ,iR,iS).skip_bounds_update() += newT_ij(iR,iS);
    }
  occ_pairs_timer.exit();

//imbn 05/10/06  delta_T.parallel_accumulate(msg_);  // don't need this now that delta_T is distributed
}

void
LMP2::compute_LMP2_residual_SFT(sma2::Array<4> &R)
{
  Timer overalltim("SFT");
  Timer tim("init");
  sma2::Index i("i"), j("j"), k("k"), r("r"), s("s"), p("p");

  for (my_occ_pairs_t::iterator iter = my_occ_pairs_.begin();
       iter != my_occ_pairs_.end();
       iter++) {
      sma2::Index ifix(iter->first - nfzc_);
      sma2::Index jfix(iter->second - nfzc_);

      // Allocate tmp
      tim.change("tmp in first loop");
      sma2::Array<2> tmp(vir_,vir_);
      std::pair<int,int> current_pair = *iter;
      const std::set<int> &current_pair_domain = domainmap_[current_pair];
      std::set<int>::const_iterator domain_iter;
      for (domain_iter = current_pair_domain.begin(); domain_iter != current_pair_domain.end();
           domain_iter++) {
          sma2::BlockInfo<2> tmp_blockinfo;
          int tmp_index = *domain_iter;
          tmp_blockinfo.block(0) = tmp_index;
          std::set<int>::const_iterator domain_iter2;
          for (domain_iter2 = current_pair_domain.begin(); domain_iter2 != current_pair_domain.end();
               domain_iter2++) {
              int tmp_index2 = *domain_iter2;
              tmp_blockinfo.block(1) = tmp_index2;
              tmp.add_unallocated_block(tmp_blockinfo);
            }
        }
      
      tmp(r,s) = T_local_(ifix,jfix,r,p) * S_(p,s);
      R(ifix,jfix,r,s).skip_bounds_update() += tmp(p,s) * F_vir_(r,p);

      tmp(r,s) = T_local_(ifix,jfix,r,p) * F_vir_(p,s);
      R(ifix,jfix,r,s).skip_bounds_update() += tmp(p,s) * S_(r,p);
    }
}

void
LMP2::compute_LMP2_residual_SFTS(sma2::Array<4> &R)
{
  if (always_use_dist_t_) {
      compute_LMP2_residual_SFTS_i(R);
    }
  else {
      compute_LMP2_residual_SFTS_ij(R);
    }
}

void
LMP2::compute_LMP2_residual_SFTS_ij(sma2::Array<4> &R)
{
  Timer overalltim("SFTS_ij");
  Timer tim("init");
  sma2::Index i("i"), j("j"), k("k"), r("r"), s("s"), p("p");

  T_.zero();
  T_(i,j,r,s) = T_local_(i,j,r,s);
  T_.parallel_accumulate(msg_);

  for (my_occ_pairs_t::iterator iter = my_occ_pairs_.begin();
       iter != my_occ_pairs_.end();
       iter++) {
      sma2::Index ifix(iter->first - nfzc_);
      sma2::Index jfix(iter->second - nfzc_);

      tim.change("tmp1 in second loop - alloc");
      sma2::Array<2> tmp1(vir_,vir_);
      tmp1(r,s) |= T_(ifix,k,r,s)*F_occ_(jfix,k);
      if (i_ge_j_) {
          tmp1(r,s) |= T_jirs_(ifix,k,s,r)*F_occ_(jfix,k);
        }

      tim.change("tmp1 in second loop - assign");
      tmp1(r,s) = T_(ifix,k,r,s)*F_occ_(jfix,k);
      if (i_ge_j_) {
          tmp1(r,s) += T_jirs_(ifix,k,s,r)*F_occ_(jfix,k);
          tmp1(r,s) -= T_(ifix,ifix,r,s)*F_occ_(jfix,ifix);
        }

#if 0
      tim.change("tmp2 in second loop");
      sma2::Array<2> tmp2(vir_,vir_);
      tmp2(r,s) |= tmp1(r,p) * S_(p,s);
      tmp2(r,s) = ~tmp1(r,p) * S_(p,s);
#endif

      tim.change("tmp3 in second loop - alloc");
      sma2::Array<2> tmp3(vir_,vir_);
      std::pair<int,int> current_pair = *iter;
      const std::set<int> &current_pair_domain = domainmap_[current_pair];
      std::set<int>::const_iterator domain_iter;
      for (domain_iter = current_pair_domain.begin(); domain_iter != current_pair_domain.end();
           domain_iter++) {
          sma2::BlockInfo<2> tmp_blockinfo;
          int tmp_index = *domain_iter;
          tmp_blockinfo.block(0) = tmp_index;
          std::set<int>::const_iterator domain_iter2;
          for (domain_iter2 = current_pair_domain.begin(); domain_iter2 != current_pair_domain.end();
               domain_iter2++) {
              int tmp_index2 = *domain_iter2;
              tmp_blockinfo.block(1) = tmp_index2;
              tmp3.add_unallocated_block(tmp_blockinfo);
            }
        }
//imbn test code
      tim.change("tmp2 in second loop - alloc");
      sma2::Array<2> tmp2(vir_,vir_);
//imbn: this is slow  tmp2(p,s) |= S_(p,r) * tmp3(r,s);
//imbn: this is slow  tmp2(r,s) |= tmp1(r,p) * S_(p,s);
//imbn: this is ~same speed  tmp2(p,s) |= R(ifix,jfix,r,s) * S_(r,p);
//imbn: this is a little slower  tmp2(p,s) |= tmp3(r,s) * S_(p,r);
      tmp2(p,s) |= tmp3(r,s) * S_(r,p);
      tim.change("tmp2 in second loop - assign");
      tmp2(r,s) = ~tmp1(r,p) * S_(p,s);
      tim.change("tmp3 in second loop - assign");
//imbn end of test code
      tmp3(r,s) = tmp2(p,s) * S_(r,p);

      tim.change("R update in second loop");
      R(ifix,jfix,r,s).skip_bounds_update() += -1.0 * tmp3(r,s);

      tim.change("tmp1 in second loop");
      // Reallocate and assign tmp1
      tmp1(r,s) |= T_(jfix,k,s,r)*F_occ_(ifix,k);
      if (i_ge_j_) {
          tmp1(r,s) |= T_jirs_(jfix,k,r,s)*F_occ_(ifix,k);
        }

      tmp1(r,s) = T_(jfix,k,s,r)*F_occ_(ifix,k);
      if (i_ge_j_) {
          tmp1(r,s) += T_jirs_(jfix,k,r,s)*F_occ_(ifix,k);
          tmp1(r,s) -= T_(jfix,jfix,s,r)*F_occ_(ifix,jfix);
        }

      tim.change("tmp2 in second loop - B");
      // Reallocate and assign tmp2
#if 0
      tmp2(r,s) |= tmp1(r,p) * S_(p,s);
#endif
      tmp2(r,s) = ~tmp1(r,p) * S_(p,s);

      tim.change("tmp3 in second loop - assign");
      tmp3(r,s) = ~tmp2(p,s) * S_(r,p);

      tim.change("R update in second loop");
      R(ifix,jfix,r,s).skip_bounds_update() += -1.0 * tmp3(r,s);
    }
}

void
LMP2::compute_LMP2_residual_SFTS_i(sma2::Array<4> &R)
{
  Timer overalltim("SFTS_i");
  Timer tim("0: init");

  sma2::Index i("i"), j("j"), k("k"), r("r"), s("s"), p("p"), q("q");

  sma2::Array<4> T_local_all_ij(occ_act_,occ_act_,vir_,vir_);
  T_local_all_ij(i,j,r,s) |= T_local_(i,j,r,s);
  T_local_all_ij(j,i,r,s) |= T_local_(i,j,r,s);
  T_local_all_ij.zero();
  for (my_occ_pairs_t::iterator iter = my_occ_pairs_.begin();
       iter != my_occ_pairs_.end();
       iter++) {
      sma2::Index ifix(iter->first - nfzc_);
      sma2::Index jfix(iter->second - nfzc_);
      T_local_all_ij(ifix,jfix,r,s) += T_local_(ifix,jfix,r,s);
      if (i_ge_j_ && ifix.value() != jfix.value()) {
          T_local_all_ij(jfix,ifix,r,s) += T_local_(ifix,jfix,s,r);
        }
    }

  tim.change("1: T_local_i");
  sma2::Array<4> T_local_i(occ_act_,occ_act_,vir_,vir_);
  sma2::CompleteBlockDistrib<4> distrib(T_local_i,msg_,0);
  T_local_i.distributed_from_distributed(msg_,distrib,T_local_all_ij,
                                         false,false);
  
  sma2::Array<4> TS_a(occ_act_,occ_act_,vir_,vir_);
  sma2::Array<4> STS_a(occ_act_,occ_act_,vir_,vir_);
  sma2::Array<4> TS_b(occ_act_,occ_act_,vir_,vir_);
  sma2::Array<4> STS_b(occ_act_,occ_act_,vir_,vir_);

  tim.change("2: R_local_i alloc");
  sma2::Array<4> R_local_i(occ_act_,occ_act_,vir_,vir_);
  R_local_i(i,j,r,s) |= T_local_i(i,j,r,s);
  R_local_i.zero();
  
  TS_a(i,k,p,s) |= T_local_i(i,k,p,q) * S_(q,s);
  STS_a(i,k,r,s) |= R_local_i(i,j,r,s) * F_occ_(k,j);

  tim.change("3: R_local_i compute A");
  TS_a(i,k,p,s) = T_local_i(i,k,p,q) * S_(q,s);
  tim.change("3: R_local_i compute B");
  STS_a(i,k,r,s) = TS_a(i,k,p,s) * S_(r,p);
  tim.change("3: R_local_i compute C");
  R_local_i(i,j,r,s) -= STS_a(i,k,r,s) * F_occ_(k,j);

  tim.change("4: R_local_j alloc");
  sma2::Array<4> R_local_j(occ_act_,occ_act_,vir_,vir_);
  R_local_j(i,j,r,s) |= T_local_i(j,i,s,r);
  R_local_j.zero();

  TS_b(j,k,p,s) |= T_local_i(j,k,q,p) * S_(q,s);
  STS_b(j,k,r,s) |= R_local_j(i,j,r,s) * F_occ_(i,k);

  tim.change("5: R_local_j compute A");
  TS_b(j,k,p,s) = T_local_i(j,k,q,p) * S_(q,s);
  tim.change("5: R_local_j compute B");
  STS_b(j,k,r,s) = TS_b(j,k,p,s) * S_(r,p);
  tim.change("5: R_local_j compute C");
  R_local_j(i,j,r,s) -= STS_b(j,k,r,s) * F_occ_(i,k);

  tim.change("6: R_contrib alloc");
  sma2::PairBlockDistrib<4> R_distrib(0,1,pair_mapping_);

  sma2::Array<4> R_contrib(occ_act_,occ_act_,vir_,vir_);
  R_contrib(i,j,r,s) |= R(i,j,r,s);

  tim.change("7: R_contrib dist");
  R_contrib.distributed_from_distributed(msg_,R_distrib,R_local_i,
                                         false,true);
  tim.change("8: R_contrib assign");
  R(i,j,r,s) += R_contrib(i,j,r,s);

  tim.change("9: R_contrib dist");
  R_contrib.distributed_from_distributed(msg_,R_distrib,R_local_j,
                                         false,true);
  tim.change("7: R");
  R(i,j,r,s) += R_contrib(i,j,r,s);
}

void
LMP2::compute_LMP2_residual(sma2::Array<4> &R)
{

  sma2::Index i("i"), j("j"), r("r"), s("s");

  R.zero();

  Timer tim("LMP2 residual");

  compute_LMP2_residual_SFT(R);

  compute_LMP2_residual_SFTS(R);

  tim.change("K_2occ update");

  R(i,j,r,s) += K_2occ_(i,j,r,s);
}

double
LMP2::iterate_LMP2_equations(double energy_tolerance, double rms_tolerance)
{
  Timer tim("LMP2");

  if (extrap_T_) extrap_T_->reinitialize();

  sma2::Array<4> Res_computed(occ_act_,occ_act_,vir_,vir_);
  allocate_R(Res_computed,domainmap_,nfzc_,my_occ_pairs_);
  Res_computed.zero();
                              
  sma2::Array<4> delta_T_computed(occ_act_,occ_act_,vir_,vir_);
  delta_T_computed("i","j","r","s") |= T_local_("i","j","r","s");
  
  ExEnv::out0() << indent << "Begin iterations for LMP2"
                << std::endl;
  ExEnv::out0() << indent
                << scprintf("%4s %16s %16s %13s %13s",
                            "iter", "E_corr   ", " delta(E_corr)",
                            "rms(delta(T))", "rms(R)")
                << std::endl;
  double ecorr_lmp2 = 0.0;
  bool converged = false;
  for (int iter=0; iter<max_iter_; iter++) {

      compute_LMP2_residual(Res_computed);
            
      double RdotR = Res_computed("i","j","r","s") * Res_computed("i","j","r","s");
      msg_->sum(RdotR);
      double Res_computed_rms
          = (T_n_element_>0?sqrt(RdotR/T_n_element_):0.0);

      tim.enter("compute deltas");
      compute_delta_T(Res_computed, delta_T_computed);
      tim.exit("compute deltas");

      T_local_("i","j","r","s") += delta_T_computed("i","j","r","s");

      if (extrap_T_) {
          Timer extim("extrap");
          Ref<sma2::Array4SCExtrapError> err_T
              = new sma2::Array4SCExtrapError(delta_T_computed,true,msg_);
          Ref<sma2::Array4SCExtrapData> dat_T
              = new sma2::Array4SCExtrapData(T_local_,true,msg_);
          extrap_T_->extrapolate(dat_T.pointer(), err_T.pointer());
          dat_T->update(T_local_);
        }
      
      // Compute new energy and check for convergence
      double ecorr_old = ecorr_lmp2;
      ecorr_lmp2 = compute_ecorr_lmp2();
      double delta_ecorr_lmp2 = ecorr_lmp2 - ecorr_old;
            
      double TdotT = delta_T_computed("i","j","r","s") * delta_T_computed("i","j","r","s");
      msg_->sum(TdotT);
      double delta_T_computed_rms
          = (T_n_element_>0?sqrt(TdotT/T_n_element_)
             :0.0);

      ExEnv::out0() << indent
                    << scprintf("%4d %16.10f %16.10f %13.10f %13.10f",
                                iter+1, ecorr_lmp2, delta_ecorr_lmp2,
                                delta_T_computed_rms,
                                Res_computed_rms)
                    << std::endl;
      if (fabs(delta_ecorr_lmp2) < energy_tolerance
          && delta_T_computed_rms < rms_tolerance) {
          converged = true;
          break;
        }
    }

  if (converged) {
      ExEnv::out0() << indent << "Converged" << std::endl;
    }
  else {
      ExEnv::out0() << indent << "Reached maximum iteration count"
                    << std::endl;
    }

  return ecorr_lmp2;
}

void
LMP2::analyze_occ_orbs(const RefSCMatrix& scf_local) {

  assert(analyze_occ_orbs_);

  const int nocc = nelectron() / 2;
  const int nocc_act = nocc - nfzc_;
  const int nao = basis()->nbasis();

  // compute dipole moment integrals
  // compute center of mass of each active occ orbital
  {
    Ref<GaussianBasisSet> bas = ref_->basis();
    Ref<PetiteList> pl = ref_->integral()->petite_list();
    Ref<OneBodyInt> m1_ints = ref_->integral()->dipole(0);

    // form skeleton mu_i in AO basis
    RefSymmSCMatrix mu_x_ao(bas->basisdim(), bas->matrixkit()); mu_x_ao.assign(0.0);
    RefSymmSCMatrix mu_y_ao(bas->basisdim(), bas->matrixkit()); mu_y_ao.assign(0.0);
    RefSymmSCMatrix mu_z_ao(bas->basisdim(), bas->matrixkit()); mu_z_ao.assign(0.0);

    const int nshell = bas->nshell();
    for(int sh1=0; sh1<nshell; sh1++) {
      int bf1_offset = bas->shell_to_function(sh1);
      int nbf1 = bas->shell(sh1).nfunction();

      int sh2max = sh1;
      for(int sh2=0; sh2<=sh2max; sh2++) {
        int bf2_offset = bas->shell_to_function(sh2);
        int nbf2 = bas->shell(sh2).nfunction();

        m1_ints->compute_shell(sh1,sh2);
        const double *m1intsptr = m1_ints->buffer();

        int bf1_index = bf1_offset;
        for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, m1intsptr+=3*nbf2) {
          int bf2_index = bf2_offset;
          const double *ptr1 = m1intsptr;
          int bf2max;
          if (sh1 == sh2)
            bf2max = bf1;
          else
            bf2max = nbf2-1;
          for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++) {

            // the negative charge of the electron is not included
            mu_x_ao.set_element(bf1_index, bf2_index, ptr1[0]);
            mu_y_ao.set_element(bf1_index, bf2_index, ptr1[1]);
            mu_z_ao.set_element(bf1_index, bf2_index, ptr1[2]);
            ptr1 += 3;

          }
        }
      }
    }
    m1_ints = 0;

    const int nbasis = bas->nbasis();
    for(int bf1=0; bf1<nbasis; bf1++) {
      for(int bf2=0; bf2<=bf1; bf2++) {
        mu_x_ao(bf2,bf1) = mu_x_ao(bf1,bf2);
        mu_y_ao(bf2,bf1) = mu_y_ao(bf1,bf2);
        mu_z_ao(bf2,bf1) = mu_z_ao(bf1,bf2);
      }
    }

    ExEnv::out0() << indent << "center-of-mass of active occupied orbitals:" << std::endl;
    ExEnv::out0() << indent << "  i              X                      Y                      Z          " << std::endl;
    ExEnv::out0() << indent << "===== ====================== ====================== ======================" << std::endl;
    assert(nao == mu_x_ao.n());
    std::vector<double> mu_x_ri(nao);
    std::vector<double> mu_y_ri(nao);
    std::vector<double> mu_z_ri(nao);
    for(int i=0; i<nocc_act; ++i) {
      for(int r=0; r<nao; ++r) {
        mu_x_ri[r] = 0.0;
        mu_y_ri[r] = 0.0;
        mu_z_ri[r] = 0.0;
        for(int c=0; c<nao; ++c) {
          mu_x_ri[r] += scf_local(c,i) * mu_x_ao(r,c);
          mu_y_ri[r] += scf_local(c,i) * mu_y_ao(r,c);
          mu_z_ri[r] += scf_local(c,i) * mu_z_ao(r,c);
        }
      }
      SCVector3 R(0.0, 0.0, 0.0);
      for(int r=0; r<nao; ++r) R(0) += scf_local(r,i) * mu_x_ri[r];
      for(int r=0; r<nao; ++r) R(1) += scf_local(r,i) * mu_y_ri[r];
      for(int r=0; r<nao; ++r) R(2) += scf_local(r,i) * mu_z_ri[r];
      r_i_[i] = R;

      ExEnv::out0() << indent << scprintf("%5d %22.15lf %22.15lf %22.15lf", i,
                                          R(0), R(1), R(2))
                    << std::endl;
    }
  }

  // print out all the pair data
  ExEnv::out0() << indent << "pair data" << std::endl;
  ExEnv::out0() << indent << "  i     j             |r|                   e(MP2)              |r|_domain      " << std::endl;
  ExEnv::out0() << indent << "===== ===== ====================== ====================== ======================" << std::endl;
  for(int i=0; i<nocc_act; ++i) {
    for(int j=0; j<=i; ++j) {
      const int ij = i*nocc_act + j;
      const double r_ij = (r_i_[i] - r_i_[j]).norm();
      ExEnv::out0() << indent << scprintf("%5d %5d %22.15lf %22.15lf %22.15lf", i, j,
                                          r_ij, emp2_ij_[ij], dist_ij_[ij]
                                         )
                          << std::endl;
    }
  }

}

}
