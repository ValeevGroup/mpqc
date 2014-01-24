//
// ccr12_info.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: TS & EFV
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

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <math/scmat/abstract.h>
#include <chemistry/qc/ccr12/ccr12_info.h>

using namespace sc;
using namespace std;

static ClassDesc CCR12_Info_cd(
  typeid(CCR12_Info),"CCR12_Info",1,"virtual public RefCount",
  0,0,0);

CCR12_Info::CCR12_Info(const Ref<R12WavefunctionWorld>& r12world,
                       const Ref<MemoryGrp>& mem, size_t memorysize,
                       const Ref<SCF> reference, int nfc, int nfv, int nirr,
                       long workmem, long memsize, int nnode, int ndiis,
                       string theory, string per, int tilef):
r12world_(r12world), mem_(mem), ref_(reference), nfzc_(nfc),
nfzv_(nfv), nirrep_(nirr), workmemsize_(workmem), theory_(theory), perturbative_(per), maxtilesize_(tilef)
{
  MPQC_ASSERT(r12world_->r12tech()->corrfactor()->nfunctions() == 1);
  restricted_ = !ref()->spin_polarized();

  needs();

  set_naocc(r12world_->refwfn()->occ_act_sb(Alpha)->rank(),
            r12world_->refwfn()->occ_act_sb(Beta )->rank());

//MPQC_ASSERT(r12world_->r12tech()->ansatz()->amplitudes() != R12Technology::GeminalAmplitudeAnsatz_fullopt);
  MPQC_ASSERT(r12world_->r12tech()->ansatz()->amplitudes() != R12Technology::GeminalAmplitudeAnsatz_scaledfixed);
  set_fixed(r12world_->r12tech()->ansatz()->amplitudes() == R12Technology::GeminalAmplitudeAnsatz_fixed);
  set_navir(r12world_->refwfn()->orbs(Alpha)->rank() - naoa() - nfzc() - nfzv(),
            r12world_->refwfn()->orbs(Beta )->rank() - naob() - nfzc() - nfzv());
  set_mosym(r12world_->refwfn()->orbs(Alpha)->orbsym(),
            r12world_->refwfn()->orbs(Beta )->orbsym());

  if (need_w1_) {
    set_ncabs(r12world_->cabs_space(Alpha)->rank(),
              r12world_->cabs_space(Beta )->rank());
    set_cabssym(r12world_->cabs_space(Alpha)->orbsym(),
                r12world_->cabs_space(Beta )->orbsym());
  } else {
    set_ncabs(0L, 0L);
  }

  determine_tilesizes();
  maxtilesize_ = *max_element(range_.begin(), range_.end());
  if (mem_->me() == 0) print_tile_info();

  irrep_e_ = 0L;
  irrep_f_ = 0L;
  irrep_v_ = 0L;
  irrep_t_ = 0L;
  irrep_y_ = 0L;

  // More to come... e.g. dipole,...
  bool do_lambda
    = (perturbative_ == "(2)T" || perturbative_ == "(2)TQ" || perturbative_ == "(2)R12FULL" || perturbative_ == "(2)TQR12");

  // compute the correlated spin-orbital space used by SMITH
  compute_corr_space();

  // now get MemoryGrp ready
  mem_->set_localsize(memorysize);
  // compute all source integrals now, before MemoryGrp is used by Tensors
  if (restricted_)
    compute_source_integrals_rhf();
  else
    compute_source_integrals_uhf();

  // Retrieve B and X intermediate
  if (need_w1()) {
    retrieve_B_and_X_ii();
  }

  // retrieve and sort eigenvalues of Fock matrix
  orbital_energies();

  long static_size = 0L;

  d_f1 = new Tensor("f1", mem_);  // for later convenience (in disk based algorithm),
                                  // we are naming all the tensors uniquely.
  offset_f1();
  static_size += d_f1->get_filesize() * sizeof(double);
  fill_in_f1();
  mem_->sync();

  // now need_w1 and need_w2 are treated within a single function
  d_v2 = new Tensor("v2", mem_);
  offset_v2();
  static_size += d_v2->get_filesize() * sizeof(double);
  mem_->sync();
  fill_in_v2();

  if (need_t1()) {
    d_t1 = new Tensor("t1",mem_);
    offset_t1(d_t1, mem_->me() == 0);
    static_size += d_t1->get_filesize() * sizeof(double) * (ndiis + 1) * 2;
    mem_->sync();
    if (do_lambda) {
      d_lambda1 = new Tensor("lambda1", mem_);
      offset_l1(d_lambda1);
      static_size += d_lambda1->get_filesize() * sizeof(double);
      mem_->sync();
    }
  }

  if (need_gt2()){
    MPQC_ASSERT(need_w1());

    d_gt2 = new Tensor("gt2", mem_);
    offset_gt2(d_gt2, true);
    static_size += d_gt2->get_filesize() * sizeof(double) * (ndiis + 1) * 2;
    mem_->sync();
    // needs to be initialized here (with fixed amplitudes or full-opt MP2-R12 amplitudes)...
    if (do_lambda || r12world_->r12tech()->ansatz()->amplitudes() != R12Technology::GeminalAmplitudeAnsatz_fullopt) {
      d_glambda2 = new Tensor("glambda2", mem_);
      offset_gt2(d_glambda2, false);
      static_size += d_glambda2->get_filesize() * sizeof(double);
      mem_->sync();
    }
  }

  if (need_t2()){
    d_t2 = new Tensor("t2", mem_);
    offset_t2(d_t2, mem_->me() == 0);
    static_size += d_t2->get_filesize() * sizeof(double) * (ndiis + 1) * 2;
    mem_->sync();

    if (do_lambda) {
      d_lambda2 = new Tensor("lambda2", mem_);
      offset_l2(d_lambda2);
      static_size += d_lambda2->get_filesize() * sizeof(double);
      mem_->sync();
    }
  }

  if (need_t3()){
    d_t3 = new Tensor("t3", mem_);
    offset_t3(d_t3, mem_->me() == 0);
    static_size += d_t3->get_filesize() * sizeof(double) * (ndiis + 1) * 2;
    mem_->sync();
  }

  if (need_t4()){
    d_t4 = new Tensor("t4", mem_);
    offset_t4(d_t4, mem_->me() == 0);
    static_size += d_t4->get_filesize() * sizeof(double) * (ndiis + 1) * 2;
    mem_->sync();
  }


  if (need_gt2() || need_w1()){
    // need_gt2, then R12 tensors are initialized.
    // (2)R12 methods do not need gt2 (hence need_gt2_ = false), but they do need w1 tensors.

    // V intermediate (and its conjugate)
    // note that V^pA_ii type integrals are not needed in CCSD(R12) calculations.
    d_vr2 = new Tensor("vr2", mem_);
    offset_vr2();
    static_size += d_vr2->get_filesize() * sizeof(double);
    mem_->sync();

    d_vd2 = new Tensor("vd2", mem_);
    offset_vd2();
    static_size += d_vd2->get_filesize() * sizeof(double);
    mem_->sync();

    fill_in_vr_and_vd();

    if (need_F_) {
      d_fr2 = new Tensor("fr2", mem_);
      offset_fr2();
      static_size += d_fr2->get_filesize() * sizeof(double);
      mem_->sync();

      d_fd2 = new Tensor("fd2", mem_);
      offset_fd2();
      static_size += d_fd2->get_filesize() * sizeof(double);
      mem_->sync();

      fill_in_fr_and_fd();
    }

    // X intermediate
    d_xs2 = new Tensor("xs2", mem_);
    offset_gt2(d_xs2, true);
    static_size += d_xs2->get_filesize() * sizeof(double);
    mem_->sync();

    // B intermediate
    d_bs2 = new Tensor("bs2", mem_);
    offset_gt2(d_bs2, true);
    static_size += d_bs2->get_filesize() * sizeof(double);
    mem_->sync();

    // F * gt2 (\tilde{t} in our papers)
    // filled at run-time
    // those are not needed in (2)R12 type theories
    if (need_gt2() && need_F_) {
      d_qy = new Tensor("qy", mem_);
      offset_qy();
      static_size += d_qy->get_filesize() * sizeof(double);
      mem_->sync();

      if (need_w2()) {
        d_qx = new Tensor("qx", mem_);
        offset_qx();
        static_size += d_qx->get_filesize() * sizeof(double);
        mem_->sync();
      }

      if (do_lambda || r12world_->r12tech()->ansatz()->amplitudes() != R12Technology::GeminalAmplitudeAnsatz_fullopt) {
        // gl2 * F^dagger (\tilde{l} in our papers)
        // filled at run-time
        d_ly = new Tensor("ly", mem_);
        offset_ly();
        static_size += d_ly->get_filesize() * sizeof(double);
        mem_->sync();

        if (need_w2()) {
          d_lx = new Tensor("lx", mem_);
          offset_lx();
          static_size += d_lx->get_filesize() * sizeof(double);
          mem_->sync();
        }
      }
    }

    if (need_w2()){ // only for full CC-R12 methods
      // P intermediate
      d_ps2 = new Tensor("ps2", mem_);
      offset_gt2(d_ps2, true);
      static_size += d_ps2->get_filesize() * sizeof(double);
      mem_->sync();
    }

    fill_in_iiii();
  }

  if (perturbative_ == "(T)R12[DT]") {
    if (r12world_->r12tech()->ansatz()->orbital_product_GG() != R12Technology::OrbProdGG_pq) {
      throw std::runtime_error("(T)R12[DT] requires GG space to be pq.");
    }
    // generalized V intermediate.
    d_vd2_gen = new Tensor("vd2_gen", mem_);
    // Setting need_CABS = false, need_xx = true.
    offset_vd2_gen(false, false);
    fill_in_vd2_gen(false, false);
    // Creates B and X intermediates of ip space.
    retrieve_B_and_X_ip();
  }

  // prediagonalization set-up for certain class of methods.
  // EFV: 06/19/2011 seem to need this even in CCSD(R12)?
//  if (need_w1() && theory_ != "CCSD(R12)") {
  if (need_w1())
    prediagon(bdiag_, lmatrix_);
//  }

  // Making initial guess for t2.
  if (!need_gt2() || (need_gt2() && (perturbative_ == "(T)R12" || perturbative_ == "(2)R12"))) {
    guess_t2(d_t2);
  } else if (r12world_->r12tech()->ansatz()->amplitudes() != R12Technology::GeminalAmplitudeAnsatz_fullopt) {
    // assuming d_gt2 is already filled in.
    guess_t2_r12(d_t2, d_gt2);
  } else {
    // now using MP1-R12/SP amplitude even for full-opt calculations
    guess_t2_r12(d_t2, d_gt2);
  }


/// perhaps we need somehow to obtain 6-index B for Jacobi iterations.
/// The tensors would better be inverted before they are transformed
/// into block-wise structure.

#ifndef DISK_BASED_SMITH
  const long input_tensors = static_size / 1000000L / nnode;
#else
  const long input_tensors = 0L;
#endif
  const long work_space    = workmem / 1000000L;
  const long intermediates = memsize / 1000000L - work_space - input_tensors;

  if(intermediates < 0){
      throw ProgrammingError("CCR12_Info::CCR12_Info --- needs more memory", __FILE__, __LINE__);
  }
  if (need_w1()) ExEnv::out0() << endl;
  ExEnv::out0() << indent << "input tensors (total)    : " << setw(6) << input_tensors*nnode << " MB" << endl;
  ExEnv::out0() << indent << "input tensors (per node) : " << setw(6) << input_tensors       << " MB" << endl << endl;
  ExEnv::out0() << indent << "work space    (per node) : " << setw(6) << work_space          << " MB" << endl << endl;
  ExEnv::out0() << indent << "intermediates (total)    : " << setw(6) << intermediates*nnode << " MB" << endl;
  ExEnv::out0() << indent << "intermediates (per node) : " << setw(6) << intermediates       << " MB" << endl;

  mem_->sync();

}

CCR12_Info::~CCR12_Info(){
  mem_->set_localsize(0);
}


void CCR12_Info::determine_tilesizes(){

  double workdouble = static_cast<double>(workmemsize_) / sizeof(double);
  determine_maxtilesize(workdouble);

  ExEnv::out0() << indent << "max tile size: " << maxtilesize_ << endl;

  int offset = 0;
  int alpha  = 0;
  // alpha hole
  noab_ = determine_tilesize_each(nfzc(), nfzc() + naoa(), 1, mosyma_, true, offset, alpha);
  // beta  hole
  if (restricted_) {
    noab_ += determine_tilesize_each(nfzc(), nfzc() + naob(), 2, mosymb_, true, offset, alpha);
  } else {
    alpha = range_.size();
    noab_ += determine_tilesize_each(nfzc(), nfzc() + naob(), 2, mosymb_, true, offset, alpha);
  }
  // alpha particle
  alpha = range_.size();
  nvab_ = determine_tilesize_each(nfzc() + naoa(), nfzc() + naoa() + nava(), 1, mosyma_, true, offset, alpha);
  // beta  particle
  if (restricted_) {
    nvab_ += determine_tilesize_each(nfzc() + naob(), nfzc() + naob() + navb(), 2, mosymb_, true, offset, alpha);
  } else {
    alpha = range_.size();
    nvab_ += determine_tilesize_each(nfzc() + naob(), nfzc() + naob() + navb(), 2, mosymb_, true, offset, alpha);
  }
  // alpha CABS
  alpha = range_.size();
  ncab_ = determine_tilesize_each(0, nria(), 1, cabssyma_, false, offset, alpha);
  // beta  CABS
  if (restricted_) {
    ncab_ += determine_tilesize_each(0, nrib(), 2, cabssymb_, false, offset, alpha);
  } else {
    alpha = range_.size();
    ncab_ += determine_tilesize_each(0, nrib(), 2, cabssymb_, false, offset, alpha);
  }

}


int CCR12_Info::determine_tilesize_each(int istart, int iend, int spin,
           vector<unsigned int> &orbsym, bool mo, int &offset, int alpha){

// Integers in smith-generated codes must be "long int."
  int blocks = 0;
  for (int irrep = 0; irrep < nirrep_int(); ++irrep){
    int num_orb = 0;
    for (int orb = istart; orb < iend; ++orb){
      if (static_cast<int>(orbsym[orb]) == irrep){
        orbital_sym_sorted_.push_back(static_cast<long>(irrep));
        orbital_spin_sorted_.push_back(static_cast<long>(spin));
        if (mo) {
//        orbital_evl_sorted_.push_back((double)orb); // eivenvalue needed here
          momap_.push_back(static_cast<long>(orb - nfzc()));
        }
        num_orb++;
      }
    }
    int acc_range = 0;
    int nblock = num_orb >0 ? (num_orb - 1) / maxtilesize_ + 1 : 0;
    for (int ib = 0; ib < nblock; ++ib){
      sym_.push_back(   static_cast<long>(irrep));
      range_.push_back( static_cast<long>((num_orb * (ib + 1)) / nblock - acc_range));
      spin_.push_back(  static_cast<long>(spin));
      offset_.push_back(static_cast<long>(offset));
      alpha_.push_back( static_cast<long>(alpha));
      acc_range += range_[range_.size() - 1];
      offset    += range_[range_.size() - 1];
      alpha++;
    }
    blocks += nblock;
  }
  return blocks;
}

void CCR12_Info::determine_maxtilesize(double memory){

  if (maxtilesize_ > 0) {
    ExEnv::out0() << indent << ":: caution :: max tile size is manually set to " << maxtilesize_ << endl; 
  } else {

    string t = theory_;
    // transform(t.begin(),t.end(),t.begin(),(int (*)(int))tolower);
    const double max_blocks = 7.0;

    if (t == "CCSD" || t == "CCSD(R12)" || t == "CCSD-R12"){
      maxtilesize_=static_cast<int>(::pow(memory / max_blocks, 0.25));
    } else if (t == "CCSDT" || t == "CCSDT(R12)" || t == "CCSDT-R12") {
      maxtilesize_ = static_cast<int>(::pow(memory / max_blocks, 1.0 / 6.0));
    } else if (t == "CCSDTQ" || t == "CCSDTQ(R12)" || t == "CCSDTQ-R12") {
      maxtilesize_ = static_cast<int>(::pow(memory / max_blocks, 0.125));
    } else {
      ExEnv::out0() << indent << "theory: " << theory_ << endl;
      throw ProgrammingError("CCR12_Info::tilesize -- not yet implemented", __FILE__, __LINE__);
    }

    // Perturbative correction needs an additional memory area.
    // There is a room to let tilesize larger than this by calculating explicitly the memory demands...
    if (perturbative_ == "(T)" || perturbative_ == "(T)R12[DT]" || perturbative_ == "(T)R12"){
      const int p_maxtilesize = static_cast<int>(::pow(memory / 5.0, 1.0 / 6.0));
      if (p_maxtilesize < maxtilesize_) maxtilesize_ = p_maxtilesize;
    } else if (perturbative_ == "(2)T") {
      const int p_maxtilesize = static_cast<int>(::pow(memory / 5.0, 1.0 / 6.0));
      if (p_maxtilesize < maxtilesize_) maxtilesize_ = p_maxtilesize;
    } else if (perturbative_ == "(2)TQ" || perturbative_ == "(2)TQR12") {
      const int p_maxtilesize = static_cast<int>(::pow(memory / 5.0, 1.0 / 8.0));
      if (p_maxtilesize < maxtilesize_) maxtilesize_ = p_maxtilesize;
    } else if (perturbative_ == "(2)R12" || perturbative_ == "(2)R12FULL") {
      /// nothing happens
    }
  }

}

void CCR12_Info::print_tile_info(){
  string bar(40, '-');
  ExEnv::out0() << endl << indent << setw(33) << ">>>> tile infomation <<<<" << endl << endl;
  ExEnv::out0() << indent << bar     << endl;
  ExEnv::out0() << indent << setw(8) << "spin"    << setw(8) << "irrep"
                          << setw(8) << "length"  << setw(8) << "offset"
                          << setw(8) << "alpha"   << endl;
  for (int i = 0; i < range_.size(); ++i)
        ExEnv::out0() << indent << setw(8) << spin_[i]   << setw(8) << sym_[i]
                                << setw(8) << range_[i]  << setw(8) << offset_[i]
                                << setw(8) << alpha_[i]  << endl;
  ExEnv::out0() << indent <<  bar    << endl << endl;
}

void CCR12_Info::needs(){
  if (theory_ == "CCSD") {
    need_w1_ = false; need_w2_ = false; need_t1_ = true; need_t2_ = true;
    need_gt2_ = false; need_t3_ = false; need_t4_ = false; need_FAA_ = false;
    need_VpA_ = false; need_F_ = false;
  } else if(theory_ == "CCSD-R12") {
    need_w1_ = true; need_w2_ = true; need_t1_ = true; need_t2_ = true;
    need_gt2_ = true; need_t3_ = false; need_t4_ = false; need_FAA_ = true;
    need_VpA_ = true; need_F_ = true;
  } else if(theory_ == "CCSD(R12)") {
    need_w1_ = true; need_w2_ = false; need_t1_ = true; need_t2_ = true;
    need_gt2_ = true; need_t3_ = false; need_t4_ = false; need_FAA_ = false;
    need_VpA_ = false; need_F_ = true;
  } else if(theory_ == "CCSDT"){
    need_w1_ = false; need_w2_ = false; need_t1_ = true; need_t2_ = true;
    need_gt2_ = false; need_t3_ = true; need_t4_ = false; need_FAA_ = false;
    need_VpA_ = false; need_F_ = false;
  } else if(theory_ == "CCSDTQ"){
    need_w1_ = false; need_w2_ = false; need_t1_ = true; need_t2_ = true;
    need_gt2_ = false; need_t3_ = true; need_t4_ = true; need_FAA_ = false;
    need_VpA_ = false; need_F_ = false;
  } else {
    throw ProgrammingError("CCR12_Info::needs", __FILE__, __LINE__);
  }

  // unscreened version
  if (perturbative_ == "(2)R12FULL" || perturbative_ == "(2)TQR12") {
    need_w1_ = true;
    need_F_ = true;
  }
  // screened version
  if (perturbative_ == "(T)R12[DT]" || perturbative_ == "(T)R12" || perturbative_ == "(2)R12") {
    need_w1_ = true;
    if (this->r12world()->r12tech()->ansatz()->diag()) {
      need_gt2_ = true;
    }
  }

}


double CCR12_Info::get_e(const Ref<Tensor>& d_e_){
  double out;
  d_e_->get_block(0L, &out);
  return out;
}


void CCR12_Info::print(ostream& o){
  o << indent << "theory = " << theory_ << endl;
  o << indent << "perturbative = " << perturbative_ << endl;
  o << indent << "reference wave function = " << endl;
  o << incindent;
  ref_->print(o);
  o << decindent;
}

void CCR12_Info::orbital_energies(){

  Ref<OrbitalSpace> aobs_space = r12world()->refwfn()->orbs_sb(Alpha);
  Ref<OrbitalSpace> bobs_space = r12world()->refwfn()->orbs_sb(Beta);

  // compute map from indices in full spin-orbital space to indices in the respective spin spaces
  vector<long> amap;
  {
    vector<int> smith_to_obs = sc::map(*aobs_space, *corr_space_, false);
    amap.resize(smith_to_obs.size());
    std::copy(smith_to_obs.begin(), smith_to_obs.end(), amap.begin());
  }
  vector<long> bmap;
  {
    vector<int> smith_to_obs = map(*bobs_space, *corr_space_, false);
    bmap.resize(smith_to_obs.size());
    std::copy(smith_to_obs.begin(), smith_to_obs.end(), bmap.begin());
  }

  orbital_evl_sorted_.resize(naoa() + naob() + nava() + navb());
  int count = 0;
  for (long g1b = 0L; g1b < nab(); ++g1b) {
    if (g1b < noab() + nvab()) {
      const long size1 = get_range(g1b);
      const int spin1 = get_spin(g1b);
      const int offset1 = (spin1 == 1) ? amap[get_offset(g1b)] : bmap[get_offset(g1b)];

      for(int i1=0; i1<size1; ++i1, ++count) {
        const int ii1 = i1 + offset1;
        if (spin1 == 1)  // Alpha
          orbital_evl_sorted_[count] = F_[Alpha].get_element(ii1,ii1);
        else
          orbital_evl_sorted_[count] = F_[Beta].get_element(ii1,ii1);

      }
    }
  }

#if 0
  ExEnv::out0() << "orbital_evl_sorted_ = " << endl;
  for(std::vector<double>::iterator i=orbital_evl_sorted_.begin(); i!=orbital_evl_sorted_.end(); ++i)
    ExEnv::out0() << *i << endl;
#endif
}


void CCR12_Info::retrieve_B_and_X_ii() {

  Ref<OrbitalSpace> Gspace = r12eval()->GGspace(Alpha);
  Ref<OrbitalSpace> ispace = r12eval()->occ_act(Alpha);

  RefSymmSCMatrix b_tmp = r12int_eval_->B(AlphaBeta);
  RefSymmSCMatrix x_tmp = r12int_eval_->X(AlphaBeta);

  const unsigned int ni = ispace->rank();
  const unsigned int nG = Gspace->rank();
  const int target_pairsize = ni * ni;
  RefSCDimension dim_occpair(new SCDimension(target_pairsize));
  Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
  B_ = kit->symmmatrix(dim_occpair);
  X_ = kit->symmmatrix(dim_occpair);

  // map ispace to Gspace
  vector<unsigned int> imap = (*Gspace << *ispace);
  unsigned int r01 = 0;
  for(unsigned int r0=0; r0<ni; ++r0) {
    const unsigned int rr0 = imap[r0];
    for(unsigned int r1=0; r1<ni; ++r1, ++r01) {
      const unsigned int rr1 = imap[r1];
      const unsigned int rr01 = rr0 * nG + rr1;

      unsigned int c01 = 0;
      for(unsigned int c0=0; c0<ni; ++c0) {
        const unsigned int cc0 = imap[c0];
        for(unsigned int c1=0; c1<ni; ++c1, ++c01) {
          const unsigned int cc1 = imap[c1];
          const unsigned int cc01 = cc0 * nG + cc1;

          const double xval = x_tmp(rr01,cc01);
          const double bval = b_tmp(rr01,cc01);
          X_(r01,c01) = xval;
          B_(r01,c01) = bval;
        }
      }
    }
  }

}


void CCR12_Info::retrieve_B_and_X_ip() {

  Ref<OrbitalSpace> Gspace = r12eval()->GGspace(Alpha);
  Ref<OrbitalSpace> aspace = r12eval()->vir_act(Alpha);
  Ref<OrbitalSpace> ispace = r12eval()->occ_act(Alpha);

  // extract ip subblock
  RefSymmSCMatrix b_tmp = r12int_eval_->B(AlphaBeta);
  RefSymmSCMatrix x_tmp = r12int_eval_->X(AlphaBeta);

  const unsigned int ni = ispace->rank();
  const unsigned int na = aspace->rank();
  const unsigned int nG = Gspace->rank();
  const int target_pairsize = ni * ni + na * ni * 2;
  RefSCDimension dim_pair(new SCDimension(target_pairsize));
  Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
  B_ip_ = kit->symmmatrix(dim_pair);
  X_ip_ = kit->symmmatrix(dim_pair);

  // map ispace and pspace to Gspace
  vector<unsigned int> imap = (*Gspace << *ispace);
  vector<unsigned int> amap = (*Gspace << *aspace);

  typedef pair<unsigned int, unsigned int> PairUint;
  vector<PairUint> mapping(target_pairsize);
  vector<PairUint>::iterator iter = mapping.begin();

  unsigned int r01 = 0;
  for(unsigned int r0=0; r0<ni; ++r0) {
    const unsigned int rr0 = imap[r0];
    for(unsigned int r1=0; r1<ni; ++r1, ++r01, ++iter)
      *iter = make_pair(r01, rr0*nG+imap[r1]);
  }
  for(unsigned int r0=0; r0<ni; ++r0) {
    const unsigned int rr0 = imap[r0];
    for(unsigned int r1=0; r1<na; ++r1, ++r01, ++iter)
      *iter = make_pair(r01, rr0*nG+amap[r1]);
  }
  for(unsigned int r0=0; r0<na; ++r0) {
    const unsigned int rr0 = amap[r0];
    for(unsigned int r1=0; r1<ni; ++r1, ++r01, ++iter)
      *iter = make_pair(r01, rr0*nG+imap[r1]);
  }

  for (vector<PairUint>::const_iterator riter = mapping.begin(); riter != mapping.end(); ++riter) {
    const unsigned int r01 = riter->first;
    const unsigned int rr01 = riter->second;
    for (vector<PairUint>::const_iterator citer = mapping.begin(); citer != mapping.end(); ++citer) {
      const unsigned int c01 = citer->first;
      const unsigned int cc01 = citer->second;
      const double xval = x_tmp(rr01,cc01);
      const double bval = b_tmp(rr01,cc01);
      X_ip_(r01,c01) = xval;
      B_ip_(r01,c01) = bval;
    }
  }

}


