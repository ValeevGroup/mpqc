//
// cadfclhf.cc
//
// Copyright (C) 2013 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Nov 13, 2013
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

// Standard library includes
#include <memory>
#include <math.h>

// Boost includes
#include <boost/tuple/tuple_io.hpp>

// MPQC includes
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <util/group/messmpi.h>
#include <util/misc/regtime.h>
#include <util/misc/scexception.h>
#include <util/misc/xmlwriter.h>
#include <math/scmat/blas.h>
#include <util/container/conc_cache.h>

#include "cadfclhf.h"

#define EIGEN_NO_AUTOMATIC_RESIZING 1

#define CADF_USE_BLAS 0


using namespace sc;
using namespace std;
using boost::thread;
using boost::thread_group;
using namespace sc::parameter;
static constexpr bool xml_debug_ = false;

typedef std::pair<int, int> IntPair;
typedef CADFCLHF::CoefContainer CoefContainer;
typedef CADFCLHF::Decomposition Decomposition;
typedef std::pair<CoefContainer, CoefContainer> CoefPair;

static boost::mutex debug_print_mutex;

////////////////////////////////////////////////////////////////////////////////
// Debugging asserts and outputs


////////////////////////////////////////////////////////////////////////////////

ClassDesc CADFCLHF::cd_(
    typeid(CADFCLHF),
    "CADFCLHF",
    1, // verion number
    "public CLHF",
    0, // default constructor pointer
    create<CADFCLHF>, // KeyVal constructor pointer
    create<CADFCLHF> // StateIn constructor pointer
);

//////////////////////////////////////////////////////////////////////////////////

CADFCLHF::CADFCLHF(StateIn& s) :
    SavableState(s),
    CLHF(s),
    stats_(),
    iter_stats_(0)
{
  //----------------------------------------------------------------------------//
  // Nothing to do yet
  throw FeatureNotImplemented("SavableState construction of CADFCLHF",
      __FILE__, __LINE__, class_desc());
  //----------------------------------------------------------------------------//
}

//////////////////////////////////////////////////////////////////////////////////

CADFCLHF::CADFCLHF(const Ref<KeyVal>& keyval) :
    CLHF(keyval),
#if USE_INTEGRAL_CACHE
    // TODO Initialize the number of bins to a more reasonable size, check for enough memory, etc
    ints3_(make_shared<ThreeCenterIntCache>(
        gbs_->nshell() * gbs_->nshell() * gbs_->nshell() / scf_grp_->n()
    )),
    // TODO Initialize the number of bins to a more reasonable size, check for enough memory, etc
    ints2_(make_shared<TwoCenterIntCache>(
        gbs_->nshell() * gbs_->nshell() / scf_grp_->n()
    )),
#endif
    decomps_(make_shared<DecompositionCache>(
        molecule()->natom() * molecule()->natom()
    )),
    ints4maxes_(make_shared<FourCenterMaxIntCache>(
        gbs_->nshell() * gbs_->nshell()
    )),
    local_pairs_spot_(0)
{
  //----------------------------------------------------------------------------//
  // Get the auxiliary basis set
  dfbs_ << keyval->describedclassvalue("df_basis", KeyValValueRefDescribedClass(0));
  if(dfbs_.null()){
    throw InputError("CADFCLHF requires a density fitting basis set",
        __FILE__, __LINE__, "df_basis");
  }
  //----------------------------------------------------------------------------//
  // get the bound for filtering pairs.  This should be smaller than the general
  //   Schwarz bound.
  // TODO Figure out a reasonable default value for this
  pair_screening_thresh_ = keyval->doublevalue("pair_screening_thresh", KeyValValuedouble(1e-8));
  density_screening_thresh_ = keyval->doublevalue("density_screening_thresh", KeyValValuedouble(1e-8));
  full_screening_thresh_ = keyval->doublevalue("full_screening_thresh", KeyValValuedouble(1e-8));
  distance_screening_thresh_ = keyval->doublevalue("distance_screening_thresh", KeyValValuedouble(full_screening_thresh_));
  coef_screening_thresh_ = keyval->doublevalue("coef_screening_thresh", KeyValValuedouble(1e-8));
  full_screening_expon_ = keyval->doublevalue("full_screening_expon", KeyValValuedouble(1.0));
  distance_damping_factor_ = keyval->doublevalue("distance_damping_factor", KeyValValuedouble(1.0));
  //----------------------------------------------------------------------------//
  // For now, use coulomb metric.  We can easily make this a keyword later
  metric_oper_type_ = TwoBodyOper::eri;
  //----------------------------------------------------------------------------//
  do_linK_ = keyval->booleanvalue("do_linK", KeyValValueboolean(false));
  linK_block_rho_ = keyval->booleanvalue("linK_block_rho", KeyValValueboolean(false));
  linK_sorted_B_contraction_ = keyval->booleanvalue("linK_sorted_B_contraction", KeyValValueboolean(false));
  linK_use_distance_ = keyval->booleanvalue("linK_use_distance", KeyValValueboolean(false));
  //----------------------------------------------------------------------------//
  print_screening_stats_ = keyval->intvalue("print_screening_stats", KeyValValueint(0));
  //----------------------------------------------------------------------------//
  print_iteration_timings_ = keyval->booleanvalue("print_iteration_timings", KeyValValueboolean(false));
  //----------------------------------------------------------------------------//
  use_norms_nu_ = keyval->booleanvalue("use_norms_nu", KeyValValueboolean(true));
  use_norms_sigma_ = keyval->booleanvalue("use_norms_sigma", KeyValValueboolean(true));
  old_two_body_ = keyval->booleanvalue("old_two_body", KeyValValueboolean(old_two_body_));
  xml_screening_data_ = keyval->booleanvalue("xml_screening_data", KeyValValueboolean(xml_screening_data_));
  //----------------------------------------------------------------------------//
  xml_debug_ = keyval->booleanvalue("xml_debug", KeyValValueboolean(false));
  //----------------------------------------------------------------------------//
  initialize();
  //----------------------------------------------------------------------------//
}

//////////////////////////////////////////////////////////////////////////////////

CADFCLHF::~CADFCLHF()
{
  if(have_coefficients_){
    // Clean up the coefficient data
    deallocate(coefficients_data_);
  }
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::initialize()
{
  //----------------------------------------------------------------------------//
  // Check that the density is local
  if(!local_dens_){
    throw FeatureNotImplemented("Can't handle density matrices that don't fit on one node",
        __FILE__, __LINE__, class_desc());
  }
  //----------------------------------------------------------------------------//
  // need a nonblocked cl_gmat_ in this method
  Ref<PetiteList> pl = integral()->petite_list();
  gmat_ = gbs_->so_matrixkit()->symmmatrix(pl->SO_basisdim());
  gmat_.assign(0.0);
  //----------------------------------------------------------------------------//
  // Determine if the message group is an instance of MPIMessageGrp
  using_mpi_ = dynamic_cast<MPIMessageGrp*>(scf_grp_.pointer()) ? true : false;
  //----------------------------------------------------------------------------//
  have_coefficients_ = false;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::init_threads()
{
  Timer timer("init threads");
  assert(not threads_initialized_);

  //----------------------------------------------------------------------------//

  const int nthread = threadgrp_->nthread();
  const int n_node = scf_grp_->n();

  //----------------------------------------------------------------------------//
  // initialize the two electron integral classes

  // ThreeCenter versions
  integral()->set_basis(gbs_, gbs_, dfbs_);
  size_t storage_avail = integral()->storage_unused();
  eris_3c_.resize(nthread);
  do_threaded(nthread, [&](int ithr){
    eris_3c_[ithr] = (integral()->coulomb<3>());
    eris_3c_[ithr]->set_integral_storage(storage_avail/nthread);
  });
  for (int i=0; i < nthread; i++) {
    // TODO different fitting metrics
    if(metric_oper_type_ == coulomb_oper_type_){
      metric_ints_3c_.push_back(eris_3c_[i]);
    }
    else{
      throw FeatureNotImplemented("non-coulomb metrics in CADFCLHF", __FILE__, __LINE__, class_desc());
    }
  }

  // TwoCenter versions
  integral()->set_basis(dfbs_, dfbs_);
  eris_2c_.resize(nthread);
  do_threaded(nthread, [&](int ithr){
    eris_2c_[ithr] = integral()->coulomb<2>();
  });
  for (int i=0; i < nthread; i++) {
    if(metric_oper_type_ == coulomb_oper_type_){
      metric_ints_2c_.push_back(eris_2c_[i]);
    }
    else{
      throw FeatureNotImplemented("non-coulomb metrics in CADFCLHF", __FILE__, __LINE__, class_desc());
    }
  }

  // Reset to normal setup
  integral()->set_basis(gbs_, gbs_, gbs_, gbs_);

  //----------------------------------------------------------------------------//
  // TODO fix this so that deallocating the tbis_ array doesn't cause a seg fault when this isn't called (we don't need it)
  SCF::init_threads();

  //----------------------------------------------------------------------------//
  // Set up the all pairs vector, needed to prescreen Schwarz bounds
  if(not dynamic_){
    const int nshell = gbs_->nshell();
    for(int ish=0, inode=0; ish < nshell; ++ish){
      for(int jsh=0; jsh <= ish; ++jsh, ++inode){
        IntPair ij(ish, jsh);
        pair_assignments_[AllPairs][ij] = inode % n_node;
      }
    }
    // Make the backwards mapping for the current node
    const int me = scf_grp_->me();
    for(auto it : pair_assignments_[AllPairs]){
      if(it.second == me){
        local_pairs_[AllPairs].push_back(it.first);
      }
    }
  }

  //----------------------------------------------------------------------------//

  init_significant_pairs();

  //----------------------------------------------------------------------------//
  if(not dynamic_){

    const int me = scf_grp_->me();
    int inode = 0;
    for(auto&& sig_pair : sig_pairs_) {
      const int assignment = inode % n_node;
      pair_assignments_[SignificantPairs][sig_pair] = assignment;
      if(assignment == me){
        local_pairs_[SignificantPairs].push_back(sig_pair);
      }
      ++inode;
    }

    // Make the assignments for the mu, X pairs in K
    for(int mu_set = 0; mu_set < sig_blocks_.size(); ++mu_set) {
      for(auto&& sig_block : sig_blocks_[mu_set]) {

        const int assignment = inode % n_node; ++inode;
        auto pair = std::make_pair(mu_set, sig_block);
        pair_assignments_k_[SignificantPairs][pair] = assignment;
        if(assignment == me) {
          local_pairs_k_[SignificantPairs].push_back(pair);
        }

      } // end loop over blocks for mu
    } // end loop over mu sets

    for(auto&& ish : shell_range(gbs_)) {
      for(auto&& Yblk : shell_block_range(dfbs_, 0, 0, NoLastIndex, NoRestrictions)) {
        const int assignment = inode % n_node; ++inode;
        auto pair = std::make_pair((int)ish, Yblk);
        pair_assignments_k_[AllPairs][pair] = assignment;
        if(assignment == me) {
          local_pairs_k_[AllPairs].push_back(pair);
        }
      }
    }
  }

  //----------------------------------------------------------------------------//
  threads_initialized_ = true;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::print(ostream&o) const
{
  o << indent << "Closed Shell Hartree-Fock (CLHF):" << endl;
  o << incindent;
  CLHF::print(o);
  if(print_screening_stats_) {
    stats_.print_summary(o, gbs_, dfbs_, print_screening_stats_);
  }
  o << decindent;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::done_threads(){
  CLHF::done_threads();
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::init_significant_pairs()
{
  Timer timer("init significant pairs");
  ExEnv::out0() << "  Computing significant shell pairs" << endl;
  std::atomic_int n_significant_pairs(0);
  const int nthread = threadgrp_->nthread();
  boost::mutex pair_mutex;
  std::vector<std::pair<double, IntPair>> pair_values;
  //----------------------------------------//
  schwarz_frob_.resize(gbs_->nshell(), gbs_->nshell());
  schwarz_frob_ = Eigen::MatrixXd::Zero(gbs_->nshell(), gbs_->nshell());
  local_pairs_spot_ = 0;
  do_threaded(nthread, [&](int ithr){

    ShellData ish, jsh;
    std::vector<std::pair<double, IntPair>> my_pair_vals;

    while(get_shell_pair(ish, jsh, AllPairs)){
      const double norm_val = ints4maxes_->get(ish, jsh, ish, jsh, coulomb_oper_type_, [&]() -> double {
        tbis_[ithr]->compute_shell(ish, jsh, ish, jsh);
        const double* buffer = tbis_[ithr]->buffer(coulomb_oper_type_);
        double frob_val = 0.0;
        for(int i = 0; i < ish.nbf*jsh.nbf*ish.nbf*jsh.nbf; ++i) {
          frob_val += fabs(buffer[i]);
        }
        return sqrt(frob_val);
      });
      schwarz_frob_(ish, jsh) = norm_val;
      schwarz_frob_(jsh, ish) = norm_val;
      my_pair_vals.push_back({norm_val, IntPair(ish, jsh)});
    } // end while get shell pair
    //----------------------------------------//
    // put our values on the node-level vector
    boost::lock_guard<boost::mutex> lg(pair_mutex);
    for(auto item : my_pair_vals){
      pair_values.push_back(item);
    }

  });
  //----------------------------------------//
  // All-to-all the shell-wise Frobenius norms of the Schwarz matrix
  scf_grp_->sum(schwarz_frob_.data(), gbs_->nshell() * gbs_->nshell());
  //const double schwarz_norm = schwarz_frob_.norm();
  // In this case we do actually want the max norm
  const double schwarz_norm = schwarz_frob_.maxCoeff();
  //----------------------------------------//
  // Now go through the list and figure out which ones are significant
  shell_to_sig_shells_.resize(gbs_->nshell());
  std::vector<double> sig_values;
  do_threaded(nthread, [&](int ithr){
    std::vector<IntPair> my_sig_pairs;
    std::vector<double> my_sig_values;
    for(int i = ithr; i < pair_values.size(); i += nthread){
      auto item = pair_values[i];
      if(item.first * schwarz_norm > pair_screening_thresh_){
        my_sig_pairs.push_back(item.second);
        my_sig_values.push_back(item.first);
        ++n_significant_pairs;
        ++stats_.sig_pairs;
        const int nfxn = gbs_->shell(item.second.first).nfunction() * gbs_->shell(item.second.second).nfunction();
        stats_.sig_pairs_fxn += nfxn;
        if(item.second.first != item.second.second) {
          ++stats_.sig_pairs;
          stats_.sig_pairs_fxn += nfxn;
        }

      }
    }
    //----------------------------------------//
    // put our values on the node-level vector
    boost::lock_guard<boost::mutex> lg(pair_mutex);
    for(auto item : my_sig_pairs){
      sig_pairs_.push_back(item);
    }
    for(auto item : my_sig_values){
      sig_values.push_back(item);
    }
  });
  //----------------------------------------//
  // Now all-to-all the significant pairs
  // This should be done with an MPI_alltoall_v or something like that
  int n_sig = n_significant_pairs;
  int my_n_sig = n_significant_pairs;
  int n_sig_pairs[scf_grp_->n()];
  for(int inode = 0; inode < scf_grp_->n(); n_sig_pairs[inode++] = 0);
  n_sig_pairs[scf_grp_->me()] = n_sig;
  scf_grp_->sum(n_sig_pairs, scf_grp_->n());
  scf_grp_->sum(n_sig);
  {
    int sig_data[n_sig*2];
    double sig_vals[n_sig];
    for(int inode = 0; inode < n_sig*2; sig_vals[inode] = 0.0, sig_data[inode++] = 0);
    int data_offset = 0;
    for(int inode = 0; inode < scf_grp_->me(); inode++){
      data_offset += 2 * n_sig_pairs[inode];
    }
    for(int ipair = 0; ipair < my_n_sig; ++ipair) {
      sig_data[data_offset + 2 * ipair] = sig_pairs_[ipair].first;
      sig_data[data_offset + 2 * ipair + 1] = sig_pairs_[ipair].second;
      sig_vals[data_offset/2 + ipair] = sig_values[ipair];
    }
    scf_grp_->sum(sig_data, n_sig*2);
    scf_grp_->sum(sig_vals, n_sig);
    sig_pairs_.clear();
    sig_values.clear();
    for(int ipair = 0; ipair < n_sig; ++ipair){
      sig_pairs_.push_back(IntPair(sig_data[2*ipair], sig_data[2*ipair + 1]));
      sig_values.push_back(sig_vals[ipair]);
    }
  } // get rid of sig_data and sig_vals
  //----------------------------------------//
  // Now compute the significant pairs for the outer loop of the exchange and the sig_partners_ array
  sig_partners_.resize(gbs_->nbasis());
  sig_blocks_.resize(gbs_->nbasis());
  int pair_index = 0;
  for(auto&& pair : sig_pairs_) {
    ShellData ish(pair.first, gbs_), jsh(pair.second, gbs_);
    sig_partners_[ish].insert(jsh);
    L_schwarz[ish].insert(jsh, schwarz_frob_(ish, jsh));
    if(ish != jsh) {
      sig_partners_[jsh].insert(ish);
      L_schwarz[jsh].insert(ish, schwarz_frob_(ish, jsh));
    }
    for(auto&& Xblk : iter_shell_blocks_on_center(dfbs_, ish.center)) {
      sig_blocks_[ish].insert(Xblk);
      sig_blocks_[jsh].insert(Xblk);
    }
    for(auto&& Xblk : iter_shell_blocks_on_center(dfbs_, jsh.center)) {
      sig_blocks_[ish].insert(Xblk);
      sig_blocks_[jsh].insert(Xblk);
    }
    //----------------------------------------//
    ++pair_index;
  }
  out_assert(L_schwarz.size(), >, 0);
  //----------------------------------------//
  do_threaded(nthread, [&](int ithr){
    auto L_schwarz_iter = L_schwarz.begin();
    const auto& L_schwarz_end = L_schwarz.end();
    L_schwarz_iter.advance(ithr);
    while(L_schwarz_iter != L_schwarz_end) {
      L_schwarz_iter->second.sort();
      L_schwarz_iter.advance(nthread);
    }
  });
  out_assert(L_schwarz[0].size(), >, 0);
  //----------------------------------------//
  max_schwarz_.resize(gbs_->nshell());
  for(auto&& ish : shell_range(gbs_)) {
    double max_val = 0.0;
    for(auto&& jsh : L_schwarz[ish]) {
      if(jsh.value > max_val) {
        max_val = jsh.value;
      }
      break;
    }
    max_schwarz_[ish] = max_val;
  }
  //----------------------------------------//
  centers_.resize(molecule()->natom());
  for(int iatom = 0; iatom < molecule()->natom(); ++iatom) {
    const double* r = molecule()->r(iatom);
    centers_[iatom] << r[0], r[1], r[2];
  }
  for(auto&& ish : shell_range(gbs_)) {
    const auto& ishell = gbs_->shell((int)ish);
    const std::vector<double>& i_exps = ishell.exponents();
    assert(ishell.ncontraction() == 1);
    for(auto&& jsh : iter_significant_partners(ish)) {
      const auto& jshell = gbs_->shell((int)jsh);
      const std::vector<double>& j_exps = jshell.exponents();
      //----------------------------------------//
      Eigen::Vector3d weighted_center;
      weighted_center << 0.0, 0.0, 0.0;
      double coef_tot = 0.0;
      for(int i_prim = 0; i_prim < ishell.nprimitive(); ++i_prim) {
        for(int j_prim = 0; j_prim < jshell.nprimitive(); ++j_prim) {
          const double coef_prod = fabs(ishell.coefficient_unnorm(0, i_prim)
              * jshell.coefficient_unnorm(0, j_prim));
          weighted_center += (coef_prod / (i_exps[i_prim] + j_exps[j_prim])) *
              (i_exps[i_prim] * centers_[ish.center] + j_exps[j_prim] * centers_[jsh.center]);
          coef_tot += coef_prod;
        }
      }
      //----------------------------------------//
      pair_centers_[{(int)ish, (int)jsh}] = (1.0 / coef_tot) * weighted_center;
    }
  }
  //----------------------------------------//
  ExEnv::out0() << "  Number of significant shell pairs:  " << n_sig << endl;
  ExEnv::out0() << "  Number of total basis pairs:  " << (gbs_->nshell() * (gbs_->nshell() + 1) / 2) << endl;
  ExEnv::out0() << "  Schwarz Norm = " << schwarz_norm << endl;

}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::save_data_state(StateOut& so)
{
  //----------------------------------------------------------------------------//
  // Nothing to do yet
  throw FeatureNotImplemented("SavableState writing in CADFCLHF",
      __FILE__, __LINE__, class_desc());
  //----------------------------------------------------------------------------//
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::ao_fock(double accuracy)
{
  /*=======================================================================================*/
  /* Setup                                                 		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  Timer timer("ao_fock");
  //---------------------------------------------------------------------------------------//
  if(not have_coefficients_) {
    ints_computed_locally_ = 0;
    compute_coefficients();
    ints_computed_ = ints_computed_locally_;
    scf_grp_->sum(ints_computed_);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << "  Computed " << ints_computed_ << " integrals to determine coefficients." << endl;
    }
  }
  //---------------------------------------------------------------------------------------//
  timer.enter("misc");
  int nthread = threadgrp_->nthread();
  iter_stats_ = &(stats_.next_iteration());
  //----------------------------------------//
  // transform the density difference to the AO basis
  RefSymmSCMatrix dd = cl_dens_diff_;
  Ref<PetiteList> pl = integral()->petite_list();
  cl_dens_diff_ = pl->to_AO_basis(dd);
  //----------------------------------------//
  double gmat_accuracy = accuracy;
  if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }
  //----------------------------------------//
  // copy over the density
  // cl_dens_diff_ includes total density right now, so halve it
  D_ = cl_dens_diff_.copy().convert2RefSCMat(); D_.scale(0.5);
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form G                                                		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  // compute J and K
  timer.change("build");
  if(xml_debug_){
    begin_xml_context("compute_fock", "compute_fock.xml");
  }
  RefSCMatrix G;
  {
    ints_computed_locally_ = 0;
    if(xml_debug_) begin_xml_context("compute_J");
    RefSCMatrix J = compute_J();
    if(xml_debug_) write_as_xml("J", J), end_xml_context("compute_J");
    G = J.copy();
    ints_computed_ = ints_computed_locally_;
    scf_grp_->sum(ints_computed_);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << "        Computed " << ints_computed_ << " integrals for J part" << endl;
    }
  }
  {
    ints_computed_locally_ = 0;
    if(xml_debug_) begin_xml_context("compute_K");
    RefSCMatrix K = compute_K();
    if(xml_debug_) write_as_xml("K", K), end_xml_context("compute_K");
    G.accumulate( -1.0 * K);
    ints_computed_ = ints_computed_locally_;
    scf_grp_->sum(ints_computed_);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << "        Computed " << ints_computed_ << " integrals for K part" << endl;
    }
  }
  if(xml_debug_) end_xml_context("compute_fock"), assert(false);
  density_reset_ = false;
  //---------------------------------------------------------------------------------------//
  // Move data back to a RefSymmSCMatrix, transform back to the SO basis
  Ref<SCElementOp> accum_G_op = new SCElementAccumulateSCMatrix(G.pointer());
  RefSymmSCMatrix G_symm = G.kit()->symmmatrix(G.coldim()); G_symm.assign(0.0);
  G_symm.element_op(accum_G_op); G = 0;
  G_symm = pl->to_SO_basis(G_symm);
  //---------------------------------------------------------------------------------------//
  // Accumulate difference back into gmat_
  gmat_.accumulate(G_symm); G_symm = 0;
  //---------------------------------------------------------------------------------------//
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Clean up                                             		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  timer.change("misc");
  //----------------------------------------//
  // restore the SO version of the density difference
  cl_dens_diff_ = dd;
  //----------------------------------------//
  // F = H+G
  cl_fock_.result_noupdate().assign(hcore_);
  cl_fock_.result_noupdate().accumulate(gmat_);
  accumddh_->accum(cl_fock_.result_noupdate());
  cl_fock_.computed()=1;
  //---------------------------------------------------------------------------------------//
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::reset_density()
{
  CLHF::reset_density();
  gmat_.assign(0.0);
  density_reset_ = true;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::loop_shell_pairs_threaded(
    PairSet pset,
    const std::function<void(int, const ShellData&, const ShellData&)>& f
)
{
  local_pairs_spot_ = 0;
  boost::thread_group compute_threads;
  const int nthread = threadgrp_->nthread();
  // Loop over number of threads
  for(int ithr = 0; ithr < nthread; ++ithr) {
    // ...and create each thread that computes pairs
    compute_threads.create_thread([&,ithr](){
      ShellData ish, jsh;
      //----------------------------------------//
      while(get_shell_pair(ish, jsh, pset)){
        f(ithr, ish, jsh);
      }

    });
  }
  compute_threads.join_all();
}

//////////////////////////////////////////////////////////////////////////////////

bool
CADFCLHF::get_shell_pair(ShellData& mu, ShellData& nu, PairSet pset)
{
  // Atomicly access and increment
  int spot = local_pairs_spot_++;
  if(spot < local_pairs_[pset].size()) {
    IntPair& next_pair = local_pairs_[pset][spot];
    //----------------------------------------//
    if(dynamic_) {
      // Here's where we'd need to check if we're running low on pairs and prefetch some more
      // When implemented, this should use a std::async or something like that
      throw FeatureNotImplemented("dynamic load balancing", __FILE__, __LINE__, class_desc());
    }
    //----------------------------------------//
    mu = ShellData(next_pair.first, gbs_.pointer(), dfbs_.pointer());
    nu = ShellData(next_pair.second, gbs_.pointer(), dfbs_.pointer());
  }
  else{
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////////


