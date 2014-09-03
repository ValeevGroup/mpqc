//
// init.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Mar 5, 2014
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

#include <random>
#include <algorithm>

// Boost includes
// NOTE:  THIS CAUSES VALGRIND ERRORS
#include <boost/math/special_functions/erf.hpp>

// MPQC includes
#include <chemistry/qc/basis/petite.h>
#include <util/group/messmpi.h>
#include <util/container/conc_cache.h>

#include "cadfclhf.h"
#include "assignments.h"
#include "ordered_shells.h"

using namespace sc;
using std::endl;

typedef std::pair<int, int> IntPair;

//////////////////////////////////////////////////////////////////////////////////

double my_erfc_inv(double val) {
  //if(well_separated_thresh_ == 0.1) {
  //  erfcinv_thr = 1.1630871536766740867262542605629475934779325500020816;
  //}
  //else {
  //  throw FeatureNotImplemented("Erfc_inv of number other than 0.1", __FILE__, __LINE__, class_desc());
  //}
  return boost::math::erfc_inv(val);
}

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
  have_coefficients_ = false;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::init_vector()
{
  CLHF::init_vector();
  //core_evals_ = basis_matrixkit()->diagmatrix(oso_dimension());
  //core_evecs_ = basis_matrixkit()->matrix(oso_dimension(), oso_dimension());
  //core_hamiltonian().diagonalize(core_evals_, core_evecs_);
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::init_threads()
{
  /*=======================================================================================*/
  /* Setup                                                		                        {{{1 */ #if 1 // begin fold

  Timer timer("init threads");
  assert(not threads_initialized_);

  ExEnv::out0() << indent << "Initializing CADFCLHF" << std::endl;
  ExEnv::out0() << incindent;
  ExEnv::out0() << indent << "nbf: " << gbs_->nbasis() << std::endl;
  ExEnv::out0() << indent << "dfnbf: " << dfbs_->nbasis() << std::endl;

  //----------------------------------------------------------------------------//
  // convenience variables

  const int me = scf_grp_->me();
  const int n_node = scf_grp_->n();
  const int nbf = gbs_->nbasis();
  const int nsh = gbs_->nshell();
  const int dfnbf = dfbs_->nbasis();
  const int dfnsh = dfbs_->nshell();

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/


  /*=======================================================================================*/
  /* Initialize 2-, 3-, and 4-center integral evaluators 		                          {{{1 */ #if 1 // begin fold
  //----------------------------------------------------------------------------//

  /*=======================================================================================*/
  /* Initialize two center TEI evaluators and compute (X|Y)                           {{{1 */ #if 1 // begin fold
  //----------------------------------------------------------------------------//
  // initialize the two electron integral classes

  ExEnv::out0() << indent << "Initializing 2 center integral evaluators" << std::endl;

  integral()->set_basis(dfbs_, dfbs_);

  size_t storage_required_2c = 0;
  try {
    // Note:  This overestimates, since we're using clone()
    storage_required_2c = integral()->storage_required(
        coulomb_oper_type_, TwoBodyIntShape::value::_1_O_2, 0,
        dfbs_, dfbs_
    ) * nthread_;
    if(coulomb_oper_type_ != metric_oper_type_) {
      storage_required_2c += integral()->storage_required(
          metric_oper_type_, TwoBodyIntShape::value::_1_O_2, 0,
          dfbs_, dfbs_
      ) * nthread_;

    }
    ExEnv::out0() << incindent << indent << "Integral object reports " << data_size_to_string(storage_required_2c)
                  << " required for 2 center integral evaluators." << decindent << endl;
  }
  catch(sc::Exception& e) {
    ExEnv::out0() << incindent << indent << "Integral object is not reporting the amount of"
                  << " storage needed for 2 center integral evaluators." << decindent << endl;
  }

  eris_2c_.resize(nthread_);
  eris_2c_[0] = integral()->coulomb<2>();
  for(int ithr = 1; ithr < nthread_; ++ithr) {
    eris_2c_[ithr] = eris_2c_[0]->clone();
  }
  for (int i=0; i < nthread_; i++) {
    if(metric_oper_type_ == coulomb_oper_type_){
      metric_ints_2c_.push_back(eris_2c_[i]);
    }
    else{
      throw FeatureNotImplemented("non-coulomb metrics in CADFCLHF", __FILE__, __LINE__, class_desc());
    }
  }
  consume_memory(storage_required_2c);

  //----------------------------------------------------------------------------//
  // Compute the two center integrals, then dispense with the evaluators
  // TODO this will need to be changed for different metric kernels

  ExEnv::out0() << indent << "Computing two center integrals" << std::endl;
  g2_full_ptr_ = ints_to_eigen_threaded(
      ShellBlockData<>(dfbs_),
      ShellBlockData<>(dfbs_),
      eris_2c_, coulomb_oper_type_
  );
  consume_memory(sizeof(TwoCenterIntContainer) + dfnbf*dfnbf*sizeof(double));
  const auto& g2 = *g2_full_ptr_;

  // Compute the (X|X)^1/2 schwarz matrix
  schwarz_df_.resize(dfnsh);
  consume_memory(dfnsh*sizeof(double));
  for(auto&& Xsh : shell_range(dfbs_)) {
    schwarz_df_(Xsh) = g2.block(Xsh.bfoff, Xsh.bfoff, Xsh.nbf, Xsh.nbf).norm();
  }

  // Release the integral evaluators
  for(int ithr = 0; ithr < nthread_; ++ithr) {
    eris_2c_[ithr] = 0;
    metric_ints_2c_[ithr] = 0;
  }
  release_memory(storage_required_2c);

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/


  /*=======================================================================================*/
  /* Initializes three center TEI evaluators		                                      {{{1 */ #if 1 // begin fold

  ExEnv::out0() << indent << "Initializing 3 center integral evaluators" << std::endl;

  // ThreeCenter versions
  integral()->set_basis(gbs_, gbs_, dfbs_);

  size_t storage_required_3c = 0;
  try {
    storage_required_3c = integral()->storage_required(
        coulomb_oper_type_, TwoBodyIntShape::value::_11_O_2, 0,
        gbs_, gbs_, dfbs_
    ) * nthread_;
    if(coulomb_oper_type_ != metric_oper_type_) {
      storage_required_3c += integral()->storage_required(metric_oper_type_, TwoBodyIntShape::value::_11_O_2, 0,
          gbs_, gbs_, dfbs_
      ) * nthread_;
    }
    ExEnv::out0() << incindent << indent << "Integral object reports " << data_size_to_string(storage_required_3c)
                  << " required for 3 center integral evaluators." << decindent << endl;
  }
  catch(sc::Exception& e) {
    ExEnv::out0() << incindent << indent << "Integral object is not reporting the amount of"
                  << " storage needed for 3 center integral evaluators." << decindent << endl;
  }

  size_t storage_avail = integral()->storage_unused();
  eris_3c_.resize(nthread_);

  eris_3c_[0] = integral()->coulomb<3>();
  for(int ithr = 1; ithr < nthread_; ++ithr) {
    eris_3c_[ithr] = eris_3c_[0]->clone();
    //eris_3c_[ithr] = integral()->coulomb<3>();
    //eris_3c_[ithr]->set_integral_storage(storage_avail/nthread_);
  }

  for (int i=0; i < nthread_; i++) {
    // TODO different fitting metrics
    if(metric_oper_type_ == coulomb_oper_type_){
      metric_ints_3c_.push_back(eris_3c_[i]);
    }
    else{
      throw FeatureNotImplemented("non-coulomb metrics in CADFCLHF", __FILE__, __LINE__, class_desc());
    }
  }
  memory_used_ += storage_required_3c;

  // Reset to normal setup
  integral()->set_basis(gbs_, gbs_, gbs_, gbs_);

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/


  /*=======================================================================================*/
  /* Initialize four center integral evaluators                 		                  {{{1 */ #if 1 // begin fold
  //----------------------------------------------------------------------------//
  ExEnv::out0() << indent << "Initializing 4 center integral evaluators" << endl;
  size_t storage_required_4c = 0;
  try {
    storage_required_4c = integral()->storage_required_eri(gbs_, gbs_, gbs_, gbs_);

    ExEnv::out0() << incindent << indent << "Integral object reports "
                  << data_size_to_string(storage_required_4c)
                  << " required for 4 center integral evaluators." << decindent << endl;
  }
  catch(sc::Exception& e) {
    ExEnv::out0() << incindent << indent << "Integral object is not reporting the amount of"
                  << " storage needed for 4 center integral evaluators." << decindent << endl;
  }
  bool need_nthr_4c_ints = thread_4c_ints_ or exact_diagonal_J_ or exact_diagonal_K_;
  tbis_ = new Ref<TwoBodyInt>[need_nthr_4c_ints ? nthread_ : 1];
  tbis_[0] = integral()->electron_repulsion();
  if(thread_4c_ints_) {
    for (int i=1; i < (need_nthr_4c_ints ? nthread_ : 1); i++) {
      tbis_[i] = tbis_[0]->clone();
    }
  }
  consume_memory(storage_required_4c);
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/


  /*=======================================================================================*/
  /* Compute static distributions of work and the significant pairs list              {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//

  /*-----------------------------------------------------*/
  /* Set up the all pairs vector                    {{{2 */ #if 2 // begin fold

  // Set up the all pairs vector, needed to prescreen Schwarz bounds
  ExEnv::out0() << indent << "Computing static distribution of all pairs" << endl;
  const int nshell = gbs_->nshell();
  for(int ish=0, inode=0; ish < nshell; ++ish){
    for(int jsh=0; jsh <= ish; ++jsh, ++inode){
      IntPair ij(ish, jsh);
      pair_assignments_[AllPairs][ij] = inode % n_node;
    }
  }
  // Make the backwards mapping for the current node
  for(auto it : pair_assignments_[AllPairs]){
    if(it.second == me){
      local_pairs_all_.push_back(it.first);
    }
  }

  /********************************************************/ #endif //2}}}
  /*-----------------------------------------------------*/

  //----------------------------------------------------------------------------//

  boost::shared_ptr<cadf::Node> my_part_ptr;
  if(new_exchange_algorithm_) {
    ExEnv::out0() << indent << "Computing static distribution of (obs, dfbs) pairs for new exchange" << endl;
    assignments_new_k_ = make_shared<cadf::assignments::Assignments>(
        gbs_, dfbs_, min_atoms_per_node_, n_node, me
    );
    assignments_new_k_->print_detail(ExEnv::out0());
    my_part_ptr = assignments_new_k_->nodes[me];
  }
  else {
    ExEnv::out0() << indent << "Computing static distribution of (obs, dfbs) pairs for exchange" << endl;
    atom_pair_assignments_k_ = make_shared<cadf::AssignmentGrid>(
        gbs_, dfbs_, scf_grp_->n(), scf_grp_->me()
    );
    atom_pair_assignments_k_->print_detail(ExEnv::out0(), !distribute_coefficients_);
    my_part_ptr = atom_pair_assignments_k_->my_assignments_ptr(me);
  }
  auto& my_part = *my_part_ptr;

  //----------------------------------------------------------------------------//

  ExEnv::out0() << indent << "Initializing significant basis function pairs" << endl;
  init_significant_pairs();

  // 4c int objects are destroyed in init_significant_pairs()
  if(not exact_diagonal_J_ and not exact_diagonal_K_) {
    release_memory(storage_required_4c);
  }

  //----------------------------------------------------------------------------//

  ExEnv::out0() << indent << "Computing static distribution of significant pairs for Coulomb" << endl;
  int inode = 0;
  // TODO More efficient distribution based on load balancing and minimizing thread collisions
  if(shuffle_J_assignments_) {
    std::random_shuffle(sig_pairs_.begin(), sig_pairs_.end());
  }
  for(auto&& sig_pair : sig_pairs_) {
    const int assignment = inode % n_node;
    pair_assignments_[SignificantPairs][sig_pair] = assignment;
    if(assignment == me){
      ShellData ish(sig_pair.first, gbs_);
      ShellData jsh(sig_pair.second, gbs_);
      if(ish.nbf > max_fxn_obs_j_ish_) max_fxn_obs_j_ish_ = ish.nbf;
      if(jsh.nbf > max_fxn_obs_j_jsh_) max_fxn_obs_j_jsh_ = jsh.nbf;
      local_pairs_sig_.push_back(sig_pair);
    }
    ++inode;
  }

  schwarz_df_mine_.resize(my_part.dfnsh());
  out_assert(my_part.dfnsh(), ==, my_part.assigned_dfbs_shells().size());
  int Xsh_off = 0;
  for(auto&& Xsh_index : my_part.assigned_dfbs_shells()) {
    ShellData Xsh(Xsh_index, dfbs_, gbs_);
    schwarz_df_mine_(Xsh_off) = g2.block(Xsh.bfoff, Xsh.bfoff, Xsh.nbf, Xsh.nbf).norm();
    ++Xsh_off;
  }

  // Get maximum sizes for intermediate data blocks
  for(auto&& pair : my_part.pairs) {
    {
      ShellBlockData<> Xblk = ShellBlockData<>::atom_block(pair.Xatom, dfbs_, gbs_);
      if(Xblk.nbf > max_fxn_atom_dfbs_todo_) max_fxn_atom_dfbs_todo_ = Xblk.nbf;
      if(Xblk.atom_obsnbf > max_obs_atom_fxn_on_dfbs_center_todo_) {
        max_obs_atom_fxn_on_dfbs_center_todo_ = Xblk.atom_obsnbf;
      }
      ShellData ish(pair.ish, gbs_, dfbs_);
      if(ish.nbf > max_fxn_obs_todo_) max_fxn_obs_todo_ = ish.nbf;
      if(do_linK_) {
        for(auto&& Xsh : shell_range(Xblk)) {
          if(Xsh.nbf > max_fxn_dfbs_todo_) max_fxn_dfbs_todo_ = Xsh.nbf;
          local_pairs_linK_.emplace((int)ish, (int)Xsh);
          linK_local_map_[Xsh].push_back((int)ish);
          linK_local_map_ish_Xsh_[ish].push_back((int)Xsh);
        }
      }
      else {
        local_pairs_k_.emplace_back(pair.ish, Xblk);
      }


    }
  }
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/


  /*=======================================================================================*/
  /* Summary and cleanup                                  		                        {{{1 */ #if 1 // begin fold

  threads_initialized_ = true;

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Done initializing CADFCLHF" << endl;

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::init_significant_pairs()
{

  /*=======================================================================================*/
  /* Setup                                                  		                      {{{1 */ #if 1 // begin fold
  Timer timer("init significant pairs");

  boost::shared_ptr<cadf::Node> my_part_ptr;
  if(new_exchange_algorithm_) {
    my_part_ptr = assignments_new_k_->nodes[scf_grp_->me()];
  }
  else {
    my_part_ptr = atom_pair_assignments_k_->my_assignments_ptr(scf_grp_->me());
  }
  auto& my_part = *my_part_ptr;

  ExEnv::out0() << incindent;
  ExEnv::out0() << indent << "Computing Schwarz matrix" << endl;

  //----------------------------------------------------------------------------//
  // convenience variables

  const int me = scf_grp_->me();
  const int n_node = scf_grp_->n();
  const int nbf = gbs_->nbasis();
  const int nsh = gbs_->nshell();
  const int dfnbf = dfbs_->nbasis();
  const int dfnsh = dfbs_->nshell();

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/

  /*=======================================================================================*/
  /* Compute the Schwarz matrix in threads and in a distributed manner                {{{1 */ #if 1 // begin fold

  std::vector<std::pair<double, IntPair>> pair_values;
  boost::mutex pair_mutex;
  schwarz_frob_.resize(gbs_->nshell(), gbs_->nshell());
  memory_used_ += gbs_->nshell() * gbs_->nshell() * sizeof(double);
  schwarz_frob_ = Eigen::MatrixXd::Zero(gbs_->nshell(), gbs_->nshell());
  local_pairs_spot_ = 0;

  typedef Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>> ConstVectorMap;
  bool need_nthr_4c_ints = thread_4c_ints_ or exact_diagonal_J_ or exact_diagonal_K_;
  do_threaded((need_nthr_4c_ints ? nthread_ : 1), [&](int ithr){

    ShellData ish, jsh;
    std::vector<std::pair<double, IntPair>> my_pair_vals;

    while(get_shell_pair(ish, jsh, AllPairs)){
      tbis_[ithr]->compute_shell(ish, jsh, ish, jsh);
      const double* buffer = tbis_[ithr]->buffer(coulomb_oper_type_);
      const ConstVectorMap buff_map(buffer, ish.nbf*jsh.nbf*ish.nbf*jsh.nbf);
      const double norm_val = sqrt(buff_map.cwiseAbs().sum());
      schwarz_frob_(ish, jsh) = norm_val;
      if(ish != jsh) schwarz_frob_(jsh, ish) = norm_val;
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
  // Delete the tbis_ if we can
  if(not exact_diagonal_J_ and not exact_diagonal_K_) {
    // At this point, we're done with the tbis_
    for (int i=0; i < (thread_4c_ints_ ? nthread_ : 1); i++) tbis_[i] = 0;
    delete[] tbis_;
    tbis_ = 0;
  }

  //----------------------------------------//
  // All-to-all the shell-wise Frobenius norms of the Schwarz matrix
  ExEnv::out0() << indent << "Distributing Schwarz matrix" << endl;
  scf_grp_->sum(schwarz_frob_.data(), nsh * nsh);

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/

  /*=======================================================================================*/
  /* Form the significant pair lists                      		                        {{{1 */ #if 1 // begin fold

  //----------------------------------------//
  // Get the Schwarz norm
  // In this case we do actually want the max norm, since we
  //   want to determine which pairs have ANY significant quartets
  //   that they could be a part of.
  const double schwarz_norm = schwarz_frob_.maxCoeff();

  /*-----------------------------------------------------*/
  /* Figure out which ones are significant          {{{2 */ #if 2 // begin fold

  // Go through the list and figure out which ones are significant
  std::vector<double> sig_values;
  std::atomic_int n_significant_pairs(0);
  do_threaded(nthread_, [&](int ithr){

    std::vector<IntPair> my_sig_pairs;
    std::vector<double> my_sig_values;
    for(int i = ithr; i < pair_values.size(); i += nthread_){
      auto item = pair_values[i];
      if(item.first * schwarz_norm > pair_screening_thresh_ or count_ints_only_){
        my_sig_pairs.push_back(item.second);
        my_sig_values.push_back(item.first);
        ++n_significant_pairs;
        if(print_screening_stats_) {
          ++(stats_->sig_pairs);
          const int nfxn = gbs_->shell(item.second.first).nfunction() * gbs_->shell(item.second.second).nfunction();
          stats_->sig_pairs_fxn += nfxn;
          if(item.second.first != item.second.second) {
            ++(stats_->sig_pairs);
            stats_->sig_pairs_fxn += nfxn;
          }
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

  /********************************************************/ #endif //2}}}
  /*-----------------------------------------------------*/

  /*-----------------------------------------------------*/
  /* All-to-all the significant pairs               {{{2 */ #if 2 // begin fold
  // This should be done with an MPI_alltoall_v or something like that

  ExEnv::out0() << indent << "Distributing significant pairs list" << endl;

  int n_sig = n_significant_pairs;
  int my_n_sig = n_significant_pairs;
  int n_sig_pairs[scf_grp_->n()];
  for(int inode = 0; inode < scf_grp_->n(); n_sig_pairs[inode++] = 0);

  n_sig_pairs[scf_grp_->me()] = n_sig;
  scf_grp_->sum(n_sig_pairs, scf_grp_->n());
  scf_grp_->sum(n_sig);
  {

    int* sig_data = new int[n_sig*2];
    double* sig_vals = new double[n_sig];
    for(int inode = 0; inode < n_sig*2; sig_vals[inode/2] = 0.0, sig_data[inode++] = 0);
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
    delete[] sig_data;
    delete[] sig_vals;

  } // get rid of sig_data and sig_vals
  /********************************************************/ #endif //2}}}
  /*-----------------------------------------------------*/

  /*-----------------------------------------------------*/
  /* Build the significant partners arrays          {{{2 */ #if 2 // begin fold

  ExEnv::out0() << indent << "Computing the sig partners array" << endl;
  sig_partners_.resize(gbs_->nshell());

  for(auto&& pair : sig_pairs_) {
    ShellData ish(pair.first, gbs_), jsh(pair.second, gbs_);
    sig_partners_[ish].insert(jsh);
    L_schwarz[ish].insert(jsh, schwarz_frob_(ish, jsh));
    if(ish != jsh) {
      sig_partners_[jsh].insert(ish);
      L_schwarz[jsh].insert(ish, schwarz_frob_(ish, jsh));
    }
  }

  for(auto&& idxlist : sig_partners_) {
    sig_partner_blocks_.emplace_back(
        idxlist, gbs_, dfbs_,
        std::set<int>{NoRestrictions, SameCenter}
    );
  }

  //----------------------------------------//
  // Sort the Schwarz lists
  do_threaded(nthread_, [&](int ithr){
    auto L_schwarz_iter = L_schwarz.begin();
    const auto& L_schwarz_end = L_schwarz.end();
    L_schwarz_iter.advance(ithr);
    while(L_schwarz_iter != L_schwarz_end) {
      L_schwarz_iter->second.sort();
      L_schwarz_iter.advance(nthread_);
    }
  });

  /********************************************************/ #endif //2}}}
  /*-----------------------------------------------------*/

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/

  /*=======================================================================================*/
  /* Compute the centers, the pair centers, overlap, and the pair extents             {{{1 */ #if 1 // begin fold

  ExEnv::out0() << indent << "Computing overlap" << endl;
  // Get the overlap integrals
  if(dist_factor_use_overlap_){
    resize_and_zero_matrix(S_frob_, nsh, nsh);
    integral()->set_basis(gbs_, gbs_, gbs_, gbs_);
    Ref<OneBodyInt> sint = integral()->overlap();
    for(auto&& ish : shell_range(gbs_)) {
      for(auto&& jsh : shell_range(gbs_)) {
        sint->compute_shell((int)ish, (int)jsh);
        const Eigen::Map<const Eigen::VectorXd> buffmap(sint->buffer(), ish.nbf*jsh.nbf);
        S_frob_(ish, jsh) = buffmap.norm();
      }
    }
  }

  ExEnv::out0() << indent << "Computing pair centers and extents" << endl;
  centers_.resize(molecule()->natom());
  for(int iatom = 0; iatom < molecule()->natom(); ++iatom) {
    const double* r = molecule()->r(iatom);
    centers_[iatom] << r[0], r[1], r[2];
  }


  double erfcinv_thr = my_erfc_inv(well_separated_thresh_);

  for(auto&& ish : shell_range(my_part.obs_shells_to_do, gbs_, dfbs_)) {

    const auto& ishell = gbs_->shell((int)ish);
    const std::vector<double>& i_exps = ishell.exponents();

    if(ishell.ncontraction() != 1) {
      throw FeatureNotImplemented("Generally contracted basis sets", __FILE__, __LINE__, class_desc());
    }

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
      auto& r_ij = (1.0 / coef_tot) * weighted_center;
      pair_centers_[{(int)ish, (int)jsh}] = r_ij;

      if(use_extents_) {
        double min_exponent = std::numeric_limits<double>::infinity();
        double extent_max = 0.0;
        double extent_sum = 0.0;
        double avg_expon = 0.0;

        for(int i_prim = 0; i_prim < ishell.nprimitive(); ++i_prim) {
          for(int j_prim = 0; j_prim < jshell.nprimitive(); ++j_prim) {
            const double zeta_p = i_exps[i_prim];
            const double zeta_q = j_exps[j_prim];
            //const double reduced_exponent = (zeta_p * zeta_q) / (zeta_p + zeta_q);
            const double reduced_exponent = fabs(zeta_p + zeta_q);
            if(reduced_exponent < min_exponent) min_exponent = reduced_exponent;
            auto& r_pq = (1.0 / (i_exps[i_prim] + j_exps[j_prim])) *
                (i_exps[i_prim] * centers_[ish.center] + j_exps[j_prim] * centers_[jsh.center]);
            const double rdiff = (r_pq - r_ij).norm();
            const double ext_pq = sqrt(2.0 / (zeta_p + zeta_q)) + rdiff;
            if(ext_pq > extent_max) {
              extent_max = ext_pq;
            }
            const double coef_prod = fabs(ishell.coefficient_unnorm(0, i_prim)
                * jshell.coefficient_unnorm(0, j_prim));
            avg_expon += reduced_exponent * coef_prod;
            extent_sum += coef_prod * ext_pq;
          }
        }

        effective_pair_exponents_[{(int)ish, (int)jsh}] = min_exponent;
        //effective_pair_exponents_[{(int)ish, (int)jsh}] = avg_expon / coef_tot;

        if(use_max_extents_) {
          pair_extents_[{(int)ish, (int)jsh}] = extent_max * erfcinv_thr;
        }
        else {
          pair_extents_[{(int)ish, (int)jsh}] = extent_sum/coef_tot * erfcinv_thr;
        }
      }

    }
  }

  if(use_extents_) {
    df_extents_.reserve(dfbs_->nshell());
    for(auto&& Xsh : shell_range(dfbs_)) {

      if(safe_extents_) {
        erfcinv_thr = my_erfc_inv(pow(well_separated_thresh_, double(1+Xsh.am)));
      }

      const auto& Xshell = dfbs_->shell((int)Xsh);
      if(Xshell.ncontraction() != 1) {
        throw FeatureNotImplemented("Generally contracted basis sets", __FILE__, __LINE__, class_desc());
      }

      const std::vector<double>& x_exps = Xshell.exponents();

      double ext_max = 0.0;
      double ext_sum = 0.0;
      double coef_sum = 0.0;
      double min_exponent = std::numeric_limits<double>::infinity();

      for(int x_prim = 0; x_prim < Xshell.nprimitive(); ++x_prim) {
        const double zeta_x = fabs(x_exps[x_prim]);
        if(zeta_x < min_exponent) min_exponent = zeta_x;
        const double coef = fabs(Xshell.coefficient_unnorm(0, x_prim));
        double ext = sqrt(2.0 / zeta_x);
        if(ext > ext_max) {
          ext_max = ext;
        }
        coef_sum += coef;
        ext_sum += coef * ext;
      }
      if(use_max_extents_) {
        df_extents_[(int)Xsh] = ext_max * erfcinv_thr;
        //const double l_X = Xsh.am;
        //df_extents_[(int)Xsh] = ext_max * erfcinv_thr *
        //    (1.0 + 0.023215658687093542 * l_X - 0.00008543584437817054 * l_X * l_X +
        //     0.24989969347856592 * log(1 + 0.5114639422780547 * l_X));
      }
      else {
        df_extents_[(int)Xsh] = ext_sum/coef_sum * erfcinv_thr;
      }
      effective_df_exponents_.push_back(min_exponent);
    }
  }


  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/

  /*=======================================================================================*/
  /* Summary and cleanup                                  		                        {{{1 */ #if 1 // begin fold

  ExEnv::out0() << indent << "Number of significant shell pairs:  " << n_sig << endl;
  ExEnv::out0() << indent << "  Number of total basis pairs:  " << (gbs_->nshell() * (gbs_->nshell() + 1) / 2) << endl;
  ExEnv::out0() << indent << "  Schwarz Norm = " << schwarz_norm << endl;
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Done computing significant pairs" << endl;

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/

}


