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

// Boost includes
#include <boost/math/special_functions/erf.hpp>

// MPQC includes
#include <chemistry/qc/basis/petite.h>
#include <util/group/messmpi.h>
#include <util/container/conc_cache.h>

#include "cadfclhf.h"
#include "assignments.h"

using namespace sc;
using std::endl;

typedef std::pair<int, int> IntPair;

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

  ExEnv::out0() << indent << "Initializing CADFCLHF" << std::endl;
  ExEnv::out0() << incindent;

  //----------------------------------------------------------------------------//

  const int n_node = scf_grp_->n();

  //----------------------------------------------------------------------------//
  // initialize the two electron integral classes

  ExEnv::out0() << indent << "Initializing 3 center integral evaluators" << std::endl;

  // ThreeCenter versions
  integral()->set_basis(gbs_, gbs_, dfbs_);
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

  ExEnv::out0() << indent << "Initializing 2 center integral evaluators" << std::endl;

  // TwoCenter versions
  integral()->set_basis(dfbs_, dfbs_);
  eris_2c_.resize(nthread_);
  eris_2c_[0] = integral()->coulomb<2>();
  for(int ithr = 1; ithr < nthread_; ++ithr) {
    eris_2c_[ithr] = eris_2c_[0]->clone();
    //eris_2c_[ithr] = integral()->coulomb<2>();
  }
  for (int i=0; i < nthread_; i++) {
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
  ExEnv::out0() << indent << "Initializing 4 center integral evaluators" << endl;
  //SCF::init_threads();
  tbis_ = new Ref<TwoBodyInt>[nthread_];
  tbis_[0] = integral()->electron_repulsion();
  for (int i=1; i < nthread_; i++) {
    tbis_[i] = tbis_[0]->clone();
  }

  //----------------------------------------------------------------------------//
  // Set up the all pairs vector, needed to prescreen Schwarz bounds
  if(not dynamic_){
    ExEnv::out0() << indent << "Computing static distribution of all pairs" << endl;
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
        local_pairs_all_.push_back(it.first);
      }
    }
  }

  //----------------------------------------------------------------------------//

  ExEnv::out0() << indent << "Initializing significant basis function pairs" << endl;
  init_significant_pairs();

  //----------------------------------------------------------------------------//

  ExEnv::out0() << indent << "Computing static distribution of significant pairs for Coulomb" << endl;
  const int me = scf_grp_->me();
  int inode = 0;
  for(auto&& sig_pair : sig_pairs_) {
    const int assignment = inode % n_node;
    pair_assignments_[SignificantPairs][sig_pair] = assignment;
    if(assignment == me){
      local_pairs_sig_.push_back(sig_pair);
    }
    ++inode;
  }

  ExEnv::out0() << indent << "Computing static distribution of (obs, dfbs) pairs for exchange" << endl;
  atom_pair_assignments_k_ = make_shared<cadf::AssignmentGrid>(
      gbs_, dfbs_, scf_grp_->n()
  );
  atom_pair_assignments_k_->print_detail();
  auto& my_part = atom_pair_assignments_k_->my_assignments(me);

  for(auto&& ish_ptr : my_part.bin->assigned_obs_shells) {
    ShellData ish(ish_ptr->index, gbs_, dfbs_);
    for(auto&& Xatom_ptr : my_part.bin->assigned_dfbs_atoms) {
      ShellBlockData<> Xblk = ShellBlockData<>::atom_block(Xatom_ptr->index, dfbs_, gbs_);
      if(do_linK_) {
        for(auto&& Xsh : shell_range(Xblk)) {
          local_pairs_linK_.emplace(ish, (int)Xsh);
          linK_local_map_[Xsh].push_back(ish);
        }
      }
      else {
        local_pairs_k_.emplace_back(ish, Xblk);
      }


    }
  }


  //----------------------------------------------------------------------------//
  threads_initialized_ = true;

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Done initializing CADFCLHF" << endl;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::init_significant_pairs()
{
  Timer timer("init significant pairs");

  ExEnv::out0() << incindent;
  ExEnv::out0() << indent << "Computing Schwarz matrix" << endl;

  std::atomic_int n_significant_pairs(0);
  boost::mutex pair_mutex;
  std::vector<std::pair<double, IntPair>> pair_values;
  //----------------------------------------//
  schwarz_frob_.resize(gbs_->nshell(), gbs_->nshell());
  memory_used_ += gbs_->nshell() * gbs_->nshell() * sizeof(double);
  schwarz_frob_ = Eigen::MatrixXd::Zero(gbs_->nshell(), gbs_->nshell());
  local_pairs_spot_ = 0;
  do_threaded(nthread_, [&](int ithr){

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

  ExEnv::out0() << indent << "Distributing Schwarz matrix" << endl;

  scf_grp_->sum(schwarz_frob_.data(), gbs_->nshell() * gbs_->nshell());
  //const double schwarz_norm = schwarz_frob_.norm();
  // In this case we do actually want the max norm
  const double schwarz_norm = schwarz_frob_.maxCoeff();
  //----------------------------------------//
  // Now go through the list and figure out which ones are significant
  shell_to_sig_shells_.resize(gbs_->nshell());
  std::vector<double> sig_values;
  do_threaded(nthread_, [&](int ithr){

    std::vector<IntPair> my_sig_pairs;
    std::vector<double> my_sig_values;
    for(int i = ithr; i < pair_values.size(); i += nthread_){
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

  //----------------------------------------//

  //============================================================================//
  // Now compute the significant pairs for the outer loop of the exchange and the sig_partners_ array
  ExEnv::out0() << indent << "Computing the sig partners array" << endl;
  sig_partners_.resize(gbs_->nshell());
  int pair_index = 0;
  for(auto&& pair : sig_pairs_) {
    ShellData ish(pair.first, gbs_), jsh(pair.second, gbs_);
    sig_partners_[ish].insert(jsh);
    L_schwarz[ish].insert(jsh, schwarz_frob_(ish, jsh));
    if(ish != jsh) {
      sig_partners_[jsh].insert(ish);
      L_schwarz[jsh].insert(ish, schwarz_frob_(ish, jsh));
    }
    //----------------------------------------//
    ++pair_index;
  }
  //----------------------------------------//
  do_threaded(nthread_, [&](int ithr){
    auto L_schwarz_iter = L_schwarz.begin();
    const auto& L_schwarz_end = L_schwarz.end();
    L_schwarz_iter.advance(ithr);
    while(L_schwarz_iter != L_schwarz_end) {
      L_schwarz_iter->second.sort();
      L_schwarz_iter.advance(nthread_);
    }
  });
  //----------------------------------------//
  // compute the max Schwarz frobnorm ( mu rho | mu rho ) for each ish
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

  //============================================================================//
  // Compute the centers, the pair centers, and the pair extents
  ExEnv::out0() << indent << "Computing pair centers and extents" << endl;
  centers_.resize(molecule()->natom());
  for(int iatom = 0; iatom < molecule()->natom(); ++iatom) {
    const double* r = molecule()->r(iatom);
    centers_[iatom] << r[0], r[1], r[2];
  }

  const double erfcinv_thr = boost::math::erfc_inv(well_separated_thresh_);

  for(auto&& ish : shell_range(gbs_)) {

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
        double extent_max = 0.0;
        double extent_sum = 0.0;

        for(int i_prim = 0; i_prim < ishell.nprimitive(); ++i_prim) {
          for(int j_prim = 0; j_prim < jshell.nprimitive(); ++j_prim) {
            const double zeta_p = i_exps[i_prim];
            const double zeta_q = j_exps[j_prim];
            auto& r_pq = (1.0 / (i_exps[i_prim] + j_exps[j_prim])) *
                (i_exps[i_prim] * centers_[ish.center] + j_exps[j_prim] * centers_[jsh.center]);
            const double rdiff = (r_pq - r_ij).norm();
            const double ext_pq = sqrt(2.0 / (zeta_p + zeta_q)) + rdiff;
            if(ext_pq > extent_max) {
              extent_max = ext_pq;
            }
            const double coef_prod = fabs(ishell.coefficient_unnorm(0, i_prim)
                * jshell.coefficient_unnorm(0, j_prim));
            extent_sum += coef_prod * ext_pq;
          }
        }

        if(use_max_extents_) {
          pair_extents_[{(int)ish, (int)jsh}] = extent_max * erfcinv_thr;
        }
        else {
          pair_extents_[{(int)ish, (int)jsh}] = extent_sum/coef_tot * erfcinv_thr;
        }
      }

    }

    if(use_extents_) {
      df_extents_.reserve(dfbs_->nshell());
      for(auto&& Xsh : shell_range(dfbs_)) {

        const auto& Xshell = dfbs_->shell((int)Xsh);
        if(Xshell.ncontraction() != 1) {
          throw FeatureNotImplemented("Generally contracted basis sets", __FILE__, __LINE__, class_desc());
        }

        const std::vector<double>& x_exps = Xshell.exponents();

        double ext_max = 0.0;
        double ext_sum = 0.0;
        double coef_sum = 0.0;

        for(int x_prim = 0; x_prim < Xshell.nprimitive(); ++x_prim) {
          const double zeta_x = x_exps[x_prim];
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
        }
        else {
          df_extents_[(int)Xsh] = ext_sum/coef_sum * erfcinv_thr;
        }
      }
    }

  }


  //----------------------------------------//

  ExEnv::out0() << indent << "Number of significant shell pairs:  " << n_sig << endl;
  ExEnv::out0() << indent << "  Number of total basis pairs:  " << (gbs_->nshell() * (gbs_->nshell() + 1) / 2) << endl;
  ExEnv::out0() << indent << "  Schwarz Norm = " << schwarz_norm << endl;
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Done computing significant pairs" << endl;

}


