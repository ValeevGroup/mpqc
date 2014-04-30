//
// compute_k.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Feb 13, 2014
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

#include <numeric>
#include <cassert>

#include <chemistry/qc/basis/petite.h>
#include <util/misc/xmlwriter.h>

#include "cadfclhf.h"
#include "assignments.h"

using namespace sc;
using std::cout;
using std::endl;


#define DEBUG_K_INTERMEDIATES 0


double CADFCLHF::get_distance_factor(
    const ShellData& ish, const ShellData& jsh, const ShellData& Xsh
) const
{
  double dist_factor;
  if(linK_use_distance_){
    const double R = get_R(ish, jsh, Xsh);
    const double l_X = double(dfbs_->shell(Xsh).am(0));
    double r_expon = l_X + 1.0;
    if(ish.center == jsh.center) {
      if(ish.center == Xsh.center) {
        r_expon = 0.0;
      }
      else {
        const double l_i = double(gbs_->shell(ish).am(0));
        const double l_j = double(gbs_->shell(jsh).am(0));
        r_expon += abs(l_i-l_j);
      }
    }
    dist_factor = 1.0 / (pow(R, pow(r_expon, distance_damping_factor_)));            //latex `\label{sc:link:dist_damp}`
    // If the distance factor actually makes the bound larger, then ignore it.
    dist_factor = std::min(1.0, dist_factor);
  }
  else {
    dist_factor = 1.0;
  }
  return dist_factor;
};


double CADFCLHF::get_R(
    const ShellData& ish,
    const ShellData& jsh,
    const ShellData& Xsh
) const
{
  double rv = (pair_centers_.at({(int)ish, (int)jsh}) - centers_[Xsh.center]).norm();
  if(use_extents_) {
    const double ext_a = pair_extents_.at({(int)ish, (int)jsh});
    const double ext_b = df_extents_[(int)Xsh];
    if(subtract_extents_) {
      rv -= ext_a + ext_b;
    }
    else if(ext_a + ext_b >= rv){
      rv = 1.0; // Don't do distance screening
    }
  }
  if(rv < 1.0) {
    rv = 1.0;
  }
  return rv;
};


typedef std::pair<int, int> IntPair;


RefSCMatrix
CADFCLHF::compute_K()
{
  /*=======================================================================================*/
  /* Setup                                                 		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // Convenience variables
  Timer timer("compute K");
  const int me = scf_grp_->me();
  const int n_node = scf_grp_->n();
  const Ref<GaussianBasisSet>& obs = gbs_;
  const int nbf = obs->nbasis();
  const int dfnbf = dfbs_->nbasis();
  const int natom = obs->ncenter();
  auto& my_part = atom_pair_assignments_k_->my_assignments(me);

  //----------------------------------------//
  // Get the density in an Eigen::Map form
  //double* D_data = allocate<double>(nbf*nbf);
  double* __restrict__ D_data = new double[nbf*nbf];
  //double* D_data = new double[nbf*nbf];
  D_.convert(D_data);
  typedef Eigen::Map<Eigen::VectorXd> VectorMap;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrix;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMatrix;
  //typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> MatrixMap;
  typedef Eigen::Map<ColMatrix> MatrixMap;
  // Matrix and vector wrappers, for convenience
  VectorMap d(D_data, nbf*nbf);
  MatrixMap D(D_data, nbf, nbf);
  // Match density scaling in old code:
  D *= 0.5;

  //----------------------------------------//
  Eigen::MatrixXd Kt(nbf, nbf);
  Kt = Eigen::MatrixXd::Zero(nbf, nbf);

  //----------------------------------------//
  const auto& g2 = *g2_full_ptr_;

  //----------------------------------------//
  // if we're doing the exact diagonal, form the Z intermediate
  std::vector<RowMatrix> Z;
  std::vector<RowMatrix> Z_tilde;
  std::vector<RowMatrix> Z_tilde_bar;
  size_t Z_size = 0;
  if(exact_diagonal_K_) {
    timer.enter("compute Z for exact diagonal");

    for(int iatom = 0; iatom < natom; ++iatom) {

      const int atom_nbf = gbs_->nbasis_on_center(iatom);

      Z.emplace_back(atom_nbf*atom_nbf, dfnbf);
      auto& Z_iatom = Z.back();
      Z_iatom = RowMatrix::Zero(atom_nbf*atom_nbf, dfnbf);
      Z_size += atom_nbf*atom_nbf*dfnbf*sizeof(double) + sizeof(RowMatrix);

      for(auto&& X : function_range(dfbs_, gbs_)) {
        for(auto&& nu : iter_functions_on_center(obs, iatom)) {
          Z_iatom.col(X).segment(nu.bfoff_in_atom*nu.atom_nbf, nu.atom_nbf) += 2.0 *
              D.block(nu.atom_bfoff, X.atom_obsbfoff, nu.atom_nbf, X.atom_obsnbf)
                * coefs_transpose_[X].col(nu);
        }
      }

    } // end loop over atoms


    for(auto&& X : function_range(dfbs_, gbs_)) {

      // Use the orbital basis as the "auxiliary" for the dfbs, just for convenience
      const int atom_nbf = X.atom_obsnbf;

      Z_tilde.emplace_back(atom_nbf, nbf);
      Z_tilde_bar.emplace_back(atom_nbf, nbf);
      auto& Zt_iatom = Z_tilde.back();
      auto& Ztb_iatom = Z_tilde_bar.back();
      Zt_iatom = RowMatrix::Zero(atom_nbf, nbf);
      Ztb_iatom = RowMatrix::Zero(atom_nbf, nbf);
      Z_size += 2 * X.atom_dfnbf*nbf*sizeof(double) + sizeof(RowMatrix);

      for(auto&& jblk : shell_block_range(obs, dfbs_, 0, NoLastIndex, SameCenter)) {
        Zt_iatom.middleCols(jblk.bfoff, jblk.nbf) += 2.0 *
            coefs_transpose_[X].middleCols(jblk.atom_bfoff, jblk.atom_nbf)
            * D.block(jblk.bfoff, jblk.atom_bfoff, jblk.nbf, jblk.atom_nbf);
      }
      Ztb_iatom += 2.0 * D.block(X.atom_obsbfoff, X.atom_obsbfoff, X.atom_obsnbf, X.atom_obsnbf) * coefs_transpose_[X];
    } // end loop over X in dfbs

    if(xml_debug_) {
      for(int iatom = 0; iatom < natom; ++iatom) {
        for(auto&& nu : iter_functions_on_center(obs, iatom)) {
          for(auto&& rho : iter_functions_on_center(obs, iatom)) {
            write_as_xml("Z", Z[iatom].row(nu.bfoff_in_atom*rho.atom_nbf + rho.bfoff_in_atom), attrs<int>{
              {"ao_index1", nu},
              {"ao_index2", rho}
            });
          }
        }
      }
      for(auto&& X : function_range(dfbs_)) {
        for(auto&& nu : iter_functions_on_center(obs, X.center)) {
          write_as_xml("Z_tilde", Z_tilde[X].row(nu.bfoff_in_atom), attrs<int>{
              {"ao_index1", X},
              {"ao_index2", nu}
          });
        }
      }
      for(auto&& X : function_range(dfbs_)) {
        for(auto&& rho : iter_functions_on_center(obs, X.center)) {
          write_as_xml("Z_tilde_bar", Z_tilde_bar[X].row(rho.bfoff_in_atom), attrs<int>{
              {"ao_index1", X},
              {"ao_index2", rho}
          });
        }
      }

    }

    timer.exit();
  }

  memory_used_ += Z_size;
  if(exact_diagonal_K_) {
    ExEnv::out0() << indent
        << "Z intermediate requires " << data_size_to_string(Z_size)
        << std::endl << indent << "Total memory usage is now at least "
        << data_size_to_string(memory_used_)
        << std::endl;
  }

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Make the CADF-LinK lists                                                         {{{1 */ #if 1 //latex `\label{sc:link}`

  std::vector<std::tuple<int, int, int>> L_3_keys;
  std::unordered_map<int, std::unordered_set<int>> L_3_keys_X_to_mu;

  // Now make the linK lists if we're doing linK
  if(do_linK_){
    timer.enter("LinK lists");

    //============================================================================//
    // Get the Frobenius norms of the density matrix shell blocks

    Eigen::MatrixXd D_frob(obs->nshell(), obs->nshell());
    do_threaded(nthread_, [&](int ithr){
      for(int lsh_index = ithr; lsh_index < obs->nshell(); lsh_index += nthread_) {           //latex `\label{sc:link:ld}`
        ShellData lsh(lsh_index, obs, dfbs_);
        for(auto&& jsh : shell_range(obs)) {
          double dnorm = D.block(lsh.bfoff, jsh.bfoff, lsh.nbf, jsh.nbf).norm();
          D_frob(lsh, jsh) = dnorm;
        }
      }
    });
    Eigen::MatrixXd D_frob_sq(obs->nshell(), obs->nshell());
    D_frob_sq = D_frob.array().square();

    //----------------------------------------//                                             //latex `\label{sc:link:setupend}`
    // Form L_DC

    timer.enter("build L_DC");
    // TODO Do this only for my Xsh
    RowMatrix dtilde(gbs_->nshell(), dfbs_->nshell());
    dtilde = D_frob * C_bar_ * schwarz_df_.asDiagonal();
    for(auto&& jsh : shell_range(gbs_)) {
      auto& L_DC_jsh = L_DC[jsh];
      L_DC_jsh.acquire_and_sort(dtilde.data() + dfbs_->nshell() * jsh, dfbs_->nshell(), 0.0);
      L_DC_jsh.set_basis(dfbs_, gbs_);
    }

    //----------------------------------------//                                             //latex `\label{sc:link:ldc:end}`
    // Form L_3

    timer.change("build L_3");                                                               //latex `\label{sc:link:l3}`

    // TODO Scale the B_screening_thresh also

    double epsilon = full_screening_thresh_;
    double epsilon_dist = distance_screening_thresh_;
    double epsilon_B = B_screening_thresh_;

    if(density_reset_){
      prev_density_frob_ = D_frob.norm();
      prev_epsilon_ = epsilon;
      prev_epsilon_dist_ = epsilon_dist;
      prev_epsilon_B_ = epsilon_B;
    }
    else{
      if(scale_screening_thresh_) {
        const double ratio = D_frob.norm() / prev_density_frob_;
        epsilon = prev_epsilon_ * ratio;
        epsilon_dist = prev_epsilon_dist_ * ratio;
        epsilon_B = prev_epsilon_B_ * ratio;
      }

      prev_epsilon_ = epsilon;
      prev_epsilon_dist_ = epsilon_dist;

      if(full_screening_expon_ != 1.0) {
        epsilon = pow(epsilon, full_screening_expon_);                                      //latex `\label{sc:link:expon}`
        epsilon_dist = pow(epsilon_dist, full_screening_expon_);
      }

      epsilon = std::max(full_screening_thresh_min_, epsilon);
      epsilon_dist = std::max(full_screening_thresh_min_, epsilon_dist);
    }


    // Update the user on the effective thresholds
    if(scale_screening_thresh_) {
      if(linK_use_distance_ and epsilon != epsilon_dist) {
        if(screen_B_) {
          ExEnv::out0() << indent << "  Effective LinK screening, LinK distance screening, and B screening thresholds are now "
                        << scprintf("%5.2e, %5.2e, and %5.2e", epsilon, epsilon_dist, epsilon_B) << endl;
        }
        else {
          ExEnv::out0() << indent << "  Effective LinK screening and LinK distance screening thresholds are now "
                        << scprintf("%5.2e and %5.2e", epsilon, epsilon_dist) << endl;
        }
      }
      else {
        if(screen_B_) {
          if(epsilon_B != epsilon) {
            ExEnv::out0() << indent << "  Effective LinK screening and B screening thresholds are now "
                          << scprintf("%5.2e and %5.2e", epsilon, epsilon_B) << endl;
          }
          else {
            ExEnv::out0() << indent << "  Effective LinK screening (and B screening) threshold is now "
                          << scprintf("%5.2e", epsilon) << endl;
          }
        }
        else {
          ExEnv::out0() << indent << "  Effective LinK screening threshold is now "
                        << scprintf("%5.2e", epsilon) << endl;
        }
      }
    }


    // Loop over all jsh
    const auto& const_loc_pairs = local_pairs_linK_;
    const auto& const_loc_pairs_map = linK_local_map_;
    const auto& loc_pairs_end = const_loc_pairs.end();
    do_threaded(nthread_, [&](int ithr) {
      for(auto&& jsh : thread_over_range(shell_range(gbs_, dfbs_), ithr, nthread_)) {
        auto& L_sch_jsh = L_schwarz[jsh];
        // TODO take this as the intersection with local Xsh list
        for(auto&& Xsh : L_DC[jsh]) {

          auto found = const_loc_pairs_map.find(Xsh);
          if(found != const_loc_pairs_map.end()) {
            const double pf = Xsh.value;
            const double eps_prime = epsilon / pf;
            const double eps_prime_dist = epsilon_dist / pf;
            const double Xsh_schwarz = schwarz_df_[Xsh];
            // TODO Figure out a CORRECT way to exit the Xsh loop early in parallel context
            auto local_schwarz_jsh = L_sch_jsh.intersection_with(found->second);
            bool stop_mu = false;
            bool stop_rho = not distribute_coefficients_;
            for(auto&& ish : local_schwarz_jsh) {

              const double dist_factor = get_distance_factor(ish, jsh, Xsh);

              if(!stop_mu) {
                if(ish.value > eps_prime) {
                  if(!linK_use_distance_ or ish.value * dist_factor > eps_prime_dist) {
                    auto& L_3_ish_Xsh = L_3[{ish, Xsh}];
                    L_3_ish_Xsh.insert(jsh,
                        ish.value * dist_factor
                    );
                    if(screen_B_) {
                      double B_scr_contrib = ish.value * ish.value;
                      if(screen_B_use_distance_) B_scr_contrib *= dist_factor * dist_factor;
                      L_3_ish_Xsh.add_to_aux_value_vector(B_scr_contrib, D_frob_sq.col(jsh));
                    }

                    if(print_screening_stats_) {
                      ++iter_stats_->K_3c_needed;
                      iter_stats_->K_3c_needed_fxn += ish.nbf * jsh.nbf * Xsh.nbf;
                    }

                  }
                  else if(print_screening_stats_ and linK_use_distance_) {
                    ++iter_stats_->K_3c_dist_screened;
                    iter_stats_->K_3c_dist_screened_fxn += ish.nbf * jsh.nbf * Xsh.nbf;
                  }

                }
                else {
                  stop_mu = true;
                }
              }

              if(!stop_rho and ish.center != jsh.center) {
                const double ish_value_rho = ish.value * dtilde(ish, Xsh);
                if(ish_value_rho > epsilon) {
                  if(!linK_use_distance_ or ish_value_rho * dist_factor > epsilon_dist) {
                    auto& L_3_ish_Xsh = L_3[{ish, Xsh}];
                    L_3_ish_Xsh.insert_in_aux_set(jsh);
                  }
                }
                // We can't stop with rho because we aren't sorted by dtilde(ish, Xsh)
                // TODO sort and do seperate loop so that rho can exit early
              }

              if(stop_mu and stop_rho) {
                break;
              }

            } // end loop over ish

          }
        }

      }
    });

    L_DC.clear();

    //----------------------------------------//                                             //latex `\label{sc:link:lB}`
    // Form L_B

    timer.change("sort L_3 and build L_B");


    {
      std::mutex L_3_key_merge_mutex;
      do_threaded(nthread_, [&](int ithr){
        std::unordered_map<int, std::set<int>> my_L_3_keys_X_to_mu;

        auto L_3_iter = L_3.begin();
        const auto& L_3_end = L_3.end();
        L_3_iter.advance(ithr);
        while(L_3_iter != L_3_end) {
          //if(not linK_block_rho_) L_3_iter->second.set_sort_by_value(false);
          // Do this even if we're blocking by rho, since it maximizes block sizes
          L_3_iter->second.set_sort_by_value(false);
          L_3_iter->second.sort();
          L_3_iter->second.set_aux_value(sqrt(L_3_iter->second.get_aux_value()));
          if(distribute_coefficients_) {
            L_3_iter->second.sort_aux_set();
          }

          // Build the B screening list

          if(screen_B_) {
            int ish, Xsh;
            ish = L_3_iter->first.first;
            Xsh = L_3_iter->first.second;
            const int Xsh_center = dfbs_->shell_to_center(Xsh);
            auto& L_B_ish_Xsh = L_B[{ish, Xsh}];
            // Note that we could probably avoid a copy here by doing a const cast
            Eigen::VectorXd aux_vect = L_3_iter->second.get_aux_vector();
            const auto& L_3_ish_Xsh = L_3_iter->second;
            const auto& aux_val = L_3_ish_Xsh.get_aux_value();
            aux_vect = aux_vect.array().cwiseSqrt() * aux_val * schwarz_df_[Xsh];

            // TODO we can further restrict this loop by prescreening it
            for(auto&& lsh : shell_range(obs)) {
              if(lsh.center == Xsh_center) continue;
              if(fabs(C_bar_(lsh, Xsh) * aux_vect(lsh)) > epsilon_B) {
                L_B_ish_Xsh.insert(lsh);
                L_B_ish_Xsh.set_aux_value(lsh.nbf + L_B_ish_Xsh.get_aux_value());
              }
            }

            L_B_ish_Xsh.set_sort_by_value(false);
            L_B_ish_Xsh.sort();
            if(distribute_coefficients_) {
              throw FeatureNotImplemented("distribute coefficients with B screening", __FILE__, __LINE__, class_desc());
            }
            if(L_B_ish_Xsh.size() > 0) {
              my_L_3_keys_X_to_mu[Xsh].insert(ish);
            }
          }
          else {
            my_L_3_keys_X_to_mu[L_3_iter->first.second].insert(L_3_iter->first.first);
          }

          L_3_iter.advance(nthread_);
        }

        //----------------------------------------//
        // Merge in our part of the L_3 keys map
        std::lock_guard<std::mutex> lg(L_3_key_merge_mutex);
        for(auto&& kvpair : my_L_3_keys_X_to_mu) {
          L_3_keys_X_to_mu[kvpair.first].insert(
              kvpair.second.begin(), kvpair.second.end()
          );
        }
        //----------------------------------------//
      });                                                                                      //latex `\label{sc:link:l3:end}`
    }


    for(auto&& pair : L_3) {
      if(screen_B_) {
        const auto& L_B_part = L_B[pair.first];
        if(L_B_part.size() == 0) continue;
      }
      L_3_keys.emplace_back(
          pair.first.first,
          pair.first.second,
          pair.second.size()
      );
    }

    std::sort(L_3_keys.begin(), L_3_keys.end(), [](std::tuple<int, int, int> const& A, std::tuple<int, int, int> const& B){
      int ia, ib, Xa, Xb, sza, szb;
      std::tie(ia, Xa, sza) = A;
      std::tie(ib, Xb, szb) = B;
      if(sza > szb) return true;
      else if(sza == szb) {
        if(ia < ib) return true;
        else if(ia == ib) {
          if(Xa < Xb) return true;
        }
      }
      return false;
    });


    timer.exit();
  }

  timer.exit("LinK lists");


  /*****************************************************************************************/ #endif //1}}} //latex `\label{sc:link:end}`
  /*=======================================================================================*/
  /* Loop over local shell pairs for three body contributions                         {{{1 */ #if 1 //latex `\label{sc:k3b:begin}`

  {

    Timer timer("three body contributions");
    boost::mutex tmp_mutex;
    std::mutex L_3_mutex;

    auto L_3_key_iter = L_3_keys.begin();
    //auto L_3_prime_key_iter = L_3_prime.begin();
    boost::thread_group compute_threads;

    // reset the iteration over local pairs
    local_pairs_spot_ = 0;

    MultiThreadTimer mt_timer("threaded part", nthread_);

    auto get_ish_Xblk_3 = [&](ShellData& ish, ShellBlockData<>& Xblk) -> bool {                                                //latex `\label{sc:k3b:getiX}`
      if(do_linK_) {
        std::lock_guard<std::mutex> lg(L_3_mutex);
        if(L_3_key_iter == L_3_keys.end()) {
          return false;
        }
        else {

          int ishidx, Xshidx, list_size;
          std::tie(ishidx, Xshidx, list_size) = *L_3_key_iter;
          ish = ShellData(ishidx, gbs_, dfbs_);
          Xblk = ShellBlockData<>(ShellData(Xshidx, dfbs_, gbs_));
          ++L_3_key_iter;
          return true;
        }
      }
      else {
        int spot = local_pairs_spot_++;
        if(spot < local_pairs_k_.size()){
          auto& sig_pair = local_pairs_k_[spot];
          ish = ShellData(sig_pair.first, gbs_, dfbs_);
          Xblk = sig_pair.second;
          return true;
        }
        else{
          return false;
        }
      }
    };                                                                                   //latex `\label{sc:k3b:getiXend}`

    // Loop over number of threads
    for(int ithr = 0; ithr < nthread_; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){                                              //latex `\label{sc:k3b:thrpatend}`

        // thread-local portion of K_tilde
        Eigen::MatrixXd Kt_part(nbf, nbf);
        Kt_part = Eigen::MatrixXd::Zero(nbf, nbf);

        ShellData ish;
        ShellData jsh;
        ShellBlockData<> Xblk;

        // Allocate buffers beforehand so that we don't have to do so inside the main loop
        long b_buff_offset = 0;
        size_t actual_size = 0;
        actual_size = B_buffer_size_/sizeof(double);
        double* __restrict__ b_buffer = new double[actual_size];
        double* __restrict__ B_ish_data = new double[max_fxn_obs_todo_*max_fxn_dfbs_todo_*nbf];
        double* __restrict__ D_B_buff_data;
        double* __restrict__ D_sd_data;
        double* __restrict__ C_X_diff_data;
        double* __restrict__ dt_ish_X_data;
        double* __restrict__ dt_prime_data;
        double* __restrict__ gt_ish_X_data;
        double* __restrict__ C_X_block_data;
        double* __restrict__ D_ordered_data;
        if(distribute_coefficients_) {
           dt_ish_X_data = new double[max_fxn_obs_todo_*max_fxn_dfbs_todo_*nbf];
           dt_prime_data = new double[max_fxn_obs_todo_*max_fxn_dfbs_todo_*max_obs_atom_fxn_on_dfbs_center_todo_];
           gt_ish_X_data = new double[max_fxn_obs_todo_*max_fxn_dfbs_todo_*nbf];
        }
        int D_B_buff_n = 0;
        if(B_use_buffer_) {
          D_B_buff_data = new double[nbf*nbf];
          D_B_buff_n = nbf;
        }
        if(linK_block_rho_) {
          D_ordered_data = new double[nbf*nbf];
        }
        if(screen_B_) {
          D_sd_data = new double[nbf*nbf];
          C_X_diff_data = new double[max_obs_atom_fxn_on_dfbs_center_todo_*max_fxn_dfbs_todo_*nbf];
          if(screen_B_transfer_as_transpose_) {
            C_X_block_data = new double[max_obs_atom_fxn_on_dfbs_center_todo_*max_fxn_dfbs_todo_*nbf];
          }
        }
        Eigen::Map<ColMatrix> D_B_buff(D_B_buff_data, D_B_buff_n, D_B_buff_n);
        Eigen::Map<ColMatrix> D_sd(D_sd_data, 0, 0);
        Eigen::Map<RowMatrix> C_X_diff(C_X_diff_data, 0, 0);
        Eigen::Map<RowMatrix> B_sd(B_ish_data, 0, 0);
        Eigen::Map<RowMatrix> dt_ish_X(dt_ish_X_data, 0, 0);
        Eigen::Map<RowMatrix> dt_prime(dt_prime_data, 0, 0);
        Eigen::Map<RowMatrix> gt_ish_X(gt_ish_X_data, 0, 0);
        Eigen::Map<ColMatrix> D_ordered(D_ordered_data, 0, 0);

        //============================================================================//
        // Main loop
        //============================================================================//

        while(get_ish_Xblk_3(ish, Xblk)) {                                                            //latex `\label{sc:k3b:while}`

          /*-----------------------------------------------------*/
          /* Compute B intermediate                         {{{2 */ #if 2 // begin fold      //latex `\label{sc:k3b:b}`

          // Timer stuff
          mt_timer.enter("compute B", ithr);
          auto ints_timer = mt_timer.get_subtimer("compute ints", ithr);
          auto k2_part_timer = mt_timer.get_subtimer("k2 part", ithr);
          auto contract_timer = mt_timer.get_subtimer("contract", ithr);
          auto ex_timer = mt_timer.get_subtimer("exact diagonal", ithr);

          // Because of the way the L_3 lists are made, we can only do one X shell
          //   at a time if we're doing linK.  Otherwise, we want to block by atoms.
          assert(!do_linK_ or Xblk.nshell == 1);
          auto&& Xsh = Xblk.first_shell;


          // What list of J are we using?
          OrderedShellList* jlist;
          if(not do_linK_) {
            jlist = new OrderedShellList(sig_partners_[ish], gbs_, dfbs_);
          }
          else {
            jlist = &(L_3[{ish, Xsh}]);
          }
          const int jlist_size = jlist->size();


          // Reordering of D, only if linK_block_rho_ is true
          int block_offset = 0;
          if(do_linK_ and linK_block_rho_) {

            mt_timer.enter("rearrange D", ithr);
            new (&D_ordered) Eigen::Map<ColMatrix>(D_ordered_data, nbf, jlist_size);
            D_ordered.resize(nbf, nbf);
            for(auto jblk : shell_block_range(*jlist, Contiguous)){
              D_ordered.middleCols(block_offset, jblk.nbf) = D.middleCols(jblk.bfoff, jblk.nbf);
              block_offset += jblk.nbf;
            }
            mt_timer.exit(ithr);
          }

          /*-----------------------------------------------------*/
          /* Initialize memory for B screening if screen_B_ {{{2 */ #if 2 // begin fold

          // Form the B_sd, D_sd, and C_X matrices
          int l_b_size;

          if(screen_B_) {

            mt_timer.enter("build screened B part", ithr);

            const auto& L_B_ish_Xsh = L_B[{ish, Xsh}];
            l_b_size = int(L_B_ish_Xsh.get_aux_value());
            if(l_b_size == 0){
              mt_timer.exit(ithr); // build screened B part
              mt_timer.exit(ithr); // compute B
              continue;
            }

            new (&B_sd) Eigen::Map<RowMatrix>(B_ish_data, ish.nbf*Xblk.nbf, Xblk.atom_obsnbf + l_b_size);

            if(linK_block_rho_) {
              new (&D_sd) Eigen::Map<ColMatrix>(D_sd_data, Xblk.atom_obsnbf + l_b_size, jlist_size);
              D_sd.topRows(Xblk.atom_obsnbf) = D_ordered.middleRows(Xblk.atom_obsbfoff, Xblk.atom_obsnbf);
            }
            else {
              new (&D_sd) Eigen::Map<ColMatrix>(D_sd_data, Xblk.atom_obsnbf + l_b_size, nbf);
              D_sd.topRows(Xblk.atom_obsnbf) = D.middleCols(Xblk.atom_obsbfoff, Xblk.atom_obsnbf).transpose();
            }

            if(screen_B_transfer_as_transpose_) {

              Eigen::Map<ColMatrix> C_X_diff_cols(C_X_diff_data, Xblk.nbf*Xblk.atom_obsnbf, l_b_size);

              block_offset = 0;
              Eigen::Map<ColMatrix> C_X_block(C_X_block_data, Xblk.nbf*Xblk.atom_obsnbf, nbf);
              C_X_block = coefs_transpose_blocked_other_[Xblk.center].middleRows(
                    Xblk.bfoff_in_atom*Xblk.atom_obsnbf, Xblk.nbf*Xblk.atom_obsnbf
              );
              for(auto&& lblk : shell_block_range(L_B_ish_Xsh, Contiguous)) {
                C_X_diff_cols.middleCols(block_offset, lblk.nbf) = C_X_block.middleCols(lblk.bfoff, lblk.nbf);
                if(linK_block_rho_) {
                  D_sd.middleRows(Xblk.atom_obsnbf + block_offset, lblk.nbf) = D_ordered.middleRows(lblk.bfoff, lblk.nbf);
                }
                else {
                  D_sd.middleRows(Xblk.atom_obsnbf + block_offset, lblk.nbf) = D.middleCols(lblk.bfoff, lblk.nbf).transpose();
                }
                block_offset += lblk.nbf;
              }

              new (&C_X_diff) Eigen::Map<RowMatrix>(C_X_block_data, Xblk.nbf*Xblk.atom_obsnbf, l_b_size);
              C_X_diff = C_X_diff_cols;

            }
            else {

              new (&C_X_diff) Eigen::Map<RowMatrix>(C_X_diff_data, Xblk.nbf*Xblk.atom_obsnbf, l_b_size);

              block_offset = 0;
              const auto& C_X_block = coefs_transpose_blocked_other_[Xblk.center].middleRows(
                    Xblk.bfoff_in_atom*Xblk.atom_obsnbf, Xblk.nbf*Xblk.atom_obsnbf
              );
              for(auto&& lblk : shell_block_range(L_B_ish_Xsh, Contiguous)) {
                C_X_diff.middleCols(block_offset, lblk.nbf) = C_X_block.middleCols(lblk.bfoff, lblk.nbf);
                if(linK_block_rho_) {
                  D_sd.middleRows(Xblk.atom_obsnbf + block_offset, lblk.nbf) = D_ordered.middleRows(lblk.bfoff, lblk.nbf);
                }
                else {
                  D_sd.middleRows(Xblk.atom_obsnbf + block_offset, lblk.nbf) = D.middleCols(lblk.bfoff, lblk.nbf).transpose();
                }
                block_offset += lblk.nbf;
              }

            }

            mt_timer.exit(ithr);
          }

          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/

          // Create B_ish and the B buffer
          Eigen::Map<RowMatrix> B_ish(B_ish_data, ish.nbf * Xblk.nbf, nbf);

          // Set the portion of the B_ish_data buffer we're going to use to zero
          if(not screen_B_) std::memset(B_ish_data, 0, ish.nbf*Xblk.nbf*nbf*sizeof(double));
          else std::memset(B_ish_data, 0, ish.nbf * Xblk.nbf * (Xblk.atom_obsnbf + l_b_size) * sizeof(double));

          // We can also buffer the B mat to do fewer contractions, even if
          //   we're not screening B.  This doesn't really seem to help things...
          int b_buff_nrows, b_buff_ncols;
          if(B_use_buffer_) {
            b_buff_nrows = std::min(size_t(nbf), actual_size / ish.nbf / Xblk.nbf);
            b_buff_ncols = ish.nbf * Xblk.nbf;
          }
          else {
            b_buff_ncols = b_buff_nrows = 0;
          }
          Eigen::Map<RowMatrix> B_buffer_mat(b_buffer, b_buff_nrows, b_buff_ncols);
          b_buff_offset = 0;

          // Exact diagonal itermediate storage
          RowMatrix M_mu_X, W_mu_X, W_mu_X_bar;
          if(exact_diagonal_K_) {
            M_mu_X.resize(ish.nbf * Xblk.nbf, ish.atom_nbf);
            M_mu_X = RowMatrix::Zero(ish.nbf*Xblk.nbf, ish.atom_nbf);
            W_mu_X.resize(ish.nbf * Xblk.nbf, nbf);
            W_mu_X = RowMatrix::Zero(ish.nbf*Xblk.nbf, nbf);
            W_mu_X_bar.resize(ish.nbf * Xblk.nbf, nbf);
            W_mu_X_bar = RowMatrix::Zero(ish.nbf*Xblk.nbf, nbf);
          }

          // TODO figure out how to take advantage of L_3 sorting


          /*-----------------------------------------------------*/
          /* Transpose part if distribute_coefficients_     {{{2 */ #if 2 // begin fold
          // Build dt for the transpose, only if we're distributing coefficients
          if(distribute_coefficients_) {
            if(do_linK_) {
              mt_timer.change("transpose part", ithr);

              mt_timer.enter("compute d_tilde", ithr);
              new (&dt_ish_X) Eigen::Map<RowMatrix>(dt_ish_X_data, ish.nbf * Xsh.nbf, nbf);
              dt_ish_X  = RowMatrix::Zero(ish.nbf * Xsh.nbf, nbf);
              new (&dt_prime) Eigen::Map<RowMatrix>(dt_prime_data, ish.nbf * Xsh.nbf, Xsh.atom_obsnbf);
              dt_prime  = RowMatrix::Zero(ish.nbf * Xsh.nbf, Xsh.atom_obsnbf);
              {
                Eigen::Map<RowMatrix, Eigen::Default, Eigen::OuterStride<Eigen::Dynamic>> dt_ish_part(
                    NULL, 0, 0, Eigen::OuterStride<Eigen::Dynamic>(Xsh.nbf*nbf)
                );
                for(auto&& X : function_range(Xsh)) {
                  new (&dt_ish_part) Eigen::Map<RowMatrix, Eigen::Unaligned, Eigen::OuterStride<Eigen::Dynamic>>(
                      dt_ish_X.data() + X.off*nbf, ish.nbf, nbf, Eigen::OuterStride<Eigen::Dynamic>(Xsh.nbf*nbf)
                  );
                  dt_ish_part += 2.0 *
                      D.middleRows(ish.bfoff, ish.nbf).middleCols(Xsh.atom_obsbfoff, Xsh.atom_obsnbf)
                      * coefs_X_nu_other.at(Xsh.center).middleRows(
                          X.bfoff_in_atom*Xsh.atom_obsnbf, Xsh.atom_obsnbf
                      );
                }
              }
              {
                Eigen::Map<RowMatrix, Eigen::Default, Eigen::OuterStride<Eigen::Dynamic>> dt_p_part(
                    NULL, 0, 0, Eigen::OuterStride<Eigen::Dynamic>(Xsh.nbf*nbf)
                );
                for(auto&& X : function_range(Xsh)) {
                  new (&dt_p_part) Eigen::Map<RowMatrix, Eigen::Unaligned, Eigen::OuterStride<Eigen::Dynamic>>(
                      dt_prime.data() + X.off*Xsh.atom_obsnbf, ish.nbf, Xsh.atom_obsnbf,
                      Eigen::OuterStride<Eigen::Dynamic>(Xsh.nbf*Xsh.atom_obsnbf)
                  );
                  dt_p_part += 2.0 *
                      D.middleRows(ish.bfoff, ish.nbf)
                      * coefs_X_nu_other.at(Xsh.center).middleRows(
                          X.bfoff_in_atom*Xsh.atom_obsnbf, Xsh.atom_obsnbf
                      ).transpose();
                }
              }

              mt_timer.change("form g and K contributions", ithr);
              auto form_g_timer = mt_timer.get_subtimer("form g", ithr);
              auto k_contrib_timer = mt_timer.get_subtimer("K contrib", ithr);

              /// Now loop over rho in L_3 auxiliary list and compute k contributions
              for(const auto&& jblk : shell_block_range(
                  L_3[{ish, Xsh}].get_aux_indices(), gbs_.pointer(), dfbs_.pointer(), Contiguous
              )){

                TimerHolder holder(form_g_timer);
                new (&gt_ish_X) Eigen::Map<RowMatrix>(gt_ish_X_data, ish.nbf*Xblk.nbf, jblk.nbf);
                gt_ish_X = RowMatrix::Zero(ish.nbf*Xblk.nbf, jblk.nbf);
                for(auto&& mu : function_range(ish)) {
                  gt_ish_X.middleRows(mu.off*Xblk.nbf, Xblk.nbf) -= 0.5 *
                      g2.middleCols(ish.atom_dfbfoff, ish.atom_dfnbf).middleRows(Xblk.bfoff, Xblk.nbf)
                      * coefs_mu_X_other.at(ish).middleRows(mu.off*nbf + jblk.bfoff, jblk.nbf).transpose();
                }


                holder.change(k_contrib_timer);
                for(auto&& mu : function_range(ish)) {
                  Kt_part.middleCols(jblk.bfoff, jblk.nbf).transpose().noalias() +=
                      gt_ish_X.middleRows(mu.off*Xblk.nbf, Xblk.nbf).transpose()
                      * dt_ish_X.middleRows(mu.off*Xblk.nbf, Xblk.nbf);
                  Kt_part.middleCols(jblk.bfoff, jblk.nbf).middleRows(Xblk.atom_obsbfoff, Xblk.atom_obsnbf).transpose().noalias() +=
                      gt_ish_X.middleRows(mu.off*Xblk.nbf, Xblk.nbf).transpose()
                      * dt_prime.middleRows(mu.off*Xblk.nbf, Xblk.nbf);
                  Kt_part.middleCols(jblk.bfoff, jblk.nbf).middleRows(Xblk.atom_obsbfoff, Xblk.atom_obsnbf).transpose().noalias() -=
                      gt_ish_X.middleRows(mu.off*Xblk.nbf, Xblk.nbf).transpose()
                      * dt_ish_X.middleRows(mu.off*Xsh.nbf, Xblk.nbf).middleCols(Xblk.atom_obsbfoff, Xblk.atom_obsnbf);
                }

              }

              mt_timer.exit(ithr);

            }
            else {
              throw FeatureNotImplemented("non-linK", __FILE__, __LINE__, class_desc());
            }
          }
          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/


          mt_timer.change("compute B", ithr);

          //============================================================================//
          // Loop over the largest blocks of J at once that we can

          int restrictions = linK_block_rho_ ? NoRestrictions : Contiguous;
          // Reset the block offset for the jblk loop
          block_offset = 0;
          for(const auto&& jblk : shell_block_range(*jlist, restrictions)){             //latex `\label{sc:k3b:noblk:loop}`
            TimerHolder subtimer(ints_timer);

            auto g3_in = ints_to_eigen_map(
                jblk, ish, Xblk,
                eris_3c_[ithr], coulomb_oper_type_,
                b_buffer + (B_use_buffer_ ? (b_buff_offset * ish.nbf * Xblk.nbf) : 0)
            );

            //----------------------------------------//

            // Now view the integrals as a jblk.nbf x (ish.nbf*Xsh.nbf) matrix, which makes
            //   the contraction more convenient.  Doesn't require any movement of data

            Eigen::Map<ThreeCenterIntContainer> g3(g3_in.data(), jblk.nbf, ish.nbf*Xblk.nbf);

            //----------------------------------------//
            /* Two-body part                     {{{3 */ #if 3 // begin fold

            // TODO This breaks integral caching (if I ever use it again)

            subtimer.change(k2_part_timer);

            int subblock_offset = 0;
            for(const auto&& jsblk : shell_block_range(jblk, Contiguous|SameCenter)) {
              int inner_size = ish.atom_dfnbf;
              if(ish.center != jsblk.center) {
                inner_size += jsblk.atom_dfnbf;
              }

              if(distribute_coefficients_) {

                for(int mu_off = 0; mu_off < ish.nbf; ++mu_off) {
                  g3.middleCols(mu_off*Xblk.nbf, Xblk.nbf).middleRows(subblock_offset, jsblk.nbf) -= 0.5 *
                      coefs_mu_X_other.at(ish).middleRows(mu_off*nbf + jsblk.bfoff, jsblk.nbf)
                      * g2.middleRows(ish.atom_dfbfoff, ish.atom_dfnbf).middleCols(Xblk.bfoff, Xblk.nbf);
                }


              }
              else {
                if(store_coefs_transpose_) {

                  const int tot_cols = coefs_blocked_[jsblk.center].cols();
                  const int col_offset = coef_block_offsets_[jsblk.center][ish.center]
                      + ish.bfoff_in_atom*inner_size;
                  double* data_start = coefs_blocked_[jsblk.center].data() +
                      jsblk.bfoff_in_atom * tot_cols + col_offset;

                  StridedRowMap Ctmp(data_start, jsblk.nbf, ish.nbf*inner_size,
                      Eigen::OuterStride<>(tot_cols)
                  );

                  RowMatrix C(Ctmp.nestByValue());
                  C.resize(jsblk.nbf*ish.nbf, inner_size);

                  g3_in.middleRows(subblock_offset*ish.nbf, jsblk.nbf*ish.nbf) -= 0.5
                      * C.rightCols(ish.atom_dfnbf) * g2.block(
                          ish.atom_dfbfoff, Xblk.bfoff,
                          ish.atom_dfnbf, Xblk.nbf
                  );
                  if(ish.center != jsblk.center) {
                    g3_in.middleRows(subblock_offset*ish.nbf, jsblk.nbf*ish.nbf) -= 0.5
                        * C.leftCols(jsblk.atom_dfnbf) * g2.block(
                            jsblk.atom_dfbfoff, Xblk.bfoff,
                            jsblk.atom_dfnbf, Xblk.nbf
                    );
                  }

                }
                else {

                  int rho_block_offset = 0;
                  for(auto&& mu : function_range(ish)) {
                    g3.middleRows(subblock_offset, jsblk.nbf).middleCols(mu.off*Xblk.nbf, Xblk.nbf) -= 0.5 *
                        coefs_transpose_blocked_[ish.center].middleCols(mu.bfoff_in_atom*nbf + jsblk.bfoff, jsblk.nbf).transpose()
                          * g2.middleRows(ish.atom_dfbfoff, ish.atom_dfnbf).middleCols(Xblk.bfoff, Xblk.nbf);
                  }

                  if(ish.center != jsblk.center) {
                    int rho_block_offset = 0;
                    for(auto&& rho : function_range(gbs_, dfbs_, jsblk.bfoff, jsblk.last_function)) {
                      g3_in.middleRows(
                          (subblock_offset + rho_block_offset)*ish.nbf, ish.nbf
                      ) -= 0.5 *
                          coefs_transpose_blocked_[jsblk.center].middleCols(
                              rho.bfoff_in_atom*nbf + ish.bfoff, ish.nbf
                          ).transpose() * g2.middleRows(
                              jsblk.atom_dfbfoff, jsblk.atom_dfnbf
                          ).middleCols(
                                  Xblk.bfoff, Xblk.nbf
                          );
                      ++rho_block_offset;
                    }
                  }

                }
              }

              if(exact_diagonal_K_) {
                if(distribute_coefficients_) {
                  throw FeatureNotImplemented("exact diagonal with distributed coefficients", __FILE__, __LINE__, class_desc());
                }
                subtimer.change(ex_timer);

                if(jsblk.center == Xblk.center or ish.center == Xblk.center) {
                  // Build W and Wbar
                  for(auto&& mu : function_range(ish)) {
                    for(auto&& Y : iter_functions_on_center(dfbs_, ish.center)) {
                      // TODO remove one loop here
                      W_mu_X.middleRows(mu.off*Xblk.nbf, Xblk.nbf).middleCols(jsblk.bfoff, jsblk.nbf) +=
                          g2.col(Y).segment(Xblk.bfoff, Xblk.nbf)
                          * coefs_transpose_[Y].row(mu.bfoff_in_atom).segment(jsblk.bfoff, jsblk.nbf);
                    }
                    if(ish.center != jsblk.center){
                      for(auto&& Y : iter_functions_on_center(dfbs_, jsblk.center)) {
                        W_mu_X_bar.middleRows(mu.off*Xblk.nbf, Xblk.nbf).middleCols(jsblk.bfoff, jsblk.nbf) +=
                            g2.col(Y).segment(Xblk.bfoff, Xblk.nbf)
                            * coefs_transpose_[Y].col(mu).segment(jsblk.bfoff-jsblk.atom_bfoff, jsblk.nbf).transpose();
                      }
                    }
                  }
                }

                // Build M
                if(jsblk.center == Xblk.center) {
                  M_mu_X += (
                      4.0 * g3.middleRows(subblock_offset, jsblk.nbf).transpose()
                      - W_mu_X.middleCols(jsblk.bfoff, jsblk.nbf)
                      - W_mu_X_bar.middleCols(jsblk.bfoff, jsblk.nbf)
                      ) * D.middleCols(ish.atom_bfoff, ish.atom_nbf).middleRows(jsblk.bfoff, jsblk.nbf);
                }

                if(Xblk.center == ish.center) {
                  for(auto&& mu : function_range(ish)) {
                    for(auto&& nu : iter_functions_on_center(obs, jsblk.center)) {
                      Kt_part(mu, nu) -= (Z[nu.center].middleCols(Xblk.bfoff, Xblk.nbf).middleRows(
                          nu.bfoff_in_atom*jsblk.atom_nbf + jsblk.bfoff_in_atom, jsblk.nbf
                        ).array() * (
                             2.0 * g3.middleCols(mu.off*Xblk.nbf, Xblk.nbf).middleRows(subblock_offset, jsblk.nbf)
                              - 0.5 * W_mu_X.middleRows(mu.off*Xblk.nbf, Xblk.nbf).middleCols(jsblk.bfoff, jsblk.nbf).transpose()
                      ).array()).sum();
                      if(ish.center != jsblk.center) {
                        Kt_part(mu, nu) += (Z[nu.center].middleCols(Xblk.bfoff, Xblk.nbf).middleRows(
                            nu.bfoff_in_atom*jsblk.atom_nbf + jsblk.bfoff_in_atom, jsblk.nbf
                          ).array()
                          * 0.5 * W_mu_X_bar.middleRows(mu.off*Xblk.nbf, Xblk.nbf).middleCols(jsblk.bfoff, jsblk.nbf).transpose().array()
                        ).sum();

                      }
                    }
                  }
                }

                if(Xblk.center == ish.center) {
                  if(ish.center != jsblk.center) {
                    for(auto&& X : function_range(Xblk)) {
                      for(auto&& mu : function_range(ish)) {
                        Kt_part.row(mu).segment(ish.atom_bfoff, ish.atom_nbf).transpose() -=
                            Z_tilde[X].middleCols(jsblk.bfoff, jsblk.nbf) * (
                             2.0 * g3.middleRows(subblock_offset, jsblk.nbf).col(mu.off*Xblk.nbf + X-Xblk.bfoff)
                             - 0.5 * W_mu_X.row(mu.off*Xblk.nbf + X-Xblk.bfoff).segment(jsblk.bfoff, jsblk.nbf).transpose()
                             - 0.5 * W_mu_X_bar.row(mu.off*Xblk.nbf + X-Xblk.bfoff).segment(jsblk.bfoff, jsblk.nbf).transpose()
                        );
                      }
                    }
                  }
                }

                if(ish.center != Xblk.center) {
                  if(Xblk.center == jsblk.center) {
                    for(auto&& X : function_range(Xblk)) {
                      for(auto&& mu : function_range(ish)) {
                        Kt_part.row(mu).segment(ish.atom_bfoff, ish.atom_nbf).transpose() -=
                            Z_tilde_bar[X].middleRows(jsblk.bfoff_in_atom, jsblk.nbf).middleCols(ish.atom_bfoff, ish.atom_nbf).transpose() * (
                             2.0 * g3.middleRows(subblock_offset, jsblk.nbf).col(mu.off*Xblk.nbf + X-Xblk.bfoff)
                             - 0.5 * W_mu_X.row(mu.off*Xblk.nbf + X-Xblk.bfoff).segment(jsblk.bfoff, jsblk.nbf).transpose()
                             - 0.5 * W_mu_X_bar.row(mu.off*Xblk.nbf + X-Xblk.bfoff).segment(jsblk.bfoff, jsblk.nbf).transpose()
                        );
                      }
                    }

                  }
                }

                subtimer.change(k2_part_timer);

              } // end if exact_diagonal_k

              //----------------------------------------//

              subblock_offset += jsblk.nbf;

            }

            /******************************************/ #endif //3}}}
            //----------------------------------------//

            //----------------------------------------//
            /* Screening stats                   {{{3 */ #if 3 // begin fold
            if(do_linK_ and print_screening_stats_ > 2) {
              mt_timer.enter("count underestimated ints", ithr);

              double epsilon;
              if(density_reset_){ epsilon = full_screening_thresh_; }
              else{ epsilon = pow(full_screening_thresh_, full_screening_expon_); }

              int offset_in_block = 0;
              for(const auto&& jsh : shell_range(jblk)) {
                const double g3_norm = g3.middleRows(offset_in_block, jsh.nbf).norm();
                offset_in_block += jsh.nbf;
                const int nfxn = jsh.nbf*ish.nbf*Xblk.nbf;
                if(L_3[{ish, Xsh}].value_for_index(jsh) < g3_norm) {
                  ++iter_stats_->K_3c_underestimated;
                  iter_stats_->K_3c_underestimated_fxn += jsh.nbf*ish.nbf*Xsh.nbf;
                }
                if(g3_norm * L_DC[jsh].value_for_index(Xsh) > epsilon) {
                  ++iter_stats_->K_3c_perfect;
                  iter_stats_->K_3c_perfect_fxn += nfxn;
                }
                if(xml_screening_data_ and iter_stats_->is_first) {
                  iter_stats_->int_screening_values.mine(ithr).push_back(
                      L_3[{ish, Xsh}].value_for_index(jsh)
                  );
                  iter_stats_->int_actual_values.mine(ithr).push_back(g3_norm);
                  iter_stats_->int_distance_factors.mine(ithr).push_back(
                      get_distance_factor(ish, jsh, Xsh)
                  );
                  iter_stats_->int_distances.mine(ithr).push_back(
                      get_R(ish, jsh, Xsh)
                  );
                  iter_stats_->int_indices.mine(ithr).push_back(
                      std::make_tuple(ish, jsh, Xsh)
                  );
                  iter_stats_->int_ams.mine(ithr).push_back(
                      std::make_tuple(ish.am, jsh.am, Xsh.am)
                  );
                }
              }
              mt_timer.exit(ithr);
            }

            /******************************************/ #endif //3}}}
            //----------------------------------------//

            subtimer.change(contract_timer);

            if(B_use_buffer_) {

              if(b_buff_offset + jblk.nbf > b_buff_nrows) {
                if(b_buff_offset == 0) {
                  throw SCException("B_buffer_size smaller than single contiguous block.  Set B_use_buffer to no and try again.");
                }
                if(screen_B_) {
                  B_sd.noalias() += 2.0 * B_buffer_mat.topRows(b_buff_offset).transpose() * D_B_buff.leftCols(b_buff_offset).transpose();
                }
                else {
                  B_ish.noalias() += 2.0 * B_buffer_mat.topRows(b_buff_offset).transpose() * D_B_buff.leftCols(b_buff_offset).transpose();
                }
                b_buff_offset = 0;
              }

              if(screen_B_) {
                throw FeatureNotImplemented("screen_B with B buffer", __FILE__, __LINE__, class_desc());
              }
              else if (linK_block_rho_) {
                throw FeatureNotImplemented("linK_block_rho_ with B buffer", __FILE__, __LINE__, class_desc());
              }
              else {
                std::copy(D_data + jblk.bfoff, D_data + jblk.bfoff + jblk.nbf, D_B_buff_data + b_buff_offset);
                b_buff_offset += jblk.nbf;
              }

              b_buff_offset += jblk.nbf;

            }
            else {

              if(linK_block_rho_) {
                if(screen_B_) {
                  B_sd.noalias() += 2.0 * g3.transpose() * D_sd.middleCols(block_offset, jblk.nbf).transpose();
                }
                else {
                  B_ish.noalias() += 2.0 * g3.transpose() * D_ordered.middleCols(block_offset, jblk.nbf).transpose();
                }
              }
              else {
                if(screen_B_) {
                  B_sd.noalias() += 2.0 * g3.transpose() * D_sd.middleCols(jblk.bfoff, jblk.nbf).transpose();
                }
                else {
                  B_ish.noalias() += 2.0 * g3.transpose() * D.middleCols(jblk.bfoff, jblk.nbf).transpose();
                }
              }

            }

            //#if !CADF_USE_BLAS
            //#else
            //// BLAS version
            //const char notrans = 'n', trans = 't';
            //const blasint M = ish.nbf*Xblk.nbf;
            //const blasint K = jblk.nbf;
            //const double alpha = 2.0, one = 1.0;
            //double* D_data;
            //if(linK_block_rho_) {
            //  D_data = D_ordered.data() + block_offset;
            //}
            //else {
            //  D_data = D.data() + jblk.bfoff;
            //}
            //F77_DGEMM(&notrans, &notrans,
            //    &M, &nbf, &K,
            //    &alpha, g3.data(), &M,
            //    D_data, &nbf,
            //    &one, B_ish.data(), &M
            //);
            //#endif

            block_offset += jblk.nbf;

          } // end loop over jsh

          if(not do_linK_) {
            delete jlist;
          }

          if(B_use_buffer_ and b_buff_offset > 0) {
            TimerHolder subtimer(contract_timer);
            B_ish.noalias() += 2.0 * B_buffer_mat.topRows(b_buff_offset).transpose() * D_B_buff.leftCols(b_buff_offset).transpose();
            b_buff_offset = 0;
          }

          if(xml_debug_ and exact_diagonal_K_) {
            Eigen::MatrixXd Mout(ish.nbf*Xblk.nbf, nbf);
            Mout = Eigen::MatrixXd::Zero(ish.nbf*Xblk.nbf, nbf);
            Mout.middleCols(ish.atom_bfoff, ish.atom_nbf) = M_mu_X;
            for(auto&& mu : function_range(ish)) {
              for(auto&& X : function_range(Xblk)) {
                write_as_xml("M", Mout.row(mu.off*Xblk.nbf + X-Xblk.bfoff), attrs<int>{
                  {"ao_index1", mu},
                  {"ao_index2", X},
                  {"partial", 1}
                });
                write_as_xml("W_k", W_mu_X.row(mu.off*Xblk.nbf + X-Xblk.bfoff), attrs<int>{
                  {"ao_index1", mu},
                  {"ao_index2", X}
                });
                write_as_xml("Wbar", W_mu_X_bar.row(mu.off*Xblk.nbf + X-Xblk.bfoff), attrs<int>{
                  {"ao_index1", mu},
                  {"ao_index2", X}
                });
              }
            }

          }

          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/

          /*-----------------------------------------------------*/
          /* Compute K contributions                        {{{2 */ #if 2 // begin fold      //latex `\label{sc:k3b:kcontrib}`
          mt_timer.change("K contributions", ithr);
          const int obs_atom_bfoff = obs->shell_to_function(obs->shell_on_center(Xblk.center, 0));
          const int obs_atom_nbf = obs->nbasis_on_center(Xblk.center);
          for(auto&& X : function_range(Xblk)) {
            for(auto&& mu : function_range(ish)) {
              // TODO get rid of one of these loops

              // B_mus[mu.bfoff_in_shell] is (nbf x Ysh.nbf)
              // C_Y is (Y.{obs_}atom_nbf x nbf)
              // result should be (Y.{obs_}atom_nbf x 1)

              if(distribute_coefficients_) {
                if(screen_B_) {
                  throw FeatureNotImplemented("B screening with distributed coefficients", __FILE__, __LINE__, class_desc());
                }
                else {
                  Eigen::Map<RowMatrix> C_X_view(
                      coefs_X_nu.at(Xsh.center).row(X.bfoff_in_atom).data(),
                      Xblk.atom_obsnbf, nbf
                  );
                  Kt_part.col(mu).segment(obs_atom_bfoff, obs_atom_nbf).noalias() +=
                      C_X_view * B_ish.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).transpose();

                  Kt_part.col(mu).noalias() += C_X_view.transpose()
                      * B_ish.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf).transpose();

                  Kt_part.col(mu).segment(obs_atom_bfoff, obs_atom_nbf).noalias() -=
                      C_X_view.middleCols(obs_atom_bfoff, obs_atom_nbf).transpose()
                      * B_ish.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf).transpose();
                }
              }
              else {
                const auto& C_X = coefs_transpose_[X];
                const auto& C_X_diff_X = C_X_diff.middleRows(
                    screen_B_ ? X.bfoff_in_block*Xblk.atom_obsnbf : 0,
                    screen_B_ ? Xblk.atom_obsnbf : 0
                );
                if(screen_B_) {
                  Kt_part.col(mu).segment(obs_atom_bfoff, obs_atom_nbf).noalias() +=
                      C_X_diff_X
                      * B_sd.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).tail(l_b_size).transpose();

                  Kt_part.col(mu).noalias() += C_X.transpose()
                      * B_sd.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).head(obs_atom_nbf).transpose();
                }
                else {
                  Kt_part.col(mu).segment(obs_atom_bfoff, obs_atom_nbf).noalias() +=
                      C_X * B_ish.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).transpose();

                  Kt_part.col(mu).noalias() += C_X.transpose()
                      * B_ish.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf).transpose();

                  Kt_part.col(mu).segment(obs_atom_bfoff, obs_atom_nbf).noalias() -=
                      C_X.middleCols(obs_atom_bfoff, obs_atom_nbf).transpose()
                      * B_ish.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf).transpose();
                }
              }

              //----------------------------------------//
            }
          }
          if(exact_diagonal_K_) {
            mt_timer.enter("exact diagonal subtract", ithr);
            assert(Xblk.atom_dfbfoff != NotAssigned and Xblk.atom_dfnbf != NotAssigned);

            if(Xblk.center != ish.center) {
              for(auto&& X : function_range(Xblk)) {
                for(auto&& mu : function_range(ish)) {
                  Kt_part.row(mu).middleCols(Xblk.atom_dfbfoff, Xblk.atom_dfnbf) -=
                      M_mu_X.row(mu.off*Xblk.nbf + X-Xblk.bfoff)
                      * coefs_transpose_[X].middleCols(ish.atom_bfoff, ish.atom_nbf).transpose();
                }
              }
            }

            mt_timer.exit(ithr);
          }
          mt_timer.exit(ithr);
          /*******************************************************/ #endif //2}}}            //latex `\label{sc:k3b:kcontrib:end}`
          /*-----------------------------------------------------*/
        } // end while get ish and Xblk

        delete[] b_buffer;
        delete[] B_ish_data;
        if(B_use_buffer_) {
          delete[] D_B_buff_data;
        }
        if(screen_B_) {
          delete[] D_sd_data;
          delete[] C_X_diff_data;
          if(screen_B_transfer_as_transpose_) {
            delete[] C_X_block_data;
          }
        }
        if(linK_block_rho_) {
          delete[] D_ordered_data;
        }
        if(distribute_coefficients_) {
          delete[] dt_ish_X_data;
          delete[] dt_prime_data;
          delete[] gt_ish_X_data;
        }

        //============================================================================//
        //----------------------------------------//
        // Sum Kt parts within node
        boost::lock_guard<boost::mutex> lg(tmp_mutex);
        Kt += Kt_part;
      }); // end create_thread
    } // end enumeration of threads

    compute_threads.join_all();
    mt_timer.exit();
    timer.insert(mt_timer);
    if(print_iteration_timings_) mt_timer.print(ExEnv::out0(), 12, 45);
    L_3.clear();
    L_B.clear();

  } // compute_threads is destroyed here
  /*****************************************************************************************/ #endif //1}}} //latex `\label{sc:k3b:end}`
  /*=======================================================================================*/
  /* Add back in the exact diagonal                       		                        {{{1 */ #if 1 //latex `\label{sc:kglobalsum}`
  if(exact_diagonal_K_) {
    timer.enter("exact diagonal contributions");
    boost::mutex tmp_mutex;
    boost::thread_group compute_threads;
    MultiThreadTimer mt_timer("threaded part", nthread_);
    //----------------------------------------//
    // reset the iteration over local pairs
    local_pairs_spot_ = 0;
    // Loop over number of threads
    for(int ithr = 0; ithr < nthread_; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){
        Eigen::MatrixXd Kt_part(nbf, nbf);
        Kt_part = Eigen::MatrixXd::Zero(nbf, nbf);
        //----------------------------------------//
        ShellData ish, jsh;
        while(get_shell_pair(ish, jsh, SignificantPairs)){
          // TODO Permutational symmetry (the 14 cases thing...)
          for(auto&& ksh : iter_shells_on_center(obs, ish.center)) {
            for(auto&& lsh : iter_shells_on_center(obs, jsh.center)) {
              if(not is_sig_pair(ksh, lsh)) continue;

              auto g4_ptr = ints_to_eigen(ish, jsh, ksh, lsh, tbis_[ithr], coulomb_oper_type_);
              const auto& g4 = *g4_ptr;

              for(auto&& mu : function_range(ish)) {
                for(auto&& rho : function_range(jsh)) {
                  if(rho > mu) continue;
                  // TODO more vectorization
                  for(auto&& nu : function_range(ksh)) {
                    Kt_part(mu, nu) +=
                      g4.row(mu.bfoff_in_shell*jsh.nbf + rho.bfoff_in_shell).segment(nu.bfoff_in_shell*lsh.nbf, lsh.nbf)
                        * D.col(rho).segment(lsh.bfoff, lsh.nbf);
                    if(mu != rho) {
                      Kt_part(rho, nu) +=
                        g4.row(mu.bfoff_in_shell*jsh.nbf + rho.bfoff_in_shell).segment(nu.bfoff_in_shell*lsh.nbf, lsh.nbf)
                          * D.col(mu).segment(lsh.bfoff, lsh.nbf);
                      if(ish.center != jsh.center) {
                        Kt_part.row(mu).segment(lsh.bfoff, lsh.nbf) +=
                            g4.row(mu.bfoff_in_shell*jsh.nbf + rho.bfoff_in_shell).segment(nu.bfoff_in_shell*lsh.nbf, lsh.nbf)
                            * D(rho, nu);
                        Kt_part.row(rho).segment(lsh.bfoff, lsh.nbf) +=
                            g4.row(mu.bfoff_in_shell*jsh.nbf + rho.bfoff_in_shell).segment(nu.bfoff_in_shell*lsh.nbf, lsh.nbf)
                            * D(mu, nu);
                      }
                    }
                  }
                }
              }

            }
          }
        }
        // Sum the thread's contributions to the node-level J
        boost::lock_guard<boost::mutex> tmp_lock(tmp_mutex);
        Kt += Kt_part;
      }); // end create_thread
    } // end enumeration of threads
    compute_threads.join_all();
    mt_timer.exit();
    timer.insert(mt_timer);
    timer.exit();
    if(xml_debug_) {
      #if DEBUG_K_INTERMEDIATES
      write_as_xml("Ktex", Ktex);
      #endif
    }
  } // compute_threads is destroyed here
  /*****************************************************************************************/ #endif //1}}} //latex `\label{sc:k3b:end}`
  /*=======================================================================================*/
  /* Global sum K                                         		                        {{{1 */ #if 1 //latex `\label{sc:kglobalsum}`
  //----------------------------------------//
  scf_grp_->sum(Kt.data(), nbf*nbf);
  //----------------------------------------//
  // Symmetrize K
  Eigen::MatrixXd K(nbf, nbf);
  K = Kt + Kt.transpose();
  //ExEnv::out0() << indent << "K checksum: " << scprintf("%20.15f", K.sum()) << std::endl;
  //----------------------------------------//
  /*****************************************************************************************/ #endif //1}}} //latex `\label{sc:kglobalsumend}`
  /*=======================================================================================*/
  /* Transfer K to a RefSCMatrix                           		                        {{{1 */ #if 1 // begin fold
  Ref<Integral> localints = integral()->clone();
  localints->set_basis(obs);
  Ref<PetiteList> pl = localints->petite_list();
  RefSCDimension obsdim = pl->AO_basisdim();
  RefSCMatrix result(
      obsdim,
      obsdim,
      obs->so_matrixkit()
  );
  result.assign(K.data());
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Clean up                                             		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  memory_used_ -= Z_size;
  //deallocate(D_data);
  delete[] D_data;
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  return result;
}

//////////////////////////////////////////////////////////////////////////////////


