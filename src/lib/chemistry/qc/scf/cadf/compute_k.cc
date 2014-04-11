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

#include <chemistry/qc/basis/petite.h>
#include <util/misc/xmlwriter.h>

#include "cadfclhf.h"

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


  #if DEBUG_K_INTERMEDIATES
  Eigen::MatrixXd Ktex(nbf, nbf);
  Ktex = Eigen::MatrixXd::Zero(nbf, nbf);
  #endif

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

  // Now make the linK lists if we're doing linK
  if(do_linK_){
    timer.enter("LinK lists");

    // First clear all of the lists
    //L_D.clear();
    L_DC.clear();
    L_3.clear();
    L_B.clear();

    //============================================================================//
    // Get the Frobenius norms of the density matrix shell blocks
    Eigen::MatrixXd D_frob(obs->nshell(), obs->nshell());
    do_threaded(nthread_, [&](int ithr){
      for(int lsh_index = ithr; lsh_index < obs->nshell(); lsh_index += nthread_) {           //latex `\label{sc:link:ld}`
        ShellData lsh(lsh_index, obs, dfbs_);
        for(auto&& jsh : shell_range(obs)) {
          double dnorm = D.block(lsh.bfoff, jsh.bfoff, lsh.nbf, jsh.nbf).norm();
          D_frob(lsh, jsh) = dnorm;
          //L_D[lsh].insert(jsh, dnorm);
        }
        //L_D[lsh].sort();
      }
    });
    Eigen::MatrixXd D_frob_sq(obs->nshell(), obs->nshell());
    D_frob_sq = D_frob.array().square();

    //----------------------------------------//                                             //latex `\label{sc:link:setupend}`
    // Form L_DC

    // TODO Distribute this over threads and MPI processes
    timer.enter("build L_DC");
    do_threaded(nthread_, [&](int ithr) {
      for(auto&& shpair : thread_over_range(
          product_range(shell_range(gbs_, dfbs_), shell_range(dfbs_)),
          ithr, nthread_)
      ) {                                                       //latex `\label{sc:link:ldc}`

        ShellData jsh, Xsh;
        boost::tie(jsh, Xsh) = shpair;
        const auto& Drho = D_frob.row(jsh);

        const auto& Cmaxes_X = Cmaxes_[Xsh];

        // TODO Optimize this to not be N^3

        double max_val = 0.0;

        if(use_norms_sigma_) {
          for(auto&& lsh : shell_range(obs)) {
            //const double this_val = Drho(lsh) * Cmaxes_X[lsh].value;
            //max_val += this_val * this_val;
            max_val += Drho(lsh) * Cmaxes_X[lsh].value;
          } // end loop over lsh
          //max_val = sqrt(max_val);
        }
        else {
          for(auto&& lsh : shell_range(obs)) {
            const double this_val = Drho(lsh) * Cmaxes_X[lsh].value;
            if(this_val > max_val) max_val = this_val;
          } // end loop over lsh
        }
        L_DC[jsh].insert(Xsh,                                                                //latex `\label{sc:link:ldcstore}`
            max_val * schwarz_df_[Xsh]
        );
      } // end loop over jsh, Xsh
    });

    do_threaded(nthread_, [&](int ithr){
      for(auto&& kvpair : thread_over_range(L_DC, ithr, nthread_)) {
        kvpair.second.sort();
      }
    });

    //----------------------------------------//                                             //latex `\label{sc:link:ldc:end}`
    // Form L_3
    timer.change("build L_3");                                                               //latex `\label{sc:link:l3}`

    // TODO Scale the B_screening_thresh also

    double epsilon = full_screening_thresh_;
    double epsilon_dist = distance_screening_thresh_;

    if(density_reset_){
      prev_density_frob_ = D_frob.norm();
      prev_epsilon_ = epsilon;
      prev_epsilon_dist_ = epsilon_dist;
    }
    else{
      if(scale_screening_thresh_) {
        const double ratio = D_frob.norm() / prev_density_frob_;
        epsilon = prev_epsilon_ * ratio;
        epsilon_dist = prev_epsilon_dist_ * ratio;
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


    if(scale_screening_thresh_) {
      if(linK_use_distance_ and epsilon != epsilon_dist) {
        ExEnv::out0() << indent << "  Effective LinK screening and LinK distance screening thresholds are now "
                      << scprintf("%5.2e and %5.2e", epsilon, epsilon_dist) << endl;
      }
      else {
        ExEnv::out0() << indent << "  Effective LinK screening threshold is now "
                      << scprintf("%5.2e", epsilon) << endl;
      }
    }


    // Loop over all jsh
    const auto& const_loc_pairs = local_pairs_linK_;
    const auto& const_loc_pairs_map = linK_local_map_;
    const auto& loc_pairs_end = const_loc_pairs.end();
    do_threaded(nthread_, [&](int ithr) {
      for(auto&& jsh : thread_over_range(shell_range(gbs_, dfbs_), ithr, nthread_)) {
        auto& L_sch_jsh = L_schwarz[jsh];
        for(auto&& Xsh : L_DC[jsh]) {

          const double pf = Xsh.value;
          const double eps_prime = epsilon / pf;
          const double eps_prime_dist = epsilon_dist / pf;
          const double Xsh_schwarz = schwarz_df_[Xsh];
          //bool jsh_added = false;
          //bool ish_found = false;
          auto found = const_loc_pairs_map.find(Xsh);
          if(found != const_loc_pairs_map.end()) {
            auto local_schwarz_jsh = L_sch_jsh.intersection_with(found->second);
            for(auto&& ish : local_schwarz_jsh) {

              double dist_factor = get_distance_factor(ish, jsh, Xsh);
              assert(ish.value == schwarz_frob_(ish, jsh));

              if(ish.value > eps_prime) {
                //jsh_added = true;
                if(!linK_use_distance_ or ish.value * dist_factor > eps_prime_dist) {
                  auto& L_3_ish_Xsh = L_3[{ish, Xsh}];
                  L_3_ish_Xsh.insert(jsh,
                      ish.value * dist_factor * Xsh_schwarz
                  );
                  if(screen_B_) {
                    L_3_ish_Xsh.add_to_aux_value_vector(ish.value * ish.value, D_frob_sq.col(jsh));
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
                break;
              }

            } // end loop over ish

          }
        }

      }
    });

    for(auto&& pair : L_3) {
      L_3_keys.emplace_back(
          pair.first.first,
          pair.first.second,
          pair.second.size()
      );
    }

    //----------------------------------------//                                             //latex `\label{sc:link:lB}`
    // Form L_B

    timer.change("sort L_3 and build L_B");


    do_threaded(nthread_, [&](int ithr){
      auto L_3_iter = L_3.begin();
      const auto& L_3_end = L_3.end();
      L_3_iter.advance(ithr);
      while(L_3_iter != L_3_end) {
        if(not linK_block_rho_) L_3_iter->second.set_sort_by_value(false);
        L_3_iter->second.sort();
        L_3_iter->second.set_aux_value(sqrt(L_3_iter->second.get_aux_value()));

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
            if(fabs(Cmaxes_[Xsh][lsh] * aux_vect(lsh)) > B_screening_thresh_) {
              L_B_ish_Xsh.insert(lsh);
              L_B_ish_Xsh.set_aux_value(lsh.nbf + L_B_ish_Xsh.get_aux_value());
            }
          }

          L_B_ish_Xsh.set_sort_by_value(false);
          L_B_ish_Xsh.sort();
        }

        L_3_iter.advance(nthread_);
      }
    });                                                                                      //latex `\label{sc:link:l3:end}`

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
    boost::thread_group compute_threads;
    // reset the iteration over local pairs
    local_pairs_spot_ = 0;
    // Loop over number of threads
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

    for(int ithr = 0; ithr < nthread_; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){                                              //latex `\label{sc:k3b:thrpatend}`

        Eigen::MatrixXd Kt_part(nbf, nbf);
        Kt_part = Eigen::MatrixXd::Zero(nbf, nbf);
        //----------------------------------------//
        ShellData ish;
        ShellBlockData<> Xblk;

        long b_buff_offset = 0;
        size_t actual_size = 0;
        actual_size = B_buffer_size_/sizeof(double);
        double* __restrict__ b_buffer = new double[actual_size];
        double* __restrict__ B_ish_data = new double[max_fxn_obs_*max_fxn_dfbs_*nbf];
        double* __restrict__ D_B_buff_data;
        double* __restrict__ D_sd_data;
        int D_B_buff_n = 0;
        if(B_use_buffer_) {
          D_B_buff_data = new double[nbf*nbf];
          D_B_buff_n = nbf;
        }
        if(screen_B_) {
          D_sd_data = new double[nbf*nbf];
        }
        Eigen::Map<ColMatrix> D_B_buff(D_B_buff_data, D_B_buff_n, D_B_buff_n);

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

          assert(!do_linK_ or Xblk.nshell == 1);
          auto&& Xsh = Xblk.first_shell;

          // Form the B_same, B_diff, D_same, D_diff, and C_X matrices
          std::vector<Eigen::MatrixXd> C_X_diff;
          //ColMatrix D_diff, D_same;
          //RowMatrix B_diff, B_same;
          //ColMatrix D_sd;
          //RowMatrix B_sd;
          int l_b_size;
          int B_sd_rows = 0, B_sd_cols = 0;
          int D_sd_rows = 0, D_sd_cols = 0;

          if(screen_B_) {

            mt_timer.enter("build screened B parts", ithr);

            const auto& L_B_ish_Xsh = L_B[{ish, Xsh}];
            l_b_size = int(L_B_ish_Xsh.get_aux_value());
            if(l_b_size == 0){
              mt_timer.exit(ithr);
              mt_timer.exit(ithr); // compute B
              continue;
            }

            C_X_diff.resize(Xsh.nbf);
            //D_diff.resize(l_b_size, nbf);
            //D_same.resize(Xblk.atom_obsnbf, nbf);
            //D_same = D.middleCols(Xblk.atom_obsbfoff, Xblk.atom_obsnbf).transpose();

            //B_diff.resize(ish.nbf*Xblk.nbf, l_b_size);
            //B_diff = RowMatrix::Zero(ish.nbf*Xblk.nbf, l_b_size);
            //B_same.resize(ish.nbf*Xblk.nbf, Xblk.atom_obsnbf);
            //B_same = RowMatrix::Zero(ish.nbf*Xblk.nbf, Xblk.atom_obsnbf);

            //B_sd.resize(ish.nbf*Xblk.nbf, Xblk.atom_obsnbf + l_b_size);
            //B_sd = RowMatrix::Zero(ish.nbf*Xblk.nbf, Xblk.atom_obsnbf + l_b_size);
            B_sd_rows = ish.nbf*Xblk.nbf;
            B_sd_cols = Xblk.atom_obsnbf + l_b_size;
            D_sd_rows = Xblk.atom_obsnbf + l_b_size;
            D_sd_cols = nbf;
          }

          Eigen::Map<ColMatrix> D_sd(D_sd_data, D_sd_rows, D_sd_cols);
          //D_sd = Eigen::::Zero

          if(screen_B_) {
            D_sd.topRows(Xblk.atom_obsnbf) = D.middleCols(Xblk.atom_obsbfoff, Xblk.atom_obsnbf).transpose();


            int Xblk_offset = 0;
            for(auto&& X : function_range(Xblk)) {
              assert(Xblk.atom_obsnbf > 0);
              C_X_diff[Xblk_offset].resize(Xblk.atom_obsnbf, l_b_size);
              Xblk_offset++;
            }

            int block_offset = 0;
            const auto& L_B_ish_Xsh = L_B[{ish, Xsh}];
            for(auto&& lblk : shell_block_range(L_B_ish_Xsh, Contiguous)) {
              Xblk_offset = 0;
              for(auto&& X : function_range(Xblk)) {
                C_X_diff[Xblk_offset].middleCols(block_offset, lblk.nbf) = coefs_transpose_[X].middleCols(lblk.bfoff, lblk.nbf);
                Xblk_offset++;
              }

              // Would this be more efficient if we did it the other way and then did a transpose in place?
              //D_diff.middleRows(block_offset, lblk.nbf) = D.middleCols(lblk.bfoff, lblk.nbf).transpose();
              D_sd.middleRows(Xblk.atom_obsnbf + block_offset, lblk.nbf) = D.middleCols(lblk.bfoff, lblk.nbf).transpose();

              block_offset += lblk.nbf;
            }

            mt_timer.exit(ithr);

          }

          Eigen::Map<RowMatrix> B_sd(B_ish_data, B_sd_rows, B_sd_cols);


          // Create B_ish and the B buffer
          Eigen::Map<RowMatrix> B_ish(B_ish_data, ish.nbf * Xblk.nbf, nbf);
          if(not screen_B_) {
            ::memset(B_ish_data, 0, ish.nbf*Xblk.nbf*nbf*sizeof(double));
          }
          else{
            ::memset(B_ish_data, 0, B_sd_rows*B_sd_cols*sizeof(double));
          }
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


          // What list of J are we using?
          OrderedShellList* jlist;
          if(not do_linK_) {
            jlist = new OrderedShellList(sig_partners_[ish], gbs_, dfbs_);
          }
          else {
            jlist = &(L_3[{ish, Xsh}]);
          }

          // Reordering of D, only if linK_block_rho_ is true
          int block_offset = 0;
          Eigen::MatrixXd D_ordered;                                                   //latex `\label{sc:k3b:reD}`
          if(do_linK_ and linK_block_rho_) {

            mt_timer.enter("rearrange D", ithr);
            D_ordered.resize(nbf, nbf);
            for(auto jblk : shell_block_range(L_3[{ish, Xsh}], Contiguous)){

              D_ordered.middleCols(block_offset, jblk.nbf) = D.middleCols(jblk.bfoff, jblk.nbf);
              block_offset += jblk.nbf;

            }                                                                            //latex `\label{sc:k3b:reD:end}`
            mt_timer.exit(ithr);

          }

          //============================================================================//
          // Loop over the largest blocks of J at once that we can

          int restrictions = linK_block_rho_ ? NoRestrictions : Contiguous;
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

              if(exact_diagonal_K_) {
                subtimer.change(ex_timer);

                if(jsblk.center == Xblk.center or ish.center == Xblk.center) {
                  // Build W and Wbar
                  for(auto&& mu : function_range(ish)) {
                    for(auto&& Y : iter_functions_on_center(dfbs_, ish.center)) {
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

            // Eigen version
            if(B_use_buffer_) {

              if(b_buff_offset + jblk.nbf > b_buff_nrows) {
                if(b_buff_offset == 0) {
                  throw SCException("B_buffer_size smaller than single contiguous block.  Set B_use_buffer to no and try again.");
                }
                B_ish.noalias() += 2.0 * B_buffer_mat.topRows(b_buff_offset).transpose() * D_B_buff.leftCols(b_buff_offset).transpose();
                b_buff_offset = 0;
              }


              std::copy(D_data + jblk.bfoff, D_data + jblk.bfoff + jblk.nbf, D_B_buff_data + b_buff_offset);

              //D_B_buff.middleCols(b_buff_offset, jblk.nbf) = D.middleCols(jblk.bfoff, jblk.nbf);
              b_buff_offset += jblk.nbf;

            }
            else {
              if(linK_block_rho_) {
                B_ish.noalias() += 2.0 * g3.transpose() * D_ordered.middleCols(block_offset, jblk.nbf).transpose();
              }
              else {
                if(screen_B_) {
                  B_sd.noalias() += 2.0 * g3.transpose() * D_sd.middleCols(jblk.bfoff, jblk.nbf).transpose();
                  //B_diff.noalias() += 2.0 * g3.transpose() * D_diff.middleCols(jblk.bfoff, jblk.nbf).transpose();
                  //B_same.noalias() += 2.0 * g3.transpose() * D_same.middleCols(jblk.bfoff, jblk.nbf).transpose();
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
            const auto& C_X = coefs_transpose_[X];
            //const auto& C_X_diff_X = C_X_diff[X.bfoff_in_block];
            for(auto&& mu : function_range(ish)) {

              // B_mus[mu.bfoff_in_shell] is (nbf x Ysh.nbf)
              // C_Y is (Y.{obs_}atom_nbf x nbf)
              // result should be (Y.{obs_}atom_nbf x 1)

              if(screen_B_) {
                // TODO Offset C_X_diff outside of the mu loop
                Kt_part.col(mu).segment(obs_atom_bfoff, obs_atom_nbf).noalias() +=
                    C_X_diff[X.bfoff_in_block]
                    * B_sd.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).tail(l_b_size).transpose();
                //Kt_part.col(mu).segment(obs_atom_bfoff, obs_atom_nbf).noalias() +=
                //    C_X_diff[X.bfoff_in_block]
                //    * B_diff.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).transpose();

                Kt_part.col(mu).noalias() += C_X.transpose()
                    * B_sd.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).head(obs_atom_nbf).transpose();
                //Kt_part.col(mu).noalias() += C_X.transpose()
                //    * B_same.row(mu.bfoff_in_shell*Xblk.nbf + X.bfoff_in_block).transpose();
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
        } // end while get ish Xblk pair
        //deallocate(b_buffer);
        delete[] b_buffer;
        delete[] B_ish_data;
        if(B_use_buffer_) {
          delete[] D_B_buff_data;
        }
        if(screen_B_) {
          delete[] D_sd_data;
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


