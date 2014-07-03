//
// new_exchange.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Jun 30, 2014
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

#include <tuple>

#include <chemistry/qc/basis/petite.h>

#include "cadfclhf.h"
#include "assignments.h"

using namespace sc;
using std::cout;
using std::endl;

typedef uint64_t uli;

// Help out the Eclipse parser without affecting performance
#ifndef SH
#  if ECLIPSE_PARSER_ONLY
#    define SH ShellData
#    define SHV ShellDataWithValue
#    define BF BasisFunctionData
#    define BLK ShellBlockData<>
#  else
#    define SH auto&&
#    define SHV auto&&
#    define BF auto&&
#    define BLK auto&&
#  endif
#endif


RefSCMatrix
CADFCLHF::new_compute_K()
{
  /*=======================================================================================*/
  /* Setup                                                 		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // Convenience variables
  Timer timer("compute K");
  const int me = scf_grp_->me();
  const int n_node = scf_grp_->n();
  const int nbf = gbs_->nbasis();
  const int nsh = gbs_->nshell();
  const int dfnbf = dfbs_->nbasis();
  const int dfnsh = dfbs_->nshell();
  const int natom = gbs_->ncenter();

  auto& my_part = *(assignments_new_k_->nodes[me]);
  const auto& g2 = *g2_full_ptr_;

  //----------------------------------------//
  // Get the density in an Eigen::Map form
  double* __restrict__ D_data = new double[nbf*nbf];
  D_.convert(D_data);
  Eigen::Map<ColMatrix> D(D_data, nbf, nbf);
  auto D_tracker = hold_memory(nbf*nbf*sizeof(double) + sizeof(decltype(D)));
  // Match density scaling in old code:
  D *= 0.5;

  //----------------------------------------//
  // Initialize the K tilde matrix
  ColMatrix Kt(nbf, nbf);
  Kt.setZero();
  auto Kt_tracker = hold_memory(nbf*nbf*sizeof(double) + sizeof(ColMatrix));


  //----------------------------------------//
  // Lots of other stuff not implemented also,
  //   but a couple deserve specific attention
  if(exact_diagonal_K_) {
    throw FeatureNotImplemented("Exact diagonal with new exchange algorithm", __FILE__, __LINE__, class_desc());
  }
  if(not do_linK_) {
    throw FeatureNotImplemented("new exchange algorithm always does linK, so do_linK_ = false is not implemented", __FILE__, __LINE__, class_desc());
  }
  if(not screen_B_) {
    throw FeatureNotImplemented("new exchange algorithm always screens the B intermediate, so screen_B_ = false is not implemented", __FILE__, __LINE__, class_desc());
  }

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/


  /*=======================================================================================*/
  /* Make the CADF-LinK lists                                                         {{{1 */ #if 1 // begin fold

  timer.enter("LinK lists");

  typedef std::unordered_map<std::pair<int, int>, ContiguousShellBlockList, sc::hash<std::pair<int, int>>> FastBlockListMap2;
  FastBlockListMap2 L_3_blocked;
  FastBlockListMap2 L_B_blocked;
  IndexListMap L_C;
  IndexListMap L_DF;
  std::vector<int> mu_to_do;
  std::atomic<uli> max_L_B_size(0);
  std::atomic<uli> max_L_3_size(0);

  {
    /*-----------------------------------------------------------------*/
    /* Get the Frobenius norms of the density matrix shell blocks {{{2 */ #if 2 // begin fold
    timer.enter("build frob D");
    MAKE_MATRIX(RowMatrix, D_frob_sq, nsh, nsh);
    do_threaded(nthread_, [&](int ithr){
      for(SH lsh : thread_over_range(shell_range(gbs_, dfbs_), ithr, nthread_)) {
        for(SH jsh : shell_range(gbs_)) {
          double dnorm = D.block(lsh.bfoff, jsh.bfoff, lsh.nbf, jsh.nbf).squaredNorm();
          D_frob_sq(lsh, jsh) = dnorm;
        }
      }
    });
    MAKE_MATRIX(RowMatrix, D_frob, nsh, nsh);
    D_frob = D_frob_sq.cwiseSqrt();
    timer.exit();
    /*******************************************************************/ #endif //2}}}
    /*-----------------------------------------------------------------*/

    /*-----------------------------------------------------*/
    /* Adjust the thresholds and tell user about it   {{{2 */ #if 2 // begin fold

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
      prev_epsilon_B_ = epsilon_B;

      if(full_screening_expon_ != 1.0) {
        epsilon = pow(epsilon, full_screening_expon_);                                      //latex `\label{sc:link:expon}`
        epsilon_dist = pow(epsilon_dist, full_screening_expon_);
      }

      epsilon = std::max(full_screening_thresh_min_, epsilon);
      epsilon_dist = std::max(full_screening_thresh_min_, epsilon_dist);
      // TODO Make this an option
      epsilon_B = std::max(full_screening_thresh_min_*(B_screening_thresh_/full_screening_thresh_), epsilon_B);

      epsilon = std::min(epsilon, full_screening_thresh_);
      epsilon_dist = std::min(epsilon_dist, distance_screening_thresh_);
      epsilon_B = std::min(epsilon_B, B_screening_thresh_);
    }


    // Update the user on the effective thresholds
    if(scale_screening_thresh_) {
      if(linK_use_distance_ and epsilon != epsilon_dist) {
        ExEnv::out0() << indent << "  Effective LinK screening, LinK distance screening, and B screening thresholds are now "
                      << scprintf("%5.2e, %5.2e, and %5.2e", epsilon, epsilon_dist, epsilon_B) << endl;
      }
      else {
        if(epsilon_B != epsilon) {
          ExEnv::out0() << indent << "  Effective LinK screening and B screening thresholds are now "
                        << scprintf("%5.2e and %5.2e", epsilon, epsilon_B) << endl;
        }
        else {
          ExEnv::out0() << indent << "  Effective LinK screening (and B screening) threshold is now "
                        << scprintf("%5.2e", epsilon) << endl;
        }
      }
    }

    /********************************************************/ #endif //2}}}
    /*-----------------------------------------------------*/

    /*-----------------------------------------------------*/
    /* Build d_bar and form L_DC                      {{{2 */ #if 2 // begin fold

    timer.change("build L_DC");

    {
      Ref<MessageGrp> X_grp = scf_grp_->split(my_part.cluster->index);

      // Build d_bar
      MAKE_MATRIX(RowMatrix, d_bar, nsh, my_part.dfnsh());
      d_bar.setZero();
      int my_begin, my_size;
      get_split_range_part(X_grp, 0, nsh, my_begin, my_size);
      do_threaded(nthread_, [&](int ithr) {
         int my_begin_obs, my_size_obs;
         get_split_range_part(ithr, nthread_, 0, gbs_->nshell(), my_begin_obs, my_size_obs);
         d_bar.middleRows(my_begin_obs, my_size_obs).noalias() =
             D_frob.middleRows(my_begin_obs, my_size_obs).middleCols(my_begin, my_size)
             * C_bar_mine_.middleRows(my_begin, my_size) * schwarz_df_mine_.asDiagonal();
      });
      X_grp->sum(d_bar.data(), d_bar.rows() * d_bar.cols());

      // Form the L_DC lists
      do_threaded(nthread_, [&](int ithr) {
        for(auto&& jsh : thread_over_range(shell_range(gbs_), ithr, nthread_)) {
          auto& L_DC_jsh = L_DC[jsh];
          L_DC_jsh.acquire_and_sort(
              my_part.assigned_dfbs_shells(),
              d_bar.data() + my_part.dfnsh() * jsh,
              0.0, true // sort by value
          );
          L_DC_jsh.set_basis(dfbs_, gbs_);
        }
      });
    } // d_bar and X_grp deleted

    // For whatever reason, splitting the scf_grp_ messes up the SCFormIO processor 0 label
    sc::SCFormIO::init_mp(scf_grp_->me());


    /*******************************************************/ #endif //2}}}
    /*-----------------------------------------------------*/

    /*-----------------------------------------------------*/
    /* Build L_3                                      {{{2 */ #if 2 // begin fold

    const auto& const_loc_pairs = local_pairs_linK_;
    const auto& const_loc_pairs_map = linK_local_map_;
    const auto& loc_pairs_end = const_loc_pairs.end();
    do_threaded(nthread_, [&](int ithr) {
      for(SH jsh : thread_over_range(shell_range(gbs_, dfbs_), ithr, nthread_)) {
        auto& L_sch_jsh = L_schwarz[jsh];

        for(SHV Xsh : L_DC[jsh]) {

          auto found = const_loc_pairs_map.find(Xsh);

          if(found != const_loc_pairs_map.end()) {

            const double pf = Xsh.value;
            const double eps_prime = epsilon / pf;
            const double eps_prime_dist = epsilon_dist / pf;

            auto local_schwarz_jsh = L_sch_jsh.intersection_with(found->second);

            bool found_ish = false;

            for(SHV ish : local_schwarz_jsh) {

              const double dist_factor = get_distance_factor(ish, jsh, Xsh);

              if(ish.value > eps_prime) {
                found_ish = true;
                if(!linK_use_distance_ or ish.value * dist_factor > eps_prime_dist) {
                  auto& L_3_ish_Xsh = L_3[{ish, Xsh}];
                  L_3_ish_Xsh.insert(jsh,
                      ish.value * dist_factor,
                      jsh.nbf
                  );

                  L_C[ish].insert(jsh);

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

            // We know we can exit early if this Xsh iteration contains all local ish
            if(!found_ish and local_schwarz_jsh.size() == my_part.obs_shells_to_do.size()) {
              break;
            }

          }

        }


      }

    });

    L_DC.clear();

    /*******************************************************/ #endif //2}}}
    /*-----------------------------------------------------*/

    /*-----------------------------------------------------*/
    /* Sort L_3 and build L_B                         {{{2 */ #if 2 // begin fold

    timer.change("sort L_3 and build L_B");

    {
      std::mutex L_3_mutex;
      std::mutex L_B_mutex;

      do_threaded(nthread_, [&](int ithr){

        auto L_3_iter = L_3.begin();
        const auto& L_3_end = L_3.end();
        L_3_iter.advance(ithr);
        for(; L_3_iter != L_3_end; L_3_iter.advance(nthread_)) {
          L_3_iter->second.set_sort_by_value(false);
          L_3_iter->second.sort();

          // Build the B screening list

          // Convenience variables
          const ShellData ish(L_3_iter->first.first, gbs_, dfbs_);
          const ShellData Xsh(L_3_iter->first.second, dfbs_, gbs_);
          auto& L_B_ish_Xsh = L_B[{ish, Xsh}];
          auto& L_3_ish_Xsh = L_3_iter->second;

          // Find the number of basis functions in L_3_ish_Xsh
          int l3_nbf = 0;
          for(auto&& jblk : shell_block_range(L_3_ish_Xsh, NoRestrictions)) {
            l3_nbf += jblk.nbf;
          }
          L_3[{ish, Xsh}].nbf = l3_nbf;
          if(l3_nbf == 0) { continue; }

          // Build up the distance factors for use in b_bar build
          int jsh_offset = 0;
          Eigen::VectorXd dist_facts;
          if(screen_B_use_distance_) {
            dist_facts.resize(L_3_ish_Xsh.size());
            for(auto&& jsh : L_3_ish_Xsh) {
              const double dist_fact = get_distance_factor(ish, jsh, Xsh);
              if(use_norms_B_) {
                dist_facts(jsh_offset)  = dist_fact * dist_fact;
              }
              else {
                dist_facts(jsh_offset)  = dist_fact;
              }
              ++jsh_offset;
            }
          }

          // Build b_bar
          Eigen::VectorXd b_bar(nsh);
          b_bar.setZero();
          jsh_offset = 0;
          for(auto&& jblk : shell_block_range(L_3_ish_Xsh, Contiguous)) {
            if(use_norms_B_) {
              if(screen_B_use_distance_) {
                b_bar += dist_facts.segment(jsh_offset, jblk.nshell).asDiagonal()
                    * D_frob.middleCols(jblk.first_shell, jblk.nshell).rowwise().squaredNorm()
                    * schwarz_frob_.col(ish).segment(jblk.first_shell, jblk.nshell).squaredNorm();
                jsh_offset += jblk.nshell;
              }
              else{
                b_bar += D_frob.middleCols(jblk.first_shell, jblk.nshell).rowwise().squaredNorm()
                    * schwarz_frob_.col(ish).segment(jblk.first_shell, jblk.nshell).squaredNorm();
              }
            }
            else {
              if(screen_B_use_distance_) {
                b_bar += D_frob.middleCols(jblk.first_shell, jblk.nshell)
                    * dist_facts.segment(jsh_offset, jblk.nshell).asDiagonal()
                    * schwarz_frob_.col(ish).segment(jblk.first_shell, jblk.nshell);
                jsh_offset += jblk.nshell;
              }
              else {
                b_bar += D_frob.middleCols(jblk.first_shell, jblk.nshell)
                    * schwarz_frob_.col(ish).segment(jblk.first_shell, jblk.nshell);
              }
            }
          }
          if(use_norms_B_) {
            b_bar = b_bar.cwiseSqrt();
          }
          b_bar *= schwarz_df_[Xsh];
          b_bar.array() *= C_bar_mine_.col(my_part.cluster->dfbs_shell_map[Xsh]).array();
          b_bar = b_bar.cwiseAbs();

          // Acquire L_B from b_bar
          // Avoid double-counting by zeroing the same-atom block
          b_bar.segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh).setZero();
          L_B_ish_Xsh.acquire_and_sort(
              b_bar.data(), nsh, epsilon_B, false
          );
          L_B_ish_Xsh.set_basis(gbs_, dfbs_);

          // Find the number of basis functions in L_B_ish_Xsh
          int lb_nbf = 0;
          for(auto&& lblk : shell_block_range(L_B_ish_Xsh, NoRestrictions)) {
            lb_nbf += lblk.nbf;
          }
          L_B_ish_Xsh.nbf = lb_nbf;


          // We can't exclude empty L_B lists here because the same
          //  atom terms will still need to be done

          // Add Xsh to the list of df shells to do for a given ish
          L_DF[ish].insert(Xsh);

          {
            // Pre-block L_3
            ContiguousShellBlockList L3_mu_X_list(
                L_3_ish_Xsh, gbs_, dfbs_, { NoRestrictions }
            );
            std::lock_guard<std::mutex> lg(L_3_mutex);
            L_3_blocked[{ish, Xsh}] = L3_mu_X_list;
            if(l3_nbf > max_L_3_size) max_L_3_size = l3_nbf;
          }

          {
            // Pre-block L_B
            ContiguousShellBlockList LB_mu_X_list(
                L_B_ish_Xsh, gbs_, dfbs_, { NoRestrictions }
            );
            std::lock_guard<std::mutex> lg(L_B_mutex);
            L_B_blocked[{ish, Xsh}] = LB_mu_X_list;
            if(lb_nbf > max_L_B_size) max_L_B_size = lb_nbf;
          }


        }

      });
    }

    timer.change("misc");

    // Sort L_C(mu) for each mu
    do_threaded(nthread_, [&](int ithr){

      auto L_C_iter = L_C.begin();
      const auto& L_C_end = L_C.end();
      L_C_iter.advance(ithr);
      for(; L_C_iter != L_C_end; L_C_iter.advance(nthread_)) {
        L_C_iter->second.set_sort_by_value(false);
        L_C_iter->second.sort();
        L_C_iter->second.set_basis(gbs_, dfbs_);
      }

    });

    // Sort L_DF(mu) for each mu
    do_threaded(nthread_, [&](int ithr){

      auto L_DF_iter = L_DF.begin();
      const auto& L_DF_end = L_DF.end();
      L_DF_iter.advance(ithr);
      for(; L_DF_iter != L_DF_end; L_DF_iter.advance(nthread_)) {
        L_DF_iter->second.set_sort_by_value(false);
        L_DF_iter->second.sort();
        L_DF_iter->second.set_basis(dfbs_, gbs_);
      }

    });

    // Make a list of mu where L_DF and L_C are not empty
    for(auto&& pair : L_C) {
      if(L_DF.find(pair.first) != L_DF.end()) {
        mu_to_do.push_back(pair.first);
      }
    }



    /*******************************************************/ #endif //2}}}
    /*-----------------------------------------------------*/

  }

  timer.exit();

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/


  /*=======================================================================================*/
  /* Main loop over mu and X to compute K tilde contributions                         {{{1 */ #if 1 // begin fold
  /*---------------------------------------------------------------------------------------*/

  timer.enter("main loop");

  std::mutex mu_list_mutex;
  auto mu_iter = mu_to_do.begin();

  auto get_next_ish = [&](ShellData& ishell) -> bool {
    if(mu_iter == mu_to_do.end()) return false;
    else {
      std::lock_guard<std::mutex> lg(mu_list_mutex);
      if(mu_iter == mu_to_do.end()) return false;
      ishell = ShellData(*mu_iter, gbs_, dfbs_);
      ++mu_iter;
      return true;
    }
  };

  MultiThreadTimer mt_timer("threaded part", nthread_);

  do_threaded(nthread_, [&](int ithr) {

    /*-----------------------------------------------------*/
    /* Thread-local setup                             {{{2 */ #if 2 // begin fold

    // Allocate memory buffers outside of loop to avoid malloc overhead
    double* __restrict__ ints_buffer = new double[B_buffer_size_/sizeof(double)];
    double* __restrict__ rec_coefs_data = new double[max_fxn_obs_todo_ * nbf * max_fxn_atom_dfbs_*2];
    double* __restrict__ B_sd_data = new double[max_fxn_obs_todo_ * max_fxn_dfbs_todo_ * (max_L_B_size + max_fxn_atom_obs_)];
    double* __restrict__ D_sd_data = new double[max_L_3_size * (max_L_B_size + max_fxn_atom_obs_)];
    double* __restrict__ C_X_diff_data = new double[max_obs_atom_fxn_on_dfbs_center_todo_*max_fxn_dfbs_todo_*max_L_B_size];

    Eigen::Map<RowMatrix> B_sd(B_sd_data, 0, 0);
    Eigen::Map<RowMatrix> B_sd_other(B_sd_data, 0, 0);
    Eigen::Map<ColMatrix> D_sd(D_sd_data, 0, 0);
    Eigen::Map<RowMatrix> C_X_diff(D_sd_data, 0, 0);

    ShellData ish;

    /*******************************************************/ #endif //2}}}
    /*-----------------------------------------------------*/

    while(get_next_ish(ish)) {

      mt_timer.enter("recompute coefs", ithr);

      /*------------------------------------------------------*/
      /* Compute the part of C that we need to recompute {{{2 */ #if 2 // begin fold

      std::unordered_map<int, Eigen::Map<RowMatrix>> rec_coefs;
      uli coef_spot_offset = 0;

      // Compute the coefficients that we need for the current ish
      for(BLK jblk : shell_block_range(L_C[ish], SameCenter, NoMaximumBlockSize)) {
        const int df_size = (ish.center == jblk.center) ? ish.atom_dfnbf : ish.atom_dfnbf + jblk.atom_dfnbf;
        rec_coefs.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(jblk.center),
            std::forward_as_tuple(
                rec_coefs_data + coef_spot_offset,
                jblk.atom_nbf * ish.nbf, df_size
            )
        );
        coef_spot_offset += jblk.atom_nbf * ish.nbf * df_size;
        auto& rec_C = rec_coefs.at(jblk.center);
        rec_C = RowMatrix::Constant(jblk.atom_nbf * ish.nbf, df_size, std::numeric_limits<double>::infinity());

        // Get the decomposed two center ints
        auto decomp = get_decomposition(ish, jblk.first_shell, metric_ints_2c_[ithr]);

        if(ish.center == jblk.center) {

          auto Zblk = ShellBlockData<>::atom_block(ish.center, dfbs_, gbs_);

          // Loop over jsh in the block
          for(SH jsh : shell_range(jblk)) {

            // Compute the integrals we need
            auto g3_part = ints_to_eigen_map(
                jsh, ish, Zblk,
                eris_3c_[ithr], coulomb_oper_type_,
                ints_buffer
            );

            // Get the coefficients for each basis function pair
            for(BF mu : function_range(ish)) {
              for(BF rho : function_range(jsh)) {
                rec_C.row(rho.bfoff_in_atom*ish.nbf + mu.off).transpose() = decomp->solve(
                    g3_part.row(rho.off*ish.nbf + mu.off).transpose()
                );

              } // end nu loop
            } // end mu loop

          } // end loop over jblk

        } // end same center case
        else {

          // Same thing for the different center case, but joined blocks this time
          auto Zblk = ShellBlockData<>::atom_block(ish.center, dfbs_, gbs_) + ShellBlockData<>::atom_block(jblk.center, dfbs_, gbs_) ;

          // Loop over jsh in the block
          for(SH jsh : shell_range(jblk)) {

            // Compute the integrals we need
            auto g3_part = ints_to_eigen_map(
                jsh, ish, Zblk,
                eris_3c_[ithr], coulomb_oper_type_,
                ints_buffer
            );

            // Get the coefficients for each basis function pair
            for(BF mu : function_range(ish)) {
              for(BF rho : function_range(jsh)) {
                rec_C.row(rho.bfoff_in_atom*ish.nbf + mu.off).transpose() = decomp->solve(
                    g3_part.row(rho.off*ish.nbf + mu.off).transpose()
                );

              } // end nu loop
            } // end mu loop

          } // end loop over jblk

        } // end different center case

      } // end loop over jblk in L_C

      /********************************************************/ #endif //2}}}
      /*------------------------------------------------------*/


      // Now loop over the local Xsh connected to the current ish
      for(SH Xsh : L_DF[ish]) {

        mt_timer.change("Compute B", ithr);

        /*-----------------------------------------------------*/
        /* Setup                                          {{{2 */ #if 2 // begin fold

        // Convenience variables
        const auto& L_B_ish_Xsh = L_B[{ish, Xsh}];
        const auto& L3iter = L_3.find({ish, Xsh});
        const auto& L_3_ish_Xsh = L3iter->second;
        bool L_B_empty = (L_B_blocked.find({ish, Xsh}) == L_B_blocked.end());
        const auto& L_B_blocks = L_B_blocked[{ish, Xsh}];
        const auto& L_3_blocks = L_3_blocked[{ish, Xsh}];
        const int l_b_size = L_B_ish_Xsh.nbf;
        const int l_3_size = L_3_ish_Xsh.nbf;
        int jblk_offset, block_offset;

        /*******************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/

        /*-----------------------------------------------------*/
        /* Initialize memory for B screening              {{{2 */ #if 2 // begin fold

        // Form the B_sd, D_sd, and C_X matrices

        new (&B_sd) Eigen::Map<RowMatrix>(B_sd_data, ish.nbf*Xsh.nbf, Xsh.atom_obsnbf + l_b_size);
        new (&B_sd_other) Eigen::Map<RowMatrix>(B_sd_data, ish.nbf, Xsh.nbf*(Xsh.atom_obsnbf + l_b_size));
        new (&D_sd) Eigen::Map<ColMatrix>(D_sd_data, Xsh.atom_obsnbf + l_b_size, l_3_size);
        new (&C_X_diff) Eigen::Map<RowMatrix>(C_X_diff_data, Xsh.nbf*Xsh.atom_obsnbf, l_b_size);
        B_sd.setZero();

        jblk_offset = 0;
        for(ShellBlockData<> jblk : L_3_blocks[NoRestrictions]) {
          D_sd.block(
              0, jblk_offset,
              Xsh.atom_obsnbf, jblk.nbf
          ) = D.block(
              Xsh.atom_obsbfoff, jblk.bfoff,
              Xsh.atom_obsnbf, jblk.nbf
          );
          jblk_offset += jblk.nbf;
        }

        block_offset = 0;
        const auto& C_X_block = coefs_X_nu_other.at(Xsh.center).middleRows(
              Xsh.bfoff_in_atom*Xsh.atom_obsnbf, Xsh.nbf*Xsh.atom_obsnbf
        );
        if(not L_B_empty) {
          for(ShellBlockData<> lblk : L_B_blocks[NoRestrictions]) {
            C_X_diff.middleCols(block_offset, lblk.nbf) = C_X_block.middleCols(lblk.bfoff, lblk.nbf);
            block_offset += lblk.nbf;
          }
        }
        block_offset = 0;

        // TODO combine this loop with the one above it when we no longer need seperate times
        if(not L_B_empty) {
          for(ShellBlockData<> lblk : L_B_blocks[NoRestrictions]) {
            jblk_offset = 0;
            for(ShellBlockData<> jblk : L_3_blocks[NoRestrictions]) {
              D_sd.block(Xsh.atom_obsnbf + block_offset, jblk_offset,
                  lblk.nbf, jblk.nbf
              ) = D.block(
                  lblk.bfoff, jblk.bfoff,
                  lblk.nbf, jblk.nbf
              );
              jblk_offset += jblk.nbf;
            }
            block_offset += lblk.nbf;
          }
        }


        /*******************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/

        /*-----------------------------------------------------*/
        /* Build B                                        {{{2 */ #if 2 // begin fold

        block_offset = 0;

        for(auto&& jblk : shell_block_range(L_3_ish_Xsh, NoRestrictions)){

          //----------------------------------------//

          // Compute the integrals

          auto g3_in = ints_to_eigen_map(
              jblk, ish, Xsh,
              eris_3c_[ithr], coulomb_oper_type_,
              ints_buffer
          );

          //----------------------------------------//

          // Two body part

          jblk_offset = 0;
          for(const auto&& jsblk : shell_block_range(jblk, Contiguous|SameCenter)) {

            g3_in.middleRows(jblk_offset*ish.nbf, jsblk.nbf*ish.nbf) -= 0.5 *
                rec_coefs.at(jsblk.center).middleRows(jsblk.bfoff_in_atom*ish.nbf, jsblk.nbf*ish.nbf).leftCols(ish.atom_dfnbf)
                * g2.middleRows(ish.atom_dfbfoff, ish.atom_dfnbf).middleCols(Xsh.bfoff, Xsh.nbf);
            if(ish.center != jsblk.center) {
              g3_in.middleRows(jblk_offset*ish.nbf, jsblk.nbf*ish.nbf) -= 0.5 *
                  rec_coefs.at(jsblk.center).middleRows(jsblk.bfoff_in_atom*ish.nbf, jsblk.nbf*ish.nbf).rightCols(jsblk.atom_dfnbf)
                  * g2.middleRows(jsblk.atom_dfbfoff, jsblk.atom_dfnbf).middleCols(Xsh.bfoff, Xsh.nbf);
            }
            jblk_offset += jsblk.nbf;

          }

          //----------------------------------------//

          // Now view the integrals as a jblk.nbf x (ish.nbf*Xsh.nbf) matrix, which makes
          //   the contraction more convenient.  Doesn't require any movement of data

          Eigen::Map<ThreeCenterIntContainer> g3(g3_in.data(), jblk.nbf, ish.nbf*Xsh.nbf);

          //----------------------------------------//

          // Contract to get contribution to B

          B_sd.noalias() += 2.0 * g3.transpose() * D_sd.middleCols(block_offset, jblk.nbf).transpose();

          block_offset += jblk.nbf;

        } // end loop over jblk

        /*******************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/

        mt_timer.change("K contributions", ithr);

        /*-----------------------------------------------------*/
        /* K contributions                                {{{2 */ #if 2 // begin fold

        for(BF X : function_range(Xsh)) {

          Kt.middleCols(ish.bfoff, ish.nbf).middleRows(Xsh.atom_obsbfoff, Xsh.atom_obsnbf).noalias() +=
              C_X_diff.middleRows(X.off*Xsh.atom_obsnbf, Xsh.atom_obsnbf)
              * B_sd_other.middleCols(X.off*(l_b_size+Xsh.atom_obsnbf) + Xsh.atom_obsnbf, l_b_size).transpose();

          Eigen::Map<RowMatrix> C_X_view(
              coefs_X_nu.at(Xsh.center).row(X.bfoff_in_atom).data(),
              Xsh.atom_obsnbf, nbf
          );
          for(SH lsh : iter_shells_on_center(gbs_, Xsh.center, dfbs_)) {
            block_offset = 0;
            for(ShellBlockData<> kblk : sig_partner_blocks_[lsh][NoRestrictions]) {
              Kt.middleCols(ish.bfoff, ish.nbf).middleRows(kblk.bfoff, kblk.nbf).noalias() +=
                  C_X_view.middleRows(lsh.bfoff_in_atom, lsh.nbf).middleCols(kblk.bfoff, kblk.nbf).transpose()
                  * B_sd_other.middleCols(X.off*(l_b_size+Xsh.atom_obsnbf) + lsh.bfoff_in_atom, lsh.nbf).transpose();
              block_offset += kblk.nbf;
            } // end loop over significant partners
          } // end loop over lsh

        } // end loop over X functions

        /*******************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/


      } // end loop over Xsh

      mt_timer.exit(ithr);

    } // end loop over ish

    // Deallocate buffers
    delete[] ints_buffer;
    delete[] rec_coefs_data;
    delete[] B_sd_data;
    delete[] D_sd_data;
    delete[] C_X_diff_data;

  }); // end of do_threaded

  mt_timer.exit();
  timer.insert(mt_timer);

  timer.exit();

  L_3.clear();
  L_B.clear();

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/


  /*=======================================================================================*/
  /* Global sum K tilde and transfer to back a RefSCMatrix                            {{{1 */ #if 1 // begin fold
  /*---------------------------------------------------------------------------------------*/

  // Perform the global sum
  timer.enter("global sum");
  scf_grp_->sum(Kt.data(), nbf*nbf);
  timer.exit();


  //----------------------------------------//
  // Symmetrize K
  Eigen::MatrixXd K(nbf, nbf);
  K = Kt + Kt.transpose();

  //----------------------------------------//
  // Transfer K back into a RefSCMatrix object
  Ref<Integral> localints = integral()->clone();
  localints->set_basis(gbs_);
  Ref<PetiteList> pl = localints->petite_list();
  RefSCDimension obsdim = pl->AO_basisdim();
  RefSCMatrix result(
      obsdim,
      obsdim,
      gbs_->so_matrixkit()
  );
  result.assign(K.data());

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/


  /*=======================================================================================*/
  /* Clean up                                             		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  delete[] D_data;
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/

  return result;
}

#undef SH
#undef SHV
#undef BF
#undef BLK

