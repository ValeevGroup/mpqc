//
// count_ints.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Jul 16, 2014
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

#include "cadfclhf.h"

using namespace sc;
using namespace std;

void
CADFCLHF::count_ints()
{
  /*=======================================================================================*/
  /* Setup                                                 		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  const int me = scf_grp_->me();
  const int n_node = scf_grp_->n();
  const int nbf = gbs_->nbasis();
  const int nsh = gbs_->nshell();
  const int dfnbf = dfbs_->nbasis();
  const int dfnsh = dfbs_->nshell();
  const int natom = gbs_->ncenter();
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/

  if(n_node > 1) {
    throw FeatureNotImplemented("Parallel integral counting", __FILE__, __LINE__, class_desc());
  }

  std::mutex histogram_mtx;

  do_threaded(nthread_, [&](int ithr) {

    // Make local copies
    std::vector<cadf::Histogram2d> my_dist_hists = iter_stats_->distance_hists;
    std::vector<cadf::Histogram2d> my_dist_ns_hists = iter_stats_->distance_noschwarz_hists;
    std::vector<cadf::Histogram2d> my_int_hists = iter_stats_->values_hists;
    std::vector<cadf::Histogram2d> my_exp_hists = iter_stats_->exponent_ratio_hists;
    std::vector<ull> counts = iter_stats_->int_am_counts;
    std::vector<double> sums = iter_stats_->int_am_ratio_sums;
    std::vector<double> lsums = iter_stats_->int_am_ratio_log_sums;



    for(auto&& ish : shell_range(gbs_, dfbs_)) {
      for(auto&& jsh : shell_range(gbs_, dfbs_)) {
        for(auto&& Xsh : shell_range(dfbs_, gbs_)) {
          if((ish*nsh*dfnsh+jsh*nsh+Xsh) % nthread_ != ithr) continue;
          const double distance_factor = get_distance_factor(ish, jsh, Xsh);
          const double estimate = schwarz_frob_(ish, jsh) * schwarz_df_(Xsh) * distance_factor;
          auto actual_ints = ints_to_eigen(ish, jsh, Xsh,
              eris_3c_[ithr], coulomb_oper_type_,
              gbs_, gbs_, dfbs_
          );
          const double actual = actual_ints->norm();
          const double R = get_R(ish, jsh, Xsh, true);
          //const double min_pair_expon = effective_pair_exponents_.at({ish, jsh});
          const double min_pair_expon = effective_df_exponents_[Xsh];

          if(estimate > count_ints_exclude_thresh_ or actual > count_ints_exclude_thresh_) {
            if(count_ints_use_norms_) {
              const double ratio = fabs(estimate/actual);
              if(ratio > count_ints_hist_min_ratio_ and ratio < count_ints_hist_max_ratio_) {
                sums[Xsh.am] += ratio;
                lsums[Xsh.am] += log10(ratio);
                ++counts[Xsh.am];
              }
              if(count_ints_histogram_) {
                if(actual > count_ints_exclude_thresh_) {
                  my_dist_hists[Xsh.am].insert_value(get_R(ish, jsh, Xsh), ratio);
                  if(distance_factor < 1.0) {
                    my_dist_ns_hists[Xsh.am].insert_value(R, ratio);
                    my_exp_hists[Xsh.am].insert_value(min_pair_expon, ratio);
                  }
                }
                my_int_hists[Xsh.am].insert_value(estimate, actual);
              }
              else {
                iter_stats_->int_screening_values.mine(ithr).push_back(estimate);
                iter_stats_->int_actual_values.mine(ithr).push_back(actual);
                iter_stats_->int_distance_factors.mine(ithr).push_back(distance_factor);
                iter_stats_->int_distances.mine(ithr).push_back(R);
                iter_stats_->int_ams.mine(ithr).push_back(std::make_tuple(ish.am, jsh.am, Xsh.am));
                iter_stats_->int_indices.mine(ithr).push_back(std::make_tuple(ish, jsh, Xsh));
              }
            }
            else {

              for(auto&& ibf : function_range(ish)) {
                for(auto&& jbf : function_range(jsh)) {
                  for(auto&& Xbf : function_range(Xsh)) {
                    const double act = (*actual_ints)(ibf.off*jsh.nbf + jbf.off, Xbf.off);
                    if(count_ints_histogram_) {
                      if(fabs(act) > count_ints_exclude_thresh_) {
                        my_dist_hists[Xsh.am].insert_value(R, fabs(estimate/act));
                        my_int_hists[Xsh.am].insert_value(estimate, fabs(act));
                        if(distance_factor < 1.0) {
                          my_dist_ns_hists[Xsh.am].insert_value(R, fabs(estimate/act));
                          my_exp_hists[Xsh.am].insert_value(min_pair_expon, fabs(estimate/act));
                        }

                      }
                    }
                    else {
                      iter_stats_->int_screening_values.mine(ithr).push_back(estimate);
                      iter_stats_->int_actual_values.mine(ithr).push_back(act);
                      iter_stats_->int_distance_factors.mine(ithr).push_back(distance_factor);
                      iter_stats_->int_distances.mine(ithr).push_back(R);
                      iter_stats_->int_ams.mine(ithr).push_back(std::make_tuple(ish.am, jsh.am, Xsh.am));
                    }
                  }
                }
              }

            }


          }
        }
      }
    }
    //----------------------------------------//
    if(count_ints_histogram_) {
      std::lock_guard<std::mutex> lg(histogram_mtx);
      for(int am = 0; am < dfbs_->max_angular_momentum() + 1; ++am) {
        iter_stats_->distance_hists[am].accumulate(my_dist_hists[am]);
        iter_stats_->distance_noschwarz_hists[am].accumulate(my_dist_ns_hists[am]);
        iter_stats_->values_hists[am].accumulate(my_int_hists[am]);
        iter_stats_->exponent_ratio_hists[am].accumulate(my_exp_hists[am]);
        iter_stats_->int_am_counts[am] += counts[am];
        iter_stats_->int_am_ratio_sums[am] += sums[am];
        iter_stats_->int_am_ratio_log_sums[am] += lsums[am];

      }
    }


  });




}


