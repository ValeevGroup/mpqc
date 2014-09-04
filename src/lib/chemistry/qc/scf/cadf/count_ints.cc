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

namespace {
  struct EstimatedIntegralValue {
      double ratio = 0.0;
      double est_value = 0.0;
      double act_value = 0.0;
      int ish = -1;
      int jsh = -1;
      int Xsh = -1;
      double R = 0.0;

      EstimatedIntegralValue() = delete;

      EstimatedIntegralValue(
          double ratio
      ) : ish(-1), jsh(-1), Xsh(-1), est_value(0.0), act_value(0.0), ratio(ratio), R(0.0)
      { }

      EstimatedIntegralValue(
          int ish, int jsh, int Xsh, double est_value, double act_value, double ratio, double R
      ) : ish(ish), jsh(jsh), Xsh(Xsh), est_value(est_value), act_value(act_value), ratio(ratio), R(R)
      { }

      bool operator< (const EstimatedIntegralValue& other) const {
        if(ratio < other.ratio) return true;
        else if(ratio > other.ratio) return false;
        else if(est_value < other.est_value) return true;
        else if(est_value > other.est_value) return false;
        else if(R < other.R) return true;
        else return false;
      }

      bool operator== (const EstimatedIntegralValue& other) const {
        return ratio == other.ratio
            and est_value == other.est_value
            and act_value == other.act_value
            and R == other.R;
      }

      bool operator> (const EstimatedIntegralValue& other) const {
        return not (*this < other) and not (*this == other);
      }
  };
}

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

  //if(n_node > 1) {
  //  throw FeatureNotImplemented("Parallel integral counting", __FILE__, __LINE__, class_desc());
  //}

  std::mutex histogram_mtx;


  const int n_min_max = count_ints_n_integral_extrema_;
  std::vector<std::vector<EstimatedIntegralValue>> all_mins(
      dfbs_->max_angular_momentum()+1, std::vector<EstimatedIntegralValue>(n_min_max*nthread_, std::numeric_limits<double>::infinity()));
  std::vector<std::vector<EstimatedIntegralValue>> all_maxes(
      dfbs_->max_angular_momentum()+1, std::vector<EstimatedIntegralValue>(n_min_max*nthread_, 0));
  std::vector<double> all_sqsums(dfbs_->max_angular_momentum()+1, 0.0);
  std::vector<double> all_sqlsums(dfbs_->max_angular_momentum()+1, 0.0);

  do_threaded(nthread_, [&](int ithr) {

    // Make local copies
    std::vector<cadf::Histogram2d> my_dist_hists = iter_stats_->distance_hists;
    std::vector<cadf::Histogram2d> my_dist_ns_hists = iter_stats_->distance_noschwarz_hists;
    std::vector<cadf::Histogram2d> my_int_hists = iter_stats_->values_hists;
    std::vector<cadf::Histogram2d> my_exp_hists = iter_stats_->exponent_ratio_hists;

    std::vector<ull> counts(dfbs_->max_angular_momentum()+1, 0);
    std::vector<double> sums(dfbs_->max_angular_momentum()+1, 0.0);
    std::vector<double> lsums(dfbs_->max_angular_momentum()+1, 0.0);
    std::vector<double> sqsums(dfbs_->max_angular_momentum()+1, 0.0);
    std::vector<double> sqlsums(dfbs_->max_angular_momentum()+1, 0.0);

    std::vector<std::set<EstimatedIntegralValue>> max_ratios(dfbs_->max_angular_momentum()+1);
    std::vector<std::set<EstimatedIntegralValue, std::greater<EstimatedIntegralValue>>> min_ratios(dfbs_->max_angular_momentum()+1);
    //    std::vector<EstimatedIntegralValue>(n_min_max, std::numeric_limits<double>::infinity())
    //);
    for(int iam = 0; iam < dfbs_->max_angular_momentum()+1; ++iam) {
      max_ratios[iam].emplace(0);
      min_ratios[iam].emplace(std::numeric_limits<double>::infinity());
      //counts[iam] = 0;
      //sums[iam] = 0.0;
      //lsums[iam] = 0.0;
    }

    for(auto&& ish : shell_range(gbs_, dfbs_)) {
      for(auto&& jsh : shell_range(gbs_, dfbs_)) {
        if(count_ints_exclude_two_center_ and ish.center == jsh.center) continue;
        for(auto&& Xsh : shell_range(dfbs_, gbs_)) {
          if((ish*nsh*dfnsh+jsh*nsh+Xsh) % nthread_ != ithr) continue;
          if(((ish*nsh*dfnsh+jsh*nsh+Xsh) / nthread_) % n_node != me) continue;
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
              if(distance_factor < 1.0
                  and ratio < count_ints_hist_max_ratio_ and ratio > count_ints_hist_min_ratio_) {
                sums[Xsh.am] += ratio;
                const double lratio = log10(ratio);
                lsums[Xsh.am] += lratio;
                counts[Xsh.am] += 1;
                sqsums[Xsh.am] += ratio*ratio;
                sqlsums[Xsh.am] += lratio*lratio;
              }
              if(count_ints_histogram_) {
                if(actual > count_ints_exclude_thresh_) {
                  my_dist_hists[Xsh.am].insert_value(R, ratio);
                  if(distance_factor < 1.0) {
                    my_dist_ns_hists[Xsh.am].insert_value(R, ratio);
                    my_exp_hists[Xsh.am].insert_value(min_pair_expon, ratio);
                  }
                }
                my_int_hists[Xsh.am].insert_value(estimate, actual);

                if(distance_factor < 1.0) {
                  // Figure out if it goes in the min list
                  EstimatedIntegralValue ival(ish, jsh, Xsh, estimate, actual, ratio, R);
                  if(ival < *(min_ratios[(int)Xsh.am].begin())) {
                    if(min_ratios[(int)Xsh.am].size() >= n_min_max) {
                      min_ratios[(int)Xsh.am].erase(min_ratios[(int)Xsh.am].begin());
                    }
                    min_ratios[(int)Xsh.am].insert(ival);
                  }
                  if(ival > *(max_ratios[(int)Xsh.am].begin())) {
                    if(max_ratios[(int)Xsh.am].size() >= n_min_max) {
                      max_ratios[(int)Xsh.am].erase(max_ratios[(int)Xsh.am].begin());
                    }
                    max_ratios[(int)Xsh.am].insert(ival);
                  }
                }

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
    {
      std::lock_guard<std::mutex> lg(histogram_mtx);
      if(count_ints_histogram_) {
        for(int am = 0; am < dfbs_->max_angular_momentum() + 1; ++am) {
          iter_stats_->distance_hists[am].accumulate(my_dist_hists[am]);
          iter_stats_->distance_noschwarz_hists[am].accumulate(my_dist_ns_hists[am]);
          iter_stats_->values_hists[am].accumulate(my_int_hists[am]);
          iter_stats_->exponent_ratio_hists[am].accumulate(my_exp_hists[am]);
          iter_stats_->int_am_counts[am] += counts[am];
          iter_stats_->int_am_ratio_sums[am] += sums[am];
          iter_stats_->int_am_ratio_log_sums[am] += lsums[am];
          all_sqsums[am] += sqsums[am];
          all_sqlsums[am] += sqlsums[am];
        }
      }

      for(int iam = 0; iam < dfbs_->max_angular_momentum()+1; ++iam) {
        std::copy(max_ratios[iam].begin(), max_ratios[iam].end(), all_maxes[iam].begin() + ithr*n_min_max);
        std::copy(min_ratios[iam].begin(), min_ratios[iam].end(), all_mins[iam].begin() + ithr*n_min_max);
      }
    }

  });

  auto& o = ExEnv::out0();
  auto print_valdata = [&](const EstimatedIntegralValue& v, int ilist) {
      o << setw(2) << std::right << ilist << ". (" << v.ish << " " << v.jsh << " | " << v.Xsh << " )" << endl;
      o << "    Ratio: " << setw(10) << v.ratio
        << ",  Actual: " << v.act_value << ",  Estimate: " << v.est_value << endl;
      o << "    Separation: " << v.R << " ( > " << pair_extents_.at({v.ish, v.jsh}) + df_extents_[v.Xsh] << " )" << endl;
      ShellData ish(v.ish, gbs_);
      o << "    Shell " << v.ish << " on center " << ish.center << " (Z = " << molecule()->atom(ish.center).Z() << "), am = " << ish.am << endl;
      for(int i = 0; i < gbs_->shell(ish).nprimitive(); ++i) {
        o << "       coef: "
          << gbs_->shell(ish).coefficient_unnorm(0, i)
          << ", expon: "
          << gbs_->shell(ish).exponent(i)
          << endl;
      }
      ShellData jsh(v.jsh, gbs_);
      o << "    Shell " << v.jsh << " on center " << jsh.center << " (Z = " << molecule()->atom(jsh.center).Z() << "), am = " << jsh.am << endl;
      for(int i = 0; i < gbs_->shell(jsh).nprimitive(); ++i) {
        o << "       coef: "
          << gbs_->shell(jsh).coefficient_unnorm(0, i)
          << ", expon: "
          << gbs_->shell(jsh).exponent(i)
          << endl;
      }
      ShellData Xsh(v.Xsh, dfbs_);
      o << "    Shell " << v.Xsh << " on center " << Xsh.center << " (Z = " << molecule()->atom(Xsh.center).Z() << ")" << endl;
      for(int i = 0; i < dfbs_->shell(Xsh).nprimitive(); ++i) {
        o << "       coef: "
          << dfbs_->shell(Xsh).coefficient_unnorm(0, i)
          << ", expon: "
          << dfbs_->shell(Xsh).exponent(i)
          << endl;
      }
      o << "    Pair ( " << v.ish << " " << v.jsh << " | extent: " << pair_extents_.at({v.ish, v.jsh}) << endl;
      o << "    | " << v.Xsh << " ) extent: " << df_extents_[Xsh] << endl;
      o << "    Schwarz ( " << v.ish << " " << v.jsh << " | " << v.ish << " " << v.jsh << " )^(1/2): " << double(schwarz_frob_(v.ish, v.jsh)) << endl;
      if(dist_factor_use_overlap_) {
        o << "    Overlap S_{ " << v.ish << " " << v.jsh << " } " << double(S_frob_(v.ish, v.jsh)) << endl;
        o << "    Ratio S/Schwarz: " << double(S_frob_(v.ish, v.jsh) / schwarz_frob_(v.ish, v.jsh)) << endl;
      }
  };

  if(n_node == 1) {
    o << "==================================" << endl;
    o << "= Integral estimation statistics =" << endl;
    o << "==================================" << endl;
    for(int iam = 0; iam < dfbs_->max_angular_momentum()+1; ++iam) {
      o << "For l_X = " << iam << ":" << endl;
      o << n_min_max << " most overestimated integral shell triplets:" << endl;
      int ilist = 1;
      std::sort(all_maxes[iam].begin(), all_maxes[iam].end(), std::greater<EstimatedIntegralValue>());
      for(auto&& v : all_maxes[iam]) {
        if(v.ish == -1) break;
        if(ilist > n_min_max) break;
        print_valdata(v, ilist);
        ++ilist;
      }
      o << n_min_max << " most underestimated integral shell triplets:" << endl;
      ilist = 1;
      std::sort(all_mins[iam].begin(), all_mins[iam].end());
      for(auto&& v : all_mins[iam]) {
        if(v.ish == -1) break;
        if(ilist > n_min_max) break;
        print_valdata(v, ilist);
        ++ilist;
      }
      o << "----------------------------------" << endl;

    }
    o << "==================================" << endl;
  }
  else {
    for(int iam = 0; iam < dfbs_->max_angular_momentum()+1; ++iam) {
      std::sort(all_maxes[iam].begin(), all_maxes[iam].end(), std::greater<EstimatedIntegralValue>());
      std::sort(all_mins[iam].begin(), all_mins[iam].end());
    }
  }

  o << "\n\n===============================" << endl;
  o << "= Integral statistics summary =" << endl;
  o << "===============================" << endl;

  auto& out = ExEnv::out0();

  double r_sum_tot = 0.0;
  double r_lsum_tot = 0.0;
  double r_sqsum_tot = 0.0;
  double r_sqlsum_tot = 0.0;
  ull r_count_tot = 0.0;
  double max_ratio_tot = 0.0;
  double min_ratio_tot = std::numeric_limits<double>::infinity();
  for(int l = 0; l < dfbs_->max_angular_momentum() + 1; ++l) {
    double r_min = all_mins[l].begin()->ratio;
    scf_grp_->min(r_min);
    double r_max = all_maxes[l].begin()->ratio;
    scf_grp_->max(r_max);
    scf_grp_->sum(iter_stats_->int_am_ratio_sums[l]);
    const double r_sum = iter_stats_->int_am_ratio_sums[l];
    scf_grp_->sum(iter_stats_->int_am_ratio_log_sums[l]);
    const double r_lsum = iter_stats_->int_am_ratio_log_sums[l];
    scf_grp_->sum(all_sqsums[l]);
    const double r_sqsum = all_sqsums[l];
    scf_grp_->sum(all_sqlsums[l]);
    const double r_sqlsum = all_sqlsums[l];
    long r_count = iter_stats_->int_am_counts[l];
    scf_grp_->sum(&r_count, 1);
    out << indent << "Ratio statistics for l_X = " << l << ":" << endl;
    out << indent << "  Average: " << r_sum / double(r_count) << endl;
    out << indent << "  Average log: " << r_lsum / double(r_count) << endl;
    out << indent << "  Count: " << r_count << endl;
    out << indent << "  Min ratio: " << double(r_min) << endl;
    out << indent << "  Max ratio: " << double(r_max) << endl;
    out << indent << "  Sum: " << setw(20) << setprecision(15) << double(r_sum) << endl;
    out << indent << "  Sum of Squares: " << setw(20) << setprecision(15) << double(r_sqsum) << endl;
    out << indent << "  Sum of Logs: " << setw(20) << setprecision(15) << double(r_lsum) << endl;
    out << indent << "  Sum of Square Logs: " << setw(20) << setprecision(15) << double(r_sqlsum) << endl;
    out << indent << "  Standard Deviation: "
        << sqrt((r_sqsum - (r_sum*r_sum/r_count))/double(r_count - 1)) << endl;
    out << indent << "  Log Standard Deviation: "
        << sqrt((r_sqlsum - (r_lsum*r_lsum/r_count))/double(r_count - 1)) << endl;
    r_sum_tot += r_sum;
    r_lsum_tot += r_lsum;
    r_sqsum_tot += r_sqsum;
    r_sqlsum_tot += r_sqlsum;
    r_count_tot += r_count;
    if(r_min < min_ratio_tot) min_ratio_tot = r_min;
    if(r_max > max_ratio_tot) max_ratio_tot = r_max;
  }
  o << "-------------------------------" << endl;
  out << indent << "Ratio statistics for all 3-center integrals:" << endl;
  out << indent << "  Average: " << r_sum_tot / double(r_count_tot) << endl;
  out << indent << "  Average log: " << r_lsum_tot / double(r_count_tot) << endl;
  out << indent << "  Count: " << r_count_tot << endl;
  out << indent << "  Min ratio: " << double(min_ratio_tot) << endl;
  out << indent << "  Max ratio: " << double(max_ratio_tot) << endl;
  out << indent << "  Sum: " << setw(20) << setprecision(15) << double(r_sum_tot) << endl;
  out << indent << "  Sum of Squares: " << setw(20) << setprecision(15) << double(r_sqsum_tot) << endl;
  out << indent << "  Sum of Logs: " << setw(20) << setprecision(15) << double(r_lsum_tot) << endl;
  out << indent << "  Sum of Square Logs: " << setw(20) << setprecision(15) << double(r_sqlsum_tot) << endl;
  out << indent << "  Standard Deviation: "
      << sqrt((r_sqsum_tot - (r_sum_tot*r_sum_tot/r_count_tot))/double(r_count_tot - 1)) << endl;
  out << indent << "  Log Standard Deviation: "
      << sqrt((r_sqlsum_tot - (r_lsum_tot*r_lsum_tot/r_count_tot))/double(r_count_tot - 1)) << endl;
  o << "===============================" << endl;

}


