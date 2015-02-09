//
// statistics.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Feb 25, 2014
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

#include <util/misc/xmlwriter.h>

using namespace sc;
using namespace std;
using std::endl;
using std::setw;

cadf::Histogram2d::Histogram2d(
    int nbins_h, double min_h, double max_h,
    int nbins_v, double min_v, double max_v,
    bool log_h, bool log_v,
    bool clip_edges
) : nbins_h_(nbins_h), nbins_v_(nbins_v),
    //hist_mat_(nbins_v_, nbins_h_),
    min_h_(min_h), max_h_(max_h), min_v_(min_v), max_v_(max_v),
    log_h_(log_h), log_v_(log_v),
    clip_edges_(clip_edges)
{
  hist_mat_.resize(nbins_v_, nbins_h_);
  if(log_h_) {
    min_h_ = log10(min_h_);
    max_h_ = log10(max_h_);
  }
  interval_h_ = (max_h_ - min_h_) / nbins_h_;

  if(log_v_) {
    min_v_ = log10(min_v_);
    max_v_ = log10(max_v_);
  }
  interval_v_ = (max_v_ - min_v_) / nbins_v_;

  hist_mat_.setZero();
}

void
cadf::Histogram2d::insert_value(
    double value_h,
    double value_v
)
{
  double hval = value_h;
  if(log_h_) hval = log10(hval);
  double vval = value_v;
  if(log_v_) vval = log10(vval);
  if(clip_edges_ and (hval < min_h_ or hval >= max_h_ or vval < min_v_ or vval >= max_v_)) {
    return;
  }

  int bin_h = int((hval - min_h_) / interval_h_);
  int bin_v = int((vval - min_v_) / interval_v_);

  // This should only happen if clip_edges_ is false
  if(bin_h >= nbins_h_) bin_h = nbins_h_ - 1;
  if(bin_v >= nbins_v_) bin_v = nbins_v_ - 1;
  if(bin_h < 0) bin_h = 0;
  if(bin_v < 0) bin_v = 0;

  ++hist_mat_(bin_v, bin_h);
}


void CADFCLHF::ScreeningStatistics::print_summary(
    std::ostream& out,
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbs,
    const CADFCLHF* parent,
    int print_level_in, bool new_exchange
) const
{
  int print_lvl = print_level;
  if(print_level_in != -1) {
    print_lvl = print_level_in;
  }

  //if(parent->count_ints_only_) {
  //  for(int l = 0; l < parent->dfbs_->max_angular_momentum() + 1; ++l) {
  //    const double r_sum = iterations[0].int_am_ratio_sums[l];
  //    const double r_lsum = iterations[0].int_am_ratio_log_sums[l];
  //    const ull r_count = iterations[0].int_am_counts[l];
  //    out << indent << "Ratio statistics for l_X = " << l << ":" << endl;
  //    out << indent << "  Average: " << r_sum / double(r_count) << endl;
  //    out << indent << "  Average log: " << r_lsum / double(r_count) << endl;
  //    out << indent << "  Count: " << double(r_count) << endl;
  //  }
  //}

  out << indent << "CADFCLHF Screening Statistics" << endl;
  out << indent << "-----------------------------" << endl;
  const count_t total_3c = ((uint64_t)basis->nshell())
      * ((uint64_t)basis->nshell())
      * ((uint64_t)dfbs->nshell());
  const count_t total_3c_fxn = ((uint64_t)basis->nbasis())
      * ((uint64_t)basis->nbasis()) * ((uint64_t)dfbs->nbasis());
  out << indent << "Total shell triplets: " << total_3c << endl
      << indent << "Total function triplets: " << total_3c_fxn << endl
      << indent << "Total shell triplets after Schwarz screening: "
      << sig_pairs.load() * count_t(dfbs->nshell()) << endl
      << indent << "Total function triplets after Schwarz screening: "
      << sig_pairs_fxn.load() * count_t(dfbs->nbasis()) << endl;
  int iteration = 1;
  out << incindent;
  out << indent << setw(38) << " "
      << setw(23) << std::internal <<  "Shell-wise"
      << setw(23) << std::internal <<  "Function-wise"
      << endl;
  out << decindent;

  auto pr_iter_stat = [&](
      const std::string& title,
      const count_t sh_num,
      const count_t fxn_num
  ) {
    out << indent << setw(38) << std::left << title
        << setw(14) << std::right << sh_num
        << setw(7) << std::right << scprintf("(%3.1f", 100.0 * double(sh_num)/double(total_3c)) << "%)"
        << setw(14) << std::right << fxn_num
        << setw(7) << std::right << scprintf("(%3.1f", 100.0 * double(fxn_num)/double(total_3c_fxn)) << "%)"
        << endl;
  };

  for(auto& iter : iterations) {
    out << indent << "Iteration " << iteration++ << ":" << endl;
    out << incindent;
    if(new_exchange) {
      pr_iter_stat("K: 3c ints computed", iter.K_3c_needed, iter.K_3c_needed_fxn);
      out << indent << setw(38) << std::left << "L3 build comparisons"
          << setw(14) << std::right << iter.L3_build_compares
          << setw(9) << std::right << "  (N/A)"
          << setw(14) << std::right << "-"
          << setw(9) << std::right << "  (N/A)"
          << endl;
      out << indent << setw(38) << std::left << "B contraction multiplies"
          << setw(14) << std::right << "-"
          << setw(9) << std::right << "  (N/A)"
          << setw(14) << std::right << iter.K_3c_contract_fxn
          << setw(9) << std::right << "  (N/A)"
          << endl;
      out << indent << setw(38) << std::left << "g2 contraction multiplies"
          << setw(14) << std::right << "-"
          << setw(9) << std::right << "  (N/A)"
          << setw(14) << std::right << iter.K_2c_contract_fxn
          << setw(9) << std::right << "  (N/A)"
          << endl;
      out << indent << setw(38) << std::left << "Kt contraction 1 multiplies"
          << setw(14) << std::right << "-"
          << setw(9) << std::right << "  (N/A)"
          << setw(14) << std::right << iter.Kt_contract1_fxn
          << setw(9) << std::right << "  (N/A)"
          << endl;
      out << indent << setw(38) << std::left << "Kt contraction 2 multiplies"
          << setw(14) << std::right << "-"
          << setw(9) << std::right << "  (N/A)"
          << setw(14) << std::right << iter.Kt_contract2_fxn
          << setw(9) << std::right << "  (N/A)"
          << endl;
    }
    else {
      if(print_lvl > 1) {
        pr_iter_stat("K: 3c ints needed", iter.K_3c_needed, iter.K_3c_needed_fxn);
        pr_iter_stat("K: 3c ints screened by distance", iter.K_3c_dist_screened, iter.K_3c_dist_screened_fxn);
        if(print_lvl > 2) {
          pr_iter_stat("K: 3c ints underestimated", iter.K_3c_underestimated, iter.K_3c_underestimated_fxn);
          pr_iter_stat("K: 3c ints needed, \"perfect screening\"", iter.K_3c_perfect, iter.K_3c_perfect_fxn);
          pr_iter_stat("K: \"extra\" ints computed",
              iter.K_3c_needed - iter.K_3c_perfect,
              iter.K_3c_needed_fxn - iter.K_3c_perfect_fxn
          );
        }
      }
    }

    out << decindent;

  }
}


using boost::property_tree::ptree;

ptree&
sc::write_xml(
    const CADFCLHF::ScreeningStatistics& obj,
    ptree& parent,
    const XMLWriter& writer
)
{
  ptree* child_ptr;
  if(writer.fold_in_class_name()) {
    parent.put("<xmlattr>.type", "CADFCLHF::ScreeningStatistics");
    child_ptr = &(parent);
  }
  else{
    ptree& tmp = parent.add_child("ScreeningStatistics", ptree());
    child_ptr = &tmp;
  }
  ptree& child = *child_ptr;

  int i = 0;
  bool tmp = obj.xml_stats_saved;
  for(auto&& iter : obj.iterations) {
    writer.insert_child(child, iter, "iteration", std::map<std::string, int>{
      {"number", i}
    });
    // only print xml stats for the first iteration, for now
    // TODO option to turn this behavior on or off
    obj.xml_stats_saved = false;
  }
  obj.xml_stats_saved = tmp;

  return child;
}

boost::property_tree::ptree&
sc::write_xml(
    const CADFCLHF::ScreeningStatistics::Iteration& obj,
    ptree& parent,
    const XMLWriter& writer
)
{
  typedef CADFCLHF::ScreeningStatistics::accumulate_t accumulate_t;
  ptree* child_ptr;
  if(writer.fold_in_class_name()) {
    parent.put("<xmlattr>.type", "CADFCLHF::ScreeningStatistics::Iteration");
    child_ptr = &(parent);
  }
  else{
    ptree& tmp = parent.add_child("Iteration", ptree());
    child_ptr = &tmp;
  }
  ptree& child = *child_ptr;

  auto add_stat = [](
      ptree& p,
      const std::string& name,
      const accumulate_t& sh_count,
      const accumulate_t& fxn_count
  ){
    ptree& sub = p.add_child(name, ptree());
    sub.put("shell_tuples", sh_count.load());
    sub.put("function_tuples", fxn_count.load());
  };

  auto am_stat = [](ptree& p, const std::string& name, int am, double value) {
    ptree& sub = p.add_child(name, ptree());
    sub.put("<xmlattr>.am", am);
    sub.put_value(value);
  };
  auto am_stat_ull = [](ptree& p, const std::string& name, int am, ull value) {
    ptree& sub = p.add_child(name, ptree());
    sub.put("<xmlattr>.am", am);
    sub.put_value(value);
  };

  std::string fmt = "%16.12g";
  std::string ifmt = "%5d";
  auto double_str = [&fmt](double val) -> std::string {
     return std::string(scprintf(fmt.c_str(), val).str());
  };
  auto int_str = [&ifmt](int val) -> std::string {
     return std::string(scprintf(ifmt.c_str(), val).str());
  };

  ptree& kstats = child.add_child("compute_k", ptree());
  if(obj.parent->print_level > 0) {
    add_stat(kstats, "needed", obj.K_3c_needed, obj.K_3c_needed_fxn);
    add_stat(kstats, "distance_screened", obj.K_3c_dist_screened, obj.K_3c_dist_screened_fxn);
    if(obj.parent->print_level > 2) {
      add_stat(kstats, "underestimated", obj.K_3c_underestimated, obj.K_3c_underestimated_fxn);
      add_stat(kstats, "needed_perfect", obj.K_3c_perfect, obj.K_3c_perfect_fxn);
      if(obj.parent->xml_stats_saved) {
        writer.insert_child(child, obj.int_screening_values, "int_estimates");
        writer.insert_child(child, obj.int_actual_values, "int_values");
        writer.insert_child(child, obj.int_distance_factors, "int_distance_factors");
        writer.insert_child(child, obj.int_distances, "int_distances");
        writer.insert_child(child, obj.int_indices.merged(), "int_indices");
        writer.insert_child(child, obj.int_ams.merged(), "int_ams");
      }
      for(int l = 0; l <= obj.int_am_counts.size(); ++l) {
        const double r_sum = obj.int_am_ratio_sums[l];
        const double r_lsum = obj.int_am_ratio_log_sums[l];
        const ull r_count = obj.int_am_counts[l];
        am_stat(kstats, "ratio_sum", l, r_sum);
        am_stat(kstats, "ratio_log_sum", l, r_lsum);
        am_stat(kstats, "average_ratio", l, r_sum / double(r_count));
        am_stat(kstats, "average_log_ratio", l, r_lsum / double(r_count));
        am_stat_ull(kstats, "int_count", l, r_count);
      }
      if(obj.parent->histogram_mode) {
        int am = 0;
        for(auto&& hist : obj.distance_hists) {
          writer.insert_child(child, hist.matrix(), "ratio_vs_distance_histogram",
              std::map<std::string, std::string>{
                { "am", int_str(am) },
                { "min_ratio", double_str(hist.min_v_) },
                { "max_ratio", double_str(hist.max_v_) },
                { "min_distance", double_str(hist.min_h_) },
                { "max_distance", double_str(hist.max_h_) }
          });
          ++am;
        }
        am = 0;
        for(auto&& hist : obj.distance_noschwarz_hists) {
          writer.insert_child(child, hist.matrix(), "ratio_vs_distance_noschwarz_histogram",
              std::map<std::string, std::string>{
                { "am", int_str(am) },
                { "min_ratio", double_str(hist.min_v_) },
                { "max_ratio", double_str(hist.max_v_) },
                { "min_distance", double_str(hist.min_h_) },
                { "max_distance", double_str(hist.max_h_) }
          });
          ++am;
        }
        am = 0;
        for(auto&& hist : obj.values_hists) {
          writer.insert_child(child, hist.matrix(), "actual_vs_estimate_histogram",
              std::map<std::string, std::string>{
                { "am", int_str(am) },
                { "min_value", double_str(hist.min_h_) },
                { "max_value", double_str(hist.max_h_) },
              }
          );
          ++am;
        }
        am = 0;
        for(auto&& hist : obj.exponent_ratio_hists) {
          writer.insert_child(child, hist.matrix(), "ratio_vs_exponent_histogram",
              std::map<std::string, std::string>{
                { "am", int_str(am) },
                { "min_exponent", double_str(hist.min_h_) },
                { "max_exponent", double_str(hist.max_h_) },
                { "min_ratio", double_str(hist.min_v_) },
                { "max_ratio", double_str(hist.max_v_) },
              }
          );
          ++am;
        }
      }
    }
  }



  return parent;
}




