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


void CADFCLHF::ScreeningStatistics::print_summary(
    std::ostream& out,
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbs,
    int print_level_in
) const
{
  using std::endl;
  using std::setw;
  int print_lvl = print_level;
  if(print_level_in != -1) {
    print_lvl = print_level_in;
  }
  const auto& old_loc = out.getloc(); out.imbue(std::locale(""));
  out << indent << "CADFCLHF Screening Statistics" << endl;
  out << indent << "-----------------------------" << endl;
  const count_t total_3c = basis->nshell() * basis->nshell() * dfbs->nshell();
  const count_t total_3c_fxn = basis->nbasis() * basis->nbasis() * dfbs->nbasis();
  out << indent << "Total shell triplets: " << total_3c << endl
      << indent << "Total function triplets: " << total_3c_fxn << endl
      << indent << "Total shell triplets after Schwarz screening: "
      << sig_pairs.load() * count_t(dfbs->nshell()) << endl
      << indent << "Total function triplets after Schwarz screening: "
      << sig_pairs_fxn.load() * count_t(dfbs->nbasis()) << endl;
  int iteration = 1;
  out << incindent;
  out << indent << setw(38) << " "
      << setw(22) << std::internal <<  "Shell-wise"
      << setw(22) << std::internal <<  "Function-wise"
      << endl;
  out << decindent;

  for(auto& iter : iterations) {
    if(print_lvl > 1) {
      out << indent << "Iteration " << iteration++ << ":" << endl;
      auto pr_iter_stat = [&](
          const std::string& title,
          const count_t sh_num,
          const count_t fxn_num
      ) {
        out << indent << setw(38) << std::left << title
            << setw(14) << std::right << sh_num
            << setw(7) << scprintf(" (%3.1f", 100.0 * double(sh_num)/double(total_3c)) << "%)"
            << setw(14) << std::right << fxn_num
            << setw(7) << scprintf(" (%3.1f", 100.0 * double(fxn_num)/double(total_3c_fxn)) << "%)"
            << endl;
      };
      out << incindent;
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
      out << decindent;
    }
  }
  out.imbue(old_loc);
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
      }
    }
  }

  return parent;
}




