//
// thread_timer.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Feb 7, 2014
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

#include <vector>

#include <util/misc/thread_timer.h>
#include <util/misc/regtime.h>

using namespace sc;

TimedRegion* MultiThreadTimer::collect_regions_recursive(
    const std::vector<const ThreadTimer*>& subtimers,
    const std::string& my_name,
    TimedRegion* parent
) const
{
  TimedRegion* rv = new TimedRegion(my_name.c_str());
  rv->up_ = parent;
  rv->cpu_time_ = 0;
  rv->wall_time_ = 0;

  std::vector<std::vector<bool>> dones;

  for(auto st : subtimers) {
    rv->cpu_time_ += st->read_seconds();
    rv->wall_time_ += st->read_seconds();
    dones.emplace_back(st->section_names_.size(), false);
  }
  // For now, treat wall time as an average, since it isn't well defined for a nested, threaded region
  rv->wall_time_ /= subtimers.size();

  // A couple of useful lambdas
  auto all = [](const std::vector<bool> v) -> bool {
    for(const auto& i : v){
      if(!i) return false;
    }
    return true;
  };
  auto first_false_index = [](const std::vector<bool> v) -> int {
    int idx = 0;
    for(const auto& i : v){
      if(!i) return idx;
      else ++idx;
    }
    return -1;
  };

  TimedRegion* prev_subregion = 0;
  // Now loop and get the subsections
  while(true){

    std::vector<const ThreadTimer*> next_subs;
    std::string curr_name;
    bool name_found = false;

    for(int i = 0; i < subtimers.size(); ++i) {
      const ThreadTimer* sub = subtimers[i];
      if(not all(dones[i])){
        if(not name_found){
          const int idx = first_false_index(dones[i]);
          curr_name = sub->section_names_[idx];
          dones[i][idx] = true;
          name_found = true;
          next_subs.push_back(&(sub->subtimers_.at(curr_name)));
        }
        else {
          if(sub->subtimers_.find(curr_name) != sub->subtimers_.end()) {
            int iname = 0;
            for(const auto& subname : sub->section_names_) {
              if(subname == curr_name) {
                dones[i][iname] = true;
                next_subs.push_back(&(sub->subtimers_.at(curr_name)));
                break;
              }
              else ++iname;
            }
          }
        }
      }
    } // end loop over subtimers

    if(name_found){
      TimedRegion* next_sub = collect_regions_recursive(next_subs, curr_name, rv);
      if(prev_subregion) {
        prev_subregion->next_ = next_sub;
        next_sub->prev_ = prev_subregion;
        prev_subregion = next_sub;
      }
      else {
        rv->subregions_ = next_sub;
        prev_subregion = next_sub;
      }
    }
    else{
      break;
    }
  } // end while all_done

  return rv;
}


TimedRegion*
MultiThreadTimer::make_timed_region() const
{
  std::vector<const ThreadTimer*> subtimers;
  for(const auto& tim : thread_timers_) subtimers.push_back(&tim);
  return collect_regions_recursive(subtimers, name_, 0);
}
