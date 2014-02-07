//
// thread_timer.h
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

#ifndef _util_misc_thread_timer_h
#define _util_misc_thread_timer_h

// standard library includes
#include <atomic>
#include <chrono>

// boost includes
#include <boost/thread/thread.hpp>

// MPQC includes
#include <util/misc/formio.h>

namespace sc {

template<
  typename DurationType = std::chrono::nanoseconds,
  typename ClockType = std::chrono::high_resolution_clock,
  typename AccumulateToType = std::atomic_uint_fast64_t
>
class auto_time_accumulator {

  public:

    typedef DurationType duration_type;
    typedef ClockType clock_type;
    typedef AccumulateToType accumulate_to_type;
    typedef auto_time_accumulator<
        duration_type, clock_type, accumulate_to_type
    > self_type;

    auto_time_accumulator() = delete;
    auto_time_accumulator(self_type&) = delete;
    auto_time_accumulator(const self_type&) = delete;
    void* operator new(size_t) = delete;

    auto_time_accumulator(accumulate_to_type& dest)
      : dest_(dest)
    {
      start_ = clock_type::now();
    }

    auto_time_accumulator(self_type&& other) = default;

    ~auto_time_accumulator()
    {
      auto end = clock_type::now();
      dest_ += std::chrono::duration_cast<DurationType>(end - start_).count();
    }

  private:
    typename clock_type::time_point start_;
    accumulate_to_type& dest_;
};

template<
  typename AccumulateToType = std::atomic_uint_fast64_t
>
auto_time_accumulator<
  std::chrono::nanoseconds,
  std::chrono::high_resolution_clock,
  AccumulateToType
>
make_auto_timer(AccumulateToType& dest) {
  return auto_time_accumulator<
    std::chrono::nanoseconds,
    std::chrono::high_resolution_clock,
    AccumulateToType
  >(dest);
}

// TODO doesn't work with nesting!!!
template<
  typename DurationType = std::chrono::nanoseconds,
  typename ClockType = std::chrono::high_resolution_clock,
  typename AccumulateToType = std::atomic_uint_fast64_t
>
class time_accumulator_factory {

  public:

    typedef DurationType duration_type;
    typedef ClockType clock_type;
    typedef AccumulateToType accumulate_to_type;
    typedef time_accumulator_factory<
        duration_type, clock_type, accumulate_to_type
    > self_type;
    typedef auto_time_accumulator<
        duration_type, clock_type, accumulate_to_type
    > generated_type;
    typedef decltype(accumulate_to_type().load()) accumulated_value_type;

    time_accumulator_factory() = delete;

    explicit time_accumulator_factory(accumulate_to_type& dest)
      : dest_(dest)
    { }

    generated_type create() const {
      return generated_type(dest_);
    }

    accumulated_value_type total_time() const {
      return dest_.load();
    }

  private:

    accumulate_to_type& dest_;

};

class MultiThreadTimer;

class ThreadTimer {

  public:

    typedef std::unordered_map<std::string, ThreadTimer> section_map;
    typedef typename time_accumulator_factory<>::clock_type clock_type;
    typedef std::chrono::time_point<clock_type> time_type;
    typedef std::chrono::nanoseconds duration_type;
    typedef std::chrono::duration<double> fp_seconds;

  private:

    time_type begin_time_;
    duration_type accum_time_{ 0 };

    std::vector<std::string> section_names_;
    section_map subtimers_;

    ThreadTimer* active_subsection_;
    std::string active_subname_ = "";
    bool stopped_{ true };

    int depth_;

    void start() {
      assert(stopped_);
      stopped_ = false;
      begin_time_ = clock_type::now();
    }

    void stop() {
      assert(!stopped_);
      accum_time_ += std::chrono::duration_cast<duration_type>(
          clock_type::now() - begin_time_
      );
      stopped_ = true;
    }

    struct Holdable {
        ThreadTimer* to_hold;
        ThreadTimer* parent;
        Holdable(ThreadTimer* to_hold, ThreadTimer* parent)
          : to_hold(to_hold), parent(parent)
        { }
    };

  public:

    ThreadTimer() = delete;

    explicit ThreadTimer(int depth, bool start=true)
      : section_names_(0),
        subtimers_(0),
        active_subsection_(0),
        depth_(depth)
    {
      if(start) this->start();
    }

    Holdable get_subtimer(const std::string& subname, bool start=false) {
      if(active_subsection_) {
        return active_subsection_->get_subtimer(subname, start);
      }
      else {
        ThreadTimer* rv_ptr;
        active_subname_ = subname;
        auto subspot = subtimers_.find(subname);
        if(subspot != subtimers_.end()) {
          rv_ptr = &(subspot->second);
          if(start) {
            active_subsection_ = rv_ptr;
            rv_ptr->start();
          }
        }
        else {
          auto insertion_pair = subtimers_.emplace(
              std::piecewise_construct,
              std::forward_as_tuple(subname),
              std::forward_as_tuple(depth_+1, start)
          );
          section_names_.push_back(subname);
          assert(insertion_pair.second);
          rv_ptr = &(insertion_pair.first->second);
          if(start) {
            active_subsection_ = rv_ptr;
          }
        }
        return ThreadTimer::Holdable(rv_ptr, this);
      }
    }

    void enter(const std::string& subname) {
      get_subtimer(subname, true);
    }

    void exit() {
      // TODO This should throw exceptions on failure rather than just asserting
      if(active_subsection_) {
        active_subsection_->exit();
        if(active_subsection_->stopped_) {
          active_subsection_ = 0;
        }
      }
      else{
        this->stop();
      }
    }

    void change(const std::string& newsub) {
      this->exit();
      enter(newsub);
    }

    bool is_stopped() const { return stopped_; }

    double read_seconds() const {
      assert(stopped_);
      return fp_seconds(accum_time_).count();
    }

    friend class MultiThreadTimer;
    friend class TimerHolder;

};

class TimerHolder {

    ThreadTimer* held;
    ThreadTimer* parent;

  public:

    explicit TimerHolder(const ThreadTimer::Holdable& to_hold)
      : held(to_hold.to_hold), parent(to_hold.parent) {
      held->start();
      parent->active_subsection_ = held;
    }

    ~TimerHolder()
    {
      held->stop();
      parent->active_subsection_ = 0;
    }

    void change(ThreadTimer::Holdable& other) {
      assert(other.parent == parent);
      held->stop();
      held = other.to_hold;
      parent->active_subsection_ = held;
      held->start();
    }

};

class MultiThreadTimer {

    std::vector<ThreadTimer> thread_timers_;
    int nthreads_;
    std::string name_;

    typename time_accumulator_factory<>::accumulate_to_type overhead_nanos_{ 0 };
    time_accumulator_factory<> overhead_factory_{ overhead_nanos_ };

    boost::thread::id creator_id_;

    void print_sub(
        std::ostream& out,
        int indent_size,
        const std::vector<const ThreadTimer*>& subtimers,
        const std::string& name,
        int label_width
    ) {
      double sum = 0.0;
      double min = std::numeric_limits<double>::infinity();
      double max = 0.0;
      for(auto timer : subtimers) {
        // TODO throw exception rather than just asserting
        assert(timer->is_stopped());
        const double time = timer->read_seconds();
        sum += time;
        if(time < min) min = time;
        if(time > max) max = time;
      }
      const double avg = sum / (double)subtimers.size();
      const std::string indent(indent_size, ' ');
      out << std::setw(label_width) << std::left << (indent + name + ":")
          << scprintf("%7.2f %7.2f %7.2f",
               avg, min, max
             )
          << std::endl;
      //----------------------------------------//
      std::vector<std::vector<bool>> dones;
      for(auto st : subtimers) { dones.emplace_back(st->section_names_.size(), false); }
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
          print_sub(out, indent_size+2, next_subs, curr_name, label_width);
        }
        else{
          break;
        }
      } // end while all_done

    }

  public:

    MultiThreadTimer(const std::string& name, int nthreads)
      : name_(name),
        nthreads_(nthreads),
        creator_id_(boost::this_thread::get_id()),
        overhead_nanos_(0),
        overhead_factory_(overhead_nanos_)
    {
      for(int i = 0; i < nthreads_; ++i) {
        thread_timers_.emplace_back(0);
      }
    }

    void enter(const std::string& subname, int ithr) {
      auto overtime = overhead_factory_.create();
      thread_timers_[ithr].enter(subname);
    }

    void exit(int ithr) {
      auto overtime = overhead_factory_.create();
      thread_timers_[ithr].exit();
    }

    ThreadTimer::Holdable
    get_subtimer(const std::string& subname, int ithr) {
      auto overtime = overhead_factory_.create();
      return thread_timers_[ithr].get_subtimer(subname);
    }

    void exit() {
      auto overtime = overhead_factory_.create();
      const boost::thread::id& my_id = boost::this_thread::get_id();
      assert(my_id == creator_id_);
      for(auto& tim : thread_timers_) tim.exit();
    }

    void change(const std::string& subname, int ithr) {
      auto overtime = overhead_factory_.create();
      thread_timers_[ithr].change(subname);
    }

    void print(
        std::ostream& out=ExEnv::out0(),
        int indent_size = 0,
        int label_width=50,
        const std::string& title = ""
    ) {

      std::vector<const ThreadTimer*> tim_ptrs;
      for(const auto& tim : thread_timers_) tim_ptrs.push_back(&tim);
      const std::string indent = std::string(indent_size, ' ');
      out << std::setw(label_width) << std::left << (indent + title)
          << std::setw(8) << std::internal << "avg"
          << std::setw(8) << std::internal << "min"
          << std::setw(8) << std::internal << "max"
          << std::endl;
      print_sub(out, indent_size, tim_ptrs, name_, label_width);
      out << indent << "Timer overhead: " << std::setprecision(3)
          << (double)((unsigned long long)overhead_nanos_)/1.e9
          << std::endl;
      // TODO walltime/thread time ratio and efficiency
    }

};

//############################################################################//

} // end namespace sc

#endif /* _util_misc_thread_timer_h */
