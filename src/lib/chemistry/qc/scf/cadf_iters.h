
// cadf_iters.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Jan 10, 2014
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

#ifndef _chemistry_qc_scf_cadf_iters_h
#define _chemistry_qc_scf_cadf_iters_h

#include <chemistry/qc/scf/clhf.h>
#include <util/misc/property.h>
#include <boost/thread/thread.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
//#include <boost/format.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/join.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <utility>
#include <type_traits>
#include <memory>
#include <iterator>
#include <atomic>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <limits>

#define DEFAULT_TARGET_BLOCK_SIZE 100  // functions
#define LINK_TIME_LIST_INSERTIONS 1
#define LINK_SORTED_INSERTION 0

// Relatively trivial boost extensions that promote readability without loss of performance
namespace boost_ext {

// From http://stackoverflow.com/questions/8511035/sequence-zip-function-for-c11
// Make sure the containers are the same length.  Behaviour is undefined otherwise
template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
    auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}

}

namespace sc {

#define DUMP(expr) std::cout << #expr << " = " << (expr) << std::endl;
#define DUMP2(expr1, expr2) std::cout << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << std::endl;
#define DUMP3(expr1, expr2, expr3) std::cout << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << ", " << #expr3 << " = " << (expr3) << std::endl;
#define DUMP4(expr1, expr2, expr3, expr4) std::cout << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << ", " << #expr3 << " = " << (expr3) << ", " << #expr4 << " = " << (expr4) << std::endl;
#define DUMP5(expr1, expr2, expr3, expr4, expr5) std::cout << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << ", " << #expr3 << " = " << (expr3) << ", " << #expr4 << " = " << (expr4) << ", " << #expr5 << " = " << (expr5) << std::endl;
#define out_assert(a, op, b) assert(a op b || ((std::cout << "Failed assertion output: " << #a << " ( = " << a << ") " << #op << " " << #b <<  " ( = " << b << ")" << std::endl), false))

enum {
  NotAssigned = -998,
  NoLastIndex = -999,
  NoMaximumBlockSize = -1001,
  IndicesEnd = -1002
};

// Forward declarations
class ShellData;
class BasisFunctionData;
class ShellIndexWithValue;
class OrderedShellList;
template<typename range_type> class ShellBlockIterator;
template<typename Iterator, typename ParentContainer> class IterableBasisElementData;

using int_range = decltype(boost::counting_range(1, 2));

struct ignored_argument { };

template <typename T>
class OptionalRefParameter {
    T* p_ = 0;

  public:

    OptionalRefParameter(T* p) : p_(p) { }

    OptionalRefParameter(const Ref<T>& rp) : p_(rp.nonnull() ? rp.pointer() : 0) { }

    operator T*() const { return p_; }
    operator const T*() const { return p_; }

};

typedef enum {
  NoRestrictions = 0,
  SameCenter = 1,
  SameAngularMomentum = 2
} BlockCompositionRequirement;

//############################################################################//

template<int size=1024, typename... Args>
std::string
stream_printf(const char* fmt, Args&&... args)
{
    char str[size];
    sprintf(str, fmt, std::forward<Args>(args)...);
    std::string strobj(str);
    return strobj;
}


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

    void enter(const std::string& subname) {
      if(active_subsection_) {
        active_subsection_->enter(subname);
      }
      else {
        active_subname_ = subname;
        auto subspot = subtimers_.find(subname);
        if(subspot != subtimers_.end()) {
          subspot->second.start();
          active_subsection_ = &(subspot->second);
        }
        else {
          auto insertion_pair = subtimers_.emplace(
              std::piecewise_construct,
              std::forward_as_tuple(subname),
              std::forward_as_tuple(depth_+1)
          );
          section_names_.push_back(subname);
          assert(insertion_pair.second);
          active_subsection_ = &(insertion_pair.first->second);
        }
      }
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
          << stream_printf("%7.2f %7.2f %7.2f",
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

template <typename Iterator, typename DataContainer> struct BasisElementIteratorDereference;

class BasisElementData {

  public:

    BasisElementData(
        int index,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis
    ) : index(index),
        basis(basis),
        dfbasis(dfbasis)
    { }

    BasisElementData() = delete;

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;

    int index = NotAssigned;

    bool operator!=(const BasisElementData& other) const
    {
      return index != other.index;
    }

    virtual ~BasisElementData() { }

    int center = NotAssigned;
    int bfoff_in_atom = NotAssigned;

};

struct ShellData : public BasisElementData {

    template<typename Iterator> using with_iterator =
        BasisElementIteratorDereference<ShellData, Iterator>;

    ShellData(
        int index,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0,
        int block_offset = NotAssigned
    ) : BasisElementData(index, basis, dfbasis),
        block_offset(NotAssigned)
    {
      init();
    }

    template <typename Iterator>
    ShellData(
        const with_iterator<Iterator>& deref
    ) : ShellData(deref.index, deref.basis, deref.dfbasis)
    { }

    ShellData() : BasisElementData(NotAssigned, 0, 0) { }


    int bfoff = NotAssigned;
    int nbf = NotAssigned;
    int atom_bfoff = NotAssigned;
    int atom_shoff = NotAssigned;
    int atom_nsh = NotAssigned;
    int atom_nbf = NotAssigned;
    int shoff_in_atom = NotAssigned;
    int atom_last_function = NotAssigned;
    int atom_last_shell = NotAssigned;
    int last_function = NotAssigned;

    // Used when an auxiliary basis is set in the parent ShellIter.  Otherwise, set to -1
    int atom_dfshoff = NotAssigned;
    int atom_dfbfoff = NotAssigned;
    int atom_dfnbf = NotAssigned;
    int atom_dfnsh = NotAssigned;
    int atom_df_last_function = NotAssigned;
    int atom_df_last_shell = NotAssigned;


    // Reserved for future use
    int block_offset = NotAssigned;

    operator const int() const {
      out_assert(basis, !=, 0);
      out_assert(index, <, basis->nshell());
      out_assert(index, >=, 0);
      return index;
    }

    static int max_index(GaussianBasisSet* basis) {
      return basis->nshell();
    }

  protected:

    void init(){

      assert(basis);

      if(index == NotAssigned || index == basis->nshell()) return;
      if(index >= basis->nshell() || index < 0) return;

      //out_assert(index, <, basis->nshell());
      //out_assert(index, >=, 0);

      nbf = basis->shell(index).nfunction(); 
      bfoff = basis->shell_to_function(index); 
      center = basis->shell_to_center(index); 
      atom_shoff = basis->shell_on_center(center, 0);
      atom_bfoff = basis->shell_to_function(atom_shoff); 
      atom_nsh = basis->nshell_on_center(center); 
      atom_nbf = basis->nbasis_on_center(center); 
      shoff_in_atom = index - atom_shoff; 
      bfoff_in_atom = bfoff - atom_bfoff; 
      atom_last_function = atom_bfoff + atom_nbf - 1; 
      last_function = bfoff + nbf - 1; 
      atom_last_shell = atom_shoff + atom_nsh - 1; 
      if(dfbasis != 0) {
        atom_dfshoff = dfbasis->shell_on_center(center, 0); 
        atom_dfbfoff = dfbasis->shell_to_function(atom_dfshoff); 
        atom_dfnsh = dfbasis->nshell_on_center(center); 
        atom_dfnbf = dfbasis->nbasis_on_center(center); 
        atom_df_last_function = atom_dfbfoff + atom_dfnbf - 1; 
        atom_df_last_shell = atom_dfshoff + atom_dfnsh - 1; 
      }

    }
};

struct BasisFunctionData : public BasisElementData {

    template<typename Iterator> using with_iterator =
        BasisElementIteratorDereference<BasisFunctionData, Iterator>;

    BasisFunctionData(
        int idx,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        int block_offset = NotAssigned
    ) : BasisElementData(idx, basis, dfbasis),
        block_offset(block_offset)
    {
      init();
    }

    template <typename Iterator>
    BasisFunctionData(
        const with_iterator<Iterator>& deref
    ) : BasisFunctionData(deref.index, deref.basis, deref.dfbasis, deref.block_offset)
    { }

    BasisFunctionData()
      : BasisFunctionData(NotAssigned, 0, 0)
    { }

    int shell_index = NotAssigned;
    int shell_bfoff = NotAssigned;
    int bfoff_in_shell = NotAssigned;
    int atom_dfshoff = NotAssigned;
    int atom_dfbfoff = NotAssigned;
    int atom_dfnbf = NotAssigned;
    int atom_shoff = NotAssigned;
    int atom_bfoff = NotAssigned;
    int block_offset = NotAssigned;
    int bfoff_in_block = NotAssigned;

    operator const int() const { return index; }

    static int max_index(GaussianBasisSet* basis) {
      return basis->nbasis();
    }

  protected:

    void init()
    {
      if(index == NotAssigned || index == basis->nbasis()) return;

      out_assert(basis, !=, 0);
      out_assert(index, <, basis->nbasis());
      out_assert(index, >=, 0);

      shell_index = basis->function_to_shell(index);
      shell_bfoff = basis->shell_to_function(shell_index);
      center = basis->shell_to_center(shell_index);
      bfoff_in_shell = index - shell_bfoff;
      atom_shoff = basis->shell_on_center(center, 0);
      atom_bfoff =  basis->shell_to_function(atom_shoff);
      bfoff_in_atom = index - atom_bfoff;
      if(dfbasis != 0){
        atom_dfshoff = dfbasis->shell_on_center(center, 0);
        atom_dfbfoff = dfbasis->shell_to_function(atom_dfshoff);
        atom_dfnbf = dfbasis->nbasis_on_center(center);
      }
      if(block_offset != NotAssigned) {
        bfoff_in_block = index - block_offset;
      }
    }
};

template <typename DataContainer, typename Iterator>
struct BasisElementIteratorDereference : DataContainer {

    const Iterator iterator;

    template<typename... Args>
    BasisElementIteratorDereference(
        Iterator iter,
        Args&&... args
    ) : DataContainer(std::forward<Args>(args)...), iterator(iter)
    { }

};

//############################################################################//

template <typename DataContainer, typename Iterator=int_range::iterator>
class basis_element_iterator
  : public boost::iterator_adaptor<
      basis_element_iterator<DataContainer, Iterator>,
      Iterator,
      boost::use_default,
      boost::bidirectional_traversal_tag,
      BasisElementIteratorDereference<DataContainer, Iterator>
    >
{

  public:

    typedef basis_element_iterator<DataContainer, Iterator> self_t;
    typedef boost::iterator_adaptor<
         self_t, Iterator, boost::use_default, boost::bidirectional_traversal_tag,
         BasisElementIteratorDereference<DataContainer, Iterator>
     > super_t;
    typedef typename DataContainer::template with_iterator<Iterator> value_reference_t;


  protected:

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;
    int block_offset = NotAssigned;

    friend class boost::iterator_core_access;

    value_reference_t dereference() const {
      auto base_spot = super_t::base();
      return value_reference_t(
        base_spot,
        int(*base_spot),
        basis, dfbasis,
        block_offset
      );
    }

  public:

    basis_element_iterator() : basis(0), dfbasis(0) { }

    basis_element_iterator(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        Iterator iter,
        int block_offset = NotAssigned
    ) : basis_element_iterator::iterator_adaptor_(iter),
        basis(basis),
        dfbasis(dfbasis),
        block_offset(block_offset)
    { }


};

template <typename DataContainer, typename Iterator=int_range::iterator>
class basis_element_with_value_iterator
  : public boost::iterator_adaptor<
      basis_element_with_value_iterator<DataContainer, Iterator>,
      Iterator,
      boost::use_default,
      boost::bidirectional_traversal_tag,
      BasisElementIteratorDereference<DataContainer, Iterator>
    >
{

  public:

    typedef basis_element_with_value_iterator<DataContainer, Iterator> self_t;
    typedef boost::iterator_adaptor<
         self_t, Iterator, boost::use_default, boost::bidirectional_traversal_tag,
         BasisElementIteratorDereference<DataContainer, Iterator>
     > super_t;
    typedef BasisElementIteratorDereference<DataContainer, Iterator> value_reference_t;

    basis_element_with_value_iterator() : basis(0), dfbasis(0) { }

    basis_element_with_value_iterator(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        Iterator iter,
        int block_offset = NotAssigned
    ) : super_t(iter),
        basis(basis),
        dfbasis(dfbasis),
        block_offset(block_offset)
    { }

    operator basis_element_iterator<DataContainer, Iterator>()
    {
      typedef basis_element_iterator<DataContainer, Iterator> result_type;
      return result_type(basis, dfbasis, super_t::base(), block_offset);
    }

  private:

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;
    int block_offset = NotAssigned;

    friend class boost::iterator_core_access;

    value_reference_t dereference() const {
      auto base_spot = super_t::base();
      return value_reference_t(
        base_spot,
        (*base_spot).index,
        (*base_spot).value,
        basis, dfbasis,
        block_offset
      );
    }
};

template<typename Iterator=int_range::iterator>
    using shell_iterator = basis_element_iterator<ShellData, Iterator>;

template<typename Iterator=int_range::iterator>
    using function_iterator = basis_element_iterator<BasisFunctionData, Iterator>;


//############################################################################//

template<typename DataContainer, typename Iterator=int_range::iterator>
    using range_of = decltype(boost::make_iterator_range(
        basis_element_iterator<DataContainer, Iterator>(),
        basis_element_iterator<DataContainer, Iterator>()
    ));

template<typename DataContainer, typename Iterator=int_range::iterator>
    using valued_range_of = decltype(boost::make_iterator_range(
        basis_element_with_value_iterator<DataContainer, Iterator>(),
        basis_element_with_value_iterator<DataContainer, Iterator>()
    ));

template <typename DataContainer>
inline range_of<DataContainer>
basis_element_range(
    GaussianBasisSet* basis,
    OptionalRefParameter<GaussianBasisSet> dfbasis,
    int first_index,
    int last_index,
    int block_offset = NotAssigned
)
{
  const int actual_last_index = last_index == NoLastIndex ?
      DataContainer::max_index(basis) - 1 : last_index;
  auto base_range = boost::counting_range(first_index, actual_last_index + 1);
  return boost::make_iterator_range(
      basis_element_iterator<DataContainer>(basis, dfbasis, base_range.begin(), block_offset),
      basis_element_iterator<DataContainer>(basis, dfbasis, base_range.end(), block_offset)
  );

}

template <typename DataContainer>
inline range_of<DataContainer>
basis_element_range(
    GaussianBasisSet* basis,
    int last_index
)
{
  return basis_element_range<DataContainer>(basis, 0, 0, last_index);
}

template <typename DataContainer>
inline range_of<DataContainer>
basis_element_range(
    GaussianBasisSet* basis,
    GaussianBasisSet* dfbasis = 0,
    int last_index = NoLastIndex
)
{
  return basis_element_range<DataContainer>(basis, dfbasis, 0, last_index);
}

template <typename DataContainer>
inline range_of<DataContainer>
basis_element_range(
    GaussianBasisSet* basis,
    int first_index,
    int last_index
)
{
  return basis_element_range<DataContainer>(basis, 0, first_index, last_index);
}

template <typename Iterator>
inline range_of<ShellData, Iterator>
shell_range(
    const ShellData::with_iterator<Iterator>& begin,
    const ShellData::with_iterator<Iterator>& end
)
{
  return boost::make_iterator_range(
      basis_element_iterator<ShellData, Iterator>(begin.basis, begin.dfbasis, begin.iterator),
      basis_element_iterator<ShellData, Iterator>(begin.basis, begin.dfbasis, end.iterator)
  );
}

template <typename Iterator, typename Iterator2>
inline range_of<ShellData, Iterator>
shell_range(
    const ShellData::with_iterator<Iterator>& begin,
    const Iterator2& end
)
{
  return boost::make_iterator_range(
      basis_element_iterator<ShellData, Iterator>(begin.basis, begin.dfbasis, begin.iterator),
      basis_element_iterator<ShellData, Iterator>(begin.basis, begin.dfbasis, Iterator(end))
  );
}

template <typename Iterator>
inline range_of<ShellData, Iterator>
shell_range(
    const Iterator& begin,
    const Iterator& end,
    GaussianBasisSet* basis,
    GaussianBasisSet* dfbasis
)
{
  return boost::make_iterator_range(
      basis_element_iterator<ShellData, Iterator>(basis, dfbasis, begin),
      basis_element_iterator<ShellData, Iterator>(basis, dfbasis, end)
  );
}

template <typename Iterator>
inline range_of<BasisFunctionData, Iterator>
function_range(
    const BasisFunctionData::with_iterator<Iterator>& begin,
    const BasisFunctionData::with_iterator<Iterator>& end
)
{
  return boost::make_iterator_range(
      basis_element_iterator<BasisFunctionData, Iterator>(begin.basis, begin.dfbasis, begin.iterator),
      basis_element_iterator<BasisFunctionData, Iterator>(begin.basis, begin.dfbasis, end.iterator)
  );
}

template <typename... Args>
inline range_of<ShellData>
shell_range(GaussianBasisSet* basis, Args&&... args) {
  return basis_element_range<ShellData>(basis, std::forward<Args>(args)...);
}

template <typename... Args>
inline range_of<BasisFunctionData>
function_range(GaussianBasisSet* basis, Args&&... args) {
  return basis_element_range<BasisFunctionData>(basis, std::forward<Args>(args)...);
}

inline range_of<BasisFunctionData>
function_range(const ShellData& ish) {
  return basis_element_range<BasisFunctionData>(
      ish.basis, ish.dfbasis,
      ish.bfoff, ish.last_function
  );
}

//############################################################################//

template<typename Iterator> class ShellBlockSkeleton;

template<typename Range = range_of<ShellData>>
class ShellBlockData {

    void init();

  public:

    Range shell_range;
    ShellData first_shell;
    ShellData last_shell;
    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;
    int restrictions;

    ShellBlockData() : restrictions(NotAssigned) { }

    ShellBlockData(const ShellBlockSkeleton<Range>&);
    
    explicit ShellBlockData(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0
    ) : ShellBlockData(
          sc::shell_range(basis, dfbasis, 0, basis->nshell() - 1),
          basis->nshell(), basis->nbasis(),
          NoRestrictions
        )
    { }

    // Construct a one-shell shell block
    ShellBlockData(const ShellData& ish)
      : ShellBlockData(
          sc::shell_range(ish.basis, ish.dfbasis, ish, ish),
          1, ish.nbf, SameCenter|SameAngularMomentum
        )
    { }

    ShellBlockData(
        Range sh_range,
        int nshell, int nbf,
        int requirements
    ) : shell_range(sh_range),
        first_shell(*shell_range.begin()),
        last_shell(shell_range.begin() == shell_range.end() ?
            *shell_range.end() : *(--shell_range.end())
        ),
        basis(first_shell.basis),
        dfbasis(first_shell.dfbasis),
        nshell(nshell), nbf(nbf),
        restrictions(requirements)
    { init(); }

    int nbf;
    int bfoff;
    int nshell;
    int last_function;
    int center = NotAssigned;
    int atom_bfoff = NotAssigned;
    int atom_shoff = NotAssigned;
    int atom_nsh = NotAssigned;
    int atom_nbf = NotAssigned;
    int bfoff_in_atom = NotAssigned;
    int shoff_in_atom = NotAssigned;
    int atom_last_function = NotAssigned;
    int atom_last_shell = NotAssigned;
    int atom_dfshoff = NotAssigned;
    int atom_dfbfoff = NotAssigned;
    int atom_dfnbf = NotAssigned;
    int atom_dfnsh = NotAssigned;
    int atom_df_last_function = NotAssigned;
    int atom_df_last_shell = NotAssigned;

    static int max_index(GaussianBasisSet* basis) { return basis->nshell(); }

};

template<typename Range = range_of<ShellData>>
class ShellBlockSkeleton {

  private:

    int nshell;
    int nbf;


    friend class ShellBlockData<Range>;

  public:

    Range shell_range;
    int first_index;
    int restrictions;

    ShellBlockSkeleton() : nshell(0), nbf(0), first_index(NotAssigned) { }

    static ShellBlockSkeleton<Range> end_skeleton() { return ShellBlockSkeleton(); }

    ShellBlockSkeleton(const ShellBlockData<Range>& shbd)
      : shell_range(shbd.shell_range), nshell(shbd.nshell), nbf(shbd.nbf),
        restrictions(shbd.restrictions),
        first_index((*shell_range.begin()).index)
    { }

    ShellBlockSkeleton(Range range, int nsh, int nbas, int restrictions)
      : shell_range(range), nshell(nsh), nbf(nbas),
        restrictions(restrictions),
        first_index((*shell_range.begin()).index)
    { }


    template<typename OtherRange>
    bool operator==(const ShellBlockSkeleton<OtherRange>& other) const
    {
      return !(this->operator!=(other));
    }

    template<typename OtherRange>
    bool operator!=(const ShellBlockSkeleton<OtherRange>& other) const
    {
      return first_index != other.first_index
          or nshell != other.nshell
          or nbf != other.nbf;
    }

    template<typename OtherRange>
    bool operator<(const ShellBlockSkeleton<OtherRange>& other) const
    {
      const int my_index = (*shell_range.begin()).index;
      const int other_index = (*other.shell_range.begin()).index;
      if(my_index < other_index) return true;
      else if(my_index > other_index) return false;
      else if(nshell < other.nshell) return true;
      else if(nshell > other.nshell) return false;
      else if(nbf < other.nbf) return true;
      else return false;
    }

};

template<typename ShellIterator, typename ShellRange=range_of<ShellData, ShellIterator>>
class shell_block_iterator
  : public boost::iterator_facade<
      shell_block_iterator<ShellIterator, ShellRange>,
      ShellBlockSkeleton<range_of<ShellData, ShellIterator>>,
      boost::forward_traversal_tag,
      ShellBlockData<range_of<ShellData, ShellIterator>>
    >
{
  public:

    typedef shell_block_iterator<ShellIterator> self_type;
    typedef ShellBlockData<ShellRange> value_reference;

  private:

    int target_size;
    int restrictions;

    ShellRange all_shells;
    ShellBlockSkeleton<ShellRange> current_skeleton;

    void init();
    void init_from_spot(const decltype(all_shells.begin())& start_spot);


    friend class boost::iterator_core_access;

    value_reference dereference() const {
      return value_reference(
          current_skeleton
      );
    }

    template <typename OtherIterator>
    bool equal(shell_block_iterator<OtherIterator> const& other) const
    {
      return current_skeleton == other.current_skeleton;
    }

    void increment()
    {
      // Don't allow incrementing of end iterator
      assert(current_skeleton.first_index != NotAssigned);
      all_shells = shell_range(*current_skeleton.shell_range.end(), *all_shells.end());
      init();
    }

  public:

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;

    shell_block_iterator() : basis(0), dfbasis(0), restrictions(0), target_size(0) { }

    shell_block_iterator(
        const ShellRange& all_shells_in,
        int requirements = SameCenter,
        int target_size = DEFAULT_TARGET_BLOCK_SIZE
    ) : basis(all_shells_in.begin()->basis),
        dfbasis(all_shells_in.begin()->dfbasis),
        restrictions(requirements),
        target_size(target_size),
        all_shells(all_shells_in)
    {
      init();
    }

    shell_block_iterator(
        const ShellIterator& index_begin,
        const ShellIterator& index_end,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0,
        int requirements = SameCenter,
        int target_size = DEFAULT_TARGET_BLOCK_SIZE
    ) : basis(basis),
        dfbasis(dfbasis),
        restrictions(requirements),
        target_size(target_size),
        all_shells(shell_range(index_begin, index_end, basis, dfbasis))
    {
      init();
    }

    shell_block_iterator(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0,
        int first_index = 0,
        int last_index = NoLastIndex,
        int requirements = SameCenter,
        int target_size = DEFAULT_TARGET_BLOCK_SIZE
    ) : basis(basis),
        dfbasis(dfbasis),
        restrictions(requirements),
        target_size(target_size),
        all_shells(shell_range(
            basis, dfbasis, first_index,
            last_index == NoLastIndex ? ShellData::max_index(basis) - 1 : last_index
        ))
    {
      init();
    }

    static shell_block_iterator end_of_basis(GaussianBasisSet* basis) {
      return shell_block_iterator(basis, 0, basis->nshell(), basis->nshell()-1);
    }

    static shell_block_iterator end_with_last_index(int index, GaussianBasisSet* basis) {
      // The end iterator needs to start after the "last_index" that is passed in.  But
      //   because last_index gets incremented before being made into a counting range,
      //   the last_index parameter of the end iterator should be one less than the first.
      //   This is counterintuitive, but not meant to be exposed to the outside world anyway.
      return shell_block_iterator(basis, 0, index + 1, index);
    }

    const self_type& requiring(int restr) {
      restrictions = restr;
      init();
      return *this;
    }

};

// Type alias specialization for shell blocks
template<typename Iterator=int_range::iterator, typename Range=range_of<ShellData, Iterator>>
    using range_of_shell_blocks = decltype(boost::make_iterator_range(
        shell_block_iterator<Iterator, Range>(),
        shell_block_iterator<Iterator, Range>()
    ));

inline range_of_shell_blocks<>
shell_block_range(
    GaussianBasisSet* basis,
    OptionalRefParameter<GaussianBasisSet> dfbasis = 0,
    int first_index = 0,
    int last_index = NoLastIndex,
    int reqs = SameCenter,
    int target_size = DEFAULT_TARGET_BLOCK_SIZE
)
{
  typedef shell_block_iterator<int_range::iterator> result_t;
  const int actual_last_index = last_index == NoLastIndex ?
      basis->nshell() - 1 : last_index;
  return boost::make_iterator_range(
      result_t(basis, dfbasis, first_index,
          actual_last_index, reqs, target_size
      ),
      result_t::end_with_last_index(actual_last_index, basis)
  );

}

inline std::ostream&
operator << (std::ostream& out, const ShellData& ish)
{
  out << "ShellData: {"
      << "index = " << ish.index << ", ";
  out << "center = ";
  if(ish.center == NotAssigned) out << "NotAssigned, ";
  else out << ish.center << ", ";
  out << "bfoff = ";
  if(ish.bfoff == NotAssigned) out << "NotAssigned, ";
  else out << ish.bfoff << ", ";
  out << "nbf = ";
  if(ish.nbf == NotAssigned) out << "NotAssigned, ";
  else out << ish.nbf << ", ";
  out << "last_function = ";
  if(ish.nbf == NotAssigned) out << "NotAssigned";
  else out << ish.last_function;
  out << "}";
  return out;
}

template<typename Range>
std::ostream&
operator << (std::ostream& out, const ShellBlockData<Range>& blk)
{
  auto write_if_assigned = [&](int val) {
    if(val == NotAssigned) out << "NotAssigned";
    else out << val;
  };
  out << "\nShellBlockData: {"
      << "\n  nshell: " << blk.nshell
      << ", nbf: " << blk.nbf;
  out << "\n  shoff: ";
  write_if_assigned(blk.first_shell.index);
  out << ", bfoff: ";
  write_if_assigned(blk.bfoff);
  out << "\n  (basis nshell = "
      << blk.basis->nshell()
      << ", nbf = "
      << blk.basis->nbasis()
      << ")"
      << "\n  center: ";
  write_if_assigned(blk.center);
  out << "\n  shells:";
  for(auto sh : shell_range(blk)) {
    out << "\n    " << sh;
  }
  out << "\n}" << std::endl;
  return out;
}

//############################################################################//

inline range_of<BasisFunctionData>
iter_functions_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0
)
{
  const int shoff = basis->shell_on_center(center, 0);
  const int bfoff = basis->shell_to_function(shoff);
  return function_range(
      basis, dfbasis,
      bfoff, bfoff + basis->nbasis_on_center(center) - 1
  );
}

inline range_of<ShellData>
iter_shells_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0
)
{
  const int shoff = basis->shell_on_center(center, 0);
  return shell_range(
      basis, dfbasis,
      shoff, shoff + basis->nshell_on_center(center) - 1
  );
}

template<typename DataContainer, typename Iterator=int_range::iterator>
using joined_range_of = decltype(
    boost::join(
        range_of<DataContainer, Iterator>(),
        range_of<DataContainer, Iterator>()
    )
);

inline joined_range_of<ShellData>
iter_shells_on_centers(
    const Ref<GaussianBasisSet>& basis,
    int center1,
    int center2,
    const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0
)
{
  const int shoff1 = basis->shell_on_center(center1, 0);
  const int shoff2 = basis->shell_on_center(center2, 0);
  if(center1 == center2) {
    return boost::join(
      shell_range(
        basis, dfbasis,
        shoff1, shoff1 + basis->nshell_on_center(center1) - 1
      ),
      // An empty range
      shell_range(
        basis, dfbasis,
        shoff1, shoff1 - 1
      )
    );
  }
  else {
    return boost::join(
      shell_range(
        basis, dfbasis,
        shoff1, shoff1 + basis->nshell_on_center(center1) - 1
      ),
      shell_range(
        basis, dfbasis,
        shoff2, shoff2 + basis->nshell_on_center(center2) - 1
      )
    );
  }
}

inline range_of_shell_blocks<>
iter_shell_blocks_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0,
    int reqs = SameCenter,
    int target_size = DEFAULT_TARGET_BLOCK_SIZE
)
{
  const int shoff = basis->shell_on_center(center, 0);
  return shell_block_range(
      basis, dfbasis,
      shoff,
      shoff + basis->nshell_on_center(center) - 1,
      reqs|SameCenter,
      target_size
  );
}

//############################################################################//

struct ShellDataWithValue : public ShellData {

    template<typename Iterator> using with_iterator =
        BasisElementIteratorDereference<ShellDataWithValue, Iterator>;

    double value;

    ShellDataWithValue(
        int index,
        double value,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        int block_offset = NotAssigned
    ) : ShellData(index, basis, dfbasis), value(value)
    { }

    ShellDataWithValue(
        const ShellData& ish,
        double value
    ) : ShellData(ish), value(value)
    { }

};

template <typename Iterator>
inline valued_range_of<ShellDataWithValue, Iterator>
shell_range(
    const ShellDataWithValue::with_iterator<Iterator>& begin,
    const ShellDataWithValue::with_iterator<Iterator>& end
)
{
  return boost::make_iterator_range(
      basis_element_with_value_iterator<ShellDataWithValue, Iterator>(begin.basis, begin.dfbasis, begin.iterator),
      basis_element_with_value_iterator<ShellDataWithValue, Iterator>(begin.basis, begin.dfbasis, end.iterator)
  );
}

template <typename Iterator, typename Iterator2>
inline valued_range_of<ShellData, Iterator>
shell_range(
    const ShellDataWithValue::with_iterator<Iterator>& begin,
    const Iterator2& end
)
{
  return boost::make_iterator_range(
      basis_element_with_value_iterator<ShellDataWithValue, Iterator>(begin.basis, begin.dfbasis, begin.iterator),
      basis_element_with_value_iterator<ShellDataWithValue, Iterator>(begin.basis, begin.dfbasis, Iterator(end))
  );
}

template <typename Iterator>
inline valued_range_of<ShellData, Iterator>
shell_range(
    const basis_element_with_value_iterator<ShellData, Iterator>& begin,
    const basis_element_with_value_iterator<ShellData, Iterator>& end
)
{
  return boost::make_iterator_range(begin, end);
}

struct ShellIndexWithValue{

    int index;
    mutable double value;

    ShellIndexWithValue() : index(NotAssigned) { }

    ShellIndexWithValue(int index)
      : index(index), value(0.0)
    { }

    ShellIndexWithValue(int index, double value)
      : index(index), value(value)
    { }

    bool operator==(const ShellIndexWithValue& other) {
      return index == other.index;
    }

    operator int() const {
      return index;
    }
};

namespace {
  template<typename T>
  struct hash_;

  template<>
  struct hash_<ShellIndexWithValue> {

      std::size_t operator()(const ShellIndexWithValue& v) const {
        return std::hash<int>()(v.index);
      }

  };

  struct index_equal_ {
      bool operator()(const ShellIndexWithValue& a, const ShellIndexWithValue& b) {
        return a.index == b.index;
      }
  };
}

class OrderedShellList {

    struct index_compare {
        bool operator()(const ShellIndexWithValue& a, const ShellIndexWithValue& b) const {
          return a.index < b.index;
        }
    };

    struct value_compare {
        bool operator()(const ShellIndexWithValue& a, const ShellIndexWithValue& b) const {
          return a.value > b.value or (a.value == b.value && a.index > b.index);
        }
    };

  public:
    #if LINK_SORTED_INSERTION
    typedef std::set<ShellIndexWithValue, value_compare> index_list;
    #else
    typedef std::vector<ShellIndexWithValue> index_list;
    typedef std::unordered_set<
        ShellIndexWithValue,
        hash_<ShellIndexWithValue>,
        index_equal_
    > index_set;
    #endif
    typedef index_list::const_iterator index_iterator;
    typedef basis_element_with_value_iterator<ShellDataWithValue, index_iterator> iterator;

  private:

    boost::mutex get_next_mtx_;
    boost::mutex insert_mtx_;

    #if LINK_SORTED_INSERTION
    index_list indices_{ value_compare() };
    std::set<ShellIndexWithValue, index_compare> index_ordered_{ index_compare() };
    #else
    index_list indices_;
    index_set idx_set_;
    bool sorted_ = false;
    #endif

  public:

    GaussianBasisSet* basis_ = 0;
    GaussianBasisSet* dfbasis_ = 0;

    #if LINK_TIME_LIST_INSERTIONS
    static std::atomic_uint_fast64_t insertion_time_nanos;
    static std::atomic_uint_fast64_t n_insertions;
    time_accumulator_factory<> insertion_timer{ insertion_time_nanos };
    #endif

    OrderedShellList()
    { }

    OrderedShellList(const OrderedShellList& other)
      : indices_(other.indices_)
        #if LINK_SORTED_INSERTION
        , index_ordered_(other.index_ordered_)
        #else
        , idx_set_(other.idx_set_)
        , sorted_(other.sorted_)
        #endif

    { }
    
    void set_basis(GaussianBasisSet* basis, GaussianBasisSet* dfbasis = 0)
    {
      basis_ = basis;
      if(dfbasis) dfbasis_ = dfbasis;
    }

    void insert(const ShellData& ish, double value) {
      #if LINK_TIME_LIST_INSERTIONS
      ++n_insertions;
      auto timer = insertion_timer.create();
      #endif
      //----------------------------------------//
      assert(basis_ == 0 || ish.basis == basis_);
      assert(ish.basis);
      if(basis_ == 0) basis_ = ish.basis;
      assert(dfbasis_ == 0 || ish.dfbasis == dfbasis_);
      if(dfbasis_ == 0) dfbasis_ = ish.dfbasis;
      //----------------------------------------//
      boost::lock_guard<boost::mutex> lg(insert_mtx_);
      #if LINK_SORTED_INSERTION
      ShellIndexWithValue insert_val((int)ish, value);
      const auto& found = index_ordered_.find(insert_val);
      if(found != index_ordered_.end()) {
        const auto& found_old = indices_.find(*found);
        assert(found_old != indices_.end());
        const double old_value = found_old->value;
        if(value > old_value) {
          found->value = value;
          indices_.erase(found_old);
          indices_.insert(insert_val);
        }
      }
      else{
        index_ordered_.insert(insert_val);
        indices_.insert(insert_val);
      }
      #else
      sorted_ = false;
      ShellIndexWithValue insert_val(ish, value);
      const auto& found = idx_set_.find(insert_val);
      if(found != idx_set_.end()) {
        const double old_val = found->value;
        if(value > old_val) found->value = value;
      }
      else {
        idx_set_.insert(insert_val);
      }
      #endif
    }

    #if !LINK_SORTED_INSERTION
    void sort() {
      std::move(idx_set_.begin(), idx_set_.end(), std::back_inserter(indices_));
      std::sort(indices_.begin(), indices_.end(),
          [](const ShellIndexWithValue& a, const ShellIndexWithValue& b){
            return a.value > b.value;
          }
      );
      sorted_ = true;
    }
    #endif

    iterator begin() const
    {
      #if !LINK_SORTED_INSERTION
      assert(sorted_);
      #endif
      return iterator(
          basis_, dfbasis_,
          indices_.cbegin()
      );
    }

    index_iterator index_begin() const
    {
      #if !LINK_SORTED_INSERTION
      assert(sorted_);
      #endif
      return indices_.cbegin();
    }

    iterator end() const
    {
      #if !LINK_SORTED_INSERTION
      assert(sorted_);
      #endif
      return iterator(
          basis_, dfbasis_,
          indices_.cend()
      );
    }

    index_iterator index_end() const
    {
      #if !LINK_SORTED_INSERTION
      assert(sorted_);
      #endif
      return indices_.cend();
    }

    auto size() -> decltype(indices_.size())
    {
      #if !LINK_SORTED_INSERTION
      assert(sorted_);
      #endif
      return indices_.size();
    }

    void clear()
    {
      indices_.clear();
      #if LINK_SORTED_INSERTION
      index_ordered_.clear();
      #else
      idx_set_.clear();
      #endif
    }

};

inline range_of_shell_blocks<OrderedShellList::index_iterator>
shell_block_range(
    const OrderedShellList& shlist,
    int requirements=SameCenter,
    int target_size=DEFAULT_TARGET_BLOCK_SIZE
)
{
  return boost::make_iterator_range(
      shell_block_iterator<OrderedShellList::index_iterator>(
          shlist.index_begin(), shlist.index_end(),
          shlist.basis_, shlist.dfbasis_, requirements, target_size
      ),
      shell_block_iterator<OrderedShellList::index_iterator>(
          shlist.index_end(), shlist.index_end(),
          shlist.basis_, shlist.dfbasis_, requirements, target_size
      )
  );
}


inline
std::ostream&
operator <<(std::ostream& out, const OrderedShellList& list) {
  out << "\nOrderedShellList: {";
  for(auto ish : list) {
    out << "\n  " << ish << ", with value " << ish.value;
  }
  out << "\n}" << std::endl;
  return out;
}

//############################################################################//

} // end namespace sc

#include <chemistry/qc/scf/cadf_iters_impl.h>

#endif /* _chemistry_qc_scf_cadf_iters_h */
