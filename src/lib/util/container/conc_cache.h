//
// conc_cache.h
//
// Copyright (C) 2013 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Nov 19, 2013
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

#include <utility>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/pop_back.hpp>
#include <boost/mpl/back.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
//#include <boost/fusion/algorithm/transformation/zip.hpp>
//#include <boost/fusion/include/zip.hpp>
//#include <boost/fusion/algorithm/iteration/for_each.hpp>
//#include <boost/fusion/include/for_each.hpp>
//#include <boost/fusion/adapted/boost_tuple.hpp>
//#include <boost/fusion/include/boost_tuple.hpp>
//#include <boost/fusion/sequence/intrinsic/at.hpp>
//#include <boost/fusion/include/at.hpp>
#include <boost/utility.hpp>


#ifndef _util_container_conc_cache_h
#define _util_container_conc_cache_h


//#define LOCK_SEQ(ltype, k1, k2) \
//  /*mutex_type& m1 = std::get<0>(mutexes_)[k1];*/ \
//  /*mutex_type& m2 = std::get<1>(mutexes_)[k2];*/ \
//  mutex_type& m1 = get_mutex_1(k1); \
//  mutex_type& m2 = get_mutex_2(k2); \
//  ltype l1(m1, defer_lock); \
//  ltype l2(m2, defer_lock); \
//  lock(l1, l2);
//#define READ_LOCK(k1, k2) LOCK_SEQ(read_lock, k1, k2)
//#define WRITE_LOCK(k1, k2) LOCK_SEQ(write_lock, k1, k2)


namespace mpl = boost::mpl;

// These could be switched to std library versions someday
#define USE_BOOST_THREAD 1
#if USE_BOOST_THREAD
#include <boost/function.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/locks.hpp>
using boost::shared_mutex;
using boost::shared_lock;
using boost::unique_lock;
using boost::function;
using boost::lock;
using boost::defer_lock;
#else
#include <mutex>
#include <boost/thread/shared_mutex.hpp>
// This requires C++14 support
using std::shared_lock;
using std::unique_lock;
using std::shared_mutex;
using std::function;
using std::lock;
using std::defer_lock;
#endif

using std::make_pair;

// A copy-constructable Lockable type
//template<typename lock_type, typename mutex_type>
//struct lock_wrapper {
//  public:
//    lock_wrapper() : lock_(mutex_type(), defer_lock) { }
//    lock_wrapper(lock_wrapper& lock_w) : lock_(std::move(lock_w.lock_)) { }
//    lock_wrapper(mutex_type&& mtx) : lock_(std::forward<mutex_type>(mtx)) { }
//    lock_wrapper& operator=(lock_wrapper& other) { lock_.swap(other.lock_); return this; }
//    auto lock() -> decltype(this->lock_.lock()) { return lock_.lock(); }
//    auto unlock() -> decltype(this->lock_.unlock()) { return lock_.unlock(); }
//    auto try_lock() -> decltype(this->lock_.try_lock()) { return lock_.try_lock(); }
//    void set_mutex(mutex_type&& mtx) { lock_.swap(mtx); }
//  private:
//    lock_type lock_;
//};
/*
 *
    using boost::fusion::begin; \
    using boost::fusion::end; \
    using boost::fusion::next; \
    auto key_it = begin(key_tup);\
    auto mut_it = begin(mutexes_); \
    const int length = boost::tuples::length<key_tuple>::value; \
    for(int __ikey = 0; __ikey < length; ++__ikey) {\
      auto& keyi = *key_it; \
      auto& mtx_map = *mut_it; \
      auto found_spot = mtx_map.find(keyi); \
      if(found_spot == mtx_map.end()){ \
        mtx_map.emplace( \
            std::piecewise_construct, \
            std::forward_as_tuple(keyi), \
            std::forward_as_tuple() \
        );\
      }\
      __locks.emplace_back(\
          mtx_map[keyi],\
          defer_lock\
      );\
      if(__ikey != length-1) { \
        key_it = next(key_it); \
        mut_it = next(mut_it); \
      } \
    } \
 */

#define LOCK_IMPL(ltype, key_tup) \
  std::vector<ltype> __locks; \
  { \
    boost_ext::for_each(do_lock<decltype(__locks)>(&__locks), key_tup, mutexes_); \
    boost::lock(__locks.begin(), __locks.end()); \
  }

#define READ_LOCK(key_tup) LOCK_IMPL(read_lock, key_tup)
#define WRITE_LOCK(key_tup) LOCK_IMPL(write_lock, key_tup)


namespace sc {

// TODO arbitrary number/types of keys using the variable length template arguments in C++11
// TODO a number of other generalizations need to be made
/**
 * A cache of objects that can be safely accessed concurrently by threads that share memory.
 * For now, a simple two-key version is implemented
 */
//template <typename T, typename key_type_1, typename key_type_2>
//class ConcurrentCache2 {

//  public:
//    typedef T value_type;
//    typedef shared_mutex mutex_type;
//    typedef std::map<key_type_1, mutex_type> mutex_map_1;
//    typedef std::map<key_type_2, mutex_type> mutex_map_2;
//    template<typename... Types> using tuple_type = std::pair<Types...>;
//    typedef tuple_type<key_type_1, key_type_2> key_tuple;
//    typedef shared_lock<mutex_type> read_lock;
//    typedef unique_lock<mutex_type> write_lock;

//    // constness may change if LRU capabilities are implemented
//    bool contains(const key_type_1& k1, const key_type_2& k2) const {
//      READ_LOCK(k1, k2);
//      key_tuple k(k1, k2);
//      auto found_spot = cached_values_.find(k);
//      return found_spot != cached_values_.end();
//    }

//    // constness may change if LRU capabilities are implemented
//    boost::optional<value_type> get(const key_type_1& k1, const key_type_2& k2) const {
//      READ_LOCK(k1, k2);
//      key_tuple k(k1, k2);
//      auto found_spot = cached_values_.find(k);
//      if(found_spot == cached_values_.end()){
//        return boost::optional<value_type>();
//      }
//      else{
//        return boost::optional<value_type>(found_spot->second);
//      }
//    }

//    value_type get(
//        const key_type_1& k1, const key_type_2& k2,
//        const function<value_type()>& compute_fxn
//    )
//    {
//      key_tuple k(k1, k2);
//      {
//        READ_LOCK(k1, k2);
//        auto found_spot = cached_values_.find(k);
//        if(found_spot != cached_values_.end()){
//          return found_spot->second;
//        }
//      } // Release the read lock
//      //----------------------------------------//
//      // If we haven't returned at this point, we need to
//      //   compute the value
//      {
//        // All threads that get here will have
//        //   to wait for the one thread doing
//        //   the computing.  Threads that call
//        //   get() after the call_once() call
//        //   has completed will exit after the
//        //   read lock section.
//        WRITE_LOCK(k1, k2);
//        // Check again, in case another thread already
//        //   did the compute
//        auto found_spot = cached_values_.find(k);
//        if(found_spot != cached_values_.end()){
//          return found_spot->second;
//        }
//        else {
//          // Looks like we're the lucky ones who get
//          //   to do the compute.
//          cached_values_[k] = compute_fxn();
//        }
//      }
//      //----------------------------------------//
//      return cached_values_[k];
//    }

//    void set(
//        const key_type_1& k1, const key_type_2& k2,
//        const value_type& val
//    )
//    {
//      key_tuple k(k1, k2);
//      WRITE_LOCK(k1, k2);
//      cached_values_[k] = val;
//    }



//  private:

//    inline mutex_type& get_mutex_1(const key_type_1& k1){
//      auto found_spot = mutexes_.first.find(k1);
//      if(found_spot == mutexes_.first.end()){
//        // This is a mess, and probably could be
//        //   written less confusingly, but it works
//        mutexes_.first.emplace(
//            std::piecewise_construct,
//            std::forward_as_tuple(k1),
//            std::forward_as_tuple()
//        );
//      }
//      return mutexes_.first[k1];
//    }

//    inline mutex_type& get_mutex_2(const key_type_2& k2){
//      auto found_spot = mutexes_.second.find(k2);
//      if(found_spot == mutexes_.second.end()){
//        // This is a mess, and probably could be
//        //   written less confusingly, but it works
//        mutexes_.second.emplace(
//            std::piecewise_construct,
//            std::forward_as_tuple(k2),
//            std::forward_as_tuple()
//        );
//      }
//      return mutexes_.second[k2];
//    }

//    tuple_type<mutex_map_1, mutex_map_2> mutexes_;
//    std::map<key_tuple, value_type> cached_values_;


//};

namespace mpl_ext{
// mpl_ext::apply "splats" an mpl::vector (or any mpl back extensible sequence) into C++11 variadic
//   template arguments.
// credit: https://github.com/scientific-coder/Computer-Languages/blob/master/interpreting/apply.hxx
template<template<typename...> class T, bool empty, typename C, typename... Types> struct apply_helper {
  typedef typename boost::mpl::pop_back<C>::type rest;
  typedef typename apply_helper<T, boost::mpl::empty<rest>::value, rest, typename boost::mpl::back<C>::type, Types...>::type type;
};
template<template<typename...> class T,typename C, typename... Types> struct apply_helper<T, true, C, Types...> {
  typedef T<Types...> type;
};

template<template<typename...> class T,typename C> struct apply : apply_helper<T, boost::mpl::empty<C>::value, C> {};
} // end namespace mpl_ext

namespace boost_ext{

template<typename FuncT, int I = 0, typename... tuple_types>
inline typename boost::enable_if_c<
    I == boost::tuples::length<
        typename mpl::at<
            typename mpl::vector<tuple_types...>::type, mpl::int_<0>
        >::type
    >::value,
    void
>::type
for_each(FuncT, const tuple_types&...)
{ }

template<typename FuncT, int I = 0, typename... tuple_types>
inline typename boost::enable_if_c<
    I <= boost::tuples::length<
        typename mpl::at<
            typename mpl::vector<tuple_types...>::type, mpl::int_<0>
        >::type
    >::value - 1,
    void
>::type
for_each(FuncT f, tuple_types&... tups)
{
  f(boost::get<I>(tups)...);
  for_each<FuncT, I + 1, tuple_types...>(f, tups...);
}

} // end namespace boost_ext

template <typename val_type_param, typename... key_types>
class ConcurrentCache {

  public:
    typedef val_type_param value_type;
    typedef shared_mutex mutex_type;
    template<typename K, typename V> using map_type = std::map<K, V>;
    template<int keynum>
    struct mutex_map {
      typedef typename std::map<
        mpl::at<mpl::vector<key_types...>, mpl::int_<keynum>>,
        mutex_type
      > type;
    };
    template<typename... Types> using tuple_type = boost::tuple<Types...>;
    typedef tuple_type<key_types...> key_tuple;
    typedef shared_lock<mutex_type> read_lock;
    typedef unique_lock<mutex_type> write_lock;
    typedef typename mpl_ext::apply<
        tuple_type, typename mpl::transform<mpl::vector<key_types...>, std::map<mpl::_, mutex_type>>::type
    >::type mutex_map_tuple;

    typedef map_type<key_tuple, value_type> value_map;

    // constness may change if LRU capabilities are implemented
    bool contains(key_types... keys) const {
      key_tuple k(boost::forward<key_types>(keys)...);
      READ_LOCK(k);
      auto found_spot = cached_values_.find(k);
      return found_spot != cached_values_.end();
    }

    // constness may change if LRU capabilities are implemented
    boost::optional<value_type> get(key_types... keys) const {
      key_tuple k(boost::forward<key_types>(keys)...);
      READ_LOCK(k);
      auto found_spot = cached_values_.find(k);
      if(found_spot == cached_values_.end()){
        return boost::optional<value_type>();
      }
      else{
        return boost::optional<value_type>(found_spot->second);
      }
    }

    value_type get(
        key_types... keys,
        const function<value_type()>& compute_fxn
    )
    {
      key_tuple k(boost::forward<key_types>(keys)...);
      {
        READ_LOCK(k);
        auto found_spot = cached_values_.find(k);
        if(found_spot != cached_values_.end()){
          return found_spot->second;
        }
      } // Release the read lock
      //----------------------------------------//
      // If we haven't returned at this point, we need to
      //   compute the value
      {
        // All threads that get here will have
        //   to wait for the one thread doing
        //   the computing.  Threads that call
        //   get() after the call_once() call
        //   has completed will exit after the
        //   read lock section.
        WRITE_LOCK(k);
        // Check again, in case another thread already
        //   did the compute
        auto found_spot = cached_values_.find(k);
        if(found_spot != cached_values_.end()){
          return found_spot->second;
        }
        else {
          // Looks like we're the lucky ones who get
          //   to do the compute.
          cached_values_[k] = compute_fxn();
        }
      }
      //----------------------------------------//
      return cached_values_[k];
    }

    void set(
        key_types... keys,
        const value_type& val
    )
    {
      key_tuple k(boost::forward<key_types>(keys)...);
      WRITE_LOCK(k);
      cached_values_[k] = val;
    }

  private:

    mutex_map_tuple mutexes_;
    value_map cached_values_;

    template <typename lock_vector>
    struct do_lock {

      private:
        lock_vector* lockv_;

      public:

        do_lock(lock_vector* v) : lockv_(v) { }

        template<typename key_type, typename mutex_map_type>
        void operator()(key_type& key, mutex_map_type& mtx_map) const
        {
          auto found_spot = mtx_map.find(key);
          if(found_spot == mtx_map.end()){
            // This is a mess, and probably could be
            //   written less confusingly, but it works
            mtx_map.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(key),
                std::forward_as_tuple()
            );
          }

          (*lockv_).emplace_back(
              mtx_map[key],
              defer_lock
          );
        }
    };

};


} // end namespace sc
#endif /* _util_container_conc_cache_h */
