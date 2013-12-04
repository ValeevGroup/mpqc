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

#include <cassert>
#include <utility>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/pop_back.hpp>
#include <boost/mpl/back.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/size.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/utility.hpp>
#include <boost/type_traits.hpp>


#ifndef _util_container_conc_cache_h
#define _util_container_conc_cache_h

namespace mpl = boost::mpl;

#define ALLOW_DYNAMIC_MUTEX_CREATION 0

#define USE_BOOST_THREAD 1
#if USE_BOOST_THREAD
#include <boost/function.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/locks.hpp>
using boost::mutex;
using boost::shared_mutex;
using boost::upgrade_mutex;
using boost::shared_lock;
using boost::unique_lock;
using boost::function;
using boost::lock;
using boost::defer_lock;
using boost::upgrade_lock;
using boost::upgrade_to_unique_lock;
#else
#include <mutex>
#include <shared_mutex>
using std::mutex;
// This requires C++14 support
using std::shared_lock;
using std::unique_lock;
using std::shared_mutex;
using std::function;
using std::lock;
using std::defer_lock;
using std::upgrade_lock;
using std::upgrade_to_unique_lock;
#endif

namespace sc {

////////////////////////////////////////////////////////////////////////////////

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

template <typename T, int N>
struct repeated_vector{
  typedef typename mpl::copy<
      mpl::range_c<int, 0, N>,
      mpl::inserter<
        mpl::vector<>,
        mpl::push_back<mpl::_1, T>
      >
  >::type type;
};

} // end namespace mpl_ext

////////////////////////////////////////////////////////////////////////////////

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
for_each(FuncT, tuple_types&...)
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

////////////////////////////////////////////////////////////////////////////////

// TODO This has a serious bug.  STL maps do not support simultaneous writes, even of different objects.
//        Thus, this class needs to lock the whole value map when writing, killing efficiency.
//        A workaround would be to add a boolean template parameter that only allows adding or changing
//        of keys in a setup phase, similarly to how mutex_map_mutexes_ functions now.  A more
//        perminant solution is going to be using Intel's Thread Building Blocks (TBB) library,
//        with the current implementation as a fallback when the compilation doesn't have TBB
//        available.

/**
 * A cache of objects that can be safely accessed concurrently by threads that share memory.
 */
template <
    typename val_type,
    bool row_locked,
    bool shared_read,
    typename... key_types
>
class ConcurrentCacheBase {

  public:
    template<typename K, typename V> using map_type = std::map<K, V>;
    //----------------------------------------//
    // types for keys
    template<typename... Types> using tuple_type = boost::tuple<Types...>;
    typedef tuple_type<key_types...> key_tuple;
    //----------------------------------------//
    // Types for the mutexes and locks on mutexes
    typedef typename mpl::if_c<
        shared_read,
        shared_mutex,
        mutex
    >::type mutex_type;
    template<int keynum>
    struct mutex_map {
      typedef typename std::map<
        mpl::at<mpl::vector<key_types...>, mpl::int_<keynum>>,
        mutex_type
      > type;
    };
    typedef typename mpl::if_c<
        shared_read,
        shared_lock<mutex_type>,
        unique_lock<mutex_type>
    >::type read_lock;
    typedef unique_lock<mutex_type> write_lock;
    template<typename lock_type> using lock_set_container = std::vector<lock_type>;
    template<typename lock_type> using lock_set = typename mpl::if_c<
        row_locked,
        lock_set_container<lock_type>,
        lock_type
    >::type;
    typedef lock_set<read_lock> read_lock_set;
    typedef lock_set<write_lock> write_lock_set;
    typedef typename mpl::if_c<
      row_locked,
      typename mpl_ext::apply<
        tuple_type,
        typename mpl::transform<mpl::vector<key_types...>, map_type<mpl::_, mutex_type>>::type
      >::type,
      map_type<key_tuple, mutex_type>
    >::type mutex_map_set;
    //----------------------------------------//
    // types for mutexes to lock the mutex_map
    typedef upgrade_mutex mutex_mutex_type;
    typedef typename mpl::if_c<
      row_locked,
      typename mpl_ext::apply<
        tuple_type,
        typename mpl_ext::repeated_vector<
          mutex_mutex_type,
          mpl::size<mpl::vector<key_types...>>::value
        >::type
      >::type,
      mutex_mutex_type
    >::type mutex_set;
    //----------------------------------------//
    // types for map of values
    typedef val_type value_type;
    typedef map_type<key_tuple, value_type> value_map;


    // constness may change if LRU capabilities are implemented
    bool contains(key_types... keys) const {
      key_tuple k(boost::forward<key_types>(keys)...);
      read_lock_set locks = get_read_locks(k);
      auto found_spot = cached_values_.find(k);
      return found_spot != cached_values_.end();
    }

    // constness may change if LRU capabilities are implemented
    boost::optional<value_type> get(key_types... keys) const {
      key_tuple k(boost::forward<key_types>(keys)...);
      read_lock_set locks = get_read_locks(k);
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
      // TODO We should probably use upgrade locks and upgrade mutexes in this method,
      //   but that would require implementing a separate get_locks function
      //   for the upgrade_to_unique_locks, which takes a little bit of effort
      key_tuple k(boost::forward<key_types>(keys)...);
      {
        read_lock_set locks = get_read_locks(k);
        auto found_spot = cached_values_.find(k);
        if(found_spot != cached_values_.end()){
          return found_spot->second;
        }
      } // Release the read lock
      //----------------------------------------//
      // If we haven't returned at this point, we need to
      //   compute the value
      {
        // All threads that get here will have to wait for the one thread doing
        //   the computing.  Threads that call get() after the compute has
        //   completed will exit after the read lock section above.
        write_lock_set locks = get_write_locks(k);
        // Check again to make sure another thread didn't get here before us.
        //   This is known as the "singleton with double-checked locking pattern"
        //   or simply the "double-checked locking pattern"
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
      write_lock_set locks = get_write_locks(k);
      cached_values_[k] = val;
    }

    void
    init_mutex_impl(
        key_types... keys
    )
    {
      key_tuple key_tup(boost::forward<key_types>(keys)...);
      write_lock mtxs_lock(mutex_map_mutexes_);
      auto found_spot = mutexes_.find(key_tup);
      if(found_spot == mutexes_.end()){
        mutexes_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(key_tup),
            std::forward_as_tuple()
        );
      }
    }

    void init_mutex(
        key_types... keys
    )
    {
      init_mutex_impl(keys..., boost::mpl::bool_<row_locked>());
    }

  protected:

    // map of keys to mutexes.  If row_locked is true, this is a
    //   tuple of maps to mutexes, one for each key type
    mutex_map_set mutexes_;

    // Creation locks for the mutexes in the (tuple of) mutex maps
    mutex_set mutex_map_mutexes_;

    // The actual map from key tuples to values
    value_map cached_values_;

    template<typename row_locked_type>
    inline typename boost::enable_if<row_locked_type, void>::type
    init_mutex_impl(key_types... keys, row_locked_type row_locked_true)
    {
      key_tuple key_tup(boost::forward<key_types>(keys)...);
      boost_ext::for_each(mutex_initializer(), key_tup, mutex_map_mutexes_, mutexes_);
    }

    template<typename row_locked_type>
    inline typename boost::enable_if<mpl::not_<row_locked_type>, void>::type
    init_mutex_impl(key_types... keys, row_locked_type row_locked_false)
    {
      key_tuple key_tup(boost::forward<key_types>(keys)...);
      unique_lock<mutex_mutex_type> mtxs_lock(mutex_map_mutexes_);
      auto found_spot = mutexes_.find(key_tup);
      if(found_spot == mutexes_.end()){
        mutexes_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(key_tup),
            std::forward_as_tuple()
        );
      }
    }

    template<typename lock_set_type>
    inline typename boost::enable_if_c<row_locked, lock_set_type>::type
    get_locks(const key_tuple& key_tup){
      lock_set_type locks;
      {
        boost_ext::for_each(
            key_row_locker<decltype(locks)>(&locks),
            key_tup,
            mutex_map_mutexes_,
            mutexes_
        );
        boost::lock(locks.begin(), locks.end());
      }
      // This return will use move semantics, since lock_set_type won't be copyable
      return locks;
    }

    template<typename lock_set_type>
    inline typename boost::enable_if_c<not row_locked, lock_set_type>::type
    get_locks(const key_tuple& key_tup){
      upgrade_lock<mutex_mutex_type> mtxs_read_lock(mutex_map_mutexes_);
      auto found_spot = mutexes_.find(key_tup);
      if(found_spot == mutexes_.end()){
        #if ALLOW_DYNAMIC_MUTEX_CREATION == 0
        assert(false);
        #else
        upgrade_to_unique_lock<mutex_mutex_type> mtxs_write_lock(mtxs_read_lock);
        // Check again to make sure another thread didn't get here before us.
        //   This is known as the "singleton with double-checked locking pattern"
        //   or simply the "double-checked locking pattern"
        auto found_spot = mutexes_.find(key_tup);
        if(found_spot == mutexes_.end()){
          mutexes_.emplace(
              std::piecewise_construct,
              std::forward_as_tuple(key_tup),
              std::forward_as_tuple()
          );
        }
        #endif
      }
      // No need to defer locking, since we're only locking one thing
      lock_set_type lock(mutexes_[key_tup]);
      // This return will use move semantics, since lock_set_type won't be copyable
      return lock;
    }

    inline read_lock_set get_read_locks(const key_tuple& key_tup){
      return get_locks<read_lock_set>(key_tup);
    }

    inline write_lock_set get_write_locks(const key_tuple& key_tup){
      return get_locks<write_lock_set>(key_tup);
    }

    struct mutex_initializer {
      template<typename key_type, typename mutex_map_mtx_type, typename mutex_map_type>
      void operator()(
          key_type& key,
          mutex_map_mtx_type& mtx_map_mtx,
          mutex_map_type& mtx_map
      ) const
      {
        unique_lock<mutex_mutex_type> mtxs_lock(mtx_map_mtx);
        auto found_spot = mtx_map.find(key);
        if(found_spot == mtx_map.end()){
          mtx_map.emplace(
              std::piecewise_construct,
              std::forward_as_tuple(key),
              std::forward_as_tuple()
          );
        }
      }
    };

    template <typename lock_vector>
    struct key_row_locker {
      private:
        lock_vector* lockv_;
      public:
        key_row_locker(lock_vector* v) : lockv_(v) { }

        template<typename key_type, typename mutex_map_mtx_type, typename mutex_map_type>
        void operator()(
            key_type& key,
            mutex_map_mtx_type& mtx_map_mtx,
            mutex_map_type& mtx_map
        ) const
        {
          upgrade_lock<mutex_mutex_type> mtxs_read_lock(mtx_map_mtx);
          auto found_spot = mtx_map.find(key);
          if(found_spot == mtx_map.end()){
            #if ALLOW_DYNAMIC_MUTEX_CREATION == 0
            assert(false);
            #else
            upgrade_to_unique_lock<mutex_mutex_type> mtxs_write_lock(mtxs_read_lock);
            // Check again to make sure another thread didn't get here before us.
            //   This is known as the "singleton with double-checked locking pattern"
            //   or simply the "double-checked locking pattern"
            auto found_spot = mtx_map.find(key);
            if(found_spot == mtx_map.end()){
              mtx_map.emplace(
                  std::piecewise_construct,
                  std::forward_as_tuple(key),
                  std::forward_as_tuple()
              );
            }
            #endif
          }

          (*lockv_).emplace_back(
              mtx_map[key],
              defer_lock
          );
        }
    };

};

/*
 * Different ways of dealing with the problem that mutexes can be kind of
 *   large (shared_mutex is 408 bytes in my boost implementation).  Either
 *   we can lock a mutex for each key ("row locking") and thus have
 *   k*n mutexes, or we can lock a mutex for each tuple of keys
 *   ("element locking").  We can further reduce memory cost by
 *   using a plain mutex ("exclusive read") rather than a shared_mutex
 *   (boost::mutex is only 64 bytes on my machine), which may
 *   speed things up.
 */

template <typename val_type, typename... key_types>
class RowLockedSharedReadCache
    : public ConcurrentCacheBase<val_type, true, true, key_types...>
{ };

template <typename val_type, typename... key_types>
class ElementLockedSharedReadCache
    : public ConcurrentCacheBase<val_type, false, true, key_types...>
{ };

template <typename val_type, typename... key_types>
class RowLockedExclusiveReadCache
    : public ConcurrentCacheBase<val_type, true, false, key_types...>
{ };

template <typename val_type, typename... key_types>
class ElementLockedExclusiveReadCache
    : public ConcurrentCacheBase<val_type, false, false, key_types...>
{ };



} // end namespace sc
#endif /* _util_container_conc_cache_h */
