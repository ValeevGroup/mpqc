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

#ifndef _util_container_conc_cache_h
#define _util_container_conc_cache_h

// Standard library includes
#include <cassert>
#include <utility>
#include <functional>

// Boost includes
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/utility.hpp>
#include <boost/type_traits.hpp>
#include <boost/thread/future.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/replace_if.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/back_inserter.hpp>

// Madness includes
#include <madness/world/worldhashmap.h>
#include <madness/world/future.h>

// MPQC imcludes
#include <util/misc/iterators.h>
#include <util/misc/meta.h>
#include <util/container/conc_cache_fwd.h>

namespace {

  template<typename T>
  struct _hash { };

  template<typename... key_types>
  struct _hash<boost::tuple<key_types...>> {
      inline madness::hashT
      operator()(const boost::tuple<key_types...>& key) const {
        madness::hashT seed = 0;
        hash_recursive(key, seed);
        return seed;
      }
    private:

      template<typename H, typename T>
      inline typename boost::disable_if<
        boost::is_same<T, boost::tuples::null_type>>::type
      hash_recursive(const boost::tuples::cons<H, T>& key_part, madness::hashT& seed) const
      {
        madness::hash_combine(seed, make_hashable(key_part.get_head()));
        hash_recursive(key_part.get_tail(), seed);
      }

      template<typename H, typename T>
      inline typename boost::enable_if<boost::is_same<T, boost::tuples::null_type>>::type
      hash_recursive(const boost::tuples::cons<H, T>& key_part, madness::hashT& seed) const
      {
        madness::hash_combine(seed, make_hashable(key_part.get_head()));
      }

      template<typename T>
      inline typename boost::enable_if<
        boost::is_enum<T>,
        int
      >::type
      make_hashable(const T& val) const {
        return (int)val;
      }

      template<typename T>
      inline typename boost::disable_if<
        boost::is_enum<T>,
        const T&
      >::type
      make_hashable(const T& val) const {
        return val;
      }
  };

}

//////////////////////////////////////////////////////////////////////////////////

namespace sc {

template<int n_indices, int... PermutedIndices>
struct KeyPermutation {
    static constexpr int order = sizeof...(PermutedIndices);
};

template<int n_indices, int idx_1, int idx_2>
struct KeyTransposition
  : public meta::splat_values<int,
      KeyPermutation,
      typename boost::mpl::replace_if<
        typename boost::mpl::replace_if<
          typename boost::mpl::copy<
            boost::mpl::range_c<int, 0, n_indices>,
            boost::mpl::back_inserter<
              boost::mpl::vector_c<int, n_indices>
            >
          >::type,
          boost::mpl::equal_to<
            boost::mpl::_,
            boost::mpl::int_<idx_1>
          >,
          boost::mpl::int_<idx_2>
        >::type,
        boost::mpl::equal_to<
          boost::mpl::_,
          boost::mpl::int_<idx_2>
        >,
        boost::mpl::int_<idx_1>
      >::type
    >::type
{
    static constexpr int first_index = idx_1;
    static constexpr int second_index = idx_2;
};

template <int nkeys>
struct IdentityKeyPermutation
  : public meta::splat_values<int,
      KeyPermutation,
      typename boost::mpl::copy<
        boost::mpl::range_c<int, 0, nkeys>,
        boost::mpl::back_inserter< boost::mpl::vector_c<int, nkeys> >
      >::type
    >::type
{ };

template<typename... Permutations>
struct KeySymmetry {
    static constexpr int n_permutations = sizeof...(Permutations);
};


template<
  typename ValueType,
  typename... key_types
>
class ConcurrentCacheBase {

  public:

    /// types for keys
    //@{
    /// tuple type alias (may use std::tuple in the future)
    template<typename... Types> using tuple_type = boost::tuple<Types...>;
    /// The type used as a key into the map
    typedef tuple_type<key_types...> key_tuple;
    //@}

    /// types for map of values
    //@{
    /// The type of the cached value; returned by the get() functions
    typedef ValueType value_type;
    /// A future of a value, used to prevent recomputation
    typedef boost::shared_future<value_type> future_value;
    /// A map from keys to futures
    typedef madness::ConcurrentHashMap<key_tuple, future_value, _hash<key_tuple>> future_map;
    /// The accessor type used by the madness ConcurrentHashMap
    typedef typename future_map::accessor future_map_accessor;
    /// The const accessor used by the madness ConcurrentHashMap
    typedef typename future_map::const_accessor future_map_const_accessor;
    //@}


    ConcurrentCacheBase() : cached_values_()
    { }

    ConcurrentCacheBase(int nbins) : cached_values_(nbins)
    { }

    // Old version, provided for backwards compatibility
    value_type get(
        key_types... keys,
        const std::function<value_type()>& compute_fxn
    )
    {

      key_tuple k(std::forward<key_types>(keys)...);
      future_map_accessor a;

      if(cached_values_.insert(a, std::make_pair(k, future_value()))){
        boost::promise<value_type> p;
        p.set_value(compute_fxn());
        a->second = p.get_future().share();
      }
      return a->second.get();

    }

    value_type get(
        key_types... keys,
        const std::function<value_type(const key_tuple&)>& compute_fxn
    )
    {

      key_tuple k(std::forward<key_types>(keys)...);
      future_map_accessor a;

      if(cached_values_.insert(a, std::make_pair(k, future_value()))){
        boost::promise<value_type> p;
        p.set_value(compute_fxn(k));
        a->second = p.get_future().share();
      }
      return a->second.get();

    }

    void clear()
    {
      cached_values_.clear();
    }

  protected:

    // The actual map from key tuples to values
    future_map cached_values_;

};

/**
 * A cache of objects that can be safely accessed concurrently by threads that share memory.
 * Some cached objects may be transformable into objects corresponding to other keys via
 * a transformation.
 */
template <
    typename val_type,
    typename symmetry,
    typename... key_types
>
class ConcurrentCacheWithSymmetry
  : public ConcurrentCacheBase<val_type, key_types...>
{
  public:
    typedef ConcurrentCacheBase<val_type, key_types...> super_t;
    using ConcurrentCacheBase<val_type, key_types...>::ConcurrentCacheBase;

    typename super_t::value_type get_or_permute(
        key_types... keys,
        const std::function<typename super_t::value_type(
            const typename super_t::key_tuple&
        )>& compute_fxn,
        const std::function<typename super_t::value_type(
            const typename super_t::value_type&,
            const typename super_t::key_tuple&
        )>& permute_fxn
    );
};

/**
 * Specialization for the identity
 */
template <typename val_type, typename... key_types>
class ConcurrentCacheWithSymmetry<
  val_type,
  KeySymmetry<IdentityKeyPermutation<sizeof...(key_types)>>,
  key_types...
> : public ConcurrentCacheBase<val_type, key_types...>
{
  public:
    typedef ConcurrentCacheBase<val_type, key_types...> super_t;
    using super_t::ConcurrentCacheBase;

    typename super_t::value_type get_or_permute(
        key_types... keys,
        const std::function<typename super_t::value_type(
            const typename super_t::key_tuple&
        )>& compute_fxn,
        const std::function<typename super_t::value_type(
            const typename super_t::value_type&,
            const typename super_t::key_tuple&
        )>& permute_fxn
    )
    {
      return this->get(keys..., compute_fxn);
    }
};

/**
 * Specialization with only one transposition other than the identity
 */
template <typename val_type, int n_keys, int idx1, int idx2, typename... key_types>
class ConcurrentCacheWithSymmetry<
  val_type,
  KeySymmetry<
    IdentityKeyPermutation<sizeof...(key_types)>,
    KeyTransposition<n_keys, idx1, idx2>
  >,
  key_types...
> : public ConcurrentCacheBase<val_type, key_types...>
{
  public:
    typedef ConcurrentCacheBase<val_type, key_types...> super_t;

    using ConcurrentCacheBase<val_type, key_types...>::ConcurrentCacheBase;

    typedef KeyTransposition<n_keys, idx1, idx2> transposition_type;
    typedef typename super_t::key_tuple key_tuple;
    typedef typename super_t::value_type value_type;

    static_assert(
        boost::is_convertible<
          typename boost::mpl::at_c<boost::mpl::vector<key_types...>, idx1>::type,
          typename boost::mpl::at_c<boost::mpl::vector<key_types...>, idx2>::type
        >::value
        and
        boost::is_convertible<
          typename boost::mpl::at_c<boost::mpl::vector<key_types...>, idx2>::type,
          typename boost::mpl::at_c<boost::mpl::vector<key_types...>, idx1>::type
        >::value
        and
        boost::has_less<
          typename boost::mpl::at_c<boost::mpl::vector<key_types...>, idx1>::type,
          typename boost::mpl::at_c<boost::mpl::vector<key_types...>, idx2>::type
        >::value
        ,
        "Types must be compatible to have symmetry in ConcurrentCacheWithSymmetry"
    );

    inline value_type
    get_or_permute(
        key_types... keys,
        const std::function<value_type(
            const key_tuple&
        )>& compute_fxn,
        const std::function<value_type(
            const value_type&,
            const key_tuple&
        )>& permute_fxn
    )
    {
      key_tuple k(std::forward<key_types>(keys)...);

      if(not (boost::tuples::get<idx1>(k) < boost::tuples::get<idx2>(k))) {

        // regular case
        typename super_t::future_map_accessor a;
        if(this->cached_values_.insert(a, std::make_pair(k, typename super_t::future_value()))){
          boost::promise<value_type> p;
          p.set_value(
              compute_fxn(k)
          );
          a->second = p.get_future().share();
        }
        return a->second.get();

      }
      else {

        // permuted case
        typename super_t::future_map_accessor a_perm, a_canon;
        if(this->cached_values_.insert(a_perm, std::make_pair(k, typename super_t::future_value()))){

          key_tuple k_canon = make_permuted_key_tuple(k);

          // Now see if we need to comptue the canonical ordered value
          if(this->cached_values_.insert(a_canon, std::make_pair(k_canon, typename super_t::future_value()))){

            boost::promise<value_type> p;
            p.set_value(
                compute_fxn(k_canon)
            );
            a_canon->second = p.get_future().share();

          }

          // Now call the permute function on the result
          boost::promise<value_type> perm_p;
          perm_p.set_value(
              permute_fxn(a_canon->second.get(), k_canon)
          );
          a_perm->second = perm_p.get_future().share();

        }
        return a_perm->second.get();

      }
    }

  private:

    inline key_tuple
    make_permuted_key_tuple(const key_tuple& orig) {
      return make_permuted_key_tuple_recursive(
          static_cast<typename key_tuple::inherited>(orig),
          orig, boost::mpl::int_<0>()
      );
    }

    template<typename I, typename H, typename T>
    inline typename boost::enable_if_c<
      I::value == idx1,
      boost::tuples::cons<H, T>
    >::type
    make_permuted_key_tuple_recursive(
        const boost::tuples::cons<H, T>& kpart,
        const key_tuple& orig, I not_used
    )
    {
      return boost::tuples::cons<H, T>(
          boost::tuples::get<idx2>(orig),
          make_permuted_key_tuple_recursive(
              kpart.get_tail(), orig,
              typename boost::mpl::plus<I, boost::mpl::int_<1>>::type()
          )
      );
    }

    template<typename I, typename H, typename T>
    inline typename boost::enable_if_c<
      I::value == idx2,
      boost::tuples::cons<H, T>
    >::type
    make_permuted_key_tuple_recursive(
        const boost::tuples::cons<H, T>& kpart,
        const key_tuple& orig, I not_used
    )
    {
      return boost::tuples::cons<H, T>(
          boost::tuples::get<idx1>(orig),
          make_permuted_key_tuple_recursive(
              kpart.get_tail(), orig,
              typename boost::mpl::plus<I, boost::mpl::int_<1>>::type()
          )
      );
    }

    template<typename I, typename H, typename T>
    inline typename boost::disable_if_c<
      I::value == idx1 or I::value == idx2,
      boost::tuples::cons<H, T>
    >::type
    make_permuted_key_tuple_recursive(
        const boost::tuples::cons<H, T>& kpart,
        const key_tuple& orig, I not_used
    )
    {
      return boost::tuples::cons<H, T>(
          kpart.get_head(),
          make_permuted_key_tuple_recursive(
              kpart.get_tail(), orig,
              typename boost::mpl::plus<I, boost::mpl::int_<1>>::type()
          )
      );
    }

    template<typename I>
    inline const boost::tuples::null_type&
    make_permuted_key_tuple_recursive(
        const boost::tuples::null_type& nt,
        const key_tuple& orig, I not_used
    )
    {
      return nt;
    }

};


} // end namespace sc
#endif /* _util_container_conc_cache_h */
