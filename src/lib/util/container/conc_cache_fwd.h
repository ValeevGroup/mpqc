//
// conc_cache_fwd.h
//
// Forward declarations for conc_cache.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Feb 13, 2014
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

#ifndef _util_misc_conc_cache_fwd_h
#define _util_misc_conc_cache_fwd_h

namespace sc {

template<int n_indices, int... PermutedIndices>
struct KeyPermutation;

template<int n_indices, int idx_1, int idx_2>
struct KeyTransposition;

template<typename... Permutations>
struct KeySymmetry;

template <int nkeys>
struct IdentityKeyPermutation;

template <int nkeys>
using IdentityKeySymmetry =
    KeySymmetry<
      IdentityKeyPermutation<nkeys>
    >;

template <int nkeys, int idx1, int idx2>
using SingleTranspositionKeySymmetry =
    KeySymmetry<
      IdentityKeyPermutation<nkeys>,
      KeyTransposition<nkeys, idx1, idx2>
    >;

template <typename val_type, typename... key_types>
class ConcurrentCacheBase;

template <typename val_type, typename symmetry, typename... key_types>
class ConcurrentCacheWithSymmetry;

template <typename val_type, typename... key_types>
class ConcurrentCacheWithSymmetry<
  val_type,
  KeySymmetry<IdentityKeyPermutation<sizeof...(key_types)>>,
  key_types...
>;

template <typename val_type, int n_keys, int idx1, int idx2, typename... key_types>
class ConcurrentCacheWithSymmetry<
  val_type,
  KeySymmetry<
    IdentityKeyPermutation<sizeof...(key_types)>,
    KeyTransposition<n_keys, idx1, idx2>
  >,
  key_types...
>;

/**
 * A cache of objects that can be safely accessed concurrently by threads that share memory.
 */
template <typename val_type, typename... key_types>
using ConcurrentCache = ConcurrentCacheBase<val_type, key_types...>;


} // end namespace sc

#endif /* _util_misc_conc_cache_fwd_h */
