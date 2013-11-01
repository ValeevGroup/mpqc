//
// registry.timpl.h
//
// Copyright (C) 2008 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_registrytimpl_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_registrytimpl_h

#include <cassert>
#include <iterator>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/misc/registry.h>

namespace sc {

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  Ref< Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual> >
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::instance()
  {
    return CreationPolicy< Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual> >::instance();
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  void
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::save_instance(const Ref<Registry>& obj, StateOut& so)
  {
    CreationPolicy< Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual> >::save_instance(obj, so);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  Ref< Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual> >
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::restore_instance(StateIn& si)
  {
    return CreationPolicy< Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual> >::instance(si);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::Registry() :
    lock_(ThreadGrp::get_default_threadgrp()->new_lock())
  {
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::Registry(StateIn& si) :
    lock_(ThreadGrp::get_default_threadgrp()->new_lock())
  {
    si.get(map_);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  void
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::save_data_state(StateOut& so)
  {
    so.put(map_);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  void
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::clear()
  {
    ThreadLockHolder lh(lock_);
    map_.clear();
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  typename Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::const_iterator
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::find_by_key(const Key& key) const
  {
    // if KeyEqual is std::equal_to<Key> then can use the fast find function
    if (SameType<KeyEqual,std::equal_to<Key> >::result) {
      const_iterator result = map_.find(key);
      return result;
    }
    else {
      KeyEqual keyeq;
      for(const_iterator v=map_.begin(); v!= map_.end(); ++v)
        if (keyeq(v->first,key))
          return v;
      return map_.end();
    }
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  typename Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::iterator
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::find_by_key(const Key& key)
  {
    // if KeyEqual is std::equal_to<Key> then can use the fast find function
    if (SameType<KeyEqual,std::equal_to<Key> >::result) {
      iterator result = map_.find(key);
      return result;
    }
    else {
      KeyEqual keyeq;
      for(iterator v=map_.begin(); v!= map_.end(); ++v)
        if (keyeq(v->first,key))
          return v;
      return map_.end();
    }
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  typename Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::const_iterator
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::find_by_value(const Value& value) const
  {
    const_iterator result;

    ValueEqual valeq;
    for(const_iterator v=map_.begin(); v!= map_.end(); ++v)
      if (valeq(v->second,value))
        return v;

    return map_.end();
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  bool
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::key_exists(const Key& key) const
  {
    bool result = false;

    // although this does not modify the map, cannot search map while someone else is changing it
    const_iterator v = find_by_key(key);
    ThreadLockHolder lh(lock_);
    if (v != map_.end())
      result = true;

    return result;
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  bool
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::value_exists(const Value& value) const
  {
    bool result = false;

    // although this does not modify the map, cannot search map while someone else is changing it
    ThreadLockHolder lh(lock_);
    const_iterator v = find_by_value(value);
    if (v != map_.end())
      result = true;

    return result;
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  const Key&
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::key(const Value& value) const
  {
    // although this does not modify the map, cannot search map while someone else is changing it
    ThreadLockHolder lh(lock_);
    const_iterator v = find_by_value(value);
    if (v != map_.end()) {
      return v->first;
    }
    else {
      throw not_found("value not found");
    }
    // unreachable
    assert(false);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  const Value&
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::value(const Key& key) const
  {
    // although this does not modify the map, cannot search map while someone else is changing it
    ThreadLockHolder lh(lock_);
    const_iterator v = find_by_key(key);
    if (v != map_.end()) {
      return v->second;
    }
    else {
      lh.unlock();
      this->print(ExEnv::out0());
      throw not_found("key not found");
    }
    // unreachable
    assert(false);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  void
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::add(const Key& key,
                                          const Value& value)
  {
    // check if key already exists
    if (key_exists(key)) {
      this->print(ExEnv::out0());
      throw std::logic_error("key already exists");
    }
    ThreadLockHolder lh(lock_);
    map_[key] = value;
  }
  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  void
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::add(const std::pair<Key,Value>& keyval_pair)
  {
    this->add(keyval_pair.first,keyval_pair.second);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  void
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::remove(const Key& key)
  {
    ThreadLockHolder lh(lock_);
    iterator v = find_by_key(key);
    if (v != map_.end())
      map_.erase(v);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy, typename KeyEqual, typename ValueEqual >
  void
  Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual>::print(std::ostream& os) const
  {
    ThreadLockHolder lh(lock_);
    for(typename Map::const_iterator iter = map_.begin();
        iter != map_.end();
        ++iter) {
      os << iter->first << " " << iter->second << std::endl;
    }
  }

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
