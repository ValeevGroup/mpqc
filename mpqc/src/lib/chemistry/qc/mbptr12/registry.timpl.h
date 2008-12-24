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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_registrytimpl_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_registrytimpl_h

#include <cassert>
#include <chemistry/qc/mbptr12/registry.h>

namespace sc {

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  Ref< Registry<Key,Value,CreationPolicy> >
  Registry<Key,Value,CreationPolicy>::instance()
  {
    return CreationPolicy< Registry<Key,Value,CreationPolicy> >::instance();
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  void
  Registry<Key,Value,CreationPolicy>::save_instance(const Ref<Registry>& obj, StateOut& so)
  {
    CreationPolicy< Registry<Key,Value,CreationPolicy> >::save_instance(obj, so);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  Ref< Registry<Key,Value,CreationPolicy> >
  Registry<Key,Value,CreationPolicy>::restore_instance(StateIn& si)
  {
    return CreationPolicy< Registry<Key,Value,CreationPolicy> >::instance(si);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  Registry<Key,Value,CreationPolicy>::Registry() :
    lock_(ThreadGrp::get_default_threadgrp()->new_lock())
  {
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  Registry<Key,Value,CreationPolicy>::Registry(StateIn& si) :
    lock_(ThreadGrp::get_default_threadgrp()->new_lock())
  {
    si.get(map_);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  void
  Registry<Key,Value,CreationPolicy>::save_data_state(StateOut& so)
  {
    so.put(map_);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  typename Registry<Key,Value,CreationPolicy>::const_iterator
  Registry<Key,Value,CreationPolicy>::find_by_key(const Key& key) const
  {
    const_iterator result = map_.find(key);
    return result;
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  typename Registry<Key,Value,CreationPolicy>::const_iterator
  Registry<Key,Value,CreationPolicy>::find_by_value(const Value& value) const
  {
    const_iterator result;

    for(const_iterator v=map_.begin(); v!= map_.end(); ++v)
      if (v->second == value)
        return v;

    return map_.end();
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  bool
  Registry<Key,Value,CreationPolicy>::key_exists(const Key& key) const
  {
    bool result = false;

    // although this does not modify the map, cannot search map while someone else is changing it
    lock_->lock();
    const_iterator v = find_by_key(key);
    if (v != map_.end())
      result = true;
    lock_->unlock();

    return result;
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  bool
  Registry<Key,Value,CreationPolicy>::value_exists(const Value& value) const
  {
    bool result = false;

    // although this does not modify the map, cannot search map while someone else is changing it
    lock_->lock();
    const_iterator v = find_by_value(value);
    if (v != map_.end())
      result = true;
    lock_->unlock();

    return result;
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  const Key&
  Registry<Key,Value,CreationPolicy>::key(const Value& value) const
  {
    // although this does not modify the map, cannot search map while someone else is changing it
    lock_->lock();
    const_iterator v = find_by_value(value);
    if (v != map_.end()) {
      lock_->unlock();
      return v->first;
    }
    else {
      lock_->unlock();
      throw not_found("value not found");
    }
    // unreachable
    assert(false);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  const Value&
  Registry<Key,Value,CreationPolicy>::value(const Key& key) const
  {
    // although this does not modify the map, cannot search map while someone else is changing it
    lock_->lock();
    const_iterator v = find_by_key(key);
    if (v != map_.end()) {
      lock_->unlock();
      return v->second;
    }
    else {
      lock_->unlock();
      throw not_found("key not found");
    }
    // unreachable
    assert(false);
  }

  template <typename Key, typename Value, template <typename> class CreationPolicy >
  void
  Registry<Key,Value,CreationPolicy>::add(const Key& key,
                                          const Value& value)
  {
    lock_->lock();
    map_[key] = value;
    lock_->unlock();
  }
  template <typename Key, typename Value, template <typename> class CreationPolicy >
  void
  Registry<Key,Value,CreationPolicy>::add(const std::pair<Key,Value>& keyval_pair)
  {
    this->add(keyval_pair.first,keyval_pair.second);
  }

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
