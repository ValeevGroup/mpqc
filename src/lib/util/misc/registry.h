//
// registry.h
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_registry_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_registry_h

#include <map>
#include <stdexcept>
#include <util/group/thread.h>
#include <util/state/state.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>

namespace sc {

  namespace detail {
    /** SingletonCreationPolicy is used to create Singletons.
     */
    template <typename T>
    class SingletonCreationPolicy {
      protected:
        static Ref<T> instance() {
          return instance_;
        }
        static Ref<T> instance(StateIn& si) {
          if (!instance_restored_) {
            instance_ = new T(si);
            instance_restored_ = true;
          }
          return instance_;
        }
        static void save_instance(const Ref<T>& instance, StateOut& so) {
          if (!instance_saved_) {
            instance_->save_data_state(so);
            instance_saved_ = true;
          }
        }
        static Ref<T> instance_;
        static bool instance_restored_;
        static bool instance_saved_;
    };

    template <typename T> Ref<T> SingletonCreationPolicy<T>::instance_(new T);
    template <typename T> bool SingletonCreationPolicy<T>::instance_restored_ = false;
    template <typename T> bool SingletonCreationPolicy<T>::instance_saved_ = false;

    /// NonsingletonCreationPolicy is used to create non-Singletons on heap.
    template<typename T>
    class NonsingletonCreationPolicy {
      protected:
        static Ref<T> instance() {
          return new T;
        }
        static Ref<T> instance(StateIn& si) {
          bool nonnull; si.get(nonnull);
          if (nonnull)
            return new T(si);
          else
            return 0;
        }
        static void save_instance(const Ref<T>& instance, StateOut& so) {
          if (instance) {
            so.put(true);
            instance->save_data_state(so);
          }
          else
            so.put(false);
        }
    };

  } // end of namespace detail

  /** Registry wraps std::map and can be policy-configured to act as a Singleton or a regular object.
      Registry is not a SavableState, but it can be serialized (see save_instance and restore_instance methods).
      Therefore both Key and Value are assumed to be serializable with StateIn and StateOut.

      KeyEqual and ValueEqual are Functor types that define the
      equivalence of two Key (or Value) objects. By default std::equal_to is used.
      For pointer types std::equal_to may not be appropriate since it will not distinguish
      2 equivalent but distinct objects.

      \sa sc::detail::SingletonCreationPolicy \sa sc::detail::NonsingletonCreationPolicy
    */
  template <typename Key, typename Value, template <typename> class CreationPolicy,
            typename KeyEqual = std::equal_to<Key>,
            typename ValueEqual = std::equal_to<Value>
           >
    class Registry : public RefCount, public CreationPolicy< Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual> > {
      public:
        static Ref<Registry> instance();
        static void save_instance(const Ref<Registry>&, StateOut&);
        static Ref<Registry> restore_instance(StateIn&);

        /// erases all entries
        void clear();
        /// removes the object corresponding to key
        void remove(const Key& key);
        /// removes all objects whose keys evaluate predicate to true:  p(key) == true
        template <typename Pred> void remove_if(const Pred& p) {
          ThreadLockHolder lh(lock_);
          typename Map::iterator i = map_.begin();
          for(; i != map_.end(); ) {
                      if (p(*i)) {
                        map_.erase(i++);
                      } else {
                        ++i;
                      }
          }
        }
        /// key exists?
        bool key_exists(const Key& key) const;
        /// value exists?
        bool value_exists(const Value& value) const;
        /// returns object that corresponds to this key. If key is not known, throws
        const Value& value(const Key& key) const;
        /// returns key that corresponds to this object. If obj is not known, throws
        const Key& key(const Value& obj) const;
        /// registers this object
        void add(const Key& key,
                 const Value& obj);
        /// a useful variation of Registry::add(const Key&, const Value&)
        void add(const std::pair<Key,Value>& keyval_pair);

        void print(std::ostream& os = ExEnv::out0()) const;

        class not_found : public std::logic_error {
          public:
            not_found(const char* what) : std::logic_error(what) {}
        };

      private:
        // creation policy must be able to construct Registry
        friend class CreationPolicy< Registry<Key,Value,CreationPolicy,KeyEqual,ValueEqual> >;

        // access only through instance() and related methods
        Registry();
        // use restore_instance
        Registry(StateIn&);
        // use save_instance
        void save_data_state(StateOut&);

        typedef std::map<Key,Value> Map;
        typedef typename Map::const_iterator const_iterator;
        typedef typename Map::iterator iterator;
        Map map_;

        // assumes that map is already locked
        const_iterator find_by_key(const Key& key) const;
        // assumes that map is already locked
        iterator find_by_key(const Key& key);
        // assumes that map is already locked
        const_iterator find_by_value(const Value& value) const;

        // std::map's operations are not reentrant, hence lock the map every time
        Ref<ThreadLock> lock_;

        /// used to compare types
        template <typename T1, typename T2> struct SameType {
          static const bool result = false;
        };
        template <typename T> struct SameType<T,T> {
          static const bool result = true;
        };

    };

} // end of namespace sc

#include <util/misc/registry.timpl.h>

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
