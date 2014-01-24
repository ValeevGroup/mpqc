//
// consumableresources.h
//
// Copyright (C) 2010 Edward Valeev
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

#ifndef _mpqc_src_lib_util_misc_consumableresources_h
#define _mpqc_src_lib_util_misc_consumableresources_h

#include <map>
#include <limits>
#include <assert.h>
#include <iterator>
#include <mpqc_config.h>
#include <util/keyval/keyval.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/misc/scexception.h>
#include <util/group/thread.h>
#include <util/misc/bug.h>

namespace sc {

  /// ConsumableResources keeps track of consumable resources (memory, disk).
  class ConsumableResources : virtual public SavableState {
    public:
      /** A KeyVal constructor is used to generate a ConsumableResources
          object from the input. The full list of keywords
          that are accepted is below.

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

          <tr><td><tt>memory</tt><td>integer<td>1GB<td>number of bytes; user is allowed to use KB/MB/GB abbreviations \sa KeyValValuestring::sizevalue()

          <tr><td><tt>disk</tt><td>[string integer] pair<td>["./" 0]<td>specifies location of scratch files and available storage in bytes ("0" means unlimited)

          </table>
       */
      ConsumableResources(const Ref<KeyVal>& kv);
      ConsumableResources(StateIn&);
      ConsumableResources();
      ~ConsumableResources();
      void save_data_state(StateOut&);

      //@{
      /// how much resource was given
      size_t max_memory() const;
      size_t max_disk() const;
      //@}

      //@{
      /// how much resource is currently available
      size_t memory() const;
      size_t disk() const;
      //@}

      //@{
      /// consume resource, may throw LimitExceeded<size_t> if not enough available
      void consume_memory(size_t value);
      void consume_disk(size_t value);
      //@}

      //@{
      /// release resouce, may throw ProgrammingError if releasing more resource than how much has been consumed to this point
      void release_memory(size_t value);
      void release_disk(size_t value);
      //@}

      /// UNIX path (absolute or relative) to the disk resource
      const std::string& disk_location() const;

      /** Create a ConsumableResources object.  This routine looks for a -resource
          argument, then the environmental variable MPQC_RESOURCES.
          The argument to -resources should be a string for the KeyVal constructor. */
      static Ref<ConsumableResources> initial_instance(int &argc, char **argv);
      /// Specifies a new default ConsumableResources
      static void set_default_instance(const Ref<ConsumableResources>&);
      /// Returns the default ConsumableResources object
      static const Ref<ConsumableResources>& get_default_instance();

      /// @brief prints ConsumableResources
      ///
      /// @param os output stream
      /// @param print_state set to true to print out currently available resources
      /// @param print_stats set to true to print out maximum resource usage
      void print_summary(std::ostream& os = ExEnv::out0(),
                         bool print_state = true,
                         bool print_stats = false) const;
      /// prints short definition to a string \sa print_summary()
      std::string sprint() const;

      /// allocate array of T size elements long using operator new[] (keeps track of memory)
      template <typename T> T* allocate(std::size_t size) {
        if (size == 0)  return 0;

        ThreadLockHolder lh(lock_);

        T* array = 0;
        try {
          array = (size > 1) ? new T[size] : new T;
        }
        catch (std::bad_alloc&) {
          std::ostringstream oss;
          oss << "ConsumableResources::allocate(size=" << size << "): allocation failed";
          this->print_summary(ExEnv::out0(), true, true);
          throw MemAllocFailed(oss.str().c_str(),__FILE__,__LINE__,size*sizeof(T));
        }

        if (do_profile()) {
          size *= sizeof(T);
          try {
            consume_memory_(size);
          }
          catch (LimitExceeded<size_t>&) {
            if (size > sizeof(T)) delete[] array;
            else delete array;
            this->print_summary(ExEnv::out0(), true, true);
            throw;
          }

          void* array_ptr = static_cast<void*>(array);
          if (debug()) {
            ExEnv::out0() << indent << "ConsumableResources::allocate(size="
                << size << ") => array=" << array_ptr << std::endl;
            // make sure the pointer is not managed (may happen if delete was called directly, not deallocate on a managed buffer before)
            std::map<void*, ResourceAttribites>::iterator pos =
                managed_arrays_.find(array_ptr);
            if (pos != managed_arrays_.end()) {
              ExEnv::out0() << indent << "WARNING: " << array_ptr
                  << " is on the list of managed buffers. size="
                  << managed_arrays_[array_ptr].size << std::endl;
            }
          }
          ResourceAttribites attr(size);
          managed_arrays_[array_ptr] = attr;
        }

        return array;
      }
      //@{
      /// deallocate array of T that was allocated using ConsumableResources::allocate() using operator delete[] (keeps track of memory)
      /// will throw ProgrammingError if this array is not managed by ConsumableResources (i.e. not allocated using allocate() )
      template <typename T> void deallocate(T* const & array) {
        if (object_is_gone()) return;

        if (not do_profile()) {
          delete[] array;
          return;
        }

        if (array != 0) {
          ThreadLockHolder lh(lock_);
          void* array_ptr = static_cast<void*>(array);
          // make sure it's managed by me
          std::map<void*, ResourceAttribites>::iterator pos = managed_arrays_.find(array_ptr);
          if (pos != managed_arrays_.end()) {
            const size_t size = pos->second.size;
            if (size > sizeof(T)) delete[] array;
            else delete array;
            release_memory_(size);
            if (debug()) {
              ExEnv::out0() << indent << "ConsumableResources::deallocate(array=" << array_ptr << "): size=" << size << std::endl;
            }
            managed_arrays_.erase(pos);
          }
          else
            throw ProgrammingError("ConsumableResources::deallocate() -- non-managed array given",
                                   __FILE__, __LINE__, class_desc());
        }
      }
      /// same as  before, but will set array to 0 after deallocation
      template <typename T> void deallocate(T*& array) {
        if (object_is_gone()) return;

        if (array != 0) {
          deallocate<T>(const_cast<T* const &>(array));
          array = 0;
        }
      }
      //@}


      /// adds array to the list of managed arrays and decrements the memory counter

      /// this is useful if need to do special memory allocation, e.g. in pinned memory buffer, etc.
      template <typename T> void manage_array(T* const & array, std::size_t size) {
        if (not do_profile()) return;

        if (array != 0) {
          ThreadLockHolder lh(lock_);
          size *= sizeof(T);
          void* array_ptr = static_cast<void*>(array);
          // make sure it's NOT managed by me
          std::map<void*, ResourceAttribites>::iterator pos = managed_arrays_.find(array_ptr);
          if (pos != managed_arrays_.end()) {
            std::ostringstream oss;
            oss << indent << "ConsumableResources::manage_array() -- given managed array (ptr="
                << array_ptr << " size=" << size << ")";
            throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__, class_desc());
          }

          try {
            consume_memory_(size);
          }
          catch (LimitExceeded<size_t>&) {
            this->print_summary(ExEnv::out0(), true, true);
            throw;
          }

          ResourceAttribites attr(size);
          managed_arrays_[array_ptr] = attr;
          if (debug()) {
            ExEnv::out0() << indent << "ConsumableResources::manage_array(array=" << array_ptr << ", size=" << size << ")" << std::endl;
          }
        }
      }
      /// removes array to the list of managed arrays and increments the memory counter
      template <typename T> void unmanage_array(T* const & array) {
        if (object_is_gone()) return;
        if (not do_profile()) return;

        if (array != 0) {
          ThreadLockHolder lh(lock_);
          void* array_ptr = static_cast<void*>(array);
          // make sure it's managed by me
          std::map<void*, ResourceAttribites>::iterator pos = managed_arrays_.find(array_ptr);
          if (pos != managed_arrays_.end()) {
            const size_t size = pos->second.size;
            release_memory_(size);
            if (debug()) {
              ExEnv::out0() << indent << "ConsumableResources::unmanage_array(array=" << array_ptr << ": size=" << size << ")" << std::endl;
            }
            managed_arrays_.erase(pos);
          }
          else {
            void* array_ptr = static_cast<void*>(array);
            std::ostringstream oss;
            oss << "ConsumableResources::unmanage_array() -- given non-managed array (ptr="
                << array_ptr << " size=" << pos->second.size << ")";
            throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__, class_desc());
          }
        }
      }

      // need to check and report any anomalous events?
      static bool debug() { return MPQC_MEMORY_CHECK >= 3; }
      // need to profile resource usage?
      static bool do_profile() { return MPQC_MEMORY_CHECK >= 1; }

    private:
      static ClassDesc class_desc_;

      //@{ these do not lock
      /// consume resource, may throw LimitExceeded<size_t> if not enough available
      void consume_memory_(size_t value);
      void consume_disk_(size_t value);
      //@}

      //@{ these do not lock
      /// release resource, may throw ProgrammingError if releasing more resource than how much has been consumed to this point
      void release_memory_(size_t value);
      void release_disk_(size_t value);
      //@}

      /// return true, and print a warning to stderr, if this object is gone
      bool object_is_gone();
      /// return true, and print a warning to stderr, if the default instance of this class is gone
      static bool default_object_is_gone();

      /**
       * Checks that there are no outstanding debts
       * @param os output stream to which the info is written (default = sc::ExEnv::err0())
       */
      void summarize_unreleased_resources(std::ostream& os = sc::ExEnv::err0()) const;

      struct defaults {
          static size_t memory;
          static std::pair<std::string,size_t> disk;
      };

      template <typename T> class ResourceCounter {
        public:
          ResourceCounter() :
            max_value_(),
            value_(),
            lowest_value_()
            {
            }
          ResourceCounter(const T& max_value) :
            max_value_(max_value),
            value_(max_value),
            lowest_value_(max_value)
            {
            }
          ResourceCounter(const T& max_value, const T& value) :
            max_value_(max_value),
            value_(value),
            lowest_value_(value)
            {
             MPQC_ASSERT(value <= max_value);
            }
          ResourceCounter(const ResourceCounter& other) :
            max_value_(other.max_value_),
            value_(other.value_),
            lowest_value_(other.lowest_value_)
          {
          }
          ResourceCounter& operator=(const ResourceCounter& other) {
            max_value_ = other.max_value_;
            value_ = other.value_;
            lowest_value_ = other.lowest_value_;
            return *this;
          }
          operator T() const { return value_; }
          ResourceCounter& operator+=(const T& val) { value_ = std::min(max_value_, value_ + val); return *this; }
          // nonthrowing
          ResourceCounter& operator-=(const T& val) {
            value_ = std::max(T(0), value_ - val);
            lowest_value_ = std::min(value_,lowest_value_);
            return *this;
          }

          const T& max_value() const { return max_value_; }
          const T& value() const { return value_; }
          const T& lowest_value() const { return lowest_value_; }

          static std::string value_to_string(T value) {
            return to_string(value, true);
          }
          static std::string difference_to_string(T value) {
            return to_string(value, false);
          }

          void operator &(StateIn& s) {
            s.get(max_value_);
            s.get(value_);
            s.get(lowest_value_);
          }
          void operator &(StateOut& s) {
            s.put(max_value_);
            s.put(value_);
            s.put(lowest_value_);
          }

        private:
          T max_value_;
          T value_;
          T lowest_value_;  //< keeps track of the lowest value of the resource

          static std::string to_string(T t, bool may_be_unlimited) {
            const int prec = 3; // print this many digits

            if (may_be_unlimited && t == std::numeric_limits<T>::max())
              return "unlimited";

            // determine m such that 1000^m <= t <= 1000^(m+1)
            char m = 0;
            double thousand_m = 1;
            double thousand_mp1 = 1000;
            while (t >= thousand_mp1) {
              ++m;
              thousand_m = thousand_mp1;
              thousand_mp1 *= 1000;
            }

            // determine units
            std::string unit;
            switch (m) {
              case 0: unit = "B"; break;
              case 1: unit = "kB"; break;
              case 2: unit = "MB"; break;
              case 3: unit = "GB"; break;
              case 4: unit = "TB"; break;
              case 5: unit = "PB"; break;
              case 6: unit = "EB"; break;
              case 7: unit = "ZB"; break;
              case 8: unit = "YB"; break;
              default: MPQC_ASSERT(false); break;
            }

            // compute normalized mantissa
            std::ostringstream oss;
            oss.precision(prec);
            oss << (double)t/(double)thousand_m << unit;
            return oss.str();
          }
      };

      typedef ResourceCounter<size_t> rsize;
      rsize memory_;
      std::pair<std::string, rsize> disk_;

      /// Describes attributes of a resource (size, where allocated, etc.)
      struct ResourceAttribites {
          ResourceAttribites() : size(0) {
          }
          ResourceAttribites(std::size_t s) : size(s) {
          }
          std::size_t size;
#if MPQC_MEMORY_CHECK >= 2
          Debugger::Backtrace backtrace;
#endif
          operator std::string() const {
            std::ostringstream oss;
            oss << "size=" << size;
#if MPQC_MEMORY_CHECK >= 2
            const size_t nframes_to_skip = 1;
            oss << " allocated at:" << std::endl << backtrace.str(nframes_to_skip);
#endif
            return oss.str();
          }
      };

      /// this keeps track of arrays of data explicitly managed by ConsumableResources
      /// value is the size of array in bytes
      std::map<void*, ResourceAttribites> managed_arrays_;
      Ref<ThreadLock> lock_; // the lock used to protect the map and resource counters

      static Ref<ConsumableResources> default_instance_;

  };

  //@{
  /// allocate and deallocate array of data using new or new[] (delete or delete[]) and using default ConsumableResources object
  template <typename T> T* allocate(std::size_t size) {
    return ConsumableResources::get_default_instance()->allocate<T>(size);
  }
  /// this version will set array to 0 upon return \sa deallocate
  template <typename T> void deallocate(T*& array) {
    // deallocation after this is gone is OK, .pointer() ensures that dereferencing does not abort
    ConsumableResources::get_default_instance().pointer()->deallocate(array);
  }
  template <typename T> void deallocate(T* const & array) {
    // deallocation after this is gone is OK, .pointer() ensures that dereferencing does not abort
    ConsumableResources::get_default_instance().pointer()->deallocate(array);
  }
  //@}

  //@{
  /// manage or unmanaged array of data using default ConsumableResources object
  template <typename T> void manage_array(T* const & array, std::size_t size) {
    ConsumableResources::get_default_instance()->manage_array(array, size);
  }
  template <typename T> void unmanage_array(T* const & array) {
    ConsumableResources::get_default_instance()->unmanage_array(array);
  }


} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
