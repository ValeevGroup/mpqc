#ifndef MPQC_UTILITY_MUTEX_HPP
#define MPQC_UTILITY_MUTEX_HPP

#include <boost/thread/mutex.hpp>
#include <util/misc/exenv.h>

namespace mpqc {

    /// @addtogroup CoreUtility
    /// @{

    /// Static mutex factory
    template<typename T>
    struct static_mutex {
        static const bool debug = false; // change to true
        static void lock() {
          if (debug) {
            std::ostringstream oss;
            oss << "entered mpqc::static_mutex<>::lock(): count = " << lock_count << std::endl;
            sc::ExEnv::out0() << oss.str();
          }
          get().lock();
          if (debug) {
            ++lock_count;
            std::ostringstream oss;
            oss << "completed mpqc::static_mutex<>::lock(): count = " << lock_count << std::endl;
            sc::ExEnv::out0() << oss.str();
          }
        }
        static void unlock() {
          if (debug) {
            std::ostringstream oss;
            oss << "entered mpqc::static_mutex<>::unlock(): count = " << lock_count << std::endl;
            sc::ExEnv::out0() << oss.str();
          }
          get().unlock();
          if (debug) {
            std::ostringstream oss;
            oss << "completed mpqc::static_mutex<>::unlock(): count = " << lock_count << std::endl;
            sc::ExEnv::out0() << oss.str();
          }
        }
        static boost::mutex& get() {
            return mutex;
        }

        static boost::mutex mutex;
        static int64_t lock_count;
    };

    template<typename T>
    boost::mutex static_mutex<T>::mutex;
    template<typename T>
    int64_t static_mutex<T>::lock_count = int64_t(0);
    
    /// Static mutex instances.
    /// Example:
    /// @code
    /// mutex::global::lock();
    /// // critical code
    /// mutex::global::unlock()
    /// @endcode
    struct mutex {
        struct global_mutex_tag {};
        typedef static_mutex<global_mutex_tag> global;
    };

    /// @} mpqc::Utility

}

#endif // MPQC_UTILITY_MUTEX_HPP
