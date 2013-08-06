#ifndef MPQC_UTILITY_MUTEX_HPP
#define MPQC_UTILITY_MUTEX_HPP

#include <boost/thread/mutex.hpp>

namespace mpqc {

    /// @addtogroup Utility
    /// @{

    /// Static mutex factory
    template<typename>
    struct static_mutex {
        static void lock() { get().lock(); }
        static void unlock() { get().unlock(); }
        static boost::mutex& get() {
            return mutex;
        }
	static boost::mutex mutex;
    };

    template<typename T>
    boost::mutex static_mutex<T>::mutex;
    
    /// Static mutex instance container
    /// Use as mutex::global::lock()/mutex::global::unlock()
    struct mutex {
        struct global_mutex_tag;
        typedef static_mutex<global_mutex_tag> global;
    };

    /// @} mpqc::Utility

}

#endif // MPQC_UTILITY_MUTEX_HPP
