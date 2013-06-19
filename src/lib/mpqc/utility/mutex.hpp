#ifndef MPQC_MUTEX_HPP
#define MPQC_MUTEX_HPP

#include <boost/thread/mutex.hpp>

namespace mpqc {

    template<typename>
    struct static_mutex {
        static void lock() { get().lock(); }
        static void unlock() { get().unlock(); }
        static boost::mutex& get() {
            static boost::mutex mutex;
            return mutex;
        }
    };
    
    struct mutex {
        struct global_mutex_tag;
        typedef static_mutex<global_mutex_tag> global;
    };

}

#endif // MPQC_MUTEX_HPP
