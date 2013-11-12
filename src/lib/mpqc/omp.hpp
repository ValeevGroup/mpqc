#ifndef MPQC_OMP_HPP
#define MPQC_OMP_HPP

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

namespace mpqc {
namespace omp {

    inline bool master() {
	int master = false;
#pragma omp master
	master = true;
	return master;
    }

    inline int max_threads() {
#ifndef _OPENMP
        return 1;
#else
        return omp_get_max_threads();
#endif
    };

    template <typename T>
    struct task : boost::noncopyable {
        task() : value_() {}
        T operator++(int) {
            T v;
#pragma omp critical(mpqc_omp_task)
            v = value_++;
            return v;
        }
    private:
        T value_;
    };

    struct mutex : boost::noncopyable {
#ifndef _OPENMP
        mutex() {}
        void lock() {}
        void unlock() {}
#else // _OPENMP
        mutex() {
            omp_init_lock(&lock_);
        }
        void lock() {
            omp_set_lock(&lock_);
        }
        void unlock() {
            omp_unset_lock(&lock_);
        }
    private:
        omp_lock_t lock_;
#endif // _OPENMP
    };
    
}
}

#endif /* MPQC_OMP_HPP */
