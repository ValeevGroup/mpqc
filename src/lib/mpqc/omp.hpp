#ifndef MPQC_OMP_HPP
#define MPQC_OMP_HPP

namespace mpqc {
namespace omp {

    inline bool master() {
	int master = false;
#pragma omp master
	master = true;
	return master;
    }

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
    
}
}

#endif /* MPQC_OMP_HPP */
