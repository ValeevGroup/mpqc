#ifndef MPQC_MPI_TASK_HPP
#define MPQC_MPI_TASK_HPP

#include "mpqc/mpi.hpp"
#include "mpqc/utility/mutex.hpp"

#ifdef MPQC_PARALLEL
extern "C" {
#include <armci.h>
}
#endif

namespace mpqc {
namespace MPI {

    /// Distributed task
    /// @ingroup MPI
    struct Task : boost::noncopyable {

        typedef int T;

        /// Construct new task
        /// @warning NOT threadsafe
        explicit Task(MPI::Comm comm)
            : comm_(comm), data_(0)
        {
#ifdef MPQC_PARALLEL
            assert(comm == MPI_COMM_WORLD);
            ARMCI_Init();
            data_.resize(comm_.size());
            ARMCI_Malloc(&data_[0], sizeof(T));
#endif
            reset(0);
        }

        /// Destructor
        /// @warning NOT threadsafe
        ~Task() {
#ifdef MPQC_PARALLEL
            ARMCI_Free(data_[comm_.rank()]);
#endif
        }

        /// Reset task
        void reset(const T &value = T(0)) {
            mutex::global::lock();
#ifdef MPQC_PARALLEL
            comm_.barrier();
            if (comm_.rank() == 0) {
                ARMCI_PutValueInt(value, this->value(), 0);
                ARMCI_Fence(0);
            }
            comm_.barrier();
#else
            data_ = 0;
#endif
            mutex::global::unlock();
        }

        /// Get next task
        T operator++(int) {
            int next;
            mutex::global::lock();
#ifdef MPQC_PARALLEL
            ARMCI_Rmw(ARMCI_FETCH_AND_ADD, &next, this->value(), 1, 0);
#else
            next = this->data_++;
#endif
            mutex::global::unlock();
            return next;
        }

        /// Get next task range
        mpqc::range next(range r, int block = 1) {            
            int i = block*((*this)++);
            return (r & mpqc::range(i, i+block));
        }

    private:

        MPI::Comm comm_;
#ifdef MPQC_PARALLEL
        std::vector< void* > data_;
        T* value() {
            return (T*)data_[0];
        }
#else
        T data_;
#endif

    };


} // namespace mpi
} // namespace mpqc

#endif /* MPQC_MPI_TASK_HPP */
