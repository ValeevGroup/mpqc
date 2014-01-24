#ifndef MPQC_MPI_TASK_HPP
#define MPQC_MPI_TASK_HPP

#include "mpqc/mpi.hpp"
#include "mpqc/utility/mutex.hpp"

#ifdef HAVE_ARMCI
extern "C" {
#include <armci.h>
}
#endif

#if (defined HAVE_MPI) && !(defined HAVE_ARMCI)
#error mpqc::MPI::Task requires ARMCI if using MPI
#endif

namespace mpqc {
namespace MPI {

    /// Distributed task
    /// @ingroup CoreMPI
    struct Task : boost::noncopyable {

        typedef int T;

        /// Construct new task
        /// @warning NOT threadsafe
        explicit Task(const MPI::Comm &comm)
            : comm_(comm), data_(0)
        {
#ifdef HAVE_ARMCI
            MPQC_ASSERT(comm == MPI_COMM_WORLD);
            ARMCI_Init();
            data_.resize(comm_.size());
            ARMCI_Malloc(&data_[0], sizeof(T));
#endif
            reset(0);
        }

        /// Destructor
        /// @warning NOT threadsafe
        ~Task() {
#ifdef HAVE_ARMCI
            ARMCI_Free(data_[comm_.rank()]);
#endif
        }

        /// Reset task
        void reset(const T &value = T(0)) {
            mutex::global::lock();
#ifdef HAVE_ARMCI
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
#ifdef HAVE_ARMCI
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

        template<typename Iterator>
        Iterator next(Iterator begin, Iterator end) {
            Iterator it = begin;
            int n = (*this)++;
            for (int i = 0; i < n; ++i) {
                if (it == end) break;
                ++it;
            }
            return it;
        }

    private:

        MPI::Comm comm_;
#ifdef HAVE_ARMCI
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
