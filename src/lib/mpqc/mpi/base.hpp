#ifndef MPQC_MPI_BASE_HPP
#define MPQC_MPI_BASE_HPP

#include "mpqc_config.h"

#ifdef HAVE_MPI
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif

#include "mpqc/utility/mutex.hpp"

#include <stdio.h>
#include <stdexcept>

namespace mpqc {
namespace MPI {

/// @addtogroup CoreMPI
/// @{

    void initialize(int thread_level);
    void finalize();
    std::string get_processor_name();

    // MPI stubs
#ifndef HAVE_MPI

    inline void initialize(int thread_level) {}
    inline void finalize() {}
    inline std::string get_processor_name() { return "localhost"; }

#endif // HAVE_MPI
    
#ifdef HAVE_MPI

#ifndef DOXYGEN
    /// RAII mutex::global lock
    /// Will only lock if MPI is not MPI_THREAD_MULTIPLE
    struct threadsafe : boost::noncopyable {
        threadsafe() {
            locked_ = false;
            int thread = MPI_THREAD_SINGLE;
            MPI_Query_thread(&thread);
            if (thread != MPI_THREAD_MULTIPLE) {
                mutex::global::lock();
                locked_ = true;
            }
        }
        ~threadsafe() {
            if (locked_) mutex::global::unlock();
        }
    private:
        bool locked_;
    };
#define MPQC_MPI_THREADSAFE threadsafe _threadsafe;
#endif // DOXYGEN

    template<typename T>
    MPI_Datatype type();
    // {
    // 	BOOST_STATIC_ASSERT_MSG(false, "MPI Type mapping is not implemented");
    // }

#define MPQC_MPI_TYPE(T, MPI_DATATYPE)                                  \
    template<> inline MPI_Datatype type<T>() { return MPI_DATATYPE; }

    MPQC_MPI_TYPE(int, MPI_INT)
    MPQC_MPI_TYPE(double, MPI_DOUBLE)

#undef MPQC_MPI_TYPE

    inline void initialize(int thread_level) {
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (!initialized) {
            int argc = 0;
            char **argv = NULL;
            int thread;
            MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread);
        }
        int thread = MPI_THREAD_SINGLE;
        MPI_Query_thread(&thread);
        if (thread < MPI_THREAD_SERIALIZED) {
	    throw std::runtime_error("thread < MPI_THREAD_SERIALIZED");
	}
    }

    inline void finalize() {
	MPI_Finalize();
    }

    inline std::string get_processor_name() {
        char name[MPI_MAX_PROCESSOR_NAME];
        int len;
        MPI_Get_processor_name(name, &len);
        return std::string(name, len);
    }

    inline MPI_Status wait(MPI_Request request, size_t us = 0) {
        MPQC_MPI_THREADSAFE {
            MPI_Status status;
            if (!us) {
                MPI_Wait(&request, &status);
                return status;
            }
            while (true) {
                int flag = 0;
                MPI_Test(&request, &flag, &status);
                if (flag) return status;
                boost::this_thread::sleep(boost::posix_time::microseconds(us));
            }
        }
    }

#endif // HAVE_MPI

    /// @}

}
}

#endif // MPQC_MPI_BASE_HPP
