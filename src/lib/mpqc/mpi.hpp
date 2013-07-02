#ifndef MPQC_MPI_HPP
#define MPQC_MPI_HPP

#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX

#include "mpqc/config.h"

#ifndef HAVE_MPI
#error Missing MPI (HAVE_MPI not defined)
#endif

#include <mpi.h>

extern "C" {
#include <armci.h>
}

#include <stdio.h>
#include <stdarg.h>

#include <stdexcept>

#include "mpqc/utility/string.hpp"
#include "mpqc/utility/mutex.hpp"

namespace mpqc {
namespace mpi {

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

    template<typename T>
    MPI_Datatype type();
    // {
    // 	BOOST_STATIC_ASSERT_MSG(false, "MPI Type mapping is not implemented");
    // }

#define MPQC_MPI_TYPE(T, MPI_DATATYPE)                                  \
    template<> inline MPI_Datatype type<T>() { return MPI_DATATYPE; }

    MPQC_MPI_TYPE(int, MPI_INT)
    MPQC_MPI_TYPE(double, MPI_DOUBLE)


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

    struct Stream {
	explicit Stream(MPI_Comm comm) {
	    MPI_Comm_rank(comm, &rank_);
	}
	std::ostream& stream() const {
	    std::cout << rank_ << ": ";
	    return std::cout;
	}
    private:
	int rank_;
    };


    struct Comm {

	const Stream cout;

        explicit Comm(MPI_Comm comm)
            : cout(comm), comm_(comm) {}

        operator MPI_Comm() const {
            return comm_;
        }

        int rank() const {
            MPQC_MPI_THREADSAFE {
                int rank;
                MPI_Comm_rank(comm_, &rank);
                return rank;
            }
        }

        int size() const {
            MPQC_MPI_THREADSAFE {
                int size;
                MPI_Comm_size(comm_, &size);
                return size;
            }
        }

        void barrier() const {
            //printf("%i: Comm::barrier()\n", this->rank());
            MPI_Barrier(this->comm_);
        }

        MPI_Status recv(void *data, int count, MPI_Datatype type,
                        int src, int tag) const {
            MPQC_MPI_THREADSAFE {
                // printf("mpi::recv(data=%p, count=%i, src=%i, tag=%i)\n",
                // 	   data, count, src, tag);
                MPI_Status status;
                MPI_Recv(data, count, type, src, tag, comm_, &status);
                return status;
            }
        }

        template <typename T>
        T recv(int src, int tag) const {
            T data;
            mpi::Comm::recv(&data, sizeof(T), MPI_BYTE, src, tag);
            return data;
        }


        void send(const void *data, int count, MPI_Datatype type,
                  int dst, int tag) const {
            threadsafe lock;
            // printf("mpi::send(data=%p, count=%i, dst=%i, tag=%i)\n",
            // 	   data, count, dst, tag);
            MPI_Send((void*)data, count, type, dst, tag, comm_);
        }

        template <typename T>
        void send(const T &data, int dst, int tag) const {
            mpi::Comm::send((void*)&data, sizeof(T), MPI_BYTE, dst, tag);
        }


        void ssend(const void *data, int count, MPI_Datatype type,
                   int dst, int tag) const {
            threadsafe lock;
            // printf("mpi::ssend(data=%p, count=%i, dst=%i, tag=%i)\n",
            // 	   data, count, dst, tag);
            MPI_Ssend((void*)data, count, type, dst, tag, comm_);
        }

        template <typename T>
        void ssend(const T &data, int dst, int tag) const {
            mpi::Comm::ssend((void*)&data, sizeof(T), MPI_BYTE, dst, tag);
        }


        MPI_Request irecv(void *data, int count, MPI_Datatype type,
                          int src, int tag) const {
            threadsafe lock;
            //printf("Thread::recv(count=%i, src=%i, tag=%i)\n", count, src, tag);
            MPI_Request request;
            MPI_Irecv(data, count, type, src, tag, comm_, &request);
            return request;
        }

        MPI_Request isend(const void *data, int count, MPI_Datatype type,
                          int src, int tag) const {
            threadsafe lock;
            //printf("Thread::recv(count=%i, src=%i, tag=%i)\n", count, src, tag);
            MPI_Request request;
            MPI_Isend((void*)data, count, type, src, tag, comm_, &request);
            return request;
        }


        template <typename T>
        void broadcast(T &value, int root) const {
            MPI_Bcast(&value, sizeof(T), MPI_BYTE, root, comm_);
        }

        template <typename T>
        void broadcast(T *data, int count, int root) const {
            MPI_Bcast(data, sizeof(T)*count, MPI_BYTE, root, comm_);
        }

        bool any(const bool &value) const {
            int result = value;
            MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_LOR, comm_);
            return result;
        }

        template <typename T>
        void sum(T *value, int count) const {
            MPI_Allreduce(MPI_IN_PLACE, value, count, mpi::type<T>(), MPI_SUM, comm_);
        }

        template <typename T>
        void sum(T &value) const {
            sum(&value, 1);
        }

        template <typename T>
        std::vector<T> allgather(T value) const {
            int bytes = sizeof(T);
            std::vector<T> data(this->size(), T());
            MPI_Allgather(&value, bytes, MPI_BYTE, &data[0], bytes, MPI_BYTE, comm_);
            return data;
        }



        void printf(std::string fmt, ...) const {
            fmt = string_cast(this->rank()) + ": " + fmt;
            va_list args;
            va_start(args, fmt);
            vfprintf(stdout, fmt.c_str(), args);
            va_end(args);
        }    

    protected:

        MPI_Comm comm_;

    };


    template <typename T>
    std::ostream& operator<<(const Stream &s, const T &t) {
	s.stream() << t;
	return s.stream();
    }

    struct Task : boost::noncopyable {
        
        typedef int T;

        explicit Task(Comm comm)
            : comm_(MPI_COMM_WORLD)
        {
            assert(comm == MPI_COMM_WORLD);
            ARMCI_Init();
            data_.resize(comm_.size());
            ARMCI_Malloc(&data_[0], sizeof(T));
            reset(0);
        }

        ~Task() {
            ARMCI_Free(data_[comm_.rank()]);
        }

        void reset(const T &value = T(0)) {
            comm_.barrier();
            if (comm_.rank() == 0) {
                ARMCI_PutValueInt(value, this->value(), 0);
                ARMCI_Fence(0);
            }
            comm_.barrier();
        }

        T operator++(int) {
            int data;
            ARMCI_Rmw(ARMCI_FETCH_AND_ADD, &data, this->value(), 1, 0);
            return data;
        }

        range next(range r, int block = 1) {            
            int i = block*((*this)++);
            return (r & range(i, i+block));
        }

    private:

        std::vector< void* > data_;
        mpi::Comm comm_;

        T* value() {
            return (T*)data_[0];
        }

    };


}
}

#endif // MPQC_MPI_HPP
