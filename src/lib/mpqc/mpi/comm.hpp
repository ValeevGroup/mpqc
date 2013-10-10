#ifndef MPQC_MPI_COMM_HPP
#define MPQC_MPI_COMM_HPP

#include "mpqc/mpi.hpp"
#include "mpqc/utility/string.hpp"

#include <stdarg.h>

namespace mpqc {
namespace MPI {

    /// MPI_Comm object wrapper/stub
    /// @ingroup CoreMPI
    struct Comm {

    protected:
#ifdef HAVE_MPI
        // warning: comm_ must be initialized prior to cout
        MPI_Comm comm_;
#endif

    public:
        struct OStream {
            explicit OStream(const MPI::Comm &comm) {
                this->rank_ = comm.rank();
            }
            std::ostream& stream() const {
                std::cout << rank_ << ": ";
                return std::cout;
            }
        private:
            int rank_;
        };

	const OStream cout;

        void printf(std::string fmt, ...) const {
            fmt = string_cast(this->rank()) + ": " + fmt;
            va_list args;
            va_start(args, fmt);
            vfprintf(stdout, fmt.c_str(), args);
            va_end(args);
        }

        // serial "MPI" stubs, PTP aren't available in serial
#ifndef HAVE_MPI
    public:
        static MPI::Comm Self() { return Comm(); }
        static MPI::Comm World() { return Comm(); }
        static MPI::Comm dup(MPI::Comm comm) { return Comm(); }
        bool operator==(const Comm &comm) const {
            return true;
        }
        void free() {}
        int rank() const { return 0; }
        int size() const { return 1; }
        void barrier() const {}
    private:
        Comm() : cout(*this) {}
#endif

        // following functions are used if MPI is available
#ifdef HAVE_MPI

    public:

        static MPI::Comm Self() {
            return Comm(MPI_COMM_SELF);
        }

        static MPI::Comm World() {
            return Comm(MPI_COMM_WORLD);
        }

        /// Duplicate communicator
        static MPI::Comm dup(MPI::Comm comm) {
            MPI_Comm dup;
            MPI_Comm_dup(comm, &dup);
            return Comm(dup);
        }

        /// Free communicator
        /// @warning currently no mechanism to ensure shared comms are freed properly
        void free() {
            MPI_Comm_free(&this->comm_);
        }

        explicit Comm(MPI_Comm comm)
            : comm_(comm), cout(*this) {}

        operator MPI_Comm() const {
            return comm_;
        }

        bool operator==(const Comm &comm) const {
            return this->comm_ == comm.comm_;
        }

        int rank() const {
            MPQC_MPI_THREADSAFE;
            int rank;
            MPI_Comm_rank(comm_, &rank);
            return rank;
        }

        int size() const {
            MPQC_MPI_THREADSAFE;
            int size;
            MPI_Comm_size(comm_, &size);
            return size;
        }

        void barrier() const {
            //printf("%i: Comm::barrier()\n", this->rank());
            MPI_Barrier(this->comm_);
        }

        MPI_Status recv(void *data, int count, MPI_Datatype type,
                        int src, int tag) const {
            MPQC_MPI_THREADSAFE;
            // printf("MPI::recv(data=%p, count=%i, src=%i, tag=%i)\n",
            // 	   data, count, src, tag);
            MPI_Status status;
            MPI_Recv(data, count, type, src, tag, comm_, &status);
            return status;
        }

        template <typename T>
        T recv(int src, int tag) const {
            T data;
            MPI::Comm::recv(&data, sizeof(T), MPI_BYTE, src, tag);
            return data;
        }

        void send(const void *data, int count, MPI_Datatype type,
                  int dst, int tag) const {
            MPQC_MPI_THREADSAFE;
            // printf("MPI::send(data=%p, count=%i, dst=%i, tag=%i)\n",
            // 	   data, count, dst, tag);
            MPI_Send((void*)data, count, type, dst, tag, comm_);
        }

        template <typename T>
        void send(const T &data, int dst, int tag) const {
            MPI::Comm::send((void*)&data, sizeof(T), MPI_BYTE, dst, tag);
        }


        void ssend(const void *data, int count, MPI_Datatype type,
                   int dst, int tag) const {
            MPQC_MPI_THREADSAFE;
            // printf("MPI::ssend(data=%p, count=%i, dst=%i, tag=%i)\n",
            // 	   data, count, dst, tag);
            MPI_Ssend((void*)data, count, type, dst, tag, comm_);
        }

        template <typename T>
        void ssend(const T &data, int dst, int tag) const {
            MPI::Comm::ssend((void*)&data, sizeof(T), MPI_BYTE, dst, tag);
        }

        MPI_Request irecv(void *data, int count, MPI_Datatype type,
                          int src, int tag) const {
            MPQC_MPI_THREADSAFE;
            //printf("Thread::recv(count=%i, src=%i, tag=%i)\n", count, src, tag);
            MPI_Request request;
            MPI_Irecv(data, count, type, src, tag, comm_, &request);
            return request;
        }

        MPI_Request isend(const void *data, int count, MPI_Datatype type,
                          int src, int tag) const {
            MPQC_MPI_THREADSAFE;
            //printf("Thread::recv(count=%i, src=%i, tag=%i)\n", count, src, tag);
            MPI_Request request;
            MPI_Isend((void*)data, count, type, src, tag, comm_, &request);
            return request;
        }

#endif // HAVE_MPI

        // Collective calls, either MPI or serial stubs

    public:

        template <typename T>
        void broadcast(T &value, int root) const {
#ifdef HAVE_MPI
            MPI_Bcast(&value, sizeof(T), MPI_BYTE, root, comm_);
#endif
        }

        template <typename T>
        void broadcast(T *data, int count, int root) const {
#ifdef HAVE_MPI
            MPI_Bcast(data, sizeof(T)*count, MPI_BYTE, root, comm_);
#endif
        }

        bool any(const bool &value) const {
            int result = value;
#ifdef HAVE_MPI
            MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_LOR, comm_);
#endif
            return result;
        }

        template <typename T>
        void sum(T *value, int count) const {
#ifdef HAVE_MPI
            MPI_Allreduce(MPI_IN_PLACE, value, count, MPI::type<T>(), MPI_SUM, comm_);
#endif
        }

        template <typename T>
        void sum(T &value) const {
            sum(&value, 1);
        }

        template <typename T>
        std::vector<T> allgather(T value) const {
            int bytes = sizeof(T);
            std::vector<T> data(this->size(), T());
#ifdef HAVE_MPI
            MPI_Allgather(&value, bytes, MPI_BYTE, &data[0], bytes, MPI_BYTE, comm_);
#else
            data[0] = value;
#endif
            return data;
        }

    };


    template <typename T>
    std::ostream& operator<<(const MPI::Comm::OStream &s, const T &t) {
	s.stream() << t;
	return s.stream();
    }


}
}

#endif // MPQC_MPI_COMM_HPP
