#ifndef MPQC_ARRAY_CORE_HPP
#define MPQC_ARRAY_CORE_HPP

#include "mpqc/array/forward.hpp"

#include "mpqc_config.h"
#include "mpqc/mpi.hpp"
#ifdef HAVE_MPI
#include "mpqc/array/parallel.hpp"
#ifdef HAVE_ARMCI
extern "C" {
#include <armci.h>
}
#endif
#endif // HAVE_MPI

#include "mpqc/range.hpp"
#include "mpqc/utility/foreach.hpp"
#include "mpqc/utility/mutex.hpp"
#include "mpqc/utility/exception.hpp"

#include <boost/noncopyable.hpp>

namespace mpqc {
namespace detail {

    struct array_core_driver {
        array_core_driver() {}
    };

    template<typename T>
    struct array_impl<T, array_core_driver>
	: ArrayBase, boost::noncopyable
    {

        template<typename Extent>
	array_impl(const std::string &name,
                   const std::vector<Extent> &extents,
                   const MPI::Comm &comm)
            : ArrayBase(name, extents, comm)
        {
            if (!(ArrayBase::comm_ == MPI::Comm::Self()))
                throw MPQC_EXCEPTION("Serial Array implementation must use MPI_COMM_SELF");
            data_.resize(this->size());
	}

        void sync() {}

    protected:

	void _put(const std::vector<range> &r, const void *buffer) {
            apply(putv, r, this->dims_, &this->data_[0], (const T*)buffer);
	}

	void _get(const std::vector<range> &r, void *buffer) const {
            apply(getv, r, this->dims_, &this->data_[0], (T*)buffer);
	}

    private:

        std::vector<T> data_;

        static size_t putv(size_t size, T *data, const T *buffer) {
            std::copy(buffer, buffer+size, data);
            return size;
        }

        static size_t getv(size_t size, const T *data, T *buffer) {
            std::copy(data, data+size, buffer);
            return size;
        }

        template<class F, typename data_ptr, typename buffer_ptr>
        static size_t apply(F f,
                            const std::vector<range> &r,
                            const std::vector<size_t> &dims,
                            data_ptr data, buffer_ptr buffer) {
            size_t begin = 0;
            std::vector<size_t> counts;
            std::vector<size_t> strides;
            counts.push_back(1);
            strides.push_back(1);
            for (size_t i = 0, stride = 1; i < dims.size(); ++i) {
                counts.back() *= r[i].size();
                begin += *r[i].begin()*stride;
                stride *= dims[i];
                // std::cout << "counts = " << counts.back() << std::endl;
                // std::cout << "strides = " << strides.back() << std::endl;
                if (dims[i] != r[i].size()) {
                    counts.push_back(1);
                    strides.push_back(stride);
                }
            }
            return apply(f, counts, strides, data+begin, buffer, strides.size());
        }

        template<class F, typename data_ptr, typename buffer_ptr>
        static size_t apply(F f,
                            const std::vector<size_t> &counts,
                            const std::vector<size_t> &strides,
                            data_ptr data, buffer_ptr buffer,
                            int level) {
            if (level == 0) 
                return 0;
            if (level == 1) {
                //std::cout << "count1 = " << counts[0] << std::endl;
                return f(counts[0], data, buffer);
            }
            size_t size = 0;
            --level;
            for (size_t i = 0; i < counts[level]; ++i) {
                //std::cout << size << " " << strides[level] << std::endl;
                size_t n = apply(f, counts, strides, data, buffer+size, level);
                data += strides[level];
                size += n;
            }
            return size;
        }


    };

#ifdef HAVE_ARMCI

    template<typename T>
    struct array_parallel_impl<T, array_core_driver>
	: ArrayBase, boost::noncopyable
    {

	array_parallel_impl(const std::string &name,
			    const std::vector<size_t> &dims,
			    const MPI::Comm &comm)
            : ArrayBase(name, dims, comm)
        {
            initialize(ArrayBase::dims_, ArrayBase::comm_);
        }

    private:

        void initialize(const std::vector<size_t> &dims, const MPI::Comm &comm) {

            // dont have code to use ARMCI groups yet
	    assert(ArrayBase::comm_ == MPI::Comm::World());

	    size_t size = 1;
	    std::vector<range> extents;
	    for (int i = 0; i < int(dims.size())-1; ++i) {
		size *= dims[i];
		extents.push_back(range(0, dims[i]));
	    }
	    size_t block = (dims.back() + comm.size() - 1)/comm.size();
	    size *= block;

	    check(ARMCI_Init(), "ARMCI_Init", comm);
	    
	    std::vector<void*> data(comm.size(), NULL);
	    check(ARMCI_Malloc(&data[0], size*sizeof(T)), "ARMCI_Malloc", comm);
            data_ = data[comm.rank()];

	    for (size_t i = 0; i < comm.size(); ++i) {
		Tile tile;
		tile.data = data[i];
		tile.proc = i;
		tile.local = (i == comm.rank());
		size_t begin = std::min<size_t>(dims.back(), i*block);
		size_t end = std::min<size_t>(dims.back(), begin+block);
		tile.extents = extents;
		tile.extents.push_back(range(begin, end));
		tiles_.push_back(tile);
		//printf("tile[%i].ptr=%p\n", i, tile.ptr);
	    }

	}

	~array_parallel_impl() {
            ARMCI_Free(this->data_);
	}

	void sync() {
	    ARMCI_Barrier();
	}

        const MPI::Comm& comm() const {
            return this->comm_;
        }

    protected:

	void _put(const std::vector<range> &r, const void *buffer) {
	    apply<PUT>(this->tiles_, r, (void*)buffer);
	}

	void _get(const std::vector<range> &r, void *buffer) const {
	    apply<GET>(this->tiles_, r, buffer);
	}

    private:

	enum OP { PUT, GET };

        typedef array_tile<void*> Tile;
        void *data_;
	std::vector<Tile> tiles_;

	static void check(int err, const std::string &func,
			  MPI::Comm comm = MPI::Comm::Self()) {
	    if (comm.any((err != 0))) {
		throw std::runtime_error(func + " failed");
	    }
	}

	template<OP op>
	static void apply(const std::vector<Tile> &tiles, 
			  const std::vector<range> &r,
			  void *buffer) {
	    
	    T *local = (T*)buffer;

	    foreach (Tile t, tiles) {

		std::vector<range> x = t.subset(r);
		//std::cout << t.extents << " ^ " << r << " = " << x << std::endl;
		if (x.empty()) continue;

		size_t total = 1;
		size_t N = x.size();
		int remote_strides[N];
		int local_strides[N];
		int count[N];

		for (size_t i = 0; i < N; ++i) {
		    total *= x[i].size();
		    count[i] = x[i].size();
		    local_strides[i] =
			x[i].size()*((i > 0) ? local_strides[i-1] : 1);
		    remote_strides[i] =
			t.extents[i].size()*((i > 0) ? remote_strides[i - 1] : 1);
		    local_strides[i] *= sizeof(T);
		    remote_strides[i] *= sizeof(T);
		}
                count[0] *= sizeof(T);

                T *remote = (T*)t.data;
                for (size_t i = 0, stride = 1; i < N; ++i) {
                    remote += (x[i].begin() - t.extents[i].begin())*stride;
                    stride *= t.extents[i].size();
                }

		// for (size_t i = 0; i < N; ++i) {
		//     printf("%i count=%i, remote_strides=%i, local_strides=%i, total=%lu\n",
		// 	   i, count[i], remote_strides[i], local_strides[i], total);
		// }
		// printf("tile.ptr=%p\n", t.ptr);

                mutex::global::lock();
		if (op == GET) {
		    int err;
		    err = ARMCI_GetS(remote, remote_strides,
				     local, local_strides,
				     count, N-1, t.proc);
		    check(err, "ARMCI_GetS failed");
		}
		if (op == PUT) {
		    int err;
		    err = ARMCI_PutS(local, local_strides,
				     remote, remote_strides,
				     count, N-1, t.proc);
		    check(err, "ARMCI_PutS failed");
		}
                mutex::global::unlock();

		local += total;

	    }
	}

    };
#endif // HAVE_ARMCI

} // namespace detail
} // namespace mpqc


namespace mpqc {
    static const detail::array_core_driver ARRAY_CORE;
}


#endif /* MPQC_ARRAY_CORE_HPP */
