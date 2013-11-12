#ifndef MPQC_ARRAY_PARALLEL_HPP
#define MPQC_ARRAY_PARALLEL_HPP

#include "mpqc/range.hpp"
#include "mpqc/array/forward.hpp"
#include "mpqc/array/thread.hpp"
#include "mpqc/utility/string.hpp"

namespace mpqc {
namespace detail {

    template<class Data>
    struct array_tile {
	std::vector<range> extents;
	int proc, local;
	Data data;
	std::vector<range> subset(const std::vector<range> &ranges) const {
	    std::vector<range> s;
	    for (int i = 0; i < ranges.size(); ++i) {
		range r = (ranges[i] & this->extents[i]);
		if (!r.size()) return std::vector<range>();
		s.push_back(r);
	    }
	    return s;
	}
    };

    template<typename T, typename Driver>
    struct array_parallel_impl
	: ArrayBase, boost::noncopyable
    {

        typedef array_impl<T,Driver> Impl;
        typedef array_tile<Impl*> Tile;

	array_parallel_impl(const std::string &name,
			    const std::vector<size_t> &dims,
			    MPI::Comm comm)
            : ArrayBase(name, dims), thread_comm_(comm)
        {
            initialize(ArrayBase::dims_, comm);
        }

        void initialize(const std::vector<size_t> &dims, const MPI::Comm &comm) {

            int np = comm.size();
            size_t block = (dims.back() + np - 1)/np;

            for (int i = 0; i < dims.back(); i += block) {
                Tile tile;

                int n = std::min(dims.back()-i, block);

                foreach (size_t dim, dims) {
                    tile.extents.push_back(range(0, dim));
                }
                tile.extents.back() = range(i, i+n);

                tile.proc = (i/block)%np;
                tile.local = (tile.proc == comm.rank());

		tile.data = NULL;
                if (tile.local) {
                    std::string suffix = ".part" + string_cast(comm.rank());
		    try {
			tile.data = new Impl(this->name() + suffix, tile.extents);
		    }
		    catch (std::exception e) {
			comm.cout << e.what() << std::endl;
		    }
                }

                comm.broadcast(tile.data, tile.proc);
		if (!tile.data) {
		    throw std::runtime_error("failed to create parallel array segment");
		}
                tiles_.push_back(tile);
            }
	}

	~array_parallel_impl() {
	    foreach (Tile t, tiles_) {
	      if (t.local) delete t.data;
	    }
	    thread_comm_.comm().free();
	}

        void sync() {
            thread_comm_.sync();
        }

    protected:

        void _put(const std::vector<range> &r, const void *buffer) {
            _put(r, (const T*)buffer);
        }

        void _get(const std::vector<range> &r, void *buffer) const {
            _get(r, (T*)buffer);
        }

    private:

	void _put(const std::vector<range> &r, const T *buffer) {
            size_t total = size(r);
            size_t count = 0;
            foreach (const auto &tile, tiles_) {
                auto x = tile.subset(r);
                if (!x.empty()) {
                    // if (tile.local) {
                    //     tile.object->put(x, buffer);
                    // }
                    // // else {
                    // MPQC_PROFILE_LINE;
                    // printf("comm::write\n");
                    //if (tile.proc == 0)
                    thread_comm_.write(buffer + count, tile.data, x, tile.proc);
                    // //}
                    count += size(x);
                }
            }
            assert(total == count);
	}

	void _get(const std::vector<range> &r, T *buffer) const {
            size_t total = size(r);
            size_t count = 0;
            foreach (const auto &tile, tiles_) {
                auto x = tile.subset(r);
                if (!x.empty()) {
                    // if (tile.local) {
                    //     tile.object->get(x, buffer);
                    // }
                    // // else {
                    // MPQC_PROFILE_LINE;
                    // printf("comm::read\n");
                    //if (tile.proc == 0)
                    // std::cout << "read " << x << " from " << tile.proc << std::endl;
                    thread_comm_.read(buffer + count, tile.data, x, tile.proc);
                    count += size(x);
                    // //}
                }
            }
            assert(total == count);
	}

    private:

        static size_t size(const std::vector<range> &R) {
            size_t size = R.empty() ? 0 : 1;
            foreach (range r, R) {
                size *= r.size();
            }
            //printf("size = %lu\n", size);
            return size;
        }

	std::vector<Tile> tiles_;
        array_thread_comm thread_comm_;

    };


} // namespace detail
} // namespace mpqc

#endif /* MPQC_ARRAY_PARALLEL_HPP */
