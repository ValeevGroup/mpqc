#ifndef MPQC_ARRAY_HPP
#define MPQC_ARRAY_HPP

#include "mpqc/range.hpp"
#include "mpqc/range/operator.hpp"
#include "mpqc/mpi.hpp"

#include "mpqc/array/core.hpp"
#include "mpqc/array/file.hpp"
#include "mpqc/array/parallel.hpp"

#include "mpqc/utility/foreach.hpp"

namespace mpqc {

    template<typename T>
    struct Array {

        Array() {}

        template<typename Extent>
        Array(const std::string &name,
              const std::vector<Extent> &extents,
              mpi::Comm comm = mpi::Comm(MPI_COMM_SELF))
        {
	    initialize(name, extents, ARRAY_CORE, comm);
        }

        template<typename Extent, class Driver>
        Array(const std::string &name,
              const std::vector<Extent> &extents,
              Driver driver,
              mpi::Comm comm = mpi::Comm(MPI_COMM_SELF))
        {
	    initialize(name, extents, driver, comm);
        }

        Array(const Array &a)
        {
            range_ = a.range_;
            dims_ = a.dims_;
            impl_ = a.impl_;
        }

	void put(const T *buffer) {
	    impl_->put(this->range_, buffer);
	}

	void get(T *buffer) const {
	    impl_->get(this->range_, buffer);
        }

        void sync() {
            impl_->sync();
        }

        const std::vector<size_t>& dims() const {
            return dims_;
        }

        size_t rank() const {
            return dims().size();
        }

        operator vector<T>() const {
            vector<T> p(range_[0].size());
            *this >> p;
            return p;
        }

        operator matrix<T>() const {
            matrix<T> p(range_[0].size(), range_[1].size());
            *this >> p;
            return p;
        }

        Array operator()(const std::vector<range> &r) {
            return Array(*this, r);
        }

        Array operator()(const std::vector<range> &r) const {
            return Array(*this, r);
        }

        // generate range operators up to rank 4
        MPQC_RANGE_OPERATORS(4, Array, this->operator())
        MPQC_RANGE_CONST_OPERATORS(4, Array, this->operator())

        void write(File::Dataspace<T> f) const {
            assert(this->rank() == 2);
            size_t B = block();
            range ri = range_[0];
            vector<T> block(ri.size()*B);
            foreach (auto rj, range_[1].block(B)) {
                this->operator()(ri,rj) >> block;
                f(ri,rj) << block;
            }
        }

        void read(File::Dataspace<T> f) {
            assert(this->rank() == 2);
            size_t B = block();
            range ri = range_[0];
            vector<T> block(ri.size()*B);
            foreach (auto rj, range_[1].block(B)) {
                f(ri,rj) >> block;
		//printf("read %lu*%lu, data=%p\n", ri.size(), rj.size(), block.data());
                this->operator()(ri,rj) << block;
            }
        }

    protected:

	template<class Extent, class Driver>
	void initialize(const std::string &name,
			const std::vector<Extent> &extents,
			Driver driver, mpi::Comm comm) {
	    
            foreach (auto e, extents) {
                range_.push_back(detail::ArrayBase::extent(e));
                dims_.push_back(range_.back().size());
            }

	    detail::ArrayBase *impl = NULL;
	    if (comm == MPI_COMM_SELF) {
		impl = new detail::array_impl<T, Driver>(name, dims_);
	    }
	    else {
		impl = new detail::array_parallel_impl<T, Driver>(name, dims_, comm);
	    }
	    this->impl_.reset(impl);
	}

        Array(const Array &A, const std::vector<range> &r)
        {
            range_ = r;
            foreach (range r, range_) {
		//std::cout << r << ",";
                dims_.push_back(r.size());
            }
	    //std::cout << std::endl;
            impl_ = A.impl_;
        }

        size_t block() const { 
            assert(this->rank() == 2);
            range ri = range_[0];
            return std::max<size_t>((8 << 20)/(sizeof(T)*ri.size()), 1);
        }

    private:
        std::vector<range> range_;
        std::vector<size_t> dims_;
	std::shared_ptr<detail::ArrayBase> impl_;
    };


    template<typename T, class V>
    void operator<<(Array<T> A, const V &v) {
        A.put(v.data());
    }

    template<typename T, class V>
    void operator>>(Array<T> A, V &v) {
        A.get(v.data());
    }

    template<typename T, class E>
    void operator>>(Array<T> A, Eigen::Block<E> b) {
        //A >> b;
    }

}

#endif /* MPQC_ARRAY_HPP */
