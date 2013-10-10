#ifndef MPQC_ARRAY_HPP
#define MPQC_ARRAY_HPP

#include "mpqc/range.hpp"
#include "mpqc/range/operator.hpp"
#include "mpqc/mpi.hpp"

#include "mpqc/array/core.hpp"
#include "mpqc/array/file.hpp"
#ifdef MPQC_PARALLEL
#include "mpqc/array/parallel.hpp"
#endif

#include "mpqc/assert.hpp"
#include "mpqc/utility/foreach.hpp"

namespace mpqc {

    /// @addtogroup MathArray
    /// @{


    /// Array implementation.
    /// Array can either be serial or parallel and reside in memory or filesystem.
    template<typename T>
    struct Array {

        /// Create new empty array
        Array() {}

        /// Create new array using CORE driver (in-memory)
        template<typename Extent>
        Array(const std::string &name,
              const std::vector<Extent> &extents,
              MPI::Comm comm = MPI::Comm::Self())
        {
	    initialize(name, extents, ARRAY_CORE, comm);
        }

        /// Create new array using an arbitrary driver
        template<typename Extent, class Driver>
        Array(const std::string &name,
              const std::vector<Extent> &extents,
              Driver driver,
              MPI::Comm comm = MPI::Comm(MPI::Comm::Self()))
        {
	    initialize(name, extents, driver, comm);
        }

        /// Shallow copy constructor
        /// @warning Shallow copy
        Array(const Array &a)
        {
            range_ = a.range_;
            dims_ = a.dims_;
            impl_ = a.impl_;
        }

        /// Put buffer data into array
        /// @warning buffer must be contigous
	void put(const T *buffer) {
	    impl_->put(this->range_, buffer);
	}

        /// Get array data into buffer
        /// @warning buffer must be contigous
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
            MPQC_CHECK(this->rank() == r.size());
            return Array(*this, r);
        }

        Array operator()(const std::vector<range> &r) const {
            MPQC_CHECK(this->rank() == r.size());
            return Array(*this, r);
        }


#ifdef DOXYGEN
        /**
           N-ary sub-array access operators.
           The parameters R should be either integral types (a single element)
           or of type mpqc::range (a range of elements)
           The method packs arguments into <c>std::vector<range></c>
           and calls the equivalent operator.
        */
        template<class R, ...>
        Array operator()(const R &r, ...);
#else
        template<class S>
        Array operator()(const range::tie<S> &t) {
            return this->operator()(std::vector<range>(t));
        }
        template<class S>
        Array operator()(const range::tie<S> &t) const {
            return this->operator()(std::vector<range>(t));
        }
        MPQC_RANGE_OPERATORS(4, Array, this->operator())
        MPQC_RANGE_CONST_OPERATORS(4, Array, this->operator())
#endif

        /// Write array data to File::Dataspace
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

        /// Read array data from File::Dataspace
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
			Driver driver, MPI::Comm comm) {
	    
            foreach (auto e, extents) {
                range_.push_back(detail::ArrayBase::extent(e));
                dims_.push_back(range_.back().size());
                //std::cout << "extent=" << range_.back() << std::endl;
            }

	    detail::ArrayBase *impl = NULL;
	    if (comm == MPI::Comm::Self()) {
		impl = new detail::array_impl<T, Driver>(name, dims_);
	    }
	    else {
#ifdef HAVE_MPI
		impl = new detail::array_parallel_impl<T, Driver>(name, dims_, comm);
#else // HAVE_MPI
                throw std::runtime_error("Parallel array not implemented: "
                                         "library was compiled without MPI");
#endif // HAVE_MPI
	    }
	    this->impl_.reset(impl);
	}

        Array(const Array &A, const std::vector<range> &r)
        {
            MPQC_CHECK(A.rank() == r.size());
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


    /**
       Write to Array from a generic vector V.
       @tparam V Vector with member <c>const T* V::data() const</c>
       @param A Array to write to
       @param v Vector to read from.
       @warning The pointer returned by V::data() must be contigous
    */
    template<typename T, class V>
    void operator<<(Array<T> A, const V &v) {
        A.put(v.data());
    }


    /**
       Read from Array to a generic vector V.
       @tparam V Vector with member <c>T* V::data()</c>
       @param A Array to read from
       @param v Vector to write to.
       @warning The pointer returned by V::data() must be contigous
    */
    template<typename T, class V>
    void operator>>(Array<T> A, V &v) {
        A.get(v.data());
    }

    // template<typename T, class E>
    // void operator>>(Array<T> A, Eigen::Block<E> b) {
    //     //A >> b;
    // }

    /// @} // group Array

}

#endif /* MPQC_ARRAY_HPP */
