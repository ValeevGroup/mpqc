#ifndef MPQC_ARRAY_FILE_HPP
#define MPQC_ARRAY_FILE_HPP

#include "mpqc/range.hpp"
#include "mpqc/array/forward.hpp"
#include "mpqc/file.hpp"

namespace mpqc {
namespace detail {

    struct array_file_driver {
        array_file_driver() {}
    };

    template<typename T>
    struct array_impl<T, array_file_driver>
	: ArrayBase, boost::noncopyable
    {
        template<typename Extent>
	array_impl(const std::string &name,
                   const std::vector<Extent> &extents)
            : ArrayBase(name, extents)
        {
            mpqc::File file(name + ".h5");
            data_.reset(new mpqc::File::Dataset<T>(file, "data", this->dims_));
	}
        void sync() {
            //data_->flush();
        }
	void _put(const std::vector<range> &r, const void *buffer) {
            //std::cout << "write: " << r << std::endl;
            (*data_)(r).write((const T*)buffer);
	}
	void _get(const std::vector<range> &r, void *buffer) const {
            //std::cout << "read: " << r << std::endl;
            (*data_)(r).read((T*)buffer);
	}
    private:
        std::auto_ptr< mpqc::File::Dataset<T> > data_;
    };

#ifdef H5_HAVE_PARALLEL
#if 0 // TODO Andrey, do we need this?
    template<typename T>
    struct array_parallel_impl<T, array_file_driver>
	: ArrayBase, boost::noncopyable
    {
        array_parallel_impl(const std::string &name,
		     const std::vector<size_t> &dims,
		     mpi::Comm comm)
            : ArrayBase(dims)
        {
	    File file = (comm == MPI_COMM_NULL) ?
		File(name + ".h5") :
		File(name + ".h5", File::Driver::MPIIO(comm));
	    dataset_ = new File::Dataset<T>(file, "data", this->dims_);
	}
	void put(const std::vector<range> &r, const void *buffer) {
	    dataset_(r).put((const void*)buffer);
	}
	void get(const std::vector<range> &r, void *buffer) const {
	    dataset_(r).get((void*)buffer);
	}
    private:
	File::Dataset<T> dataset_;
    };
#endif
#endif // H5_HAVE_PARALLEL


} // namespace detail
} // namespace mpqc

namespace mpqc {

    static const detail::array_file_driver ARRAY_FILE;

}

#endif /* MPQC_ARRAY_FILE_HPP */
