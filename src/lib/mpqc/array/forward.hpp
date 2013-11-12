#ifndef MPQC_ARRAY_FORWARD_HPP
#define MPQC_ARRAY_FORWARD_HPP

#include "mpqc/range.hpp"
 #include "mpqc/utility/debug.hpp"

namespace mpqc {

    template<typename T>
    struct Array;

}

namespace mpqc {
namespace detail {

    template<typename T, typename Driver>
    struct array_impl;
    
    template<typename T, class Driver>
    struct array_parallel_impl;

    struct ArrayBase {

    public:

	template<typename U>
	static range extent(U n) {
	    return range(0,n);
	}

	static range extent(range r) {
	    return r;
	}

    protected:

	virtual void _put(const std::vector<range> &r, const void *buffer) = 0;
	virtual void _get(const std::vector<range> &r, void *buffer) const = 0;

    public:

	virtual ~ArrayBase() {}
        virtual void sync() = 0;

        template<typename Extent>
	ArrayBase(const std::string &name,
                  const std::vector<Extent> &extents)
        {
            name_ = name;
            for (size_t i = 0; i < extents.size(); ++i) {
                range r = extent(extents.at(i));
                base_.push_back(*r.begin());
                dims_.push_back(r.size());
            }
	}

        const std::string name() const {
            return this->name_;
        }

        size_t rank() const {
            assert(dims_.size() == base_.size());
            return dims_.size();
        }

        size_t size() const {
            size_t size = 1;
            for (size_t i = 0; i < dims_.size(); ++i) {
                size *= dims_.at(i);
            }
            return size;
        }

	void put(const std::vector<range> &r, const void *buffer) {
            _put(rebase(r), buffer);
	}

	void get(const std::vector<range> &r, void *buffer) const {
            _get(rebase(r), buffer);
	}

    protected:

        std::string name_;
        std::vector<size_t> dims_;
        std::vector<size_t> base_;

        void check_range(const std::vector<range> &R) const {
            assert(this->rank() == R.size());
            for (size_t i = 0; i < this->rank(); ++i) {
                range r = R.at(i);
                //std::cout << "range = " << r << "base = " << base_.at(i) << std::endl;
                assert(base_.at(i) <= *r.begin());
                assert(*r.end() <= base_.at(i) + dims_.at(i));
            }
        }

        std::vector<range> rebase(std::vector<range> R) const {
            check_range(R);
            for (size_t i = 0; i < R.size(); ++i) {
                range r = R.at(i);
                size_t base = base_.at(i);
                R.at(i) = range(*r.begin()-base, *r.end()-base);
            }
            return R;
        }

    };    


} // namespace detail
} // namespace mpqc


#endif /* MPQC_ARRAY_FORWARD_HPP */
