#ifndef MPQC_RANGE_HPP
#define MPQC_RANGE_HPP

#include <iostream>
#include <vector>

#include <boost/range/irange.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/preprocessor/repetition.hpp>

#include <boost/fusion/include/size.hpp>
#include <boost/fusion/include/at_c.hpp>

/// @defgroup Range mpqc.Math.Range
/// Range objects for iterating or accessing slices of data

namespace mpqc {

/// @addtogroup Range
/// @{

    struct range_block_iterator;

    struct range : 
	boost::iterator_range<boost::range_detail::integer_iterator<int> >
    {
	typedef boost::iterator_range<
	    boost::range_detail::integer_iterator<int> > iterator_range;

        template<class S>
        struct tie;

	explicit range(int size = 0) : iterator_range(0, size) {}

	range(int begin, int end) : iterator_range(begin, end) {}

	int size() const { return end() - begin(); }

        bool operator==(const range &r) const {
            return (r.begin() == begin() && r.end() == end());
        }

        static range intersection(const range &a, const range &b) {
            int begin = std::max(*a.begin(), *b.begin());
            int end = std::min(*a.end(), *b.end());
            return (begin < end) ? range(begin, end) : range();
        }

        boost::iterator_range<range_block_iterator> block(size_t N) const;

        static std::vector<range> block(range r, size_t N) {
            std::vector<range> blocks;
            for (int i = *r.begin(); i < *r.end(); i += N) {
                blocks.push_back(range(i, std::min<int>(i+N, *r.end())));
            }
            return blocks;
        }

        std::vector<range> split2(size_t n) const {
            std::vector<range> v;
            int m = this->size()%n;
            int b = this->size()/n;
            int r = 0;
            for (int i = 0; i < m; ++i) {
                v.push_back(range(r, r+b+1));
                r += b+1;
            }
            for (int i = m; i < n; ++i) {
                v.push_back(range(r, r+b));
                r += b;
            }
            assert(r == this->size());
            return v;
        }

    };

    /// Range intersection
    inline range operator&(const range &a, const range &b) {
        return range::intersection(a,b);
    }

    /// print range as "begin:end"
    inline std::ostream& operator<<(std::ostream &os, const range &r) {
        os << *r.begin() << ":" << *r.end();
        return os;
    }

    /// print range vector as "[ begin:end, ... ]"
    inline std::ostream& operator<<(std::ostream &os, const std::vector<range> &r) {
        os << "[ ";
	for (int i = 0; i < r.size(); ++i) {
	    os << r[i] << ",";
	}
        os << " ]";
        return os;
    }

    inline range range0(const range &r) {
        return range(0, r.size());
    }

    struct range_block_iterator
        : boost::iterator_facade<range_block_iterator,
                                 range,
                                 boost::random_access_traversal_tag,
                                 range // reference
                                 >
    {
        range_block_iterator(range::iterator it, range r, int block)
            : data_(it), range_(r), block_(block) {}

    private:
        friend class boost::iterator_core_access;

        void increment() {
            advance(1);
        }

        bool equal(const range_block_iterator &it) const {
            assert(this->block_ == it.block_);
            assert(this->range_ == it.range_);
            return this->data_ == it.data_;
        }

        void advance(ptrdiff_t n) {
            data_ = std::min(data_ + n*block_, range_.end());;
        }

        ptrdiff_t distance_to(const range_block_iterator &it) const {
            assert(this->block_ == it.block_);
	    assert(this->range_ == it.range_);
	    ptrdiff_t n = (this->data_ - it.data_);
            return -(n + block_ - 1)/block_;
        }

        range dereference() const {
            range r(*data_, std::min<int>(*data_ + block_, *range_.end()));
            // std::cout << "*" << r << std::endl;
            return r;
        }

    private:
        range::iterator data_;
        range range_;
        int block_;
    };


    inline
    boost::iterator_range<range_block_iterator> range::block(size_t N) const {
        typedef range_block_iterator It;
        typedef boost::iterator_range<It> Range;
        return Range(It(this->begin(), *this, N), It(this->end(), *this, N));
    }

    /// Cast range to range, return argument unchanged
    inline range range_cast(const range &r) {
        return r;
    }

    /// Cast integral argument r to range [r:r+1)
    template<class R>
    typename boost::enable_if<boost::is_integral<R>, range>::type
    range_cast(const R &r) {
        //std::cout << "range_cast<R>" << r << std::endl;
        return range(r,r+1);
    }

    /// @} // Range group


}

namespace mpqc {

    template<class It>
    struct Range : boost::iterator_range<It>
    {
	typedef boost::iterator_range<It> iterator_range;
	Range(It begin, It end) : iterator_range(begin, end) {}
	//size_t size() const { return iterator_range::end() - begin(); }
    };

}

namespace sc {
    using mpqc::range;
}

#endif // MPQC_RANGE_HPP
