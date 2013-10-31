#ifndef MPQC_RANGE_HPP
#define MPQC_RANGE_HPP

#include <iostream>
#include <vector>
#include <stdint.h>

#include <boost/range/irange.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/preprocessor/repetition.hpp>

#include <boost/fusion/include/size.hpp>
#include <boost/fusion/include/at_c.hpp>

namespace mpqc {

    /// @addtogroup MathRange
    /// @{

    struct range : 
	boost::iterator_range<boost::range_detail::integer_iterator<int64_t> >
    {

        // typedef boost::iterator_range<
        //     boost::range_detail::integer_iterator<int>
        //     >::iterator_category iterator_category;

	typedef boost::iterator_range<
	    boost::range_detail::integer_iterator<int64_t> > iterator_range;

        template<class S>
        struct tie;

	explicit range(int64_t size = 0)
            : iterator_range(int64_t(0), size) {}

	range(int64_t begin, int64_t end)
            : iterator_range(begin, end) {}

	int64_t size() const {
            return end() - begin();
        }

        bool operator==(const range &r) const {
            return (r.begin() == begin() && r.end() == end());
        }

        bool test(int64_t value) const {
            return ((*this->begin() <= value) && (value < *this->end()));
        }

        static range intersection(const range &a, const range &b) {
            int64_t begin = std::max(*a.begin(), *b.begin());
            int64_t end = std::min(*a.end(), *b.end());
            return (begin < end) ? range(begin, end) : range();
        }

        static std::vector<range> split(range r, size_t N) {
            if (N >= r.size())
                return std::vector<range>(1, r);
            std::vector<range> blocks;
            for (auto i = *r.begin(); i < *r.end(); i += N) {
                blocks.push_back(range(i, std::min<int64_t>(i+N, *r.end())));
            }
            return blocks;
        }

        static std::vector<range> split2(range R, size_t N) {
            std::vector<range> v;
            int64_t m = R.size()%N;
            int64_t b = R.size()/N;
            int64_t r = 0;
            for (int64_t i = 0; i < m; ++i) {
                v.push_back(range(r, r+b+1));
                r += b+1;
            }
            for (int64_t i = m; i < N; ++i) {
                v.push_back(range(r, r+b));
                r += b;
            }
            assert(r == R.size());
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
