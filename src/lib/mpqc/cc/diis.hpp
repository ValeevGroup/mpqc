#ifndef CC_DIIS_HPP
#define CC_DIIS_HPP

#include "cc/cc.hpp"
#include "exception.hpp"
#include "utility/array.hpp"

#include <cstdio>
#include <boost/numeric/ublas/io.hpp>
#include "boost/utility/profiler.hpp"

#include "blas.hpp"
#include "lapack.hpp"

namespace cchem {
namespace cc {

    struct DIIS {
    private:
	size_t no_, nv_;
	File::Group data_;
	bool enabled_;
    public:
	DIIS(size_t no, size_t nv, const File::Group &data, size_t max = 5)
	    : no_(no), nv_(nv), data_(data), enabled_(max > 0)
	{
	    if (!this->enabled_) return;
	    size_t dims[] = { no*no*nv, nv+1, max };
	    data_.create_dataset<double,3>("r", dims);
	    data_.create_dataset<double,3>("t", dims);
	    File::Dataset t = data_.dataset("t");
	    t.create_attribute<int>("max");
	    t.create_attribute<int>("size");
	    t.create_attribute<int>("iteration");
	    t.attribute("max") = int(max);
	    t.attribute("size") = int(-1);
	    t.attribute("iteration") = int(-1);
	}

	operator bool() const { return enabled_; }

	void initialize(Array &t1, Array &t2) {
	    using utility::make_array;

	    CCHEM_ASSERT(this->enabled_);

	    size_t no = this->no_;
	    size_t nv = this->nv_;
	    size_t no2v = no*no*nv;

	    File::Dataset t = data_.dataset("t");
	    File::Dataset r = data_.dataset("r");
	    int iter = -1;

	    Vector vj(no*no*nv);

	    // t2
	    for (size_t v = 0; v < nv; ++v) {
		t2.get(vj.data().begin(),
		       make_array<size_t>(0,0,0,v),
		       make_array<size_t>(no,no,nv,v+1));
		t.put(vj.data().begin(),
		      make_array<size_t>(0, v, iter+1),
		      make_array<size_t>(no2v, v+1, iter+2));
	    }
	    // t1
	    t1.get(vj.data().begin(),
		   make_array<size_t>(0, 0),
		   make_array<size_t>(no, nv));
	    t.put(vj.data().begin(),
		  make_array<size_t>(0, nv, iter+1),
		  make_array<size_t>(no*nv, nv+1, iter+2));

	    t.attribute("iteration") = int(iter+1);
	    t.attribute("size") = 0;

	}

	/**
	   Apply DIIS update to t1 and t2.
	 */
	double apply(Array &t1, Array &t2) {
	    using utility::make_array;

	    CCHEM_ASSERT(this->enabled_);

	    double rmax = 0;
	    size_t no = this->no_;
	    size_t nv = this->nv_;
	    size_t no2v = no*no*nv;

	    File::Dataset t = data_.dataset("t");
	    File::Dataset r = data_.dataset("r");
	    int iter = t.attribute("iteration");
	    int M = t.attribute("max"); // iter%M implement circular buffer
	    int N = std::min(iter+1, M)+1;

	    assert(iter > -1);

	    Vector vi(no*no*nv);
	    Vector vj(no*no*nv);

	    // t1 ranges
	    size_t start1[] = { 0, 0 };
	    size_t stop1[] = { no, nv };

	    Matrix xr(N,N,0);;
	    for (int i = 1; i < N; ++i) {
		xr(0,i) = -1;
		xr(i,0) = -1;
	    }

	    for (size_t v = 0; v < nv+1; ++v) {
		BOOST_PROFILE_LINE;
		size_t start[] = { 0, v, iter%M };
		size_t stop[] = { no2v, v+1, iter%M + 1 };
		bool v1 = (v == nv); // t1 iteration
		// r(i) = p(i+1) - p(i)
		if (v1) {
		    stop[0] = no*nv;
		    vj.clear();
		    vi.clear();
		    t1.get(vj.data().begin(), start1, stop1);
		}
		else { 
		    t2.get(vj.data().begin(),
			   make_array<size_t>(0,0,0,v),
			   make_array<size_t>(no,no,nv,v+1));
		}
		t.get(vi.data().begin(), start, stop);
		vj -= vi;
		rmax = std::max(rmax, norm_inf(vj));
		r.put(vj.data().begin(), start, stop);
		for (size_t j = N-1, J = iter; j > 0; --j, --J) {
		    start[2] = J%M;
		    stop[2] = start[2]+1;
		    r.get(vj.data().begin(), start, stop);
		    xr(j,j) += inner_prod(vj,vj);
		    for (size_t i = j-1, I = J-1; i > 0; --i, --I) {
			start[2] = I%M;
			stop[2] = start[2]+1;
			r.get(vi.data().begin(), start, stop);
			double d = inner_prod(vi,vj);
			xr(i,j) += d;
			xr(j,i) += d;
		    }
		}
	    }

	    // for (int i = 0; i < N; ++i) {
	    // 	std::cout << "A "<< row(xr,i) << std::endl;
	    // }

	    // e = (0,0,..,-1)
	    Vector e(N,0);
	    e(0) = -1;
	    Vector c = solve(xr,e);

	    // std::cout << "b* " << e << std::endl;
	    // std::cout << "c* " << c << std::endl;

	    for (size_t v = 0; v < nv+1; ++v) {
		BOOST_PROFILE_LINE;
		size_t start[] = { 0, v, 0 };
		size_t stop[] = { no2v, v+1, 0 };
		bool v1 = (v == nv); // t1 iteration
		if (v1) stop[0] = no*nv;
		vj.clear();
		// t(i+1) = c*t(i) + cr(i);
		for (int i = N-1, I = iter; i > 0; --i, --I) {
		    start[2] = I%M;
		    stop[2] = start[2]+1;
		    t.get(vi.data().begin(), start, stop);
		    vj += c(i)*vi;
		    r.get(vi.data().begin(), start, stop);
		    vj += c(i)*vi;
		}
		start[2] = (iter+1)%M;
		stop[2] = start[2]+1;
		t.put(vj.data().begin(), start, stop);
		if (v1) {
		    t1.put(vj.data().begin(), start1, stop1);
		}
		else {
		    t2.put(vj.data().begin(),
			   make_array<size_t>(0,0,0,v),
			   make_array<size_t>(no,no,nv,v+1));
		}
	    }
	    t.attribute("iteration") = int(iter+1);

	    return rmax;
	}

	template<class E>
    	static Vector solve(const E &ae, Vector b) {
	    BOOST_PROFILE_LINE;
	    Matrix A(ae);
	    boost::numeric::ublas::vector<lapack_int> ipiv(A.size1());
	    int status = lapack::gesv(A, ipiv, b);
	    if (status != 0) throw CCHEM_EXCEPTION("DIIS::solve failed");
	    return b;
	}
    };

} // namespace cc
} // namespace cchem

#endif // CC_DIIS_HPP
