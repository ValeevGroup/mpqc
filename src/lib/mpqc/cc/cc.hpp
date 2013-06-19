#ifndef MPQC_CC_CC_HPP
#define MPQC_CC_CC_HPP

// #include "runtime.hpp"
// #include "core/wavefunction.hpp"
// #include "integrals/screening.hpp"

#include "mpqc/array.hpp"

#include <cstdlib>
#include <map>
#include <stdexcept>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/ref.hpp>

namespace cchem {
namespace cc {

    template<class P>
    struct Map {
	struct no_such_key : std::runtime_error {
	    no_such_key(const char *what) : std::runtime_error(what) {}
	};

	typedef typename std::map<std::string, P>::iterator iterator;
	typedef typename std::map<std::string, P>::const_iterator const_iterator;

	Map() {}
	template<class U>
	Map(const Map<U> &other)
	    : data_(other.data().begin(), other.data().end()) {}

	void clear() { data_.clear(); }

	iterator begin() { return data_.begin(); }
	iterator end() { return data_.end(); }

	const_iterator begin() const { return data_.begin(); }
	const_iterator end() const { return data_.end(); }

	P& operator[](const char *key) {
	    // std::cout << "[" << key << "]" << std::endl;
	    return data_[key];
	}
	const P& operator[](const char *key) const {
	    // std::cout << "(" << key << ")" << std::endl;
	    if (!has(key)) throw no_such_key(key);
	    return data_.find(key)->second;
	}
	bool has(const char *key) const {
	    return (data_.find(key) != data_.end());
	}
	const std::map<std::string, P>& data() const {
	    return data_;
	}
    private:
	std::map<std::string, P> data_;
    };

    template<class P>
    std::ostream& operator<<(std::ostream &os, const Map<P> &m) {
	BOOST_AUTO(it, m.data().begin());
	while (it != m.data().end()) {
	    os << "[" << it->first << "]";
	    ++it;
	}
	return os;
    }

    typedef boost::numeric::ublas::vector<double> Vector;
    typedef boost::numeric::ublas::column_major Layout;
    typedef boost::numeric::ublas::matrix<double, Layout> Matrix;

    struct DIIS;

    struct Denominator {
	template<class E>
	Denominator(const boost::numeric::ublas::matrix_expression<E> &f) {
	    const E &e = f();
	    assert((e.size1() == e.size2()));
	    e_.resize(e.size1());
	    for (size_t i = 0; i < e.size1(); ++i) {
		e_(i) = e(i,i);
	    }
	}
	double operator()(int i) const {
	    return e_(i);
	}
	double operator()(int i, int j) const {
	    return (e_(i) + e_(j));
	}
	const Vector& data() const { return e_; }
    private:
	Vector e_;
    };


    struct Energy {
	double& operator[](const std::string &key) {
	    //assert(data_.count(key));
	    //return data_.find(key)->second;
	    return data_[key];
	}
	const double& operator[](const std::string &key) const {
	    assert(data_.count(key));
	    return data_.find(key)->second;
	}
    private:
	std::map<std::string, double> data_;
    };

    typedef ::Array<double> Array;

    Energy energy(Wavefunction wf, Runtime &rt, const std::string &method);


    void integrals(const Parallel &pe, const Wavefunction &wf,
		   const Map<Array*> &v,
		   const ::integrals::Screening &screening);

    double mp2(const Parallel &pe, const Wavefunction &wf,
	       const Array &V, Array &t);

    struct sd {
	bool ccd;
	double cutoff;
	double convergence;
	size_t max_iter;
	int iteration_;

	sd(const Runtime &rt, bool ccd = false) :
	    ccd(ccd),
	    cutoff(rt.get<double>("/cc/integrals/cutoff", 1e-10)), 
	    convergence(rt.get<double>("/cc/convergence", 1e-10)),
	    max_iter(rt.get<int>("/cc/max_iter", 10))
	{
	    iteration_ = 0;
	}

	double energy(const Parallel &pe,
		      const Wavefunction &wf,
		      const Map<const Array*> &V,
		      const Map<Array*> &A,
		      DIIS *diis);

	static void direct(const Parallel &pe,
			   const Wavefunction &wf,
			   const Matrix &t1, const Array &t2,
			   Map<Matrix> &O1, const Map<Array*> &O2,
			   const integrals::Screening &screening);

	void energy1(const Parallel &pe,
		     const Wavefunction &wf,
		     const Map<const Array*> &V,
		     const Map<Array*> &A,
		     double *E, double *rmax);

	double energy1(const Parallel &pe,
		       const Wavefunction &wf,
		       const Array &V,
		       const Matrix &t1, const Array &t2);

    };

    namespace triples {
	Energy energy(const Parallel &pe,
		      const Wavefunction &wf,
		      const Map<const Array*> &V,
		      const Map<const Array*> &t);
    }

    struct Correction {
	double ots, otd, ets, etd;
	Correction() : ots(0), otd(0), ets(0), etd(0) {}
	Correction& operator+=(const Correction &c) {
	    ots += c.ots;
	    otd += c.otd;
	    ets += c.ets;
	    etd += c.etd;
	    return *this;
	}
    };

    struct Triples {

	typedef cc::Correction Correction;

	typedef std::map<
	    std::string,
	    boost::reference_wrapper<const Array>
	    > Data;


	Correction operator()(size_t no, size_t nv,
			      Map<const Array*> t,
			      Map<const Array*> V,
			      const Map<Vector> &e);

	// Correction operator()(size_t no, size_t nv,
	// 		      const Data &data,
	// 		      const vector &eh, const vector &ep,
	// 		      int i, int j, int k,
	// 		      double* u1, double* t3);
	
    };


} // namespace cc
} // namespace cchem

#endif // MPQC_CC_CC_HPP
