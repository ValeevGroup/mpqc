#include "mpqc/cc/cc.hpp"
#include "mpqc/cc/tensor.hpp"
#include "mpqc/cc/symmetrize.hpp"
#include "mpqc/cc/utility.hpp"
#include "mpqc/cc/diis.hpp"

#include "mpqc/thread.hpp"
#include "mpqc/parallel/environment.hpp"
#include "mpqc/blas.hpp"
#include "mpqc/foreach.hpp"
#include "mpqc/utility/timer.hpp"
#include "mpqc/utility/array.hpp"

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/adaptors.hpp>
#include <stdio.h>

// #define BOOST_PROFILE_ENABLE
//#include "boost/utility/profiler.hpp"

namespace cchem {
namespace cc {

    void singles(const Parallel &pe,
		 const Wavefunction &wf,
		 const Map<const Array*> &V,
		 const Matrix &t1, Map<Matrix> &O,
		 const Map<const Array*> &O2) {

	BOOST_PROFILE_REGISTER_THREAD;

	namespace ublas = boost::numeric::ublas;
	using ublas::make_matrix;
	using ublas::make_vector;
	using ublas::column_major;

	using tensor::as_matrix;
	using tensor::as_vector;

	const Array &t2 = *O2["t(ijab)"];

	ublas::range ro(0, wf.active().size());
	ublas::range rv(0, wf.virtuals().size());

	const size_t no = ro.size();
	const size_t nv = rv.size();

	// I(ab) += (1 - \delta(a,b))*f(ab)
	{
	    O["I(ab)"] = trans(project(wf.F(), wf.virtuals(), wf.virtuals()));
	    for (size_t i = 0; i < nv; ++i) {
		O["I(ab)"](i,i) = 0;
	    }
	}

	O["'ij"] = project(wf.F(), wf.active(), wf.active());
	O["ia"] = project(wf.F(), wf.active(), wf.virtuals());
	for (size_t i = 0; i < no; ++i) {
	    O["'ij"](i,i) = 0;
	}
	O["ij"] = Matrix(no,no,0);

	Symbol< tensor_reference<4> > S;
	Buffer<double> buffer;
	buffer[0].resize(no*no*nv);
	buffer[1].resize(no*no*nv);

	// t(ia) += vt2
	{
	    BOOST_PROFILE_LINE;
	    tensor_reference<3> vt(buffer[0], boost::extents[nv][no][no]);
	    for (size_t k = 0; k < no; ++k) {
		size_t start[] = { 0, 0, 0, nv+k };
		size_t finish[] = { no, no, nv, nv+k+1 };
		O2["vt2"]->get(vt.data(), start, finish);
		symmetrize2<0,1>(vt, 2, -1);
		ublas::range rk(k*no, (k+1)*no);
		O["t(ia)"] += project(as_matrix<2,1>(vt), rk, rv);
	    }
	}

	// I(i,a) += (2*V(m,i,e,a) - V(i,m,e,a))*t(m,e)
	for (size_t u = 0; u < nv; ++u) {
	    BOOST_PROFILE_LINE;
	    S.load("v(ijab)", buffer[0], ro, ro, rv, u, *V["ijab"]);
	    BOOST_AUTO(v, S["v(ijab)"][0]);
	    symmetrize2<0,1>(v, -1, 2);
	    BOOST_AUTO(ov, column(O["ia"], u));
	    blas::gemv(1,
	    	       as_matrix<1,2>(S["v(ijab)"][0]),
	    	       as_vector(t1),
	    	       1, ov);
	}

	// I(b,a) = I(a,b)
	O["I(ab)"] = trans(Matrix(O["I(ab)"]));

	if (pe.rank() != 0) {
	    O["I(ab)"].clear();
	    O["'ij"].clear();
	    O["t(ia)"].clear();
	}

	Thread::Task<Parallel::Task&> task(pe.task());
	task.reset();

	//for (size_t u = 0; u < nv; ++u) {
	while (task++ < nv) {
	    size_t u = task;

	    // I(b,a) += (2*vt1' - vt1')
	    {
		BOOST_PROFILE_LINE;
		boost::multi_array_ref<double,3>
		    vt0(buffer[0], boost::extents[nv][no][no]),
		    vt1(buffer[1], boost::extents[nv][no][no]);
		size_t start[] = { 0, 0, 0, u };
		size_t finish[] = { no, no, nv, u+1 };
		O2["vt1\""]->get(vt0.data(), start, finish);
		O2["vt1\'"]->get(vt1.data(), start, finish);
		for (size_t v = 0; v < nv; ++v) {
		    for (size_t k = 0; k < no; ++k) {
			O["I(ab)"](u,v) += 2*vt0[v][k][k];
			O["I(ab)"](v,u) -= vt1[v][k][k];
		    }
		}
	    }

	    {
		BOOST_PROFILE_LINE;
		S.load("t(ijab)", buffer[1], ro, ro, rv, u, t2);
	    }

	    {
		BOOST_AUTO(t, S["t(ijab)"][0]);
		symmetrize2<0,1>(t, -1, 2);
	    }

	    // o(i,a) += I(m,e)*(2*t(i,m,e,a) - t(m,i,e,a))
	    {
	    	BOOST_AUTO(oa, column(O["t(ia)"], u));
	    	blas::gemv(1, as_matrix<1,2>(S["t(ijab)"][0]), as_vector(O["ia"]),
	    		   1, oa);
	    }

	    {
		BOOST_PROFILE_LINE;
		S.load("v(ijab)", buffer[0], ro, ro, rv, u, *V["ijab"]);
		S.load("v(iajb)", buffer[1], ro, rv, ro, u, *V["iajb"]);
	    }

	    // o(i,a) += 2*V(m,i,e,a)*t(m,e)
	    // o(i,a) -= V(i,e,m,a)*t(m,e)
	    {
	    	BOOST_AUTO(oa, column(O["t(ia)"],u));
	    	Matrix tt(trans(t1));
	    	S["v(ijab)"][0].permute<1,0,2>();
	    	blas::gemv(2, as_matrix<1,2>(S["v(ijab)"][0]), as_vector(t1),
	    	 	   1, oa);
	    	blas::gemv(-1, as_matrix<1,2>(S["v(iajb)"][0]), as_vector(tt),
	    	 	   1, oa);
	    	S["v(ijab)"][0].permute<1,0,2>();
	    }

	    S.load("t(ijab)", buffer[1], ro, ro, rv, u, t2);

	    // I(b,a) += (2*V(m,n,a,e) - V(n,m,a,e))*t(m,n,b,e)
	    // I'(i,j) += (2*V(i,m,f,e) - V(m,i,f,e))*t(j,m,f,e)
	    {
		BOOST_AUTO(v, S["v(ijab)"][0]);
	    	symmetrize2<0,1>(v, 2, -1);
	    	blas::gemm(-1,
	    		   trans(as_matrix<2,1>(v)),
			   (as_matrix<2,1>(S["t(ijab)"][0])),
	    		   1, O["I(ab)"]);
	    	blas::gemm(1,
	    		   as_matrix<1,2>(v),
	    		   trans(as_matrix<1,2>(S["t(ijab)"][0])),
	    		   1, O["'ij"]);
	    }
	    
	    //S.load("t(ijab)", t2, ro, ro, rv, u);
	    S.load("v(ijka)", buffer[0], ro, ro, ro, u, *V["ijka"]);

	    // o(i,a) -= v(n,m,i,e)*(2*t(n,m,a,e) - t(m,n,a,e)) 
	    {
		BOOST_AUTO(t, S["t(ijab)"][0]);
	    	symmetrize2<0,1>(t, 2, -1);
	    	blas::gemm(-1,
	    		   trans(as_matrix<2,1>(S["v(ijka)"][0])),
	    		   as_matrix<2,1>(t),
	    		   1, O["t(ia)"]);
	    }

	    // I'(i,j) += (2*V(i,m,j,e) - V(m,i,j,e))*t(m,e)
	    {
	    	BOOST_AUTO(v, S["v(ijka)"][0]);
	    	symmetrize2<0,1>(v, -1, 2);
	    	BOOST_AUTO(ov, as_vector(O["'ij"]));
	    	blas::gemv(1,
	    		   trans(as_matrix<1,2>(S["v(ijka)"][0])),
	    		   column(t1, u),
	    		   1, ov);
	    }

	}

	pe.reduce("+", O["I(ab)"].data());
	pe.reduce("+", O["t(ia)"].data());
	pe.reduce("+", O["'ij"].data());

	// I(i,j) = I'(i,j) + I(i,e)*t(j,e)
	for (size_t j = 0; j < no; ++j) {
	    BOOST_AUTO(oj, column(O["ij"], j));
	    oj = column(O["'ij"], j);
	    blas::gemv(1, O["ia"], column(trans(t1),j), 1, oj);
	}

	Matrix f = project(wf.F(), wf.active(), wf.virtuals()); 

	// I(b,a) += I(m,b)*t(m,a) !  papers says f(m,b)*t(m,a), code says I(m,b)*t(m,a)
	// o(i,a) += f(i,a)
	// o(i,a) += t(i,e)*I(a,e)'
	// o(i,a) += I'(m,i)*t(m,a)'
	for (size_t a = 0; a < nv; ++a) {
	    BOOST_AUTO(Iab, column(O["I(ab)"],a));
	    blas::gemv(-1, trans(O["ia"]), column(t1,a), 1, Iab);
	    BOOST_AUTO(oa, column(O["t(ia)"], a));
	    oa += column(f,a);
	    blas::gemv(1, t1, Iab, 1, oa);
	    blas::gemv(-1, trans(O["'ij"]), column(t1,a), 1, oa);
	}

	// I(a,b) = I(b,a)
	O["I(ab)"] = trans(Matrix(O["I(ab)"]));

	BOOST_AUTO(&t, O["t(ia)"]);
	Denominator eh(project(wf.F(), wf.active(), wf.active()));
	Denominator ep(project(wf.F(), wf.virtuals(), wf.virtuals()));
	for (size_t a = 0; a < nv; ++a) {
	    for (size_t i = 0; i < no; ++i) {
		t(i,a) /= (eh(i) - ep(a));
	    }
	}

	BOOST_PROFILE_DUMP(std::cout);    

	pe.barrier();

    }    

    // implemented in doubles.cpp
    void doubles(const Parallel &pe,
		 const Wavefunction &wf,
		 const Map<const Array*> &V,
		 const Matrix &t1,
		 const Map<Matrix> &A1,
		 const Map<Array*> &A,
		 double *rmax);

} // namespace cc
} // namespace cchem

double cchem::cc::sd::energy1(const Parallel &pe,
			      const Wavefunction &wf,
			      const Array &V,
			      const Matrix &t1, const Array &t2) {

    namespace ublas = boost::numeric::ublas;
    using tensor::as_vector;
    using tensor::as_matrix;

    ublas::range ro(0, wf.active().size());
    ublas::range rv(0, wf.virtuals().size());

    const size_t no = ro.size();
    const size_t nv = rv.size();

    Symbol< tensor_reference<4> > S;
    Buffer<double> buffer;
    buffer[0].resize(no*no*nv);
    buffer[1].resize(no*no*nv);

    Thread::Task<Parallel::Task&> task(pe.task());

    double E = 0;

    // E += (2*v(ijab) - v(jiab))*c(ijab)
    {
	Matrix c, u;
	task.reset();
	while (++task < nv) {
	    size_t b = task;
	    S.load("v(ijab)", buffer[0], ro, ro, rv, b, V);
	    S.load("t(ijab)", buffer[1], ro, ro, rv, b, t2);
	    for (size_t a = 0; a < nv; ++a) {
		BOOST_AUTO(v, (as_matrix<1,1>(S["v(ijab)"][0][a])));
		BOOST_AUTO(t, (as_matrix<1,1>(S["t(ijab)"][0][a])));
		c = (t + outer_prod(column(t1,a), column(t1,b)));
		u = (2*v - trans(v));
		E += inner_prod(as_vector(u), as_vector(c));
	    }
	}
    }
    pe.reduce("+", E);

    // E += 2*f(i,a)*t(i,a)
    Matrix f = project(wf.F(), wf.active(), wf.virtuals());
    E += 2*inner_prod(as_vector(f), as_vector(t1));

    //std::cout << "E(CCSD) " << E << std::endl;
    return E;
}


void cchem::cc::sd::energy1(const Parallel &pe,
			    const Wavefunction &wf,
			    const Map<const Array*> &V,
			    const Map<Array*> &A,
			    double *E, double *rmax) {

    namespace ublas = boost::numeric::ublas;
    using ublas::make_matrix;
    using ublas::make_vector;
    using ublas::column_major;

    using utility::make_array;

    const size_t no = wf.active().size();
    const size_t nv = wf.virtuals().size();

    ublas::range ro(0, no);
    ublas::range rv(0, nv);

    Matrix t1(no,nv,0);
    A["t(ia)"]->get(t1.data().begin(),
		    make_array<size_t>(0,0),
		    make_array<size_t>(no,nv));

    bool damp1 = 0;//(this->iteration_ < 3);
    if (damp1) t1.clear();

    Map<Matrix> O1;

    // partially evaluated directly
    O1["t(ia)"] = Matrix(no, nv, 0);
    O1["I(ab)"] = Matrix(nv, nv, 0);

    {
	utility::timer timer;
	integrals::Screening screening(wf.basis(), this->cutoff);
	sd::direct(pe, wf, t1, *A["t(ijab)"], O1, A, screening);
	pe.cout() << "ccsd vt time: " << timer << std::endl;
    }

    {
	utility::timer timer;
	singles(pe, wf, V, t1, O1, A);
	pe.cout() << "ccsd t1 time: " << timer << std::endl;
    }

    {
	utility::timer timer;
	*rmax = 0;
	doubles(pe, wf, V, t1, O1, A, rmax);
	pe.cout() << "ccsd t2 time: " << timer << std::endl;
    }
    if (damp1) O1["t(ia)"].clear();

    if (pe.rank() == 0) {
	A["t(ia)"]->put(O1["t(ia)"].data().begin(),
			make_array<size_t>(0,0),
			make_array<size_t>(no,nv));
	for (size_t a = 0; a < nv; ++a) {
	    for (size_t i = 0; i < no; ++i) {
		double r = t1(i,a) - O1["t(ia)"](i,a);
		*rmax = std::max(*rmax, fabs(r));
	    }
	}
    }
    pe.broadcast(*rmax, 0);
    A["t(ia)"]->flush();
    pe.barrier();

    *E = energy1(pe, wf, *V["ijab"], O1["t(ia)"], *A["t(ijab)"]);

    return;
}

double cchem::cc::sd::energy(const Parallel &pe,
			     const Wavefunction &wf,
			     const Map<const Array*> &V,
			     const Map<Array*> &A,
			     DIIS *diis) {

    double E = 0;
    size_t num_iter = this->max_iter;

    // std::cout << "xxx" << std::endl;
    // std::cout << num_iter << std::endl;
    // std::cout << convergence << std::endl;

    while (num_iter--) {

	utility::timer t;
	double Ei;
	double rmax;

	energy1(pe, wf, V, A, &Ei, &rmax);

	//double dE = (E - Ei);
	E = Ei;

	if (diis) {
	    utility::timer t;
	    if (this->iteration_ == 0) {
		diis->initialize(*A["t(ia)"], *A["t(ijab)"]);
		pe.cout() << "ccsd diis initialized" << t << std::endl;
	    }
	    else {
		rmax = diis->apply(*A["t(ia)"], *A["t(ijab)"]);
		pe.cout() << "ccsd diis time: " << t << std::endl;
	    }
	}

	A["t(ia)"]->flush();
	A["t(ijab)"]->flush();
	pe.broadcast(rmax, 0);
	pe.barrier();

	pe.cout() << "ccsd iter. " << this->iteration_
	   << "\t" << std::setprecision(10) << E
	   << "\t" << std::setprecision(10) << rmax
	   << "\t" << t
	   << std::endl;
	// printf("ccsd iter. %i: \t%13.10f\t%13.10f\n",
	//        this->iteration_, E, rmax);

	if (rmax <= this->convergence && this->iteration_)
	    goto converged;
	++this->iteration_;
    }

    pe.barrier();
    pe.cout() << "ccsd not converged" << E << std::endl;
    throw CCHEM_EXCEPTION("ccsd not converged");

 converged:
    pe.barrier();
    pe.cout() << "ccsd energy:" << "\t"
	      << std::setprecision(10) << E << std::endl;
    return E;

}




