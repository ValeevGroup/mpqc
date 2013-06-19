#include "mpqc/cc/cc.hpp"
#include "mpqc/cc/tensor.hpp"
#include "mpqc/cc/symmetrize.hpp"
#include "mpqc/cc/utility.hpp"

#include "mpqc/thread.hpp"
#include "mpqc/blas.hpp"
#include "mpqc/omp.hpp"
#include "mpqc/foreach.hpp"
#include "mpqc/utility/timer.hpp"
#include "mpqc/utility/array.hpp"

#ifdef HAVE_CUBLAS
//#warning HAVE_CUBLAS
#include "cublas.hpp"
#include <boost/numeric/bindings/cublas/cublas.hpp>
#endif


#include <boost/typeof/typeof.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/adaptors.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <set>
#include <stdio.h>

//#define BOOST_PROFILE_ENABLE
#include "boost/utility/profiler.hpp"

namespace cchem {

namespace cc {
namespace detail {


    template<typename T, class T1, class V>
    void add_t1(tensor_reference<3> t2,
		const T &alpha, const T1 &t1, const V &v) {
	for (size_t u = 0; u < t2.shape()[0]; ++u) {
	    BOOST_AUTO(c, (tensor::as_matrix<1,1>(t2[u])));
	    blas::ger(alpha, column(t1,u), v, c);
	}
    }

    template<typename T, class T1, class V>
    void add_t1(const T &alpha, const T1 &t1, const V &v,
		tensor_reference<3> t2) {
	for (size_t u = 0; u < t2.shape()[0]; ++u) {
	    BOOST_AUTO(c, (tensor::as_matrix<1,1>(t2[u])));
	    blas::ger(alpha, column(t1,u), v, c);
	}
    }

    template<typename T, class T1, class V>
    void add_t1(const T &alpha, const T1 &t1, const V &v,
		tensor_reference<4> &t2) {
	for (size_t u = 0; u < t2.shape()[0]; ++u) {
	    BOOST_AUTO(c, t2[u]);
	    add_t1(alpha, t1, column(v,u), c); 
	}
    }
    

    struct doubles {
    private:
	boost::ptr_map<int,Thread> threads_;
	omp::task<int> task_;
    public:
	void terms_with_t2(boost::numeric::ublas::range ro,
			   boost::numeric::ublas::range rv,
			   boost::numeric::ublas::range rb,
			   Symbol< tensor_reference<4> > &S,
			   const Array &t2);
	template<class T1>
	void terms_with_v2(boost::numeric::ublas::range ro,
			   boost::numeric::ublas::range rv,
			   boost::numeric::ublas::range rb,
			   Symbol< tensor_reference<4> > &S,
			   const Map<const Array*> &V,
			   const T1 &t1);
    };

	
    template<class T1>
    void doubles::terms_with_v2(boost::numeric::ublas::range ro,
				boost::numeric::ublas::range rv,
				boost::numeric::ublas::range rb,
				Symbol< tensor_reference<4> > &S,
				const Map<const Array*> &V,
				// column(t1,b)
				const T1 &t1) {

	using tensor::as_vector;
	using tensor::as_matrix;

#pragma omp critical
	{
	    // allocate thread data
	    this->threads_[omp::thread()];
	}
	Thread &my = this->threads_[omp::thread()];

	{
	    size_t no = ro.size();
	    size_t nv = rv.size();

	    my.buffer[0].resize(no*no*nv);
	    my.S.clear();
	}

	this->task_.reset();

	// terms with v(i,j,k,a)
	for (size_t u = 0, next = this->task_++; u < rv.size(); ++u) {
	    BOOST_PROFILE_LINE;
	    
	    if (u != next) continue;
	    next = this->task_++;

	    {
		my.S.load("v(ijka)", my.buffer[0], ro, ro, ro, u, *V["ijka"]);
	    }
	    BOOST_AUTO(v, my.S["v(ijka)"][0]);


	    // I(i,j,a,b) -= v(m,i,j,b)*t(m,a)
	    for (size_t b = 0; b < rb.size(); ++b) {
		BOOST_AUTO(I, (as_vector(S["I(ijab)"][b][u])));
		blas::gemv(-1, trans(as_matrix<1,2>(v)),
			   column(t1,b),
			   1, I);
	    }

	    tensor::permute<0,2,1>(v);

	    // I(i,a,j,b) -= v(i,m,j,b)*t(m,a)
	    //            -= v(j,m,i,b)*t(m,a)
	    for (size_t b = 0; b < rb.size(); ++b) {
		BOOST_AUTO(I, (as_vector(S["I(iajb)"][b][u])));
		blas::gemv(-1, (as_matrix<2,1 >(v)),
			   column(t1,b),
			   1, I);
	    }

	    my.S.erase("v(ijka)");
	}

#pragma omp barrier
#pragma omp single
	{	
	    tensor::permute<1,0,2,3>(S["I(ijab)"]);
	}

	this->task_.reset();
	for (size_t u = 0, next = this->task_++; u < rv.size(); ++u) {
	    BOOST_PROFILE_LINE;

	    if (u != next) continue;
	    next = this->task_++;

	    {
		BOOST_PROFILE_LINE;
		my.S.load("v(ijab)", my.buffer[0], ro, ro, rv, u, *V["ijab"]);
	    }
	    BOOST_AUTO(v, (my.S["v(ijab)"][0]));

	    for (size_t b = 0; b < rb.size(); ++b) {

		// I(i,a,j,b) -= 0.5*V(i,m,e,b)*c(m,j,a,e)
		//            -= 0.5*V(m,i,b,e)*c(j,m,e,a)
		{
		    BOOST_AUTO(I, (as_matrix<1,1>(S["I(iajb)"][b][u])));
		    blas::gemm(-0.5,
			       (as_matrix<1,2>(S["c(ijab)"][b])),
			       trans(as_matrix<1,2>(v)),
			       1, I);
		}
	    }

	    tensor::permute<1,0,2>(v);

	    // I(i,j,a,b) -= 0.5*v(i,m,b,e)*c(j,m,e,a)
	    for (size_t b = 0; b < rb.size(); ++b) {
		BOOST_AUTO(I, (as_matrix<1,1>(S["I(ijab)"][b][u])));
		blas::gemm(-0.5,
			   (as_matrix<1,2>(S["c(ijab)"][b])),
			   trans(as_matrix<1,2>(v)),
			   1, I);
	    }
	    
	    symmetrize2<0,1>(v, 1.0, -0.5);

	    // I(i,j,a,b) += (v(i,m,b,e) - 0.5*v(i,m,e,b))*t(j,m,a,e)
	    //            += -0.5*(-2*v(i,m,b,e) + v(m,i,b,e))*t(m,j,e,a)
	    for (size_t b = 0; b < rb.size(); ++b) {
		// S["t(ijab)"][0].permute<1,0,2>();
		BOOST_AUTO(I, (as_matrix<1,1>(S["I(ijab)"][b][u])));
		blas::gemm(1,
			   (as_matrix<1,2>(S["t(jiab)"][b])),
			   trans(as_matrix<1,2>(v)),
			   1, I);
		// S["t(ijab)"][0].permute<1,0,2>();
	    }

	} // u index loop

#pragma omp barrier
#pragma omp single
	{	
	    tensor::permute<1,0,2,3>(S["I(ijab)"]);
	}

    }


    void doubles::terms_with_t2(boost::numeric::ublas::range ro,
				boost::numeric::ublas::range rv,
				boost::numeric::ublas::range rb,
				Symbol< tensor_reference<4> > &S,
				const Array &t2) {

	// -0.7940365155

	using tensor::as_vector;
	using tensor::as_matrix;

#pragma omp barrier
#pragma omp single
	{
	    tensor::permute<1,0,2,3>(S["I(iajb)"]);
	    tensor::permute<1,0,2,3>(S["I(ijab)"]);
	    tensor::permute<1,0,2,3>(S["Dt"]);
	}

#pragma omp critical
	{
	    // allocate thread data
	    this->threads_[omp::thread()];
	}
	Thread &my = this->threads_[omp::thread()];

	size_t no = ro.size();
	size_t nv = rv.size();
	// size_t nb = rb.size();

	{
	    my.S.clear();

	    my.buffer[-1].resize(no*no*nv);
	    my.S.set("host.t", my.buffer[-1], ro, ro, rv, 0);

	    my.S.set("I(ijab)", S["I(ijab)"].data(), ro, ro, rv, rb);
	    my.S.set("I(iajb)", S["I(iajb)"].data(), ro, ro, rv, rb);
	}


	this->task_.reset();
	size_t next = this->task_++;

	// terms with t2
	for (size_t u = 0; u < rv.size(); ++u) {
	    if (u != next) continue;
	    next = this->task_++;

	    {
		BOOST_PROFILE_LINE;
		my.S.load("host.t", my.S["host.t"].data(), ro, ro, rv, u, t2);
	    }

	    BOOST_AUTO(t, my.S["host.t"][0]);	

	    // Dt(j,i,a,b) -= t(m,j,a,e)*I(i,e,m,b)
	    for (size_t b = 0; b < rb.size(); ++b) {
		BOOST_PROFILE_LINE;
		BOOST_AUTO(Dt, (as_matrix<1,1>(S["Dt"][b][u])));
		blas::gemm(-1,
			   (as_matrix<1,2>(t)),
			   trans(as_matrix<1,2>(my.S["I(iajb)"][b])),
			   1, Dt);
	    }

	    tensor::permute<1,0,2>(t);
	
	    for (size_t b = 0; b < rb.size(); ++b) {

		BOOST_AUTO(Dt, (as_matrix<1,1>(S["Dt"][b][u])));

		// Dt(i,j,b,a) -= I(m,a,i,e)*t(m,j,e,b)
		//             -= I(m,a,i,e)*t(j,m,b,e)
		{
		    BOOST_PROFILE_LINE;
		    blas::gemm(-1,
			       (as_matrix<1,2>(my.S["I(iajb)"][b])),
			       trans(as_matrix<1,2>(t)),
			       1, Dt);
		}

	    }

	    // t(m,i,e,a) = (2*t(m,i,e,a) - t(i,m,e,a))
	    {
		BOOST_PROFILE_LINE;
		symmetrize2<0,1>(t, 2.0, -1.0);
	    }

	    // Dt(j,i,a,b) += (2*t(m,i,e,a) - t(i,m,e,a))*I(m,j,b,e)
	    for (size_t b = 0; b < rb.size(); ++b) {
		BOOST_PROFILE_LINE;
		BOOST_AUTO(Dt, (as_matrix<1,1>(S["Dt"][b][u])));
		blas::gemm(1,
			   (as_matrix<1,2>(my.S["I(ijab)"][b])),
			   trans(as_matrix<1,2>(t)),
			   1, Dt);
	    }


	} // u index loop

#pragma omp barrier

#pragma omp single
	{
	    tensor::permute<1,0,2,3>(S["I(iajb)"]);
	    tensor::permute<1,0,2,3>(S["I(ijab)"]);
	    tensor::permute<1,0,2,3>(S["Dt"]);
	}

    }


} // namespace detail
} // namespace cc


namespace cc {

    // declared in sd.cpp
    void doubles(const Parallel &pe,
		 const Wavefunction &wf,
		 const Map<const Array*> &V,
		 const Matrix &t1,
		 const Map<Matrix> &A1,
		 const Map<Array*> &A,
		 double *rmax) {

	using cc::detail::add_t1;

	namespace ublas = boost::numeric::ublas;
	using ublas::make_matrix;
	using ublas::make_vector;
	using ublas::column_major;

	using tensor::as_vector;
	using tensor::as_matrix;

	using utility::make_array;

	ublas::range ro(0, wf.active().size());
	ublas::range rv(0, wf.virtuals().size());

	const size_t no = ro.size();
	const size_t nv = rv.size();

	BOOST_PROFILE_REGISTER_THREAD;

	//std::cout << "# CC doubles" << std::endl;

	struct {
	    utility::timer total;
	    utility::timer::value_type term[2];
	} time;

	Buffer<double> buffer;
	Symbol< tensor_reference<4> > S;

	utility::timer timer;

	Thread::Task<Parallel::Task&> task(pe.task());

	// t(i,j,a,b) += 0.5(vt2 + vt1)
	{
	    BOOST_PROFILE_LINE;
	    Tensor<3> vt(boost::extents[nv][no][no]);
	    Tensor<3> t(boost::extents[nv][no][no]);
	    task.reset();
	    while (++task < nv) {
		BOOST_PROFILE_LINE;
		int v = task;
		size_t start[] = { 0, 0, 0, v };
		size_t finish[] = { no, no, nv, v+1 };
		A["vt2"]->get(t.data(), start, finish);

		A["vt1"]->get(vt.data(), start, finish);
		as_vector(t).plus_assign(as_vector(vt));
		as_vector(t) *= 0.5;
		//as_vector(t).plus_assign(as_vector(vt));
		A["u(ijab)"]->put(t.data(), start, finish);
	    }
	}
	A["u(ijab)"]->flush();
	pe.barrier();

	{
	    BOOST_PROFILE_LINE;

	    buffer[0].resize(no*no*nv);
	    buffer[1].resize(no*no*nv);

	    Tensor<4> I(boost::extents[no][no][no][no]);
	    tensor::fill(I, 0);

	    // I(i,j,k,l) = V(i,j,k,l) + v(i,j,e,f)*c(e,f,k,l) + P(t(k,e)*V(i,j,e,l))
	    //            = V(i,j,k,l) + v(i,j,e,f)*c(e,f,k,l) + P(t(k,e)*V(j,i,l,e))

	    // I(i,j,k,l) = t(k,e)*V(i,j,e,l)
	    //              t(k,e)*V(j,i,l,e)
	    task.reset();
	    while (++task < nv) {
		size_t v = task;
		S.load("v(ijka)", buffer[0], ro, ro, ro, v, *V["ijka"]);
		BOOST_AUTO(Im, (as_matrix<3,1>(I)));
		blas::ger(1, as_vector(S["v(ijka)"]), column(t1,v), Im);
	    }

	    // I(i,j,k,l) = P(I(i,j,k,l))
	    BOOST_AUTO(tr1, make_matrix<column_major>(no, no, buffer[0].data()));
	    for (size_t l = 0; l < no; ++l) {
	    	for (size_t k = 0; k <= l; ++k) {
		    BOOST_AUTO(lk, (as_matrix<1,1>(I[l][k])));
		    BOOST_AUTO(kl, (as_matrix<1,1>(I[k][l])));
		    tr1.assign(trans(kl));
	    	    lk.plus_assign(tr1);
		    kl.assign(trans(lk));
	    	}
	    }

	    task.reset();
	    while (++task < no) {
		// I(ijkl) += vt(klij) + vt2(klij)
		size_t j = task;
		ublas::range ri(nv, nv+no);
		S.load("vt", buffer[0], ro, ro, ri, j+nv, *A["vt1"]);
		S.load("vt2", buffer[1], ro, ro, ri, j+nv, *A["vt2"]);
		as_vector(S["vt"]).plus_assign(as_vector(S["vt2"]));
		tensor::permute<0,2,1,3>(S["vt"]);
		tensor::permute<1,0,2,3>(S["vt"]);
		for (size_t l = 0; l < no; ++l) {
		    for (size_t k = 0; k < no; ++k) {
			BOOST_AUTO(Ij, as_vector(I[l][k][j]));
			Ij.plus_assign(as_vector(S["vt"][0][l][k]));
		    }
		}
		// I(ijkl) += v(ijkl)
		size_t l = task;
	    	S.load("v(ijkl)", buffer[0], ro, ro, ro, l, *V["ijkl"]);
	    	as_vector(I[l]).plus_assign(as_vector(S["v(ijkl)"]));
	    }

	    pe.reduce("+", I.data(), I.num_elements());

	    task.reset();
	    while (++task < nv) {
		size_t v = task;

		S.load("t(ijab)", buffer[0], ro, ro, rv, v, *A["t(ijab)"]);
		S.load("Dt", buffer[1], ro, ro, rv, v, *A["u(ijab)"]);
		//tensor::fill(S["Dt"], 0); // %%%

		{
		    //S["Dt"][0].permute<1,0,2>();
		    BOOST_AUTO(Dt, (as_matrix<1,2>(S["Dt"][0])));
		    // t(j,i,b,a) += t(i,m,a,b)*I(m,j)
		    //            += t(m,i,b,a)*I(m,j)
		    blas::gemm(-1,
		    	       trans(A1["ij"]),
		    	       as_matrix<1,2>(S["t(ijab)"][0]),
		    	       1, Dt);
		    //S["Dt"][0].permute<1,0,2>();
		}

		{
		    BOOST_AUTO(Dt, (as_matrix<2,1>(S["Dt"][0])));
		    // t(j,i,b,a) += t(j,i,e,a)*I(e,b)
		    blas::gemm(1,
		    	       as_matrix<2,1>(S["t(ijab)"][0]),
		    	       trans(A1["I(ab)"]),
		    	       1, Dt);
		    // t(j,i,b,a) += 0.5*c(n,m,b,a)*I(n,m,j,i)
		    add_t1(S["t(ijab)"][0], 1.0, t1, column(t1,v));
		    blas::gemm(0.5,
			       trans(as_matrix<2,2>(I)),
			       as_matrix<2,1>(S["t(ijab)"][0]),
			       1, Dt);
		}
		S.store("Dt", ro, ro, rv, v, *A["u(ijab)"]);	    
	    }
	}
	A["u(ijab)"]->flush();
	pe.barrier();
	//std::cout << "    v(ijkl) term: " << timer << std::endl;

	{
	    using std::max;

	    S.clear();

	    size_t B = 10;
	    detail::doubles doubles;

	    buffer[0].resize(B*no*no*nv);
	    buffer[1].resize(B*2*no*no*nv);
	    buffer[2].resize(B*no*no*nv);

	    task.reset();
	    while (++task < (nv+B-1)/B) {
		
		ublas::range rb(task*B, std::min(nv, task*B+B));

		{
		    double *ptr = buffer[1].data();
		    S.set("I(iajb)", ptr, ro, ro, rv, rb);
		    ptr += B*no*no*nv;
		    S.set("I(ijab)", ptr, ro, ro, rv, rb);
		}

		tensor::fill(S["I(iajb)"], 0);
		tensor::fill(S["I(ijab)"], 0);

		{
		    S.load("t(ijab)", buffer[0], ro, ro, rv, rb, *A["t(ijab)"]);
		    S.set("c(ijab)", buffer[2], ro, ro, rv, rb);

		    tensor::copy(S["t(ijab)"], S["c(ijab)"]);
		    add_t1(2.0, t1, project(t1,ro,rb), S["c(ijab)"]);

		    tensor::permute<1,0,2,3>(S["t(ijab)"]);
		    S.set("t(jiab)", S["t(ijab)"].data(), ro, ro, rv, rb);
		    S.erase("t(ijab)");

		    BOOST_AUTO(tm, project(t1, ro, rb));

		    timer = utility::timer();
#pragma omp parallel
		    {
			// terms with v(i,j,a,b)
			doubles.terms_with_v2(ro, rv, rb, S, V, tm);
		    }
		    tensor::permute<1,0,2,3>(S["I(ijab)"]);
		    time.term[0] += timer;
		    
		    S.erase("t(jiab)");
		    S.erase("c(ijab)");
		}

		S.set("vt", buffer[2], ro, ro, rv, rb);

		// I(i,j,a,b) += v(i,j,b,a) + vt'
		{
		    BOOST_PROFILE_LINE;
		    S.load("v(ijab)", buffer[0], ro, ro, rv, rb, *V["ijab"]);

		    BOOST_AUTO(I, as_vector(S["I(ijab)"]));
		    tensor::permute<1,0,2,3>(S["v(ijab)"]);
		    blas::axpy(1, as_vector(S["v(ijab)"]), I);
		    S.load("vt1", S["vt"].data(),
			   ro, ro, rv, rb, *A["vt1\'"]);
		    blas::axpy(1, as_vector(S["vt1"]), I);
		
		    S.erase("v(ijab)");
		}

		// I(i,a,j,b) += vt1" + v(i,a,j,b)
		{
		    BOOST_PROFILE_LINE;
		    S.load("v(iajb)", buffer[0], ro, rv, ro, rb, *V["iajb"]);

		    BOOST_AUTO(I, as_vector(S["I(iajb)"]));
		    S.load("vt1", S["vt"].data(),
			   ro, ro, rv, rb, *A["vt1\""]);
		    blas::axpy(1, as_vector(S["vt1"]), I);
		    tensor::permute<1,0,2,3>(S["I(iajb)"]);
		    for (size_t b = 0; b < rb.size(); ++b) {
			S["I(iajb)"][b].plus_assign<0,2,1>(S["v(iajb)"][b]);
		    }

		    S.erase("v(iajb)");
		}

		S.load("Dt", buffer[0], ro, ro, rv, rb, *A["u(ijab)"]);

		{
		    BOOST_PROFILE_LINE;
		    ublas::range rm(nv,nv+no);
		    // t(i,j,a,b) -= t(m,b)*vt'(i,j,a,m)
		    // t(j,i,b,a) -= t(m,b)*vt'(j,i,m,a)
		    for (size_t b = 0; b < rb.size(); ++b) {   
		    	S.load("vt", S["vt"].data(),
		    	       ro, ro, rm, b+rb.start(), *A["vt1\'"]);
			tensor::permute<1,0,2,3>(S["vt"]);
		    	BOOST_AUTO(Dt, (as_matrix<2,1>(S["Dt"][b])));
		    	blas::gemm(-1,
		    		   as_matrix<2,2>(S["vt"]), t1,
		    		   1, Dt);
		    }
		}

		{
		    BOOST_PROFILE_LINE;
		    ublas::range rm(nv,nv+no);
		    // t(j,i,b,a) -= t(m,a)*vt'(j,i,m,b)
		    for (size_t b = 0; b < rb.size(); ++b) {   
		    	S.load("vt", S["vt"].data(),
		    	       ro, ro, b+rb.start(), rm, *A["vt1\'"]);
		    	BOOST_AUTO(Dt, (as_matrix<2,1>(S["Dt"][b])));
		    	blas::gemm(-1,
		    		   as_matrix<2,2>(S["vt"]), t1,
		    		   1, Dt);
		    }
		}

		// Dt(j,i,b,a) += vt1'
		{
		    BOOST_PROFILE_LINE;
		    S.load("vt1'", S["vt"].data(),
			   ro, ro, rv, rb, *A["vt1\'"]);
		    tensor::permute<1,0,2,3>(S["vt1'"]);
		    as_vector(S["Dt"]).plus_assign(as_vector(S["vt1'"]));
		}

		tensor::permute<1,0,2,3>(S["I(ijab)"]);

		// terms with t2
		{
		    timer = utility::timer();
#pragma omp parallel
		    {
			// terms with t2
			doubles.terms_with_t2(ro, rv, rb, S, *A["t(ijab)"]);
		    }
		    //std::cout << "t2 " << b << ": " << timer << std::endl;
		    time.term[1] += timer;
		}

		{
		    S.set("I'(jkia)", buffer[1], ro, ro, ro, rb);
		    S.load("I'(jkia)", S["I'(jkia)"].data(),
		    	   ro, ro, ro, rb, *V["ijka"]);
		}

		// I'(j,k,i,a) += vt1 // v(i,a,e,f)*t(j,e)*t(k,f)
		{
		    BOOST_PROFILE_LINE;
		    S.load("vt", S["vt"].data(),
			   ro, ro, ublas::range(nv, nv+no), rb, *A["vt1"]);
		    BOOST_AUTO(I, as_vector(S["I'(jkia)"]));
		    blas::axpy(1, as_vector(S["vt"]), I);
		}

		// I'(j,k,i,a) += vt2
		{
		    BOOST_PROFILE_LINE;
		    S.load("vt", S["vt"].data(),
			   ro, ro, ublas::range(nv, nv+no), rb, *A["vt2"]);
		    BOOST_AUTO(I, as_vector(S["I'(jkia)"]));
		    blas::axpy(1, as_vector(S["vt"]), I);
		}

		
		// Dt(i,j,a,b) -= t(a,m)*I'(i,j,m,b)
		for (size_t b = 0; b < rb.size(); ++b) {
		    BOOST_AUTO(Dt, (as_matrix<2,1>(S["Dt"][b])));
		    blas::gemm(-1,
		    	       as_matrix<2,1>(S["I'(jkia)"][b]), t1,
		    	       1, Dt);
		}

		S.store("Dt", ro, ro, rv, rb, *A["u(ijab)"]);

	    }

	}
	A["u(ijab)"]->flush();
	pe.barrier();

	S.clear();
	buffer.clear();

	// t(ijab) = P(Dt(ijab))/D
	{
	    BOOST_PROFILE_LINE;

	    *rmax = 0;
	    Denominator eh(project(wf.F(), wf.active(), wf.active()));
	    Denominator ep(project(wf.F(), wf.virtuals(), wf.virtuals()));

	    buffer[0].resize(no*no*nv);
	    buffer[1].resize(no*no*nv);
	    buffer[2].resize(no*no*nv);

	    task.reset();
	    while (++task < nv) {
		size_t b = task;
		S.load("ija", buffer[0], ro, ro, rv, b, *A["u(ijab)"]);
		S.load("jia", buffer[1], ro, ro, b, rv, *A["u(ijab)"]);
		S.load("v(ijab)", buffer[2], ro, ro, rv, b, *V["ijab"]);
		for (size_t a = 0; a < nv; ++a) {
		    BOOST_AUTO(tij, (as_matrix<1,1>(S["ija"][0][a])));
		    BOOST_AUTO(tji, (as_matrix<1,1>(S["jia"][a][0])));
		    BOOST_AUTO(v, (as_matrix<1,1>(S["v(ijab)"][0][a])));
		    tij.plus_assign(trans(tji) + v);
		    for (size_t j = 0; j < no; ++j) {
			for (size_t i = 0; i < no; ++i) {
			    double e = (eh(i,j) - ep(a,b));
			    assert(fabs(e) > 1e-100);
			    tij(i,j) /= e;
			}
		    }
		}
		S.load("t", buffer[1], ro, ro, rv, b, *A["t(ijab)"]);
		{
		    double *t = S["t"].data();
		    double *u = S["ija"].data();
		    double r = 0;
		    for (size_t i = 0; i < no*no*nv; ++i) {
			r = std::max(r, fabs(t[i] - u[i]));
		    }
		    *rmax = std::max(*rmax, r);
		}
		S.store("ija", ro, ro, rv, b, *A["t(ijab)"]);
	    }
	}

	A["t(ijab)"]->flush();
	pe.reduce("max", *rmax);
	pe.barrier();

	BOOST_PROFILE_DUMP(std::cout);

    }    

} // namespace cc
} // namespace cchem

