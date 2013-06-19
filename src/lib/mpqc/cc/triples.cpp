//#include "mpqc/config.hpp"

#include "mpqc/cc/cc.hpp"
#include "mpqc/cc/tensor.hpp"
#include "mpqc/cc/utility.hpp"
#include "mpqc/cc/detail/map.hpp"
#include "mpqc/cc/detail/permute.hpp"

#include "mpqc/foreach.hpp"
#include "mpqc/array/hdf5.hpp"

#include "mpqc/thread.hpp"
#include "mpqc/blas.hpp"
#include "mpqc/omp.hpp"
#include "mpqc/runtime.hpp"
#include "mpqc/utility/timer.hpp"
#include "mpqc/utility/progress.hpp"


#include <utility>
#include <stdexcept>
#include <boost/typeof/typeof.hpp>

#include "boost/utility/profiler.hpp" 

namespace cchem {

namespace cc {
namespace triples {
namespace detail {
    
    struct T3 {
    private:
	tensor_reference<6> t_;
    public:
	typedef Symbol< tensor_reference<4> >::range range;
	T3(size_t n, size_t N, double *data)
	    : t_(data, boost::extents[N][N][N][n][n][n]) {}

	static int min(int a, range r) {
	    return std::min<int>(a, r.finish()-1);
	}

	template<class F>
	size_t apply(F &f, Symbol< tensor_reference<4> > &S,
		     range ra, range rb, range rc) {
	    size_t ops = 0;
	    omp::task<int> task;
#pragma omp parallel reduction(+:ops)
	    for (int it = 0,
		     next = task++,
		     c = rc.start(); c < rc.finish(); ++c) {
		for (int b = rb.start(); b <= min(c,rb); ++b) {
		    for (int a = ra.start(); a <= min(b,ra); ++a) {
			if (it++ != next) continue;
			next = task++;	
			BOOST_AUTO(t, t_[c-rc.start()][b-rb.start()][a-ra.start()]);
			ops += f(t, S, ra, rb, rc, a, b, c);
		    }
		}
	    }
	    return ops;
	}

	size_t evaluate(range ra, range rb, range rc,
			Symbol< tensor_reference<4> > &S) {
	    size_t ops = 0;
	    omp::task<int> task;
#pragma omp parallel reduction(+:ops)
	    for (int it = 0,
		     next = task++,
		     c = rc.start(); c < rc.finish(); ++c) {
		for (int b = rb.start(); b <= min(c,rb); ++b) {
		    for (int a = ra.start(); a <= min(b,ra); ++a) {
			if (it++ != next) continue;
			next = task++;
			BOOST_AUTO(t, t_[c-rc.start()][b-rb.start()][a-ra.start()]);
			tensor::fill(t, 0);
			ops += evaluate(ra, rb, rc, a, b, c, S, t);
		    }
		}
	    }
	    return ops;
	}

    private:

	size_t evaluate(range ra, range rb, range rc,
			int A, int B, int C,
			Symbol< tensor_reference<4> > &S,
			tensor_reference<3> &t) {

	    int a = A - ra.start();
	    int b = B - rb.start();
	    int c = C - rc.start();

	    size_t ops = 0;

	    // t(i,j,k) = t(i,j,e,a) V(e,k,b,c) - t(i,m,a,b) V(j,k,m,c)
	    // t(i,k,j) = t(i,k,e,a) V(e,j,c,b) - t(i,m,a,c) V(k,j,m,b)
	    ops +=
		evaluate(S["t(oova)"][a], S["V(ovbc)"][c][b],
			 S["t(oova)"][a], S["V(ovcb)"][b][c],
			 S["t(oova)"][a][B], S["V(oooc)"][c],
			 S["t(oova)"][a][C], S["V(ooob)"][b],
			 t);

 	    // t(k,i,j) = t(k,i,e,c) V(e,j,a,b) - t(k,m,c,a) V(i,j,m,b)
	    // t(k,j,i) = t(k,j,e,c) V(e,i,b,a) - t(k,m,c,b) V(j,i,m,a)
	    ops +=
	    	evaluate(S["t(oovc)"][c], S["V(ovab)"][b][a],
	    		 S["t(oovc)"][c], S["V(ovba)"][a][b],
	    		 S["t(oovc)"][c][A], S["V(ooob)"][b],
	    		 S["t(oovc)"][c][B], S["V(oooa)"][a],
	    		 t);

	    // // t(j,k,i) = t(j,k,e,b) V(e,i,a,c) - t(j,m,b,c) V(k,i,m,a)
	    // // t(j,i,k) = t(j,i,e,b) V(e,k,c,a) - t(j,m,b,a) V(i,k,m,c)
	    ops +=
	    	evaluate(S["t(oovb)"][b], S["V(ovca)"][a][c],
	    		 S["t(oovb)"][b], S["V(ovac)"][c][a],
	    		 S["t(oovb)"][b][C], S["V(oooa)"][a],
	    		 S["t(oovb)"][b][A], S["V(oooc)"][c],
	    		 t);

	    return ops;
	}

	size_t evaluate(tensor_reference<3> t_ijae, tensor_reference<2> V_kebc,
			tensor_reference<3> t_ikae, tensor_reference<2> V_jecb,
			tensor_reference<2> t_imab, tensor_reference<3> V_jkmc,
			tensor_reference<2> t_imac, tensor_reference<3> V_kjmb,
			tensor_reference<3> t) {
	    using boost::numeric::ublas::trans;
	    using tensor::as_matrix;
	    size_t ops = 0;

	    ops += contract(1, as_matrix<2,1>(t_ijae), trans(as_matrix<1,1>(V_kebc)),
	    		    as_matrix<2,1>(t));
	    ops += contract(-1, as_matrix<1,1>(t_imab), trans(as_matrix<2,1>(V_jkmc)),
	    		    as_matrix<1,2>(t));
	    tensor::permute<0,2,1>(t);

	    ops += contract(1, as_matrix<2,1>(t_ikae), trans(as_matrix<1,1>(V_jecb)),
	    		    as_matrix<2,1>(t));
	    ops += contract(-1, as_matrix<1,1>(t_imac), trans(as_matrix<2,1>(V_kjmb)),
	    		    as_matrix<1,2>(t));
	    tensor::permute<1,0,2>(t);

	    return ops;
	}	
 

	template<class A, class B, class C>
	static size_t contract(double alpha, const A &a, const B &b, C c) {
	    assert(a.size1() == c.size1());
	    assert(b.size2() == c.size2());
	    assert(a.size2() == b.size1());
	    blas::set_num_threads(1);
	    blas::gemm(alpha, a, b, 1, c);
	    return 2*(a.size1()*a.size2()*b.size2());
	}

    };


    struct Energy  {
	const size_t no, nv;
	const Vector eh, ep;
    private:
	Matrix t1_, u_;
	Correction corr_;
	struct Partial {
	    Correction corr;
	    Vector u[3];
	    explicit Partial(size_t n) {
		foreach (Vector& v, u) {
		    v.resize(n);
		    v.clear();
		}
	    }
	};
	template<class T>
	struct V { T ab, ac, ba, bc, ca, cb; };

    public:
	Energy(const Array &t1, const Vector &eh, const Vector &ep) :
	    no(eh.size()), nv(ep.size()), eh(eh), ep(ep)
	{
	    t1_.resize(no, nv, false);
	    u_.resize(no, nv, false);
	    u_.clear();
	    size_t start[] = { 0, 0 };
	    size_t stop[] = { no, nv };
	    t1.get(t1_.data().begin(), start, stop);
	}

	Correction corr(const Parallel &pe) const {
	    double ets = 0;
	    Matrix u = u_;
	    for (size_t i = 0; i < (nv*no); ++i) {
	    	ets += t1_.data()[i]*u.data()[i];
	    }
	    Correction corr = corr_;
	    corr.ets += 2*ets;
	    pe.reduce("+", corr.ets);
	    pe.reduce("+", corr.etd);
	    pe.reduce("+", corr.ots);
	    pe.reduce("+", corr.otd);
	    return corr;
	}

	template<class T3>
	size_t operator()(const T3 &t3, Symbol< tensor_reference<4> > &cache,
			  Symbol< tensor_reference<4> >::range ra,
			  Symbol< tensor_reference<4> >::range rb,
			  Symbol< tensor_reference<4> >::range rc,
			  int a, int b, int c) {
	    double S = 1.0/(1 + (a == b || b == c));
	    S *= !(a == b && b == c); // zero diagonal
	    double D = ep(a) + ep(b) + ep(c);

	    Partial p = evaluate(no, S, D,
				 cache["V(oovb)"][b-rb.start()],
				 cache["V(oovc)"][c-rc.start()],
				 t3, eh, a, b, c);
#pragma omp critical
	    {
		using boost::numeric::ublas::column;
		column(u_, a) += p.u[0];
		column(u_, b) += p.u[1];
		column(u_, c) += p.u[2];
		corr_ += p.corr;
	    }
	    return 0;
	}

    private:
	template<class V2, class T3>
	static
	Partial evaluate(size_t no, double S, double Dp,
			 const V2 &Vb, const V2 &Vc,
			 const T3 &t3, const Vector &eh,
			 int a, int b, int c) {


#define t3(i,j,k) t3[k][j][i]
#define t(i,j,k) t ## i ## j ## k
#define V(i,j,a,b) (V ## b [a][j][i])

	    Partial p(no);
	    int n = no;
	    for (int k = 0; k < n; ++k) {
		for (int j = 0; j < n; ++j) {
		    double djk = eh(j) + eh(k);

		    struct { double bc, ac, ab; }
		    Vjk = { V(j,k,b,c), V(j,k,a,c), V(j,k,a,b) };
		    struct { double bc, ac, ab; }
		    Vkj = { V(k,j,b,c), V(k,j,a,c), V(k,j,a,b) };

		    for (int i = 0; i < n; ++i) {
			double D = S/(eh(i) + djk - Dp);

			double t(i,j,k) = t3(i,j,k);
			double t(i,k,j) = t3(i,k,j);
			double t(j,i,k) = t3(j,i,k);
			double t(j,k,i) = t3(j,k,i);
			double t(k,i,j) = t3(k,i,j);
			double t(k,j,i) = t3(k,j,i);

			double f = (8*t(i,j,k) -
				    4*(t(i,k,j)+t(k,j,i)+t(j,i,k)) +
				    2*(t(j,k,i)+t(k,i,j)));
			p.corr.etd += D*f*t(i,j,k);

			// ET[S]

			double Z[] = {
			    (2*t(i,j,k) + t(k,i,j)) - (2*t(j,i,k) + t(i,k,j)), // abc
			    (2*t(i,k,j) + t(k,j,i)) - (2*t(j,k,i) + t(i,j,k)), // acb
			    (2*t(j,i,k) + t(i,k,j)) - (2*t(i,j,k) + t(k,i,j)), // bac
			    (2*t(k,i,j) + t(j,k,i)) - (2*t(k,j,i) + t(j,i,k)), // bca
			    (2*t(j,k,i) + t(i,j,k)) - (2*t(i,k,j) + t(k,j,i)), // caa
			    (2*t(k,j,i) + t(j,i,k)) - (2*t(k,i,j) + t(j,k,i))  // cba
			};

			p.u[0](i) += D*(Z[0]*Vjk.bc);
			p.u[0](i) += D*(Z[1]*Vkj.bc);
			p.u[1](i) += D*(Z[2]*Vjk.ac);
			p.u[1](i) += D*(Z[3]*Vkj.ac);
			p.u[2](i) += D*(Z[4]*Vjk.ab);
			p.u[2](i) += D*(Z[5]*Vkj.ab);

		    }
		}
	    }
#undef t3
#undef t
#undef V
	    return p;
	}
    };

    struct Threads {

	Correction operator()(Runtime &rt, const Parallel &pe,
			      const Wavefunction &wf,
			      const Map<const Array*> &V,
			      const Map<const Array*> &t) {

	    namespace ublas = boost::numeric::ublas;

	    ublas::range ro(0, wf.active().size());
	    ublas::range rv(0, wf.virtuals().size());

	    const size_t no = ro.size();
	    const size_t nv = rv.size();

	    Denominator eh(project(wf.F(), wf.active(), wf.active()));
	    Denominator ep(project(wf.F(), wf.virtuals(), wf.virtuals()));

	    struct { size_t a, b, c; } last = { -1, -1, -1 };

	    Runtime::Memory& memory = rt.memory();

	    int BLOCK = 1;
	    if (pe.rank() == 0) {
		struct name {
		    size_t operator()(size_t nb, size_t no, size_t nv) {
			size_t n1 = 3*(no*no*nv + no*no*no) + 2*(no*no*nv);
			size_t n2 = 6*(no*nv);
			size_t n3 = (no*no*no);
			return (nb*n1 + nb*nb*n2 + nb*nb*nb*n3)*sizeof(double);
		    }
		} mem;
		while (mem(BLOCK+1, no, nv) <= memory.available()) ++BLOCK;
		BLOCK = std::min<int>(BLOCK, nv);
		pe.cout() << "Block = " << BLOCK << ", "
			  << "memory needed: " << mem(BLOCK, no, nv)/1e6 << " MB"
			  << std::endl;
	    }
	    pe.broadcast(BLOCK, 0);

	    T3 t3(no, BLOCK, memory.malloc<double>(no*no*no*BLOCK*BLOCK*BLOCK));

	    Symbol< tensor_reference<4> > S;
	    // tijab
	    S.set("t(oova)", memory.malloc<double>(no*no*nv*BLOCK));
	    S.set("t(oovb)", memory.malloc<double>(no*no*nv*BLOCK));
	    S.set("t(oovc)", memory.malloc<double>(no*no*nv*BLOCK));
	    // Vijka
	    S.set("V(oooa)", memory.malloc<double>(no*no*no*BLOCK));
	    S.set("V(ooob)", memory.malloc<double>(no*no*no*BLOCK));
	    S.set("V(oooc)", memory.malloc<double>(no*no*no*BLOCK));
	    // Vijab
	    S.set("V(oovb)", memory.malloc<double>(no*no*nv*BLOCK));
	    S.set("V(oovc)", memory.malloc<double>(no*no*nv*BLOCK));
	    // Viabc
	    {
		size_t size = no*nv*BLOCK*BLOCK;
		S.set("V(ovab)", memory.malloc<double>(size));
		S.set("V(ovac)", memory.malloc<double>(size));
		S.set("V(ovba)", memory.malloc<double>(size));
		S.set("V(ovbc)", memory.malloc<double>(size));
		S.set("V(ovca)", memory.malloc<double>(size));
		S.set("V(ovcb)", memory.malloc<double>(size));
	    }

	    utility::Progress progress;
	    if (pe.rank() == 0) {
		size_t nb = (nv + BLOCK - 1)/BLOCK;
		progress.reset((nb*(nb+1)*(nb+2))/6);
	    }

	    utility::timer time;
	    utility::timer::value_type io;

	    Energy energy(*t["ia"], eh.data(), ep.data());

	    size_t ops = 0;

	    pe.task().reset();
	    size_t task = pe.task()++;

	    for (size_t c = 0, it = 0; c < nv; c += BLOCK) {
	    	for (size_t b = 0; b <= c; b += BLOCK) {
	    	    for (size_t a = 0; a <= b; a += BLOCK, ++progress) {

	    		if (it++ != task) continue;
	    		task = pe.task()++;

	    		typedef Symbol< tensor_reference<4> >::range range;
	    		range v(0,nv), o(0,no);
	    		range ra(a, std::min(a + BLOCK, nv));
	    		range rb(b, std::min(b + BLOCK, nv));
	    		range rc(c, std::min(c + BLOCK, nv));
	    		range v2(0, (nv*nv+nv)/2);

	    		utility::timer timer;

	    		if (a != last.a) {
	    		    S.load("t(oova)", *t["ijab"], o, o, v, ra);
	    		    tensor::permute<1,0,2,3>(S["t(oova)"]);
	    		    S.load("V(oooa)", *V["ijka"], o, o, o, ra);
	    		}

	    		if (b != last.b) {
	    		    S.load("t(oovb)", *t["ijab"], o, o, v, rb);
	    		    tensor::permute<1,0,2,3>(S["t(oovb)"]);
	    		    S.load("V(ooob)", *V["ijka"], o, o, o, rb);
	    		    S.load("V(oovb)", *V["ijab"], o, o, v, rb);
	    		}

	    		if (c != last.c) {
	    		    S.load("t(oovc)", *t["ijab"], o, o, v, rc);
	    		    tensor::permute<1,0,2,3>(S["t(oovc)"]);
	    		    S.load("V(oooc)", *V["ijka"], o, o, o, rc);
	    		    S.load("V(oovc)", *V["ijab"], o, o, v, rc);
	    		}

	    		// load V(iabc) as needed
	    		if (a != last.a || b != last.b) {
	    		    S.load("V(ovba)", *V["iabc"], o, v, rb, ra);
	    		    S.load("V(ovab)", *V["iabc"], o, v, ra, rb);
	    		}
	    		if (a != last.a || c != last.c) {
	    		    S.load("V(ovca)", *V["iabc"], o, v, rc, ra);
	    		    S.load("V(ovac)", *V["iabc"], o, v, ra, rc);
	    		}
	    		if (b != last.b || c != last.c) {
	    		    S.load("V(ovcb)", *V["iabc"], o, v, rc, rb);
	    		    S.load("V(ovbc)", *V["iabc"], o, v, rb, rc);
	    		}

	    		io += timer;

	    		ops += t3.evaluate(ra, rb, rc, S);
	    		t3.apply(energy, S, ra, rb, rc);

	    		last.a = a;
	    		last.b = b;
	    		last.c = c;

	    	    }
	    	}
	    }

	    Correction corr = energy.corr(pe);

	    {
		pe.cout() << "ccsd[t]: " << corr.etd << std::endl;
		pe.cout() << "ccsd(t): " << corr.ets << std::endl;
		pe.cout() <<  "memory used: " << memory.used()/1e6 << " MB" << std::endl;
		pe.cout() << "I/O time: " << io << std::endl;
		pe.cout() << "Total time: " << time << std::endl;
		pe.cout() << "GFLOP/s: " << ops/1e9 << "/" << double(time)
			  << " = " << (ops/1e9)/time << std::endl;
	    }

	    memory.clear();

	    return corr;
	}
    };

} // namespace triples
} // namespace detail
} // namespace cc
} // namespace cchem


namespace cchem {

    cc::Energy cc::triples::energy(const Parallel &pe,
				   const Wavefunction &wf,
				   const Map<const Array*> &V,
				   const Map<const Array*> &t) {

	triples::detail::Threads threads;
	Correction C = threads(Runtime::rt(), pe, wf, V, t);

	Energy E;
	E["[t]"] = C.etd;
	E["(t)"] = C.ets;    

	return E;

    }

} // namespace cchem
