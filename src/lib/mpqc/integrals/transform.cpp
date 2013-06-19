// These headers are first - otherwise bizarre errors ith gcc < 4.3 and IMPI
//#include "runtime.hpp"
#include "mpqc/omp.hpp"
#include "mpqc/blas.hpp"

#include "mpqc/cc/cc.hpp"
//#include "mpqc/cc/tensor.hpp"
#include "mpqc/cc/utility.hpp"

//#include "mpqc/thread.hpp"
#include "mpqc/core/wavefunction.hpp"
#include "mpqc/basis/basis.hpp"
#include "mpqc/integrals/eri.hpp"
#include "mpqc/foreach.hpp"
#include "mpqc/utility/timer.hpp"
#include "mpqc/utility/progress.hpp"

#include <boost/bind.hpp>

#include "boost/utility/profiler.hpp"

// namespace mpqc {
// namespace cc {
// namespace detail {

//     template<class R>
//     void copy(const double *shell, boost::multi_array_ref<double,4> &a,
// 	      const R &p, const R &q, const R &r, const R &s) {
// 	for (size_t l = s.start(); l < s.stop(); ++l) {
// 	    for (size_t k = r.start(); k < r.stop(); ++k) {
// 		for (size_t j = q.start(); j < q.stop(); ++j) {
// 		    std::copy(shell, shell+p.size(), &a[j][k][l][p.start()]);
// 		    shell += p.size();
// 		}
// 	    }
// 	}
//     }


//     template<class C, class T>
//     void transform12(const Basis &basis,
// 		     const Basis::Shell &R, const Basis::Block &S,
// 		     integrals::Eri eri, const C &C1, const C &C2,
// 		     T &t2, omp::lock &lock) {

// 	namespace ublas = boost::numeric::ublas;
// 	using matrix_cast;
// 	typedef tensor_reference<4> T4;
// 	typedef typename T4::extent_range range;
// 	using boost::extents;

// 	const size_t MAX = std::min<size_t>(basis.max_block(), 128);
// 	BOOST_AUTO(const &blocks, basis.blocks(MAX));

// 	range rs = S.range();
// 	range rr = R.range();

// 	Thread thread;
// 	thread.data[0] = thread.malloc(MAX*(MAX*rr.size()*rs.size()));
// 	thread.data[1] = thread.malloc(C1.size2()*(MAX*rr.size()*rs.size()));

// #pragma omp for schedule(dynamic,1)
// 	for (int q = 0; q < blocks.size(); ++q) {
// 	    BOOST_PROFILE_LINE;
// 	    utility::timer timer;

// 	    const Basis::Block &Q = blocks.at(q);

// 	    range rq = Q.range();
// 	    T4 t1(thread.data[1], extents[rq.size()][rr.size()][rs.size()][C1.size2()]);
// 	    std::fill_n(t1.data(), t1.num_elements(), 0);

// 	    foreach (const Basis::Block &P, blocks) {
// 		range rp = P.range();

// 		T4 G(thread.data[0], extents[rq][rr][rs][rp]);
// 		std::fill_n(G.data(), G.num_elements(), 0);

// 		eri(P, Q, R, S);
// 		size_t size = (P.shell().size()*Q.shell().size()*
// 			       R.size()*S.shell().size());
// 		for (size_t i = 0; i < eri.quartets().size(); ++i) {
// 		    BOOST_AUTO(const &q, eri.quartets()[i]);
// 		    detail::copy(eri.data() + i*size, G,
// 				 basis[q[0]], basis[q[1]],
// 				 basis[q[2]], basis[q[3]]);
// 		}
// 		ublas::range r1(rp.start(), rp.finish());
// 		ublas::range r2(0, C1.size2());
// 		BOOST_AUTO(t, (as_matrix<1,3>(t1)));
// 		blas::gemm(1,
// 			   project(trans(C1), r2, r1),
// 			   (as_matrix<1,3>(G)),
// 			   1, t);
// 	    }

// 	    // 2nd transformation
// 	    {
// 		BOOST_PROFILE_LINE;
// 		BOOST_AUTO(t, (as_matrix<3,1>(t2)));
// 		ublas::range r1(rq.start(), rq.finish());
// 		ublas::range r2(0, C2.size2());
// 		lock.set();
// 		blas::gemm(1,
// 			   matrix_cast<3,1>(t1),
// 			   project(C2, r1, r2),
// 			   1, t);
// 		lock.unset();
// 	    }
// 	} // for q

//     }


//     // (O,N,s)
//     template<class R, class C, class T2>
//     void transform3(int q, R rs, const C &c, const T2 &t2,
// 		    Array &v, double *tmp) {
// 	BOOST_PROFILE_LINE;
// 	tensor_reference<3>
// 	    t3(tmp, boost::extents[rs.size()][c.size2()][t2.shape()[2]]);
// 	for (int i = 0; i < rs.size(); ++i) {
// 	    using matrix_cast;
// 	    BOOST_AUTO(t, (as_matrix<1,1>(t3[i])));
// 	    blas::gemm(1, (as_matrix<1,1>(t2[i])), c, 0, t);
// 	}
// 	size_t start[4] = { 0, 0, q, rs.start() };
// 	size_t finish[4] = { t3.shape()[2], t3.shape()[1],
// 			     q+1, rs.start()+rs.size() };
// 	BOOST_PROFILE_LINE;
// 	v.put(t3.data(), start, finish);
//     }


//     template<class U, class C, class T>
//     void transform(const U &u, const C &c, T &t) {
// 	// blas::gemm(1, u, c, 0, t);
// 	const int B = 128;
// 	using boost::numeric::ublas::range;
// 	int m = u.rows();
// #pragma omp parallel for schedule(dynamic,1)
// 	for (int i = 0; i < m; i += B) {
// 	    range K(0, u.size2());
// 	    range M(i, std::min(i+B, m));
// 	    range N(0, c.size2());
// 	    BOOST_AUTO(ti, project(t, M, N));
// 	    blas::gemm(1, project(u, M, K), c, 0, ti);
// 	}
//     }


//     template<class T, class C>
//     void transform4(int r, const C &c, const T &t, Array &V, double *buffer) {
// 	BOOST_PROFILE_LINE;
// 	{
// 	    using matrix_cast;
// 	    using boost::numeric::ublas::make_matrix;
// 	    using boost::numeric::ublas::column_major;
// 	    int m = t.shape()[3]*t.shape()[2];
// 	    BOOST_AUTO(t4, make_matrix<column_major>(m, c.size2(), buffer));
// 	    transform(as_matrix<2,2>(t), c, t4);
// 	}
// 	size_t start[] = { 0, 0, r, 0 };
// 	size_t finish[] = { t.shape()[3], t.shape()[2], r+1, c.size2() };
// 	V.put(buffer, start, finish);
//     }


// }
// }
// } // namespace mpqc


void mpqc::integrals::transform(const Comm &comm,
                                const Wavefunction &wf,
                                Array::Map<double> v,
                                Array::Array<double> *v2 = NULL) {

    BOOST_PROFILE_REGISTER_THREAD;

    const Basis &basis = wf.basis();
    size_t N = basis.size();
    size_t no = wf.active().size();
    size_t nv = wf.virtuals().size();
    Matrix Co = wf.C(wf.active()); // C(p,i)
    Matrix Cv = wf.C(wf.virtuals()); // C(p,a)

    const size_t MAXS = basis.max().size();
    const size_t MAXR = basis.max().size();
    const size_t MAXP = std::min<size_t>(basis.max_block(), 128);
    const size_t MAXQ = MAXP;

    typedef mpqc::TensorBase<double,3> T3;
    typedef mpqc::TensorBase<double,4> T4;

    const bool IABC = v.has("iabc");

    typedef Symbol< tensor_reference<4> >::range range;
    range ro(0,no), rv(0,nv), ra(0,N);
    Matrix C(N,no+nv);
    C(ra,rv) = Cv;
    C(ra,ro) = Co;

    struct {
	std::vector<Basis::Block> P, Q, R, S;
    } blocks;

    // block shells into segments
    blocks.P = basis.blocks(MAXP);
    blocks.Q = basis.blocks(MAXQ);
    blocks.S = basis.blocks(MAXS);

    int max_threads = omp_get_max_threads();

    {
	comm.cout << "Integral transformation: ";
	comm.cout << ((IABC) ? "(iabc)" : "(ija*), (iajb)") << std::endl;
	comm.cout << "Threads: " << max_threads << std::endl;
    }

    utility::Progress progress;
    if (comm.rank() == 0) progress.reset(basis.size());

    Comm::Task task(comm);

    foreach (const Basis::Block &S, blocks.S) {

	BOOST_PROFILE_LINE;

	comm.barrier();
        task.reset();
	omp_set_num_threads(max_threads);

	range rs = S.range();

#pragma omp parallel
	{

            double *buffer[4] = { 0 }:
            buffer[2] = memory.malloc<double>(no*C.cols()*MAXR*MAXS);

            for (int r = 0; r < basis.shells().size(); ++r) {
                if (r != task++) continue;

		const auto R = task.next(basis.shells());
                if (r == basis.shells().end())

                range ri = no;
                range rq = C.cols();
		range rr = R.range();
                range rs = S.range();

		T4 t2(thread.buffer[2], range::dimensions(ri, rq, rr, rs));
                t2.zero();

// 		omp::lock lock;
// #pragma omp parallel
// 		{
// 		    integrals::Eri eri(basis, &screening);
// 		    detail::transform12(basis, R, S, eri, Co, C, t2, lock);
// 		} // parallel

                v2(ri,rq,rr,rs) << t2;

	    }

            memory.free(buffer[2]);

	} // parallel

	progress += rs.size();

    } // foreach S
	
}


