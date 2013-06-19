#include "mpqc/cc/cc.hpp"
#include "mpqc/cc/tensor.hpp"
#include "mpqc/cc/utility.hpp"

// #include "mpqc/core/wavefunction.hpp"
// #include "mpqc/basis/basis.hpp"
// #include "mpqc/integrals/eri.hpp"

#include "mpqc/thread.hpp"
#include "mpqc/omp.hpp"
#include "mpqc/blas.hpp"
#include "mpqc/foreach.hpp"
#include "mpqc/utility/timer.hpp"

#include <list>
#include <boost/progress.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/adaptors.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/preprocessor/cat.hpp>

#include <boost/math/special_functions/fpclassify.hpp>

#include "boost/utility/profiler.hpp"

#ifdef CCHEM_CC_DIRECT_CUDA
#define HOST_REGISTER(ptr, size) cuda::host_register((ptr), (size)*sizeof(double));
#define HOST_UNREGISTER(ptr) cuda::host_unregister((ptr));
#else
#define HOST_REGISTER(ptr, size)
#define HOST_UNREGISTER(ptr)
#endif

namespace cchem {
namespace cc {
namespace detail {
namespace direct {


    typedef enum {
	COULOMB, EXCHANGE
    } Tag;


    struct QS {
    public:
	struct Pair {
	    int q, s;
	    size_t size;
	};
	struct Block {
	    Block() : size_(0) {}
	    const std::list<Pair>& list() const { return data_; }
	    size_t size() const { return size_; }
	private:
	    friend class QS;
	    std::list<Pair> data_;
	    size_t size_;
	};
	// split QS loop into blocks
	QS(const Basis &basis, size_t max) {
	    data_.push_back(QS::Block());
	    for (size_t s = 0; s < basis.shells().size(); ++s) {
		for (size_t q = 0; q <= s; ++q) {
		    //for (size_t q = 0; q < basis.shells().size(); ++q) {
		    Pair qs = {
			q, s,
			(basis.shells()[q].size()*
			 basis.shells()[s].size())
		    };
		    CCHEM_ASSERT(qs.size <= max);
		    if ((data_.back().size_ + qs.size) > max) {
			data_.push_back(QS::Block());
		    }
		    data_.back().data_.push_back(qs);
		    data_.back().size_ += qs.size;
		}
	    }
	}
	size_t size() const { return data_.size(); }
	const Block& operator[](int i) const { return data_.at(i); }
    private:
	std::vector<Block> data_;
    };

    static const size_t alignment = 64;
    size_t align(size_t n) {
	return alignment*((n + alignment-1)/alignment);
    }

    static void transform(const Parallel &pe,
			  const Matrix &C3, const Matrix &C4,
			  const Array &t, Array &u) {
	size_t n1 = t.shape()[0];
	size_t n2 = t.shape()[1];
	size_t n3 = C3.size1();
	Thread::Task<Parallel::Task&> task(pe.task());
	{
	    BOOST_PROFILE_LINE;
	    size_t n4 = C4.size1();
	    Matrix u0(n1*n2, n4);
	    Matrix ut(n1*n2, C4.size2());
	    task.reset();
	    while (task++ < n3) {
		int i = task;
		//std::cout << pe.rank() << " " << i << std::endl;
		size_t start[] = { 0, 0, i, 0 };
		size_t finish[] = { n1, n2, i+1, n4 };
		t.get(u0.data().begin(), start, finish);
		blas::gemm(1, u0, C4, 0, ut);
		finish[3] = ut.size2();
		u.put(ut.data().begin(), start, finish);
	    }
	    u.flush();
	    pe.barrier();
	}
	{
	    BOOST_PROFILE_LINE;
	    size_t n4 = C4.size2();
	    Matrix u0(n1*n2, n3);
	    Matrix ut(n1*n2, C3.size2());
	    task.reset();
	    while (task++ < n4) {
		int i = task;
		size_t start[] = { 0, 0, 0, i };
		size_t finish[] = { n1, n2, n3, i+1 };
		u.get(u0.data().begin(), start, finish);
		blas::gemm(1, u0, C3, 0, ut);
		finish[2] = ut.size2();
		u.put(ut.data().begin(), start, finish);
	    }
	    u.flush();
	    pe.barrier();
	}
    }


    template<size_t Order>
    static void permute(const Basis &basis, const QS::Block &block,
			tensor_reference<3> &t, double *tmp) {
	BOOST_STATIC_ASSERT((Order == 1032 || Order == 0132));
	typedef tensor_reference<3> T;
	size_t nqs = 0;
	size_t n1 = t.shape()[2];
	size_t n2 = t.shape()[1];
	foreach (QS::Pair qs, block.list()) {
	    size_t n3 = basis.shells().at(qs.q).size();
	    size_t n4 = basis.shells().at(qs.s).size();
	    T tqs(t[nqs].data(), boost::extents[n4][n3][n2*n1]);
	    T tsq(tmp, boost::extents[n3][n4][n2*n1]);
	    tensor::copy<0,2,1>(tqs, tsq);
	    T tij(tsq.data(), boost::extents[n4*n3][n2][n1]);
	    T tji(tqs.data(), boost::extents[n4*n3][n1][n2]);
	    if (Order == 0132) tensor::copy(tij, tji);
	    if (Order == 1032) tensor::copy<1,0,2>(tij, tji);
	    nqs += n3*n4;
	}
	    
    }

    template<bool Trans>
    static void put(const Basis &basis, const QS::Block &block,
		    const tensor_reference<3> &t, Array &a) {
	size_t nqs = 0;
	foreach (QS::Pair qs, block.list()) {
	    if (Trans) std::swap(qs.q, qs.s);
	    const Basis::Shell &Q = basis.shells().at(qs.q);
	    const Basis::Shell &S = basis.shells().at(qs.s);
	    size_t start[] = { 0, 0, Q.start(), S.start() };
	    size_t finish[] = { t.shape()[2], t.shape()[1],
				Q.stop(), S.stop() };
	    a.put(t[nqs].data(), start, finish);
	    nqs += (Q.size()*S.size());
	}
    }


    template<Tag tag, class P, class Q, class R, class S>
    void copy(const double *shell, boost::multi_array_ref<double,3> &a,
	      const P &p, const Q &q, const R &r, const S &s) {
	for (size_t l = 0; l < s.size(); ++l) {
	    for (size_t k = r.start(); k < r.finish(); ++k) {
		for (size_t j = q.start(); j < q.finish(); ++j) {
		    //for (size_t j = 0; j < q.size(); ++j) {
		    BOOST_AUTO(ai, (tag == COULOMB
				    ? a[k][j+l*q.size()]
				    : a[j][k+l*r.size()]));
		    std::copy(shell, shell+p.size(), &ai[p.start()]);
		    shell += p.size();
		}
	    }
	}
    }


    template<class F>
    void apply(Tag tag,
	       const Basis &basis, const QS::Block &block,
	       const integrals::Screening &screening,
	       Thread &my, F &f) {

	BOOST_PROFILE_LINE;
	
	using cc::tensor::as_matrix;
	
	namespace ublas = boost::numeric::ublas;
	using ublas::make_vector;
	using ublas::make_matrix;
	using ublas::column_major;

	const size_t N = basis.size();

	typedef tensor_reference<3> T3;
	typedef T3::extent_range range;

	integrals::Eri eri(basis, &screening);
	my.data[0] = my.realloc(my.data[0], N*block.size()*basis.max().size());
	if (my.device)
	    HOST_REGISTER(my.data[0], N*block.size()*basis.max().size());

#pragma omp for schedule(dynamic,1)
	for (int r = 0; r < (int)basis.shells().size(); ++r) {
	    BOOST_PROFILE_LINE;

	    const Basis::Shell &R = basis.shells().at(r); 

	    range rr = R.range();
	    T3 t0(my.data[0], boost::extents[rr][block.size()][N]);
	    tensor::fill(t0, 0);

	    foreach (const Basis::Block &P, basis.blocks()) {

		BOOST_PROFILE_LINE;

		size_t nqs = 0;
		foreach (QS::Pair qs, block.list()) {

		    const Basis::Shell &Q = basis.shells().at(qs.q);
		    const Basis::Shell &S = basis.shells().at(qs.s);
	
		    range rs = S.range();
		    range rq = Q.range();

		    {
			// BOOST_PROFILE_LINE;
			if (tag == COULOMB) eri(P, Q, R, S);
			if (tag == EXCHANGE) eri(P, R, Q, S);
		    }

		    // BOOST_PROFILE_LINE;
		    for (size_t i = 0; i < eri.quartets().size(); ++i) {
		    	BOOST_AUTO(const &t, eri.quartets()[i]);
			BOOST_AUTO(const &P, basis.shells().at(t[0]));
			size_t size = (P.size()*R.size()*Q.size()*S.size());
			range q(nqs, nqs+Q.size());
			range r(R.start(), R.start()+R.size());
			if (tag == COULOMB) {
			    copy<COULOMB>(eri.data() + i*size, t0, P, q, r, S);
			}
			else {
			    copy<EXCHANGE>(eri.data() + i*size, t0, P, r, q, S);
			}
		    }

		    nqs += Q.size()*S.size();
		}
	    }

	    BOOST_PROFILE_LINE;
	    t0.reindex(0);
	    f(my, R, block, t0);

	} // for R

	if (my.device)
	    HOST_UNREGISTER(my.data[0]);

    }


    struct Transform {
	typedef QS::Block Block;
	Transform(size_t N, size_t no,
		  boost::reference_wrapper<const Matrix> Co,
		  boost::reference_wrapper<const Matrix> t1,
		  boost::reference_wrapper<const Array> t2,
		  boost::reference_wrapper<Symbol< tensor_reference<3> > > S)
	    : rn_(0,N), ro_(0,no),
	      Co_(Co), t1_(t1), t2_(t2), S_(S) {}
	template<class R>
	void operator()(Thread &my,
			const R &r, const Block &block,
			const tensor_reference<3> &t0);
    private:
	boost::numeric::ublas::range rn_, ro_;
	boost::reference_wrapper<const Matrix> Co_;
	boost::reference_wrapper<const Matrix> t1_;
	boost::reference_wrapper<const Array> t2_;
	boost::reference_wrapper<Symbol< tensor_reference<3> > > S_;

	template<class C1, class C2, class Context>
	static void transform1(const C1 &c1, const C2 &c2,
			       tensor_const_reference<2> t0,
			       tensor_reference<2> t1,
			       tensor_reference<3> vt1,
			       const Context &context) {
	    BOOST_PROFILE_LINE;
	    using cc::tensor::as_matrix;
	    using cc::tensor::as_vector;
	    BOOST_AUTO(t, (as_matrix<1,1>(t1)));
	    // V(p[rs]q)*c1(pi) -> x(i[rs]q)
	    //t.clear();
	    blas::gemm(1, c1, (as_matrix<1,1>(t0)), 0, t, context);
	    // x(i[rs]q)*c2(qj) -> x(ji[rs])
	    BOOST_AUTO(vt, (as_matrix<1,2>(vt1)));
	    blas::ger(1, c2, (as_vector(t)), vt, context);
	}


	template<class R, class C, class Context>
	static void transform1(const R &r, const C &C1, const C &C2,
			       tensor_const_reference<3> t0,
			       tensor_reference<2> t1,
			       tensor_reference<3> vt1,
			       const Context &context) {
	    BOOST_PROFILE_LINE;
	    for (size_t f = 0; f < r.size(); ++f) {
		transform1(C1, column(C2,f+r.start()), t0[f], t1, vt1, context);
	    }
	}

	template<class R, class C>
	static void transform1(Thread &my, const R &r,
			       const C &C1, const C &C2,
			       const tensor_reference<3> &t0,
			       tensor_reference<3> vt1);

    };

    template<class R, class C>
    void Transform::transform1(Thread &my, const R &r,
			       const C &C1, const C &C2,
			       const tensor_reference<3> &t0,
			       tensor_reference<3> vt1) {

	using cc::tensor::as_matrix;
	using cc::tensor::as_vector;

	size_t m1 = C1.size1();
	size_t m2 = t0.shape()[1];
	my.data[1] = my.realloc(my.data[1], m1*m2);
	tensor_reference<2> t1(my.data[1], boost::extents[m2][m1]);

	size_t n1 = C2.size1();
	my.data[2] = my.realloc(my.data[2], n1*m1*m2);
	tensor_reference<3> ut1(my.data[2], boost::extents[m2][m1][n1]);
	tensor::fill(ut1, 0);

	transform1(r, C1, C2, t0, t1, ut1, blas::Context());
#pragma omp critical(BOOST_PP_CAT(cc_direct_, __LINE__))
	as_vector(vt1).plus_assign(as_vector(ut1));
    }

    template<class R>
    void Transform::operator()(Thread &my,
			       const R &r, const Block &block,
			       const tensor_reference<3> &t0) {

	BOOST_PROFILE_LINE;

	namespace ublas = boost::numeric::ublas;
	using cc::tensor::as_matrix;
	using cc::tensor::as_vector;

	size_t N = rn_.size();
	BOOST_AUTO(const &ro, ro_);

	BOOST_AUTO(const &Co, Co_.get());
	BOOST_AUTO(const &t1, t1_.get());
	BOOST_AUTO(const &t2, t2_.get());

#ifdef CCHEM_CC_DIRECT_CUDA
	if (my.device) {
	    BOOST_PROFILE_LINE;
	    BOOST_AUTO(&device, *my.device);
	    cuda::synchronize();

	    device.S.reset("t1", boost::extents[1][1][t1.size2()][t1.size1()]);
	    cublas::set_vector(t1.data().size(),
			       t1.data().begin(),
			       device.S["t1"].data());

	    device.S.reset("Co", boost::extents[1][1][Co.size2()][Co.size1()]);
	    cublas::set_vector(Co.data().size(),
			       Co.data().begin(),
			       device.S["Co"].data());
	}
#endif

	{
	    BOOST_PROFILE_LINE;

#ifdef CCHEM_CC_DIRECT_CUDA
	    if (my.device) {
		BOOST_AUTO(&device, *my.device);
		BOOST_AUTO(&streams, my.device->streams);
		size_t ns = streams.size();

		BOOST_AUTO(const &t1, (as_matrix<1,1>(device.S["t1"][0][0])));
		BOOST_AUTO(const &Co, (as_matrix<1,1>(device.S["Co"][0][0])));

		size_t m1 = Co.size1();
		size_t m2 = t0.shape()[1];

		foreach (Device::Stream &stream, streams) {
		    stream.S.reset("t0", boost::extents[1][1][m2][N]);
		    stream.S.reset("t'", boost::extents[1][1][m2][m1]);
		}

		for (size_t f = 0; f < r.size(); ++f) {
		    BOOST_AUTO(stream, streams.at(f%ns));
		    stream.synchronize();
		    BOOST_AUTO(t, (stream.S["t'"][0][0]));
		    BOOST_AUTO(g, stream.S["t0"][0][0]);
		    cublas::set_vector(t0[f].num_elements(),
				       t0[f].data(), g.data(),
				       stream);
		    stream.wait(device.events.at((f+1)%ns));
		    cublas::set_stream(stream);
		    transform1(Co, column(t1,f+r.start()), g,
			       t, device.S["vt1\'"][0],
			       device.cublas_handle);
		    transform1(t1, column(Co,f+r.start()), g,
			       t, device.S["vt1\""][0],
			       device.cublas_handle);
		    if (device.S.contains("vt1")) {
			transform1(t1, column(t1,f+r.start()), g,
				   t, device.S["vt1"][0],
				   device.cublas_handle);
		    }
		    device.events.at(f%ns).record(stream);
		}
		cuda::synchronize();

	    }
#endif

	    if (!my.device) {
		BOOST_PROFILE_LINE;
		transform1(my, r, Co, t1, t0, S_.get()["vt1\'"]);
		transform1(my, r, t1, Co, t0, S_.get()["vt1\""]);
		if (S_.get().contains("vt1")) {
		    transform1(my, r, t1, t1, t0, S_.get()["vt1"]);
		}
	    }

	}

	if (!S_.get().contains("vt2")) {
	    return;
	}

	size_t no = ro.size();
	size_t n2 = t0.shape()[1];
	my.data[1] = my.realloc(my.data[1], no*no*n2);
	tensor_reference<3> vt2(my.data[1], boost::extents[n2][no][no]);
	tensor::fill(vt2, 0);

	{
	    BOOST_PROFILE_LINE;

	    my.data[2] = my.realloc(my.data[2], no*no*N);
	    my.S.set("t", my.data[2], boost::extents[1][no][no][N]);
	    if (my.device)
		HOST_REGISTER(my.S["t"].data(), my.S["t"].num_elements());

	    for (size_t f = 0; f < r.size(); ++f) {

		{
		    // BOOST_PROFILE_LINE;
		    my.S.load("t", t2, ro, ro, ublas::range(0,N), f+r.start());
		}

#ifdef CCHEM_CC_DIRECT_CUDA
		if (my.device) {
		    BOOST_PROFILE_LINE;

		    BOOST_AUTO(&device, *my.device);
		    size_t ns = device.streams.size();
		    BOOST_AUTO(&stream, device.streams.at(f%ns));

		    stream.S.reset("t'", boost::extents[1][N][no][no]);
		    BOOST_AUTO(t2, stream.S["t'"][0]);
		    stream.S.reset("t0", boost::extents[1][1][n2][N]);
		    BOOST_AUTO(g0, stream.S["t0"][0][0]);

		    {
			// BOOST_PROFILE_LINE;
			stream.synchronize();
			// BOOST_PROFILE_LINE;
			cublas::set_vector(my.S["t"].num_elements(),
					   my.S["t"].data(), t2.data());
			cublas::set_vector(t0[f].num_elements(),
					   t0[f].data(), g0.data(),
					   stream);
		    }

		    // ensure vt have been updated by previous stream
		    stream.wait(device.events.at((f+1)%ns));
		    cublas::set_stream(stream);
		    {
			BOOST_AUTO(&device, *my.device);
			BOOST_AUTO(vt, (as_matrix<2,2>(device.S["vt2"])));
			blas::gemm(1,
				   (as_matrix<2,1>(t2)),
				   (as_matrix<1,1>(g0)),
				   1, vt,
				   device.cublas_handle);
		    }
		    device.events.at(f%ns).record(stream);

		    continue;

		}
#endif

		{
		    BOOST_PROFILE_LINE;
		    BOOST_AUTO(t, (as_matrix<2,1>(my.S["t"][0])));
		    BOOST_AUTO(g, (as_matrix<1,1>(t0[f])));	
		    BOOST_AUTO(vt, (as_matrix<2,1>(vt2)));
		    blas::gemm(1, t, g, 1, vt);
		}

	    }

	}

#ifdef CCHEM_CC_DIRECT_CUDA
	if (my.device) {
	    BOOST_PROFILE_LINE;
	    cuda::synchronize();
	    HOST_UNREGISTER(my.S["t"].data());
	    return;
	}
#endif

#pragma omp critical(BOOST_PP_CAT(cc_direct_, __LINE__))
	S_.get()["vt2"].plus_assign(vt2);

    }



    void evaluate(Runtime &rt, const Parallel &pe,
		  const Wavefunction &wf,
		  const integrals::Screening &screening,
		  const Matrix &t1, const Array &t2,
		  const Map<Array*> &A,
		  size_t MAXQS, Tag tag) {

	BOOST_PROFILE_LINE;

	const Basis &basis = wf.basis();
	Matrix Co = trans(wf.C(wf.active())); // C(p,a)
	size_t N = basis.size();
	size_t no = Co.size1();

	typedef QS::Block Block;
	QS blocks(basis, MAXQS);

	Runtime::Memory &memory = rt.memory();
	std::map<int, double*> buffer;

	buffer[0] = memory.malloc<double>(MAXQS*no*no);
	buffer[1] = memory.malloc<double>(MAXQS*no*no);

	Symbol< tensor_reference<3> > S;
	S.set("vt1\'", buffer[0], boost::extents[0][0][0]);
	S.set("vt1\"", buffer[1], boost::extents[0][0][0]);

	if (A.has("vt1")) {
	    buffer[-1] = memory.malloc<double>(MAXQS*no*no);
	    S.set("vt1", buffer[-1], boost::extents[0][0][0]);
	}

	if (A.has("vt2")) {
	    buffer[2] = memory.malloc<double>(MAXQS*no*no);
	    S.set("vt2", buffer[2], boost::extents[0][0][0]);
	}

	Thread::Task<Parallel::Task&> task(pe.task());
	task.reset();
	task.index = 0;
	++task;

#pragma omp parallel
	{
	    Thread thread;

#ifdef CCHEM_CC_DIRECT_CUDA
	    std::vector<int> devices = rt.devices();
#pragma omp critical
	    if (devices.size()) {
		try {
		    BOOST_PROFILE_LINE;
		    cuda::set_device(devices.at(omp::thread()%devices.size()));
		    thread.device = new Device(2);
		    BOOST_AUTO(&device, *thread.device);
		    foreach (Device::Stream &s, device.streams) {
		    	s.S.set("t0", device.malloc(MAXQS*N));
		    	size_t n1 = no*MAXQS;
		    	size_t n2 = N*no*no;
		    	s.S.set("t'", device.malloc(std::max(n1, n2)));
		    }
		    foreach (const std::string &s, S.keys()) {
		    	device.S.set(s, device.malloc(MAXQS*no*no));
		    }
		    device.S.set("t1", device.malloc(t1.data().size()));
		    device.S.set("Co", device.malloc(Co.data().size()));
		}
		// failed to allocate device/resources
		catch (const cuda::error &e) {
		    delete thread.device;
		    thread.device = NULL;
		}
		
	    }
#endif // CCHEM_CC_DIRECT_CUDA

	    for (int qs = 0; qs < blocks.size(); ++qs) {

		if (qs != task) continue;
#pragma omp barrier
#pragma omp master
		++task;
#pragma omp barrier

		const Block &block = blocks[qs];
		size_t M = block.size();

#pragma omp barrier
#pragma omp master
		foreach (const std::string &s, S.keys()) {
		    S.reset(s, boost::extents[M][no][no]);
		    tensor::fill(S[s], 0);
		}
#pragma omp barrier

#ifdef CCHEM_CC_DIRECT_CUDA
		if (thread.device) {
		    BOOST_AUTO(&device, *thread.device);
		    foreach (const std::string &s, S.keys()) {
			device.S.reset(s, boost::extents[1][M][no][no]);
			cublas::clear(device.S[s].num_elements(),
				      device.S[s].data());
		    }
		}
#endif // CCHEM_CC_DIRECT_CUDA

		// apply transforms
		{
		    BOOST_PROFILE_LINE;
		    using boost::ref;
		    using boost::cref;
		    Transform f(N, no, cref(Co), cref(t1), cref(t2), ref(S));
		    apply(tag, basis, block, screening, thread, f);
		}
#pragma omp barrier


#ifdef CCHEM_CC_DIRECT_CUDA
		if (thread.device) {
		    BOOST_PROFILE_LINE;
		    using tensor::as_vector;
		    cuda::synchronize();
		    BOOST_AUTO(&device, *thread.device);
		    thread.data[0] = thread.realloc(thread.data[0], M*no*no);
		    tensor_reference<4>
			host(thread.data[0], boost::extents[1][M][no][no]);
		    foreach (const std::string &s, S.keys()) {
		    	cublas::get_vector(device.S[s].num_elements(),
		    			   device.S[s].data(), host.data());
#pragma omp critical
			as_vector(S[s]).plus_assign(as_vector(host));
		    }
		}
#pragma omp barrier
#endif

#pragma omp master
		{
		    BOOST_PROFILE_LINE;

		    put<0>(basis, block, S["vt1\'"], *A["vt1'"]);

		    // use vt1' as temp buffer
		    double *tmp = S["vt1\'"].data();

		    permute<1032>(basis, block, S["vt1\""], tmp);
		    put<1>(basis, block, S["vt1\""], *A["vt1'"]);

		    if (A.has("vt1")) {
			tensor::permute<1,0,2>(S["vt1"]);
			put<0>(basis, block, S["vt1"], *A["vt1"]);
			permute<1032>(basis, block, S["vt1"], tmp);
			put<1>(basis, block, S["vt1"], *A["vt1"]);
		    }

		    if (A.has("vt2")) {
			put<0>(basis, block, S["vt2"], *A["vt2"]);
			permute<1032>(basis, block, S["vt2"], tmp);
			put<1>(basis, block, S["vt2"], *A["vt2"]);
		    }

		}
	    }
	    
	} // parallel

	BOOST_PROFILE_LINE;

	// free memory
	{
	    BOOST_AUTO(it, buffer.begin());
	    while (it != buffer.end()) {
		memory.free(it->second);
		++it;
	    }
	}

	A["vt1'"]->flush();
	if (A.has("vt1")) A["vt1"]->flush();
	if (A.has("vt2")) A["vt2"]->flush();
	pe.barrier();
	    
    }

    /* U needs to reside in fast memory for performance */
    void transform34(const Parallel &pe,
		     const Wavefunction &wf,
		     const Matrix &C3, const Matrix &C4,
		     Array &A, Array &u) {

	namespace ublas = boost::numeric::ublas;
	using tensor::as_matrix;

	size_t N = wf.basis().size();
	ublas::range r3(0, C3.size2());
	ublas::range r4(0, C4.size2());
	ublas::range ro(0, wf.active().size());

	Thread::Task<Parallel::Task&> task(pe.task());

	Symbol< tensor_reference<4> > S;
	Runtime::Memory &memory = Runtime::rt().memory();
	double* buffer[2] = {};
	buffer[0] = memory.malloc<double>
	    (ro.size()*ro.size()*N);
	buffer[1] = memory.malloc<double>
	    (ro.size()*ro.size()*std::max(r3.size(), r4.size()));

	task.reset();
	while (++task < N) {
	    int r = task;
	    BOOST_PROFILE_LINE;
	    S.load("U", buffer[0], ro, ro, ublas::range(0,N), r, A);
	    S.set("t", buffer[1], ro, ro, r3, 0);
	    BOOST_AUTO(t, (as_matrix<2,1>(S["t"][0])));
	    blas::gemm(1, as_matrix<2,2>(S["U"]), C3, 0, t);
	    S.store("t", ro, ro, r3, r, u);
	}
	u.flush();
	pe.barrier();

	task.reset();
	while (++task < r3.size()) {
	    int b = task;
	    BOOST_PROFILE_LINE;
	    S.load("U", buffer[0], ro, ro, b, ublas::range(0,N), u);
	    S.set("t", buffer[1], ro, ro, r4, 0);
	    BOOST_AUTO(t, (as_matrix<2,1>(S["t"][0])));
	    blas::gemm(1, as_matrix<2,2>(S["U"]), C4, 0, t);
	    S.store("t", ro, ro, b, r4, u);
	}
	u.flush();
	pe.barrier();

	// transfer from U to final A
	task.reset();
	while (++task < r4.size()) {
	    int b = task;
	    S.load("U", buffer[0], ro, ro, r3, b, u);
	    S.store("U", ro, ro, r3, b, A);
	}
	A.flush();
	pe.barrier();

	memory.free(buffer[0]);
	memory.free(buffer[1]);
    }

} // namespace direct
} // namespace detail
} // namespace cc
} // namespace cchem

/** evaluate vt1, vt2 diagrams */
void cchem::cc::sd::direct(const Parallel &pe,
			   const Wavefunction &wf,
			   const Matrix &t1, const Array &t2,
			   Map<Matrix> &O1, const Map<Array*> &O2,
			   const ::integrals::Screening &screening) {

    BOOST_PROFILE_REGISTER_THREAD;

    struct {
	utility::timer total;
	utility::timer::value_type eri, trans[4], ab;
	utility::timer::value_type device;
    } time;

    namespace ublas = boost::numeric::ublas;

    using tensor::as_matrix;
    using tensor::as_vector;

    const Basis &basis = wf.basis();

    size_t N = basis.size();
    size_t no = wf.active().size();
    size_t nv = wf.virtuals().size();
  
    Symbol< tensor_reference<4> > S;	
    ublas::range ro(0, no), rv(0, nv);

    Runtime &rt = Runtime::rt();

    Matrix Cv = wf.C(wf.virtuals()); // C(p,a)
    Matrix Co = wf.C(wf.active());
    Matrix C(N,nv+no);
    subrange(C, 0,N, 0,nv) = Cv;
    subrange(C, 0,N, nv,nv+no) = Co;

    size_t MAXQS = 0;
    {
	size_t available = rt.memory().available();
	size_t threads = omp_get_max_threads();
	size_t max = basis.max().size();
	struct {
	    size_t operator()(size_t no, size_t N, size_t max,
			      size_t block, size_t threads) const {
		size_t memory = 0;
		memory += N*block*max;
		memory += no*no*std::max(N,block);
		memory += no*no*block;
		memory *= threads;
		memory += 4*no*no*block;
		return memory*sizeof(double);
	    }
	} memory;
	size_t nodes = pe.size();
	while (memory(no, N, max, MAXQS+1, threads) <= available) {
	    if (MAXQS >= (N*N+nodes-1)/nodes) break;
	    ++MAXQS;
	}
	CCHEM_ASSERT(MAXQS >= max*max);
    }

    {
	pe.cout() <<  "direct vt terms: block = " << MAXQS << std::endl;
    }

    {
	namespace Direct = cc::detail::direct;
	Matrix u1 = prod(t1,trans(Cv));
	Array &u2 = *O2["u(ijab)"];
	Direct::transform(pe, trans(Cv), trans(Cv), t2, u2);
	BOOST_PROFILE_DUMP(std::cout);
	Map<Array*> A;
	A["vt1'"] = O2["vt1\""];
	Direct::evaluate(rt, pe, wf, screening, u1, u2, A, MAXQS, Direct::EXCHANGE);
	BOOST_PROFILE_DUMP(std::cout);
	A["vt2"] = O2["vt2"];
	A["vt1"] = O2["vt1"];
	A["vt1'"] = O2["vt1\'"];
	Direct::evaluate(rt, pe, wf, screening, u1, u2, A, MAXQS, Direct::COULOMB);
	BOOST_PROFILE_DUMP(std::cout);
    }
    pe.barrier();

    // std::cout << "    eri: " << time.eri << std::endl;
    // std::cout << "    trans 1: " << time.trans[0] << std::endl;
    // std::cout << "    trans 2: " << time.trans[1] << std::endl;
    // std::cout << "    t(i,j,q,s): " << time.ab << std::endl;
    // std::cout << "    device: " << time.device << std::endl;

    {
    	BOOST_PROFILE_LINE;
	Array &u = *O2["u(ijab)"];
    	cc::detail::direct::transform34(pe, wf, C, C, *O2["vt2"], u);
    	cc::detail::direct::transform34(pe, wf, C, C, *O2["vt1"], u);
    	cc::detail::direct::transform34(pe, wf, C, C, *O2["vt1\'"], u);
    	cc::detail::direct::transform34(pe, wf, C, C, *O2["vt1\""], u);
    }
    pe.barrier();

    BOOST_PROFILE_DUMP(std::cout);

}


