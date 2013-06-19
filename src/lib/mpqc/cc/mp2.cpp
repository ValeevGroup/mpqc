#include "mpqc/cc/cc.hpp"
#include "mpqc/cc/tensor.hpp"
#include "mpqc/cc/symmetrize.hpp"
#include "mpqc/cc/utility.hpp"

#include "mpqc/thread.hpp"
#include "mpqc/parallel/environment.hpp"
#include "mpqc/blas.hpp"
#include "mpqc/foreach.hpp"
#include "mpqc/utility/timer.hpp"

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/adaptors.hpp>

#include <stdio.h>

double cchem::cc::mp2(const Parallel &pe, const Wavefunction &wf,
		      const Array &V, Array &t2) {

    namespace ublas = boost::numeric::ublas;
    using tensor::as_matrix;
    using tensor::as_vector;

    ublas::range ro(0, wf.active().size());
    ublas::range rv(0, wf.virtuals().size());

    const int no = ro.size();
    const int nv = rv.size();

    Symbol< tensor_reference<4> > S;
    Buffer<double> buffer;
    buffer[0].resize(no*no*nv);
    buffer[1].resize(no*no*nv);

    Denominator ep(ublas::project(wf.F(), wf.virtuals(), wf.virtuals()));
    Denominator eh(ublas::project(wf.F(), wf.active(), wf.active()));

    double E = 0;

    Thread::Task<Parallel::Task&> task(pe.task());
    task.reset();

    while (task++ < nv) {
	int b = task;
	S.load("t(ijab)", buffer[0], ro, ro, rv, b, V);
#pragma omp parallel for reduction(+:E)
	for (int a = 0; a < nv; ++a) {
	    BOOST_AUTO(t, (as_matrix<1,1>(S["t(ijab)"][0][a])));
	    for (int j = 0; j < no; ++j) {
		for (int i = 0; i <= j; ++i) {
		    double e = (eh(i,j) - ep(a,b));
		    //assert(fabs(e) > 1e-100);
		    double tij = t(i,j);
		    double tji = t(j,i);
		    int sym = 2 - (i == j);
		    E += sym*tij*(2*tij - tji)/e;
		    t(i,j) = tij/e;
		    t(j,i) = tji/e;
		}
	    }
	}
	// std::cout <<  task << ":" << E << std::endl;
	S.store("t(ijab)", ro, ro, rv, b, t2);
	// S["t(ijab)"][0].permute<1,0,2>();
	// S.store("t(ijab)", ro, ro, b, rv, t2);
    }

    t2.flush();
    pe.reduce("+", &E, 1);
    pe.barrier();

    return E;

}

