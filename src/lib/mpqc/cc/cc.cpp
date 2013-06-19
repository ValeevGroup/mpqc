#include "mpqc/cc/cc.hpp"
#include "mpqc/utility/array.hpp"
#include "mpqc/cc/diis.hpp"



cchem::cc::Energy cchem::cc::energy(Wavefunction wf, Runtime &rt,
				    const std::string &method) {

    using utility::make_array;

    wf.sort();
    wf.reverse();

    double cutoff = rt.get<double>("/cc/integrals/cutoff", 1e-10);
    ::integrals::Screening screening(wf.basis(), cutoff);

    size_t N = wf.basis().size();
    size_t no = wf.active().size();
    size_t nv = wf.virtuals().size();

    // std::cout <<  "atomic orbitals: " << N << std::endl;
    // std::cout <<  "occupied orbitals: " << no << std::endl;
    // std::cout <<  "virtual orbitals: " << nv << std::endl;

    Parallel pe;
    
#define RT_ARRAYS_ALLOCATE(name, dims)					\
    if (!rt.arrays().contains(name)) {					\
	rt.arrays().allocate<double>(name, make_array dims, pe);	\
	pe.cout() << *rt.arrays().find<Array>(name) << std::endl;	\
    }

    // first ones may be allocated in faster storage
    RT_ARRAYS_ALLOCATE("cc.t(ijab)", (no,no,nv,nv));
    RT_ARRAYS_ALLOCATE("cc.u(ijab)", (no,no,N,N));
    RT_ARRAYS_ALLOCATE("integrals.v(ijab)", (no,no,nv+no,N));
    
    size_t B = wf.basis().max().size();
    RT_ARRAYS_ALLOCATE("integrals.v(iqrs)", (no,B,N,no+nv));

    RT_ARRAYS_ALLOCATE("integrals.v(ijka)", (no,no,no,nv));
    RT_ARRAYS_ALLOCATE("integrals.v(ijkl)", (no,no,no,no));
    RT_ARRAYS_ALLOCATE("integrals.v(iajb)", (no,nv,no,N));

    {
    	Map<Array*> V;
    	V["ijkl"] = rt.arrays().find<Array>("integrals.v(ijkl)");
    	V["ijka"] = rt.arrays().find<Array>("integrals.v(ijka)");
    	V["ijab"] = rt.arrays().find<Array>("integrals.v(ijab)");
    	V["iajb"] = rt.arrays().find<Array>("integrals.v(iajb)");
    	// V["iabc"] = rt.arrays().find<Array>("integrals.v(iabc)");
    	V["iqrs"] = rt.arrays().find<Array>("integrals.v(iqrs)");
    	utility::timer timer;
    	double cutoff = rt.get<double>("/cc/integrals/cutoff", 1e-10);
    	integrals::Screening screening(wf.basis(), cutoff);
    	cc::integrals(pe, wf, V, screening);
    	pe.cout() << "integrals time: " << timer << std::endl;
    }

    if (method == "ccsd")
	rt.arrays().erase<Array>("integrals.v(iqrs)");

    Map<const Array*> V;
    V["ijkl"] = rt.arrays().find<Array>("integrals.v(ijkl)");
    V["ijka"] = rt.arrays().find<Array>("integrals.v(ijka)");
    V["ijab"] = rt.arrays().find<Array>("integrals.v(ijab)");
    V["iajb"] = rt.arrays().find<Array>("integrals.v(iajb)");
    // V["iabc"] = rt.arrays().find<Array>("integrals.v(iabc)");

    Map<Array*> A;
    A["t(ijab)"] = rt.arrays().find<Array>("cc.t(ijab)");
    A["u(ijab)"] = rt.arrays().find<Array>("cc.u(ijab)");

    cc::Energy E;
    E["mp2"] = mp2(pe, wf, *V["ijab"], *A["t(ijab)"]);
    pe.cout() << "mbpt(2) energy: "
       << std::setprecision(10)
       << E["mp2"] << std::endl;

    RT_ARRAYS_ALLOCATE("cc.vt2(ijab)", (no,no,N,N));
    RT_ARRAYS_ALLOCATE("cc.vt1(ijab)", (no,no,N,N));
    RT_ARRAYS_ALLOCATE("cc.vt1\'", (no,no,N,N));
    RT_ARRAYS_ALLOCATE("cc.vt1\"", (no,no,N,N));
    RT_ARRAYS_ALLOCATE("cc.t(ia)", (no, nv));

    A["t(ia)"] = rt.arrays().find<Array>("cc.t(ia)");
    A["vt1"] = rt.arrays().find<Array>("cc.vt1(ijab)");
    A["vt2"] = rt.arrays().find<Array>("cc.vt2(ijab)");
    A["vt1\'"] = rt.arrays().find<Array>("cc.vt1\'");
    A["vt1\""] = rt.arrays().find<Array>("cc.vt1\"");

    std::auto_ptr<DIIS> diis;
    if (pe.rank() == 0) {
	File::Group fg = rt.file("cc").create_group("diis");
	diis.reset(new DIIS(no, nv, fg, rt.get<int>("/cc/diis/max", 5)));
    }

    E["ccsd"] = sd(rt).energy(pe, wf, V, A, diis.get());

    rt.arrays().erase<Array>("cc.vt2(ijab)");
    rt.arrays().erase<Array>("cc.vt1(ijab)");
    rt.arrays().erase<Array>("cc.vt1\"");
    rt.arrays().erase<Array>("cc.vt1\'");

    V.clear();
    A.clear();

    if (method == "ccsd(t)") {
	Map<Array*> V;

	RT_ARRAYS_ALLOCATE("integrals.v(iabc)", (no,nv,N,nv));
	V["iqrs"] = rt.arrays().find<Array>("integrals.v(iqrs)");
	V["iabc"] = rt.arrays().find<Array>("integrals.v(iabc)");

	{
	    utility::timer timer;
	    double cutoff = rt.get<double>("/cc/integrals/cutoff", 1e-10);
	    integrals::Screening screening(wf.basis(), cutoff);
	    cc::integrals(pe, wf, V, screening);
	    pe.cout() << "integrals time: " << timer << std::endl;
	}

	V["ijka"] = rt.arrays().find<Array>("integrals.v(ijka)");
	V["ijab"] = rt.arrays().find<Array>("integrals.v(ijab)");

	Map<Array*> t;
	t["ia"] = rt.arrays().find<Array>("cc.t(ia)");
	t["ijab"] = rt.arrays().find<Array>("cc.t(ijab)");

	cc::Energy e = cc::triples::energy(pe, wf, V, t);

	E["ccsd[t]"] = E["ccsd"] + e["[t]"];
	E["ccsd(t)"] = E["ccsd[t]"] + e["(t)"];

	pe.cout() << "ccsd[t]: " << E["ccsd[t]"] << std::endl;
	pe.cout() << "ccsd(t): " << E["ccsd(t)"] << std::endl;
    }

    return E;

}

