#include "chemistry/qc/ci/ci.h"
#include "util/misc/consumableresources.h"
#include <stdexcept>

//#define MPQC_PROFILE_ENABLE
#include "mpqc/utility/profile.hpp"

#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/integrals.hpp"
#include "mpqc/ci/direct.hpp"
#include "mpqc/utility/string.hpp"

#include <memory>


std::vector<double> sc::CI::compute(const Ref<RefWavefunction> &wfn,
                                    const Ref<KeyVal> &kv) {

    using mpqc::Vector;
    using mpqc::Matrix;

    ExEnv::out0() << "Beginning CI\n";
    
    mpqc::ci::Config config;
    Matrix C = Matrix(wfn->orbs()->coefs()).transpose();

    {
        size_t ne = wfn->nelectron();
        size_t no = C.rows();
        typedef KeyValValueint Int;
        config.core = kv->intvalue("core", Int(0));
        config.orbitals = kv->intvalue("orbitals", Int(no));
        config.alpha = kv->intvalue("alpha", Int((ne+1)/2));
        config.beta = kv->intvalue("beta", Int(ne/2));
        config.level = kv->intvalue("level", Int(0));
        config.max = kv->intvalue("max", Int(30));
        config.collapse = kv->intvalue("collapse", Int(config.collapse));
        config.cutoff = kv->intvalue("cutoff", Int(config.cutoff));
        config.block = kv->intvalue("block", Int(config.block));
    }

    Vector h;
    Matrix V;

    // compute molecular integrals
    {
        range mo(config.core, config.core+config.orbitals);
	range ao(0, C.cols());

	if (mo.size() > C.rows())
	    throw std::runtime_error("Number of orbitals too great");

	const auto &basis = wfn->basis();
	config.e_ref = basis->molecule()->nuclear_repulsion_energy();

        mpqc::ci::integrals(basis, C(mo, ao), wfn->integral()->hcore(), h);
        mpqc::ci::integrals(basis, C(mo, ao), wfn->integral()->electron_repulsion(), V);
    }

    mpqc::mpi::Comm comm(MPI_COMM_WORLD);

    std::string fname =
        ConsumableResources::get_default_instance()->disk_location() +
        SCFormIO::fileext_to_filename(".h5");
    
    std::auto_ptr<mpqc::File> file;
    file.reset(new mpqc::File(fname + "." + mpqc::string_cast(comm.rank())));

    std::vector<double> E;

    if (config.level > 0) {
        mpqc::ci::CI<mpqc::ci::Truncated> ci(config, comm, file->group());
        E = mpqc::ci::direct(ci, h, V);
    }
    else {
        mpqc::ci::CI<mpqc::ci::Full> ci(config, comm, file->group());
        E = mpqc::ci::direct(ci, h, V);
    }

    // delete file;

    
 
}
