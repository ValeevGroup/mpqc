#include "mpqc/chemistry/qc/phf/phf.h"

#include <clocale>
#include <sstream>

#include <tiledarray.h>

#include "mpqc/util/external/madworld/parallel_file.h"
#include "mpqc/util/external/madworld/parallel_print.h"
#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/util/keyval/keyval.h"

namespace mpqc {
namespace phf {

PHF::PHF(const KeyVal &kv) : Wavefunction(kv), kv_(kv), pao_factory_(kv) {}

void PHF::init(const KeyVal& kv) {

    // retrieve world from pao_factory
    auto& world = pao_factory_.world();

    auto init_start = mpqc::fenced_now(world);

    // set OrbitalBasisRegistry
    pao_factory_.set_orbital_basis_registry(this->wfn_world()->basis_registry());

    //TODO: retrive unitcell from pao_factory
    auto& unitcell = *kv.keyval("wfn_world:molecule").class_ptr<UnitCell>();

    if (world.rank() == 0) {
        std::cout << pao_factory_ << std::endl;
        std::cout << unitcell << std::endl;
    }

    auto charge = unitcell.charge();
    docc_ = unitcell.occupation(charge) / 2;

    // compute nuclear-repulsion energy
    auto RJ_max = pao_factory_.RJ_max();
    repulsion_ = unitcell.nuclear_repulsion(RJ_max);
    if (world.rank() == 0)
      std::cout << "\nNuclear Repulsion: " << repulsion_ << std::endl;

    e_converge_ = kv.value<double>("converge", 1.0E-8); // 1.0E(-N)
    d_converge_ = pow(10.0, log(e_converge_) / log(10) / 2.0); // 1.0E(-N/2)
    maxiter_ = kv.value<int64_t>("max_iter", 30);

    T_ = pao_factory_.compute(L"<κ|T|λ>"); // Kinetic
    V_ = pao_factory_.compute(L"<κ|V|λ>"); // Nuclear-attraction
    S_ = pao_factory_.compute(L"<κ|λ>"); // Overlap in real space
    Sk_ = pao_factory_.transform_real2recip(S_); // Overlap in reciprocal space
    H_("mu, nu") = T_("mu, nu") + V_("mu, nu"); // One-body hamiltonian in real space
    Hk_ = pao_factory_.transform_real2recip(H_); // One-body hamiltonian in reciprocal space

    // compute density matrix using core guess
    // TODO: write SOAD guess
    D_ = pao_factory_.compute_density(Hk_, Sk_, docc_);

    auto init_end = mpqc::fenced_now(world);
    std::cout << "\nPHF init time: " << mpqc::duration_in_s(init_start, init_end) << "s\n";
}

bool PHF::solve() {

    auto& world = pao_factory_.world();

    auto iter = 0;
    auto rms = 0.0;
    TArray Ddiff;
    auto converged = false;
    auto ephf = 0.0;
    auto ediff = 0.0;

    do {
        auto iter_start = mpqc::fenced_now(world);

        ++iter;
        // Save a copy of energy and density
        auto ephf_old = ephf;
        auto D_old = D_;

        J_ = pao_factory_.compute(L"(μ ν| J|κ λ)"); // Coulomb
        K_ = pao_factory_.compute(L"(μ ν| K|κ λ)"); // Exchange

        // F = H + 2J - K
        F_ = H_;
        F_("mu, nu") += 2.0 * J_("mu, nu") - K_("mu, nu");

        Fk_ = pao_factory_.transform_real2recip(F_);
        D_ = pao_factory_.compute_density(Fk_, Sk_, docc_);

        // compute PHF energy
        F_("mu, nu") += H_("mu, nu");
        std::complex<double> e_complex = F_("mu, nu") * D_("mu, nu");
        ephf = e_complex.real();

        // compute difference with last iteration
        ediff = ephf - ephf_old;
        Ddiff("mu, nu") = D_("mu, nu") - D_old("mu, nu");
        rms = Ddiff("mu, nu").norm();
        if ((rms <= d_converge_) || fabs(ediff) <= e_converge_) converged = true;

        auto iter_end = mpqc::fenced_now(world);
        auto iter_duration = mpqc::duration_in_s(iter_start, iter_end);
        // Print out information
        std::string niter = "Iter", nEle = "E(HF)", nTot = "E(tot)",
                    nDel = "Delta(E)", nRMS = "RMS(D)", nT = "Time(s)";
        if (world.rank() == 0) {
          if (iter == 1)
            printf("\n\n %4s %20s %20s %20s %20s %20s\n", niter.c_str(),
                   nEle.c_str(), nTot.c_str(), nDel.c_str(), nRMS.c_str(),
                   nT.c_str());
          printf(" %4d %20.12f %20.12f %20.12f %20.12f %20.3f\n", iter, ephf, ephf + repulsion_,
                 ediff, rms, iter_duration);
        }

    } while ((iter < maxiter_) && (!converged));

    // save total energy to energy_ no matter if PHF converges
    energy_ = ephf + repulsion_;

    if (!converged) {
        std::cout << "\nPeriodic Hartree-Fock iterations did not converge!\n" << std::endl;
        return false;
    }
    else {
        std::cout << "\nPeriodic Hartree-Fock iterations have converged!" << std::endl;
        printf("\nTotal Periodic Hartree-Fock energy = %20.12f\n", energy_);
        return true;
    }
}

double PHF::value() {
    init(kv_);
    solve();
    return energy_;
}

void PHF::obsolete() {
    Wavefunction::obsolete();
}

void PHF::compute(qc::PropertyBase *pb) {
    throw std::logic_error("Not implemented!");
}

} // end of namespace phf
} // end of namespace mpqc

