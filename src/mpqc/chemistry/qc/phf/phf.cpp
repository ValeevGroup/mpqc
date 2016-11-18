#include "mpqc/chemistry/qc/phf/phf.h"
#include "periodic_soad.h"

#include <clocale>
#include <sstream>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/util/external/madworld/parallel_file.h"
#include "mpqc/util/external/madworld/parallel_print.h"
#include "mpqc/util/keyval/keyval.h"

namespace mpqc {
namespace phf {

PHF::PHF(const KeyVal& kv) : Wavefunction(kv), kv_(kv), pao_factory_(kv) {}

void PHF::init(const KeyVal& kv) {
  // TODO: retrive unitcell from pao_factory
  auto& unitcell = *kv.keyval("wfn_world:molecule").class_ptr<UnitCell>();

  converge_ = kv.value<double>("converge", 1.0E-8);  // 1.0E(-N)
  maxiter_ = kv.value<int64_t>("max_iter", 30);
  bool soad_guess = kv.value<bool>("soad_guess", true);
  print_detail_ = kv.value<bool>("print_detail", false);

  // retrieve world from pao_factory
  auto& world = pao_factory_.world();

  auto init_start = mpqc::fenced_now(world);

  // set OrbitalBasisRegistry
  pao_factory_.set_orbital_basis_registry(this->wfn_world()->basis_registry());

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

  T_ = pao_factory_.compute(L"<κ|T|λ>");        // Kinetic
  V_ = pao_factory_.compute(L"<κ|V|λ>");        // Nuclear-attraction
  S_ = pao_factory_.compute(L"<κ|λ>");          // Overlap in real space
  Sk_ = pao_factory_.transform_real2recip(S_);  // Overlap in reciprocal space
  H_("mu, nu") =
      T_("mu, nu") + V_("mu, nu");  // One-body hamiltonian in real space

  // compute density matrix using soad/core guess
  if (!soad_guess) {
    F_ = H_;
  } else {
    F_ = periodic_fock_soad(world, unitcell, H_, pao_factory_);
  }

  // transform Fock from real to reciprocal space
  Fk_ = pao_factory_.transform_real2recip(F_);
  // compute orthogonalizer matrix
  X_ = pao_factory_.gen_orthogonalizer(Sk_);
  // compute guess density
  D_ = pao_factory_.compute_density(Fk_, X_, docc_);

  auto init_end = mpqc::fenced_now(world);
  init_duration_ = mpqc::duration_in_s(init_start, init_end);
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

    if (print_detail_)
      if (world.rank() == 0) std::cout << "\nIteration: " << iter << "\n";

    // compute PHF energy
    F_("mu, nu") += H_("mu, nu");
    std::complex<double> e_complex = F_("mu, nu") * D_("mu, nu");
    ephf = e_complex.real();

    auto j_start = mpqc::fenced_now(world);
    J_ = pao_factory_.compute(L"(μ ν| J|κ λ)");  // Coulomb
    auto j_end = mpqc::fenced_now(world);
    j_duration_ += mpqc::duration_in_s(j_start, j_end);

    auto k_start = mpqc::fenced_now(world);
    K_ = pao_factory_.compute(L"(μ ν| K|κ λ)");  // Exchange
    auto k_end = mpqc::fenced_now(world);
    k_duration_ += mpqc::duration_in_s(k_start, k_end);

    // F = H + 2J - K
    F_ = H_;
    F_("mu, nu") += 2.0 * J_("mu, nu") - K_("mu, nu");

    // transform Fock from real to reciprocal space
    auto trans_start = mpqc::fenced_now(world);
    Fk_ = pao_factory_.transform_real2recip(F_);
    auto trans_end = mpqc::fenced_now(world);
    trans_duration_ += mpqc::duration_in_s(trans_start, trans_end);

    // compute new density
    auto d_start = mpqc::fenced_now(world);
    D_ = pao_factory_.compute_density(Fk_, X_, docc_);
    auto d_end = mpqc::fenced_now(world);
    d_duration_ += mpqc::duration_in_s(d_start, d_end);

    // compute difference with last iteration
    ediff = ephf - ephf_old;
    Ddiff("mu, nu") = D_("mu, nu") - D_old("mu, nu");
    rms = Ddiff("mu, nu").norm();
    if ((rms <= converge_) || fabs(ediff) <= converge_) converged = true;

    auto iter_end = mpqc::fenced_now(world);
    auto iter_duration = mpqc::duration_in_s(iter_start, iter_end);
    scf_duration_ += iter_duration;

    // Print out information
    if (print_detail_) {
      if (world.rank() == 0) {
        std::cout << "\nPHF Energy: " << ephf << "\n"
                  << "Total Energy: " << ephf + repulsion_ << "\n"
                  << "Delta(E): " << ediff << "\n"
                  << "RMS(D): " << rms << "\n"
                  << "Coulomb Build Time: "
                  << mpqc::duration_in_s(j_start, j_end) << " s\n"
                  << "Exchange Build Time: "
                  << mpqc::duration_in_s(k_start, k_end) << " s\n"
                  << "Transform Fock (Real->Recip) Time: "
                  << mpqc::duration_in_s(trans_start, trans_end) << " s\n"
                  << "Density Time: " << mpqc::duration_in_s(d_start, d_end)
                  << " s\n"
                  << "Iteration Time: " << iter_duration << " s\n";
      }
    } else {
      std::string niter = "Iter", nEle = "E(HF)", nTot = "E(tot)",
                  nDel = "Delta(E)", nRMS = "RMS(D)", nT = "Time(s)";
      if (world.rank() == 0) {
        if (iter == 1)
          printf("\n\n %4s %20s %20s %20s %20s %20s\n", niter.c_str(),
                 nEle.c_str(), nTot.c_str(), nDel.c_str(), nRMS.c_str(),
                 nT.c_str());
        printf(" %4d %20.12f %20.12f %20.12f %20.12f %20.3f\n", iter, ephf,
               ephf + repulsion_, ediff, rms, iter_duration);
      }
    }

  } while ((iter < maxiter_) && (!converged));

  // save total energy to energy_ no matter if PHF converges
  energy_ = ephf + repulsion_;

  if (!converged) {
    if (world.rank() == 0) {
      std::cout << "\nPeriodic Hartree-Fock iterations did not converge!\n"
                << std::endl;
    }
    return false;
  } else {
    if (world.rank() == 0) {
      std::cout << "\nPeriodic Hartree-Fock iterations have converged!"
                << std::endl;
      printf("\nTotal Periodic Hartree-Fock energy = %20.12f\n", energy_);

      if (print_detail_) {
          auto orb_e = pao_factory_.eps();
          Eigen::IOFormat fmt(5);
          std::cout << "\n k | orbital energies" << std::endl;
          for (auto k = 0; k < pao_factory_.k_size(); ++k) {
              std::cout << k << " | " << orb_e[k].real().transpose().format(fmt) << std::endl;
          }
      }

      // print out timings
      printf("\nTime(s):\n");
      printf("\tInit:                %20.3f\n", init_duration_);
      printf("\tCoulomb term:        %20.3f\n", j_duration_);
      printf("\tExchange term:       %20.3f\n", k_duration_);
      printf("\tReal->Recip trans:   %20.3f\n", trans_duration_);
      printf("\tDiag + Density:      %20.3f\n", d_duration_);
      printf("\tTotal:               %20.3f\n\n", scf_duration_);
    }
    return true;
  }
}

double PHF::value() {
  init(kv_);
  solve();
  return energy_;
}

void PHF::obsolete() { Wavefunction::obsolete(); }

void PHF::compute(qc::PropertyBase* pb) {
  throw std::logic_error("Not implemented!");
}

}  // end of namespace phf
}  // end of namespace mpqc
