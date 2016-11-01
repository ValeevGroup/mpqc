/// New MPQC Main file with KeyVal

#include "../include/tiledarray.h"

#include "../utility/parallel_file.h"
#include "../utility/parallel_print.h"

#include <mpqc/chemistry/qc/properties/energy.h>
#include <mpqc/chemistry/qc/wfn/wfn.h>
#include <mpqc/util/keyval/keyval.hpp>

// include linkage file
#include <mpqc/chemistry/molecule/linkage.h>

#include <mpqc/chemistry/qc/integrals/periodic_atomic_integral.h>

#include <sstream>
#include <clocale>

using namespace mpqc;

typedef std::vector<TA::DistArray<TA::TensorZ, TA::SparsePolicy>> TArrayVec;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic,
                      Eigen::RowMajor> Matrixc;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> Vectorc;

int try_main(int argc, char *argv[], madness::World &world) {
  if (argc != 2) {
    std::cout << "usage: " << argv[0] << " <input_file.json>" << std::endl;
    throw std::invalid_argument("no input file given");
  }

  std::stringstream ss;
  utility::parallel_read_file(world, argv[1], ss);

  KeyVal kv;
  kv.read_json(ss);
  kv.assign("world", &world);

  auto threshold = 1e-15;  // Hardcode for now.
  TiledArray::SparseShape<float>::threshold(threshold);

  libint2::initialize();

  auto mol = kv.keyval("molecule").class_ptr<molecule::Molecule>();
  auto charge = mol->charge();
  auto docc = mol->occupation(charge) / 2;
  auto enuc = mol->nuclear_repulsion();

  if (world.rank() == 0) {
    std::cout << *mol << std::endl;
    std::cout << "Nuclear Repulsion: " << enuc << std::endl;
  }

  /// Build 1-e integrals
  integrals::PeriodicAtomicIntegral<TA::TensorZ, TA::SparsePolicy> pao_int(kv);
  auto obr = std::make_shared<basis::OrbitalBasisRegistry>(kv);
  pao_int.set_orbital_basis_registry(obr);
  // Read max iter and thresholds
  auto maxiter = 50;
  auto d_conv = 1.0E-8;
  auto e_conv = 1.0E-8;

  // Kinetic matrix
  auto pT = pao_int.compute(L"<κ|T|λ>");

  // Nuclear-attraction matrix
  auto pV = pao_int.compute(L"<κ|V|λ>");

  // Overlap matrix in real space: <u0|vR>
  auto pS_R = pao_int.compute(L"<κ|λ>");

  // Overlap matrix in reciprocal space: Sum_R ( <u0|vR>exp(ikR) )
  auto pS_k = pao_int.transform_real2recip(pS_R);

  // One-body fock matrix in real space
  auto p_oneFock_real = pT;
  p_oneFock_real("mu, nu") += pV("mu, nu");

  // One-body fock matrix in reciprocal space
  auto p_oneFock_recip = pao_int.transform_real2recip(p_oneFock_real);

  // Diagonalize fock in recip space and form density matrix
  auto pD = pao_int.compute_density(p_oneFock_real, p_oneFock_recip,
                                    pS_k, docc);

  /// Main iterative loop
  auto iter = 0;
  auto rms = 0.0;
  auto converged = false;
  auto ehf = 0.0;
  auto ediff = 0.0;

  do {
      ++ iter;
      // Save a copy of energy and density
      auto ehf_old = ehf;
      auto pD_old = pD;

      // Coulomb matrix
      auto pJ = pao_int.compute(L"(μ ν| J|κ λ)");
      // Exchange matrix
      auto pK = pao_int.compute(L"(μ ν| K|κ λ)");
      // F = H + 2J - K
      auto pF = p_oneFock_real;
      pF("mu, nu") += 2.0 * pJ("mu, nu") - pK("mu, nu");

      // Transform Fock from real sapce to reciprocal space
      auto pF_recip = pao_int.transform_real2recip(pF);
      pD = pao_int.compute_density(pF, pF_recip, pS_k, docc);

      // Compute SCF energy
      pF("mu, nu") += p_oneFock_real("mu, nu");
      std::complex<double> E = pF("mu, nu") * pD("mu, nu");
      ehf = E.real();

      // Compute difference with last iteration
      ediff = ehf - ehf_old;
      auto pD_rms = pD;
      pD_rms("mu, nu") -= pD_old("mu, nu");
      rms = pD_rms("mu, nu").norm();
      if ((rms <= d_conv) || fabs(ediff) <= e_conv) converged = true;

      // Print out information
      std::string niter = "Iter", nEle = "E(HF)", nTot = "E(tot)",
                  nDel = "Delta(E)", nRMS = "RMS(D)", nT = "Time(s)";
      if (world.rank() == 0) {
          if (iter == 1)
            printf("\n\n %4s %20s %20s %20s %20s %13s\n", niter.c_str(),
                   nEle.c_str(), nTot.c_str(), nDel.c_str(), nRMS.c_str(),
                   nT.c_str());
          printf(" %4d %20.12f %20.12f %20.12f %20.12f\n", iter, ehf,
                 ehf + enuc, ediff, rms);
      }

  } while ((iter < maxiter) && (!converged));



  libint2::finalize();
  madness::finalize();

  return 0;
}

int main(int argc, char *argv[]) {
  int rc = 0;

  auto &world = madness::initialize(argc, argv);
  mpqc::utility::print_par(world, "MADNESS process total size: ", world.size(),
                           "\n");

  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::cout << std::setprecision(15);
  std::wcout.sync_with_stdio(false);
  std::wcerr.sync_with_stdio(false);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);

  try {
    try_main(argc, argv, world);

  } catch (TiledArray::Exception &e) {
    std::cerr << "!! TiledArray exception: " << e.what() << "\n";
    rc = 1;
  } catch (madness::MadnessException &e) {
    std::cerr << "!! MADNESS exception: " << e.what() << "\n";
    rc = 1;
  } catch (SafeMPI::Exception &e) {
    std::cerr << "!! SafeMPI exception: " << e.what() << "\n";
    rc = 1;
  } catch (std::exception &e) {
    std::cerr << "!! std exception: " << e.what() << "\n";
    rc = 1;
  } catch (...) {
    std::cerr << "!! exception: unknown exception\n";
    rc = 1;
  }

  return rc;
}
