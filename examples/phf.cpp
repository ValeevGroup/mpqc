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

  if (world.rank() == 0) {
    std::cout << *mol << std::endl;
    std::cout << "Nuclear Repulsion: " << mol->nuclear_repulsion() << std::endl;
  }

  auto charge = mol->charge();
  auto docc = mol->occupation(charge) / 2;

  auto obr = std::make_shared<basis::OrbitalBasisRegistry>(kv);

  // Build 1-e integrals
  //  integrals::AtomicIntegral<TA::TensorD, TA::SparsePolicy> ao_int(kv);
  //  ao_int.set_orbital_basis_registry(obr);

  integrals::PeriodicAtomicIntegral<TA::TensorZ, TA::SparsePolicy> pao_int(kv);
  pao_int.set_orbital_basis_registry(obr);
  // Kinetic matrix
  auto pT = pao_int.compute(L"<κ|T|λ>");
//  if (world.rank() == 0)
//      std::cout << "\nKinetic matrix : \n" << pT << std::endl;

  // Nuclear-attraction matrix
  auto pV = pao_int.compute(L"<κ|V|λ>");
//  if (world.rank() == 0)
//      std::cout << "\nNuclear-attraction matrix : \n" << pV << std::endl;

  // Overlap matrix in real space: <u,0|v,R>
  auto pS_R = pao_int.compute(L"<κ|λ>");
//  if (world.rank() == 0)
//      std::cout << "\nOverlap matrix : \n" << pS_R << std::endl;

  // Overlap matrix in reciprocal space: Sum_R(<u,0|v,R>exp(ikR) )
  auto pS_k = pao_int.transform_real2recip(pS_R);
//  if (world.rank() == 0)
//      std::cout << "\nOverlap matrix in k space: \n" << pS_k << std::endl;

  // S^(-1/2)
//  auto pS_inv_sqrt = pao_int.compute(L"<κ|λ>[inv_sqr]");

  // One-body fock matrix in real space
  auto p_oneFock_real = pT;
  p_oneFock_real("mu, nu") += pV("mu, nu");

  // One-body fock matrix in reciprocal space
  auto p_oneFock_recip = pao_int.transform_real2recip(p_oneFock_real);

  // Diagonalize fock in recip space and form density matrix
//  auto pD = pao_int.compute_density(p_oneFock_real, p_oneFock_recip,
//                                    pS_inv_sqrt, docc);

  /// Main iterative loop
  auto iter = 0;
  auto rms = 0;
  auto converged = false;
  auto ehf = 0.0;
  auto ediff = 0.0;

  // 2-e 4-center Coulomb matrix
//  auto pERI = pao_int.compute(L"(μ ν| G|κ λ)");

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
