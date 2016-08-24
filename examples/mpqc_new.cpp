/// New MPQC Main file with KeyVal

#include "../include/tiledarray.h"

#include "../utility/parallel_file.h"
#include "../utility/parallel_print.h"

#include <mpqc/chemistry/molecule/molecule.h>
#include <mpqc/chemistry/molecule/clustering_functions.h>
#include <mpqc/chemistry/molecule/make_clusters.h>
#include <mpqc/util/keyval/keyval.hpp>
#include <mpqc/chemistry/qc/wfn/ao_wfn.h>
#include <mpqc/chemistry/qc/properties/energy.h>
#include <mpqc/chemistry/qc/integrals/atomic_integral.h>

#include <sstream>

using namespace mpqc;

TA::TensorD ta_pass_through(TA::TensorD &&ten) { return std::move(ten); }

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

  auto threshold = 1e-11;  // Hardcode for now.
  TiledArray::SparseShape<float>::threshold(threshold);

  //
  // construct molecule
  //
  auto mol = molecule::Molecule(kv.keyval("molecule"));
  auto nclusters = 1;  // Hard Coded for now
  auto clustered_mol = std::make_shared<molecule::Molecule>(
      molecule::attach_hydrogens_and_kmeans(mol.clusterables(), nclusters));
  kv.assign("molecule", clustered_mol);

  //
  // construct basis registry
  //
  auto bs = basis::Basis(kv.keyval("obs"));
  auto dfbs = basis::Basis(kv.keyval("dfbs"));

  auto bs_registry = std::make_shared<OrbitalBasisRegistry>();
  bs_registry->add(OrbitalIndex(L"κ"), bs);
  bs_registry->add(OrbitalIndex(L"Κ"), dfbs);

  //
  // construct integrals
  //
  libint2::initialize();
  integrals::AtomicIntegral<TA::TensorD, TA::SparsePolicy> ao_int(
      world, ta_pass_through, kv.class_ptr<molecule::Molecule>("molecule"),
      bs_registry);

  kv.assign("ao_integrals", &ao_int);

  auto wfn_world = qc::WfnWorld(kv);
  kv.assign("wfn_world", &wfn_world);

  qc::Wfn* wfn = new qc::AOWfn(kv);

  auto energy_prop = qc::Energy(kv);
  auto energy_prop_ptr = &energy_prop;

  wfn->visit(energy_prop_ptr);
  std::cout << "Wfn energy is: " << *energy_prop.result() << std::endl;

  {  // Test ints from Wfn
    auto S = wfn->wfn_world()->ao_integrals().compute(L"<κ|λ>");
    // TA::DistArray<TA::TensorD, TA::SparsePolicy> S = ao_int.compute(L"<κ|λ>");
    assert(S.is_initialized());
  }

  delete wfn;

  return 0;
}

int main(int argc, char *argv[]) {
  int rc = 0;

  auto &world = madness::initialize(argc, argv);
  mpqc::utility::print_par(world, "MADNESS process total size: ", world.size(),
                           "\n");

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

  madness::finalize();
  return rc;
}
