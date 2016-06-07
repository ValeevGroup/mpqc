/// New MPQC Main file with KeyVal

#include "../include/tiledarray.h"

#include "../molecule/molecule.h"
#include "../utility/parallel_file.h"
#include "../utility/parallel_print.h"

using namespace mpqc;

int try_main(int argc, char *argv[], madness::World &world) {
  if (argc != 2) {
    std::cout << "usage: " << argv[0] << " <input_file.json>" << std::endl;
    throw std::invalid_argument("no input file given");
  }

  std::stringstream input_stream = utility::parallel_read_file(world, argv[1]);

  KeyVal kv;
  kv.read_json(input_stream);
  kv.assign("world", &world);

  /// construct Molecule

//  std::shared_ptr<molecule::Molecule> mol =
//      std::make_shared<molecule::Molecule>(kv.keyval("molecule"));

  std::shared_ptr<basis::Basis> bs =
      std::make_shared<basis::Basis>(kv.keyval("obs"));
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