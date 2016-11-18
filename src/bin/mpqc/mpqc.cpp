// Massively Parallel Quantum Chemistry program
// computes properties of a Wavefunction

#include <sstream>

#include <libint2.hpp>
#include <tiledarray.h>

#include "mpqc_task.h"
#include "mpqc/mpqc_config.h"
#include "mpqc/util/external/madworld/parallel_file.h"
#include "mpqc/util/external/madworld/parallel_print.h"

#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/misc/assert.h"
#include "mpqc/util/misc/exception.h"
#include "mpqc/util/misc/exenv.h"
#include "mpqc/util/options/GetLongOpt.h"

// linkage files to force linking in of ALL Wavefunction-based classes
// this list must be in sync with CMakeLists.txt
#include "mpqc/chemistry/qc/f12/linkage.h"
#include "mpqc/chemistry/qc/cc/linkage.h"
#include "mpqc/chemistry/qc/mbpt/linkage.h"
#include "mpqc/chemistry/qc/scf/linkage.h"
#include "mpqc_init.h"

namespace mpqc {

void announce() {
  const char title1[] = "MPQC4: Massively Parallel Quantum Chemistry (v4)";
  const char title2[] = "Version " MPQC_VERSION;
  const auto target_width = 80;
  ExEnv::out0() << std::endl;
  ExEnv::out0() << indent;
  for (auto i = 0; i < (target_width - sizeof(title1)) / 2; i++)
    ExEnv::out0() << ' ';
  ExEnv::out0() << title1 << std::endl;
  ExEnv::out0() << indent;
  for (auto i = 0; i < (target_width - sizeof(title2)) / 2; i++)
    ExEnv::out0() << ' ';
  ExEnv::out0() << title2 << std::endl << std::endl;
}

}  // namespace mpqc

int try_main(int argc, char *argv[], madness::World& world) {
  using namespace mpqc;

  // define default MPQC options
  auto options = make_options();

  // initialize MPQC
  initialize(argc, argv, options, world);

  // parse and process options
  options->parse(argc, argv);
  std::string input_filename, output_filename;
  std::tie(input_filename, output_filename) = process_options(options);

  // redirect the output to output_file
  std::ofstream output;
  if (!output_filename.empty()) output.open(output_filename);
  if (!output.good()) throw FileOperationFailed("failed to open output file",
                                                __FILE__, __LINE__, output_filename.c_str(),
                                                FileOperationFailed::OpenW);
  auto cout_streambuf_reset = [](std::streambuf *p) { std::cout.rdbuf(p); };
  std::unique_ptr<std::streambuf, decltype(cout_streambuf_reset)> cout_buffer_holder(
      std::cout.rdbuf(), cout_streambuf_reset);
  if (!output_filename.empty()) std::cout.rdbuf(output.rdbuf());

  MPQCInit::instance().set_basename(input_filename, output_filename);

  // make keyval
  std::shared_ptr<KeyVal> kv = MPQCInit::instance().make_keyval(world, input_filename);

  // redirect filenames in KeyVal to the directory given by -p cmdline option
  auto prefix_opt = options->retrieve("p");
  if (prefix_opt) { // set file prefix, if given
    kv->assign("file_prefix", *prefix_opt);
  }

  // announce ourselves
  announce();

  // run
  MPQCTask task(world, kv);
  task.run();
  ExEnv::out0() << indent << "Wfn energy is: " << kv->value<double>("wfn:energy") << std::endl;

  return 0;
}

int main(int argc, char *argv[]) {
  int rc = 0;

  madness::World* world_ptr;
  // initialize MADNESS first
  try {
    world_ptr = &madness::initialize(argc, argv);
  }
  catch (...) {
    std::cerr << "!! Failed to initialize MADWorld: " << "\n";
    return 1;
  }

  try {
    try_main(argc, argv, *world_ptr);

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
