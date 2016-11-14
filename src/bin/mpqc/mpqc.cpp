// Massively Parallel Quantum Chemistry program
// computes properties of a Wavefunction

#include <sstream>

#include <libint2.hpp>
#include <tiledarray.h>

#include "mpqcinit.h"
#include "mpqc/mpqc_config.h"
#include "mpqc/util/external/madworld/parallel_file.h"
#include "mpqc/util/external/madworld/parallel_print.h"

#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/misc/assert.h"
#include "mpqc/util/misc/exenv.h"
#include "mpqc/util/options/GetLongOpt.h"

// linkage files to force linking in of ALL Wavefunction-based classes
// this list must be in sync with CMakeLists.txt
#include "mpqc/chemistry/qc/cc/linkage.h"
#include "mpqc/chemistry/qc/f12/linkage.h"
#include "mpqc/chemistry/qc/mbpt/linkage.h"
#include "mpqc/chemistry/qc/scf/linkage.h"

using namespace mpqc;

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

int try_main(int argc, char *argv[], madness::World &world) {
  // parse commandline options
  std::shared_ptr<GetLongOpt> options = std::make_shared<GetLongOpt>();

  options->usage("[options] input_file.json");
  options->enroll("o", GetLongOpt::MandatoryValue, "the name of the output file");
  options->enroll("p", GetLongOpt::MandatoryValue, "prefix for all relative paths in KeyVal");
  options->enroll("W", GetLongOpt::MandatoryValue, "set the working directory",
                 ".");
  //options->enroll("c", GetLongOpt::NoValue, "check input then exit");
  options->enroll("v", GetLongOpt::NoValue, "print the version number");
  options->enroll("w", GetLongOpt::NoValue, "print the warranty");
  options->enroll("L", GetLongOpt::NoValue, "print the license");
  //options->enroll("d", GetLongOpt::NoValue, "debug");
  options->enroll("h", GetLongOpt::NoValue, "print this message");

  mpqc::initialize(argc, argv, options);

  const int optind = options->parse(argc, argv);

  // set the working dir
  if (*options->retrieve("W") != ".") {
    std::string dir = *options->retrieve("W");
    auto err = chdir(dir.c_str());
    if (err)
      throw FileOperationFailed("could not change directory", __FILE__,
                                __LINE__, dir.c_str(),
                                FileOperationFailed::Chdir);
  }

  // redirect the output, if needed
  auto output_opt = options->retrieve("o");
  std::ofstream output;
  std::string output_filename = output_opt ? *output_opt : std::string();
  if (!output_filename.empty()) output.open(output_filename);
  if (!output.good()) throw FileOperationFailed("failed to open output file",
                                                __FILE__, __LINE__, output_filename.c_str(),
                                                FileOperationFailed::OpenW);
  auto cout_streambuf_reset = [](std::streambuf *p) { std::cout.rdbuf(p); };
  std::unique_ptr<std::streambuf, decltype(cout_streambuf_reset)> cout_buffer_holder(
      std::cout.rdbuf(), cout_streambuf_reset);
  if (!output_filename.empty()) std::cout.rdbuf(output.rdbuf());

  if (options->retrieve("h")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << std::endl
         << indent << "compiled for " << TARGET_ARCH << std::endl
         << FormIO::copyright << std::endl;
    options->usage(ExEnv::out0());
    exit(0);
  }

  if (options->retrieve("v")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << std::endl
         << indent << "compiled for " << TARGET_ARCH << std::endl
         << FormIO::copyright;
    exit(0);
  }

  if (options->retrieve("w")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << std::endl
         << indent << "compiled for " << TARGET_ARCH << std::endl
         << FormIO::copyright << std::endl
         << FormIO::warranty;
    exit(0);
  }

  if (options->retrieve("L")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << std::endl
         << indent << "compiled for " << TARGET_ARCH << std::endl
         << FormIO::copyright << std::endl
         << FormIO::license;
    exit(0);
  }

  // get input file name
  std::string input_filename;
  if (argc - optind == 1) {
    input_filename = argv[optind];
  }
  else {
    options->usage();
    throw std::invalid_argument("input filename not given");
  }

  MPQCInit::instance().set_basename(input_filename, output_filename);
  std::shared_ptr<KeyVal> kv = MPQCInit::instance().make_keyval(world, input_filename);

  auto prefix_opt = options->retrieve("p");
  if (prefix_opt) { // set file prefix, if given
    kv->assign("file_prefix", *prefix_opt);
  }

  // announce ourselves
  announce();

  double threshold = kv->value<double>("sparse_threshold", 1e-20);
  TiledArray::SparseShape<float>::threshold(threshold);

  auto wfn = kv->keyval("wfn").class_ptr<qc::Wavefunction>();

  //  auto energy_prop = qc::Energy(kv);
  //  auto energy_prop_ptr = &energy_prop;

  double val = wfn->value();
  utility::print_par(world, "Wfn energy is: ", val, "\n");

  return 0;
}

int main(int argc, char *argv[]) {
  int rc = 0;

  // initialize MADNESS first
  try {
    madness::initialize(argc, argv);
  }
  catch (...) {
    std::cerr << "!! Failed to initialize MADWorld: " << "\n";
    return 1;
  }

  try {
    try_main(argc, argv, madness::World::get_default());

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
