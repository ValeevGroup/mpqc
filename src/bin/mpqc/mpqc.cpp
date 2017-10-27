// Massively Parallel Quantum Chemistry program
// computes properties of a Wavefunction

#include <chrono>
#include <sstream>

#include <tiledarray.h>
#include <libint2.hpp>

#include "mpqc/chemistry/qc/properties/property.h"
#include "mpqc/mpqc_config.h"
#include "mpqc/mpqc_task.h"
#include "mpqc/util/external/madworld/parallel_file.h"
#include "mpqc/util/external/madworld/parallel_print.h"

#include "mpqc/chemistry/qc/lcao/wfn/wfn.h"
#include "mpqc/chemistry/units/units.h"
#include "mpqc/util/core/exception.h"
#include "mpqc/util/core/exenv.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/misc/assert.h"
#include "mpqc/util/misc/bug.h"
#include "mpqc/util/options/GetLongOpt.h"

// linkage files to force linking in of ALL Wavefunction-based classes
// this list must be in sync with CMakeLists.txt
#include "mpqc/chemistry/qc/lcao/cc/linkage.h"
#include "mpqc/chemistry/qc/lcao/ci/linkage.h"
#include "mpqc/chemistry/qc/lcao/f12/linkage.h"
#include "mpqc/chemistry/qc/lcao/mbpt/linkage.h"
#include "mpqc/chemistry/qc/lcao/scf/linkage.h"
#include "mpqc/chemistry/qc/properties/linkage.h"
#include "mpqc/mpqc_init.h"

namespace mpqc {

void announce(madness::World &world) {
  const char title1[] = "MPQC4: Massively Parallel Quantum Chemistry (v4)";
  const char title2[] = "Version " MPQC_VERSION;
  const char title3[] = "Revision " MPQC_REVISION;
  const auto target_width = 80;
  ExEnv::out0() << std::endl;
  ExEnv::out0() << indent;
  for (auto i = 0; i < (target_width - sizeof(title1)) / 2; i++)
    ExEnv::out0() << ' ';
  ExEnv::out0() << title1 << std::endl;
  ExEnv::out0() << indent;
  for (auto i = 0; i < (target_width - sizeof(title2)) / 2; i++)
    ExEnv::out0() << ' ';
  ExEnv::out0() << title2 << std::endl;
  ExEnv::out0() << indent;
  for (auto i = 0; i < (target_width - sizeof(title3)) / 2; i++)
    ExEnv::out0() << ' ';
  ExEnv::out0() << title3 << std::endl << std::endl;

  ExEnv::out0() << indent << printf("Machine:          %s", TARGET_ARCH)
                << std::endl
                << indent
                << printf("User:             %s@%s", ExEnv::username(),
                          ExEnv::hostname())
                << std::endl;
  std::time_t tm = std::time(nullptr);
  char tm_str[256];
  std::strftime(tm_str, 256, "%c %Z", std::gmtime(&tm));
  ExEnv::out0() << indent << "Start Time:       "
                << tm_str
                // << std::put_time(std::gmtime(&tm), "%c %Z")
                << std::endl;

  auto nproc = world.size();
  ExEnv::out0() << indent << "Default World:    " << nproc << " MPI process"
                << (nproc > 1 ? "es" : "") << std::endl
                << std::endl;
}

}  // namespace mpqc

int try_main(int argc, char *argv[], madness::World &world) {
  using namespace mpqc;

  // define default MPQC options
  auto options = make_options();

  // initialize MPQC
  initialize(argc, argv, world, options);

  // parse and process options
  options->parse(argc, argv);
  std::string input_filename, output_filename;
  std::tie(input_filename, output_filename) = process_options(options);

  // redirect the output to output_file
  std::ofstream output;
  if (!output_filename.empty()) output.open(output_filename);
  if (!output.good())
    throw FileOperationFailed("failed to open output file", __FILE__, __LINE__,
                              output_filename.c_str(),
                              FileOperationFailed::OpenW);
  auto cout_streambuf_reset = [](std::streambuf *p) { std::cout.rdbuf(p); };
  std::unique_ptr<std::streambuf, decltype(cout_streambuf_reset)>
      cout_buffer_holder(std::cout.rdbuf(), cout_streambuf_reset);
  if (!output_filename.empty()) std::cout.rdbuf(output.rdbuf());

  MPQCInit::instance().set_basename(input_filename, output_filename);

  // make keyval
  std::shared_ptr<KeyVal> kv =
      MPQCInit::instance().make_keyval(world, input_filename);

  // now set up the debugger
  auto debugger = kv->class_ptr<Debugger>("debugger");
  // use default-constructed debugger if -D given and debugger keyword is
  // missing
  if (!debugger && options->retrieve("D")) {
    debugger = std::make_shared<Debugger>();
  }
  if (debugger) {
    Debugger::set_default_debugger(debugger);
    debugger->set_exec(argv[0]);
    debugger->set_prefix(world.rank());
    if (options->retrieve("d"))
      debugger->debug("Starting debugger because -d given on command line.");
  }

  // redirect filenames in KeyVal to the directory given by -p cmdline option
  auto prefix_opt = options->retrieve("p");
  if (prefix_opt) {  // set file prefix, if given
    kv->assign("file_prefix", *prefix_opt);
  }

  // configure the units system
  auto units_opt = options->retrieve("u");
  std::string units_str;
  if (units_opt) {
    units_str = *units_opt;
  } else {
    if (kv->exists("units")) {
      units_str = kv->value<std::string>("units");
    }
  }
  if (!units_str.empty()) UnitFactory::set_default(units_str);

  // announce ourselves to the world
  announce(world);

  //////////////////////////////
  // print computation metadata
  //////////////////////////////

  // the input keyval for reference
  ExEnv::out0() << indent << "Input KeyVal (format="
                << to_string(MPQCInit::instance().input_format()) << "):\n";
  switch (MPQCInit::instance().input_format()) {
    case MPQCInit::InputFormat::json:
      kv->write_json(ExEnv::out0());
      break;
    case MPQCInit::InputFormat::xml:
      kv->write_xml(ExEnv::out0());
      break;
    case MPQCInit::InputFormat::info:
      kv->write_info(ExEnv::out0());
      break;
    default:
      throw ProgrammingError(
          "unrecognized input file format returned by MPQCInit::input_format()",
          __FILE__, __LINE__);
  }
  ExEnv::out0() << std::endl;

  // units
  ExEnv::out0() << indent << "Using fundamental constants system "
                << UnitFactory::get_default()->system() << std::endl;

  // run
  // set the default World to the top-level keyword "world"
  const auto kv_has_top_world = kv->exists("world");
  if (kv_has_top_world) {
    throw InputError("top-level keyword \"world\" is reserved for use by MPQC",
                     __FILE__, __LINE__, "world");
  }
  kv->assign("world", &world);   // set "$:world" keyword to &world to define
                                 // the default execution context for this input
  TA::set_default_world(world);  // must specify default world to avoid
                                 // madness::World::get_default() getting called
  MPQCTask task(world, kv);
  task.run();
  kv->erase("world");
  if (prefix_opt) {  // unset file prefix, if did previously
    kv->erase("file_prefix");
  }

  // print the final keyval
  ExEnv::out0() << indent << "Output KeyVal (format="
                << to_string(MPQCInit::instance().input_format()) << "):\n";
  switch (MPQCInit::instance().input_format()) {
    case MPQCInit::InputFormat::json:
      kv->write_json(ExEnv::out0());
      break;
    case MPQCInit::InputFormat::xml:
      kv->write_xml(ExEnv::out0());
      break;
    case MPQCInit::InputFormat::info:
      kv->write_info(ExEnv::out0());
      break;
    default:
      throw ProgrammingError(
          "unrecognized input file format returned by MPQCInit::input_format()",
          __FILE__, __LINE__);
  }
  ExEnv::out0() << std::endl;

  return 0;
}

int main(int argc, char *argv[]) {
  int rc = 0;

  madness::World *world_ptr;
  // initialize MADNESS first
  try {
    world_ptr = &madness::initialize(argc, argv);
  } catch (madness::MadnessException &e) {
    std::cerr << "!! MADNESS exception: " << e.what() << "\n";
    return 1;
  } catch (...) {
    std::cerr << "!! Failed to initialize MADWorld: "
              << "\n";
    return 1;
  }

  try {
    try_main(argc, argv, *world_ptr);

  } catch (mpqc::Exception &e) {
    std::cerr << "!! MPQC exception: " << e.what() << "\n";
    rc = 1;
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
