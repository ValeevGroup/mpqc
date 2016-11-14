
#include "mpqc_init.h"

#include "mpqc/mpqc_config.h"

#include <libgen.h>
#include <clocale>
#include <sstream>

#include <libint2.hpp>

#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_FENV_H
#  include <fenv.h>
#endif

#include "mpqc/util/keyval/keyval.h"
//#include "mpqc/util/misc/consumableresources.h"
#include "mpqc/util/misc/exception.h"
#include "mpqc/util/misc/formio.h"
#include "mpqc/util/misc/exenv.h"
#include "mpqc/util/external/c++/memory"
#include "mpqc/util/external/madworld/parallel_file.h"

namespace mpqc {

std::unique_ptr<MPQCInit> MPQCInit::instance_;

void initialize(int &argc, char **argv, std::shared_ptr<GetLongOpt> opt)
{
  if (!madness::initialized()) {
    throw ProgrammingError(
        "MADWorld runtime must be initialized before calling mpqc::initialize()",
        __FILE__, __LINE__);
  }
  if (MPQCInit::instance_) {
    throw ProgrammingError("Only one MPQCInit object can be created!", __FILE__,
                           __LINE__);
  } else {
    MPQCInit::instance_ = std::unique_ptr<MPQCInit>(new MPQCInit(argc, argv, opt));
  }
}

void finalize() {
  MPQCInit::instance_.reset();
}

MPQCInit::MPQCInit(int &argc, char **argv, std::shared_ptr<GetLongOpt> opt)
    : opt_(opt), argv_(argv), argc_(argc) {
  if (opt_)
    opt_->enroll("max_memory", GetLongOpt::MandatoryValue,
                 "the available memory");

  atexit(mpqc::finalize);
  ExEnv::init(argc_, argv_);
  init_fp();
  init_limits();

  libint2::initialize();

  init_io();
  //init_resources(keyval);
}

MPQCInit::~MPQCInit() {
  libint2::finalize();
//  ConsumableResources::set_default_instance(0);
  FormIO::set_default_basename(0);
}

MPQCInit &MPQCInit::instance() {
  if (instance_)
    return *instance_;
  else
    throw ProgrammingError(
        "mpqc::initialize() must be called before calling MPQCInit::instance()",
        __FILE__, __LINE__);
}

void MPQCInit::init_fp() {
#ifdef HAVE_FEENABLEEXCEPT
// this uses a glibc extension to trap on individual exceptions
#ifdef FE_DIVBYZERO
  feenableexcept(FE_DIVBYZERO);
#endif
#ifdef FE_INVALID
  feenableexcept(FE_INVALID);
#endif
#ifdef FE_OVERFLOW
  feenableexcept(FE_OVERFLOW);
#endif
#endif

#ifdef HAVE_FEDISABLEEXCEPT
// this uses a glibc extension to not trap on individual exceptions
#ifdef FE_UNDERFLOW
  fedisableexcept(FE_UNDERFLOW);
#endif
#ifdef FE_INEXACT
  fedisableexcept(FE_INEXACT);
#endif
#endif
}

void MPQCInit::init_limits() {
#if defined(HAVE_SETRLIMIT)
  struct rlimit rlim;
  rlim.rlim_cur = 0;
  rlim.rlim_max = 0;
  setrlimit(RLIMIT_CORE, &rlim);
#endif
}

std::shared_ptr<mpqc::KeyVal>
MPQCInit::make_keyval(madness::World& world, const std::string &filename) {
  std::stringstream ss;
  utility::parallel_read_file(world, filename, ss);

  std::shared_ptr<mpqc::KeyVal> kv = std::make_shared<mpqc::KeyVal>();
  kv->read_json(ss);
  kv->assign("world", &world);  // set "$:world" keyword to &world to define
                                // the default execution context for this input
  return kv;
}

void MPQCInit::init_io() {
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::cout << std::setprecision(15);
  FormIO::setindent(ExEnv::outn(), 2);
  FormIO::setindent(ExEnv::errn(), 2);
  FormIO::setindent(std::cout, 2);
  FormIO::setindent(std::cerr, 2);
  FormIO::set_printnode(0);
  auto& world = madness::World::get_default();
  if (world.size() > 1) FormIO::init_mp(world.rank());
}

//void MPQCInit::init_resources(Ref<KeyVal> keyval) {
//  // get the resources object. first try commandline and environment
//  Ref<ConsumableResources> inst =
//      ConsumableResources::initial_instance(argc_, argv_);
//
//  // if we still don't have one reading from the input
//  if (inst.null() && keyval) {
//    inst << keyval->describedclassvalue("resources");
//  }
//
//  if (inst) {
//    ConsumableResources::set_default_instance(inst);
//  } else {
//    ConsumableResources::get_default_instance();
//  }
//}

void MPQCInit::set_basename(const std::string &input_filename,
                            const std::string &output_filename) {
  // get the basename for output files:
  // 1) if output filename is given, use it minus the suffix
  // 2) if output filename is not given, use input's basename minus the suffix
  const char *basename_source;
  char *input_copy = ::strdup(input_filename.c_str());
  if (!output_filename.empty())
    basename_source = output_filename.c_str();
  else {
    // get the basename(input). basename() in libgen.h is in POSIX standard
    basename_source = basename(input_copy);
  }
  // if basename_source does not contain '.' this will return a null pointer
  const char *dot_position = ::strrchr(basename_source, '.');
  const int nfilebase = (dot_position) ? (int)(dot_position - basename_source)
                                       : strlen(basename_source);
  char *basename = new char[nfilebase + 1];
  strncpy(basename, basename_source, nfilebase);
  basename[nfilebase] = '\0';
  FormIO::set_default_basename(basename);
  free(input_copy);
}

std::shared_ptr<GetLongOpt> make_options() {
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

  return options;
}

std::tuple<std::string, std::string>
process_options(const std::shared_ptr<GetLongOpt>& options) {
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
  std::string output_filename = output_opt ? *output_opt : std::string();

  if (options->retrieve("h")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << std::endl
         << indent << "compiled for " << TARGET_ARCH << std::endl
         << FormIO::copyright << std::endl;
    options->usage(ExEnv::out0());
    std::exit(0);
  }

  if (options->retrieve("v")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << std::endl
         << indent << "compiled for " << TARGET_ARCH << std::endl
         << FormIO::copyright;
    std::exit(0);
  }

  if (options->retrieve("w")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << std::endl
         << indent << "compiled for " << TARGET_ARCH << std::endl
         << FormIO::copyright << std::endl
         << FormIO::warranty;
    std::exit(0);
  }

  if (options->retrieve("L")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << std::endl
         << indent << "compiled for " << TARGET_ARCH << std::endl
         << FormIO::copyright << std::endl
         << FormIO::license;
    std::exit(0);
  }

  // get input file name
  std::string input_filename;
  if (MPQCInit::instance().argc() - options->first_unprocessed_arg() == 1) {
    input_filename = MPQCInit::instance().argv()[options->first_unprocessed_arg()];
  }
  else {
    options->usage();
    throw std::invalid_argument("input filename not given");
  }

  return std::make_tuple(input_filename, output_filename);
}

}  // namespace mpqc

