
#include "mpqc_init.h"

#include "mpqc/mpqc_config.h"

#include <libgen.h>
#include <clocale>
#include <memory>
#include <sstream>

#include <libint2.hpp>

#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_FENV_H
#include <fenv.h>
#endif

// if have c++17
#if __cplusplus > 201402L
#include <filesystem>
using std::filesystem::current_path;
using std::filesystem::exists;
using std::filesystem::is_directory;
using std::filesystem::path;
#else
#include <boost/filesystem.hpp>
using boost::filesystem::current_path;
using boost::filesystem::exists;
using boost::filesystem::is_directory;
using boost::filesystem::path;
#endif

#include "mpqc/util/keyval/keyval.h"
//#include "mpqc/util/misc/consumableresources.h"
#include "mpqc/util/core/exception.h"
#include "mpqc/util/core/exenv.h"
#include "mpqc/util/core/formio.h"
#include "mpqc/util/external/madworld/parallel_file.h"

namespace mpqc {

std::unique_ptr<MPQCInit> MPQCInit::instance_;

void initialize(int &argc, char **argv, const madness::World &top_world,
                std::shared_ptr<GetLongOpt> opt) {
  if (!madness::initialized()) {
    throw ProgrammingError(
        "MADWorld runtime must be initialized before calling "
        "mpqc::initialize()",
        __FILE__, __LINE__);
  }
  if (MPQCInit::instance_) {
    throw ProgrammingError("Only one MPQCInit object can be created!", __FILE__,
                           __LINE__);
  } else {
    MPQCInit::instance_ =
        std::make_unique<MPQCInit>(argc, argv, opt, top_world, MPQCInit::singleton_ctor_tag{});
  }
}

void finalize() { MPQCInit::instance_.reset(); }

MPQCInit::MPQCInit(int &argc, char **argv, std::shared_ptr<GetLongOpt> opt,
                   const madness::World &top_world, singleton_ctor_tag)
    : opt_(opt), argv_(argv), argc_(argc), input_format_(InputFormat::invalid) {
  if (opt_)
    opt_->enroll("max_memory", GetLongOpt::MandatoryValue,
                 "the available memory");

  atexit(mpqc::finalize);
  ExEnv::init(argc_, argv_);
  init_fp();
  init_limits();
  init_work_dir();
  const auto libint_is_verbose = false;
  libint2::initialize(!libint_is_verbose);  // turns off Libint diagnostics

  init_io(top_world);
  // init_resources(keyval);
}

MPQCInit::~MPQCInit() {
  libint2::finalize();
  //  ConsumableResources::set_default_instance(0);
  FormIO::set_default_basename({});
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

namespace {
std::shared_ptr<mpqc::KeyVal> __make_keyval(madness::World &world,
                                            const std::string &filename,
                                            MPQCInit::InputFormat format) {
  std::shared_ptr<mpqc::KeyVal> kv;
  std::stringstream ss;
  utility::parallel_read_file(world, filename, ss);
  kv = std::make_shared<mpqc::KeyVal>();
  switch (format) {
    case MPQCInit::InputFormat::xml:
      kv->read_xml(ss);
      break;
    case MPQCInit::InputFormat::json:
      kv->read_json(ss);
      break;
    case MPQCInit::InputFormat::info:
      kv->read_info(ss);
      break;
    default:
      abort();
  }
  return kv;
}
}  // namespace

std::shared_ptr<mpqc::KeyVal> MPQCInit::make_keyval(
    madness::World &world, const std::string &filename) {
  std::shared_ptr<mpqc::KeyVal> kv;
  try {
    kv = __make_keyval(world, filename, InputFormat::json);
    input_format_ = InputFormat::json;
  } catch (...) {
  }
  if (input_format_ == InputFormat::invalid) {
    try {
      kv = __make_keyval(world, filename, InputFormat::xml);
      input_format_ = InputFormat::xml;
    } catch (std::exception &e) {
    }
  }
  if (input_format_ == InputFormat::invalid) {
    try {
      kv = __make_keyval(world, filename, InputFormat::info);
      input_format_ = InputFormat::info;
    } catch (...) {
      std::cerr << "failed read_info" << std::endl;
    }
  }
  if (input_format_ == InputFormat::invalid)
    throw InputError(
        "did not recognize input file format (recognized formats: JSON, XML, "
        "INFO)",
        __FILE__, __LINE__);

  return kv;
}

void MPQCInit::init_io(const madness::World &top_world) {
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
  std::cerr << std::setprecision(std::numeric_limits<double>::max_digits10);
  FormIO::setindent(ExEnv::outn(), 2);
  FormIO::setindent(ExEnv::errn(), 2);
  FormIO::setindent(std::cout, 2);
  FormIO::setindent(std::cerr, 2);
  FormIO::set_printnode(0);
  if (top_world.size() > 1) FormIO::init_mp(top_world.rank());
}

// void MPQCInit::init_resources(Ref<KeyVal> keyval) {
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

void MPQCInit::init_work_dir() {
  char *mpqc_work_dir;
  mpqc_work_dir = std::getenv("MPQC_WORK_DIR");

  // if not user defined, use current path
  if (mpqc_work_dir == nullptr) {
    path curr_path = current_path();
    // set the work dir in FormIO
    FormIO::set_default_work_dir(curr_path.string());
  } else {
    // check the correctness of path
    path work_path(mpqc_work_dir);
    bool path_exists = exists(work_path);
    if (!path_exists) {
      throw FileOperationFailed(
          "Path ${MPQC_WORK_DIR} does not exists! Please set environment "
          "variable MPQC_WORK_DIR to a valid path.\n",
          __FILE__, __LINE__, mpqc_work_dir, FileOperationFailed::Exists);
    }
    bool is_dir = is_directory(work_path);
    if (!is_dir) {
      throw FileOperationFailed(
          "Path ${MPQC_WORK_DIR} is not a directory! Please set environment "
          "variable MPQC_WORK_DIR to an existing directory.\n",
          __FILE__, __LINE__, mpqc_work_dir, FileOperationFailed::Chdir);
    }
    // set the work dir in FormIO
    FormIO::set_default_work_dir(mpqc_work_dir);
  }
}

void MPQCInit::set_basename(const std::string &input_filename,
                            const std::string &output_filename) {
  // get the basename for output files:
  // 1) if output filename is given, use it minus the suffix
  // 2) if output filename is not given, use input's basename minus the suffix
  const char *basename_source;
  auto input_filename_copy_cstr = std::make_unique<char[]>(input_filename.size() + 1);
  input_filename_copy_cstr[input_filename.size()] = '\0';
  std::copy(cbegin(input_filename), cend(input_filename), input_filename_copy_cstr.get());
  if (!output_filename.empty())
    basename_source = output_filename.c_str();
  else {
    // get the basename(input). basename() in libgen.h is in POSIX standard
    basename_source = basename(input_filename_copy_cstr.get());
  }
  // if basename_source does not contain '.' this will return a null pointer
  const char *dot_position = ::strrchr(basename_source, '.');
  const int nfilebase = (dot_position) ? (int)(dot_position - basename_source)
                                       : strlen(basename_source);
  std::string basename(basename_source, basename_source+nfilebase);
  FormIO::set_default_basename(std::move(basename));
}

std::shared_ptr<GetLongOpt> make_options() {
  // parse commandline options
  std::shared_ptr<GetLongOpt> options = std::make_shared<GetLongOpt>();

  options->usage(
      "[options] [input_file.{json|xml|info}]\nThe input file name can be "
      "given as the last argument unless an option with omitted optional value "
      "is used (e.g. -D)");
  options->enroll("i", GetLongOpt::MandatoryValue,
                  "the name of the input file");
  options->enroll("o", GetLongOpt::MandatoryValue,
                  "the name of the output file");
  options->enroll("p", GetLongOpt::MandatoryValue,
                  "prefix for all relative paths in KeyVal");
  options->enroll("W", GetLongOpt::MandatoryValue, "set the working directory",
                  ".");
  options->enroll("u", GetLongOpt::MandatoryValue, "the units system");
  options->enroll("d", GetLongOpt::NoValue,
                  "start the program and attach a debugger");
  options->enroll("t", GetLongOpt::NoValue,
                  "throw if a deprecated keyword is read");
  options->enroll("D", GetLongOpt::OptionalValue,
                  " if \"debugger\" keyword is not given, create a default "
                  "debugger; an optional JSON string can be given to provide "
                  "KeyVal ctor input");
  // options->enroll("c", GetLongOpt::NoValue, "check input then exit");
  options->enroll("v", GetLongOpt::NoValue, "print the version number");
  options->enroll("w", GetLongOpt::NoValue, "print the warranty");
  options->enroll("L", GetLongOpt::NoValue, "print the license");
  options->enroll(
      "k", GetLongOpt::NoValue,
      "print all registered (KeyVal-constructible) DescribedClass classes");
  options->enroll("h", GetLongOpt::NoValue, "print this message");

  return options;
}

std::tuple<std::string, std::string> process_options(
    const std::shared_ptr<GetLongOpt> &options) {
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

  auto print_std_header = []() {
    ExEnv::out0() << indent << "MPQC version " << MPQC_VERSION << std::endl
                  << indent << "compiled for " << TARGET_ARCH << std::endl
                  << FormIO::copyright << std::endl;
  };

  if (options->retrieve("h")) {
    print_std_header();
    options->usage(ExEnv::out0());
    std::exit(0);
  }

  if (options->retrieve("v")) {
    print_std_header();
    std::exit(0);
  }

  if (options->retrieve("w")) {
    print_std_header();
    ExEnv::out0() << FormIO::warranty << std::endl;
    std::exit(0);
  }

  if (options->retrieve("L")) {
    print_std_header();
    ExEnv::out0() << FormIO::license << std::endl;
    std::exit(0);
  }

  if (options->retrieve("k")) {
    print_std_header();
    ExEnv::out0() << std::endl
                  << indent
                  << "DescribedClass KeyVal-ctor registry:" << incindent
                  << std::endl;
    for (const auto &entry : DescribedClass::keyval_ctor_registry()) {
      ExEnv::out0() << indent << entry.first << std::endl;
    }
    std::exit(0);
  }

  if (options->retrieve("t")) {
    KeyVal::set_throw_if_deprecated_path(true);
  }

  // get input file name ... can be given as the last option
  auto input_opt = options->retrieve("i");
  std::string input_filename;
  if (input_opt) {
    input_filename = *input_opt;
  } else if (MPQCInit::instance().argc() - options->first_unprocessed_arg() ==
             1) {
    input_filename =
        MPQCInit::instance().argv()[options->first_unprocessed_arg()];
  } else {
    options->usage();
    throw std::invalid_argument("input filename not given");
  }

  return std::make_tuple(input_filename, output_filename);
}

}  // namespace mpqc
