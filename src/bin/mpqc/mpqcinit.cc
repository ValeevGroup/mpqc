
#include "mpqcinit.h"
#include <libgen.h>

#include <sstream>

#include <mpqc_config.h>

#ifdef HAVE_SYS_RESOURCE_H
#  include <sys/resource.h>
#endif
#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif

#include <util/keyval/keyval.h>
#include <util/group/message.h>
#if defined(HAVE_MPI) && defined(ALWAYS_USE_MPI)
#  include <util/group/messmpi.h>
#endif
#include <util/group/thread.h>
#include <util/group/memory.h>
#include <util/misc/autovec.h>
#include <util/group/pregtime.h>
#include <util/misc/scexception.h>
#include <math/scmat/matrix.h>
#include <chemistry/qc/basis/integral.h>
#include <util/misc/consumableresources.h>

namespace sc {

MPQCInit* MPQCInit::instance_(0);

static void
finalize()
{
  MemoryGrp::set_default_memorygrp(0);
  MessageGrp::set_default_messagegrp(0);
  ThreadGrp::set_default_threadgrp(0);
  SCMatrixKit::set_default_matrixkit(0);
  Integral::set_default_integral(0);
  ConsumableResources::set_default_instance(0);
  RegionTimer::set_default_regiontimer(0);
  SCFormIO::set_default_basename(0);
}

MPQCInit::MPQCInit(GetLongOpt&opt, int &argc, char **argv):
  opt_(opt), argc_(argc), argv_(argv)
{
  if (instance_ != 0) {
    throw ProgrammingError("Only one MPQCInit object can be created!",
                           __FILE__, __LINE__);
  }
  else
    instance_ = this;

  opt_.enroll("messagegrp", GetLongOpt::MandatoryValue,
              "which message group to use", 0);
  opt_.enroll("threadgrp", GetLongOpt::MandatoryValue,
              "which thread group to use", 0);
  opt_.enroll("memorygrp", GetLongOpt::MandatoryValue,
              "which memory group to use", 0);
  opt_.enroll("integral", GetLongOpt::MandatoryValue,
              "which integral evaluator to use", 0);
  opt_.enroll("resources", GetLongOpt::MandatoryValue,
              "the available consumable resources (memory, disk)", 0);
}

MPQCInit::~MPQCInit()
{
  finalize();
  instance_ = 0;
}

void
MPQCInit::init_fp()
{
#ifdef HAVE_FEENABLEEXCEPT
  // this uses a glibc extension to trap on individual exceptions
# ifdef FE_DIVBYZERO
  feenableexcept(FE_DIVBYZERO);
# endif
# ifdef FE_INVALID
  feenableexcept(FE_INVALID);
# endif
# ifdef FE_OVERFLOW
  feenableexcept(FE_OVERFLOW);
# endif
#endif

#ifdef HAVE_FEDISABLEEXCEPT
  // this uses a glibc extension to not trap on individual exceptions
# ifdef FE_UNDERFLOW
  fedisableexcept(FE_UNDERFLOW);
# endif
# ifdef FE_INEXACT
  fedisableexcept(FE_INEXACT);
# endif
#endif
}

void
MPQCInit::init_limits()
{
#if defined(HAVE_SETRLIMIT)
  struct rlimit rlim;
  rlim.rlim_cur = 0;
  rlim.rlim_max = 0;
  setrlimit(RLIMIT_CORE,&rlim);
#endif
}

Ref<MessageGrp>
MPQCInit::init_messagegrp()
{
  Ref<MessageGrp> grp;
#if defined(HAVE_MPI) && defined(ALWAYS_USE_MPI)
  grp = new MPIMessageGrp(&argc_, &argv_);
#endif

  // get the message group.  first try the commandline and environment
  if (grp.null()) grp = MessageGrp::initial_messagegrp(argc_, argv_);
  if (grp.nonnull()) MessageGrp::set_default_messagegrp(grp);
  else grp = MessageGrp::get_default_messagegrp();
  return grp;
}

Ref<KeyVal>
MPQCInit::init_keyval(const Ref<MessageGrp> &grp,const std::string &filename)
{
  sc::auto_vec<char> in_char_array;
  if (grp->me() == 0) {
    std::ifstream is(filename.c_str());
    if (!is.is_open()) {
      throw sc::FileOperationFailed("MPQCInit: failed to open input file",
                                    __FILE__, __LINE__,
                                    filename.c_str(),
                                    sc::FileOperationFailed::OpenR);
    }
    std::ostringstream ostrs;
    is >> ostrs.rdbuf();
    int n = 1 + strlen(ostrs.str().c_str());
    in_char_array.reset(strcpy(new char[n],ostrs.str().c_str()));
    grp->bcast(n);
    grp->bcast(in_char_array.get(), n);
    }
  else {
    int n;
    grp->bcast(n);
    in_char_array.reset(new char[n]);
    grp->bcast(in_char_array.get(), n);
    }

  Ref<ParsedKeyVal> parsedkv = new ParsedKeyVal();
  parsedkv->parse_string(in_char_array.get());

  Ref<KeyVal> kv = parsedkv;
  return kv;
}

Ref<ThreadGrp>
MPQCInit::init_threadgrp(Ref<KeyVal> keyval)
{
  // get the thread group.  first try the commandline and environment
  Ref<ThreadGrp> thread = ThreadGrp::initial_threadgrp(argc_, argv_);

  // if we still don't have a group, try reading the thread group
  // from the input
  if (thread.null() && keyval.nonnull()) {
    thread << keyval->describedclassvalue("thread");
  }

  if (thread.nonnull())
    ThreadGrp::set_default_threadgrp(thread);
  else
    thread = ThreadGrp::get_default_threadgrp();

  return thread;
}

Ref<MemoryGrp>
MPQCInit::init_memorygrp(Ref<KeyVal> keyval)
{
  // get the memory group.  first try the commandline and environment
  Ref<MemoryGrp> memory = MemoryGrp::initial_memorygrp(argc_, argv_);

  // if we still don't have a group, try reading the memory group
  // from the input
  if (memory.null() && keyval.nonnull()) {
    memory << keyval->describedclassvalue("memory");
  }

  if (memory.nonnull())
    MemoryGrp::set_default_memorygrp(memory);
  else
    memory = MemoryGrp::get_default_memorygrp();

  return memory;
}

void
MPQCInit::init_io(const Ref<MessageGrp> &grp)
{
  SCFormIO::setindent(ExEnv::outn(), 2);
  SCFormIO::setindent(ExEnv::errn(), 2);
  SCFormIO::setindent(std::cout, 2);
  SCFormIO::setindent(std::cerr, 2);
  SCFormIO::set_printnode(0);
  if (grp->n() > 1) SCFormIO::init_mp(grp->me());
}

void
MPQCInit::init_integrals(Ref<KeyVal> keyval)
{
  // get the integral factory. first try commandline and environment
  Ref<Integral> integral = Integral::initial_integral(argc_, argv_);

  // if we still don't have a integral, try reading the integral
  // from the input
  if (integral.null() && keyval.nonnull()) {
    integral << keyval->describedclassvalue("integrals");
  }

  if (integral.nonnull()) {
    Integral::set_default_integral(integral);
  }
  else {
    Integral::get_default_integral();
  }
}

void
MPQCInit::init_resources(Ref<KeyVal> keyval)
{
  // get the resources object. first try commandline and environment
  Ref<ConsumableResources> inst = ConsumableResources::initial_instance(argc_, argv_);

  // if we still don't have one reading from the input
  if (inst.null() && keyval.nonnull()) {
    inst << keyval->describedclassvalue("resources");
  }

  if (inst.nonnull()) {
    ConsumableResources::set_default_instance(inst);
  }
  else {
    ConsumableResources::get_default_instance();
  }
}

void
MPQCInit::init_timer(const Ref<MessageGrp> &grp, Ref<KeyVal> keyval)
{
  grp->sync(); // make sure nodes are sync'ed before starting timings
  Ref<RegionTimer> tim;
  if (keyval.nonnull())
    if (keyval->exists("timer"))
      tim << keyval->describedclassvalue("timer");
  if (tim.null())
    tim = new ParallelRegionTimer(grp,"mpqc",1,1);
  RegionTimer::set_default_regiontimer(tim);
}

void
MPQCInit::init_basename(const std::string &input_filename,
                        const std::string &output_filename)
{
  // get the basename for output files:
  // 1) if output filename is given, use it minus the suffix
  // 2) if output filename is not given, use input's basename minus the suffix
  const char *basename_source;
  char *input_copy = ::strdup(input_filename.c_str());
  if (!output_filename.empty()) basename_source = output_filename.c_str();
  else {
    // get the basename(input). basename() in libgen.h is in POSIX standard
    basename_source = basename(input_copy);
  }
  // if basename_source does not contain '.' this will return a null pointer
  const char* dot_position = ::strrchr(basename_source, '.');
  const int nfilebase =  (dot_position) ?
                           (int) (dot_position - basename_source) :
                           strlen(basename_source);
  char* basename = new char[nfilebase + 1];
  strncpy(basename, basename_source, nfilebase);
  basename[nfilebase] = '\0';
  SCFormIO::set_default_basename(basename);
  free(input_copy);
}

Ref<KeyVal>
MPQCInit::init(const std::string &input_filename,
               const std::string &output_filename)
{
  atexit(sc::finalize);
  ExEnv::init(argc_, argv_);
  init_fp();
  init_limits();
  Ref<MessageGrp> grp = init_messagegrp();
  init_io(grp);
  Ref<KeyVal> keyval = init_keyval(grp,input_filename);
  init_threadgrp(keyval);
  init_memorygrp(keyval);
  init_integrals(keyval);
  init_resources(keyval);
  init_timer(grp,keyval);
  init_basename(input_filename, output_filename);

  return keyval;
}

void
MPQCInit::finalize()
{
  sc::finalize();
}

}

