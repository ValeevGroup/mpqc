
#include "mpqcinit.h"

#include <sstream>

#include <scconfig.h>

#ifdef HAVE_SYS_RESOURCE_H
#  include <sys/resource.h>
#endif
#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif

#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <util/group/memory.h>
#include <util/misc/autovec.h>
#include <util/group/pregtime.h>
#include <util/class/scexception.h>
#include <math/scmat/matrix.h>
#include <chemistry/qc/basis/integral.h>

namespace sc {

static void
finalize()
{
  MemoryGrp::set_default_memorygrp(0);
  MessageGrp::set_default_messagegrp(0);
  ThreadGrp::set_default_threadgrp(0);
  SCMatrixKit::set_default_matrixkit(0);
  Integral::set_default_integral(0);
  RegionTimer::set_default_regiontimer(0);
}

MPQCInit::MPQCInit(GetLongOpt&opt, int &argc, char **argv):
  opt_(opt), argc_(argc), argv_(argv)
{
  opt_.enroll("messagegrp", GetLongOpt::MandatoryValue,
              "which message group to use", 0);
  opt_.enroll("threadgrp", GetLongOpt::MandatoryValue,
              "which thread group to use", 0);
  opt_.enroll("memorygrp", GetLongOpt::MandatoryValue,
              "which memory group to use", 0);
}

MPQCInit::~MPQCInit()
{
  finalize();
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
MPQCInit::init_threadgrp(const Ref<KeyVal>&keyval)
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
MPQCInit::init_memorygrp(Ref<KeyVal> &keyval)
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

Ref<KeyVal>
MPQCInit::init(const std::string &filename)
{
  atexit(sc::finalize);
  ExEnv::init(argc_, argv_);
  init_fp();
  init_limits();
  Ref<MessageGrp> grp = init_messagegrp();
  init_io(grp);
  Ref<KeyVal> keyval = init_keyval(grp,filename);
  init_threadgrp(keyval);
  init_memorygrp(keyval);

  // these still need to be implemented
  //init_cca();
  //init_integrals();
  //init_timer();

  return keyval;
}

void
MPQCInit::finalize()
{
  sc::finalize();
}


}

