
#ifndef _util_group_memmpl_cc
#define _util_group_memmpl_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <unistd.h>
#include <util/group/memmpl.h>

#include <mpi.h>
#include <mpproto.h>

#define DISABLE do { fflush(stdout); } while(0)
#define ENABLE do { fflush(stdout); } while(0)


#define ACK 0

#define DEBUG 0
#define DEBREQ 0

#if DEBUG
#  undef PRINTF
#  define PRINTF(args) do { DISABLE; \
                            printf args; \
                            ENABLE; \
                           } while(0)
#else
#  undef PRINTF
#  define PRINTF(args)
#endif

///////////////////////////////////////////////////////////////////////
// The handler function and its data

static volatile int global_source, global_type, global_mid;
static MPLMemoryGrp *global_mpl_mem = 0;

static void
mpl_memory_handler(int*msgid_arg)
{
  long lmid = *msgid_arg;
  global_mpl_mem->handler(&lmid);
}

///////////////////////////////////////////////////////////////////////
// The MPLMemoryGrp class

#define CLASSNAME MPLMemoryGrp
#define PARENTS public MIDMemoryGrp
#include <util/class/classi.h>
void *
MPLMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MIDMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

long
MPLMemoryGrp::lock()
{
  int oldvalue;
  mpc_lockrnc(1, &oldvalue);
  return oldvalue;
}

void
MPLMemoryGrp::unlock(long oldvalue)
{
  int old = oldvalue;
  mpc_lockrnc(0, &old);
}

long
MPLMemoryGrp::send(void* data, int nbytes, int node, int type)
{
  int mid;
  mpc_send(data, nbytes, node, type, &mid);
  return mid;
}

long
MPLMemoryGrp::recv(void* data, int nbytes, int node, int type)
{
  int n;
  if (node == -1) n = DONTCARE;
  else n = node;
  int t = type;
  int mid;
  mpc_recv(data, nbytes, &n, &t, &mid);
  return mid;
}

long
MPLMemoryGrp::postrecv(void *data, int nbytes, int type)
{
  int global_source;
  int global_type;
  int global_mid;
  mpc_rcvncall(data, nbytes,
               &global_source, &global_type, &global_mid,
               mpl_memory_handler);
  return global_mid;
}

long
MPLMemoryGrp::wait(long mid1, long mid2)
{
  int imid;
  if (mid2 == -1) imid = (int)mid1;
  else imid = DONTCARE;
  size_t count;
  mpc_wait(&imid,&count);
  return imid;
}

MPLMemoryGrp::MPLMemoryGrp(const RefMessageGrp& msg,
                           int localsize):
  MIDMemoryGrp(msg, localsize)
{
  if (global_mpl_mem) {
      fprintf(stderr, "MPLMemoryGrp: only one allowed at a time\n");
      sleep(1);
      abort();
    }

  global_mpl_mem = this;

  use_acknowledgments_ = 0;
  use_active_messages_ = 1;

  activate();
}

MPLMemoryGrp::~MPLMemoryGrp()
{
}

void
MPLMemoryGrp::deactivate()
{
  if (!global_mpl_mem) return;
  MIDMemoryGrp::deactivate();
}

#endif
