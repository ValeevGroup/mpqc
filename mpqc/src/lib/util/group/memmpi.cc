
#ifndef _util_group_memmpi_cc
#define _util_group_memmpi_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <unistd.h>
#include <util/group/memmpi.h>

#include <mpi.h>
#include <mpproto.h>

#define DISABLE do { fflush(stdout); } while(0)
#define ENABLE do { fflush(stdout); } while(0)

#define DEBUG 0

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
static MPIMemoryGrp *global_mpi_mem = 0;

static void
mpi_memory_handler(int*msgid_arg)
{
  long lmid = *msgid_arg;
  if (!global_mpi_mem) {
      fprintf(stderr,"WARNING: Tried to call mpi_memory_handler"
              " without global_mpi_mem\n");
    }
  else {
      global_mpi_mem->handler(&lmid);
    }
}

///////////////////////////////////////////////////////////////////////
// The MPIMemoryGrp class

#define CLASSNAME MPIMemoryGrp
#define PARENTS public MIDMemoryGrp
#include <util/class/classi.h>
void *
MPIMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MIDMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

long
MPIMemoryGrp::lock()
{
  int oldvalue;
  mpc_lockrnc(1, &oldvalue);
  return oldvalue;
}

void
MPIMemoryGrp::unlock(long oldvalue)
{
  int old = oldvalue;
  mpc_lockrnc(old, &old);
}

long
MPIMemoryGrp::send(void* data, int nbytes, int node, int type)
{
  int mid;
  mpc_send(data, nbytes, node, type, &mid);
  return mid;
}

long
MPIMemoryGrp::recv(void* data, int nbytes, int node, int type)
{
  int n;
  if (node == -1) n = DONTCARE;
  else n = node;
  int t = type;
  int mid;
  mpc_recv(data, nbytes, &n, &t, &mid);
  PRINTF(("MPIMemoryGrp:: recv mid = %d\n", mid));
  return mid;
}

long
MPIMemoryGrp::postrecv(void *data, int nbytes, int type)
{
  global_type = type;
  global_source = DONTCARE; 
#if HAVE_MPI
  mpc_rcvncall(data, nbytes,
               (int*)&global_source, (int*)&global_type, (int*)&global_mid,
               mpi_memory_handler);
  PRINTF(("MPIMemoryGrp:: postrecv mid = %d\n", global_mid));
#else
  cerr << "MPIMemoryGrp::postrecv: active messages not supported\n" << endl;
  abort();
#endif
  return global_mid;
}

long
MPIMemoryGrp::wait(long mid1, long mid2)
{
  int imid;
  if (mid2 == -1) imid = (int)mid1;
  else imid = DONTCARE;
  size_t count;
  PRINTF(("MPIMemoryGrp: waiting on %d\n", imid));
  if (mpc_wait(&imid,&count)) {
      fprintf(stderr, "MPIMemoryGrp: mpc_wait failed\n");
      sleep(1);
      abort();
    }
  PRINTF(("MPIMemoryGrp: imid = %d, global_mid = %d DONTCARE = %d count = %d\n",
          imid, global_mid, DONTCARE, count));
  return (long)imid;
}

MPIMemoryGrp::MPIMemoryGrp(const RefMessageGrp& msg):
  MIDMemoryGrp(msg)
{
  PRINTF(("MPIMemoryGrp entered\n"));
  if (global_mpi_mem) {
      fprintf(stderr, "MPIMemoryGrp: only one allowed at a time\n");
      sleep(1);
      abort();
    }

  global_mpi_mem = this;

  use_acknowledgments_ = 0;
#ifdef HAVE_MPL
  use_active_messages_ = 1;
#else
  use_active_messages_ = 0;
#endif

  PRINTF(("MPIMemoryGrp activating\n"));
  activate();
  PRINTF(("MPIMemoryGrp done\n"));
}

MPIMemoryGrp::MPIMemoryGrp(const RefKeyVal& keyval):
  MIDMemoryGrp(keyval)
{
  PRINTF(("MPIMemoryGrp entered\n"));
  if (global_mpi_mem) {
      fprintf(stderr, "MPIMemoryGrp: only one allowed at a time\n");
      sleep(1);
      abort();
    }

  global_mpi_mem = this;

  use_acknowledgments_ = 0;
#ifdef HAVE_MPL
  use_active_messages_ = keyval->boolvalue("active");
#else
  use_active_messages_ = 0;
#endif

  PRINTF(("MPIMemoryGrp activating\n"));
  activate();
  PRINTF(("MPIMemoryGrp done\n"));
}

MPIMemoryGrp::~MPIMemoryGrp()
{
  PRINTF(("MPIMemoryGrp: in DTOR\n"));
  deactivate();

  int oldlock = lock();
  global_mpi_mem = 0;
  unlock(oldlock);
}

void
MPIMemoryGrp::deactivate()
{
  if (!global_mpi_mem) return;
  MIDMemoryGrp::deactivate();
}

#endif
