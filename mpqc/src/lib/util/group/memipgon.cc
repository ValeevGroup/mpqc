
#ifndef _util_group_memipgon_cc
#define _util_group_memipgon_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <unistd.h>
#include <util/group/memipgon.h>

extern "C" {
#include <nx.h>
void msgwait(long);
}

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
// The IParagonMemoryGrp class

#define CLASSNAME IParagonMemoryGrp
#define PARENTS public MIDMemoryGrp
#include <util/class/classi.h>
void *
IParagonMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MIDMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

IParagonMemoryGrp::IParagonMemoryGrp(const RefMessageGrp& msg):
  MIDMemoryGrp(msg)
{
  use_acknowledgments_ = 0;
  use_active_messages_ = 0;
}

IParagonMemoryGrp::~IParagonMemoryGrp()
{
}

long
IParagonMemoryGrp::lock()
{
}

void
IParagonMemoryGrp::unlock(long oldvalue)
{
}

long
IParagonMemoryGrp::send(void* data, int nbytes, int node, int type)
{
  return isend(type, (char*)data, nbytes, node, 0);
}

long
IParagonMemoryGrp::recv(void* data, int nbytes, int node, int type)
{
  return irecvx(type, (char*)data, nbytes, node, 0, msginfo);
}

long
IParagonMemoryGrp::postrecv(void *data, int nbytes, int type)
{
  fprintf(stderr, "IParagonMemoryGrp::postrecv: not implemented\n");
  sleep(1);
  abort();
}

long
IParagonMemoryGrp::wait(long mid1, long mid2)
{
  if (mid1 == -1) {
      fprintf(stderr, "IParagonMemoryGrp::wait: mid1 == -1\n");
      sleep(1);
      abort();
    }
  else if (mid2 == -1) {
      msgwait(mid1);
      return mid1;
    }
  else {
      while(1) {
          if (msgdone(mid1)) {
              return mid1;
            }
          if (msgdone(mid2)) {
              return mid2;
            }
        }
    }
}

#endif
