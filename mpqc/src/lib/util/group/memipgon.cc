
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

///////////////////////////////////////////////////////////////////////
// The IParagonMemoryGrp class

#define CLASSNAME IParagonMemoryGrp
#define HAVE_KEYVAL_CTOR
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

IParagonMemoryGrp::IParagonMemoryGrp(const RefKeyVal& keyval):
  MIDMemoryGrp(keyval)
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
  return 0;
}

void
IParagonMemoryGrp::unlock(long oldvalue)
{
}

long
IParagonMemoryGrp::send(void* data, int nbytes, int node, int type)
{
  long mid = isend(type, (char*)data, nbytes, node, 0);
  if (debug_) cout << me() << ": IParagonMemoryGrp::send(void*, "
                   << nbytes << ", "
                   << node << ", "
                   << type << "): mid = " << mid
                   << endl;
  return mid;
}

long
IParagonMemoryGrp::recv(void* data, int nbytes, int node, int type)
{
  long mid = irecvx(type, (char*)data, nbytes, node, 0, msginfo);
  if (debug_) cout << me() << ": IParagonMemoryGrp::recv(void*, "
                   << nbytes << ", "
                   << node << ", "
                   << type << "): mid = " << mid
                   << endl;
  return mid;
}

long
IParagonMemoryGrp::postrecv(void *data, int nbytes, int type)
{
  fprintf(stderr, "IParagonMemoryGrp::postrecv: not implemented\n");
  sleep(1);
  abort();
  return 0;
}

long
IParagonMemoryGrp::wait(long mid1, long mid2)
{
  if (debug_) {
      cout << me() << ": IParagonMemoryGrp::wait("
           << mid1 << ", "
           << mid2 << ")"
           << endl;
    }

  if (mid1 == -1) {
      cerr << me() << ": IParagonMemoryGrp::wait: mid1 == -1" << endl;
      sleep(1);
      abort();
    }
  else if (mid2 == -1) {
      msgwait(mid1);
      if (debug_)
          cout << me() << ": IParagonMemoryGrp::wait(): got " << mid1 << endl;
      return mid1;
    }
  else {
      while(1) {
          if (msgdone(mid1)) {
              if (debug_)
                  cout << me() << ": IParagonMemoryGrp::wait(): got "
                       << mid1 << endl;
              return mid1;
            }
          if (msgdone(mid2)) {
              if (debug_)
                  cout << me() << ": IParagonMemoryGrp::wait(): got "
                       << mid2 << endl;
              return mid2;
            }
        }
    }
}

#endif
