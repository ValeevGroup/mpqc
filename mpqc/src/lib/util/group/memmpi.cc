
#ifndef _util_group_memmpi_cc
#define _util_group_memmpi_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <unistd.h>
#include <util/group/memmpi.h>

#include <mpi.h>
#include <mpproto.h>

///////////////////////////////////////////////////////////////////////
// The MPIMemoryGrp class

#define CLASSNAME MPIMemoryGrp
#define HAVE_KEYVAL_CTOR
#define PARENTS public MIDMemoryGrp
#include <util/class/classi.h>
void *
MPIMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MIDMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

MPIMemoryGrp::MPIMemoryGrp(const RefMessageGrp& msg):
  MIDMemoryGrp(msg)
{
  if (debug_) cout << "MPIMemoryGrp entered" << endl;

  use_acknowledgments_ = 0;
  use_active_messages_ = 0;

  init_mid();

  activate();
}

MPIMemoryGrp::MPIMemoryGrp(const RefKeyVal& keyval):
  MIDMemoryGrp(keyval)
{
  if (debug_) cout << "MPIMemoryGrp keyval entered" << endl;

  use_acknowledgments_ = 0;
  use_active_messages_ = 0;

  init_mid();

  activate();
}

void
MPIMemoryGrp::init_mid()
{
  for (int i=0; i<max_mid; i++) mid_ready_[i] = 1;
}

long
MPIMemoryGrp::get_mid()
{
  for (int i=0; i<max_mid; i++) {
      if (!mid_ready_[i]) {
          mid_ready_[i] = 0;
          return i;
        }
    }

  cerr << "MPIMemoryGrp::get_mid(): ran out of mid's" << endl;
  abort();
  return 0;
}

void
MPIMemoryGrp::free_mid(long mid)
{
  mid_ready_[i] = 1;
}

long
MPIMemoryGrp::lock()
{
  return 0;
}

void
MPIMemoryGrp::unlock(long oldvalue)
{
}

long
MPIMemoryGrp::send(void* data, int nbytes, int node, int type)
{
  int mid = get_mid();
  //int MPI_Ibsend(void* buf, int count, MPI_Datatype datatype,
  //int dest, int tag, MPI_Comm comm, MPI_Request *request) 
  MPI_Ibsend(data, nbytes, MPI_BYTE, node, type,
             MPI_COMM_WORLD, &handles_[mid]);
  return mid;
}

long
MPIMemoryGrp::recv(void* data, int nbytes, int node, int type)
{
  int n;
  if (node == -1) n = DONTCARE;
  else n = node;
  int t = type;
  int mid = get_mid();
  // int MPI_Irecv(void* buf, int count, MPI_Datatype datatype, int
  // source, int tag, MPI_Comm comm, MPI_Request *request) 
  MPI_Irecv(data, nbytes, MPI_BYTE, n, t,
            MPI_COMM_WORLD, &handles_[mid]);
  if (debug_) cerr << "MPIMemoryGrp:: recv mid = " << mid << endl;
  return mid;
}

long
MPIMemoryGrp::postrecv(void *data, int nbytes, int type)
{
  cerr << "MPIMemoryGrp::postrecv: active messages not supported\n" << endl;
  abort();
  return 0;
}

long
MPIMemoryGrp::wait(long mid1, long mid2)
{
  MPI_Status status;
  if (debug_) {
      cout << me() << ": MPIMemoryGrp::wait("
           << mid1 << ", "
           << mid2 << ")"
           << endl;
    }

  if (mid1 == -1) {
      cerr << me() << ": MPIMemoryGrp::wait: mid1 == -1" << endl;
      sleep(1);
      abort();
    }
  else if (mid2 == -1) {
      MPI_Wait(&handle_[mid1], &status);
      free_mid(mid1);
      if (debug_)
          cout << me() << ": MPIMemoryGrp::wait(): got(1) " << mid1 << endl;
      return mid1;
    }
  else {
      while(1) {
          int flag;
          MPI_Test(&handle_[mid1], &flag, &status);
          if (flag) {
              free_mid(mid1);
              if (debug_)
                  cout << me() << ": MPIMemoryGrp::wait(): got(2a) "
                       << mid1 << endl;
              return mid1;
            }
          MPI_Test(&handle_[mid2], &flag, &status);
          if (flag) {
              free_mid(mid2);
              if (debug_)
                  cout << me() << ": MPIMemoryGrp::wait(): got(2b) "
                       << mid2 << endl;
              return mid2;
            }
        }
    }
}

MPIMemoryGrp::~MPIMemoryGrp()
{
  if (debug_) cerr << "MPIMemoryGrp: in DTOR" << endl;
  deactivate();
}

#endif
