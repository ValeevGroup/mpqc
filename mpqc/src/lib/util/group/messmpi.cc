
#include <mpi.h>
#include <util/keyval/keyval.h>
#include <util/group/messmpi.h>
#include <util/misc/formio.h>

//#define MPI_SEND_ROUTINE MPI_Ssend // hangs
//#define MPI_SEND_ROUTINE MPI_Send // hangs
#define MPI_SEND_ROUTINE MPI_Bsend // works requires the attach and detach

#define CLASSNAME MPIMessageGrp
#define PARENTS public MessageGrp
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
MPIMessageGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MessageGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

MPIMessageGrp::MPIMessageGrp()
{
  init();
}

MPIMessageGrp::MPIMessageGrp(const RefKeyVal& keyval):
  MessageGrp(keyval)
{
  init();
}

void
MPIMessageGrp::init()
{
  int me, nproc;
  int argc = 1;
  char **argv;
  argv = new char*[2];
  argv[0] = "-mpiB4"; // reduce the internal buffer since a user buffer is used
  argv[1] = 0;

  if (debug_) {
      cerr << "MPIMessageGrp::init: entered" << endl;
    }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  bufsize = 4000000;
  buf = (void*) new char[bufsize];
  MPI_Buffer_attach(buf,bufsize);
  initialize(me, nproc);

  if (debug_) {
      cerr << me << ": MPIMessageGrp::init: done" << endl;
    }
}

MPIMessageGrp::~MPIMessageGrp()
{
  MPI_Buffer_detach(&buf, &bufsize);
  delete[] (char*) buf;
  MPI_Finalize();
}

void
MPIMessageGrp::raw_send(int target, void* data, int nbyte)
{
  if (debug_) {
      cerr << scprintf("Node %d sending %d bytes to %d with tag %d\n",
                       me(), nbyte, target, 0)
           << endl;
    }
  MPI_SEND_ROUTINE(data,nbyte,MPI_BYTE,target,0,MPI_COMM_WORLD); 
  if (debug_) cerr << scprintf("Node %d sent\n", me()) << endl;
}

void
MPIMessageGrp::raw_recv(int sender, void* data, int nbyte)
{
  MPI_Status status;
  if (sender == -1) sender = MPI_ANY_SOURCE;
  if (debug_) {
      cerr << scprintf("Node %d recving %d bytes from %d with tag %d\n",
                       me(), nbyte, sender, 0)
           << endl;
    }
  MPI_Recv(data,nbyte,MPI_BYTE,sender,0,MPI_COMM_WORLD,&status);
  rnode = status.MPI_SOURCE;
  rtag = status.MPI_TAG;
  rlen = nbyte;
  if (debug_) cerr << scprintf("Node %d recvd %d bytes\n", me(), rlen) << endl;
}

void
MPIMessageGrp::raw_sendt(int target, int type, void* data, int nbyte)
{
  type = (type<<1) + 1;
  if (debug_) {
      cerr << scprintf("Node %d sending %d bytes to %d with tag %d\n",
                       me(), nbyte, target, type)
           << endl;
    }
  MPI_SEND_ROUTINE(data,nbyte,MPI_BYTE,target,type,MPI_COMM_WORLD); 
  if (debug_) cerr << scprintf("Node %d sent\n", me()) << endl;
}

void
MPIMessageGrp::raw_recvt(int type, void* data, int nbyte)
{
  MPI_Status status;
  if (type == -1) type = MPI_ANY_TAG;
  else type = (type<<1) + 1;
  if (debug_ ) {
      cerr << scprintf("Node %d recving %d bytes from %d with tag %d\n",
                       me(), nbyte, MPI_ANY_SOURCE, type)
           << endl;
    }
  MPI_Recv(data,nbyte,MPI_BYTE,MPI_ANY_SOURCE,type,MPI_COMM_WORLD,&status);
  rnode = status.MPI_SOURCE;
  rtag = status.MPI_TAG;
  rlen = nbyte;
  if (debug_) cerr << scprintf("Node %d recvd %d bytes\n", me(), rlen) << endl;
}

int
MPIMessageGrp::probet(int type)
{
  int flag;
  MPI_Status status;

  if (type == -1) type = MPI_ANY_TAG;
  else type = (type<<1) + 1;
  MPI_Iprobe(MPI_ANY_SOURCE,type,MPI_COMM_WORLD,&flag,&status);
  if (flag) {
    rnode = status.MPI_SOURCE;
    rtag = status.MPI_TAG;
    MPI_Get_count(&status, MPI_BYTE, &rlen);
    return 1;
    }
  else {
    rnode = rtag = rlen = 0;
    }
    
  return 0;
}

int
MPIMessageGrp::last_source()
{
  return rnode;
}

int
MPIMessageGrp::last_size()
{
  return rlen;
}

int
MPIMessageGrp::last_type()
{
  return rtag>>1;
}

void
MPIMessageGrp::sync()
{
  MPI_Barrier(MPI_COMM_WORLD);
}

static GrpReduce<double>* doublereduceobject;
static void
doublereduce(void*b, void*a, int*len, MPI_Datatype*datatype)
{
  doublereduceobject->reduce((double*)a, (double*)b, *len);
}
void
MPIMessageGrp::reduce(double*d, int n, GrpReduce<double>&r,
                          double*scratch, int target)
{
  doublereduceobject = &r;

  MPI_Op op;
  MPI_Op_create(doublereduce, 1, &op);

  double *work;
  if (!scratch) work = new double[n];
  else work = scratch;

  if (target == -1) {
      MPI_Allreduce(d, work, n, MPI_DOUBLE, op, MPI_COMM_WORLD);
    }
  else {
      MPI_Reduce(d, work, n, MPI_DOUBLE, op, target, MPI_COMM_WORLD);
    }

  if (target == -1 || target == me()) {
     for (int i=0; i<n; i++) d[i] = work[i];
    }

  MPI_Op_free(&op);

  if (!scratch) delete[] work;
}


void
MPIMessageGrp::raw_bcast(void* data, int nbyte, int from)
{
  if (n() == 1) return;

  if (debug_) {
      cerr << scprintf("Node %d bcast %d bytes from %d", me(), nbyte, from)
           << endl;
    }
  MPI_Bcast(data, nbyte, MPI_BYTE, from, MPI_COMM_WORLD);
  if (debug_) {
      cerr << scprintf("Node %d done with bcast", me()) << endl;
    }
}
