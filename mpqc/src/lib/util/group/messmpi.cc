//
// messmpi.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <stdio.h> // for sprintf
#include <unistd.h> // for fchdir etc.

#include <mpi.h>
extern int MPI_Initialized(int *); // missing in mpi.h

#include <util/keyval/keyval.h>
#include <util/group/messmpi.h>
#include <util/misc/formio.h>
#include <util/misc/newstring.h>

//#define MPI_SEND_ROUTINE MPI_Ssend // hangs
//#define MPI_SEND_ROUTINE MPI_Send // hangs
#define MPI_SEND_ROUTINE MPI_Bsend // works requires the attach and detach

///////////////////////////////////////////////////////////////////////

static
void
print_error_and_abort(int me, int mpierror)
{
  char msg[MPI_MAX_ERROR_STRING+1];
  int size;
  MPI_Error_string(mpierror, msg, &size);
  msg[size] = '\0';
  cerr << me << ": " << msg << endl;
  cout.flush();
  cerr.flush();
  MPI_Abort(MPI_COMM_WORLD, mpierror);
  abort();
}

///////////////////////////////////////////////////////////////////////
// The MPIMessageGrp class

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
  int argc = -1;
  char **argv = 0;
  if (keyval->exists("argv")) {
      argc = keyval->count("argv");
      argv = new char*[argc+1];
      argv[argc] = 0;
      for (int arg=0; arg<argc; arg++) {
          argv[arg] = keyval->pcharvalue("argv",arg);
        }
    }

  init(argc, argv);

  if (keyval->booleanvalue("errors_return")) {
      if (me()==0)
          cout << indent << "MPIMessageGrp: errors_return is true" << endl;
      MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    }

  SCFormIO::init_mp(me());
  if (debug_) {
      cerr << indent << "MPIMessageGrp: KeyVal CTOR: done" << endl;
    }
}

void
MPIMessageGrp::init(int argc,char **argv)
{
  int me, nproc;

  if (debug_) {
      cerr << "MPIMessageGrp::init: entered" << endl;
    }

  int flag;
  MPI_Initialized(&flag);
  if (!flag) {
      if (argc < 0) {
          argc = 1;
          argv = new char*[argc+1];
          // reduce the internal buffer since a user buffer is used
          argv[0] = "-mpiB4";
          argv[1] = 0;
        }
      // This dot business is to work around problems with some MPI
      // implementations.
      int dot = open(".",O_RDONLY);
      MPI_Init(&argc, &argv);
      fchdir(dot);
      close(dot);
    }
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
  int ret;
  if ((ret = MPI_SEND_ROUTINE(data,nbyte,MPI_BYTE,target,0,MPI_COMM_WORLD))
      != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::raw_send("
          << target << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
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
  int ret;
  if ((ret = MPI_Recv(data,nbyte,MPI_BYTE,sender,0,MPI_COMM_WORLD,&status))
      != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::raw_recv("
          << sender << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
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
  int ret;
  if ((ret = MPI_SEND_ROUTINE(data,nbyte,MPI_BYTE,target,type,MPI_COMM_WORLD))
      != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::raw_sendt("
          << target << "," << type << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
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
  int ret;
  if ((ret = MPI_Recv(data,nbyte,MPI_BYTE,MPI_ANY_SOURCE,
                      type,MPI_COMM_WORLD,&status)) != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::raw_recvt("
          << type << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
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
  int ret;
  if ((ret = MPI_Iprobe(MPI_ANY_SOURCE,type,MPI_COMM_WORLD,&flag,&status))
      != MPI_SUCCESS ) {
      cerr << me() << ": MPIMessageGrp::probet("
          << type << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
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
  int ret;
  if ((ret = MPI_Barrier(MPI_COMM_WORLD)) != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::sync(): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
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

  int ret;

  if (target == -1) {
      ret = MPI_Allreduce(d, work, n, MPI_DOUBLE, op, MPI_COMM_WORLD);
    }
  else {
      ret = MPI_Reduce(d, work, n, MPI_DOUBLE, op, target, MPI_COMM_WORLD);
    }

  if (ret != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::reduce(,"
          << n << ",,," << target << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
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
  int ret;
  if ((ret = MPI_Bcast(data, nbyte, MPI_BYTE, from, MPI_COMM_WORLD))
      != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::raw_bcast(,"
          << nbyte << "," << from << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (debug_) {
      cerr << scprintf("Node %d done with bcast", me()) << endl;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
